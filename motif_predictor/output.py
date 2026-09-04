from __future__ import annotations

import os
import shutil
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError

from .fasta import ScanRegion
from .gff3 import GeneRegion



_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".ffn", ".faa"}
_COMPRESSION_SUFFIXES = {".gz", ".bgz", ".bgzip"}


def fasta_name_tag(fasta_file: str) -> str:
    """Return a concise FASTA-derived tag such as FZX from FZX.fa or FZX.fa.gz."""
    name = Path(fasta_file).name
    path_obj = Path(name)

    if path_obj.suffix.lower() in _COMPRESSION_SUFFIXES:
        name = path_obj.stem
        path_obj = Path(name)

    if path_obj.suffix.lower() in _FASTA_SUFFIXES:
        name = path_obj.stem

    name = name.strip() or "fasta"
    # Keep dots/dashes that are often meaningful in assembly names, but avoid
    # whitespace and path-hostile separators.
    name = "_".join(name.split())
    name = name.replace("/", "_").replace("\\", "_")
    return name


def add_fasta_tag_to_output(output_value: str, fasta_file: str) -> str:
    """Insert ``.<FASTA_TAG>`` immediately before the final .tsv suffix."""
    tag = fasta_name_tag(fasta_file)
    value = os.path.expanduser(output_value)
    path_obj = Path(value)
    name = path_obj.name

    if name.lower().endswith(".motif.tsv"):
        stem = name[:-10]  # remove ".motif.tsv"
        if stem.endswith(f".{tag}"):
            tagged_name = name
        else:
            tagged_name = f"{stem}.{tag}.motif.tsv"

    elif name.lower().endswith(".tsv"):
        stem = name[:-4]
        if stem.endswith(f".{tag}"):
            tagged_name = name
        else:
            tagged_name = f"{stem}.{tag}.tsv"

    else:
        if name.endswith(f".{tag}"):
            tagged_name = f"{name}.tsv"
        else:
            tagged_name = f"{name}.{tag}.tsv"
    return str(path_obj.with_name(tagged_name))

def resolve_output_path(output_dir: str, value: str) -> str:
    value = os.path.expanduser(value)
    if os.path.isabs(value) or os.path.dirname(value):
        return value
    return os.path.join(output_dir, value)


def make_unique_path(path: str) -> str:
    """Return an unused path, appending _2, _3, ... when needed.

    For TSV outputs, the corresponding bgzip/Tabix sidecars are treated as
    part of the same result set. Thus a candidate is considered occupied if
    any of ``<name>.tsv``, ``<name>.tsv.gz``, or ``<name>.tsv.gz.tbi`` exists.
    Other file types keep the original single-path collision check.
    """
    path_obj = Path(path)

    def candidate_is_used(candidate: Path) -> bool:
        if candidate.exists():
            return True
        if candidate.name.lower().endswith(".tsv"):
            gz_path = Path(f"{candidate}.gz")
            tbi_path = Path(f"{candidate}.gz.tbi")
            return gz_path.exists() or tbi_path.exists()
        return False

    if not candidate_is_used(path_obj):
        return str(path_obj)

    lower_name = path_obj.name.lower()
    if lower_name.endswith(".tsv.gz"):
        extension = path_obj.name[-7:]
        base_name = path_obj.name[:-7]
    else:
        extension = path_obj.suffix
        base_name = path_obj.stem if extension else path_obj.name

    n = 2
    while True:
        candidate = path_obj.with_name(f"{base_name}_{n}{extension}")
        if not candidate_is_used(candidate):
            return str(candidate)
        n += 1


def _load_fimo_tsv(path: str) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
    except EmptyDataError:
        return None
    except Exception as exc:
        print(f"[WARNING] Skipping unreadable FIMO file {path}: {exc}", flush=True)
        return None

    if df.empty:
        return None

    required = [
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p-value",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"[WARNING] Skipping {path}; missing columns: {missing}", flush=True)
        return None

    return df


def _unique_gene_regions(regions: tuple[GeneRegion, ...]) -> list[GeneRegion]:
    """Keep the first occurrence of each gene-level map name."""
    out: list[GeneRegion] = []
    seen: set[str] = set()
    for region in regions:
        key = region.sequence_name
        if key in seen:
            continue
        seen.add(key)
        out.append(region)
    return out


def _map_hit_to_gene_regions(
    scan_region: ScanRegion,
    hit_start1: int,
    hit_stop1: int,
) -> list[GeneRegion]:
    """
    Map one genomic motif hit back to gene +/- flank regions.

    Coordinates passed here are 1-based inclusive genomic coordinates.

    Primary rule:
      return every gene extended region that fully contains the motif hit.

    Boundary fallback:
      if no region fully contains the motif (possible when a hit crosses the
      boundary between touching/overlapping regions after merging), return all
      gene extended regions that overlap the motif. This guarantees that a FIMO
      hit is not silently lost during gene-level remapping.
    """
    hit_start0 = hit_start1 - 1
    hit_end0 = hit_stop1  # 1-based inclusive -> 0-based half-open

    candidates = _unique_gene_regions(scan_region.gene_regions)

    contained = [
        region
        for region in candidates
        if hit_start0 >= region.start0 and hit_end0 <= region.end0
    ]
    if contained:
        return contained

    overlapping = [
        region
        for region in candidates
        if hit_start0 < region.end0 and hit_end0 > region.start0
    ]
    return overlapping


def merge_batch_fimo_results(
    fimo_tsvs: list[str],
    merged_path: str,
    region_map: dict[str, ScanRegion],
) -> str | None:
    """
    Merge FIMO batch TSVs, convert coordinates, and restore gene-level map names.

    FIMO still scans merged genomic regions for speed. After each hit is
    converted to genomic coordinates, it is mapped back to the original
    gene +/- flank region(s). The final ``sequence_name`` therefore has the
    form ``<gene_id>_extended|<chrom>``.

    If a motif belongs to multiple gene extended regions, one output row is
    emitted for each gene. Statistical columns produced by FIMO, including
    p-value and q-value, are copied without recalculation or adjustment.
    """
    expanded_rows: list[dict] = []
    fimo_columns: list[str] = []
    unmapped_count = 0

    for fimo_tsv in sorted(fimo_tsvs):
        df = _load_fimo_tsv(fimo_tsv)
        if df is None:
            continue

        if not fimo_columns:
            fimo_columns = list(df.columns)
        else:
            for col in df.columns:
                if col not in fimo_columns:
                    fimo_columns.append(col)

        sequence_names = df["sequence_name"].astype(str)
        unknown = sorted(set(sequence_names) - set(region_map))
        if unknown:
            preview = ", ".join(unknown[:5])
            raise RuntimeError(
                "FIMO output contains sequence_name values not found in the "
                f"scan-region map: {preview}"
            )

        for row_index, row in df.iterrows():
            internal_name = str(row["sequence_name"])
            scan_region = region_map[internal_name]

            try:
                rel_start = int(row["start"])
                rel_stop = int(row["stop"])
            except (TypeError, ValueError) as exc:
                raise RuntimeError(
                    f"Invalid FIMO start/stop in {fimo_tsv}, row {row_index + 1}"
                ) from exc

            rel_lo = min(rel_start, rel_stop)
            rel_hi = max(rel_start, rel_stop)

            # FIMO positions are 1-based within the merged FASTA record.
            genomic_start1 = scan_region.start0 + rel_lo
            genomic_stop1 = scan_region.start0 + rel_hi

            gene_regions = _map_hit_to_gene_regions(
                scan_region,
                genomic_start1,
                genomic_stop1,
            )

            if not gene_regions:
                # This should be extremely rare because each merged scan region
                # is the union of its source gene +/- flank regions. Do not drop
                # the FIMO hit: retain an explicit merged-region fallback name.
                unmapped_count += 1
                out_row = row.to_dict()
                out_row["#chr"] = scan_region.chrom
                out_row["start"] = genomic_start1
                out_row["stop"] = genomic_stop1
                out_row["family"] = "."
                out_row["sequence_name"] = internal_name
                expanded_rows.append(out_row)
                continue

            for gene_region in gene_regions:
                out_row = row.to_dict()
                out_row["#chr"] = scan_region.chrom
                out_row["start"] = genomic_start1
                out_row["stop"] = genomic_stop1
                out_row["family"] = "."
                out_row["sequence_name"] = gene_region.sequence_name
                expanded_rows.append(out_row)

    if not expanded_rows:
        return None

    merged_df = pd.DataFrame(expanded_rows)

    # Preferred gdViewer order. FIMO p-value/q-value values are not recalculated.
    preferred_cols = [
        "#chr",
        "start",
        "stop",
        "motif_id",
        "strand",
        "motif_alt_id",
        "family",
        "score",
        "p-value",
        "q-value",
        "matched_sequence",
        "sequence_name",
    ]
    output_cols = [c for c in preferred_cols if c in merged_df.columns]

    # Preserve additional columns emitted by the installed FIMO version.
    for col in fimo_columns:
        if col not in output_cols and col not in {"#chr", "family"}:
            output_cols.append(col)
    for col in merged_df.columns:
        if col not in output_cols:
            output_cols.append(col)

    merged_df = merged_df[output_cols].copy()
    sort_cols = [c for c in ["#chr", "start", "stop", "sequence_name"] if c in merged_df.columns]
    if sort_cols:
        merged_df.sort_values(by=sort_cols, kind="mergesort", inplace=True)

    os.makedirs(os.path.dirname(os.path.abspath(merged_path)), exist_ok=True)
    merged_df.to_csv(merged_path, sep="\t", index=False)

    qvalue_note = (
        "FIMO q-values preserved unchanged"
        if "q-value" in merged_df.columns
        else "no q-value column was emitted by FIMO"
    )
    print(
        f"[OUTPUT] Wrote {len(merged_df):,} gene-mapped motif row(s); {qvalue_note}",
        flush=True,
    )
    if unmapped_count:
        print(
            f"[WARNING] {unmapped_count:,} motif hit(s) could not be assigned to a "
            "gene extended region and retained the internal merged-region name.",
            flush=True,
        )

    return merged_path


def cleanup_work_dir(work_dir: str) -> None:
    """Remove only the temporary workspace created by the current run."""
    if os.path.isdir(work_dir):
        shutil.rmtree(work_dir)