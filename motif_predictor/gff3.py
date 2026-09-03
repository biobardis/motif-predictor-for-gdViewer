from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


@dataclass(frozen=True)
class GeneRegion:
    """One selected gene plus its flanking sequence, in 0-based half-open coordinates."""

    chrom: str
    start0: int
    end0: int
    gene_id: str
    strand: str

    @property
    def length(self) -> int:
        return max(0, self.end0 - self.start0)

    @property
    def sequence_name(self) -> str:
        """Gene-level name written to the final gdViewer TSV."""
        return f"{self.gene_id}_extended|{self.chrom}"


@dataclass(frozen=True)
class MergedRegion:
    """Union of one or more overlapping/touching GeneRegion objects."""

    chrom: str
    start0: int
    end0: int
    region_id: str
    gene_regions: tuple[GeneRegion, ...]

    @property
    def length(self) -> int:
        return max(0, self.end0 - self.start0)

    @property
    def gene_ids(self) -> tuple[str, ...]:
        """Unique gene IDs represented by this merged scan region."""
        return tuple(dict.fromkeys(region.gene_id for region in self.gene_regions))


def parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = unquote(value.strip())
        else:
            attrs[item] = ""
    return attrs


def load_gene_list(path: str) -> list[str] | None:
    """Load selected gene IDs, or return None when the literal 'all' is used."""
    if path.strip().lower() == "all":
        return None

    out: list[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.append(line.split()[0])
    gene_ids = list(dict.fromkeys(out))
    if not gene_ids:
        raise ValueError(
            "Gene list is empty. Provide at least one gene ID, or use --gene-list all "
            "to analyze all genes."
        )
    return gene_ids


def _gene_id(attrs: dict[str, str]) -> str | None:
    return attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")


def extract_gene_regions(
    gff3_file: str,
    flank: int = 2000,
    gene_ids: list[str] | None = None,
) -> list[GeneRegion]:
    """Read gene/pseudogene coordinates and return gene +/- flank regions."""
    if Path(gff3_file).suffix.lower() not in {".gff3", ".gff"}:
        raise ValueError("--gff must be a GFF3/GFF file")

    wanted = set(gene_ids) if gene_ids else None
    found: set[str] = set()
    regions: list[GeneRegion] = []

    with open(gff3_file, "r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GFF3 line {line_no}: expected 9 columns")

            chrom, _, feature_type, start, end, _, strand, _, attr_text = fields
            if feature_type.lower() not in {"gene", "pseudogene"}:
                continue

            attrs = parse_attributes(attr_text)
            gene_id = _gene_id(attrs)
            if not gene_id:
                continue
            if wanted is not None and gene_id not in wanted:
                continue

            try:
                start1 = int(start)
                end1 = int(end)
            except ValueError as exc:
                raise ValueError(f"Invalid coordinates on GFF3 line {line_no}") from exc

            if start1 < 1 or end1 < start1:
                raise ValueError(f"Invalid coordinates on GFF3 line {line_no}: {start1}-{end1}")

            start0 = max(0, start1 - 1 - flank)
            end0 = end1 + flank
            regions.append(GeneRegion(chrom, start0, end0, gene_id, strand))
            found.add(gene_id)

    if wanted is not None:
        missing = sorted(wanted - found)
        if missing:
            preview = ", ".join(missing[:20])
            suffix = " ..." if len(missing) > 20 else ""
            print(
                f"[WARNING] {len(missing)} gene ID(s) were not found in GFF3: "
                f"{preview}{suffix}",
                flush=True,
            )

    return regions


def merge_gene_regions(regions: list[GeneRegion]) -> list[MergedRegion]:
    """
    Merge overlapping or directly touching gene +/- flank intervals.

    The original GeneRegion objects are retained inside each MergedRegion. This
    lets FIMO scan each merged genomic sequence only once while allowing every
    final motif hit to be mapped back to the gene-level extended region(s) that
    contain or overlap that hit.
    """
    if not regions:
        return []

    chrom_order: list[str] = []
    by_chrom: dict[str, list[GeneRegion]] = {}
    for region in regions:
        if region.chrom not in by_chrom:
            by_chrom[region.chrom] = []
            chrom_order.append(region.chrom)
        by_chrom[region.chrom].append(region)

    raw_merged: list[tuple[str, int, int, tuple[GeneRegion, ...]]] = []

    for chrom in chrom_order:
        chrom_regions = sorted(
            by_chrom[chrom],
            key=lambda r: (r.start0, r.end0, r.gene_id),
        )

        current_start = chrom_regions[0].start0
        current_end = chrom_regions[0].end0
        current_regions: list[GeneRegion] = [chrom_regions[0]]

        for region in chrom_regions[1:]:
            # Half-open intervals overlap when next.start < current.end.
            # Equality means directly touching. Merging touching intervals adds
            # no extra genomic bases and reduces FASTA/FIMO bookkeeping.
            if region.start0 <= current_end:
                current_end = max(current_end, region.end0)
                current_regions.append(region)
            else:
                raw_merged.append(
                    (chrom, current_start, current_end, tuple(current_regions))
                )
                current_start = region.start0
                current_end = region.end0
                current_regions = [region]

        raw_merged.append(
            (chrom, current_start, current_end, tuple(current_regions))
        )

    merged: list[MergedRegion] = []
    for i, (chrom, start0, end0, gene_regions) in enumerate(raw_merged, 1):
        merged.append(
            MergedRegion(
                chrom=chrom,
                start0=start0,
                end0=end0,
                region_id=f"merged_{i:06d}",
                gene_regions=gene_regions,
            )
        )

    return merged