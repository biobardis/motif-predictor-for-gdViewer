from __future__ import annotations

import os
import shutil
import tempfile
import time

from .compression import package_large_tsv_if_needed
from .fasta import write_merged_regions_to_split_fastas
from .fimo import run_fimo_batches_parallel
from .gff3 import extract_gene_regions, load_gene_list, merge_gene_regions
from .meme import parse_meme_metadata
from .output import (
    add_fasta_tag_to_output,
    cleanup_work_dir,
    make_unique_path,
    merge_batch_fimo_results,
    resolve_output_path,
)


def _validate_input_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} does not exist: {path}")


def _resolve_executable(fimo_bin: str) -> str:
    if os.path.sep in fimo_bin:
        if not os.path.isfile(fimo_bin):
            raise FileNotFoundError(f"FIMO executable does not exist: {fimo_bin}")
        return fimo_bin
    resolved = shutil.which(fimo_bin)
    if not resolved:
        raise RuntimeError(
            f"FIMO executable {fimo_bin!r} was not found in PATH. "
            "Install MEME Suite or use --fimo-bin."
        )
    return resolved


def run_pipeline(
    *,
    meme_file: str,
    fasta_file: str,
    gff3_file: str,
    output_dir: str,
    gene_list_file: str,
    flank: int = 2000,
    num_threads: int = 8,
    num_batches: int | None = None,
    thresh: str = "1e-4",
    output_file: str ="motif_predictions.motif.tsv",
    fimo_bin: str = "fimo",
    debug: bool = False,
    debug_n_genes: int = 100,
) -> str:
    start_time = time.time()

    for path, label in (
        (meme_file, "MEME file"),
        (fasta_file, "FASTA file"),
        (gff3_file, "GFF3 file"),
    ):
        _validate_input_file(path, label)
    if gene_list_file.strip().lower() != "all":
        _validate_input_file(gene_list_file, "Gene list")

    fimo_bin = _resolve_executable(fimo_bin)
    output_dir = os.path.abspath(os.path.expanduser(output_dir))
    os.makedirs(output_dir, exist_ok=True)

    meme_metadata = parse_meme_metadata(meme_file)
    print(
        f"[MEME] Loaded {meme_metadata.motif_count} motif(s); widths="
        + ",".join(str(x) for x in meme_metadata.motif_widths),
        flush=True,
    )

    gene_ids = load_gene_list(gene_list_file)
    gene_regions = extract_gene_regions(gff3_file, flank=flank, gene_ids=gene_ids)

    if debug:
        gene_regions = gene_regions[:debug_n_genes]
        print(
            f"[DEBUG] Debug mode enabled: processing the first "
            f"{debug_n_genes} selected gene(s).",
            flush=True,
        )

    if not gene_regions:
        raise RuntimeError("No gene regions were selected from the GFF3 file")

    print(f"[MAIN] Selected {len(gene_regions):,} gene +/- flank region(s)", flush=True)

    original_bp = sum(region.length for region in gene_regions)
    merged_regions = merge_gene_regions(gene_regions)
    merged_bp = sum(region.length for region in merged_regions)
    reduction_bp = max(0, original_bp - merged_bp)
    reduction_pct = (reduction_bp / original_bp * 100.0) if original_bp else 0.0

    print(
        f"[MERGE] {len(gene_regions):,} gene region(s) -> "
        f"{len(merged_regions):,} non-overlapping merged region(s)",
        flush=True,
    )
    print(
        f"[MERGE] Raw scan length {original_bp:,} bp -> {merged_bp:,} bp "
        f"({reduction_pct:.2f}% fewer repeated bases)",
        flush=True,
    )

    num_batches = num_batches or num_threads
    num_batches = min(num_batches, len(merged_regions))

    # Each run gets its own private workspace. Existing user files in the output
    # directory are never cleared or overwritten.
    work_dir = tempfile.mkdtemp(prefix=".motif_predict_work_", dir=output_dir)

    merged_path: str | None = None
    package_result = None
    try:
        fasta_batch_dir = os.path.join(work_dir, "fasta_batches")
        fasta_result = write_merged_regions_to_split_fastas(
            fasta_file=fasta_file,
            regions=merged_regions,
            output_dir=fasta_batch_dir,
            num_batches=num_batches,
        )

        fimo_tsvs = run_fimo_batches_parallel(
            meme_file=meme_file,
            fasta_paths=fasta_result.fasta_paths,
            output_dir=work_dir,
            failure_log_dir=output_dir,
            num_threads=num_threads,
            thresh=thresh,
            fimo_bin=fimo_bin,
        )

        tagged_output_file = add_fasta_tag_to_output(output_file, fasta_file)
        requested_path = resolve_output_path(output_dir, tagged_output_file)
        merged_path = make_unique_path(requested_path)
        if os.path.abspath(merged_path) != os.path.abspath(requested_path):
            print(f"[MAIN] Output file already exists; using: {merged_path}", flush=True)

        merged_path = merge_batch_fimo_results(
            fimo_tsvs=fimo_tsvs,
            merged_path=merged_path,
            region_map=fasta_result.region_map,
        )
        if not merged_path:
            raise RuntimeError("FIMO completed, but no valid motif hits were found to merge")

        package_result = package_large_tsv_if_needed(
            merged_path,
            tmp_dir=os.path.join(work_dir, "sort_tmp"),
        )

    finally:
        cleanup_work_dir(work_dir)

    elapsed = int(time.time() - start_time)
    h, rem = divmod(elapsed, 3600)
    m, s = divmod(rem, 60)
    print(f"[MAIN] Completed in {h:02d}:{m:02d}:{s:02d}", flush=True)
    print(f"[MAIN] Output: {merged_path}", flush=True)
    if package_result is not None and package_result.packaged:
        print(f"[MAIN] bgzip: {package_result.gz_path}", flush=True)
        print(f"[MAIN] Tabix: {package_result.index_path}", flush=True)
    return merged_path