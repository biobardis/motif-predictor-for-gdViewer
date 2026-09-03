from __future__ import annotations

import heapq
import os
from dataclasses import dataclass

from Bio import SeqIO

from .gff3 import GeneRegion, MergedRegion


@dataclass(frozen=True)
class ScanRegion:
    sequence_name: str
    region_id: str
    chrom: str
    start0: int
    end0: int
    gene_regions: tuple[GeneRegion, ...]

    @property
    def length(self) -> int:
        return max(0, self.end0 - self.start0)

    @property
    def gene_ids(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(region.gene_id for region in self.gene_regions))


@dataclass(frozen=True)
class FastaBatchResult:
    fasta_paths: tuple[str, ...]
    region_map: dict[str, ScanRegion]
    total_bases: int
    batch_bases: tuple[int, ...]


def load_genome_fasta(fasta_file: str) -> dict:
    print(f"[FASTA] Loading reference genome: {fasta_file}", flush=True)
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    if not genome:
        raise RuntimeError("Reference FASTA contains no sequences")
    print(f"[FASTA] Loaded {len(genome)} sequence(s)", flush=True)
    return genome


def _prepare_regions(
    genome: dict,
    regions: list[MergedRegion],
) -> list[tuple[MergedRegion, int, int, int]]:
    prepared: list[tuple[MergedRegion, int, int, int]] = []
    missing_contigs: set[str] = set()

    for region in regions:
        record = genome.get(region.chrom)
        if record is None:
            missing_contigs.add(region.chrom)
            continue

        chrom_len = len(record.seq)
        start0 = min(max(0, region.start0), chrom_len)
        end0 = min(max(start0, region.end0), chrom_len)
        length = end0 - start0
        if length <= 0:
            continue
        prepared.append((region, start0, end0, length))

    if missing_contigs:
        preview = ", ".join(sorted(missing_contigs)[:20])
        suffix = " ..." if len(missing_contigs) > 20 else ""
        print(
            f"[WARNING] {len(missing_contigs)} GFF3 contig(s) were not found in FASTA: "
            f"{preview}{suffix}",
            flush=True,
        )

    return prepared


def _balanced_batches(
    prepared: list[tuple[MergedRegion, int, int, int]],
    num_batches: int,
) -> list[list[tuple[MergedRegion, int, int, int]]]:
    """Greedily balance batches by total sequence length rather than record count."""
    num_batches = min(num_batches, len(prepared))
    batches: list[list[tuple[MergedRegion, int, int, int]]] = [
        [] for _ in range(num_batches)
    ]
    heap: list[tuple[int, int]] = [(0, i) for i in range(num_batches)]
    heapq.heapify(heap)

    for item in sorted(prepared, key=lambda x: x[3], reverse=True):
        load, batch_idx = heapq.heappop(heap)
        batches[batch_idx].append(item)
        heapq.heappush(heap, (load + item[3], batch_idx))

    return batches


def write_merged_regions_to_split_fastas(
    fasta_file: str,
    regions: list[MergedRegion],
    output_dir: str,
    num_batches: int,
) -> FastaBatchResult:
    """Write non-overlapping merged regions into length-balanced multi-FASTA batches."""
    if num_batches < 1:
        raise ValueError("num_batches must be >= 1")

    os.makedirs(output_dir, exist_ok=True)
    genome = load_genome_fasta(fasta_file)
    prepared = _prepare_regions(genome, regions)
    if not prepared:
        raise RuntimeError("No valid merged genomic regions could be extracted from FASTA")

    batches = _balanced_batches(prepared, num_batches)
    fasta_paths: list[str] = []
    region_map: dict[str, ScanRegion] = {}
    batch_bases: list[int] = []
    total_bases = 0

    for batch_idx, batch in enumerate(batches, 1):
        if not batch:
            continue
        path = os.path.join(output_dir, f"regions_batch_{batch_idx:03d}.fa")
        current_bases = 0

        with open(path, "w", encoding="utf-8") as out:
            # Restore genomic order inside each batch for easier manual inspection.
            for region, start0, end0, _ in sorted(
                batch,
                key=lambda x: (x[0].chrom, x[1], x[2]),
            ):
                seq = str(genome[region.chrom].seq[start0:end0]).upper()
                if not seq:
                    continue

                # Internal FIMO scan name. --no-pgc keeps it unchanged in fimo.tsv.
                sequence_name = (
                    f"{region.region_id}|{region.chrom}:{start0 + 1}-{end0}"
                )
                out.write(f">{sequence_name}\n{seq}\n")

                scan_region = ScanRegion(
                    sequence_name=sequence_name,
                    region_id=region.region_id,
                    chrom=region.chrom,
                    start0=start0,
                    end0=end0,
                    gene_regions=region.gene_regions,
                )
                region_map[sequence_name] = scan_region

                seq_len = len(seq)
                current_bases += seq_len
                total_bases += seq_len

        if current_bases > 0:
            fasta_paths.append(path)
            batch_bases.append(current_bases)
        else:
            try:
                os.remove(path)
            except OSError:
                pass

    if not fasta_paths:
        raise RuntimeError("No valid merged regions were written to FASTA batches")

    print(
        f"[FASTA] Wrote {len(region_map)} merged region(s) into "
        f"{len(fasta_paths)} multi-FASTA batch(es)",
        flush=True,
    )
    print(f"[FASTA] Total sequence scanned: {total_bases:,} bp", flush=True)
    print(
        "[FASTA] Batch sizes: " + ", ".join(f"{x:,} bp" for x in batch_bases),
        flush=True,
    )

    return FastaBatchResult(
        fasta_paths=tuple(fasta_paths),
        region_map=region_map,
        total_bases=total_bases,
        batch_bases=tuple(batch_bases),
    )