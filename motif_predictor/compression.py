from __future__ import annotations

import heapq
import os
import tempfile
import time
from dataclasses import dataclass


AUTO_BGZIP_THRESHOLD_BYTES = 100 * 1024 * 1024


def human_num(n: int) -> str:
    return f"{n:,}"


def human_time(seconds: float) -> str:
    seconds = int(seconds)
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    if h > 0:
        return f"{h}h {m}m {s}s"
    if m > 0:
        return f"{m}m {s}s"
    return f"{s}s"


def human_size(num_bytes: int) -> str:
    value = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return f"{num_bytes} B"


@dataclass(frozen=True)
class PackageResult:
    tsv_path: str
    gz_path: str | None
    index_path: str | None
    packaged: bool
    size_bytes: int


def sort_large_tsv(
    input_file: str,
    output_file: str,
    *,
    chunk_lines: int = 500_000,
    progress_every: int = 200_000,
    tmp_dir: str | None = None,
    seq_col: int = 0,
    start_col: int = 1,
) -> str:
    """
    Externally sort a large TSV by chromosome/contig and numeric start.

    Lines beginning with ``#`` are treated as header/meta lines and kept at the
    beginning of the output. Data are sorted in chunks and merged with heapq,
    so the entire file does not need to be loaded into memory at once.
    """
    if chunk_lines < 1:
        raise ValueError("chunk_lines must be >= 1")
    if progress_every < 1:
        raise ValueError("progress_every must be >= 1")

    temp_files: list[str] = []
    header_lines: list[str] = []
    data_lines = 0
    chunk_idx = 0
    t0 = time.time()

    if tmp_dir is not None:
        os.makedirs(tmp_dir, exist_ok=True)

    def sort_key(line: str) -> tuple[str, int]:
        cols = line.rstrip("\n").split("\t")
        needed = max(seq_col, start_col)
        if len(cols) <= needed:
            raise ValueError(
                f"TSV row has only {len(cols)} column(s), but sorting requires "
                f"column index {needed}: {line[:200]!r}"
            )
        chrom = cols[seq_col]
        try:
            start = int(cols[start_col])
        except ValueError as exc:
            raise ValueError(
                f"Non-integer start coordinate {cols[start_col]!r} in row: "
                f"{line[:200]!r}"
            ) from exc
        return chrom, start

    def write_chunk(chunk: list[str], index: int) -> str:
        chunk.sort(key=sort_key)
        tf = tempfile.NamedTemporaryFile(
            delete=False,
            mode="w",
            encoding="utf-8",
            prefix=f"chunk_{index}_",
            suffix=".tsv",
            dir=tmp_dir,
        )
        try:
            tf.writelines(chunk)
        finally:
            tf.close()
        return tf.name

    try:
        chunk: list[str] = []
        with open(input_file, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                if line.startswith("#"):
                    header_lines.append(line)
                    continue

                chunk.append(line)
                data_lines += 1

                if data_lines % progress_every == 0:
                    print(
                        f"[PACK] read {human_num(data_lines)} data lines; "
                        f"created {chunk_idx} temp chunks; "
                        f"elapsed {human_time(time.time() - t0)}",
                        flush=True,
                    )

                if len(chunk) >= chunk_lines:
                    chunk_idx += 1
                    tmp_path = write_chunk(chunk, chunk_idx)
                    temp_files.append(tmp_path)
                    print(f"[PACK] wrote temp chunk {chunk_idx}: {tmp_path}", flush=True)
                    chunk = []

            if chunk:
                chunk_idx += 1
                tmp_path = write_chunk(chunk, chunk_idx)
                temp_files.append(tmp_path)
                print(f"[PACK] wrote temp chunk {chunk_idx}: {tmp_path}", flush=True)

        print(
            f"[PACK] Finished splitting: {human_num(data_lines)} data lines, "
            f"{chunk_idx} temp chunk(s)",
            flush=True,
        )

        # An empty result should normally never reach this function, but keep the
        # header if it does.
        if not temp_files:
            with open(output_file, "w", encoding="utf-8") as out:
                out.writelines(header_lines)
            return output_file

        def file_iter(path: str):
            with open(path, "r", encoding="utf-8") as handle:
                yield from handle

        def keyed_iter(iterator):
            for row in iterator:
                yield sort_key(row), row

        print(f"[PACK] Merging {len(temp_files)} temp chunk(s)...", flush=True)
        with open(output_file, "w", encoding="utf-8") as out:
            if header_lines:
                out.writelines(header_lines)

            iterators = [file_iter(path) for path in temp_files]
            written = 0
            for _, line in heapq.merge(*(keyed_iter(it) for it in iterators)):
                out.write(line)
                written += 1
                if written % progress_every == 0:
                    print(
                        f"[PACK] merged {human_num(written)} lines; "
                        f"elapsed {human_time(time.time() - t0)}",
                        flush=True,
                    )

        print(f"[PACK] External sort done: {output_file}", flush=True)
        return output_file

    finally:
        for path in temp_files:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass


def package_large_tsv_if_needed(
    tsv_path: str,
    *,
    threshold_bytes: int = AUTO_BGZIP_THRESHOLD_BYTES,
    chunk_lines: int = 500_000,
    progress_every: int = 200_000,
    tmp_dir: str | None = None,
) -> PackageResult:
    """
    If ``tsv_path`` is larger than the threshold, sort it, bgzip-compress it,
    and build a TBI index using columns #chr/start/stop (0/1/2).

    The motif output uses 1-based inclusive genomic coordinates, so Tabix is
    built with ``zerobased=False``. The uncompressed TSV is intentionally kept.
    """
    size_bytes = os.path.getsize(tsv_path)
    print(
        f"[PACK] Final TSV size: {human_size(size_bytes)}; auto-package threshold: "
        f"{human_size(threshold_bytes)}",
        flush=True,
    )

    if size_bytes <= threshold_bytes:
        print("[PACK] TSV does not exceed the auto-package threshold; keeping plain TSV only.", flush=True)
        return PackageResult(
            tsv_path=tsv_path,
            gz_path=None,
            index_path=None,
            packaged=False,
            size_bytes=size_bytes,
        )

    print("[PACK] TSV exceeds the auto-package threshold; sorting, bgzip-compressing, and indexing...", flush=True)

    # Sort to a temporary file first, then atomically replace the plain TSV.
    sort_parent = tmp_dir or os.path.dirname(os.path.abspath(tsv_path))
    os.makedirs(sort_parent, exist_ok=True)
    sorted_tmp = tempfile.NamedTemporaryFile(
        delete=False,
        prefix=".motif_sorted_",
        suffix=".tsv",
        dir=sort_parent,
    )
    sorted_tmp_path = sorted_tmp.name
    sorted_tmp.close()

    try:
        sort_large_tsv(
            tsv_path,
            sorted_tmp_path,
            chunk_lines=chunk_lines,
            progress_every=progress_every,
            tmp_dir=tmp_dir,
            seq_col=0,
            start_col=1,
        )
        os.replace(sorted_tmp_path, tsv_path)
    finally:
        if os.path.exists(sorted_tmp_path):
            try:
                os.unlink(sorted_tmp_path)
            except OSError:
                pass

    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError(
            "Output TSV exceeds 100 MB and automatic bgzip/Tabix packaging is required, "
            "but pysam is not installed. Install dependencies with: "
            "pip install -r requirements.txt"
        ) from exc

    gz_path = f"{tsv_path}.gz"
    index_path = f"{gz_path}.tbi"

    if os.path.exists(gz_path) or os.path.exists(index_path):
        raise FileExistsError(
            f"Automatic package target already exists: {gz_path} or {index_path}. "
            "The program will not overwrite existing results."
        )

    compress_t0 = time.time()
    print(f"[PACK] Compressing with bgzip -> {gz_path}", flush=True)
    pysam.tabix_compress(tsv_path, gz_path, force=False)
    print(
        f"[PACK] bgzip compression done in {human_time(time.time() - compress_t0)}",
        flush=True,
    )

    index_t0 = time.time()
    print(f"[PACK] Building Tabix index -> {index_path}", flush=True)
    pysam.tabix_index(
        gz_path,
        seq_col=0,
        start_col=1,
        end_col=2,
        meta_char="#",
        force=False,
        preset=None,
        zerobased=False,
    )
    print(
        f"[PACK] Tabix index done in {human_time(time.time() - index_t0)}",
        flush=True,
    )
    print(
        f"[PACK] Packaged outputs:\n"
        f" - {tsv_path}\n"
        f" - {gz_path}\n"
        f" - {index_path}",
        flush=True,
    )

    return PackageResult(
        tsv_path=tsv_path,
        gz_path=gz_path,
        index_path=index_path,
        packaged=True,
        size_bytes=size_bytes,
    )