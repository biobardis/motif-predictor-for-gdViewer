from __future__ import annotations

import multiprocessing
import os
import subprocess

from .output import make_unique_path


def prefixed_print(prefix: str, line: str) -> None:
    line = line.rstrip("\n")
    if line:
        print(f"[{prefix}] {line}", flush=True)


def run_one_fimo_batch(args: tuple) -> dict:
    """
    Run one standard FIMO process on one multi-FASTA batch.

    FIMO writes its normal output directory with --oc. ``--no-pgc`` is used
    so FIMO does not reinterpret the merged-region FASTA headers as genomic
    coordinates; this keeps sequence_name identical to the region map key.
    The q-value column in fimo.tsv is calculated by FIMO itself and is preserved
    unchanged by the downstream merge step.
    """
    meme_file, input_fasta, batch_output_dir, thresh, fimo_bin = args
    batch_name = os.path.basename(input_fasta)
    fimo_tsv = os.path.join(batch_output_dir, "fimo.tsv")

    try:
        # Do not pre-create batch_output_dir: FIMO --oc creates it itself and
        # expects the target directory not to exist.
        cmd = [
            fimo_bin,
            "--oc",
            batch_output_dir,
            "--no-pgc",
            "--thresh",
            str(thresh),
            meme_file,
            input_fasta,
        ]
        prefixed_print(batch_name, f"Running: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        if result.stdout:
            for line in result.stdout.splitlines():
                prefixed_print(batch_name, line)
        if result.stderr:
            for line in result.stderr.splitlines():
                prefixed_print(batch_name, line)

        if result.returncode != 0:
            return {
                "ok": False,
                "input_fasta": input_fasta,
                "batch_output_dir": batch_output_dir,
                "fimo_tsv": fimo_tsv,
                "cmd": " ".join(cmd),
                "stdout": result.stdout,
                "stderr": result.stderr,
                "message": f"FIMO failed with return code {result.returncode}",
            }

        if not os.path.exists(fimo_tsv):
            return {
                "ok": False,
                "input_fasta": input_fasta,
                "batch_output_dir": batch_output_dir,
                "fimo_tsv": fimo_tsv,
                "cmd": " ".join(cmd),
                "stdout": result.stdout,
                "stderr": result.stderr,
                "message": "FIMO completed but fimo.tsv was not created",
            }

        return {
            "ok": True,
            "input_fasta": input_fasta,
            "batch_output_dir": batch_output_dir,
            "fimo_tsv": fimo_tsv,
            "cmd": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "message": "success",
        }
    except Exception as exc:
        return {
            "ok": False,
            "input_fasta": input_fasta,
            "batch_output_dir": batch_output_dir,
            "fimo_tsv": fimo_tsv,
            "cmd": "",
            "stdout": "",
            "stderr": "",
            "message": str(exc),
        }


def run_fimo_batches_parallel(
    meme_file: str,
    fasta_paths: list[str] | tuple[str, ...],
    output_dir: str,
    failure_log_dir: str,
    num_threads: int,
    thresh: str,
    fimo_bin: str,
) -> list[str]:
    batch_output_root = os.path.join(output_dir, "fimo_batches")
    os.makedirs(batch_output_root, exist_ok=True)

    args_list = []
    for i, fasta_path in enumerate(fasta_paths, 1):
        args_list.append(
            (
                meme_file,
                fasta_path,
                os.path.join(batch_output_root, f"batch_{i:03d}"),
                thresh,
                fimo_bin,
            )
        )

    process_count = min(num_threads, len(args_list))
    print(
        f"[FIMO] Running {len(args_list)} multi-FASTA batch(es) with "
        f"up to {process_count} FIMO process(es)",
        flush=True,
    )

    fimo_tsvs: list[str] = []
    failed: list[dict] = []
    with multiprocessing.Pool(processes=process_count) as pool:
        for i, result in enumerate(pool.imap_unordered(run_one_fimo_batch, args_list), 1):
            fasta_name = os.path.basename(result["input_fasta"])
            if result["ok"]:
                print(f"[FIMO] [{i}/{len(args_list)}] {fasta_name}: success", flush=True)
                fimo_tsvs.append(result["fimo_tsv"])
            else:
                print(f"[FIMO] [{i}/{len(args_list)}] {fasta_name}: failed", flush=True)
                failed.append(result)

    if failed:
        failed_log = make_unique_path(os.path.join(failure_log_dir, "failed_fimo_batches.log"))
        with open(failed_log, "w", encoding="utf-8") as out:
            for item in failed:
                out.write("=" * 80 + "\n")
                for key in ("input_fasta", "batch_output_dir", "cmd", "message"):
                    out.write(f"{key}: {item.get(key, '')}\n")
                out.write(f"stdout:\n{item.get('stdout', '')}\n")
                out.write(f"stderr:\n{item.get('stderr', '')}\n")
        print(f"[WARNING] {len(failed)} FIMO batch(es) failed; see {failed_log}", flush=True)

    if not fimo_tsvs:
        raise RuntimeError("All FIMO batches failed")
    return sorted(fimo_tsvs)