from __future__ import annotations

import sys

from .cli import parse_args
from .pipeline import run_pipeline


def main(argv=None) -> None:
    args = parse_args(argv)
    try:
        run_pipeline(
            meme_file=args.meme,
            fasta_file=args.fasta,
            gff3_file=args.gff,
            output_dir=args.output_dir,
            gene_list_file=args.gene_list,
            flank=args.flank,
            num_threads=args.threads,
            num_batches=args.batches,
            thresh=args.thresh,
            output_file=args.output,
            fimo_bin=args.fimo_bin,
            debug=args.debug,
            debug_n_genes=args.debug_n_genes,
        )
    except KeyboardInterrupt:
        print("\n[ERROR] Interrupted by user", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)