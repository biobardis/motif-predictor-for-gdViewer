from __future__ import annotations

import argparse

from . import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="motif_predict.py",
        description=(
            "Merge overlapping gene +/- flank regions, scan them with FIMO in parallel, "
            "and prepare genomic-coordinate motif hits for gdViewer."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    required = parser.add_argument_group("required arguments")
    required.add_argument("--meme", required=True, help="Input motif file in MEME format.")
    required.add_argument("--fasta", required=True, help="Reference genome FASTA file.")
    required.add_argument("--gff", required=True, help="Genome annotation in GFF3 format.")
    required.add_argument(
        "--output-dir",
        required=True,
        help="Output directory. It is created automatically when absent.",
    )

    required.add_argument(
        "--gene-list",
        required=True,
        help=(
            "Gene selection. Use 'all' to analyze all genes in the GFF3 file, "
            "or provide a text file containing one gene ID per line."
        ),
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=2000,
        help="Number of flanking bases added to both sides of each gene.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Maximum number of FIMO processes run in parallel.",
    )
    parser.add_argument(
        "--batches",
        type=int,
        default=None,
        help="Number of FASTA batches. If omitted, the value of --threads is used.",
    )
    parser.add_argument(
        "--thresh",
        default="1e-4",
        help="FIMO p-value threshold passed to --thresh.",
    )
    parser.add_argument(
        "--output",
        default="motif_predictions.tsv",
        help=(
            "Base name for the final merged TSV. The FASTA sample name is inserted "
            "before .tsv automatically, e.g. motif_predictions.FZX.tsv."
        ),
    )
    parser.add_argument(
        "--fimo-bin",
        default="fimo",
        help="FIMO executable name or path.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Debug mode: only process the first --debug-n-genes genes.",
    )
    parser.add_argument(
        "--debug-n-genes",
        type=int,
        default=100,
        help="Number of genes processed in debug mode (default: 100).",
    )

    return parser


def parse_args(argv=None) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.flank < 0:
        parser.error("--flank must be >= 0")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.batches is not None and args.batches < 1:
        parser.error("--batches must be >= 1")
    if args.debug_n_genes < 1:
        parser.error("--debug-n-genes must be >= 1")
    try:
        thresh = float(args.thresh)
    except ValueError:
        parser.error("--thresh must be a numeric p-value")
    if not (0 < thresh <= 1):
        parser.error("--thresh must be > 0 and <= 1")

    return args