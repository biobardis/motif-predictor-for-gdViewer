# Development Status
>
> `motif-predictor-for-gdViewer` is currently under active development.
>
> The latest development code is available in the [`dev`](../../tree/dev) branch. The `main` branch is not yet intended for general use.

# motif-predictor-for-gdViewer

`motif-predictor-for-gdViewer` is a command-line tool for predicting motif occurrences using FIMO and preparing motif data for visualization in [gdViewer](https://github.com/biobardis/jbrowse-plugin-gdViewer).

The program integrates genome annotations, reference genome sequences, and motif matrices to identify motif occurrences in gene regions and their flanking sequences.

Gene selection must be specified with `--gene-list`. Use `--gene-list all` to analyze all genes in the GFF3 file, or provide a text file containing gene IDs to analyze only selected genes.

Genome-wide motif scanning can be computationally intensive, especially when using large motif collections. For routine analyses, using a gene list to restrict the analysis to genes of interest is recommended whenever possible.

### Requirements

* Python 3
* biopython
* pandas
* pysam
* MEME Suite with `fimo` available in `PATH`

### Tested environment

The current version has been tested with the following environment:

* Python 3.12.8
* biopython 1.85
* pandas 2.3.3
* pysam 0.23.3
* MEME Suite 5.5.6

Install the required Python dependencies using:

```bash
pip install -r requirements.txt
```

MEME Suite/FIMO can be installed separately, for example:

```bash
mamba install -c bioconda meme
```

## Command-line options

The following command-line arguments are available:

| Argument | Required | Default | Description |
|---|---|---:|---|
| `--meme` | Yes | — | Input motif file in MEME format. |
| `--fasta` | Yes | — | Reference genome FASTA file. |
| `--gff` | Yes | — | Input GFF3 annotation file. |
| `--output-dir` | Yes | — | Output directory. The directory is created automatically if it does not already exist. |
| `--gene-list` | Yes | — | Gene selection. Use `all` to analyze all genes in the GFF3 file, or provide a text file containing one gene ID per line. |
| `--flank` | No | `2000` | Length of the flanking region added to both sides of each gene, in base pairs. Must be `0` or greater. |
| `--threads` | No | `8` | Maximum number of FIMO processes that can run in parallel. |
| `--batches` | No | Same as `--threads` | Number of multi-FASTA batches used for FIMO analysis. |
| `--thresh` | No | `1e-4` | FIMO p-value threshold used for motif prediction. |
| `--output` | No | Automatically generated | Optional filename or path for the output TSV. By default, the FASTA filename is included automatically, for example `motif_predictions.FZX.tsv`. |
| `--fimo-bin` | No | `fimo` | FIMO executable name or path. |
| `--debug` | No | Disabled | Enable debug mode. Only the first `--debug-n-genes` genes are processed. |
| `--debug-n-genes` | No | `100` | Number of genes processed when `--debug` is enabled. Must be at least `1`. |

## Usage

Analyze selected genes using a gene list:

```bash
python motif_predict.py \
  --meme motifs.meme \
  --fasta genome.fa \
  --gff annotation.gff3 \
  --gene-list genes.txt \
  --output-dir motif_output
```

The gene-list file should contain one gene ID per line, for example:

```text
Gene001
Gene002
Gene003
```

To analyze all genes in the GFF3 file, use:

```bash
python motif_predict.py \
  --meme motifs.meme \
  --fasta genome.fa \
  --gff annotation.gff3 \
  --gene-list all \
  --output-dir motif_output
```

Show all available options with:

```bash
python motif_predict.py -h
```

## Output

If an output file with the same name already exists, a numeric suffix is added automatically.

For output files larger than 100 MB, the program automatically generates bgzip-compressed and Tabix-indexed files.
