# demultiplex_irepertoire

Stream paired-end FASTQ files into per-sample FASTQ files using a fixed-position inline barcode. The script reads R1 and R2 together, assigns each read pair to the sample whose barcode matches, and writes the results without loading the full files into memory.

## Requirements

- Python 3 (standard library only; no extra packages)
- Gzipped paired FASTQ input (`.fastq.gz`)

## Sample sheet

A headered, tab-separated file with one row per sample. By default the script reads the `sample_id` and `barcode` columns. Barcodes must be A/C/G/T only, the same length for every sample, and unique.

```
sample_id	barcode
S1	ACGTAC
S2	TGCATG
```

Use `--sample-column` and `--barcode-column` if your headers differ.

## Usage

```bash
python3 demultiplex_irepertoire.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --sample-sheet samples.tsv \
  --output-dir demux_out/ \
  --run-label myrun \
  --barcode-start 3
```

### Required parameters

| Parameter | Meaning |
|---|---|
| `--r1`, `--r2` | Input paired FASTQ files (gzipped) |
| `--sample-sheet` | Tab-separated sample table (see above) |
| `--output-dir` | Output directory; unrelated existing files are preserved |
| `--run-label` | Label used in output names: `<sample_id>_<run_label>_R1/R2.fastq.gz` |
| `--barcode-start` | Zero-based offset of the barcode in the selected read (`3` means bases 4 onward) |

### Optional parameters (defaults)

| Parameter | Default | Meaning |
|---|---|---|
| `--sample-column` | `sample_id` | Sample table column for output names |
| `--barcode-column` | `barcode` | Sample table column with barcode sequences |
| `--barcode-read` | `r1` | Read that carries the barcode (`r1` or `r2`) |
| `--reverse-complement` | off | Reverse-complement table barcodes before matching |
| `--max-mismatches` | `0` | Maximum substitutions allowed; tied nearest matches are ambiguous. No indels |
| `--trim-prefix` / `--keep-prefix` | trim | Trim assigned R1 reads through the end of the barcode (requires `--barcode-read r1`) |
| `--r1-trim-tail` | `1` | Trailing R1 bases to remove |
| `--r2-trim-tail` | `1` | Trailing R2 bases to remove |
| `--barcode-in-read-id` / `--no-barcode-in-read-id` | off | Append `:<observed barcode>` to read IDs. Off by default |
| `--max-pairs` | `0` | Process only this many pairs (0 = all). Useful for a quick preview |
| `--compression-level` | `1` | Gzip level 1–9 for output files |
| `--progress-every` | `1000000` | Progress interval in pairs (0 disables) |

Run `python3 demultiplex_irepertoire.py --help` for the full list.

## Outputs

In `--output-dir`:

- `<sample_id>_<run_label>_R1.fastq.gz` and `_R2.fastq.gz` for each sample
- `unassigned_<run_label>_R1/R2.fastq.gz` for pairs with no barcode match
- `ambiguous_<run_label>_R1/R2.fastq.gz` for pairs tied between samples
- `<run_label>_summary.tsv` with read-pair counts per sample
- `<run_label>_observed_barcodes.tsv` with observed barcode counts
- `<run_label>_run.json` with the parameters and totals of the run

The script prints its full parameter list at startup.

## Example preview

Check a run on the first 20,000 pairs before processing everything:

```bash
python3 demultiplex_irepertoire.py ... --max-pairs 20000
```
