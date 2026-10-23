#!/usr/bin/env python3
"""Stream paired FASTQs into sample files using a fixed-position inline barcode."""

import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    required = parser.add_argument_group("Required Parameters")
    required.add_argument("--r1", required=True, help="Input R1 FASTQ.gz")
    required.add_argument("--r2", required=True, help="Input R2 FASTQ.gz")
    required.add_argument("--sample-sheet", required=True, help="Headered tab-separated sample table")
    required.add_argument("--output-dir", required=True, help="Output directory; existing unrelated files are preserved")
    required.add_argument("--run-label", required=True, help="Run label used in output names: <sample_id>_<run_label>_R1/R2.fastq.gz")
    required.add_argument("--barcode-start", required=True, type=int, help="Zero-based barcode offset in the selected read; 3 means bases 4 onward")
    optional = parser.add_argument_group("Optional Parameters")
    optional.add_argument("--sample-column", default="sample_id", help="Sample table column used for output filenames")
    optional.add_argument("--barcode-column", default="barcode", help="Sample table column containing A/C/G/T barcode sequences")
    optional.add_argument("--barcode-read", choices=["r1", "r2"], default="r1", help="Read containing the inline barcode")
    optional.add_argument("--reverse-complement", action="store_true", help="Reverse-complement table barcodes before matching")
    optional.add_argument("--max-mismatches", type=int, default=0, help="Maximum substitutions; tied nearest matches remain ambiguous; no indels")
    optional.add_argument("--trim-prefix", dest="trim_prefix", action="store_true", default=True, help="Trim assigned R1 reads through the last base of the matched barcode")
    optional.add_argument("--keep-prefix", dest="trim_prefix", action="store_false", help="Disable R1 barcode-prefix trimming")
    optional.add_argument("--r1-trim-tail", type=int, default=1, help="Number of trailing R1 bases to remove")
    optional.add_argument("--r2-trim-tail", type=int, default=1, help="Number of trailing R2 bases to remove")
    optional.add_argument("--barcode-in-read-id", dest="barcode_in_read_id", action="store_true", default=False, help="Add :<observed barcode> to both mate IDs (off by default)")
    optional.add_argument("--no-barcode-in-read-id", dest="barcode_in_read_id", action="store_false", default=False, help="Keep read IDs unchanged (default; kept for older command lines)")
    optional.add_argument("--max-pairs", type=int, default=0, help="Process only this many pairs for preview; 0 processes all pairs")
    optional.add_argument("--compression-level", type=int, choices=range(1, 10), default=1, help="Gzip compression level")
    optional.add_argument("--progress-every", type=int, default=1000000, help="Report progress every N pairs; 0 disables progress")
    args = parser.parse_args()
    if min(args.barcode_start, args.max_mismatches, args.max_pairs, args.progress_every,
           args.r1_trim_tail, args.r2_trim_tail) < 0:
        parser.error("Offsets, trimming lengths, mismatches, pair limit, and progress interval must be nonnegative")
    if args.trim_prefix and args.barcode_read != "r1":
        parser.error("R1 barcode-prefix trimming requires --barcode-read r1; use --keep-prefix for R2 barcodes")
    return args


# Parse arguments before importing other packages or inspecting input data.
ARGS = parse_arguments()
print("Input parameters and configuration:", flush=True)
for key, value in sorted(vars(ARGS).items()):
    print(f"  {key}: {value}", flush=True)

import csv
import gzip
import json
import re
import sys
from collections import Counter
from contextlib import ExitStack
from pathlib import Path


def load_samples(args):
    samples = {}
    seen_barcodes = set()
    with open(args.sample_sheet, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {args.sample_column, args.barcode_column}.issubset(reader.fieldnames):
            raise ValueError("Sample sheet must contain the selected sample and barcode column headers")
        for row in reader:
            sample = row[args.sample_column].strip()
            barcode = row[args.barcode_column].strip().upper()
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
                raise ValueError(f"Unsafe or empty sample ID: {sample!r}")
            if sample in {"unassigned", "ambiguous"} or sample in samples:
                raise ValueError(f"Reserved or duplicate sample ID: {sample}")
            if not re.fullmatch(r"[ACGT]+", barcode) or barcode in seen_barcodes:
                raise ValueError(f"Invalid or duplicate barcode: {barcode}")
            seen_barcodes.add(barcode)
            if args.reverse_complement:
                barcode = barcode.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            samples[sample] = barcode
    lengths = {len(barcode) for barcode in samples.values()}
    if len(lengths) != 1:
        raise ValueError("Provide at least one sample; all barcodes must have the same length")
    length = lengths.pop()
    if args.max_mismatches >= length:
        raise ValueError("Maximum mismatches must be smaller than the barcode length")
    return samples, length


def read_record(handle):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().rstrip("\r\n")
    plus = handle.readline().rstrip("\r\n")
    quality_line = handle.readline()
    quality = quality_line.rstrip("\r\n")
    if not header.startswith("@") or not plus.startswith("+") or not quality_line:
        raise ValueError(f"Malformed or truncated FASTQ record: {header.strip()}")
    if not sequence or len(sequence) != len(quality):
        raise ValueError(f"Sequence/quality length mismatch or empty sequence: {header.strip()}")
    return [header.rstrip("\r\n"), sequence, plus, quality]


def pair_id(record):
    identifier = record[0].split()[0]
    if identifier.endswith(("/1", "/2")):
        identifier = identifier[:-2]
    return identifier


def classify(observed, samples, length, max_mismatches):
    if len(observed) != length or any(base not in "ACGT" for base in observed):
        return "unassigned", None
    distances = [(sum(a != b for a, b in zip(observed, barcode)), sample)
                 for sample, barcode in samples.items()]
    best_distance = min(distance for distance, sample in distances)
    if best_distance > max_mismatches:
        return "unassigned", None
    winners = [sample for distance, sample in distances if distance == best_distance]
    if len(winners) != 1:
        return "ambiguous", None
    return winners[0], best_distance


def tag_barcode(record, barcode):
    header = record[0].split(maxsplit=1)
    identifier = header[0]
    mate_suffix = ""
    if identifier.endswith(("/1", "/2")):
        mate_suffix = identifier[-2:]
        identifier = identifier[:-2]
    record[0] = f"{identifier}:{barcode}{mate_suffix}"
    if len(header) == 2:
        record[0] += " " + header[1]
    if record[2] != "+":
        record[2] = "+" + record[0][1:]


def trim_record(record, front, tail, pair_number, mate):
    if front + tail >= len(record[1]):
        raise ValueError(f"Trimming would empty {mate} at pair {pair_number}")
    end = len(record[1]) - tail
    record[1] = record[1][front:end]
    record[3] = record[3][front:end]


def main(args):
    samples, barcode_length = load_samples(args)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_label):
        raise ValueError("Run label must contain only letters, digits, underscores, dots, or hyphens and start with a letter or digit")
    if Path(args.r1).resolve() == Path(args.r2).resolve():
        raise ValueError("R1 and R2 must be different files")
    for path in (args.r1, args.r2):
        if not Path(path).is_file():
            raise ValueError(f"Input FASTQ not found: {path}")
    barcode_end = args.barcode_start + barcode_length
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    output_paths = {
        sample: [output / f"{sample}_{args.run_label}_{mate}.fastq.gz" for mate in ("R1", "R2")]
        for sample in list(samples) + ["unassigned", "ambiguous"]
    }
    metadata_path = output / f"{args.run_label}_run.json"
    summary_path = output / f"{args.run_label}_summary.tsv"
    observed_path = output / f"{args.run_label}_observed_barcodes.tsv"
    planned_paths = [metadata_path, summary_path, observed_path]
    for paths in output_paths.values():
        planned_paths.extend(paths)
    for path in planned_paths:
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Output already exists; refusing to overwrite: {path}")
    metadata = {
        "parameters": vars(args),
        "effective_barcodes": samples,
        "barcode_length": barcode_length,
        "r1_prefix_trim_method": "matched_barcode_end" if args.trim_prefix else "disabled",
        "barcode_end_offset": barcode_end,
        "status": "running",
        "processed_pairs": 0,
    }
    with metadata_path.open("x") as handle:
        handle.write(json.dumps(metadata, indent=2) + "\n")
    counts = Counter()
    exact = Counter()
    rescued = Counter()
    observed_counts = Counter()
    cache = {}
    total = 0
    prefix_trimmed = 0
    try:
        with ExitStack() as stack:
            r1 = stack.enter_context(gzip.open(args.r1, "rt", encoding="ascii"))
            r2 = stack.enter_context(gzip.open(args.r2, "rt", encoding="ascii"))
            outputs = {}
            for sample in list(samples) + ["unassigned", "ambiguous"]:
                outputs[sample] = [
                    stack.enter_context(gzip.open(
                        path, "xt",
                        encoding="ascii", compresslevel=args.compression_level,
                    ))
                    for path in output_paths[sample]
                ]
            while not args.max_pairs or total < args.max_pairs:
                records = [read_record(r1), read_record(r2)]
                if records == [None, None]:
                    break
                if any(record is None for record in records):
                    raise ValueError(f"R1/R2 have different record counts after pair {total}")
                if pair_id(records[0]) != pair_id(records[1]):
                    raise ValueError(f"R1/R2 identifiers differ at pair {total + 1}")
                selected = records[0 if args.barcode_read == "r1" else 1]
                observed = selected[1][args.barcode_start:barcode_end].upper()
                observed_counts[observed] += 1
                if observed not in cache:
                    cache[observed] = classify(observed, samples, barcode_length, args.max_mismatches)
                sample, distance = cache[observed]
                counts[sample] += 1
                if distance == 0:
                    exact[sample] += 1
                elif distance is not None:
                    rescued[sample] += 1
                if args.barcode_in_read_id and len(observed) == barcode_length:
                    for record in records:
                        tag_barcode(record, observed)
                r1_front = 0
                if args.trim_prefix and sample in samples:
                    r1_front = barcode_end
                    prefix_trimmed += 1
                trim_record(records[0], r1_front, args.r1_trim_tail, total + 1, "R1")
                trim_record(records[1], 0, args.r2_trim_tail, total + 1, "R2")
                for record, handle in zip(records, outputs[sample]):
                    handle.write("\n".join(record) + "\n")
                total += 1
                if args.progress_every and total % args.progress_every == 0:
                    print(f"Processed {total:,} pairs", flush=True)
        with open(summary_path, "x", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id", "effective_barcode", "read_pairs", "fraction", "exact_pairs", "mismatch_rescued_pairs"])
            for sample in list(samples) + ["unassigned", "ambiguous"]:
                writer.writerow([sample, samples.get(sample, ""), counts[sample],
                                 counts[sample] / total if total else 0, exact[sample], rescued[sample]])
        with open(observed_path, "x", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["observed_barcode", "read_pairs", "assignment"])
            for observed, count in observed_counts.most_common():
                writer.writerow([observed, count, cache[observed][0]])
        metadata["status"] = "completed_preview" if args.max_pairs and total == args.max_pairs else "completed"
    except BaseException as error:
        metadata["status"] = "failed"
        metadata["error"] = str(error)
        raise
    finally:
        metadata["processed_pairs"] = total
        metadata["read_pairs_by_assignment"] = dict(counts)
        metadata["r1_prefix_trimmed_pairs"] = prefix_trimmed
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Finished: {total:,} pairs; {sum(counts[sample] for sample in samples):,} assigned", flush=True)
    print(f"Summary: {summary_path}", flush=True)


if __name__ == "__main__":
    try:
        main(ARGS)
    except (OSError, ValueError, csv.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
