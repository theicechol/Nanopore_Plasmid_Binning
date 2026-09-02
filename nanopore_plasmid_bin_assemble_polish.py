#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import subprocess
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plasmid_qc_report import generate_qc_report


MARKER_WORDS = ("origin", "ori", "resistance", "ampr", "kanr", "cmr", "smr", "tetr", "beta-lactamase", "chloramphenicol", "kanamycin", "streptomycin", "tetracycline")


def open_fastq(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def canonical_kmer(sequence):
    return min(sequence, str(Seq(sequence).reverse_complement()))


def sequence_kmers(sequence, k):
    sequence = str(sequence).upper()
    return {canonical_kmer(sequence[i:i + k]) for i in range(len(sequence) - k + 1) if set(sequence[i:i + k]) <= {"A", "C", "G", "T"}}


def feature_label(feature):
    values = []
    for key in ("label", "gene", "product", "note"):
        values.extend(feature.qualifiers.get(key, []))
    return "; ".join(values)


def load_references(paths, k):
    references, marker_kmers = {}, {}
    for path_string in paths:
        path = Path(path_string)
        file_format = "genbank" if path.suffix.lower() in {".gb", ".gbk", ".genbank"} else "fasta"
        records = list(SeqIO.parse(path, file_format))
        if len(records) != 1:
            raise ValueError(f"Reference {path} must contain exactly one sequence")
        record = records[0]
        name = re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem or record.id)
        if name in references:
            raise ValueError(f"Duplicate reference name: {name}")
        references[name], marker_kmers[name] = record, {}
        if file_format == "genbank":
            for feature in record.features:
                label = feature_label(feature)
                if label and any(word in label.lower() for word in MARKER_WORDS):
                    kmers = sequence_kmers(feature.extract(record.seq), k)
                    if kmers:
                        marker_kmers[name][label] = kmers

    all_kmers = {name: sequence_kmers(record.seq, k) for name, record in references.items()}
    unique_kmers = {}
    for name, kmers in all_kmers.items():
        others = set().union(*(value for other, value in all_kmers.items() if other != name))
        unique_kmers[name] = kmers - others
    for name, markers in marker_kmers.items():
        for label in list(markers):
            markers[label] &= unique_kmers[name]
            if not markers[label]:
                del markers[label]
    return references, unique_kmers, marker_kmers


def classify_read(record, unique_kmers, marker_kmers, k, min_hits, min_ratio):
    read_kmers = sequence_kmers(record.seq, k)
    hits = {name: len(read_kmers & kmers) for name, kmers in unique_kmers.items()}
    ranked = sorted(hits, key=hits.get, reverse=True)
    best = ranked[0]
    second_hits = hits[ranked[1]] if len(ranked) > 1 else 0
    ratio = hits[best] / max(1, second_hits)
    assignment = best if hits[best] >= min_hits and ratio >= min_ratio else "ambiguous"
    markers = [] if assignment == "ambiguous" else [label for label, kmers in marker_kmers[assignment].items() if len(read_kmers & kmers) >= min_hits]
    return assignment, hits, ratio, markers


def reference_bin_reads(input_fastq, output_folder, reference_paths, k, min_hits, min_ratio):
    references, unique_kmers, marker_kmers = load_references(reference_paths, k)
    handles = {name: open(output_folder / f"{name}.fastq", "w") for name in [*references, "ambiguous"]}
    counts = Counter()
    try:
        with open(output_folder / "read_assignments.tsv", "w", newline="") as report, open_fastq(input_fastq) as source:
            writer = csv.writer(report, delimiter="\t")
            writer.writerow(["read_id", "length", "assignment", "best_to_second_ratio", "marker_evidence", *references])
            for record in SeqIO.parse(source, "fastq"):
                assignment, hits, ratio, markers = classify_read(record, unique_kmers, marker_kmers, k, min_hits, min_ratio)
                SeqIO.write(record, handles[assignment], "fastq")
                counts[assignment] += 1
                writer.writerow([record.id, len(record.seq), assignment, f"{ratio:.3f}", "; ".join(markers), *(hits[name] for name in references)])
    finally:
        for handle in handles.values():
            handle.close()
    print("Reference-specific k-mers:", ", ".join(f"{name}={len(kmers)}" for name, kmers in unique_kmers.items()))
    print("Read assignments:", ", ".join(f"{name}={counts[name]}" for name in [*references, "ambiguous"]))
    return [output_folder / f"{name}.fastq" for name in references if counts[name] > 0]


def length_bin_reads(input_fastq, output_folder):
    with open_fastq(input_fastq) as source:
        records = list(SeqIO.parse(source, "fastq"))
    if not records:
        raise ValueError("Input FASTQ contains no reads")
    lengths = np.array([len(record.seq) for record in records])
    edges = np.arange(lengths.min(), lengths.max() + 201, 200)
    counts, _ = np.histogram(lengths, bins=edges)
    cutoff = counts.mean() + 3 * counts.std()
    read_files = []
    for index in np.where(counts > cutoff)[0]:
        center = int(edges[index] + 100)
        reads = [record for record in records if center - 200 <= len(record.seq) <= center + 200]
        path = output_folder / f"length_{center}.fastq"
        SeqIO.write(reads, path, "fastq")
        read_files.append(path)
    sns.histplot(lengths, binwidth=200)
    plt.axhline(y=cutoff, color="gray", linestyle="--", label="mean + 3 SD")
    plt.xlabel("Read length")
    plt.savefig(output_folder / "read_length_hist.png", bbox_inches="tight", dpi=300)
    plt.close()
    return read_files


def run_command(command):
    print("Running:", " ".join(map(str, command)))
    subprocess.run([str(value) for value in command], check=True)


def assemble_and_polish(read_files, output_folder, canu_binary_path, num_threads):
    for read_file in read_files:
        name = read_file.stem
        lengths = [len(record.seq) for record in SeqIO.parse(read_file, "fastq")]
        if not lengths:
            continue
        genome_size = max(1, round(np.median(lengths) / 1000))
        canu_folder = output_folder / f"canu_output_{name}"
        run_command([canu_binary_path, "-p", name, "-d", canu_folder, f"genomeSize={genome_size}k", f"maxThreads={num_threads}", "-nanopore", read_file])
        contigs = list(canu_folder.glob("*.contigs.fasta"))
        if not contigs:
            raise RuntimeError(f"Canu produced no contigs for {read_file}")
        for count, record in enumerate(SeqIO.parse(contigs[0], "fasta"), start=1):
            sequence = record.seq
            if "suggestCircular=yes" in record.description and "trim=" in record.description:
                bounds = record.description.split("trim=", 1)[1].split()[0].split("-")
                sequence = sequence[int(bounds[0]):int(bounds[1])]
            trimmed = canu_folder / f"{name}_{count}.trimmed_contigs.fasta"
            SeqIO.write(SeqRecord(sequence, id=f"{name}_{count}", description=""), trimmed, "fasta")
            polished = output_folder / f"polished_output_{name}_{count}"
            run_command(["medaka_consensus", "-i", read_file, "-d", trimmed, "-o", polished, "-t", num_threads, "--bacteria"])


def main(args):
    output_folder = Path(args.output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    if args.references:
        read_files = reference_bin_reads(Path(args.input_fastq), output_folder, args.references, args.kmer_size, args.min_unique_kmers, args.min_assignment_ratio)
    else:
        read_files = length_bin_reads(Path(args.input_fastq), output_folder)
    for path in read_files:
        print(f"{path}: {sum(1 for _ in SeqIO.parse(path, 'fastq'))} reads")
    if not args.no_qc_report:
        report = generate_qc_report(
            input_fastq=args.input_fastq,
            output_folder=output_folder,
            reference_paths=args.references,
            assignment_path=output_folder / "read_assignments.tsv",
            host_reference=args.host_reference,
            host_name=args.host_name,
            host_key=args.host,
            skip_host=args.skip_host_screen,
            multimer_tolerance=args.multimer_tolerance,
        )
        print(f"QC report: {report}")
    if not args.binning_only:
        if not args.canu_binary_path:
            raise ValueError("--canu_binary_path is required unless --binning-only is used")
        assemble_and_polish(read_files, output_folder, args.canu_binary_path, args.num_threads)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bin, assemble, and polish Nanopore plasmid reads")
    parser.add_argument("--output_folder", "-o", required=True)
    parser.add_argument("--input_fastq", "-i", required=True)
    parser.add_argument("--num_threads", "-t", default=1, type=int)
    parser.add_argument("--canu_binary_path", "-c")
    parser.add_argument("--references", "-r", nargs="+", help="Candidate plasmid GenBank or FASTA files")
    parser.add_argument("--kmer-size", type=int, default=15)
    parser.add_argument("--min-unique-kmers", type=int, default=10)
    parser.add_argument("--min-assignment-ratio", type=float, default=2.0)
    parser.add_argument("--binning-only", action="store_true", help="Write read bins without assembly or polishing")
    parser.add_argument("--no-qc-report", action="store_true", help="Do not generate the graphical HTML QC report")
    parser.add_argument("--host-reference", help="Host genome FASTA (default: download E. coli K-12 MG1655 from NCBI)")
    parser.add_argument("--host", choices=("ecoli", "vibrio-natriegens", "bacillus-subtilis"), default="ecoli", help="Built-in cloning-host genome used for contamination screening")
    parser.add_argument("--host-name", default="E. coli", help="Host label shown in the report")
    parser.add_argument("--skip-host-screen", action="store_true", help="Generate QC without host-genome mapping")
    parser.add_argument("--multimer-tolerance", type=float, default=0.15, help="Relative length tolerance for monomer/dimer calls")
    main(parser.parse_args())
