#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def open_fastq(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def canonical_kmers(sequence, k=15):
    sequence = str(sequence).upper()
    reverse = str(Seq(sequence).reverse_complement())
    return {min(sequence[i:i + k], reverse[len(sequence) - i - k:len(sequence) - i]) for i in range(len(sequence) - k + 1)}


def find_length_peaks(records, bin_width=200, minimum_reads=10, minimum_plasmid_length=2000):
    lengths = np.array([len(record.seq) for record in records])
    edges = np.arange(0, lengths.max() + bin_width * 2, bin_width)
    counts, _ = np.histogram(lengths, bins=edges)
    threshold = max(minimum_reads, int(len(records) * 0.015))
    candidates = [index for index, count in enumerate(counts) if edges[index] >= minimum_plasmid_length and count >= threshold and count >= max(counts[max(0, index - 1):index + 2])]
    peaks = []
    for index in candidates:
        center = int((edges[index] + edges[index + 1]) / 2)
        members = [record for record in records if abs(len(record.seq) - center) <= max(300, center * 0.06)]
        if members:
            peaks.append({"center": int(np.median([len(record.seq) for record in members])), "records": members, "count": len(members)})
    return sorted(peaks, key=lambda peak: peak["count"], reverse=True)


def best_record(records):
    median_length = float(np.median([len(record.seq) for record in records]))
    central = [record for record in records if abs(len(record.seq) - median_length) <= max(25, median_length * 0.01)]
    return max(central or records, key=lambda record: np.mean(record.letter_annotations["phred_quality"]))


def collapse_multimers(peaks, sketch_k=15, containment_threshold=0.80):
    distinct = []
    for peak in peaks:
        representative = best_record(peak["records"])
        sketch = canonical_kmers(representative.seq, sketch_k)
        relation = None
        for accepted in distinct:
            length_ratio = peak["center"] / accepted["center"]
            containment = len(sketch & accepted["sketch"]) / max(1, min(len(sketch), len(accepted["sketch"])))
            copy_number = max(1, round(length_ratio))
            if containment >= containment_threshold and abs(length_ratio - copy_number) <= 0.15 * copy_number:
                relation = {"candidate": accepted, "copies": copy_number, "containment": containment}
                break
        if relation:
            relation["candidate"].setdefault("multimer_peaks", []).append({**peak, **relation})
        else:
            distinct.append({**peak, "representative": representative, "sketch": sketch, "multimer_peaks": []})
    return distinct


def debruijn_consensus(records, expected_length, k=41):
    counts = Counter()
    for record in records:
        sequence = str(record.seq).upper()
        for strand in (sequence, str(Seq(sequence).reverse_complement())):
            counts.update(strand[i:i + k] for i in range(len(strand) - k + 1) if "N" not in strand[i:i + k])
    representative = str(best_record(records).seq).upper()
    start_candidates = [representative[i:i + k] for i in range(len(representative) - k + 1)]
    start = max(start_candidates, key=lambda value: counts[value])
    current = start
    assembled = start
    visited = {start}
    closed = False
    edge_count = 0
    for _ in range(int(expected_length * 1.25)):
        candidates = [(counts[current[1:] + base], current[1:] + base, base) for base in "ACGT"]
        candidates.sort(reverse=True)
        selected = None
        for support, next_kmer, base in candidates:
            if not support:
                continue
            if next_kmer == start and edge_count >= expected_length * 0.75:
                selected = (support, next_kmer, base)
                closed = True
                break
            if next_kmer not in visited:
                selected = (support, next_kmer, base)
                break
        if selected is None:
            break
        _, current, base = selected
        edge_count += 1
        if closed:
            break
        visited.add(current)
        assembled += base
    if closed:
        sequence = assembled[:edge_count]
    else:
        sequence = representative
    return Seq(sequence), closed, counts


def cluster_support(sequence, records, k=15):
    reference = canonical_kmers(sequence, k)
    values = []
    for record in records:
        read = canonical_kmers(record.seq, k)
        values.append(len(read & reference) / max(1, min(len(read), len(reference))))
    return float(np.median(values))


def analyze(fastq_paths, output_folder, top_n=2):
    output_folder.mkdir(parents=True, exist_ok=True)
    sample_records = {}
    all_records = []
    for path in fastq_paths:
        with open_fastq(path) as handle:
            records = list(SeqIO.parse(handle, "fastq"))
        sample_records[Path(path).stem.replace(".fastq", "")] = records
        all_records.extend(records)
    peaks = find_length_peaks(all_records)
    candidates = collapse_multimers(peaks)
    prepared = []
    for candidate in candidates:
        sequence, circular, _ = debruijn_consensus(candidate["records"], candidate["center"])
        support = cluster_support(sequence, candidate["records"])
        prepared.append((candidate, sequence, circular, support))
    prepared.sort(key=lambda item: (item[2], item[3], item[0]["center"], item[0]["count"]), reverse=True)
    prepared = prepared[:top_n]

    rows = []
    for rank, (candidate, sequence, circular, support) in enumerate(prepared, start=1):
        record = SeqRecord(
            sequence,
            id=f"de_novo_plasmid_{rank}",
            name=f"plasmid_candidate_{rank}",
            description=f"No-reference plasmid candidate {rank}; inferred from Nanopore read clustering",
        )
        record.annotations.update({
            "molecule_type": "DNA",
            "topology": "circular" if circular else "linear",
            "data_file_division": "SYN",
            "date": "01-SEP-2026",
        })
        notes = [
            f"rank: {rank}",
            f"cluster reads: {candidate['count']}",
            f"median cluster length: {candidate['center']} bp",
            f"median sequence k-mer support: {support:.4f}",
            f"circular de Bruijn path: {'yes' if circular else 'no; highest-quality representative retained'}",
        ]
        if candidate["multimer_peaks"]:
            notes.extend(f"supporting {item['copies']}-mer peak: {item['center']} bp ({item['count']} reads; containment {item['containment']:.3f})" for item in candidate["multimer_peaks"])
        record.features.append(SeqFeature(FeatureLocation(0, len(record)), type="source", qualifiers={"organism": ["synthetic construct"], "note": notes}))
        record.features.append(SeqFeature(FeatureLocation(0, len(record)), type="misc_feature", qualifiers={"label": [f"de novo plasmid candidate {rank}"], "note": ["Sequence is unannotated and should be validated with an independent assembler/polisher before synthesis or publication."]}))
        gbk_path = output_folder / f"plasmid_candidate_{rank}.gbk"
        fasta_path = output_folder / f"plasmid_candidate_{rank}.fasta"
        SeqIO.write(record, gbk_path, "genbank")
        SeqIO.write(record, fasta_path, "fasta")
        sample_counts = {name: sum(abs(len(read.seq) - candidate["center"]) <= max(300, candidate["center"] * 0.06) for read in records) for name, records in sample_records.items()}
        rows.append({"rank": rank, "length": len(record.seq), "topology": record.annotations["topology"], "cluster_reads": candidate["count"], "median_kmer_support": f"{support:.4f}", "multimer_peaks": "; ".join(f"{item['copies']}-mer:{item['center']}bp/{item['count']}reads" for item in candidate["multimer_peaks"]), **sample_counts})

    with open(output_folder / "no_reference_summary.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(len(sample_records), 1, figsize=(12, 4 * len(sample_records)), squeeze=False)
    for axis, (name, records) in zip(axes[:, 0], sample_records.items()):
        lengths = [len(record.seq) for record in records]
        axis.hist(lengths, bins=np.arange(0, max(lengths) + 201, 200), color="#2C7FB8")
        for row in rows:
            axis.axvline(row["length"], linestyle="--", label=f"candidate {row['rank']}: {row['length']} bp")
            if 2 * row["length"] <= max(lengths) * 1.05:
                axis.axvline(2 * row["length"], linestyle=":", alpha=0.7)
        axis.set(title=name, xlabel="Read length (bp)", ylabel="Read count")
        axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_folder / "no_reference_read_lengths.png", dpi=200)
    plt.close(fig)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="No-reference plasmid discovery from high-accuracy full-length Nanopore reads")
    parser.add_argument("fastq", nargs="+")
    parser.add_argument("-o", "--output-folder", required=True, type=Path)
    parser.add_argument("--top", type=int, default=2)
    arguments = parser.parse_args()
    results = analyze(arguments.fastq, arguments.output_folder, arguments.top)
    for result in results:
        print(result)
