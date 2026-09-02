import csv
import html
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO


HOST_CATALOG = {
    "ecoli": {
        "label": "E. coli K-12 MG1655",
        "accessions": ("NC_000913.3",),
        "filename": "E_coli_K12_MG1655_NC_000913.3.fasta",
    },
    "vibrio-natriegens": {
        "label": "Vibrio natriegens ATCC 14048",
        "accessions": ("CP009977.1", "CP009978.1"),
        "filename": "Vibrio_natriegens_ATCC_14048_CP009977.1_CP009978.1.fasta",
    },
    "bacillus-subtilis": {
        "label": "Bacillus subtilis subsp. subtilis 168",
        "accessions": ("NC_000964.3",),
        "filename": "Bacillus_subtilis_168_NC_000964.3.fasta",
    },
}
COLORS = {"host": "#E69F00", "ambiguous": "#999999", "unassigned": "#56B4E9"}


def _load_assignments(path):
    with open(path, newline="") as handle:
        return {row["read_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _load_reference_lengths(reference_paths):
    lengths = {}
    for value in reference_paths or []:
        path = Path(value)
        fmt = "genbank" if path.suffix.lower() in {".gb", ".gbk", ".genbank"} else "fasta"
        record = next(SeqIO.parse(path, fmt))
        lengths[path.stem] = len(record.seq)
    return lengths


def download_catalog_host(host_key, destination):
    host = HOST_CATALOG[host_key]
    accessions = ",".join(host["accessions"])
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accessions}&rettype=fasta&retmode=text"
    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {host['label']} ({accessions}) from NCBI")
    urllib.request.urlretrieve(url, destination)
    if destination.stat().st_size < 1000:
        destination.unlink(missing_ok=True)
        raise RuntimeError(f"NCBI returned an incomplete host reference for {host_key}")
    return destination


def resolve_host_reference(output_folder, host_reference, host_name, host_key="ecoli"):
    if host_reference:
        return Path(host_reference), host_name
    host = HOST_CATALOG[host_key]
    project_reference = Path(__file__).resolve().parent / "references" / "hosts" / host["filename"]
    cache_reference = output_folder / "references" / host["filename"]
    path = project_reference if project_reference.exists() else cache_reference
    if not path.exists():
        download_catalog_host(host_key, path)
    accessions = ", ".join(host["accessions"])
    return path, f"{host['label']} ({accessions})"


def screen_host_reads(input_fastq, host_reference, min_fraction=0.50, min_mapq=20):
    try:
        import mappy
    except ImportError:
        return set(), "mappy is unavailable; host screening was not performed"
    aligner = mappy.Aligner(str(host_reference), preset="map-ont")
    if not aligner:
        return set(), f"could not index host reference {host_reference}"
    host_reads = set()
    opener = __import__("gzip").open if str(input_fastq).endswith(".gz") else open
    with opener(input_fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            for hit in aligner.map(str(record.seq)):
                aligned_fraction = (hit.q_en - hit.q_st) / max(1, len(record.seq))
                if not hit.is_secondary and hit.mapq >= min_mapq and aligned_fraction >= min_fraction:
                    host_reads.add(record.id)
                    break
    return host_reads, "complete"


def classify_multimer(length, reference_length, tolerance):
    copies = length / reference_length
    nearest = max(1, round(copies))
    if nearest <= 4 and abs(copies - nearest) <= tolerance * nearest:
        return f"{nearest}-mer"
    return "other"


def _save_length_plot(records, assignments, host_reads, reference_lengths, output):
    groups = defaultdict(list)
    for record in records:
        assignment = assignments.get(record.id, {}).get("assignment", "unassigned")
        group = "host" if record.id in host_reads else assignment
        groups[group].append(len(record.seq))
    labels = list(reference_lengths) + [name for name in ("host", "ambiguous", "unassigned") if groups[name]]
    values = [groups[label] for label in labels if groups[label]]
    labels = [label for label in labels if groups[label]]
    maximum = max(len(record.seq) for record in records)
    bins = np.arange(0, maximum + 201, 200)
    palette = [COLORS.get(label, plt.cm.tab10(index % 10)) for index, label in enumerate(labels)]
    fig, axis = plt.subplots(figsize=(12, 6))
    axis.hist(values, bins=bins, stacked=True, label=labels, color=palette)
    for ref_index, (name, length) in enumerate(reference_lengths.items()):
        axis.axvspan(length * 0.90, length * 1.10, color="#BBBBBB", alpha=0.15)
        axis.text(length, axis.get_ylim()[1] * (0.96 - 0.08 * (ref_index % 3)), f"{name}\nmonomer", ha="center", va="top", fontsize=8)
        if 2 * length <= maximum:
            axis.axvspan(length * 1.90, length * 2.10, color="#BBBBBB", alpha=0.10)
            axis.text(2 * length, axis.get_ylim()[1] * 0.96, "dimer", ha="center", va="top", fontsize=8)
    axis.set(xlabel="Read length (bp)", ylabel="Read count", title="Read-length distribution")
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def _save_assignment_plot(counts, output):
    labels, values = zip(*counts.items())
    total = sum(values)
    fig, axis = plt.subplots(figsize=(8, 5))
    bars = axis.bar(labels, values, color=[COLORS.get(label, plt.cm.tab10(i % 10)) for i, label in enumerate(labels)])
    for bar, value in zip(bars, values):
        axis.text(bar.get_x() + bar.get_width() / 2, value, f"{value}\n({value / total:.1%})", ha="center", va="bottom")
    axis.set(ylabel="Read count", title="Reads matched to candidate references")
    axis.margins(y=0.18)
    axis.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def _save_multimer_plot(multimers, output):
    plasmids = list(multimers)
    classes = ["1-mer", "2-mer", "3-mer", "4-mer", "other"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, metric, title in zip(axes, ("moles", "mass"), ("By molecule count", "By sequenced bases")):
        bottom = np.zeros(len(plasmids))
        for index, category in enumerate(classes):
            values = np.array([multimers[name][metric].get(category, 0.0) for name in plasmids])
            axis.bar(plasmids, values, bottom=bottom, label=category, color=plt.cm.viridis(index / 5))
            bottom += values
        axis.set(title=title, ylabel="Percent", ylim=(0, 100))
        axis.tick_params(axis="x", rotation=25)
    axes[1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.suptitle("Plasmid multimericity")
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def generate_qc_report(input_fastq, output_folder, reference_paths, assignment_path, host_reference=None, host_name="E. coli", host_key="ecoli", skip_host=False, multimer_tolerance=0.15):
    output_folder = Path(output_folder)
    report_folder = output_folder / "qc_report"
    report_folder.mkdir(exist_ok=True)
    assignments = _load_assignments(assignment_path) if Path(assignment_path).exists() else {}
    reference_lengths = _load_reference_lengths(reference_paths)
    opener = __import__("gzip").open if str(input_fastq).endswith(".gz") else open
    with opener(input_fastq, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
    if not records:
        raise ValueError("Cannot create QC report for an empty FASTQ")

    resolved_host, resolved_name = None, host_name
    if skip_host:
        host_reads, host_status = set(), "skipped by user"
    else:
        try:
            resolved_host, resolved_name = resolve_host_reference(output_folder, host_reference, host_name, host_key)
            host_reads, host_status = screen_host_reads(input_fastq, resolved_host)
        except Exception as error:
            host_reads, host_status = set(), f"unavailable: {error}"

    assignment_counts = Counter(row.get("assignment", "unassigned") for row in assignments.values())
    if not assignment_counts:
        assignment_counts["unassigned"] = len(records)
    multimers = {}
    for name, reference_length in reference_lengths.items():
        selected = [record for record in records if assignments.get(record.id, {}).get("assignment") == name and record.id not in host_reads]
        mole_counts, mass_counts = Counter(), Counter()
        for record in selected:
            category = classify_multimer(len(record.seq), reference_length, multimer_tolerance)
            mole_counts[category] += 1
            mass_counts[category] += len(record.seq)
        if not selected:
            continue
        multimers[name] = {
            "moles": {key: 100 * value / max(1, sum(mole_counts.values())) for key, value in mole_counts.items()},
            "mass": {key: 100 * value / max(1, sum(mass_counts.values())) for key, value in mass_counts.items()},
        }

    _save_length_plot(records, assignments, host_reads, reference_lengths, report_folder / "read_length_distribution.png")
    _save_assignment_plot(assignment_counts, report_folder / "reference_assignment_ratio.png")
    if multimers:
        _save_multimer_plot(multimers, report_folder / "multimericity.png")

    total = len(records)
    read_lengths = [len(record.seq) for record in records]
    median_length = int(np.median(read_lengths))
    cumulative = 0
    n50_target = sum(read_lengths) / 2
    for length_n50 in sorted(read_lengths, reverse=True):
        cumulative += length_n50
        if cumulative >= n50_target:
            break
    mean_quality = float(np.mean([np.mean(record.letter_annotations["phred_quality"]) for record in records]))
    host_percent = 100 * len(host_reads) / total
    host_available = host_status == "complete"
    host_display = f"{len(host_reads):,} ({host_percent:.1f}%)" if host_available else "Not assessed"
    table_rows = "".join(f"<tr><td>{html.escape(name)}</td><td>{count}</td><td>{count / total:.1%}</td><td>{reference_lengths.get(name, '—')}</td></tr>" for name, count in assignment_counts.items())
    multimer_rows = "".join(f"<tr><td>{html.escape(name)}</td><td>{data['moles'].get('1-mer', 0):.1f}%</td><td>{data['moles'].get('2-mer', 0):.1f}%</td><td>{data['mass'].get('1-mer', 0):.1f}%</td></tr>" for name, data in multimers.items())
    document = f"""<!doctype html><html><head><meta charset='utf-8'><title>Plasmid sequencing QC</title><style>
body{{font-family:Arial,sans-serif;max-width:1200px;margin:32px auto;color:#17202a}}h1,h2{{color:#173f5f}}.cards{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}}.card{{background:#eef5f9;border-radius:8px;padding:16px}}.value{{font-size:1.7rem;font-weight:bold}}img{{width:100%;border:1px solid #ddd;margin:10px 0}}table{{border-collapse:collapse;width:100%}}th,td{{padding:8px;border-bottom:1px solid #ddd;text-align:left}}.note{{color:#555;font-size:.9rem}}</style></head><body>
<h1>Plasmid sequencing QC report</h1><div class='cards'><div class='card'>Reads<div class='value'>{total:,}</div></div><div class='card'>Median length<div class='value'>{median_length:,} bp</div></div><div class='card'>Length N50<div class='value'>{length_n50:,} bp</div></div><div class='card'>Mean read quality<div class='value'>Q{mean_quality:.1f}</div></div><div class='card'>Host reads<div class='value'>{host_display}</div></div><div class='card'>Host reference<div class='value' style='font-size:1rem'>{html.escape(resolved_name)}</div></div></div>
<p class='note'>Host-screen status: {html.escape(host_status)}. A host read requires ≥50% query alignment and mapping quality ≥20.</p>
<h2>Read-length distribution</h2><img src='read_length_distribution.png' alt='Read length distribution'>
<h2>Reference assignment</h2><img src='reference_assignment_ratio.png' alt='Reference assignment ratio'><table><tr><th>Assignment</th><th>Reads</th><th>Percent</th><th>Reference length</th></tr>{table_rows}</table>
<h2>Multimericity</h2>{"<img src='multimericity.png' alt='Multimericity'>" if multimers else "<p>No reference-guided bins were available.</p>"}<table><tr><th>Plasmid</th><th>Monomer molecules</th><th>Dimer molecules</th><th>Monomer mass</th></tr>{multimer_rows}</table>
<p class='note'>Multimer classes are inferred from read length relative to each candidate reference (±{multimer_tolerance:.0%}). Results are screening metrics, not a substitute for assembly validation.</p></body></html>"""
    (report_folder / "index.html").write_text(document, encoding="utf-8")
    with open(report_folder / "summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows([("total_reads", total), ("median_read_length", median_length), ("read_length_n50", length_n50), ("mean_read_quality", f"{mean_quality:.3f}"), ("host_reads", len(host_reads) if host_available else "NA"), ("host_percent", f"{host_percent:.3f}" if host_available else "NA"), ("host_reference", resolved_name), ("host_screen_status", host_status)])
    return report_folder / "index.html"
