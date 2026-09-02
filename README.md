# Nanopore_plasmid_bin_assemble_polish

Converting whole-plasmid Nanopore reads to consensus FASTA files, including multiplexed
samples. When candidate plasmid references are available, reads can be separated using
plasmid-specific sequence rather than read length alone.

## Local installation instructions
Install the conda environment with 
```
conda env create -n Nanopore --file=Nanopore.yaml
```

Packages are listed in 
```
Nanopore.txt
```

Canu also needs to be installed from the binary distribution: https://github.com/marbl/canu/releases

To run the script:
```
conda activate Nanopore
python nanopore_plasmid_bin_assemble_polish.py -i <input_fastq_file> -o <output_folder> --canu_binary_path <path_to_canu_binary> -t <number_of_threads_for_processing>
```

### Recommended: reference-guided binning

Supply two or more candidate plasmids as GenBank or FASTA files. The pipeline derives
strand-independent k-mers that occur in only one supplied reference. A read is assigned
only when it has enough unique sequence evidence and the best reference exceeds the
runner-up by the configured ratio.

```bash
python nanopore_plasmid_bin_assemble_polish.py \
  -i mixed_reads.fastq.gz \
  -o output \
  -r plasmid_A.gbk plasmid_B.gbk \
  -c canu-2.2/bin/canu \
  -t 8
```

Use `--binning-only` to inspect assignments without running Canu or Medaka. Outputs
include one FASTQ per assigned plasmid, `ambiguous.fastq`, and
`read_assignments.tsv`. For GenBank references, the report names unique annotated
origin and antibiotic-resistance features supporting each assignment when present.

### Graphical QC report

Every run creates `qc_report/index.html` plus publication-ready PNG and TSV files:

- stacked read-length distribution, colored by plasmid assignment and host reads;
- candidate-reference read counts and percentages;
- monomer, dimer, trimer, tetramer, and off-size fractions by molecule count and mass;
- host-genome contamination summary; and
- machine-readable summary statistics.

Host screening uses minimap2 through `mappy`. By default the workflow downloads and
caches the NCBI *E. coli* K-12 MG1655 reference (`NC_000913.3`). For another cloning
host, select one of the bundled choices:

```bash
--host ecoli
--host vibrio-natriegens
--host bacillus-subtilis
```

The *V. natriegens* ATCC 14048 reference contains chromosomes `CP009977.1` and
`CP009978.1`; the *B. subtilis* 168 reference uses `NC_000964.3`. Alternatively,
provide `--host-reference host.fasta --host-name "host strain"`. A custom reference
overrides `--host`. Use
`--skip-host-screen` for offline runs or `--no-qc-report` to disable reporting.

### No-reference plasmid discovery

For high-accuracy, approximately full-length plasmid reads without candidate
references, combine technical preparations and discover the two strongest distinct
sequence families:

```bash
python no_reference_plasmid_analysis.py \
  extraction_1.fastq.gz extraction_2.fastq.gz \
  -o no_reference_output \
  --top 2
```

The workflow detects sharp length peaks, clusters them by sequence, collapses multimer
peaks into their corresponding monomer, reconstructs a de Bruijn consensus when the
circular path is unambiguous, and ranks candidates by circularity and sequence support.
It writes GenBank, FASTA, a TSV summary, and a comparative length-distribution plot.
Short degradation products below 2 kb are not considered plasmid candidates.

No-reference GenBank outputs are deliberately minimally annotated. A `linear` topology
means the sequence family is well supported but circular closure was not proven; it
must not be relabeled circular without an independent assembly or junction-spanning
evidence.

Example QC-only run:

```bash
python nanopore_plasmid_bin_assemble_polish.py \
  -i mixed_reads.fastq.gz \
  -o output \
  -r plasmid_A.gbk plasmid_B.gbk \
  --binning-only
```

Classifier controls:

- `--kmer-size 15`: signature k-mer length.
- `--min-unique-kmers 10`: minimum unique k-mer hits required.
- `--min-assignment-ratio 2.0`: minimum best/second-best hit ratio.

Reads containing only backbone shared by multiple candidates are retained in the
ambiguous bin. Candidate references must include every plasmid that should receive a
named bin; otherwise reads from an unknown plasmid may remain ambiguous or be assigned
to a similar supplied reference.

Output:
1) histogram of read lengths. The cutoff for binning different plasmids is marked by a gray bar.
2) polished assemblies are named "consensus.fasta" in the "polished_output" folders.

Pipeline description:
1) assign reads by reference-specific sequence when `--references` is supplied;
   otherwise use the legacy length histogram.
2) assemble each read bin separately with Canu.
3) polish each assembly with only the reads assigned to that bin using Medaka.

Notes:  
  
Bin width is 200 bp.  
Binning cutoff is 3 std above the mean read length.  
Low-coverage bins may not provide enough reads for Canu assembly.



## Docker instructions
Make sure "Use Virtualization framework" is enabled in Docker Settings > General

Also make sure that Rosetta emulation is turned off in Docker Settings > Features in development

1. Build the docker image using:
```
# docker build -t {name of image} .
docker build -t nanopore .
```

2. Create and run the docker container using:
```
# docker run -it --rm --platform linux/amd64 --name {container name} -v {Path to your github repo that is being mounted}:/{where it's being mounted inside} {image name}
docker run -it --rm --platform linux/amd64 --name 241019_test -v /Users/joncchen/Documents/GitHub/Nanopore_plasmid_bin_assemble_polish:/app/ nanopore
```

3. Activate the conda environment inside the container
```
conda activate Nanopore
```

4. Run the command to initiate assembly + polishing.
```
# python nanopore_plasmid_bin_assemble_polish -i {input file location} -o {output folder} --canu_binary_path canu-2.2/bin/canu -t {# of CPU threads}
python nanopore_plasmid_bin_assemble_polish.py -i input/001_347A_reads.fastq.gz -o output --canu_binary_path canu-2.2/bin/canu -t 8
```
