# Staden Sequence Analysis Package

## What is this?
The Staden Package bundles classic desktop tools that help turn raw DNA sequencing traces into trustworthy consensus sequences.  Imagine a lab bench covered in printouts and coloured electropherogram peaks: Staden was the software that let scientists view those traces, trim away poor-quality regions, align overlapping reads, and agree on the final DNA sequence before publishing it.

## Typical Workflow: Inputs and Outputs

| Stage                      | Load This                                           | Source                                           | Export Produces                        | Why It Matters                                                                 |
|----------------------------|----------------------------------------------------|--------------------------------------------------|----------------------------------------|-------------------------------------------------------------------------------|
| **Trace Viewer**<br>`trev` | Chromatogram files (`.ab1`, `.scf`, `.abi`)        | Sanger sequencing machines or core facilities    | Cleaned traces, per-read quality notes | Inspect and trim reads; flag suspicious peaks.                                 |
| **Pre-processing**<br>`pregap4` | Trace files & metadata (`.exp`, `.xml` folders)     | Lab information systems or sequencing run output | Trimmed FASTA/FASTQ, clip reports      | Produce high-quality reads for assembly and clip away contamination.           |
| **Assembly Editing**<br>`gap4`, `gap5` | Read files, trace links, assembly databases (`.gap`, `.g5d`) | Output from `pregap4` or other assemblers        | Curated contigs, consensus FASTA, annotated databases | Resolve discrepancies and produce the trusted consensus sequence.              |
| **Comparative Analysis**<br>`spin` | Consensus sequences, reference genomes            | Gap4/Gap5 assemblies or external pipelines       | Dot plots, restriction maps, reports   | Visualize differences, design experiments, clarify relationships.              |


### What do people do with the outputs?
- Submit curated consensus sequences to public DNA repositories (GenBank/ENA/DDBJ).
- Produce figures showing assembly quality, coverage, or remaining ambiguities.
- Share the project databases with collaborators so they can review the evidence supporting each base call.
- Feed the polished sequence into downstream tasks such as primer design, mutation confirmation, or cloning strategies.

## Who finds it useful today?
- **Genome centres and sequencing labs** maintaining legacy Sanger pipelines built around Gap4, Gap5, Pregap4, Trev, Spin, or related utilities.
- **Researchers and educators** demonstrating how sequencing analysis was performed before modern high-throughput platforms.
- **Bioinformatics platform engineers** who must re-run or audit projects originally analysed with the Staden toolkit.

Although modern workflows usually run in cloud environments, the Staden applications remain valuable when you need an interactive, GUI-based review of sequencing data—especially when reopening older projects with Sanger reads or hybrid assemblies.

## Getting started
1. Review the build requirements and platform-specific notes in [`BUILD.md`](BUILD.md).
2. Compile the package in a Unix-like environment (Linux, BSD, macOS, or Windows via MinGW/MSYS).
3. Launch the application that matches your task:
   - `gap4` / `gap5` for assembly review and editing
   - `pregap4` for preprocessing reads
   - `trev` for viewing electropherogram traces
   - `spin` for comparative sequence analysis

## Contributing and support
This repository preserves the historic code base.  Modernising it for new platforms often requires updates to build scripts and GUI dependencies.  If you port the package to contemporary toolchains, please consider contributing patches or notes back through pull requests or the issue tracker.

## Project status & modern context
- **Last significant upstream activity**: circa 2016. No active maintainer is known.
- **Legacy value**: reopening historic assemblies, reproducing published results, teaching classic sequencing pipelines.
- **Modern alternatives**: actively maintained tools such as SPAdes, Flye, Canu, Velvet (assembly) and IGV, Tablet, Bandage, Artemis/ACT (visualisation/curation) provide easier installs for current projects.
