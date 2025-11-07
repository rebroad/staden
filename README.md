# Staden Sequence Analysis Package

## What is this?
The Staden Package is a collection of classic desktop tools for turning raw DNA sequencing reads into polished assemblies.  Think of it as a well-stocked workbench: it shows you the original trace files, helps you join overlapping reads, highlights possible sequencing errors, and lets you tidy up the final consensus sequence before it is published.

## Who finds it useful?
- **Genome centres and sequencing labs** maintaining legacy pipelines built around Gap4, Gap5, Pregap4, Trev, Spin, or related utilities.
- **Researchers and educators** who want to demonstrate how traditional Sanger-era assembly and quality control were performed.
- **Bioinformatics platform engineers** who need to inspect or reproduce results generated with the original Staden tools.

Although modern high-throughput workflows now run in the cloud, the Staden applications remain valuable when you need an interactive, GUI-based review of sequencing data—especially when revisiting older projects with Sanger reads or hybrid assemblies.

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
- **Modern alternatives**: tools such as SPAdes, Flye, Canu, Velvet (assembly) and IGV, Tablet, Bandage, Artemis/ACT (visualisation/curation) are actively maintained and easier to install.
