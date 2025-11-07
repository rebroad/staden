# Staden Sequence Analysis Package

## What is this?
The Staden Package is a collection of desktop programs created for the era when DNA sequencing machines spat out individual read files rather than billions of reads in one go. Each read looked like a colourful graph of peaks (an electropherogram) plus a machine-guessed list of letters. Those raw files were noisy, misaligned, and often contaminated. Staden gave scientists the tools to clean up that mess, line the reads up against each other, and agree on the final DNA sequence that would go into a paper or a public database.

## Why are the inputs messy?
Modern DNA sequencers still start by recording the light coming off a dye-labelled molecule. The machine saves that raw signal as a chromatogram (`.ab1`, `.scf`, `.abi`):
- Peaks may overlap, especially near the end of a read.
- Background dye or sequencing primer can introduce junk bases.
- The machine assigns an A/C/G/T call, but with varying confidence.

On top of the trace files, labs record metadata describing the experiment (sample name, primers used, expected vector) so that analysts can keep track of which read belongs to which sample.

Without Staden—or a tool like it—researchers would have to open each raw trace in a plotting program, manually compare overlapping reads, and somehow keep track of every edit in a spreadsheet. Staden automates those comparisons, keeps audit trails, and outputs a polished consensus sequence that collaborators and public archives can trust.

## Typical workflow: inputs and outputs
| Step | What you open (lay description) | File examples | Why Staden is useful | What you export (lay description) |
|------|---------------------------------|---------------|----------------------|-----------------------------------|
| **1. Trace viewer (`trev`)** | The raw “squiggle graphs” straight from the sequencing machine. | `.ab1`, `.scf`, `.abi` | Shows the peaks, lets you trim off low-quality ends, and mark doubtful bases. | Cleaned trace files and notes about which positions look suspicious. |
| **2. Pre-processing (`pregap4`)** | Folders of traces plus the notebook information about each sample (who it belongs to, which primer was used, etc.). | Trace folders + metadata (`.exp`, `.xml`) | Batch-trims the junk (vectors, primers), rejigs filenames, and produces tidy read files ready for assembly. | FASTA/FASTQ reads with clipping reports so you know what was removed. |
| **3. Assembly editing (`gap4`, `gap5`)** | The cleaned reads lined up so overlapping regions can be compared side-by-side. | Reads + project databases (`.gap`, `.g5d`) | Highlights disagreements between reads, links back to the original trace, and records every manual decision. | A curated consensus sequence, quality scores, and an assembly project file others can review. |
| **4. Comparative analysis (`spin`)** | The agreed consensus sequence versus a reference or another clone. | Consensus FASTA, reference genomes | Generates dot plots, maps restriction sites, translates coding regions—useful for clone verification or mutation analysis. | Figures, reports, and additional annotations to share with colleagues. |

### Why the outputs are better than the inputs
- **Confidence**: Instead of trusting a single noisy read, you end up with a consensus backed by multiple traces and reviewer notes.
- **Traceability**: Every base call links back to the original raw data so collaborators can reproduce decisions.
- **Compatibility**: Public archives expect clean FASTA/FASTQ files and project metadata; Staden produces those in standard formats.

### What people do with the outputs
- Submit the cleaned consensus sequence to GenBank/ENA/DDBJ.
- Document clone verification or mutation detection in lab notebooks and reports.
- Share the project database so collaborators can audit the evidence for important base calls.
- Use the curated sequence as the starting point for primer design, vector construction, or downstream analyses.

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
