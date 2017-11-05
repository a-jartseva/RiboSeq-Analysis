# Analysis of Ribosome Profiling and RNA sequence data for virus infections of host cells.

**Overview of pipeline:**

- Quality assessment and trimming of reads.
- Map sequentially to rRNA; vRNA; mRNA; ncRNA; host genome.
- Plot the read mapping statistics.
- Plot length distribution, phasing, and positions of reads relative to start and stop codons on host RefSeq mRNAs.
- Assess contamination of samples by ribonucleoproteins (RNPs)
- For virus-infected samples, make combined plots showing length distribution of reads mapping to virus and host coding sequences.

**Pre-requisites:**

Before beginning, ensure that the following programs are installed:

- FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
- bowtie v1 (http://bowtie-bio.sourceforge.net)
- STAR (https://github.com/alexdobin/STAR)
- GNU parallel (https://www.gnu.org/software/parallel/)
- samtools (http://samtools.sourceforge.net)

Additionally, the following R packages must be installed for plotting:

- grid
- gridExtra
- reshape2
- ggplot2

**Cloning the Git respository:**

To get started, make a local copy of the Git repository:

     git clone https://github.com/adamd3/RiboSeq-Analysis.git

**Running the pipeline**

Place all of your (gzipped) Fastq format files in a single directory, and create a file called "libraries.txt" which contains a description of the libraries in your analysis. 

Follow the structure in the below example when creating your libraries.txt file:



You are now ready to start running the pipeline.
