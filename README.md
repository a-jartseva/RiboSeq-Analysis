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

If changes have been made since you last cloned the repository, you can keep your local copy up-to-date by running:

cd into/cloned/fork-repo
git remote add upstream https://github.com/adamd3/RiboSeq-Analysis.git
git fetch upstream
git pull upstream master

**Running the pipeline**

Place all of your (gzipped) Fastq format files in a single directory, and create a file called "libraries.txt" which contains a description of the libraries in your analysis. 

The libraries.txt file should follow the format of the below example:

```
Index1       MuLV    Rattus_norvegicus       mock            RiboSeq-CHX 
Index2       MuLV    Rattus_norvegicus       mock            RNASeq
Index3       MuLV    Rattus_norvegicus       infected        RiboSeq-CHX 
Index4       MuLV    Rattus_norvegicus       infected        RNASeq
```

The first column gives the library name (must be the same as the Fastq file name), the second is the virus strain (must be the same as the bowtie index name for that virus), the third is the host species name (again, must be the same as the bowtie index name for that species), fourth and fifth columns are self-explanatory. 

You are now ready to start running the pipeline.
