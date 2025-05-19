# Tutorial for bulk RNA-seq data preprocessing and analysis
#### Updated on 09 May 2025
### Table of Content
  * [Introduction](#introduction)
  * [Preparation](#preparation)
    * [1-1. Linux and Bash, your important buddies for RNA-seq data analysis](#1-1-linux-and-bash-your-important-buddies-for-rna-seq-data-analysis)
    * [1-2. Access the computing server](#1-2-access-the-computing-server)
    * [1-3. Using a computer cluster like Euler](#1-3-using-a-computer-cluster-like-euler)
    * [1-4. Install the required tools with the help from conda](#1-4-install-the-required-tools-with-the-help-from-conda)
    * [1-5. Get the public RNA-seq data from SRA](#1-5-get-the-public-rna-seq-data-from-sra)
  * [Preprocessing of RNA-seq data](#preprocessing-of-rna-seq-data)
    * [2-1 Quality control of RNA-seq data](#2-1-quality-control-of-rna-seq-data)
    * [2-2 Read mapping/pseudomapping and quantification](#2-2-read-mappingpseudomapping-and-quantification)
      * [2-2-1 Read mapping with STAR and data quantification](#2-2-1-read-mapping-with-star-and-data-quantification)
      * [2-2-2 Read pseudomapping with kallisto and data quantification](#2-2-2-read-pseudomapping-with-kallisto-and-data-quantification)


## Preparation
<sub><a href="#top">(Back to top)</a></sub></br>
In this section, I will discuss basics before we even start preprocessing the RNA-seq data, so that you can make sure that you and your computer are both ready for the following steps. These basics include the following:
Check the detailed steps of 1-5 from the notebook of QuaDBio (https://github.com/quadbio/RNAseq_tutorial/blob/main/Tutorial.md)
1. Linux and [Bash (the Bourne Again SHell)](https://www.gnu.org/software/bash/), the most commonly used command line interface in Linux
2. How to access a computing server via SSH
3. How to use an HPC like Euler
4. Conda and how to use it to install software needed for the data preprocessing

5. SRA and how to retrieve public RNA-seq data
### 1-5. Get the public RNA-seq data from SRA
<sub><a href="#top">(Back to top)</a></sub></br>
Now we have the computational environment ready for the preprocessing. We just need to data to start. As one can easily expect, there are two sources of data: the in-house data that are generated freshly by the lab for specific biological questions, and the public data which have been released and used for answer certain questions, but can be reanalyzed solely or together with other data for the same or related questions.

There are several huge repositories in the world for high-throughput sequencing data. The major players include [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) by NCBI in the US, [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) by EMBL-EBI in the UK, and [DDBJ Sequence Read Archive (DRA)](https://www.ddbj.nig.ac.jp/dra/index-e.html) by the DDBJ Center in Japan. These three repositories are also members of International Nucleotide Sequence Database Collaboration (INSDC), and cross-backup each other, meaning that data submitted to any one of the three databases are also accessible from the other two.

>**NOTE**
>While majority of the high-throughput sequencing data are archived in these three databases, there are also other emerging sequencing data repositories, though most of them are regional, mostly used by researchers in the host country. Examples include [Genome Sequencing Archive (GSA)](https://ngdc.cncb.ac.cn/gsa/) by NGDC, and [CNGA Sequence Archive (CNSA)](https://db.cngb.org/cnsa/) by CNGB, both located in China.

In this tutorial, we will retrieve the data we need from SRA. This is not only because the data we are going to use are archived at SRA, but also because of SRA-Toolkit which provides a simple way to download data given the accession numbers of the needed samples.

- SRA Toolkit: OSC gives a software usage of SRA (https://www.osc.edu/resources/available_software/software_list/sra_toolkit)
1) Open the terminal of Ascend
2) View available modules:
   ```console
   module spider sratoolkit
   ```
3) Create the bash file for downloading by SRA:
   ```console
   nano download_fastq.sh
   ```
   Paste your script (like the one in your canvas).
   
   ```console
   #!/bin/bash
   #SBATCH --job-name=download_fastq
   #SBATCH --account PAS2556
   #SBATCH --time=05:00:00
   #SBATCH --ntasks=1
   #SBATCH --output=download_fastq_%j.log
   #SBATCH --error=download_fastq_%j.err
   
   module load sratoolkit/3.0.2
   
   # List of accession numbers
   ACCESSIONS=("SRR26891264" "SRR26891265" "SRR26891266")
   
   # Download 'prefetch' and convert 'fasterq-dump' each dataset into FASTQ format in the current directory
   for ACC in "${ACCESSIONS[@]}"; do
       echo "Downloading $ACC..."
       prefetch "$ACC" --output-directory .
   
       echo "Converting $ACC to FASTQ..."
       fasterq-dump "$ACC" -O .
   done
   ```
   #Save and exit from nano: Press Ctrl+O → Enter to save. Press Ctrl+X to exit.

5) Submit the bash file:
   ```console
   sbatch download_fastq.sh
   squeue -u osc_username
   ```
      
  
   Once the download is finished, you can list the files in your working directory and see whether you have all the files as expected ('wc -l'). They should all be named as [SRR Accession].fastq
   
   ```console
   [qujia93@ascend-login02 Raw.data]$ ls -l *.fastq
   -rw-rw----+ 1 qujia93 PAS2015 21969788302 May  7 16:12 SRR26891264_1.fastq
   -rw-rw----+ 1 qujia93 PAS2015 21969788302 May  7 16:12 SRR26891264_2.fastq
   -rw-rw----+ 1 qujia93 PAS2015 19084215314 May  7 15:04 SRR26891265_1.fastq
   ...[skip the other lines]...

   [qujia93@ascend-login02 Raw.data]$ ls -l *.fastq | wc -l
   12
   ```


#### An introduction to the FASTQ and FASTA data formats
FASTQ is the standard data format used to store sequencing data. It was originally developed at the Wellcome Trust Sanger Institute to incorporate quality information into the FASTA format, which is the standard data format developed to represent nucleotide or amino acid sequences where every nucleotide or amino acid is represented by one single character (e.g. A/T/C/G for DNA sequences).

In a FASTA file, one or multiple sequences are represented, each in a unit of two components. The first component is strictly one line called "description line", which starts with the character ">". The description line gives a name and/or a unique identifier for the sequence, and potentially also additional information. The second component is the sequences in one or multiple lines, one character for one nucleotide or amino acid. In genomics, the most commonly usage of FASTA is to store the sequences of genome, annotated transcriptome and peptide products. The FASTA file looks like the following:
```
>sequence A
GGTAAGTCCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATAT
TCTGTTGCCAGAAAAAACACTTTTAGGCTATATTAGAGCCATCTTCTTTGAAGCGTTGTC
>sequence B
GGTAAGTGCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATAT
TCTGTTGCCAGATTTTACACTTTTAGGCTATATTAGAGCCATCTTCTTTGAAGCGTTGTC
TATGCATCGATCGACGACTG
```

A FASTQ file looks very relevant but different from a FASTA file. It is used to describe one or multiple sequences, but in most of the time multiple (usually tremendous amount of) sequences. Every sequence is represented in a unit of four lines:
* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

It looks like the following:
```
@SRR2815952.1 1 length=100
AGACGAGACCTACTGCATTGATAACGAAGCTCTCTACGACATTTGCTTCAGAACCCTAAAGCTGACCACGCCCACCTATGGTGACCTGAACCACCTGGTG
+SRR2815952.1 1 length=100
@?@B?@DDDDHDCGHGHIIEHGIIFEHI@?BFGABDFGHAGGIIIHIIIGGIIIBHIIFGIHACHHEEBEBCCDD@ACC:>>@CDDD?CCD(<<?A?A@C
@SRR2815952.2 2 length=100
TGGGGTTTCACCATGTTGGCCGGGCTGGCCTCGAACTCCTGACCTTGTGATGCACCCACCTCGGCCTCCCAAAGTGCTGGGATTTACAGGCGTAAGCCAC
+SRR2815952.2 2 length=100
CCCFFDFFHHHHHJJJJJJJJHGIJJJJJJJJJIJJIJIIGIJJJJJIIJJJJJJJHHHFFFFDDDDDDDDDBC@>@DDDDDDDDDEDCCBDDDDDDDDD
```

FASTQ is a great format to represent short sequences with quality values, and therefore high-throughput sequencing data where the quality values can represent the sequencing quality. In the Illumina sequencing, a Q score is used to represent the quality of each sequenced nucleotide. The Q score is defined as
$$Q = -10 \log_{10} e$$
Here, $e$ is the estimated probability of the base call being wrong. Therefore, a higher Q score means lower probability of sequencing error (or more precisely base calling error), and therefore higher quality. A Q score of 20 (Q20) means an error rate of 1 in 100, or accuracy of 99%. A Q score of 30 (Q30) means accuracy of 99.9%. At such a level, the base is considered to be accurate enough that all bases of a read is likely to be correctly called if all the bases in a read reach such a level of quality. Currently Q30 is considered a benchmark for quality in next-generation sequencing (NGS).

To represent the Q scores of all the based in a sequence in the FASTQ format, the numeric Q score of each base is encoded into a compact form based on the [ASCII codes](https://www.ascii-code.com/). Basically, the estimated Q score represented as an integer is represented by the character with the ASCII code equal to $Q+33$. It has the 33 component because Q is strictly non-negative with the minimum of 0 meaning 0% accuracy, and 33 is the smallest ASCII code that defines a one-character symbol (!). 

## Preprocessing of RNA-seq data
<sub><a href="#top">(Back to top)</a></sub></br>
Now we have the tools ready, and the data ready. It is time to move on to the next step, to preprocess the RNA-seq data. In general it contains the following steps:
1. Quality control
2. Read mapping or pseudomapping
3. Gene expression quantification for samples
4. Generate the expression matrix with sample metadata table for the following analysis

(If you want to use the specific tool for each step, please check https://github.com/quadbio/RNAseq_tutorial/blob/main/Tutorial.md?plain=1), introducing more detailed information about RNA-seq data preprocessing.)

Here, we introduce **nf-core/rnaseq**, a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. (https://nf-co.re/rnaseq/3.14.0/) 

It takes a samplesheet and FASTQ files as input, **performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report**.

# Data
RNASeq data (fastq.gz) 
>**NOTE**
> If the raw data is stored as .fastq files (e.g. uncompressed by sratoolkit), the .fastq files should be compressed to .fastq.gz files:
``` console
#!/bin/bash
#SBATCH --job-name=gzip
#SBATCH --output="%j_log.txt"
#SBATCH --account=PAS2556
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=40:00:00

gzip *.fastq
```

# Input format
Prepare a **samplesheet.csv** file with your input data (use 'auto' if you do not know the strandedness):
```console
sample,fastq_1,fastq_2,strandedness
RNAseq_naive_WT1,"	SRR26891269_R1.fastq.gz","	SRR26891269_R2.fastq.gz",auto
RNAseq_naive_WT2,"	SRR26891268_R1.fastq.gz","	SRR26891268_R2.fastq.gz",auto
RNAseq_naive_WT3,"	SRR26891267_R1.fastq.gz","	SRR26891267_R2.fastq.gz",auto
```
Each row represents a fastq file (single-end) or a pair of fastq files (paired-end). Rows with the same sample identifier are considered technical replicates and merged automatically. The strandedness refers to the library preparation and will be automatically inferred if set to auto.
>Warning: Please provide pipeline parameters via the CLI or Nextflow -params-file option. Custom config files, including those provided by the -c Nextflow option, can be used to >provide any configuration except for parameters; see docs.

2) Then, you can run the pipeline using:
First, checking the version of nf-core/rnaseq in OSC
``` console
module spider nextflow
```
Then, create a .sh file and specify the current version of nextflow (module load nextflow/24.10.4
) in your .sh file:
```console
   nano BulkRNA_Alignment.sh
```
```console
#!/bin/bash
#SBATCH --job-name=BulkRNA_Alignment
#SBATCH --output="%j_log.txt"
#SBATCH --account=PAS2556
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=40:00:00

module load nextflow/24.10.4

nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir /fs/ess/PAS2556/Bioinformatics_analysis/BulkRNA/Data/05082025/Alignment_output  \
    --genome GRCm38 \
-profile singularity
```
Submit sbatch file
```console
sbatch BulkRNA_Alignment.sh
squeue -u osc_username # check the running job 
```
For more details and further functionality, please refer to the usage documentation and the parameter documentation.

# Pipeline output
To see the results of an example test run with a full-size dataset, refer to the results tab on the nf-core website pipeline page. For more details about the output files and reports, please refer to the output documentation. (https://nf-co.re/rnaseq/output)

This pipeline quantifies RNA-sequenced reads relative to genes/transcripts in the genome and normalizes the resulting data. It does not compare the samples statistically in order to assign significance in the form of FDR or P-values. For downstream analyses, the output files from this pipeline can be analysed directly in statistical environments like R, Julia or via the nf-core/differentialabundance pipeline.
