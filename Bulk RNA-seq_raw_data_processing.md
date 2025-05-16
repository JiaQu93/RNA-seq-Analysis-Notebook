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


### 2-1 Quality control of RNA-seq data
<sub><a href="#top">(Back to top)</a></sub></br>
Before actually processing the data, it is important to make sure that the data is of high quality.

The quality of a RNA-seq data contains different aspects. The sequencing quality, represented by the quality score for each base in each read, is one of the very first thing one should consider. It varies from one base to another base, and from one read to another. We can therefore look at the distribution of the Q score per loci across all reads. Usually the sequencing reads have relatively low quality on the two sides (especially the end) and higher quality in the middle, and this is the reason why we look at the per loci quality distribution across reads. Another reason to look at that is because if it is indeed the case that the start and/or the end of reads systematically have too low quality, one can easily fix it by trimming the first and the last bases of all the reads. This will be mentioned in more details later.

Other quality metrics are more related to the sample and cDNA library quality. For instance, the adapter ligantion is a critical part during the cDNA library preparation and cases can happen that multiple adapters are ligated to the sequence and therefore become parts of the sequenced reads. This would introduce troubles later when we need to locate the transcript(s) that the read represents, as the extra adapter sequences are not a part of the transcripts and would introduce a large mismatch. Another example is the complexity of the data. Ribosomal RNA (rRNA) makes up about 80% of cellular RNA, while the rRNA genes only makes up 0.5% of the human genome, and their abundances are not very relevant to many biological processes to study. Therefore, there is usually a step of mRNA enrichment (by Oligo-T sequences) or rRNA depletion to effectively capture the more informative non-rRNA transcript fractions. However, this is not always working efficiently and the cDNA library would present a low complexity if the rRNA portion is not effectively reduced. Also, there could be reads representing pure ligation products of multiple sequencing adapter, which also dramatically reduce the library complexity.

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a tool providing a simple way to do some quality control checks on the sequencing data. It checks different aspect of data quality and provides a graphical report so that one can intuitively get the idea about the data quality. The tool contains both GUI and CLI, therefore one can choose to either run it by clicking buttons with mouse or with commands.

To run it in the command line, it is very simple. This is an example:
```console
cd [the scratch folder]
cd rawdata
mkdir fastqc
fastqc -o fastqc *.fastq.gz
```
>**NOTE**
>The `-o` option in the `fastqc` command specifies the output directory to put the output files. By default it stores the output files to the same folder as the FASTQ files. In the example script, a new folder is created to store just the FastQC reports. This is not a must, but will likely help to make it more tidy.

If you prefer GUI, you can also directly run `fastqc` to open the GUI in the server. If everything goes well you would be able to see the FastQC window poping up. Then you can open one or more FASTQ files by selecting from the menu `File` -> `Open...`. After getting the report that's shown in the window, you can save it, also two files just as the command line version, by choosing from the menu `File` -> `Save report...`.

>**NOTE**
>Please note that to see the GUI at your personal computer, you need to make the X11 forwarding possible. [X11](https://en.wikipedia.org/wiki/X_Window_System) is a windowing system for bitmap displays, proving the basic framework for a GUI environment and is commonly used in Linux. To allow the X11 to work via SSH-based remote access, there are several things one would need to do:
>* For macOS users, make sure to add the `-Y` option when running the `ssh` command. It means to use `ssh -Y <username>@bs-studentsvr04` to login the bs-studentsvr04 server, for instance. Also, you need to make sure that XQuartz, the X11 system that runs on macOS, is installed.
>* For Windows users using PuTTY, it is more complicated. First of all, you need to make sure that [Xming](http://www.straightrunning.com/XmingNotes/) or any other X11 server for Windows is installed, and being opened. And then in PuTTY, after filling in the server information at the starting page, from the option menu on the left hand side, choose "Connection" -> "SSH" -> "X11". There you shall click the "Enable X11 forwarding". Afterwards click the "Open" button to start the connection.
>* The newest Xming is not free, but needs £10 "donation" to get the download password. Its older version, however, is available at [sourceforge.net](https://sourceforge.net/projects/xming/files/Xming/6.9.0.31/Xming-6-9-0-31-setup.exe/download) for free.

For each FASTQ file there are two output files generated by default. Assume the FASTQ file is called `[filename].fastq` or `[filename].fastq.gz`, then one output file is called "[filename]_fastqc.html" and the other one called "[filename]_fastqc.zip". The HTML file is a report which can be opened with your browser, and it contains the summary of all the quality metrics, as well as a grade by the software on each perspective whether it is passed, warning, or failed. The ZIP file, once decompressed, contains also the HTML report, as well as plain text files with the same information.

This is an example screenshot of the HTML report when being opened in the browser:

<p align="center"><img src="img/fastqc.png" /></p>

It is important to mention that some of the grades are made assuming the data to be whole-genome DNA sequencing. For instance, the "Per sequence GC content" compares the distribution of G/C bases proportion per read to a theoretical distribution derived from the whole genome, which is expected to be different from the GC content of transcriptome. From my personal experience, the "Per base sequence content" and "Per sequence GC content" are the two sections that easily get the warning for failed grade for RNA-seq data, but can be ignored if other sections are fine. In addition, the "Sequence Duplication Levels" is another section that could give out warning of RNA-seq data, while it may or may not be a problem that needs to be solved later.

Meanwhile, the sections that I would suggest to pay attention to for RNA-seq data include "Per base sequence quality", "Sequence Duplication Levels", "Overrepresented sequences" and "Adapter Content". They represent potential We will need to try to fix the problem if they get a failed grade:
* If any read locus shows low quality (e.g. median <20 or even <10) in the "Per base sequence quality" section, especially at the two ends, we should try to trim them if the low-quality part is large (>10 bases), either by all reads removing the same number of bases or different number per read based on the quality scores.
* Since different transcripts have very different abundance, to make sure that very lowly expressed transcripts are also detected, it is possible that the highly expressed transcripts are over-amplified and/or over-sequenced, resulting in warning of "Sequence Duplication Levels". In this case, a de-duplication step may be wanted to collapse the identical reads into one.
* For standard mRNA-seq data with oligoT enrichment, problems of "Overrepresented sequences" and "Adapter Content" often come together and represent the adapter ligation issue mentioned above. We can try to cut the adapter sequences from reads later.

For the example shown in the screenshot above, we don't need to do anything as it looks all good.

**IMPORTANT NOTE**: It is not always necessary to do anything here even if problems were found, especially those related to base quality. For instance, many up-to-date software being used later for read mapping (e.g. STAR) has implemented a soft trimming mechanism to deal with low-quality bases at the end of a read.


### 2-2 Read mapping/pseudomapping and quantification
<sub><a href="#top">(Back to top)</a></sub></br>
#### 2-2-1 Read mapping and data quantification 
(If you want to use the STAR tool, please check https://github.com/quadbio/RNAseq_tutorial/blob/main/Tutorial.md?plain=1)

Here, we introduce nf-core/rnaseq, a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.

# Data
RNASeq data (fastq or fastq.gz)
# Input format
Prepare a samplesheet.csv file with your input data that looks as follows (you can use 'auto' if you do not know the strandedness):
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,auto



