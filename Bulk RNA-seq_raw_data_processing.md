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
