#Estimating abundances in metagenomes

Authored by Jin Choi for EDAMAME2016 

## Summary

### Overarching Goal
* This tutorial will contribute towards an understanding of quantitative analyses of **metagenome data**

### Learning Objectives
* Understanding how to estimate abundances of reads in a representative gene reference
* Understanding read mapping
* Understanding mapping file formats
* Understanding how to use a mapping program (Bowtie2, samtools, bcftools)
* Apply reference mapping to assess read abundances and quantify gene presence

### Reference mapping is useful... 
* If you want to detect SNPs
* If you want to estimate abundance of genes in metagenomic or metatranscriptomic data

## Tutorial

@Jin - add AMI information or what the assumption of the system they are on -- base Ubuntu I think
@Jin - add an overview - what is the point, objective?  Imagine you have a metagenome created to study and the reference is...etc.
@Jin - figure out dataset for assembly tutorials streamlining - this tutorial naturally goes after the assembly tutorial

### Setting up operating system
```
apt-get update
apt-get -y install screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat
```

## Install mapping software for this tutorial, Bowtie2 and SamTools
```
cd 
wget wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip
unzip bowtie2-2.2.9-linux-x86_64.zip
mv bowtie2-2.2.9 BT2_HOME
```
@Jin - why use bowtie, what is it?  link to its documentation?

Set up path (you need to set up path again if you are logged in later):
```
PATH=$PATH:~/BT2_HOME
export PATH
```

Install SamTools.  
```
apt-get -y install samtools
```
@Jin - add what is samtools and why use it

## Download data
First, download reference sequence.  @Jin add more detail here, what is the reference.

download sequencing file: this file is the first 100,000 sequences of original file
```
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098038/SRR098038.fastq.gz
```

## Do the mapping
Now let’s map all of the reads to the reference. Start by indexing the reference genome:
```
cd ~/metagenome

bowtie2-build megahit_out/final.contigs.fa reference

```
Now, do the mapping of the raw reads to the reference genome (this would be done with -1 and -2 if these were paired-end reads):
```
for x in SRR*_1.sub.fastq.gz;
  do bowtie2 -x reference -1 $x -2 ${x%_1*}_2.sub.fastq.gz -S ${x%_1*}.sam 2> ${x%_1*}.out;
done
```

This file contains all of the information about where each read hits on the reference.

Next, index the reference genome with samtools:

```
samtools faidx megahit_out/final.contigs.fa
```

Convert the SAM into a BAM file:
@Jin - add details on why these steps

```
for x in *.sam;
  do samtools import megahit_out/final.contigs.fa.fai $x $x.bam;
done
```

Sort the BAM file:
```
for x in *.bam;
  do samtools sort $x $x.sorted;
done
```

And index the sorted BAM file:
```
for x in *.sorted.bam;
  do samtools index $x;
done
```

## Visualizing alignment
@Jin - why do you want to visualize? Add details -- Also, is this on their computer or on the server...tablet may be tricky for them to get installed.
 
At this point you can visualize with samtools tview or Tablet.

‘samtools tview’ is a text interface that you use from the command line; run it like so:
```
samtools tview SRR098038.sorted.bam REL606.fa
```
The ‘.’s are places where the reads align perfectly in the forward direction, and the ‘,’s are places where the reads align perfectly in the reverse direction. Mismatches are indicated as A, T, C, G, etc.

You can scroll around using left and right arrows; to go to a specific coordinate, use ‘g’ and then type in the contig name and the position. For example, type ‘g’ and then ‘rel606:553093<ENTER>’ to go to position 553093 in the BAM file.

Use ‘q’ to quit.

For the Tablet viewer, click on the link and get it installed on your local computer. Then, start it up as an application. To open your alignments in Tablet, you’ll need three files on your local computer: REL606.fa, SRR098042.sorted.bam, and SRR098042.sorted.bam.bai. You can copy them over using Dropbox, for example.
## counting alignments
This command:
```
samtools view -c -f 4 SRR098038.bam
```
will count how many reads DID NOT align to the reference (207029).

This command:

```
samtools view -c -F 4 SRR098038.bam
```
will count how many reads DID align to the reference (6839602).

And this command:

```
gunzip -c SRR098038.fastq.gz | wc
```

will tell you how many lines there are in the FASTQ file (28186524). Reminder: there are four lines for each sequence.

## Make count file

```
for x in *.sorted.bam
do samtools idxstats $x > $x.idxstats.txt
done
```
then make one file
```
git clone https://github.com/metajinomics/mapping_tools.git
python mapping_tools/get_count_table.py *.idxstats.txt > counts.txt
```



If you are interested in SNP calling, visit [here](https://github.com/metajinomics/tutorials_en/blob/gh-pages/metagenome/counting-abundance-with-mapped-reads.md) for more tutorials
@this tutorial isn't that different - I wouldn't reference it....
