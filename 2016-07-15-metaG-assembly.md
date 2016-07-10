#Metagenome Assembly
***
Authored by Joshua Herr with contribution from Jackson Sorensen for EDAMAME2015     
[EDAMAME-2015 wiki](https://github.com/edamame-course/2015-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

##Overarching Goal  
* This tutorial will contribute towards an understanding of **metagenome analysis**

##Learning Objectives
* Understand the limitations and strengths of metagenome assembly
* Assemble a metagenome with MEGAHIt, cite alternative assemblers available
* Summarize and assess the assembly
***

## A little background on assembly

Ok, so we've just spent a little while quality checking, quality trimming, normalizing, and (possibly, but probably not) partitioning and it's time to get some results -- we're going to assemble our reads into longer contigs and (hopefully!) entire bacterial and archaeal genomes

**Disclaimer:** Assembling metagenomes is really difficult and fraught with confounding issues.  It was only a few years ago that this was first done for a very simple community that [resulted in a paper in Science](http://www.sciencemag.org/content/335/6068/587.abstract)).  You're entering treacherous territory and there will be a lot of time spent assessing your output each step of the way, along with a ton of waiting and re-running things! Enjoy!

First, there are many, many options for assembling metagenomic data.  Most assemblers ([Velvet](http://www.ebi.ac.uk/~zerbino/velvet/), [IDBA](https://code.google.com/p/hku-idba/), [SPAdes](http://bioinf.spbau.ru/spades/)) that work well on genomic data can just as easily be used for metagenomic data, but since they were designed for use on single organisms, they might not be the best choice when you have many (to potentially thousands of) organisms which share genomic signatures.  It's difficult enough to get a good genome assembly from a pure culture of a single organism -- let alone many organisms not sequenced to the depth you need.

We've gone to the trouble of installing some assembly programs to the EDAMAME [ami](), so feel free to experiment in your free time with other assemblers.  We'll be using a *newish* assembler, Megahit v0.2.1 ([program](https://github.com/voutcn/megahit/tree/v0.2.1-a) and [paper](http://bioinformatics.oxfordjournals.org/content/31/10/1674.long)), which has a couple of benefits for us.  One, the assembler is optimized for (i.e. designed to handle) metagenomic reads, and two, it's pretty fast (when compared to other assemblers (i.e. SPAdes) that provide good results on metagenomic data in addition to metagenomic data). Just to note Megahit is now up to version [v0.3-b2](https://github.com/voutcn/megahit), but we will not be using this version in the class. In v0.3-b2 Megahit can use PE information. 

## Preparing the Data

Each assembler program takes data input in slightly different formats.  Velvet is the most flexible: it can take in paired-end and/or single reads, in FASTA or FASTQ format, gzipped or not.  On the other hand, SPADes requires paired-end FASTQ, while IDBA requires paired-end FASTA.  Make sure you always spend time reading the manual of a program first.

Since we're going to use megahit, which takes single-end reads (in this case we have disregarded the paired-end metagenomic data and combined the paired reads into a single file -- yes, this throws out some of the information you get from having paired reads), we want our data in a single file as the input.  We've been running our data as single end from the start, but you should know that if you have paired end data and you want to utilize the information from the pairs you'll have to plan accordingly.

## Running megahit

First read the [megahit manual here](https://github.com/voutcn/megahit).  The paper can be found here: [Li et al. 2015 MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph](http://bioinformatics.oxfordjournals.org/content/31/10/1674.abstract).

You'll want to read the (minimal) manual first, but we're going to use a couple of flags:
  1. We have to set the memory you will use in the analysis, I suggest for our case to use `-m 0.9` which means we'll use 90% of the available CPU memory.  You don't want to use 100% or your computer will not be able to run essential operations.
  2. Megahit requires us to set the length of the reads that will be ignored.  Just to be safe I have used `-l 500` here, but change it and see if it changes your assembly.  I would not go below your average read length.

Taking that into consideration, we're going to run this code:

```
megahit -m 0.9 -l 500 --cpu-only -r {your-file}_single.fastq -o megahit_assembly
```

You should slowly see something similar to the following output:

```
23.0Gb memory in total.
Using: 21.112Gb.
MEGAHIT v0.2.1-beta
[Thu Jun  25 09:53:10 2015] Start assembly. Number of CPU threads 8.
[Thu Jun  25 09:53:10 2015] Extracting solid (k+1)-mers for k = 21
[Thu Jun  25 09:59:09 2015] Building graph for k = 21
[Thu Jun  25 10:00:49 2015] Assembling contigs from SdBG for k = 21
[Thu Jun  25 10:05:00 2015] Extracting iterative edges from k = 21 to 31
[Thu Jun  25 10:07:43 2015] Building graph for k = 31
[Thu Jun  25 10:09:39 2015] Assembling contigs from SdBG for k = 31
[Thu Jun  25 10:12:27 2015] Extracting iterative edges from k = 31 to 41
[Thu Jun  25 10:13:34 2015] Building graph for k = 41
[Thu Jun  25 10:15:35 2015] Assembling contigs from SdBG for k = 41
[Thu Jun  25 10:18:06 2015] Extracting iterative edges from k = 41 to 51
[Thu Jun  25 10:19:09 2015] Building graph for k = 51
[Thu Jun  25 10:20:53 2015] Assembling contigs from SdBG for k = 51
[Thu Jun  25 10:23:02 2015] Extracting iterative edges from k = 51 to 61
```

...and this is going to run for a while (perhaps until a k of 91 or greater) and eventually at the end you'll see something like this:

```
[Thu Jun  25 18:43:16 2015] Merging to output final contigs.
[Thu Jun  25 18:43:18 2015] ALL DONE.
```

## Assembly Assessment

So, now we want to take a look at our assembly and see how it "looks".  We'll navigate into the folder we named (`megahit_assembly`) as the output by typing `cd megahit-assembly` and then `ls -lah` the contents of the folder.  You should see something like this:

```
drwxr-xr-x 3 user group    7 Jun  4 10:43 .
drwxr-xr-x 9 user group  101 Jun  4 09:53 ..
-rw-r--r-- 1 user group 147M Jun  4 10:43 final.contigs.fa
-rw-r--r-- 1 user group 1.2M Jun  4 10:43 log
-rw-r--r-- 1 user group   71 Jun  4 09:53 opts.txt
drwxr-xr-x 2 user group   13 Jun  4 10:43 tmp
```

We're really concerned with two files, the `log` of the run and, of course, the assembly `final.contigs.fa`

Let's take a look at the `log` file first; go ahead and type `tail log`.  This will give you some of megahit's stats on the assembly.  The whole log file will give you a more in depth account of what happened than was printed to the screen.

We're interested in this part:

```
Simple path graph size: 81489
Number of unitigs removed: 113
Total length: 86752217, N50: 5464, Mean: 1482, number of contigs: 58510
Maximum length: 234281
```

Megahit actually runs numerous iterations of assemblies. Here we are looking at the last iteration and the number of contigs and total length of that iteration of the assembly. Later, we will be using QUAST to calculate all of our assembly stats.  There may be some odd terminology in some of the log file, [see this wiki for more information on terminology arising from the Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology) (the first assembly program designed for the Human Genome Project)

## Calculating summary statistics for our Metagenome assembly, using QUAST. 

Quast is a program that will calculate some statistics about our metagenome assembly to give us an idea how well it assembled. Again, we do need to quickly install Quast and one of its dependencies so we can get it running on our machine. This shouldn't take very long. Copy each of the following lines of code one line at a time to install Quast.

```
sudo apt-get install python-matplotlib
cd ~
wget https://downloads.sourceforge.net/project/quast/quast-2.3.tar.gz
tar xzvf quast-2.3.tar.gz
cd quast-2.3
```

Let's run a quick test to ensure quast is working. 
```
python metaquast.py --test
```

The test should run fairly quickly and not shoot any warning or nonfatal errors. 
Now we can take a look at our assembly using QUAST. **From the ~/quast-2.3 directory** run the following line of code. 
```
python metaquast.py -o ~/EDAMAME_MG/megahit_assembly ~/EDAMAME_MG/megahit_assembly/final.contigs.fa 
```
Once QUAST has finished running, change into the quast_output directory and use `ls` to take a look at all of the files it created. Use `less` to examine the `report.txt` file. 
```
Assembly                   final.contigs
# contigs (>= 0 bp)        33174        
# contigs (>= 1000 bp)     1842         
Total length (>= 0 bp)     22158177     
Total length (>= 1000 bp)  2184708      
# contigs                  33174        
Largest contig             2620         
Total length               22158177     
GC (%)                     58.93        
N50                        647          
N75                        560          
L50                        13442        
L75                        22678        
# N's per 100 kbp          0.00
```
  