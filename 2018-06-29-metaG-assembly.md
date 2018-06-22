
## Digital Normalization and Metagenomic Assembly

Authored by Jin Choi, Jackson Sorenson, and Joshua Herr for EDAMAME

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

## Overarching Goal  
* This tutorial will contribute towards an understanding of **metagenome analysis**
* It focuses on metagenomic assembly

## Learning Objectives
* Use digital normalization to remove redundant data
* Trim/Filter out errors from sequences by identifying low coverage kmers in high coverage areas
* Understand the limitations and strengths of metagenome assembly
* Assemble a metagenome with MEGAHIt
* Summarize and assess the assembly
***

# Getting started
We'll start using the files we generated in the previous step (quality trimming step).  Here's where we're going to be running a program for a while (how long depends on the amount of memory your computer has and how large your data set is).  

Since this process can take a while and is prone to issues with remote computing (internet cutting out, etc.) make sure you're running in `screen` or `tmux` when connecting to your EC2 instance!


## A little background on assembly

Ok, so we've just spent a little while quality checking, quality trimming, normalizing, and (possibly, but probably not) partitioning and it's time to get some results -- we're going to assemble our reads into longer contigs and (hopefully!) entire bacterial and archaeal genomes

**Disclaimer:** Assembling metagenomes is really difficult and fraught with confounding issues.  It was only a few years ago that this was first done for a very simple community that [resulted in a paper in Science](http://www.sciencemag.org/content/335/6068/587.abstract)).  You're entering treacherous territory and there will be a lot of time spent assessing your output each step of the way, along with a ton of waiting and re-running things! Enjoy!

First, there are many, many options for assembling metagenomic data.  Most assemblers ([Velvet](http://www.ebi.ac.uk/~zerbino/velvet/), [IDBA](https://code.google.com/p/hku-idba/), [SPAdes](http://bioinf.spbau.ru/spades/)) that work well on genomic data can just as easily be used for metagenomic data, but since they were designed for use on single organisms, they might not be the best choice when you have many (to potentially thousands of) organisms which share genomic signatures.  It's difficult enough to get a good genome assembly from a pure culture of a single organism -- let alone many organisms not sequenced to the depth you need.

We've gone to the trouble of installing some assembly programs to the EDAMAME [ami](), so feel free to experiment in your free time with other assemblers.  We'll be using a *newish* assembler, Megahit v1.0.6 ([program](https://github.com/voutcn/megahit) and [paper](http://www.sciencedirect.com/science/article/pii/S1046202315301183)), which has a couple of benefits for us.  One, the assembler is optimized for (i.e. designed to handle) metagenomic reads, and two, it's pretty fast (when compared to other assemblers (i.e. SPAdes) that provide good results on metagenomic data in addition to metagenomic data). 

[De Bruijn Graph](https://www.nature.com/articles/nbt.2023), [velvet](https://genome.cshlp.org/content/18/5/821.full), 
[Explanation of effect of k-mer size](https://github.com/rrwick/Bandage/wiki/Effect-of-kmer-size)

## Running megahit

#### If your data is very big, so megahit cannot handle assembly, then you may consider remove redundant seqeunces using digital normalization
[Here is the tutorial of digital normalization](https://github.com/edamame-course/Metagenome/blob/master/2016-07-15-metaG-assembly.md)


First read the [megahit manual here](https://github.com/voutcn/megahit).  The paper can be found here: [Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.](http://www.sciencedirect.com/science/article/pii/S1046202315301183).

Install Megahit - a program that we will use to assemble reads
```
cd ~
git clone https://github.com/voutcn/megahit.git
cd megahit
make
```

You'll want to read the (minimal) manual first, but we're going to use a couple of flags.  We have to set the memory you will use in the analysis, I suggest for our case to use `-m 0.9` which means we'll use 90% of the available CPU memory.  You don't want to use 100% or your computer will not be able to run essential operations. default:0.9

#### Assembly takes long time usully, so it is a good idea to use [tmux](https://github.com/edamame-course/2015-tutorials/blob/master/final/2015-06-22_tmux.md)
Taking that into consideration, we're going to run this code:
```
cd ~/metagenome
~/megahit/megahit -1 all_trim.r1.fastq -2 all_trim.r2.fastq -r all_trim.single.fastq
```
--1: forward read, --2: reverse read of paired end, [more option](https://github.com/voutcn/megahit)
You should slowly see something similar to the following output:

...and this is going to run for a while (perhaps until a k of 91 or greater) and eventually at the end you'll see something like this:

```
--- [Fri Jun 22 18:37:24 2018] Assembling contigs from SdBG for k = 119 ---
--- [Fri Jun 22 18:37:31 2018] Merging to output final contigs ---
--- [STAT] 3670 contigs, total 2977828 bp, min 211 bp, max 37252 bp, avg 811 bp, N50 942 bp
--- [Fri Jun 22 18:37:31 2018] ALL DONE. Time elapsed: 503.442851 seconds ---
```

In the end, your assembled contigs will be in the folder called `megahit_out` and the file is final.contigs.fa.  How do you take a look at that?  And now you have a de novo assembled reference for these itty bitty metagenomes, now go try it on your own datasets.  The next step is how to you use this reference to estimate gene abundances.
