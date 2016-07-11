
##Digital Normalization and Assembly

Authored by Joshua Herr with contribution from Jackson Sorensen and Jin Choi for EDAMAME2016     

[EDAMAME-2016 wiki](https://github.com/edamame-course/2016-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

##Overarching Goal  
* This tutorial will contribute towards an understanding of **metagenome analysis**


##Learning Objectives
* Use digital normalization to remove redundant data
* Trim/Filter out errors from sequences by identifying low coverage kmers in high coverage areas
* Understand the limitations and strengths of metagenome assembly
* Assemble a metagenome with MEGAHIt, cite alternative assemblers available
* Summarize and assess the assembly
***

# Getting started
We'll start using the files we generated in the previous step (quality trimming step).  Here's where we're going to be running a program for a while (how long depends on the amount of memory your computer has and how large your data set is).  

Since this process can take a while and is prone to issues with remote computing (internet cutting out, etc.) make sure you're running in `screen` or `tmux` when connecting to your EC2 instance!

# Run a Digital Normalization
Normalize everything to a coverage of 20. The normalize-by-media.py script keeps track of the number of times a particular kmer is present. The flag `-C` sets a median kmer coverage cutoff for sequence. In otherwords, if the median coverage of the kmers in a particular sequence is above this cutoff then the sequence is discarded, if it is below this cutoff then it is kept. We specify the length of kmer we want to analyze using the `-k` flag. The flags `-N` and `-x` work together to specify the amount of memory to be used by the program. As a rule of thumb, the two multiplied should be equal to the available memory(RAM) on your machine. You can check the available memory on your machine with `free -m`. For our m3.large instances we should typically have about 4GB of RAM free.    
(20 min)
```
cd ~/metagenome
normalize-by-median.py -k 20 -C 20 -N 4 -x 1e9 -s normC20k20.kh *qc.fq
```

Make sure you read the manual for this script, it's part of the [khmer](https://github.com/ged-lab/khmer) package.  This script produces a set of '.keep' files, as well as a normC20k20.kh database file.  The database file (it's a hash table in this case) can get quite large so keep in ming when you are running this script on a lot of data with not a lot of free space on your computer.

# Removing Errors from our data
We'll use the `filter-abund.py` script to trim off any k-mers that are in abundance of 1 in high-coverage reads.
(15 min)
```
filter-abund.py -V normC20k20.kh *.keep
```

The output from this step produces files ending in `.abundfilt` that contain the trimmed sequences.

If you read the manual, you see that the `-V` option is used to make this work better for variable coverage data sets, such as those you would find in metagenomic sequencing.  If you're using this tool for a genome sequencing project, you wouldn't use the `-V` flag.

This produces .abundfilt files containing the trimmed sequences.

The process of error trimming could have orphaned reads, so split the PE file into still-interleaved and non-interleaved reads:
(5 min)
```
for i in *.keep.abundfilt
do
   extract-paired-reads.py $i
done
```

Now, we'll have a file (or list of files if you're using your own data) which will have the name: `{your-file}.qc.fq.keep.abundfilt`.  We're going to check the file integrity to make sure it's not faulty and we're going to clean up the names.

Let's compress and rename your files:

```
gzip *abundfilt.pe
cat *abundfilt.pe.gz > abundfilt-all.gz
```


## A little background on assembly

Ok, so we've just spent a little while quality checking, quality trimming, normalizing, and (possibly, but probably not) partitioning and it's time to get some results -- we're going to assemble our reads into longer contigs and (hopefully!) entire bacterial and archaeal genomes

**Disclaimer:** Assembling metagenomes is really difficult and fraught with confounding issues.  It was only a few years ago that this was first done for a very simple community that [resulted in a paper in Science](http://www.sciencemag.org/content/335/6068/587.abstract)).  You're entering treacherous territory and there will be a lot of time spent assessing your output each step of the way, along with a ton of waiting and re-running things! Enjoy!

First, there are many, many options for assembling metagenomic data.  Most assemblers ([Velvet](http://www.ebi.ac.uk/~zerbino/velvet/), [IDBA](https://code.google.com/p/hku-idba/), [SPAdes](http://bioinf.spbau.ru/spades/)) that work well on genomic data can just as easily be used for metagenomic data, but since they were designed for use on single organisms, they might not be the best choice when you have many (to potentially thousands of) organisms which share genomic signatures.  It's difficult enough to get a good genome assembly from a pure culture of a single organism -- let alone many organisms not sequenced to the depth you need.

We've gone to the trouble of installing some assembly programs to the EDAMAME [ami](), so feel free to experiment in your free time with other assemblers.  We'll be using a *newish* assembler, Megahit v1.0.6 ([program](https://github.com/voutcn/megahit) and [paper](http://www.sciencedirect.com/science/article/pii/S1046202315301183)), which has a couple of benefits for us.  One, the assembler is optimized for (i.e. designed to handle) metagenomic reads, and two, it's pretty fast (when compared to other assemblers (i.e. SPAdes) that provide good results on metagenomic data in addition to metagenomic data). 


## Running megahit

First read the [megahit manual here](https://github.com/voutcn/megahit).  The paper can be found here: [Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.](http://www.sciencedirect.com/science/article/pii/S1046202315301183).

You'll want to read the (minimal) manual first, but we're going to use a couple of flags:
  1. We have to set the memory you will use in the analysis, I suggest for our case to use `-m 0.9` which means we'll use 90% of the available CPU memory.  You don't want to use 100% or your computer will not be able to run essential operations. defalut:0.9
  2. Megahit requires us to set the length of the reads that will be ignored.  Just to be safe I have used `-l 500` here, but change it and see if it changes your assembly.  I would not go below your average read length.

Taking that into consideration, we're going to run this code:
```
~/megahit/megahit --12 abundfilt-all.gz
```

You should slowly see something similar to the following output:

```
7.0Gb memory in total.
Using: 7.017Gb.
MEGAHIT v1.0.6-3-gfb1e59b
--- [Mon Jul 11 15:14:57 2016] Start assembly. Number of CPU threads 2 ---
--- [Mon Jul 11 15:14:57 2016] Available memory: 8371847168, used: 7534662451
--- [Mon Jul 11 15:14:57 2016] k list: 21,41,61,81,99 ---
--- [Mon Jul 11 15:14:57 2016] Converting reads to binaries ---
    [read_lib_functions-inl.h  : 209]     Lib 0 (abundfilt-all.gz): interleaved, 156468 reads, 100 max length
    [utils.h                   : 126]     Real: 0.3885	user: 0.0872	sys: 0.1221	maxrss: 8204
--- [Mon Jul 11 15:14:57 2016] Extracting solid (k+1)-mers for k = 21 ---
--- [Mon Jul 11 15:15:03 2016] Building graph for k = 21 ---
--- [Mon Jul 11 15:15:08 2016] Assembling contigs from SdBG for k = 21 ---
--- [Mon Jul 11 15:15:16 2016] Local assembling k = 21 ---
--- [Mon Jul 11 15:15:37 2016] Extracting iterative edges from k = 21 to 41 ---
--- [Mon Jul 11 15:15:39 2016] Building graph for k = 41 ---
--- [Mon Jul 11 15:15:42 2016] Assembling contigs from SdBG for k = 41 ---
--- [Mon Jul 11 15:15:50 2016] Local assembling k = 41 ---
--- [Mon Jul 11 15:16:21 2016] Extracting iterative edges from k = 41 to 61 ---
--- [Mon Jul 11 15:16:22 2016] Building graph for k = 61 ---
--- [Mon Jul 11 15:16:24 2016] Assembling contigs from SdBG for k = 61 ---
--- [Mon Jul 11 15:16:30 2016] Local assembling k = 61 ---

```

...and this is going to run for a while (perhaps until a k of 91 or greater) and eventually at the end you'll see something like this:

```
--- [Mon Jul 11 15:17:57 2016] Assembling contigs from SdBG for k = 99 ---
--- [Mon Jul 11 15:18:02 2016] Merging to output final contigs ---
--- [STAT] 2764 contigs, total 1377587 bp, min 216 bp, max 15822 bp, avg 498 bp, N50 491 bp
--- [Mon Jul 11 15:18:02 2016] ALL DONE. Time elapsed: 185.048481 seconds ---
```

