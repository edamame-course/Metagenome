
##Digital Normalization and Metagenomic Assembly

Authored by Jin Choi, Jackson Sorenson, and Joshua Herr for EDAMAME2016     

[EDAMAME-2016 wiki](https://github.com/edamame-course/2016-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

##Overarching Goal  
* This tutorial will contribute towards an understanding of **metagenome analysis**
* It focuses on metagenomic assembly

##Learning Objectives
* Use digital normalization to remove redundant data
* Trim/Filter out errors from sequences by identifying low coverage kmers in high coverage areas
* Understand the limitations and strengths of metagenome assembly
* Assemble a metagenome with MEGAHIt
* Summarize and assess the assembly
***

# Getting started
We'll start using the files we generated in the previous step (quality trimming step).  Here's where we're going to be running a program for a while (how long depends on the amount of memory your computer has and how large your data set is).  

Since this process can take a while and is prone to issues with remote computing (internet cutting out, etc.) make sure you're running in `screen` or `tmux` when connecting to your EC2 instance!

# Run a Digital Normalization
We have found, based on mock communities, that the coverage you need to adequately assemble a 10-20 strain community is roughly 20x coverage.  However, in a metagenome, we can have varying coverages of genes, some much more than 20x.  Our goal will be to digitaly normalize our sequencing reads to a coverage of 20, basically to remove redundant reads.

The normalize-by-media.py script keeps track of the number of times a particular kmer is present and allows us to estimate the observed coverage of each read. The flag `-C` sets a median kmer coverage cutoff for sequence. In otherwords, if the median coverage of the kmers in a particular sequence is above this cutoff then the sequence is discarded (or set aside); if it is below this cutoff then it is kept. We specify the length of kmer we want to analyze using the `-k` flag. The flags `-N` and `-x` work together to specify the amount of memory to be used by the program. As a rule of thumb, the two multiplied should be equal to the available memory(RAM) on your machine. You can check the available memory on your machine with `free -m`. For our m3.large instances we should typically have about 4GB of RAM free.    

First, let's think about what happened after we did the quality filter step.  Its similar to the adapter-trimming step.  We'll get some orphans and we need to discriminate our paired and single ended reads for assembly.  We can do this with the following script which extracts paired reads.

```
extract-paired-reads.py SRR492065.combined.qc.fq
extract-paired-reads.py SRR492066.combined.qc.fq
```

Now, we're going to run the "normalization" of our datasets to a coverage of 20, reducing the dataset size for assembly and removing extra information that could contain sequencing errors.

```
cd ~/metagenome
normalize-by-median.py --ksize 20 -R diginorm.report -C 20 --n_tables 4 --max-tablesize 1e8 -p -s normC20k20.kh SRR49206?.combined.qc.fq.pe
```
-k: k-mer size, -C: k-mer coverage level above is above this numer the read is not kept, -N: number of thread(?) or number of mailroom(?) (Actually number of table), -x: memory use, -s: save the k-mer countgraph to disk, -p: paired-end, last argument: file name

Make sure you read the manual for this script, it's part of the [khmer](https://github.com/ged-lab/khmer) package.  This script produces a set of '.keep' files, as well as a normC20k20.kh database file.  The database file (it's a hash table in this case) can get quite large so keep in ming when you are running this script on a lot of data with not a lot of free space on your computer.

# Removing Errors from our data

Just like in amplicon sequencing where you sometimes may want to remove singletons, we can also remove "low coverage" reads prior to assembly to improve our data quality as well as reduce the dataset size we need to assemble.

If you read the manual, you see that the `-V` option is used to make this work better for variable coverage data sets, such as those you would find in metagenomic sequencing.  If you're using this tool for a genome sequencing project, you wouldn't use the `-V` flag.

```
filter-abund.py -V normC20k20.kh *.keep
```

The output from this step produces files ending in `.abundfilt` that contain the trimmed sequences.

This produces .abundfilt files containing the trimmed sequences.

As for most processing, now you still have data that yet again needs to process.  When you trimmed low abundant reads, you now have more orphans that you have to deal with -- and again you'll want to extract paired end reads.  The reason for doing this is because the programs that we use require it to be so and unless we want to write our own programs, we have to learn to get our data into the shape for a specific tool.  Its nice to know a little bit about for loops now:

```
for i in *.keep.abundfilt
do
   extract-paired-reads.py $i
done
```

So, now we have a bunch of files, with similar names.  They are already pretty big - so to save size, we can compress these and combine multiple samples into one file.  Let's do this for all paired end files.

Merge single end
```
cat SRR492065.single.qc.fq SRR492066.single.qc.fq SRR492065.combined.qc.fq.se SRR492066.combined.qc.fq.se > all.single.fq
```

```
gzip *abundfilt.pe
cat *abundfilt.pe.gz > abundfilt-all.gz
```

And we did all this....so now we can FINALLY assemble.

## A little background on assembly

Ok, so we've just spent a little while quality checking, quality trimming, normalizing, and (possibly, but probably not) partitioning and it's time to get some results -- we're going to assemble our reads into longer contigs and (hopefully!) entire bacterial and archaeal genomes

**Disclaimer:** Assembling metagenomes is really difficult and fraught with confounding issues.  It was only a few years ago that this was first done for a very simple community that [resulted in a paper in Science](http://www.sciencemag.org/content/335/6068/587.abstract)).  You're entering treacherous territory and there will be a lot of time spent assessing your output each step of the way, along with a ton of waiting and re-running things! Enjoy!

First, there are many, many options for assembling metagenomic data.  Most assemblers ([Velvet](http://www.ebi.ac.uk/~zerbino/velvet/), [IDBA](https://code.google.com/p/hku-idba/), [SPAdes](http://bioinf.spbau.ru/spades/)) that work well on genomic data can just as easily be used for metagenomic data, but since they were designed for use on single organisms, they might not be the best choice when you have many (to potentially thousands of) organisms which share genomic signatures.  It's difficult enough to get a good genome assembly from a pure culture of a single organism -- let alone many organisms not sequenced to the depth you need.

We've gone to the trouble of installing some assembly programs to the EDAMAME [ami](), so feel free to experiment in your free time with other assemblers.  We'll be using a *newish* assembler, Megahit v1.0.6 ([program](https://github.com/voutcn/megahit) and [paper](http://www.sciencedirect.com/science/article/pii/S1046202315301183)), which has a couple of benefits for us.  One, the assembler is optimized for (i.e. designed to handle) metagenomic reads, and two, it's pretty fast (when compared to other assemblers (i.e. SPAdes) that provide good results on metagenomic data in addition to metagenomic data). 


## Running megahit

First read the [megahit manual here](https://github.com/voutcn/megahit).  The paper can be found here: [Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.](http://www.sciencedirect.com/science/article/pii/S1046202315301183).

Install Megahit - a program that we will use to assemble reads
```
cd ~
git clone https://github.com/voutcn/megahit.git
cd megahit
make
```

You'll want to read the (minimal) manual first, but we're going to use a couple of flags.  We have to set the memory you will use in the analysis, I suggest for our case to use `-m 0.9` which means we'll use 90% of the available CPU memory.  You don't want to use 100% or your computer will not be able to run essential operations. default:0.9


Taking that into consideration, we're going to run this code:
```
cd ~/metagenome
~/megahit/megahit --12 abundfilt-all.gz -r all.single.fq
```
--12: paired end, [more option](https://github.com/voutcn/megahit)
You should slowly see something similar to the following output:

...and this is going to run for a while (perhaps until a k of 91 or greater) and eventually at the end you'll see something like this:

```
--- [Mon Jul 11 15:17:57 2016] Assembling contigs from SdBG for k = 99 ---
--- [Mon Jul 11 15:18:02 2016] Merging to output final contigs ---
--- [STAT] 2764 contigs, total 1377587 bp, min 216 bp, max 15822 bp, avg 498 bp, N50 491 bp
--- [Mon Jul 11 15:18:02 2016] ALL DONE. Time elapsed: 185.048481 seconds ---
```

In the end, your assembled contigs will be in the folder called `megahit_out` and the file is final.contigs.fa.  How do you take a look at that?  And now you have a de novo assembled reference for these itty bitty metagenomes, now go try it on your own datasets.  The next step is how to you use this reference to estimate gene abundances.
