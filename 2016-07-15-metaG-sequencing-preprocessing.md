

#Metagenome quality trimming
Authored by Joshua Herr  
[EDAMAME-2015 wiki](https://github.com/edamame-course/2015-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

##Overarching Goal
* This tutorial will contribute towards an understanding of **microbial metagenome analysis**

##Learning Objectives
* Assess the quality of "raw" metagenome data
* Trim raw reads to meet quality standards

***


## Background On Quality Trimming

When DNA sequencing takes place, errors are inevitable.  No sequencing method is perfect and some are drastically different than others in regard to sequence length and quality.  We're going to trim the poor quality tail end sections of our Illumina reads.  Although in general Illumina reads are very high quality, this degradation at the end of the sequencing run is typical of the Illumina sequencing platforms.

Some sequencing centers will remove library adapters (our sequencing center does), but you will have to check with your data provider to know what they give you and ALWAYS check for your self to verify what you have been told.

As always, you want to make sure you read the manual of any tool to be sure you know what the tool is doing to your data and if that is the right tool for the job.  Knowing which tool to use is very important -- you wouldn't use a saw to put a nail in a piece of wood, would you?

We'll be using a tool which is not aware of paired-end reads. This is fine as the downstream metagenome assembly program we will use (megahit) only takes single end reads.  If you choose to use a different assembly program that accepts (or even requires!) paired-end reads then you will have to choose a tool to do trimming on the paired-ends.  

## Quality Trimming Your Sequence Data

1.  Start a ```m3.large``` machine from Amazon Web Services running the EDAMAME-2015 AWS ami (ami-af04f2c4).  This instance has about 8 GB of RAM, and 2 CPUs, and should be enough to complete the assembly of the example data set we will use.

**Note:** One of the issues with processing whole genome shotgun data is how long it takes for the computer to process many steps of the workflow.  This can be time consuming and you should consider using ```screen``` or ```tmux``` to ensure that an internet connection issue does not cause you to lose your workflow progress.

**Pro-Tip:** You'll also want to keep in mind that these assemblies take a lot of computer power to run which can cost you some money -- for your own benefit, you can try to optimize your scripts on a desktop or laptop before you actually fire up the AWS instance of this size.

Download the data:
```
wget https://s3.amazonaws.com/edamame/EDAMAME_MG.tar.gz
tar -zxvf EDAMAME_MG.tar.gz
```

2.  First, let's get an idea of some quality stats from our data.  We're going to first use the ```fastx_quality_stats``` [script](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_statistics_usage) from the Hannon Lab's [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html) package.

```
fastx_quality_stats -i <filename>.fastq -o quality.txt

cat quality.txt
```

This will give us some idea of what we are dealing with.  We'll want to keep this in mind when we check the quality after trimming.

Then we run this command:

```
fastq_quality_filter -i <filename>.fastq -Q33 -q 30 -p 50  -o <filename>.qc.fastq
```

This command first uses the ```fastq_quality_filter``` [script](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage) from Hannon Lab's [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html) to trim the data using Illumina-33 [Phred quality score](http://en.wikipedia.org/wiki/Phred_quality_score). 

Note that you can modify the ```fastq_quality_filter``` [script](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage) to trim to any specific length or quality level that you desire.  As always, read the [manual](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage) for information on how to use a script.

If when you are using your own data, the ```fastq_quality_filter``` complains about invalid quality scores, first try removing the -Q33 in the command.  There are numerous types of quality scores and you may have older data which did not use the Q33 output.  For more information on fastq quality scores, [this is a good overview](http://en.wikipedia.org/wiki/FASTQ_format).

For a sanity check, let's use the ```fastx_quality_stats``` script again to see what changed in our trimmed data files:

```
fastx_quality_stats -i <filename>.qc.fastq -o qc_quality.txt

cat quality.txt

cat qc_quality.txt
```

What are the differences between the raw data and the quality trimmed data?


## Other tools for quality trimming

There are other tools for quality trimming your data -- some are for specific types of data and have different features.  It's a good idea to read the manual of each tool and do a sanity check on your data after using the tools (at the very least ```head```, ```tail```, ```more```, ```less```, *etc*.) to see that you are actually removing what you think you are removing.

Some handy quality and/or adapter trimming tools you might want to investigate are:   
   1. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) - all purpose
   2. [cutadapt](https://code.google.com/p/cutadapt/) - adapter trimming
   3. [sickle](https://github.com/najoshi/sickle) - read quality trimming
   4. [scythe](https://github.com/vsbuffalo/scythe) - adapter contamination trimming


## Next step

Now we're going to take our quality trimmed reads and run digital normalization on the reads to remove redundant information and also remove some erroneous information.

## Other resources
   * [ANGUS documentation](http://angus.readthedocs.org/en/2014/short-read-quality-evaluation.html)


##Digital Normalization 

Authored by Joshua Herr with contribution from Jackson Sorensen for EDAMAME2015     

[EDAMAME-2015 wiki](https://github.com/edamame-course/2015-tutorials/wiki)

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

##Overarching Goal  
* This tutorial will contribute towards an understanding of **microbial metagenome analysis**

##Learning Objectives
* Use digital normalization to remove redundant data
* Trim/Filter out errors from sequences by identifying low coverage kmers in high coverage areas

***

# Getting started
We'll start using the files we generated in the previous step (quality trimming step).  Here's where we're going to be running a program for a while (how long depends on the amount of memory your computer has and how large your data set is).  

Since this process can take a while and is prone to issues with remote computing (internet cutting out, etc.) make sure you're running in `screen` or `tmux` when connecting to your EC2 instance!

# Run a First Round of Digital Normalization
Normalize everything to a coverage of 20. The normalize-by-media.py script keeps track of the number of times a particular kmer is present. The flag `-C` sets a median kmer coverage cutoff for sequence. In otherwords, if the median coverage of the kmers in a particular sequence is above this cutoff then the sequence is discarded, if it is below this cutoff then it is kept. We specify the length of kmer we want to analyze using the `-k` flag. The flags `-N` and `-x` work together to specify the amount of memory to be used by the program. As a rule of thumb, the two multiplied should be equal to the available memory(RAM) on your machine. You can check the available memory on your machine with `free -m`. For our m3.large instances we should typically have about 4GB of RAM free.    

```
python /usr/local/share/khmer/scripts/normalize-by-median.py -k 20 -C 20 -N 4 -x 1e9 -s normC20k20.kh *qc.fastq
```

Make sure you read the manual for this script, it's part of the [khmer](https://github.com/ged-lab/khmer) package.  This script produces a set of '.keep' files, as well as a normC20k20.kh database file.  The database file (it's a hash table in this case) can get quite large so keep in ming when you are running this script on a lot of data with not a lot of free space on your computer.

# Removing Errors from our data
We'll use the `filter-abund.py` script to trim off any k-mers that are in abundance of 1 in high-coverage reads.

```
python /usr/local/share/khmer/scripts/filter-abund.py -V normC20k20.kh *.keep
```

The output from this step produces files ending in `.abundfilt` that contain the trimmed sequences.

If you read the manual, you see that the `-V` option is used to make this work better for variable coverage data sets, such as those you would find in metagenomic sequencing.  If you're using this tool for a genome sequencing project, you wouldn't use the `-V` flag.

# Normalize down to a coverage of five
Now that we've eliminated many more erroneous k-mers from the dataset, let's ditch some more high-coverage data. Normalize down to a coverage of five using the following command. Note that here we are loading the count table from the first round of digital normalization and normalizing our error filtered data from the filter-abund.py step. 

```
python /usr/local/share/khmer/scripts/normalize-by-median.py -x 1e09 -C 5 -s normC5k20.kh -l normC20k20.kh *qc.fastq.keep.abundfilt
```

Now, we'll have a file (or list of files if you're using your own data) which will have the name: `{your-file}.qc.fastq.keep.abundfilt.keep`.  We're going to check the file integrity to make sure it's not faulty and we're going to clean up the names.

Let's rename your files:
```
mv {your-file}.qc.fastq.keep.abundfilt.keep {your-file}_single.fastq  
```
These files will be used in the next section where we assemble your metagenomic reads.

##Help and other Resources
* [khmer documentation and repo](https://github.com/dib-lab/khmer/blob/master/README.rst)
* [khmer protocols - see "Kalamazoo Protocol"](http://khmer-protocols.readthedocs.org/en/v0.8.4/)
* [khmer recipes - bite-sized tasks using khmer scripts](http://khmer-recipes.readthedocs.org/en/latest/)
* [khmer discussion group](http://lists.idyll.org/listinfo/khmer)