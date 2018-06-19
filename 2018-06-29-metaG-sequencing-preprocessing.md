

# Metagenome quality trimming
Authored by Jin Choi for EDAMAME
(with contributions from Josh Herr from EDAMAME2015)  

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

## Overarching Goal
* This tutorial will contribute towards an understanding of **microbial metagenome analysis**

## Learning Objectives
* Assess the quality of "raw" metagenome data
* Trim raw reads to meet quality standards

***


## Background On Quality Trimming

When DNA sequencing takes place, errors are inevitable.  No sequencing method is perfect and some are drastically different than others in regard to sequence length and quality.  We're going to trim the poor quality tail end sections of our Illumina reads.  Although in general Illumina reads are very high quality, this degradation at the end of the sequencing run is typical of the Illumina sequencing platforms.

Some sequencing centers will remove library adapters (our sequencing center does), but you will have to check with your data provider to know what they give you and ALWAYS check for your self to verify what you have been told.

As always, you want to make sure you read the manual of any tool to be sure you know what the tool is doing to your data and if that is the right tool for the job.  Knowing which tool to use is very important -- you wouldn't use a saw to put a nail in a piece of wood, would you?

We'll be using a tool which is aware of paired-end reads and single-end together (megahit).


## Quality Trimming Your Sequence Data

The first step to starting an assembly or metagenomic read analysis is to remove bad quality sequences, a process we call quality trimming.  We'll start a cloud compute instance and practice trimming bad quality reads.

1.  Start a ```m4.large``` machine from Amazon Web Services.  This instance has about 8 GB of RAM, and 2 CPUs, and should be enough to complete the assembly of the example data set we will use. However, this may NOT enough for your FULL-SIZE DATA. You may consider 125 GB of RAM, and 8 CPUs or bigger. [Here how to connect EC2](http://angus.readthedocs.io/en/2015/amazon/index.html)

**Note:** One of the issues with processing whole genome shotgun data is how long it takes for the computer to process many steps of the workflow.  This can be time consuming and you should consider using ```screen``` or ```tmux``` to ensure that an internet connection issue does not cause you to lose your workflow progress.

**Pro-Tip:** You'll also want to keep in mind that these assemblies take a lot of computer power to run which can cost you some money -- for your own benefit, you can try to optimize your scripts on a desktop or laptop before you actually fire up the AWS instance of this size.

Install dependencies (required software)
```
sudo apt-get update
sudo apt-get -y install python-dev python-pip fastx-toolkit unzip git zlib1g-dev default-jre
```

Install Trimmomatic - a program to quality trim reads
```
cd 
curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
```

Download the data (You did on last session): The tutorial data is from [Sharon et al. 2013](http://www.ncbi.nlm.nih.gov/pubmed/22936250); it’s two data points from an infant gut sample. And it is a subsampled file so that we can do this tutorial in some reasonable amount of time, 100,000 sequences. Go to the home directory, make a directory named `metagenome`, change directores into the folder, download data, then unzip the data or decompress the data file. 
```
cd
mkdir metagenome
cd metagenome
wget https://s3.amazonaws.com/edamame/infant_gut.sub.tar.gz
tar -zxvf infant_gut.sub.tar.gz
```
This dataset contains paired end reads.  The paired end nature of datasets can be a huge advantage for assembly -- we know the approximate insert size between these reads.  However, when we trim bad quality sequences and/or reads, sometimes a previously paired read is left as an orphan (sad face).  However, for downstream analyses, we want to keep paired end reads that are together and separate these from single ended reads.  This is a bit of a process but important to understand.    

What's the quality of the dataset?


1.  First, let's get an idea of some quality stats from our data.  We're going to first use the ```fastx_quality_stats``` [script](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_statistics_usage) from the Hannon Lab's [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html) package.

```
gunzip -c SRR492065_1.sub.fastq.gz > SRR492065_1.unzip_for_quality.fastq
head SRR492065_1.unzip_for_quality.fastq
```

It looks like 
```
@SRR492065.1 HWI-EAS385_0095_FC:2:1:6702:1434/1
TCAGCCATCGCTATGCTTGGCTTCACTGTGAAGACCACTCCAATCGCGACTTGTCACGATTGTCGTTACCATTAAANNNNNGAAAACAGGAGAACAAGTA
+
IIIIIIHIIFIIIIIIIIHIIGIIIIIIEIIIGIIIIIDGGDEEFHIIEIEHEHGEHE>EEBDFEEEECF8AA??4#####5789:/1,8B?@B<@B###
@SRR492065.2 HWI-EAS385_0095_FC:2:1:6931:1435/1
GCAATAGCAGGCTCACCGACTGTGATTTTACTAGAATTTCCAAACTCAGCGACTTGATCAACTTTGTCGGAAGAAANNNNNATCACGGCTAAATCCGTAT
+
DHHFBHDDBBGHFHHHHGHHHHHHHHHHHHHHHEHHHHHHHHHHHHHHHHHBHHGHFHH>HHHHH@HFEGCHDEH1#####=;;;??3=9B8BB=CEB@C
@SRR492065.3 HWI-EAS385_0095_FC:2:1:9984:1431/1
AAGCTCATTTTTATCATAAGCTACAAAGAACGGATAGTAAATCAAGATGGCCACTAGAATTAAAACAATATTCAAANNNNNAGCACGCCAATCGCCGCCT
```

This will give us some idea of what we are dealing with.  We'll want to keep this in mind when we check the quality after trimming.


 There are numerous types of quality scores.  For more information on fastq quality scores, [this is a good overview](http://en.wikipedia.org/wiki/FASTQ_format).
 ```
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40                                
                           -5....0........9.............................40 
                                 0........9.............................40 
                                    3.....9.............................40 
  0.2......................26...31........41                             

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
     (Note: See discussion above).
L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 ```

For a sanity check, let's use the ```fastx_quality_stats``` script again to see what changed in our trimmed data files:

What are the differences between the raw data and the quality trimmed data?
```
fastx_quality_stats -i SRR492065_1.unzip_for_quality.fastq -o before_trim.quality.txt
cat before_trim.quality.txt
```

The innerworkings:  -i: input file, -o: output file

First, let's trim the sequencing adapters from our reads (these are artificial non-biological sequences and we want them out) and then we'll place all the paired end reads into one combined file. 

```
for x in SRR*_1.sub.fastq.gz;
do java -jar ../Trimmomatic-0.36/trimmomatic-0.36.jar PE $x ${x%_1*}_2.sub.fastq.gz ${x%_1*}.r1.pe ${x%_1*}.r1.se ${x%_1*}.r2.pe ${x%_1*}.r2.se ILLUMINACLIP:../Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;
done 
```



Here's the breakdown of these commands:
java: run java program, -jar: run jar program, /usr/local/bin/trimmomatic-0.36.jar: name of the program with path, PE: paired-end, SRR492066_1.sub.fastq.gz: first pared-end, SRR492066_2.sub.fastq.gz: second paired-end, s1_pe: output of first file(paired-end), s1_se: output of first file(single-end), s2_pe: output of second file(paired-end), s2_se: output of second file(single-end), ILLUMINACLIP: use illumina clip, /usr/local/share/adapters/TruSeq2-PE.fa: adapter file with path, 2:30:10 : <seed mismatches(specifies the maximum mismatch count which will still allow a full match to be performed)>:<palindrome clip threshold(specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.)>:<simple clip threshold(specifies how accurate the match between any adapter etc. sequence must be against a read.)> [Here more detail](http://www.usadellab.org/cms/?page=trimmomatic)

Then, combine
```
cat *.r1.pe > all_trim.r1.fastq
cat *.r2.pe > all_trim.r2.fastq
cat *.se > all_trim.single.fastq
```


What's the quality of the resulting dataset?


```
fastx_quality_stats -i all_trim.r1.fastq -o after_trim.quality.txt
cat after_trim.quality.txt
```


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


##Help and other Resources
* [khmer documentation and repo](https://github.com/dib-lab/khmer/blob/master/README.rst)
* [khmer protocols - see "Kalamazoo Protocol"](http://khmer-protocols.readthedocs.org/en/v0.8.4/)
* [khmer recipes - bite-sized tasks using khmer scripts](http://khmer-recipes.readthedocs.org/en/latest/)
* [khmer discussion group](http://lists.idyll.org/listinfo/khmer)


