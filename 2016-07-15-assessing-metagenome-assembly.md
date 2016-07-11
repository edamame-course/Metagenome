# Evaluating or Assessing an Assembly 

Authored by Jin Choi for EDAMAME2016 

## Summary

### Overarching Goal
* This tutorial will contribute towards an understanding the assembly of **metagenome data**

### Learning Objectives
* Understanding methods to evaluate an assembly
* Understanding assembly metrics
* Comparing two different assemblies 

## Tutorial
@Jin add details here - what is goign on -- do we want to assemble same dataset on two assemblers and compare?  Do we want to first get assembly statistics on one assembly?  Can we add paired info on the mapping tutorial or here?  That is a great method to evaluate assemblies.


## Install software for this tutorial

Install khmer
```
cd
python -m virtualenv env
source env/bin/activate
pip install -U setuptools
pip install khmer==1.4.1
```
and download and compile the SPAdes assembler:
```
cd
curl -O http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0.tar.gz
tar xvf SPAdes-3.5.0.tar.gz
cd SPAdes-3.5.0
./spades_compile.sh
export PATH="$PATH:$(pwd)/bin"
```
as well as Quast, software for evaluating the assembly against the known reference:
```
cd
curl -L http://sourceforge.net/projects/quast/files/quast-3.0.tar.gz/download > quast-3.0.tar.gz
tar xvf quast-3.0.tar.gz
```
## Getting the data
Now, let’s create a working directory:
```
cd /mnt
mkdir assembly
cd assembly
```

Download some E. coli data. This data set (ecoli_ref-5m-trim.fastq.gz) is the trimmed data from the Chitsaz paper, E. coli reference sequencing.
```
curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoli_ref-5m-trim.fastq.gz
```
Now, pull out the paired reads:
```
extract-paired-reads.py ecoli_ref-5m-trim.fastq.gz
mv ecoli_ref-5m-trim.fastq.gz.se ecoli_ref-5m-trim.se.fq
mv ecoli_ref-5m-trim.fastq.gz.pe ecoli_ref-5m-trim.pe.fq
```

## Running an assembly
Now, let’s run an assembly:
```
spades.py --12 ecoli_ref-5m-trim.pe.fq -s ecoli_ref-5m-trim.se.fq -o spades.d
```
This will take about 15 minutes; it should end with:
```
* Corrected reads are in /mnt/assembly/spades.d/corrected/
* Assembled contigs are in /mnt/assembly/spades.d/contigs.fasta (contigs.fastg)
* Assembled scaffolds are in /mnt/assembly/spades.d/scaffolds.fasta (scaffolds.fastg)
```
##Looking at the assembly
Run QUAST:
```
~/quast-3.0/quast.py spades.d/scaffolds.fasta -o report
```
and then look at the report:
```
less report/report.txt

```
You should see:
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        160
# contigs (>= 1000 bp)     84
Total length (>= 0 bp)     4571783
Total length (>= 1000 bp)  4551354
# contigs                  93
Largest contig             264754
Total length               4557807
GC (%)                     50.75
N50                        132618
N75                        64692
L50                        12
L75                        24
# N's per 100 kbp          0.00
```

##Comparing and evaluating assemblies - QUAST
Download the true reference genome:
```
cd /mnt/assembly
curl -O https://s3.amazonaws.com/public.ged.msu.edu/ecoliMG1655.fa.gz
gunzip ecoliMG1655.fa.gz
```

and run QUAST again:
```
~/quast-3.0/quast.py -R ecoliMG1655.fa spades.d/scaffolds.fasta -o report
```
Note that here we’re looking at all the assemblies we’ve generated.

Now look at the results:
```
less report/report.txt
```
and now we have a lot more information!

## A second assembler - MEGAHIT
Let’s try out the MEGAHIT assembler. MEGAHIT is primarily intended for metagenomes but works well on microbial genomes in general.

The MEGAHIT source code is on GitHub, here: https://github.com/voutcn/megahit. Let’s go grab it and build it!

```
cd
git clone https://github.com/voutcn/megahit.git
cd megahit
make
```
Now, let’s go run an assembly –
```
cd /mnt/assembly
~/megahit/megahit --12 *.pe.fq -r *.se.fq
```
This will take about a minute, and the output will be placed in megahit_out/final.contigs.fa. Let’s evaluate it against the SPAdes assembly with QUAST:
```
cp spades.d/scaffolds.fasta spades-assembly.fa
cp megahit_out/final.contigs.fa megahit-assembly.fa
~/quast-3.0/quast.py -R ecoliMG1655.fa spades-assembly.fa \
         megahit-assembly.fa -o report
```
Let’s look at the report!
```
less report/report.txt
```

## Reference-free comparison
Above, we’ve been using the genome reference to do assembly comparisons – but often you don’t have one. What do you do to evaluate and compare assemblies without a reference?

One interesting trick is to just run QUAST with one assembly as a reference, and the other N assemblies against it. My only suggestion is to first eliminate short, fragmented contigs from the assembly you’re going to use as a reference.

Let’s try that, using extract-long-sequences.py from khmer:

```
extract-long-sequences.py -l 1000 spades-assembly.fa > spades-long.fa
```
and then re-run QUAST and put the output in report-noref/report.txt:
```
~/quast-3.0/quast.py -R spades-long.fa spades-assembly.fa \
         megahit-assembly.fa -o report-noref
```
When you look at the report,

```
less report-noref/report.txt
```

take particular note of the following –

```
Assembly                     spades-assembly  megahit-assembly
...
Misassembled contigs length  0                814643
# local misassemblies        0                9
# unaligned contigs          9 + 0 part       7 + 14 part
Unaligned length             6453             7404
Genome fraction (%)          100.000          99.833
```

#here
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
  
