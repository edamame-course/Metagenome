# Evaluating or Assessing an Assembly 

Authored by Jin Choi for EDAMAME 

## Summary

### Overarching Goal
* This tutorial will contribute towards an understanding the assembly of **metagenome data**
* This tutorial will discuss assembly evaluation

### Learning Objectives
* Understanding methods to evaluate an assembly
* Understanding assembly metrics
* Comparing two different assemblies 

## Tutorial
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
[STAT] 3669 contigs, total 2977628 bp, min 211 bp, max 37382 bp, avg 812 bp, N50 942 bp
```

Megahit actually runs numerous iterations of assemblies. Here we are looking at the last iteration and the number of contigs and total length of that iteration of the assembly. Later, we will be using QUAST to calculate all of our assembly stats.  There may be some odd terminology in some of the log file, [see this wiki for more information on terminology arising from the Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology) (the first assembly program designed for the Human Genome Project)


## Calculating summary statistics for our Metagenome assembly, using QUAST. 

Quast is a program that will calculate some statistics about our metagenome assembly to give us an idea how well it assembled. Again, we do need to quickly install Quast and one of its dependencies so we can get it running on our machine. This shouldn't take very long. Copy each of the following lines of code one line at a time to install Quast.
#### Install Quast
```
cd
curl -L https://sourceforge.net/projects/quast/files/quast-4.5.tar.gz/download > quast-4.5.tar.gz
tar xvf quast-4.5.tar.gz
```

## Getting the data
Now, let’s create a working directory and download assembled file. We are going to run assessment assembled with full data not subsampled.:
```
cd ~/metagenome
mkdir assessment
cd assessment
wget https://s3.amazonaws.com/edamame/compare_assembly.tar.gz
tar -zxvf compare_assembly.tar.gz
```

Now we can take a look at our assembly using QUAST. **From the ~/metagenomics/assessment** run the following line of code. 
##Looking at the assembly
Run QUAST:
```
~/quast-4.5/quast.py pe.se.final.contigs.fa -o report
```
Once QUAST has finished running, change into the quast_output directory and use `ls` to take a look at all of the files it created. Use `less` to examine the `report.txt` file. 
```
cat report/report.txt
```
You should see:
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    pe.se.final.contigs
# contigs (>= 0 bp)         2380               
# contigs (>= 1000 bp)      883                
# contigs (>= 5000 bp)      414                
# contigs (>= 10000 bp)     244                
# contigs (>= 25000 bp)     107                
# contigs (>= 50000 bp)     42                 
Total length (>= 0 bp)      12941251           
Total length (>= 1000 bp)   12239621           
Total length (>= 5000 bp)   11110462           
Total length (>= 10000 bp)  9873116            
Total length (>= 25000 bp)  7690252            
Total length (>= 50000 bp)  5429748            
# contigs                   1335               
Largest contig              428449             
Total length                12554273           
GC (%)                      39.44              
N50                         38587              
N75                         12783              
L50                         62                 
L75                         204                
# N's per 100 kbp           0.00      
```
N50 could be an important number. L50 (L75, LG50, LG75) is the number of contigs equal to or longer than N50 (N75, NG50, NG75)
In other words, L50, for example, is the minimal number of contigs that cover half the assembly.

##Comparing and evaluating assemblies 
Let's compare three assembly. The result from last tutorial, which used only paired end `pe.final.contigs.fa`, include single end `pe.se.final.contigs.fa`and assembly from raw read `raw.final.contigs.fa`. Those assebmly is done from full data (instead of substituted data). Note, When you compare assembly, all assembly sould be done from same dataset other than that, it dose not make sense. 
```
~/quast-4.5/quast.py pe.final.contigs.fa se.final.contigs.fa pe.se.final.contigs.fa -o report_compare
```
then, open the result:
```
cat report_compare/report.txt
```
You will see.
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    pe.final.contigs  se.final.contigs  pe.se.final.contigs
# contigs (>= 0 bp)         2367              2926              2380               
# contigs (>= 1000 bp)      896               1285              883                
# contigs (>= 5000 bp)      413               485               414                
# contigs (>= 10000 bp)     246               257               244                
# contigs (>= 25000 bp)     108               107               107                
# contigs (>= 50000 bp)     43                37                42                 
Total length (>= 0 bp)      12924337          12677989          12941251           
Total length (>= 1000 bp)   12226262          11881360          12239621           
Total length (>= 5000 bp)   11063653          9949140           11110462           
Total length (>= 10000 bp)  9861702           8340326           9873116            
Total length (>= 25000 bp)  7647216           6086492           7690252            
Total length (>= 50000 bp)  5408515           3688622           5429748            
# contigs                   1354              1867              1335               
Largest contig              427295            238482            428449             
Total length                12548148          12283090          12554273           
GC (%)                      39.44             39.54             39.44              
N50                         38068             24532             38587              
N75                         12765             7009              12783              
L50                         64                110               62                 
L75                         207               362               204                
# N's per 100 kbp           0.00              0.00              0.00 
```

Another way to evaluate an assembly is using read orientation.  Paired end sequencing are two reads that we prepare to be a specific distance apart and in opposite directions.  We can look at mapped paired end reads to a reference genome or an assembled metagenome to assess "paired end concordance."  But first, you need to know how to map reads...
