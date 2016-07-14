# Evaluating or Assessing an Assembly 

Authored by Jin Choi for EDAMAME2016 

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
[STAT] 3759 contigs, total 2935387 bp, min 206 bp, max 35461 bp, avg 781 bp, N50 885 bp
```

Megahit actually runs numerous iterations of assemblies. Here we are looking at the last iteration and the number of contigs and total length of that iteration of the assembly. Later, we will be using QUAST to calculate all of our assembly stats.  There may be some odd terminology in some of the log file, [see this wiki for more information on terminology arising from the Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology) (the first assembly program designed for the Human Genome Project)


## Calculating summary statistics for our Metagenome assembly, using QUAST. 

Quast is a program that will calculate some statistics about our metagenome assembly to give us an idea how well it assembled. Again, we do need to quickly install Quast and one of its dependencies so we can get it running on our machine. This shouldn't take very long. Copy each of the following lines of code one line at a time to install Quast.
#### Install Quast
```
cd
curl -L http://sourceforge.net/projects/quast/files/quast-3.0.tar.gz/download > quast-3.0.tar.gz
tar xvf quast-3.0.tar.gz
```

## Getting the data
Now, let’s create a working directory:
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
~/quast-3.0/quast.py pe.final.contigs.fa -o report_pe
```
Once QUAST has finished running, change into the quast_output directory and use `ls` to take a look at all of the files it created. Use `less` to examine the `report.txt` file. 
```
cat report_pe/report.txt
```
You should see:
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   pe.final.contigs
# contigs (>= 0 bp)        2707            
# contigs (>= 1000 bp)     1064            
Total length (>= 0 bp)     12946382        
Total length (>= 1000 bp)  12168690        
# contigs                  1581            
Largest contig             232526          
Total length               12532010        
GC (%)                     39.45           
N50                        30749           
N75                        10062           
L50                        95              
L75                        275             
# N's per 100 kbp          0.00   
```

##Comparing and evaluating assemblies 
Let's compare three assembly. The result from last tutorial, which used only paired end `pe.final.contigs.fa`, include single end `pe.se.final.contigs.fa`and assembly from raw read `raw.final.contigs.fa`. Those assebmly is done from full data (instead of substituted data). Note, When you compare assembly, all assembly sould be done from same dataset other than that, it dose not make sense. 
```
~/quast-3.0/quast.py pe.final.contigs.fa se.final.contigs.fa pe.se.final.contigs.fa -o report_compare

raw.final.contigs.fa -o report_compare
```
then, open the result:
```
cat report_compare/report.txt
```
You will see.
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   pe.final.contigs  se.final.contigs  pe.se.final.contigs
# contigs (>= 0 bp)        2707              3162              2670               
# contigs (>= 1000 bp)     1064              1375              1021               
Total length (>= 0 bp)     12946382          12735336          12958914           
Total length (>= 1000 bp)  12168690          11849449          12177316           
# contigs                  1581              2038              1537               
Largest contig             232526            232458            282362             
Total length               12532010          12316301          12541412           
GC (%)                     39.45             39.56             39.44              
N50                        30749             22507             32048              
N75                        10062             6464              10401              
L50                        95                113               86                 
L75                        275               392               258                
# N's per 100 kbp          0.00              0.00              0.00 
```
Let's compare asssembly using RAW read.
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   raw.final.contigs
# contigs (>= 0 bp)        2976             
# contigs (>= 1000 bp)     1149             
Total length (>= 0 bp)     13050315         
Total length (>= 1000 bp)  12189601         
# contigs                  1719             
Largest contig             232493           
Total length               12590769         
GC (%)                     39.47            
N50                        22360            
N75                        8751             
L50                        132              
L75                        357              
# N's per 100 kbp          0.00 
```

Another way to evaluate an assembly is using read orientation.  Paired end sequencing are two reads that we prepare to be a specific distance apart and in opposite directions.  We can look at mapped paired end reads to a reference genome or an assembled metagenome to assess "paired end concordance."  But first, you need to know how to map reads...
