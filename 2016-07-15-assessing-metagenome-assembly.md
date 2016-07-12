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
```

Now we can take a look at our assembly using QUAST. **From the ~/metagenomics/assessment** run the following line of code. 
##Looking at the assembly
Run QUAST:
```
~/quast-3.0/quast.py pe.final.contigs.fa -o report_pe
```
Once QUAST has finished running, change into the quast_output directory and use `ls` to take a look at all of the files it created. Use `less` to examine the `report.txt` file. 
```
less report_pe/report.txt
```
You should see:
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   pe.final.contigs
# contigs (>= 0 bp)        2982            
# contigs (>= 1000 bp)     1295            
Total length (>= 0 bp)     12952450        
Total length (>= 1000 bp)  12137548        
# contigs                  1855            
Largest contig             215667          
Total length               12533749        
GC (%)                     39.42           
N50                        17908           
N75                        7410            
L50                        151             
L75                        427             
# N's per 100 kbp          0.00 
```

##Comparing and evaluating assemblies 
Let's compare three assembly. The result from last tutorial, which used only paired end 'pe.final.contigs.fa', include single end 'pe.se.final.contigs.fa'and assembly from raw read 'raw.final.contigs.fa'. Those assebmly is done from full data (instead of substituted data). Note, When you compare assembly, all assembly sould be done from same dataset other than that, it dose not make sense. 
```
~/quast-3.0/quast.py pe.final.contigs.fa pe.se.final.contigs.fa raw.final.contigs.fa -o report_compare
```
then, open the result:
```
less report_compare/report.txt
```
You will see.
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   pe.final.contigs  pe.se.final.contigs  raw.final.contigs
# contigs (>= 0 bp)        2982              3284                 2976             
# contigs (>= 1000 bp)     1295              1447                 1149             
Total length (>= 0 bp)     12952450          12744624             13050315         
Total length (>= 1000 bp)  12137548          11835839             12189601         
# contigs                  1855              2127                 1719             
Largest contig             215667            157021               232493           
Total length               12533749          12313947             12590769         
GC (%)                     39.42             39.56                39.47            
N50                        17908             17634                22360            
N75                        7410              6220                 8751             
L50                        151               150                  132              
L75                        427               453                  357              
# N's per 100 kbp          0.00              0.00                 0.00

```

