#Estimating abundances in metagenomes

Authored by Jin Choi for EDAMAME2016 

## Summary

### Overarching Goal
* This tutorial will contribute towards an understanding of quantitative analyses of **metagenome data**
* It focuses on estimating abundances of reads to an assembled reference.

### Learning Objectives
* Understanding how to estimate abundances of reads in a representative gene reference
* Understanding read mapping
* Understanding mapping file formats
* Understanding how to use a mapping program (Bowtie2, samtools, bcftools)
* Apply reference mapping to assess read abundances and quantify gene presence

### Reference mapping is useful... 
* If you want to detect SNPs
* If you want to estimate abundance of genes in metagenomic or metatranscriptomic data

## Tutorial

### get gene calling
```
curl -o mgm4753635.3.350.genecalling.coding.faa -X GET "http://api.metagenomics.anl.gov/1/download/mgm4753635.3?file=350.1”
```

### make gtf
```
python dev/assembly_downstream/350_to_gtf.py mgm4753635.3.350.genecalling.coding.faa > assmbly.gtf
```

### count using HTSeq
```
for x in *.sorted.bam;do htseq-count -i gene_id -f bam $x assmbly.gtf > $x.htseq.count;done
```

### merge
```
python dev/assembly_downstream/htseq_count_table.py *.count > final.table
```

### download annotation
```
curl -o mgm4753635.3.function_SEED.tab -X GET "http://api.metagenomics.anl.gov/1/annotation/similarity/mgm4753635.3?type=function&source=SEED"
```
