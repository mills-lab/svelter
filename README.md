#SVelter

##Description
This software is designed to identify both simple and complex rearrangements from paired-end sequencing data. It is used by calling *SVelter.py* with listed parameters. It's also possible to run it on multiple cores by calling different sub-functions separately.

##Required third-party resources
```
R:        https://www.r-project.org/
python:   https://www.python.org/
samtools: http://samtools.sourceforge.net/
```

## Quick Start
Download and Install
```
git clone git@github.com:mills-lab/svelter.git
cd svelter
chmod +x SVelter.py
```
Index Reference genome
``` 
SVelter.py Setup --reference reference.fa --workdir /working/directory --exclude exclude.ref.bed --copyneutral CN2.ref.bed --ref-index indexed-ref/ --svelter-path SVelter/ 
```
Run SVelter with its default setting:
```
SVelter.py --sample /absolute/path/of/sample.bam --workdir /working/directory
```

##Supportive files  
`exclude.ref.bed` and `CN2.ref.bed` are available from the folder *Support* for some versions of reference genome. Users could replace with their custom version as long as both are in bed format. For more details, please see *Support*. 

Pre-indexed files of certain reference genomes have been produced and kept under folder *index-ref*. For specific reference, if not pre-indexed files provided, the optional parameter '--ref-index' could be omit and the indexed files would be produced through the setup step. 


##Attention:
reference file should have been indexed by calling samtools first:  `samtools faidx ref.fasta`

working directory is required to be writable for temporal files 

##Output
*SVelter* integrates predicted SVs in both vcf4.1 format and SVelter format. Examples of both could be found under folder *Example*

##Usage
SVelter.py  [options]  [parameters]

###Options:
```
  NullModel
  BPSearch
  BPIntegrate
  SVPredict
  SVIntegrate
```

###Parameters:

####Required:
```
  --workdir, writable working directory.
  
  --sample, input alignment file in bam format
```

####Optional:
```
--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used

--null-copyneutral-length, minimum length requirement for --copyneutral regions used to build null model (default: 2000)

--null-copyneutral-perc, percentage of regions from --copyneutral to utilize (default: 0.1)

--null-random-length, specify the length of random regions if --copyneutral parameter not used (default: 5000)

--null-random-num, specify the number of random regions if --copyneutral parameter not used (default: 10000)

--num-iteration, maximum number of iterations per structure will run in SV predicting step

--qc-map-tool, the tool extracts mappability information from a bigWig file,avaliable from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary

--qc-map-file, .bigWig file used to decide local genomic mappability, avaliable from: ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/ 

--qc-map-cutoff, the minimum mapping quality required for a breakpoint to be reported (default: 0.0)

--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)

--qc-split, minimum alighment of clipped parts of reads considered as a soft clip (default: 20)

--qc-structure, minimum quality score of a resolved structure to be considered as PASS and included in the output vcf file

--split-min-len, the minumum length of clip read considered as split; (default:10% of read length)

--prefix, output prefix for vcf and svelter files (default: input.vcf, input.svelter)

--ploidy, limit algorithm to specific zygosity (0:heterozygous only; 1:homozygous only; 2:both; default:2)
```

###Attention:

> reference file should have been indexed by calling samtools first:  `samtools faidx ref.fasta`

> working directory is required to be writable for temporal files 



###For faster processing, SVelter could run with multiple cores:

####Step1: Build null models:
```
SVelter.py NullModel --sample /absolute/path/of/sample.bam --workdir /working/directory
```

```
Optional Parameters:

--chromosome, name of chromosome to run. should match chromosome name in bam file

--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used

--null-copyneutral-length, minimum length requirement for --copyneutral regions used to build null model (default: 2000)

--null-copyneutral-perc, percentage of regions from --copyneutral to utilize (default: 0.1)

--null-random-length, specify the length of random regions if --copyneutral parameter not used (default: 5000)

--null-random-num, specify the number of random regions if --copyneutral parameter not used (default: 10000)

--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)

--qc-split, minimum alighment of clipped parts of reads considered as a soft clip (default: 20)

--split-min-len, the minumum length of clip read considered as split  (default:10% of read length)
```

####Step2: Search for Breakpoints:
```
SVelter.py BPSearch --sample /absolute/path/of/sample.bam --workdir /working/directory
```

```
Optional Parameters:

--chromosome, name of chromosome to run. should match chromosome name in bam file

--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used

--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)

--qc-split, minimum alighment of clipped parts of reads considered as a soft clip (default: 20)

--split-min-len, the minumum length of clip read considered as split; (default:10% of read length)

--qc-map-tool, the tool extracts mappability information from a bigWig file,avaliable from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary

--qc-map-file, .bigWig file used to decide local genomic mappability, avaliable from: ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/

--qc-map-cutoff, the minimum mapping quality required for a breakpoint to be reported (default: 0.0)
```

####Step3: Cluster Breakpoints:
```
SVelter.py BPIntegrate --sample /absolute/path/of/sample.bam --workdir /working/directory
```

```
Optional Parameters:

--chromosome, name of chromosome to run. should match chromosome name in bam file

--batch, specify number of structures in each separate file (if 0, output files will be calssified by chromosomes; default, all BP clustered will be integrated in one txt file)
```

####Step4: Resolve complex structural variants:
```
SVelter.py SVPredict --sample sample.bam --workdir /working/directory --bp-file input/file/containing/clustered/bps
 ```
 
 ```
Optional Parameters:

--num-iteration, maximum number of iterations per structure will run

--ploidy, limit algorithm to specific zygosity (0:heterozygous only; 1:homozygous only; 2:both; default:2)

--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used

--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)
```

####Step5: Write output in vcf and svelter format:
```
SVelter.py SVIntegrate --workdir /working/directory --prefix output  --input-path path/of/output/from/Step4
```

```
Optional Parameters:

--qc-structure, minimum quality score of a resolved structure to be considered as PASS and included in the output vcf file
```

