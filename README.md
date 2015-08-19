SVelter
======================================================

Description
-----------

This software is designed to identify both simple and complex rearrangements from paired-end sequencing data. Users could ran it easily by just alling 'SVelter.py' with proper parameters. It's also possible to ran it on multiple cores (>4) by calling a series of python scripts separately.


Required third-party resources 
------------------------------
python 
R

Usage
----------
Reference genome should be indexed first:

SVelter.py Index --Reference /absolute/path/of/reference.fa --WorkDir /working/directory --EX /absolute/path/of/exclude.ref.bed --CN /absolute/path/of/CN2.ref.bed --SVelterPath /absolute/path/of/SVelter/scripts

requirements:

1.reference file should have been indexed by calling samtools first: samtools faidx ref.fasta

2. working directory is required to be writable for temporal files 


To run SVelter with its default setting:

SVelter.py --Sample /absolute/path/of/sample.bam --WorkDir /working/directory


optional parameters:

--Output absolute path of output.vcf

--QCAlign mappling quality cutoff

--Ploidy 0 for heterozygous events only / 1 for diploid homozygous genome/ 2 for diploid heterozygous genome

--NullModel 'S' for simplier but faster models / 'C' for more complex models



To run SVelter with multiple cores:

Step1: Build null models:

SVelter.py NullModel --Sample /absolute/path/of/sample.bam --WorkDir /working/directory


Step2: Search for Breakpoints:

SVelter.py BPSearch --Sample /absolute/path/of/sample.bam --WorkDir /working/directory


Step3: Cluster Breakpoints:

SVelter.py BPIntegrate --Sample /absolute/path/of/sample.bam --WorkDir /working/directory --Batch 0 


Step4: Resolve complex structural variants:

SVelter.py SVPredict --Sample /absolute/path/of/sample.bam --WorkDir /working/directory


Step5: Write output in vcf format:

SVelter.py SVIntegrate --Sample /absolute/path/of/sample.bam --WorkDir /working/directory --Output Output.vcf  --RSPath path/of/output/from/Step4
  
                                       
