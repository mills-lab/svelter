SVelter
======================================================

Description
-----------

This software is designed to identify both simple and complex rearrangements from paired-end sequencing data. Users could ran it easily by just alling 'SVelter.py' with proper parameters. It's also possible to ran it on multiple cores (>4) by calling a series of python scripts separately.


Required third-party resources 
------------------------------
python 

Usage
----------
Reference genome should be indexed first:

SVelter0.Ref.Index.py --ref ref.fa --ppre workdir

requirements:
1.reference file should have been indexed by calling samtools first: samtools faidx ref.fasta
2. the tools/ folder should have been copied or linked under working directory
3. SVelter0.Ref.Index.py should be called separately for different working directories


To run SVelter with its default setting:
SVelter.py [options]

Required parameters:

--ppre workdir,should be writable for temporary files

--ScriptPath folder of all SVelter*.py 

-s absolute path of input sample.bam  

--ref absolute path of ref.fa


optional parameters:

-o absolute path of output.vcf

--ex genomic regions to be excluded from downstream analysis.bed

--cn copy wneutural genomic regions.bed

--QCAlign mappling quality cutoff

--MoveFlag 0 for heterozygous events only / 1 for diploid homozygous genome/ 2 for diploid heterozygous genome

--NullModel 'S' for simplier but faster models / 'C' for more complex models



To run SVelter with multiple cores:

Step1: Build null models:

SVelter1.NullModel.py --ppre wirkdir -s sample.bam --ref ref.fa 


Step2: Search for Breakpoints:

SVelter2.BP.Searching.py --ppre wirkdir -s sample.bam --ref ref.fa --chr chromosome


Step3: Cluster Breakpoints:

SVelter3.BPIntegrate.py --ppre wirkdir -s sample.bam --ref ref.fa --batch 0


Step4: Resolve complex structural variants:

SVelter4.StructureResolvation.py --ppre wirkdir -s sample.bam --ref ref.fa -f input.txt


Step5: Write output in vcf format:

SVelter5.result.integrate.py --ppre wirkdir --ref ref.fa --RSPath path/of/input.txt -o output.vcf
  
                                       
