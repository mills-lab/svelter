#Simulated Simple SVs
We used our own script *Produce.Simulated.FussyJuncs.py* to modify the germline reference genome with pre-designed simple and complex SVs in both homozygous and heterozygous status. Output files includes altered reference genomes in fasta format, and detailed record of implemented SVs. 

We also adopted the simulated cancer and matched normal references reported by Moncunill et al to evaluate performance of SVelter on somatic events. Those references are publicly available from [SMuFin](http://cg.bsc.es/smufin/)


##Usage of *Produce.Simulated.FussyJuncs.py*:
Produce.Simulated.Ref.FussyJuncs.py [options] <parameters>
 
####Options:
```
heterozygous:	simulate simple heterozygous SVs
homozygous:	simulate simple homozygous SVs
complex:		simulate complex SVs
 ```
 
####Parameters:
````
--reference: reference genme
--input-sim: input sim format,see example
--input-rec: input rec format, specially designed for complex events,see example
--output-prefix: prefix of output files
````

####Examples:
```
Produce.Simulated.Ref.FussyJuncs.py heterozygous --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --input-sim /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/het.sim --output-prefix /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/simple_het
Produce.Simulated.Ref.FussyJuncs.py homozygous --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --input-sim /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/homo.sim --output-prefix /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/simple_homo
Produce.Simulated.Ref.FussyJuncs.py complex --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --input-sim /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_het.sim --input-rec /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_het.rec --output-prefix /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_het
Produce.Simulated.Ref.FussyJuncs.py complex --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --input-sim /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_homo.sim --input-rec /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_homo.rec --output-prefix /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_homo
```

Example files in Supp1 folder; Detailed number of simulated events were integrated in Supp.table1.

Then we simulated paired end reads based on each altered genome up to different read depth (10X-50X) with [wgsim](https://github.com/lh3/wgsim):
```
wgsim -e 0.001 -d 500 -s 100 -N num.of.reads -1 101 -2 101 -r 0.001 -R 0.1 -X 0 altered.reference.fa out1.fq out2.fq
```

For the simulated tumor and matched normal references, we simulated paired end reads with [art_illumina](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/), adoping parameters used by [*Moncunill et al*] (http://www.nature.com/nbt/journal/v32/n11/full/nbt.3027.html):
```
art_illumina -sam -i altered.reference.fa -p -l 80 -ss HS20 -f 20 -m 500 -s 20 -o output 
```

and aligned reads to GRCh37 before sort and index the aligned files in bam format:
```
bwa mem human_g1k_v37.fasta out1.fq out2.fq input.sam
samtools view -h –Sb input.sam -o input.bam
samtools sort input.bam input.sorted
samtools input.sorted.bam
```


###Detailed Code Used:
<!---
####Step1: Simulate SVs according to pre-defined sizes and numbers to form altered reference genome
```
Produce.Simulated.Ref.FussyJuncs.py heterozygous --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta --input-sim /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/het.sim --output-prefix /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het
Produce.Simulated.Ref.FussyJuncs.py homozygous --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta --input-sim /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/homo.sim --output-prefix /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo
Produce.Simulated.Ref.FussyJuncs.py complex --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta --input-sim /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het.sim --input-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het.rec --output-prefix /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het
Produce.Simulated.Ref.FussyJuncs.py complex --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta --input-sim /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo.sim --input-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo.rec --output-prefix /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo
```

####Step2: Simulate each altered genome up to different read depth(10X-50X)

######For **Heterozygous Simple** events
```
wgsim  -e 0.001 -d 500 -s 100 -N 76777345 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het.het1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD10.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 76777345 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD10.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 153554690 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD20.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 153554690 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD20.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 230332035 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD30.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 230332035 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD30.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 307109380 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD40.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 307109380 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD40.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 383886725 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD50.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 383886725 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD50.2.fq
```

######For **Homozygous Simple** events
```
wgsim  -e 0.001 -d 500 -s 100 -N 76777345 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo.homo.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD10.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 153554690 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.homo.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD20.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 230332035 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.homo.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD30.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 307109380 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.homo.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD40.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 383886725 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.homo.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD50.2.fq
```

######For **Heterozygous Complex** events
```
wgsim  -e 0.001 -d 500 -s 100 -N 76777345 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD10.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 76777345 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD10.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 153554690 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD20.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 153554690 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD20.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 230332035 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD30.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 230332035 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD30.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 307109380 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD40.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 307109380 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD40.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 383886725 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD50.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 383886725 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD50.2.fq
```

######For **Homozygous Complex** events
```
wgsim  -e 0.001 -d 500 -s 100 -N 76777345 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD10.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 76777345 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD10.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 153554690 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD20.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 153554690 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD20.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 230332035 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD30.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 230332035 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD30.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 307109380 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD40.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 307109380 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD40.2.fq

wgsim  -e 0.001 -d 500 -s 100 -N 383886725 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp1.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD50.2.fq
wgsim  -e 0.001 -d 500 -s 100 -N 383886725 -1 101 -2 101 -r 0.001 -R 0.1 -X 0.01 /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_new.comp2.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD50.2.fq
```

####Step3: concatinate reads simulated from a pair of references 
######For **Heterozygous Simple** events

```
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD10.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD10.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD10.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD10.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD10.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD10.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD20.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD20.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD20.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD20.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD20.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD20.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD30.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD30.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD30.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD30.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD30.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD30.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD40.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD40.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD40.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD40.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD40.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD40.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD50.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD50.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD50.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref1.RD50.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.ref2.RD50.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD50.2.fq 
```
######For **Heterozygous Complex** events
```
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD10.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD10.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD10.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD10.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD10.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD10.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD20.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD20.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD20.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD20.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD20.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD20.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD30.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD30.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD30.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD30.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD30.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD30.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD40.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD40.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD40.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD40.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD40.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD40.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD50.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD50.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD50.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref1.RD50.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.ref2.RD50.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD50.2.fq 
```
######For **Homozygous Complex** events
```
####Step3: concatinate reads simulated from a pair of references 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD10.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD10.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD10.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD10.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD10.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD10.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD20.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD20.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD20.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD20.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD20.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD20.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD30.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD30.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD30.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD30.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD30.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD30.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD40.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD40.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD40.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD40.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD40.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD40.2.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD50.1.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD50.1.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD50.1.fq 
cat /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref1.RD50.2.fq  /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.ref2.RD50.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD50.2.fq 
```
####Step4: alignment
######For **Heterozygous Simple** events
```
bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD10.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD10.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD10.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD10.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD10.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD10.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD10.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD20.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD20.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD20.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD20.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD20.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD20.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD20.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD30.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD30.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD30.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD30.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD30.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD30.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD30.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD40.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD40.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD40.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD40.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD40.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD40.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD40.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/fastq/simple_het.RD50.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD50.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD50.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD50.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD50.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD50.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_het/alignment/simple_het.RD50.sorted.bam

```
######For **Homozygous Simple** events
```
bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD10.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD10.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD10.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD10.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD10.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD10.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD10.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD20.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD20.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD20.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD20.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD20.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD20.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD20.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD30.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD30.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD30.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD30.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD30.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD30.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD30.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD40.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD40.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD40.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD40.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD40.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD40.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD40.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/fastq/simple_homo.RD50.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD50.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD50.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD50.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD50.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD50.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/simple_homo/alignment/simple_homo.RD50.sorted.bam
```
######For **Heterozygous Complex** events
```
bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD10.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD10.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD10.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD10.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD10.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD10.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD10.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD20.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD20.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD20.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD20.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD20.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD20.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD20.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD30.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD30.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD30.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD30.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD30.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD30.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD30.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD40.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD40.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD40.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD40.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD40.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD40.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD40.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/fastq/comp_het.RD50.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD50.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD50.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD50.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD50.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD50.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_het/alignment/comp_het.RD50.sorted.bam
```
######For **Homozygous Complex** events
```
bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD10.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD10.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD10.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD10.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD10.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD10.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD10.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD10.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD20.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD20.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD20.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD20.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD20.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD20.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD20.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD20.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD30.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD30.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD30.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD30.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD30.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD30.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD30.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD30.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD40.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD40.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD40.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD40.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD40.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD40.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD40.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD40.sorted.bam

bwa mem /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD50.1.fq /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/fastq/comp_homo.RD50.2.fq > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD50.sam
samtools view -h -Sb /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD50.sam -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD50.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD50.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD50.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.rerun.test.20150901/comp_homo/alignment/comp_homo.RD50.sorted.bam
```
-->
