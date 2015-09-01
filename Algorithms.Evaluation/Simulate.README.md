Simulated Simple SVs
We used our own script to produce altered reference genome with pre-set simple and complex structural variants:


Usage:
Produce.Simulated.FussyJuncs.py [options] <parameters>

 
Options:
heterozygous:	simulate simple heterozygous SVs
homozygous:	simulate simple homozygous SVs
complex:		simulate complex SVs
 
Parameters:
--reference: reference genme
--input-sim: input sim format,see example
--input-rec: input rec format, specially designed for complex events,see example
--output-prefix: prefix of output files

Examples:
Produce.Simulated.FussyJuncs.py heterozygous \
--reference human_g1k_v37.fasta \
--input-sim het.sim  \
--output-prefix simple_het

Produce.Simulated.FussyJuncs.py heterozygous \
--reference human_g1k_v37.fasta \
--input-sim comp_homo..sim  \
--input- rec comp_homo.rec \
--output-prefix comp_homo.

Example files in Supp1 folder


Detailed number of simulated events were integrated in Supp.table1.

Then we simulated paired end reads based on each altered genome up to different read depth (10X-50X) with wgsim[]:

wgsim -e 0.001 -d 500 -s 100 -N num.of.reads -1 101 -2 101 -r 0.001 -R 0.1 -X 0 altered.reference.fa out1.fq out2.fq

and aligned reads to GRCh37 before sort and index the aligned files in bam format:

bwa mem human_g1k_v37.fasta out1.fq out2.fq input.sam
samtools view -h –Sb input.sam -o input.bam
samtools sort input.bam input.sorted
samtools input.sorted.bam


Run different algorithms on simulated data:

We compared SVelter to three other algorithms: Delly, Lumpy and Pindel  for evaluation with following settings:
Delly: 
delly -t [DEL/DUP/INV/TRA] -s 10 -x human.hg19.excl.tsv -o output.vcf -g human_g1k_v37.fasta input.sorted.bam


Lumpy:
first we extract aberrant reads through lumpy module: split_unmapped_to_fasta.pl
samtools view input.sorted.bam| ../lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > Lumpy.um.fq

then we align, sort and index it:
bwa bwasw -H -t 20 human_g1k_v37.fasta Lumpy.um.fq | samtools view -Sb -> Lumpy.sr.bam
samtools sort Lumpy.sr.bam Lumpy.sr.sorted
samtools index Lumpy.sr.sorted.bam

run lumpy first to calculate IL distribution:
samtools view input.bam | tail -n+100000 | ../lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o Lumpy.histo

turn solve the SV:
../lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x Exclude.bed -pe bam_file:input.sorted.bam,histo_file:/Lumpy.histo ,mean:ILMean ,stdev:ILStd,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 –sr bam_file:Lumpy.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > Lumpy.pesr.bedpe


Pindel:
pindel -f human_g1k_v37.fasta -i input.config.txt -c ALL -o input.bam.pindel
pindel2vcf -P input.bam.pindel -r human_g1k_v37.fasta -R human_g1k_v37 -d 20150901 -v Pindel.vcf

SVelter:
index reference genome first:
SVelter.py Index --exclude Exclude.GRCh37.bed --reference human_g1k_v37.fasta --workdir workdir/directory --copyneutral CN2.GRCh37.bed --svelter-path ../svelter-master
then run SVelter:
SVelter.py --workdir workding/directory --sample input.bam


