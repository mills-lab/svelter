#Run different algorithms on simulated data:
We compared SVelter to three other algorithms: Delly, Lumpy and Pindel  for evaluation with following settings:

##SVelter:
index reference genome first:
```
SVelter.py Index --exclude Exclude.GRCh37.bed --reference human_g1k_v37.fasta --workdir workdir/directory --copyneutral CN2.GRCh37.bed --svelter-path ../svelter-master
```
then run the main function of SVelter:
```
SVelter.py --workdir workding/directory --sample input.bam
```


##Delly: 
```
delly -t DEL -s 10 -x human.hg19.excl.tsv -o output.vcf -g human_g1k_v37.fasta input.sorted.bam
delly -t DUP -s 10 -x human.hg19.excl.tsv -o output.vcf -g human_g1k_v37.fasta input.sorted.bam
delly -t INV -s 10 -x human.hg19.excl.tsv -o output.vcf -g human_g1k_v37.fasta input.sorted.bam
delly -t TRA -s 10 -x human.hg19.excl.tsv -o output.vcf -g human_g1k_v37.fasta input.sorted.bam
```


##Lumpy:
first we extract aberrant reads through lumpy module: split_unmapped_to_fasta.pl
```
samtools view input.sorted.bam| ../lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > Lumpy.um.fq
```
then we align, sort and index it:
```
bwa bwasw -H -t 20 human_g1k_v37.fasta Lumpy.um.fq | samtools view -Sb -> Lumpy.sr.bam
samtools sort Lumpy.sr.bam Lumpy.sr.sorted
samtools index Lumpy.sr.sorted.bam
```
run lumpy first to calculate IL distribution:
```
samtools view input.bam | tail -n+100000 | ../lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o Lumpy.histo
```
then solve the SV:
```
../lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x Exclude.bed -pe bam_file:input.sorted.bam,histo_file:/Lumpy.histo ,mean:ILMean ,stdev:ILStd,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 –sr bam_file:Lumpy.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > Lumpy.pesr.bedpe
```


##Pindel:
```
pindel -f human_g1k_v37.fasta -i input.config.txt -c ALL -o input.bam.pindel
pindel2vcf -P input.bam.pindel -r human_g1k_v37.fasta -R human_g1k_v37 -d 20150901 -v Pindel.vcf
```
example of input.config.txt:
```
/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam	500	het.RD10.sorted.bam
```

#### Detailed Code Used:
######SVelter
```
```
######Delly
```
```
######Lumpy
```
```
######Pindel
```
```
