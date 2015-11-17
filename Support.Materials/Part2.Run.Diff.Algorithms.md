#Run different algorithms on simulated data:
We compared SVelter to three other algorithms: Delly, Lumpy, Pindel and SMuFin for evaluation with following settings:

##SVelter:
index reference genome first:
```
SVelter.py Setup --exclude Exclude.GRCh37.bed --reference human_g1k_v37.fasta --workdir workdir/directory --copyneutral CN2.GRCh37.bed --svelter-path ../svelter-master
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
##SMuFin
```
mpirun --np 30 SMuFin --ref hg19.fa --normal_fastq_1 SMuFin_Normal_1.txt --normal_fastq_2 SMuFin_Normal_2.txt --tumor_fastq_1 SMuFin_Tumor_1.txt --tumor_fastq_2 SMuFin_Tumor_2.txt --patient_id patient_id --cpus_per_node 30
```
#### Detailed Code Used:
<!---
For **Heterozygous Simple** events
######SVelter
```
SVelter.py Setup --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/ --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --exclude /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/Exclude.GRCh37.bed --copyneutral /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/CN2.GRCh37.bed --svelter-path /mnt/EXT/Mills-scratch2/Xuefang/svelter --ref-index /mnt/EXT/Mills-scratch2/reference/SVelter.Index/GRCh37/reference

SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam
```
######Delly
```
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD10.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD10.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD10.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD10.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD20.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD20.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD20.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD20.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD30.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD30.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD30.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD30.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD40.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD40.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD40.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD40.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD50.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD50.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD50.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Delly/het.RD50.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam
```
######Lumpy
```
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD10.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD10.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD20.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD20.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD30.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD30.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD40.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD40.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/BamFiles/het.RD50.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/Lumpy/het.RD50.pesr.bedpe
```
######Pindel
```
pindel2vcf -P het.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/pbs/Pindel/Pindel.het.RD.10.vcf

```

For **Homozygous Simple** events
######SVelter
```
SVelter.py Setup --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/ --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --exclude /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/Exclude.GRCh37.bed --copyneutral /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/CN2.GRCh37.bed --svelter-path /mnt/EXT/Mills-scratch2/Xuefang/svelter --ref-index /mnt/EXT/Mills-scratch2/reference/SVelter.Index/GRCh37/reference

SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam
```
######Delly
```
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD10.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD10.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD10.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD10.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD20.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD20.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD20.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD20.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD30.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD30.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD30.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD30.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD40.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD40.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD40.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD40.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD50.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD50.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD50.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Delly/homo.RD50.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam
```
######Lumpy
```
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD10.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD10.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD20.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD20.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD30.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD30.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD40.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD40.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/BamFiles/homo.RD50.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/Lumpy/homo.RD50.pesr.bedpe
```
######Pindel
```
pindel2vcf -P homo.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

```

For **Heterozygous Complex** events
######SVelter
```
SVelter.py Setup --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/ --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --exclude /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/Exclude.GRCh37.bed --copyneutral /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/CN2.GRCh37.bed --svelter-path /mnt/EXT/Mills-scratch2/Xuefang/svelter --ref-index /mnt/EXT/Mills-scratch2/reference/SVelter.Index/GRCh37/reference

SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam
```
######Delly
```
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD10.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD10.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD10.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD10.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD20.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD20.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD20.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD20.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD30.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD30.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD30.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD30.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD40.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD40.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD40.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD40.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD50.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD50.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD50.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Delly/comp.het.RD50.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam
```
######Lumpy
```
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD10.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD10.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD20.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD20.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD30.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD30.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD40.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD40.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/BamFiles/comp.het.RD50.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Lumpy/comp.het.RD50.pesr.bedpe
```
######Pindel
```
pindel2vcf -P het.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf

pindel2vcf -P het.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf
pindel2vcf -P het.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/pbs/Pindel/Pindel.het.RD.10.vcf
```

For **Homozygous Complex** events
######SVelter
```
SVelter.py Setup --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/ --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --exclude /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/Exclude.GRCh37.bed --copyneutral /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/CN2.GRCh37.bed --svelter-path /mnt/EXT/Mills-scratch2/Xuefang/svelter --ref-index /mnt/EXT/Mills-scratch2/reference/SVelter.Index/GRCh37/reference

SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam
SVeter.py NullModel --workdir /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/ --sample /mnt/EXT/Mills
-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam
```
######Delly
```
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD10.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD10.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD10.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD10.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD20.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD20.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD20.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD20.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD30.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD30.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD30.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD30.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD40.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD40.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD40.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD40.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD50.DEL.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam
delly -t DUP -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD50.DUP.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam
delly -t INV -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD50.INV.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam
delly -t TRA -s 10 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Delly/comp.homo.RD50.TRA.vcf -g /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam
```
######Lumpy
```
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD10.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD10.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD20.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD20.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD30.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD30.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD40.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD40.pesr.bedpe

samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta /mnt/EXT/Mills-scr
atch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-sc
ratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference.flux/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/BamFiles/comp.homo.RD50.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sorted.histo,mean:499.645945946,stdev:101.066538064,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.sr.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/Lumpy/comp.homo.RD50.pesr.bedpe
```
######Pindel
```
pindel2vcf -P homo.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD10.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD20.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD30.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD40.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf

pindel2vcf -P homo.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
pindel2vcf -P homo.RD50.sorted.bam.pindel -r /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -R human_g1k_v37 -d 20150515 -v /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.homo/pbs/Pindel/Pindel.homo.RD.10.vcf
```
-->
