#Compare different SV detecting algorithms
####Xuefang Zhao 2015-08-24

prepare a NA12878.chr21.bam file which contains only reads from chr21 of NA12878 with samtools commands:
```
samtools view -h NA12878_S1.bam chr21 > NA12878_S1.chr21.sam 
samtools view -h -Sb NA12878_S1.chr21.sam > NA12878_S1.chr21.bam
samtools sort NA12878_S1.chr21.bam NA12878_S1.chr21.sorted 
samtools index NA12878_S1.chr21.sorted.bam
```

##Run SVelter
```
Time.Tracker.py
SVelter.py index --workdir /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/ --copyneutral /mnt/EXT/Mills-scratch2/Xuefang/SVeltel-git/svelter-master/Support/CN2.hg19.bed --exclude /mnt/EXT/Mills-scratch2/Xuefang/SVeltel-git/svelter-master/Support/Exclude.hg19.bed --Reference /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa --svelter-path /mnt/EXT/Mills-scratch2/Xuefang/SVeltel-git/svelter-master
SVelter.py --workdir /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/ --sample /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.chr21.sorted.bam
Time.Tracker.py
```

##Run Delly 
```
delly -t DEL -s 10 -x /scratch/remills_flux/xuefzhao/NA12878.NGS/Delly.reference/human.hg19.excl.tsv -o 
/scratch/remills_flux/xuefzhao/NA12878.NGS/Delly/Delly.DEL.NA12878_S1.test.DEL.vcf -g /scratch/remills_flux/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam
delly -t DUP -s 10 -x /scratch/remills_flux/xuefzhao/NA12878.NGS/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/NA12878.NGS/Delly/Delly.DUP.NA12878_S1.test.DUP.vcf -g /scratch/remills_flux/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam
delly -t INV -s 10 -x /scratch/remills_flux/xuefzhao/NA12878.NGS/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/NA12878.NGS/Delly/Delly.INV.NA12878_S1.test.INV.vcf -g /scratch/remills_flux/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam
delly -t TRA -s 10 -x /scratch/remills_flux/xuefzhao/NA12878.NGS/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/NA12878.NGS/Delly/Delly.TRA.NA12878_S1.test.TRA.vcf -g /scratch/remills_flux/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam
```

#Run Lumpy
```
samtools view /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/NA12878_S1.test.sorted.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/NA12878_S1.test.sorted.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/NA12878_S1.test.sorted.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/NA12878_S1.test.sorted.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/NA12878_S1.test.sorted.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/NA12878_S1.test.sorted.bam.sr.sorted.bamsamtools view /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/Lumpy1.NA12878_S1.test.sorted.histo
samtools view /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/Lumpy1.NA12878_S1.test.sorted.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/reference/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.test.sorted.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/Lumpy1.NA12878_S1.test.sorted.histo,mean:299.189132393,stdev:103.668379205,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/header.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Lumpy/NA12878_S1.pesr.bedpe
```
#Run Pindel
```
pindel -f /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/reference/genome.fa -i /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/pbs/Pindel/NA12878_S1.bam.config.txt -c ALL -o NA12878_S1.bam.pindel
pindel2vcf -P NA12878_S1.test.sorted.bam.pindel -r /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/reference/genome.fa -R hg19 -d 20150816 -v /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/pbs/Pindel/NA12878_S1.test.Pindel.vcf
```
