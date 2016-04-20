##we compared *SVelter*, *Delly*, *Lumpy*, *Pindel* and *SMuFin* based on the somatic events SMuFin have implemented in their paper
(`Moncunill V, Gonzalez S, Beà S, et al. Comprehensive characterization of complex structural variations in cancer by directly comparing genome sequence reads[J]. Nature biotechnology, 2014, 32(11): 1106-1112.`)

In the first batch comparison, we run *SVelter*, *Delly*, *Lumpy*, *Pindel* and *SMuFin* on simulated events on Chr22 at Read Depth 30, as provided by *SMuFin*. 

we provided here the scripts we used to to apply each algorithm and interprete the results:

####Apply SVelter:
```
SVelter.py Setup --reference hg19.fa --exclude Exclude.hg19.bed --copyneutral CN2.hg19.bed --svelter-path /scratch/remills_flux/xuefzhao/svelter/ --ref-index /scratch/remills_flux/xuefzhao/svelter/Index.Reference/hg19/
SVelter.py --sample chr22_insilico_Normal.sorted.bam --workdir / --null-model C
SVelter.py --sample chr22_insilico_Tumor.sorted.bam --workdir / --null-model C
```

####Apply Delly:
```
delly -t DEL -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Normal.DEL.vcf -g hg19.fa chr22_insilico_Normal.sorted.bam
delly -t DUP -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Normal.DEL.vcf -g hg19.fa chr22_insilico_Normal.sorted.bam
delly -t INV -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Normal.DEL.vcf -g hg19.fa chr22_insilico_Normal.sorted.bam
delly -t TRA -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Normal.DEL.vcf -g hg19.fa chr22_insilico_Normal.sorted.bam

delly -t DEL -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Tumor.DEL.vcf -g hg19.fa chr22_insilico_Tumor.sorted.bam
delly -t DUP -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Tumor.DEL.vcf -g hg19.fa chr22_insilico_Tumor.sorted.bam
delly -t INV -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Tumor.DEL.vcf -g hg19.fa chr22_insilico_Tumor.sorted.bam
delly -t TRA -s 5 -x human.hg19.excl.tsv -o ../Delly/chr22_insilico_Tumor.DEL.vcf -g hg19.fa chr22_insilico_Tumor.sorted.bam
```

####Apply Lumpy :
```
samtools view chr22_insilico_Normal.sorted.bam | ../Lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > ../Lumpy/chr22_insilico_Normal.sorted.um.fq
bwa bwasw -H -t 20 hg19.fa ../Lumpy/chr22_insilico_Normal.sorted.um.fq | samtools view -Sb ->../Lumpy/chr22_insilico_Normal.sorted.um.bam
samtools sort ../Lumpy/chr22_insilico_Normal.sorted.um.bam ../Lumpy/chr22_insilico_Normal.sorted.um.sorted
samtools index ../Lumpy/chr22_insilico_Normal.sorted.um.sorted.bam
samtools view chr22_insilico_Normal.sorted.bam | tail -n+100000 | ../Lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o ../Lumpy/chr22_insilico_Normal.sorted.histo
#mean:499.2397	stdev:23.990148893
lumpy -mw 4 -tt 0.0 -x Exclude.hg19.bed -pe bam_file:chr22_insilico_Normal.sorted.bam,histo_file:../Lumpy/chr22_insilico_Normal.sorted.histo,mean:499.2397,stdev:23.990148893,read_length:80,min_non_overlap:80,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:../Lumpy/chr22_insilico_Normal.sorted.um.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > ../Lumpy/chr22_insilico_Normal.sorted.bedpe

samtools view chr22_insilico_Tumor.sorted.bam | ../Lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > ../Lumpy/chr22_insilico_Tumor.sorted.um.fq
bwa bwasw -H -t 20 hg19.fa ../Lumpy/chr22_insilico_Tumor.sorted.um.fq | samtools view -Sb ->../Lumpy/chr22_insilico_Tumor.sorted.um.bam
samtools sort ../Lumpy/chr22_insilico_Tumor.sorted.um.bam ../Lumpy/chr22_insilico_Tumor.sorted.um.sorted
samtools index ../Lumpy/chr22_insilico_Tumor.sorted.um.sorted.bam
samtools view chr22_insilico_Tumor.sorted.bam | tail -n+100000 | ../Lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o ../Lumpy/chr22_insilico_Tumor.sorted.histo
#mean:499.2397	stdev:23.990148893
lumpy -mw 4 -tt 0.0 -x Exclude.hg19.bed -pe bam_file:chr22_insilico_Tumor.sorted.bam,histo_file:../Lumpy/chr22_insilico_Tumor.sorted.histo,mean:499.5817,stdev:21.3189241077,read_length:80,min_non_overlap:80,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:../Lumpy/chr22_insilico_Tumor.sorted.um.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > ../Lumpy/chr22_insilico_Tumor.sorted.bedpe
```

####Apply Pindel
```
pindel -f hg19.fa -i ../chr22_insilico_Normal.sorted.bam.config.txt -c ALL -o chr22_insilico_Normal.sorted.bam.pindel
pindel2vcf -P chr22_insilico_Normal.sorted.bam -r hg19.fa -R hg19 -d 20151102 -v ../chr22_insilico_Normal.sorted.bam.vcf
pindel -f hg19.fa -i ../chr22_insilico_Tumor.sorted.bam.config.txt -c ALL -o chr22_insilico_Tumor.sorted.bam.pindel
pindel2vcf -P chr22_insilico_Tumor.sorted.bam.pindel -r hg19.fa -R hg19 -d 20151102 -v ../Pindel/chr22_insilico_Tumor.sorted.bam.vcf
```

####Apply SMuFin
```
mpirun --np 16 ./SMuFin --ref ref_genome/hg19.fa --normal_fastq_1 normal_fastqs_1.txt --normal_fastq_2 normal_fastqs_2.txt --tumor_fastq_1 tumor_fastqs_1.txt --tumor_fastq_2 tumor_fastqs_2.txt --patient_id chr22_insilico --cpus_per_node 16
```



###Detailed Code Used:
<!---

####Apply SVelter:
```
SVelter.py Setup --reference /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa --exclude /scratch/remills_flux/xuefzhao/svelter/Support/Exclude.hg19.bed --copyneutral /scratch/remills_flux/xuefzhao/svelter/Support/CN2.hg19.bed --svelter-path /scratch/remills_flux/xuefzhao/svelter/ --ref-index /scratch/remills_flux/xuefzhao/svelter/Index.Reference/hg19/

SVelter.py --sample /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Normal.sorted.bam --workdir /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/ --null-model C

SVelter.py --sample /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Tumor.sorted.bam --workdir /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/ --null-model C
```

####Apply Delly:
```
delly -t DEL -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Normal.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Normal.sorted.bam
delly -t DUP -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Normal.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Normal.sorted.bam
delly -t INV -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Normal.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Normal.sorted.bam
delly -t TRA -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Normal.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Normal.sorted.bam

delly -t DEL -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Tumor.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Tumor.sorted.bam
delly -t DUP -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Tumor.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Tumor.sorted.bam
delly -t INV -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Tumor.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Tumor.sorted.bam
delly -t TRA -s 5 -x /scratch/remills_flux/xuefzhao/reference/Delly.reference/human.hg19.excl.tsv -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Delly/chr22_insilico_Tumor.DEL.vcf -g /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/chr22_insilico_Tumor.sorted.bam
```

####Apply Lumpy :
```
samtools view /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/alignment/chr22_insilico_Normal.sorted.bam | /nfs/remills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.um.fq
bwa bwasw -H -t 20 /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.um.fq | samtools view -Sb ->/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.um.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.um.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.um.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.um.sorted.bam
samtools view /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/alignment/chr22_insilico_Normal.sorted.bam | tail -n+100000 | /nfs/remills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.histo
#mean:499.2397	stdev:23.990148893
lumpy -mw 4 -tt 0.0 -x /scratch/remills_flux/xuefzhao/svelter/Support/Exclude.hg19.bed -pe bam_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/alignment/chr22_insilico_Normal.sorted.bam,histo_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.histo,mean:499.2397,stdev:23.990148893,read_length:80,min_non_overlap:80,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.um.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Normal.sorted.bedpe

samtools view /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/alignment/chr22_insilico_Tumor.sorted.bam | /nfs/remills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.fq
bwa bwasw -H -t 20 /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.fq | samtools view -Sb ->/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.bam
samtools sort /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.bam /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.sorted
samtools index /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.sorted.bam
samtools view /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/alignment/chr22_insilico_Tumor.sorted.bam | tail -n+100000 | /nfs/remills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.histo
#mean:499.2397	stdev:23.990148893
lumpy -mw 4 -tt 0.0 -x /scratch/remills_flux/xuefzhao/svelter/Support/Exclude.hg19.bed -pe bam_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/alignment/chr22_insilico_Tumor.sorted.bam,histo_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.histo,mean:499.5817,stdev:21.3189241077,read_length:80,min_non_overlap:80,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.bedpe
```

#Apply Pindel
```
pindel -f /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa -i /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/pbs/Pindel/chr22_insilico_Normal.sorted.bam.config.txt -c ALL -o chr22_insilico_Normal.sorted.bam.pindel
pindel2vcf -P chr22_insilico_Normal.sorted.bam -r /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa -R hg19 -d 20151102 -v /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/pbs/Pindel/chr22_insilico_Normal.sorted.bam.vcf

pindel -f /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa -i /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/pbs/Pindel/chr22_insilico_Tumor.sorted.bam.config.txt -c ALL -o chr22_insilico_Tumor.sorted.bam.pindel
pindel2vcf -P chr22_insilico_Tumor.sorted.bam.pindel -r /scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa -R hg19 -d 20151102 -v /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/pbs/Pindel/chr22_insilico_Tumor.sorted.bam.vcf
```

-->
