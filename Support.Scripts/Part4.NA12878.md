We apply different SV detecting algorithms: SVelter, Delly, Lumpy and Pindel on two real genomes: NA12878 and CHM1. Compare the results using our in house scripts:
# On NA12878
## Run SVelter on NA12878
```
SVelter.py Setup --workdir /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/ --reference /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa --exclude /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/Exclude.hg19.bed -copyneutral /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/CN2.hg19.bed --svelter-path /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support --ref-index /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/hg19/
SVelter.py --workdir /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/ --sample /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.bam  
```
#### Process the results:
```
SV.SVelter.VCF.process.py quality-control --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/SVelter/SVelter_NA12878_S1.vcf --score -20 
  `<keep only SVs with quality score  -20>`
SV.SVelter.VCF.process.py simple-extract --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/SVelter/SVelter_NA12878_S1.cff-20.vcf.vcf  
  `<extract only predicted simple events over -20>`
SV.Simple.Output.Process.py vcf-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/SVelter/SVelter_NA12878_S1.cff-20_simple.vcf
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/SVelter/SVelter_NA12878_S1.cff-20_simple.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/SVelter/SVelter_NA12878_S1.cff-20_simple.DEL.Mappable.bed
```

## Run Delly on NA12878
```
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Delly.reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Delly/Delly_NA12878_DEL.vcf -g /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alignment/NA12878_S1.bam  
```
#### Process the results:
```
grep -v LowQual /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Delly/Delly_NA12878_DEL.vcf > /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Delly/Delly_QC_NA12878_DEL.vcf
  `<remove reports that failed quality control>`
SV.Simple.Output.Process.py vcf-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Delly/Delly_QC_NA12878_DEL.vcf
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Delly/Delly_QC_NA12878_DEL.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/hg19/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Delly/Delly_QC_NA12878_DEL.DEL.Mappable.bed
```

## Run Lumpy on NA12878
```
samtools view /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alignment/NA12878_S1.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alignment/NA12878_S1.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/Lumpy1.NA12878_S1.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/reference/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alignment/NA12878_S1.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/Lumpy1.NA12878_S1.histo,mean:299.189132393,stdev:103.668379205,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/header.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/NA12878_S1.pesr.bedpe
```
#### Process the results:
```
SV.Simple.Output.Process.py bedpe-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/Lumpy_NA12878_pesr.bedpe --reference /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/Lumpy_NA12878_pesr.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/hg19/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Lumpy/Lumpy_NA12878_pesr.DEL.Mappable.bed
```

## Run Pindel on NA12878
```
pindel -f /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa -i /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/pbs/Pindel/NA12878_S1.bam.config.txt -c ALL -o NA12878_S1.bam.pindel
pindel2vcf -P NA12878_S1.bam.pindel -r /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/reference/genome.fa -R hg19 -d 20150816 -v /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/pbs/Pindel/NA12878_S1.test.Pindel.vcf
```
#### Process the results:
```
vcf.size.filter.py -i /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Pindel/Pindel.NA12878_S1.vcf --size 100
  `<keep only SV reports with size over 100bp>`
SV.Simple.Output.Process.py vcf-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Pindel/Pindel.NA12878_S1.LargerThan100.vcf
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Pindel/Pindel.NA12878_S1.LargerThan100.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/hg19/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/Pindel/Pindel.NA12878_S1.LargerThan100.DEL.Mappable.bed
```

## Prepare reference bed from Personalis:
```
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed
change.chromo.name.py add --input Personalis_1000_Genomes_deduplicated_deletions.bed 
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Personalis_1000_Genomes_deduplicated_deletions.2.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/hg19/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Personalis_1000_Genomes_deduplicated_deletions.2.Mappable.bed
```
#### Set symbolic links:
```
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Personalis_1000_Genomes_deduplicated_deletions.2.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Delly_QC_NA12878_DEL.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Personalis_1000_Genomes_deduplicated_deletions.2.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Lumpy_NA12878_pesr.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Personalis_1000_Genomes_deduplicated_deletions.2.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Pindel.NA12878_S1.LargerThan100.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/Personalis_1000_Genomes_deduplicated_deletions.2.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/SVelter_QC_NA12878_S1.DEL.Mappable.min100.max1000000000.bed
```
## Run Comparison algorithms:
```
Produce.Pseudo.ROC.stats.py --path_ref /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/ --path_in /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alt_bed_before_pacbio/ --appdix .Mappable.min100.max1000000000.bed
```


## Prepare reference bed from Personalis and Pacbio Validation:
```
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/DEL.REF.Personalies.and.PacValied.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/hg19/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/DEL.REF.Personalies.and.PacValied.Mappable.bed
```
#### Set symbolic links:
```
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/DEL.REF.Personalies.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/Delly_QC_NA12878_DEL.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/DEL.REF.Personalies.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/Lumpy_NA12878_pesr.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/DEL.REF.Personalies.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/Pindel.NA12878_S1.LargerThan100.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/DEL.REF.Personalies.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/SVelter_QC_NA12878_S1.DEL.Mappable.min100.max1000000000.bed
```
# Run Comparison algorithms:
```
Produce.Pseudo.ROC.stats.py --path_ref /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_before_pacbio/ --path_in /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alt_bed_before_pacbio/ --appdix .Mappable.min100.max1000000000.bed
```

## Prepare alter bed from each algorithms:
```
ln -s /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alt_bed_before_pacbio/* /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alt_bed_after_pacbio/
```

## Run Comparison algorithms:
```
Produce.Pseudo.ROC.stats.py --path_ref /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/ref_bed_after_pacbio/ --path_in /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Compare_Different_algorithms/alt_bed_after_pacbio/ --appdix .Mappable.min100.max1000000000.bed
```



# On CHM1:
##Run SVelter on CHM1
```
SVelter.py Setup --workdir /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/ --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --exclude /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/Exclude.GRCh37.bed --copyneutral /mnt/EXT/Mills-scratch2/Xuefang/svelter/Support/CN2.GRCh37.bed --svelter-path /mnt/EXT/Mills-scratch2/Xuefang/svelter/ --ref-index /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/
SVelter.py --workdir /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/BamFiles/SAMN02744161.sorted.bam --sample /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.bam  
```
#### Process the results:
```
SV.SVelter.VCF.process.py quality-control --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.SVelter.vcf --score -20
SV.SVelter.VCF.process.py simple-extract --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.SVelter.cff-20.vcf 
SV.Simple.Output.Process.py vcf-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.SVelter.cff-20_simple.vcf
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.SVelter.cff-20_simple.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.SVelter.cff-20_simple.DEL.Mappable.bed
```

## Run Delly on CHM1
```
delly -t DEL -s 10 -x /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/Delly.reference/human.hg19.excl.tsv -o /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Delly/Delly_NA12878_DEL.vcf -g /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/BamFiles/SAMN02744161.sorted.bam
```
#### Process the results:
```
grep -v LowQual /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Delly/SAMN02744161.Delly.DEL.vcf > /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Delly.QC.DEL.vcf
SV.Simple.Output.Process.py vcf-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Delly.QC.DEL.vcf
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Delly.QC.DEL.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Delly.QC.DEL.DEL.Mappable.bed
```

## Run Lumpy on CHM1
```
samtools view /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/BamFiles/SAMN02744161.sorted.bam | /mnt/EXT/Mills-data/apps/lumpy-sv/scripts/split_unmapped_to_fasta.pl -b 20 > /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.um.fq
bwa bwasw -H -t 20 /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.um.fq | samtools view -Sb -> /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.bam
samtools sort /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.bam /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.sorted
samtools index /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/NA12878_S1.bam.sr.sorted.bam
samtools view /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alignment/NA12878_S1.bam | tail -n+100000 | /mnt/EXT/Mills-data/apps/temp/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/Lumpy1.NA12878_S1.histo
/mnt/EXT/Mills-data/apps/temp/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -x /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/reference/Exclude.bed -pe bam_file:/mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alignment/NA12878_S1.bam,histo_file:/mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/Lumpy1.NA12878_S1.histo,mean:299.189132393,stdev:103.668379205,read_length:101,min_non_overlap:101,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/header.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Lumpy/NA12878_S1.pesr.bedpe

```
#### Process the results:
```
mv SAMN02744161.pesr.bedpe /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Lumpy.pesr.bedpe
SV.Simple.Output.Process.py bedpe-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Lumpy.pesr.bedpe --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Lumpy.pesr.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161.Lumpy.pesr.DEL.Mappable.bed
```

## Run Pindel on CHM1
```
pindel -f /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta -i /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/pbs/Pindel/SAMN02744161.sorted.bam.config.txt -c ALL -o /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/pbs/Pindel/NA12878_S1.bam.pindel
pindel2vcf -P /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/pbs/Pindel/NA12878_S1.bam.pindel -r /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/reference/genome.fa -R GRCh37 -d 20150816 -v /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/Pindel/NA12878_S1.test.Pindel.vcf
```
#### Process the results:
```
vcf.size.filter.py -i /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161_sorted_sorted_pindel.vcf --size 100
SV.Simple.Output.Process.py vcf-to-bed --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161_sorted_sorted_pindel.LargerThan100.vcf
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161_sorted_sorted_pindel.LargerThan100.DEL.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/genome.Mappable.bed
SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/SAMN02744161_sorted_sorted_pindel.LargerThan100.DEL.Mappable.bed
```

## Prepare reference bed from Echler:
```
wget http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/deletions.bed
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/deletion.3col.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/genome.Mappable.bed

SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/deletion.3col.Mappable.bed
```
#### Set symbolic links:
```
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/deletion.3col.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/SAMN02744161_Delly_QC_DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/deletion.3col.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/SAMN02744161_Lumpy_pesr_DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/deletion.3col.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/SAMN02744161_sorted_sorted_pindel.LargerThan100.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/deletion.3col.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/SAMN02744161_SVelter_cff-20_simple_DEL.Mappable.min100.max1000000000.bed
```
## Run Comparison algorithms:
```
Produce.Pseudo.ROC.stats.py --path_ref /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_before_pacbio/ --path_in /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_before_pacbio/ --appdix .Mappable.min100.max1000000000.bed
```

## Prepare reference bed from Echler and Pacbio Validation:
```
SV.Simple.Output.Process.py Mappable-Control --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/deletion.Echler.and.PacValied.bed --ref-prefix /mnt/EXT/Mills-scratch2/Xuefang/svelter/Index.Reference/GRCh37/genome.Mappable.bed

SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/deletion.Echler.and.PacValied.Mappable.bed
```
#### Set symbolic links:
```
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/deletion.Echler.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/SAMN02744161_Delly_QC_DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/deletion.Echler.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/SAMN02744161_Lumpy_pesr_DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/deletion.Echler.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/SAMN02744161_sorted_sorted_pindel.LargerThan100.DEL.Mappable.min100.max1000000000.bed
ln -s /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/deletion.Echler.and.PacValied.Mappable.min100.max1000000000.bed /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/SAMN02744161_SVelter_cff-20_simple_DEL.Mappable.min100.max1000000000.bed
```
## Run Comparison algorithms:
```
Produce.Pseudo.ROC.stats.py --path_ref /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/ref_bed_after_pacbio/ --path_in /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Compare_Different_algorithms/alt_bed_after_pacbio/ --appdix .Mappable.min100.max1000000000.bed
```
