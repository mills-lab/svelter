#SVelter:
prefix='SV.Simple.Output.Process.py vcf-to-bed --input'
import os
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='vcf':
		os.system(r'''%s %s'''%(prefix,k1))

#Delly:
prefix='SV.Simple.Output.Process.py vcf-to-bed --input'
import os
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='vcf':
		os.system(r'''%s %s'''%(prefix,k1))

for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='vcf':
		SV=k1.split('.')[0].split('_')[-1]
		if not os.path.isfile(k1.replace('.vcf','.'+SV+'.bed')):
			fo=open(k1.replace('.vcf','.'+SV+'.bed'),'w')
			fo.close()

#Delly_RP
import os
for k1 in os.listdir('./'):
	if k1.split('_')[0]=='SmuFin':
		os.system(r'''mv %s %s'''%(k1,'DellyRP_'+k1))

prefix='SV.Simple.Output.Process.py vcf-to-bed --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='vcf':
		os.system(r'''%s %s'''%(prefix,k1))

for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='vcf':
		SV=k1.split('.')[0].split('_')[-1]
		if not os.path.isfile(k1.replace('.vcf','.'+SV+'.bed')):
			fo=open(k1.replace('.vcf','.'+SV+'.bed'),'w')
			fo.close()

#Lumpy
import os
ref='/scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa'
prefix='SV.Simple.Output.Process.py bedpe-to-bed --ref '+ref+' --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='bedpe':
		os.system(r'''%s %s'''%(prefix,k1))

#SmuFin:
cat somatic_small_SVs.SMuFin_chr*RD10 > SMuFin.somatic_small_SVs.RD10
cat somatic_small_SVs.SMuFin_chr*RD20 > SMuFin.somatic_small_SVs.RD20
cat somatic_small_SVs.SMuFin_chr*RD30 > SMuFin.somatic_small_SVs.RD30
cat somatic_small_SVs.SMuFin_chr*RD40 > SMuFin.somatic_small_SVs.RD40
cat somatic_small_SVs.SMuFin_chr*RD50 > SMuFin.somatic_small_SVs.RD50

cat somatic_large_SVs.SMuFin_chr*RD10 >SMuFin.somatic_large_SVs.RD10
cat somatic_large_SVs.SMuFin_chr*RD20 >SMuFin.somatic_large_SVs.RD20
cat somatic_large_SVs.SMuFin_chr*RD30 >SMuFin.somatic_large_SVs.RD30
cat somatic_large_SVs.SMuFin_chr*RD40 >SMuFin.somatic_large_SVs.RD40
cat somatic_large_SVs.SMuFin_chr*RD50 >SMuFin.somatic_large_SVs.RD50

cat SMuFin_chr*_RD10.bed > SMuFin_somatic.RD10.bed
cat SMuFin_chr*_RD20.bed > SMuFin_somatic.RD20.bed
cat SMuFin_chr*_RD30.bed > SMuFin_somatic.RD30.bed
cat SMuFin_chr*_RD40.bed > SMuFin_somatic.RD40.bed
cat SMuFin_chr*_RD50.bed > SMuFin_somatic.RD50.bed



Normal.ALL.min100.max1000000000.bed
Normal.DEL.min100.max1000000000.bed
Normal.INS.min100.max1000000000.bed
Normal.INV.min100.max1000000000.bed
Somatic.ALL.min100.max1000000000.bed
Somatic.DEL.min100.max1000000000.bed
Somatic.INS.min100.max1000000000.bed
Somatic.INV.min100.max1000000000.bed


#python
import os
for k1 in os.listdir('./'):
	if 'somatic_small_SVs' in k1 and not 'bed' in k1:
		fin=open(k1)
		fo=open(k1+'.bed','w')
		for line in fin:
			pin=line.strip().split()
			if not pin[1]=='Type':
				if int(pin[4])>50:
					print >>fo, ' '.join([pin[2],pin[3],str(int(pin[3])+int(pin[4]))])
		fin.close()
		fo.close()

import os
for k1 in os.listdir('./'):
	if 'somatic_large_SVs' in k1 and not 'bed' in k1:
		fin=open(k1)
		fo=open(k1+'.bed','w')
		for line in fin:
			pin=line.strip().split()
			if not pin[1]=='Type':
				if pin[2]==pin[4]:
					if int(pin[5])-int(pin[3])>50:
						print >>fo, ' '.join([pin[2],pin[3],pin[5]])
		fin.close()
		fo.close()



Pindel:
vcf.size.filter.py -i Pindel.homo.RD.10.vcf --size 100
vcf.size.filter.py -i Pindel.homo.RD.20.vcf --size 100
vcf.size.filter.py -i Pindel.homo.RD.30.vcf --size 100
vcf.size.filter.py -i Pindel.homo.RD.40.vcf --size 100
vcf.size.filter.py -i Pindel.homo.RD.50.vcf --size 100
mv Pindel.homo.RD.10.LargerThan100.vcf Pindel_homo_RD10.LargerThan100.vcf
mv Pindel.homo.RD.20.LargerThan100.vcf Pindel_homo_RD20.LargerThan100.vcf
mv Pindel.homo.RD.30.LargerThan100.vcf Pindel_homo_RD30.LargerThan100.vcf
mv Pindel.homo.RD.40.LargerThan100.vcf Pindel_homo_RD40.LargerThan100.vcf
mv Pindel.homo.RD.50.LargerThan100.vcf Pindel_homo_RD50.LargerThan100.vcf
SV.Simple.Output.Process.py vcf-to-bed --input Pindel_homo_RD10.LargerThan100.vcf
SV.Simple.Output.Process.py vcf-to-bed --input Pindel_homo_RD20.LargerThan100.vcf
SV.Simple.Output.Process.py vcf-to-bed --input Pindel_homo_RD30.LargerThan100.vcf
SV.Simple.Output.Process.py vcf-to-bed --input Pindel_homo_RD40.LargerThan100.vcf
SV.Simple.Output.Process.py vcf-to-bed --input Pindel_homo_RD50.LargerThan100.vcf

erds:
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.10.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Normal_10.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.20.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Normal_20.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.30.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Normal_30.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.40.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Normal_40.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.50.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Normal_50.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.10.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Tumor_10.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.20.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Tumor_20.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.30.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Tumor_30.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.40.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Tumor_40.DEL.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.50.RG.sorted.bam/20.del.events > erds_SmuFin_Simu_Tumor_50.DEL.bed

awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.10.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Normal_10.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.20.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Normal_20.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.30.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Normal_30.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.40.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Normal_40.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Normal.50.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Normal_50.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.10.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Tumor_10.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.20.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Tumor_20.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.30.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Tumor_30.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.40.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Tumor_40.DUP.bed
awk {'print $1,$2,$3'}  SmuFin.Simu.Tumor.50.RG.sorted.bam/20.dup.events > erds_SmuFin_Simu_Tumor_50.DUP.bed

#alt_bed:
##python
import os
def file_setup(file):
	if not os.path.isfile(file):
		fo=open(file,'w')
		fo.close()

for k1 in os.listdir('./'):
	if 'Tumor' in k1:
		if not os.path.
Extract.Normal.From.Tumor.py 




cp Delly_homo_RD10_DUP.DUP.bed Delly_homo_RD10_DUP_TANDEM.DUP_TANDEM.bed
cp Delly_homo_RD20_DUP.DUP.bed Delly_homo_RD20_DUP_TANDEM.DUP_TANDEM.bed
cp Delly_homo_RD30_DUP.DUP.bed Delly_homo_RD30_DUP_TANDEM.DUP_TANDEM.bed
cp Delly_homo_RD40_DUP.DUP.bed Delly_homo_RD40_DUP_TANDEM.DUP_TANDEM.bed
cp Delly_homo_RD50_DUP.DUP.bed Delly_homo_RD50_DUP_TANDEM.DUP_TANDEM.bed

awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD10.DEL.bed > erds_homo_RD10.DEL.geno.bed
awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD10.DUP.bed > erds_homo_RD10.DUP.geno.bed
mv erds_homo_RD10.DEL.geno.bed erds_homo_RD10.DEL.bed
mv erds_homo_RD10.DUP.geno.bed erds_homo_RD10.DUP.bed

awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD20.DEL.bed > erds_homo_RD20.DEL.geno.bed
awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD20.DUP.bed > erds_homo_RD20.DUP.geno.bed
mv erds_homo_RD20.DEL.geno.bed erds_homo_RD20.DEL.bed
mv erds_homo_RD20.DUP.geno.bed erds_homo_RD20.DUP.bed

awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD30.DEL.bed > erds_homo_RD30.DEL.geno.bed
awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD30.DUP.bed > erds_homo_RD30.DUP.geno.bed
mv erds_homo_RD30.DEL.geno.bed erds_homo_RD30.DEL.bed
mv erds_homo_RD30.DUP.geno.bed erds_homo_RD30.DUP.bed

awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD40.DEL.bed > erds_homo_RD40.DEL.geno.bed
awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD40.DUP.bed > erds_homo_RD40.DUP.geno.bed
mv erds_homo_RD40.DEL.geno.bed erds_homo_RD40.DEL.bed
mv erds_homo_RD40.DUP.geno.bed erds_homo_RD40.DUP.bed

awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD50.DEL.bed > erds_homo_RD50.DEL.geno.bed
awk '{print $1,$2,$3,"homo\t"}' erds_homo_RD50.DUP.bed > erds_homo_RD50.DUP.geno.bed
mv erds_homo_RD50.DEL.geno.bed erds_homo_RD50.DEL.bed
mv erds_homo_RD50.DUP.geno.bed erds_homo_RD50.DUP.bed

cp erds_homo_RD10.DUP.bed erds_homo_RD10.DUP_TANDEM.bed
cp erds_homo_RD20.DUP.bed erds_homo_RD20.DUP_TANDEM.bed
cp erds_homo_RD30.DUP.bed erds_homo_RD30.DUP_TANDEM.bed
cp erds_homo_RD40.DUP.bed erds_homo_RD40.DUP_TANDEM.bed
cp erds_homo_RD50.DUP.bed erds_homo_RD50.DUP_TANDEM.bed

cp Lumpy_homo_RD10.RG.pesr.DUP.bed Lumpy_homo_RD10.RG.pesr.DUP_TANDEM.bed
cp Lumpy_homo_RD20.RG.pesr.DUP.bed Lumpy_homo_RD20.RG.pesr.DUP_TANDEM.bed
cp Lumpy_homo_RD30.RG.pesr.DUP.bed Lumpy_homo_RD30.RG.pesr.DUP_TANDEM.bed
cp Lumpy_homo_RD40.RG.pesr.DUP.bed Lumpy_homo_RD40.RG.pesr.DUP_TANDEM.bed
cp Lumpy_homo_RD50.RG.pesr.DUP.bed Lumpy_homo_RD50.RG.pesr.DUP_TANDEM.bed

cp Pindel_homo_RD20.LargerThan100.DUP_TANDEM.bed Pindel_homo_RD20.LargerThan100.DUP.bed
cp Pindel_homo_RD10.LargerThan100.DUP_TANDEM.bed Pindel_homo_RD10.LargerThan100.DUP.bed
cp Pindel_homo_RD30.LargerThan100.DUP_TANDEM.bed Pindel_homo_RD30.LargerThan100.DUP.bed
cp Pindel_homo_RD40.LargerThan100.DUP_TANDEM.bed Pindel_homo_RD40.LargerThan100.DUP.bed
cp Pindel_homo_RD50.LargerThan100.DUP_TANDEM.bed Pindel_homo_RD50.LargerThan100.DUP.bed

#python
import os
ppre='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/Simulate.homo/comparison/alt_bed/'
ref_prefix='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/Simulate.homo/SVelter/reference_SVelter/genome.Mappable.bed'
for k1 in os.listdir(ppre):
	if k1.split('.')[-1]=='bed':
		if not 'Mappable.' in k1:
			os.system(r'''SV.Simple.Output.Process.py Mappable-Control --input %s --ref-prefix %s'''%(ppre+k1,ref_prefix))

TRA_rec='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/Simulate.homo/sv_rec/homo.homo.TRA.rec'
for k1 in os.listdir(ppre):
	if k1.split('.')[-1]=='bed' and k1.split('.')[-2]=='Mappable':
		os.system(r'''SV.Simple.Output.Process.py TRA-Control --TRA-rec %s --input %s'''%(TRA_rec,ppre+k1))


pre='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
for k1 in os.listdir(ppre):
	if k1.split('.')[-1]=='bed' and k1.split('.')[-3]=='Mappable':
		os.system(r'''%s %s'''%(pre,ppre+k1))




#ref_bed
awk '{print $1,$2,$3,"homo\t"}' homo.homo.DEL.bed > homo.DEL.geno.bed
awk '{print $1,$2,$3,"homo\t"}' homo.homo.DUP.bed > homo.DUP.geno.bed
awk '{print $1,$2,$3,"homo\t"}' homo.homo.DUP_TANDEM.bed > homo.DUP_TANDEM.geno.bed
awk '{print $1,$2,$3,"homo\t"}' homo.homo.INV.bed > homo.INV.geno.bed
mv homo.DEL.geno.bed  homo.DEL.bed
mv homo.DUP.geno.bed  homo.DUP.bed
mv homo.DUP_TANDEM.geno.bed  homo.DUP_TANDEM.bed
mv homo.INV.geno.bed  homo.INV.bed

import os
ppre='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/Simulate.homo/comparison/ref_bed/'
ref_prefix='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/Simulate.homo/SVelter/reference_SVelter/genome.Mappable.bed'
for k1 in os.listdir(ppre):
 if k1.split('.')[-1]=='bed':
  if not 'Mappable.' in k1:
   os.system(r'''SV.Simple.Output.Process.py Mappable-Control --input %s --ref-prefix %s'''%(ppre+k1,ref_prefix))

TRA_rec='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/Simulate.homo/sv_rec/homo.homo.TRA.rec'
for k1 in os.listdir(ppre):
 if k1.split('.')[-1]=='bed' and k1.split('.')[-2]=='Mappable':
  os.system(r'''SV.Simple.Output.Process.py TRA-Control --TRA-rec %s --input %s'''%(TRA_rec,ppre+k1))

pre='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
for k1 in os.listdir(ppre):
 if k1.split('.')[-1]=='bed' and k1.split('.')[-3]=='Mappable':
  os.system(r'''%s %s'''%(pre,ppre+k1))



import os
ref_hash={}
for k1 in os.listdir('./ref_files/'):
	if k1.split('.')[-1]=='bed':
		if not k1.split('.')[0] in ref_hash.keys():
			ref_hash[k1.split('.')[0]]={}
		if not k1.split('.')[1] in ref_hash[k1.split('.')[0]].keys():
			ref_hash[k1.split('.')[0]][k1.split('.')[1]]=[]
		ref_hash[k1.split('.')[0]][k1.split('.')[1]].append(k1)

alt_hash={}
alt_hash['Normal']={}
alt_hash['Somatic']={}
for k1 in os.listdir('./alt_beds/'):
	if 'Normal' in k1:
		key1='Normal'
	elif 'Tumor' in k1:
		key1='Somatic'
	key2=k1.split('.')[-2]
	if not key2 in alt_hash[key1].key1():
		alt_hash[key1][key2]=[]
	alt_hash[key1][key2].append(k1)









Produce.Pseudo.ROC.stats.py  --path_ref ref_bed/ --path_in alt_bed/ --appdix .Mappable.TRAFree.min100.max1000000000.bed --RO_cff 0.5

#R
data=read.table('Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Stats',header=T)
data2=data[1,]
rec=0
for(x in sort(unique(data$sample))){
for (y in sort(unique(data$sv))){
 rec=rec+1
 temp=data[data[,1]==x & data[,2]==y,]
 data2[rec,1]=x
 data2[rec,2]=y
 data2[rec,4]=sum(temp[,4])
 data2[rec,5]=sum(temp[,5])
 data2[rec,6]=sum(temp[,6])
 data2[rec,7]=data2[rec,4]/data2[rec,5]
 data2[rec,8]=data2[rec,4]/data2[rec,6]
}}
data2=data2[data2[,5]>0,]
data2=data2[,-3]
write.table(data2,'Simp.Homo.Pseudo.ROC.Mappable.TRAFree.min100.max1000000000.Barplot.Stats',quote=F,col.names=T)
