def RM_Normal_From_Tumor(path):
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bed' and 'Tumor' in k1:
			if os.path.isfile(path+k1.replace('Tumor','Normal')):
				os.system(r'''Extract.Normal.From.Tumor.py --normal %s --tumor %s --output %s'''%(path+k1.replace('Tumor','Normal'),path+k1,path+k1.replace('Tumor','Somatic')))
			else:
				os.system(r'''cp %s %s'''%(path+k1,path+k1.replace('Tumor','Somatic')))

#SVelter:
prefix='SV.Simple.Output.Process.py vcf-to-bed --input'
import os
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='vcf':
		os.system(r'''%s %s'''%(prefix,k1))

prefix1='SV.Simple.Output.Process.py Size-Control --min-size 5 --max-size 500 --input'
prefix2='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='bed':
		if not 'min' in k1:
			os.system(r'''%s %s'''%(prefix2,'./'+k1))


bed_hash={}
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='bed':
		if not k1.split('.')[2] in bed_hash.keys():
			bed_hash[k1.split('.')[2]]={}
		if not k1.split('.')[3] in bed_hash[k1.split('.')[2]].keys():
			bed_hash[k1.split('.')[2]][k1.split('.')[3]]={}
		if not k1.split('.')[-2] in bed_hash[k1.split('.')[2]][k1.split('.')[3]].keys():
			bed_hash[k1.split('.')[2]][k1.split('.')[3]][k1.split('.')[-2]]=[]
		bed_hash[k1.split('.')[2]][k1.split('.')[3]][k1.split('.')[-2]].append(k1)

for k1 in bed_hash.keys():
	for k2 in bed_hash[k1].keys():
		for k3 in bed_hash[k1][k2].keys():
			print [k1,k2,k3,bed_hash[k1][k2][k3][0],'SmuFin.Simu.'+k1+'.'+k2+'.RG.'+k3+'.bed']
			os.system(r'''cat %s > %s'''%('SmuFin.Simu.'+k1+'.'+k2+'.RG.*.'+k3+'.bed','SmuFin.Simu.'+k1+'.'+k2+'.RG.'+k3+'.bed'))


#Delly:
prefix='SV.Simple.Output.Process.py vcf-to-bed --input'
import os
for k1 in os.listdir('./Delly/'):
	if k1.split('.')[-1]=='vcf':
		os.system(r'''%s %s'''%(prefix,'./Delly/'+k1))

prefix1='SV.Simple.Output.Process.py Size-Control --min-size 5 --max-size 500 --input'
prefix2='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
for k1 in os.listdir('./Delly/'):
	if k1.split('.')[-1]=='bed':
		if not 'min' in k1:
			os.system(r'''%s %s'''%(prefix2,'./Delly/'+k1))


#Pindel:
import os
prefix='vcf.size.filter.py --size 100 -i'
all_vcfs=os.listdir('./')
for k1 in all_vcfs:
	if k1.split('.')[-1]=='vcf':
		os.system(r'''%s %s'''%(prefix,k1))

cat Integrate.SmuFin.Simu.Normal.40.RG.chr*.DEL.bed >Integrate.SmuFin.Simu.Normal.40.RG.DEL.bed
cat Integrate.SmuFin.Simu.Normal.40.RG.chr*.DUP_TANDEM.bed >Integrate.SmuFin.Simu.Normal.40.RG.DUP_TANDEM.bed
cat Integrate.SmuFin.Simu.Normal.40.RG.chr*.INV.bed >Integrate.SmuFin.Simu.Normal.40.RG.INV.bed
cat Integrate.SmuFin.Simu.Tumor.40.RG.chr*.DEL.bed >Integrate.SmuFin.Simu.Tumor.40.RG.DEL.bed
cat Integrate.SmuFin.Simu.Tumor.40.RG.chr*.DUP_TANDEM.bed >Integrate.SmuFin.Simu.Tumor.40.RG.DUP_TANDEM.bed
cat Integrate.SmuFin.Simu.Tumor.40.RG.chr*.INV.bed >Integrate.SmuFin.Simu.Tumor.40.RG.INV.bed

cat Integrate.SmuFin.Simu.Normal.50.RG.chr*.DEL.bed >Integrate.SmuFin.Simu.Normal.50.RG.DEL.bed
cat Integrate.SmuFin.Simu.Normal.50.RG.chr*.DUP_TANDEM.bed >Integrate.SmuFin.Simu.Normal.50.RG.DUP_TANDEM.bed
cat Integrate.SmuFin.Simu.Normal.50.RG.chr*.INV.bed >Integrate.SmuFin.Simu.Normal.50.RG.INV.bed
cat Integrate.SmuFin.Simu.Tumor.50.RG.chr*.DEL.bed >Integrate.SmuFin.Simu.Tumor.50.RG.DEL.bed
cat Integrate.SmuFin.Simu.Tumor.50.RG.chr*.DUP_TANDEM.bed >Integrate.SmuFin.Simu.Tumor.50.RG.DUP_TANDEM.bed
cat Integrate.SmuFin.Simu.Tumor.50.RG.chr*.INV.bed >Integrate.SmuFin.Simu.Tumor.50.RG.INV.bed

for k1 in os.listdir('./'):
    if k1.split('.')[0]=='Integrate':
            os.system(r'''mv %s %s'''%(k1,'_'.join(k1.split('.')[:5])+'.'+'.'.join(k1.split('.')[6:]))) 

prefix='SV.Simple.Output.Process.py vcf-to-bed --input'
for k1 in os.listdir('./'):
	if not k1=='vcf' or k1=='sub_chromo':
		if k1.split('.')[-1]=='vcf' and k1.split('.')[-2]=='LargerThan100':
			os.system(r'''%s %s'''%(prefix,k1))

prefix1='SV.Simple.Output.Process.py Size-Control --min-size 5 --max-size 500 --input'
prefix2='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='bed':
		if not 'min' in k1:
			os.system(r'''%s %s'''%(prefix2,k1))




#Delly_RP
import os
for k1 in os.listdir('./'):
	if k1.split('_')[0]=='SmuFin':
		os.system(r'''mv %s %s'''%(k1,'DellyRP_'+k1))

prefix='SV.Simple.Output.Process.py vcf-to-bed --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='vcf':
		os.system(r'''%s %s'''%(prefix,k1))

prefix1='SV.Simple.Output.Process.py Size-Control --min-size 5 --max-size 500 --input'
prefix2='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='bed':
		if not 'min' in k1:
			os.system(r'''%s %s'''%(prefix2,k1))


#Lumpy
import os
ref='/scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa'
prefix='SV.Simple.Output.Process.py bedpe-to-bed --ref '+ref+' --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='bedpe':
		os.system(r'''%s %s'''%(prefix,k1))


#SmuFin:
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr16_RD10
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr16_RD20
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr16_RD30
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr16_RD40
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr16_RD50
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr17_RD10
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr17_RD20
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr17_RD30
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr17_RD40
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr17_RD50
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr18_RD10
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr18_RD20
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr18_RD30
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr18_RD40
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr18_RD50
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr19_RD10
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr19_RD20
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr19_RD30
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr19_RD40
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr19_RD50
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr20_RD10
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr20_RD20
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr20_RD30
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr20_RD40
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr21_RD10
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr21_RD20
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr21_RD30
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr21_RD40
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr21_RD50
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr22_RD10
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr22_RD20
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr22_RD30
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr22_RD40
Smufin.Result.process.py large-SVs --input somatic_large_SVs.SMuFin_chr22_RD50

Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr16_RD10
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr16_RD20
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr16_RD30
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr16_RD40
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr16_RD50
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr17_RD10
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr17_RD20
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr17_RD30
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr17_RD40
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr17_RD50
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr18_RD10
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr18_RD20
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr18_RD30
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr18_RD40
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr18_RD50
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr19_RD10
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr19_RD20
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr19_RD30
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr19_RD40
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr19_RD50
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr20_RD10
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr20_RD20
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr20_RD30
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr20_RD40
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr21_RD10
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr21_RD20
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr21_RD30
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr21_RD40
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr21_RD50
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr22_RD10
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr22_RD20
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr22_RD30
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr22_RD40
Smufin.Result.process.py small-SVs --input  somatic_small_SVs.SMuFin_chr22_RD50

cat somatic_large_SVs.SMuFin_chr16_RD10.bed somatic_small_SVs.SMuFin_chr16_RD10.bed > SMuFin_chr16_RD10.bed
cat somatic_large_SVs.SMuFin_chr16_RD20.bed somatic_small_SVs.SMuFin_chr16_RD20.bed > SMuFin_chr16_RD20.bed
cat somatic_large_SVs.SMuFin_chr16_RD30.bed somatic_small_SVs.SMuFin_chr16_RD30.bed > SMuFin_chr16_RD30.bed
cat somatic_large_SVs.SMuFin_chr16_RD40.bed somatic_small_SVs.SMuFin_chr16_RD40.bed > SMuFin_chr16_RD40.bed
cat somatic_large_SVs.SMuFin_chr16_RD50.bed somatic_small_SVs.SMuFin_chr16_RD50.bed > SMuFin_chr16_RD50.bed
cat somatic_large_SVs.SMuFin_chr17_RD10.bed somatic_small_SVs.SMuFin_chr17_RD10.bed > SMuFin_chr17_RD10.bed
cat somatic_large_SVs.SMuFin_chr17_RD20.bed somatic_small_SVs.SMuFin_chr17_RD20.bed > SMuFin_chr17_RD20.bed
cat somatic_large_SVs.SMuFin_chr17_RD30.bed somatic_small_SVs.SMuFin_chr17_RD30.bed > SMuFin_chr17_RD30.bed
cat somatic_large_SVs.SMuFin_chr17_RD40.bed somatic_small_SVs.SMuFin_chr17_RD40.bed > SMuFin_chr17_RD40.bed
cat somatic_large_SVs.SMuFin_chr17_RD50.bed somatic_small_SVs.SMuFin_chr17_RD50.bed > SMuFin_chr17_RD50.bed
cat somatic_large_SVs.SMuFin_chr18_RD10.bed somatic_small_SVs.SMuFin_chr18_RD10.bed > SMuFin_chr18_RD10.bed
cat somatic_large_SVs.SMuFin_chr18_RD20.bed somatic_small_SVs.SMuFin_chr18_RD20.bed > SMuFin_chr18_RD20.bed
cat somatic_large_SVs.SMuFin_chr18_RD30.bed somatic_small_SVs.SMuFin_chr18_RD30.bed > SMuFin_chr18_RD30.bed
cat somatic_large_SVs.SMuFin_chr18_RD40.bed somatic_small_SVs.SMuFin_chr18_RD40.bed > SMuFin_chr18_RD40.bed
cat somatic_large_SVs.SMuFin_chr18_RD50.bed somatic_small_SVs.SMuFin_chr18_RD50.bed > SMuFin_chr18_RD50.bed
cat somatic_large_SVs.SMuFin_chr19_RD10.bed somatic_small_SVs.SMuFin_chr19_RD10.bed > SMuFin_chr19_RD10.bed
cat somatic_large_SVs.SMuFin_chr19_RD20.bed somatic_small_SVs.SMuFin_chr19_RD20.bed > SMuFin_chr19_RD20.bed
cat somatic_large_SVs.SMuFin_chr19_RD30.bed somatic_small_SVs.SMuFin_chr19_RD30.bed > SMuFin_chr19_RD30.bed
cat somatic_large_SVs.SMuFin_chr19_RD40.bed somatic_small_SVs.SMuFin_chr19_RD40.bed > SMuFin_chr19_RD40.bed
cat somatic_large_SVs.SMuFin_chr19_RD50.bed somatic_small_SVs.SMuFin_chr19_RD50.bed > SMuFin_chr19_RD50.bed
cat somatic_large_SVs.SMuFin_chr20_RD10.bed somatic_small_SVs.SMuFin_chr20_RD10.bed > SMuFin_chr20_RD10.bed
cat somatic_large_SVs.SMuFin_chr20_RD20.bed somatic_small_SVs.SMuFin_chr20_RD20.bed > SMuFin_chr20_RD20.bed
cat somatic_large_SVs.SMuFin_chr20_RD30.bed somatic_small_SVs.SMuFin_chr20_RD30.bed > SMuFin_chr20_RD30.bed
cat somatic_large_SVs.SMuFin_chr20_RD40.bed somatic_small_SVs.SMuFin_chr20_RD40.bed > SMuFin_chr20_RD40.bed
cat somatic_large_SVs.SMuFin_chr21_RD10.bed somatic_small_SVs.SMuFin_chr21_RD10.bed > SMuFin_chr21_RD10.bed
cat somatic_large_SVs.SMuFin_chr21_RD20.bed somatic_small_SVs.SMuFin_chr21_RD20.bed > SMuFin_chr21_RD20.bed
cat somatic_large_SVs.SMuFin_chr21_RD30.bed somatic_small_SVs.SMuFin_chr21_RD30.bed > SMuFin_chr21_RD30.bed
cat somatic_large_SVs.SMuFin_chr21_RD40.bed somatic_small_SVs.SMuFin_chr21_RD40.bed > SMuFin_chr21_RD40.bed
cat somatic_large_SVs.SMuFin_chr21_RD50.bed somatic_small_SVs.SMuFin_chr21_RD50.bed > SMuFin_chr21_RD50.bed
cat somatic_large_SVs.SMuFin_chr22_RD10.bed somatic_small_SVs.SMuFin_chr22_RD10.bed > SMuFin_chr22_RD10.bed
cat somatic_large_SVs.SMuFin_chr22_RD20.bed somatic_small_SVs.SMuFin_chr22_RD20.bed > SMuFin_chr22_RD20.bed
cat somatic_large_SVs.SMuFin_chr22_RD30.bed somatic_small_SVs.SMuFin_chr22_RD30.bed > SMuFin_chr22_RD30.bed
cat somatic_large_SVs.SMuFin_chr22_RD40.bed somatic_small_SVs.SMuFin_chr22_RD40.bed > SMuFin_chr22_RD40.bed
cat somatic_large_SVs.SMuFin_chr22_RD50.bed somatic_small_SVs.SMuFin_chr22_RD50.bed > SMuFin_chr22_RD50.bed
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

prefix1='SV.Simple.Output.Process.py Size-Control --min-size 5 --max-size 500 --input'
prefix2='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='bed':
		if not 'min' in k1:
			os.system(r'''%s %s'''%(prefix1,k1))
			os.system(r'''%s %s'''%(prefix2,k1))


#alt_bed:
mv Delly/*.bed alt_beds/
mv Delly_RP/*.bed alt_beds/
mv Lumpy/*.bed alt_beds/
mv SVelter/*.bed alt_beds/
mv Pindel/*.bed alt_beds/
mv erds/*.bed alt_beds/

##python
import os
def file_setup(file):
	if not os.path.isfile(file):
		fo=open(file,'w')
		fo.close()

for k1 in os.listdir('./'):
	if 'Tumor' in k1:
		file_setup(k1.replace('Tumor','Normal'))

for k1 in os.listdir('./'):
	if 'Tumor' in k1:
		os.system(r'''Extract.Normal.From.Tumor.py --tumor %s --normal %s --output %s --RO-cff 0.5'''%(k1,k1.replace('Tumor','Normal'),k1.replace('Tumor','Somatic')))	





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



 
ln -s Normal.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Normal_10_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Normal_20_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Normal_30_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Normal_40_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Normal_50_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Normal_10_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Normal_20_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Normal_30_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Normal_40_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Normal_50_DEL.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Normal_10.DEL_geno.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Normal_20.DEL_geno.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Normal_30.DEL_geno.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Normal_40.DEL_geno.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Normal_50.DEL_geno.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Normal_10_pesr.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Normal_20_pesr.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Normal_30_pesr.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Normal_40_pesr.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Normal_50_pesr.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Normal_10.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Normal_20.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Normal_30.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Normal_40.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Normal_50.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Normal.10.RG.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Normal.20.RG.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Normal.30.RG.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Normal.40.RG.DEL.min100.max1000000000.bed
ln -s Normal.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Normal.50.RG.DEL.min100.max1000000000.bed
 
ln -s Somatic.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_10_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_20_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_30_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_40_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_50_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_10_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_20_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_30_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_40_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_50_DEL.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Somatic_10.DEL_geno.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Somatic_20.DEL_geno.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Somatic_30.DEL_geno.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Somatic_40.DEL_geno.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed erds_erds_SmuFin_Simu_Somatic_50.DEL_geno.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_10_pesr.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_20_pesr.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_30_pesr.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_40_pesr.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_50_pesr.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_10.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_20.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_30.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_40.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_50.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.10.RG.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.20.RG.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.30.RG.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.40.RG.DEL.min100.max1000000000.bed
ln -s Somatic.DEL.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.50.RG.DEL.min100.max1000000000.bed
 
 
ln -s Somatic.INV.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_10_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_20_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_30_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_40_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_Delly_SmuFin_Simu_Somatic_50_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_10_INV.ALL.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_10_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_20_INV.ALL.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_20_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_30_INV.ALL.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_30_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_40_INV.ALL.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Delly_RP_DellyRP_SmuFin_Simu_Somatic_40_INV.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_10_pesr.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_20_pesr.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_30_pesr.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_40_pesr.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Lumpy_Lumpy_SmuFin_Simu_Somatic_50_pesr.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_10.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_20.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_30.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_40.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed Pindel_Integrate_SmuFin_Simu_Somatic_50.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.10.RG.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.20.RG.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.30.RG.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.40.RG.INV.min100.max1000000000.bed
ln -s Somatic.INV.min100.max1000000000.bed SVelter_SmuFin.Simu.Somatic.50.RG.INV.min100.max1000000000.bed
 
ln -s Somatic.ALL.min100.max1000000000.bed SMuFin_Somatic_SMuFin_RD10.ALL.min100.max1000000000.bed
ln -s Somatic.ALL.min100.max1000000000.bed SMuFin_Somatic_SMuFin_RD20.ALL.min100.max1000000000.bed
ln -s Somatic.ALL.min100.max1000000000.bed SMuFin_Somatic_SMuFin_RD30.ALL.min100.max1000000000.bed
ln -s Somatic.ALL.min100.max1000000000.bed SMuFin_Somatic_SMuFin_RD40.ALL.min100.max1000000000.bed
ln -s Somatic.ALL.min100.max1000000000.bed SMuFin_Somatic_SMuFin_RD50.ALL.min100.max1000000000.bed




Produce.Pseudo.ROC.stats.py  --path_ref ref_files/ --path_in alt_beds/ --appdix .min100.max1000000000.bed --RO_cff 0.5

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
