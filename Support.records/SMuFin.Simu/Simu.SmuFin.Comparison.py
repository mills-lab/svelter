import os
def path_mkdir(path):
    if not os.path.isdir(path):
        os.system(r'''mkdir %s'''%(path))

def path_modify(path):
	if not path[-1]=='/':
		path+='/'
	return path

def vcf_to_bed(path):
	path=path_modify(path)
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='vcf':
			os.system(r'''SV.Simple.Output.Process.py vcf-to-bed --input %s'''%(path+k1))

def bedpe_to_bed(path):
	path=path_modify(path)
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bedpe':
			os.system(r'''SV.Simple.Output.Process.py bedpe-to-bed --ref %s --input %s'''%(ref,path+k1))

def bed_size_control(path):
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bed' and not 'min' in k1 and not k1=='bed':
			os.system(r'''SV.Simple.Output.Process.py Size-Control --min-size 5 --max-size 500 --input %s'''%(path+k1))
			os.system(r'''SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input %s'''%(path+k1))

def combind_beds(path):
	data_hash={}
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bed':
			if not 'min' in k1 and not 'ALL' in k1:
				if not k1.split('.')[0] in data_hash.keys():
					data_hash[k1.split('.')[0]]=[]
				data_hash[k1.split('.')[0]].append(k1)
	for k1 in data_hash.keys():
		fout=path+k1+'.ALL.bed'
		os.system(r'''cat %s > %s'''%(' '.join([path+i for i in data_hash[k1]]),fout))

def make_alt_folder(path_in,path_to):
	path_mkdir(path_to)
	path_in=path_modify(path_in)
	file_prefix=path_in.split('/')[-2]+'_'
	for k1 in os.listdir(path_in):
		if k1.split('.')[-1]=='bed' and 'min' in k1:
			if 'Somatic' in k1 or 'Normal' in k1:
				os.system(r'''ln -s %s %s'''%(path_in+k1,path_to+file_prefix+k1))

def make_ref_folder(path_ref, path_alt):
	ref_hash={}
	for k1 in os.listdir(path_ref):
		if 'min' in k1:
			if not k1.split('.')[0] in ref_hash.keys():
				ref_hash[k1.split('.')[0]]={}
			if not k1.split('.')[1] in ref_hash[k1.split('.')[0]].keys():
				ref_hash[k1.split('.')[0]][k1.split('.')[1]]={}
			if not k1.split('.')[-3] in ref_hash[k1.split('.')[0]][k1.split('.')[1]].keys():
				ref_hash[k1.split('.')[0]][k1.split('.')[1]][k1.split('.')[-3]]=[]
			ref_hash[k1.split('.')[0]][k1.split('.')[1]][k1.split('.')[-3]].append(k1)
	alt_hash={}
	for k1 in os.listdir(path_alt):
		if 'min' in k1:
			if ref_hash.keys()[0] in k1:
				key1=ref_hash.keys()[0]
			elif ref_hash.keys()[1] in k1:
				key1=ref_hash.keys()[1]
			if not key1 in alt_hash.keys():
				alt_hash[key1]={}
			if not k1.split('.')[1] in alt_hash[key1].keys():
				alt_hash[key1][k1.split('.')[1]]={}
			if not k1.split('.')[-3] in alt_hash[key1][k1.split('.')[1]].keys():
				alt_hash[key1][k1.split('.')[1]][k1.split('.')[-3]]=[]
			alt_hash[key1][k1.split('.')[1]][k1.split('.')[-3]].append(k1)
	for k1 in ref_hash.keys():
		if k1 in alt_hash.keys():
			for k2 in ref_hash[k1].keys():
				if k2 in alt_hash[k1].keys():
					for k3 in ref_hash[k1][k2].keys():
						if k3 in alt_hash[k1][k2].keys():
							for k4 in alt_hash[k1][k2][k3]:
								os.system(r'''ln -s %s %s'''%(path_ref+ref_hash[k1][k2][k3][0],path_ref+k4))

def remove_empty_bed(path):
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bed':
			fin=os.popen(r'''wc -l %s'''%(path+k1))
			pin=fin.readline().strip().split()
			fin.close()
			if pin[0]=='0':
				os.system(r'''rm %s'''%(path+k1))

def RM_Normal_From_Tumor(path):
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bed' and 'Tumor' in k1:
			if os.path.isfile(path+k1.replace('Tumor','Normal')):
				os.system(r'''Extract.Normal.From.Tumor.py --normal %s --tumor %s --output %s'''%(path+k1.replace('Tumor','Normal'),path+k1,path+k1.replace('Tumor','Somatic')))
			else:
				os.system(r'''cp %s %s'''%(path+k1,path+k1.replace('Tumor','Somatic')))

def order_bed(path):
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bed':
			os.system(r'''SV.Simple.Output.Process.py order-bed --input %s'''%(path+k1))

def change_Lumpy_vcf_name(path):
	for k1 in os.listdir(path):
		if 'pesr' in k1:
			k2=k1.replace('.pesr','')
			os.system(r'''mv %s %s'''%(path+k1,path+k2))

def SMufin_process(path):
	for k1 in os.listdir(path):
		if k1.split('_')[0]=='somatic':
			if k1.split('_')[1]=='large':
				os.system(r'''Smufin.Result.process.py large-SVs --input %s'''%(path+k1))
			elif k1.split('_')[1]=='small':
				os.system(r'''Smufin.Result.process.py small-SVs --input %s'''%(path+k1))
	cat_hash={}
	for k1 in os.listdir(path):
		if k1.split('.')[-1]=='bed' and not k1=='bed':
			if not k1.split('.')[1] in cat_hash.keys():
				cat_hash[k1.split('.')[1]]=[]
			cat_hash[k1.split('.')[1]].append(k1)
	for k1 in cat_hash.keys():
		fout='Somatic_'+k1+'.ALL.bed'
		os.system(r'''cat %s > %s'''%(' '.join([path+i for i in cat_hash[k1]]),path+fout))

ppre='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/'
global ref
ref='/scratch/remills_flux/xuefzhao/reference/hg19/hg19.fa'
#Process SVelter:
vcf_to_bed(ppre+'SVelter/')
order_bed(ppre+'SVelter/')
combind_beds(ppre+'SVelter/')
remove_empty_bed(ppre+'SVelter/')
RM_Normal_From_Tumor(ppre+'SVelter/')
bed_size_control(ppre+'SVelter/')
make_alt_folder(ppre+'SVelter/',ppre+'alt_beds/')

#Process Delly:
vcf_to_bed(ppre+'Delly/')
order_bed(ppre+'Delly/')
combind_beds(ppre+'Delly/')
remove_empty_bed(ppre+'Delly/')
RM_Normal_From_Tumor(ppre+'Delly/')
bed_size_control(ppre+'Delly/')
make_alt_folder(ppre+'Delly/',ppre+'alt_beds/')

#Process Delly_RP:
vcf_to_bed(ppre+'Delly_RP/')
order_bed(ppre+'Delly_RP/')
combind_beds(ppre+'Delly_RP/')
remove_empty_bed(ppre+'Delly_RP/')
RM_Normal_From_Tumor(ppre+'Delly_RP/')
bed_size_control(ppre+'Delly_RP/')
make_alt_folder(ppre+'Delly_RP/',ppre+'alt_beds/')

#Process Lumpy:
change_Lumpy_vcf_name(ppre+'Lumpy/')
bedpe_to_bed(ppre+'Lumpy/')
order_bed(ppre+'Lumpy/')
combind_beds(ppre+'Lumpy/')
RM_Normal_From_Tumor(ppre+'Lumpy/')
bed_size_control(ppre+'Lumpy/')
make_alt_folder(ppre+'Lumpy/',ppre+'alt_beds/')

data_hash={}
#Process Pindel
prefix='vcf.size.filter.py --size 5 -i'
for k1 in os.listdir(ppre+'Pindel/vcf/'):
	if k1.split('.')[-1]=='vcf' and not 'LargerThan' in k1:
		os.system(r'''%s %s'''%(prefix,ppre+'Pindel/vcf/'+k1))

for k1 in os.listdir(ppre+'Pindel/vcf/'):
	if k1.split('.')[-1]=='vcf' and k1.split('.')[-2]=='LargerThan5':
		if not k1.split('.')[1] in data_hash.keys():
			data_hash[k1.split('.')[1]]={}
		if not k1.split('.')[0].split('_')[-1] in data_hash[k1.split('.')[1]].keys():
			data_hash[k1.split('.')[1]][k1.split('.')[0].split('_')[-1]]=[]
		data_hash[k1.split('.')[1]][k1.split('.')[0].split('_')[-1]].append(k1)

for k1 in data_hash.keys():
	for k2 in data_hash[k1].keys():
		fout=ppre+'Pindel/'+k1+'_'+k2+'.vcf'
		os.system(r'''cat %s > %s'''%(' '.join([ppre+'Pindel/vcf/'+x for x in data_hash[k1][k2]]),fout))

vcf_to_bed(ppre+'Pindel/')
order_bed(ppre+'Pindel/')
combind_beds(ppre+'Pindel/')
RM_Normal_From_Tumor(ppre+'Pindel/')
bed_size_control(ppre+'Pindel/')
make_alt_folder(ppre+'Pindel/',ppre+'alt_beds/')

#process erds:
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

awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_10.DEL.bed>erds_SmuFin_Simu_Normal_10.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_10.DUP.bed>erds_SmuFin_Simu_Normal_10.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_20.DEL.bed>erds_SmuFin_Simu_Normal_20.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_20.DUP.bed>erds_SmuFin_Simu_Normal_20.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_30.DEL.bed>erds_SmuFin_Simu_Normal_30.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_30.DUP.bed>erds_SmuFin_Simu_Normal_30.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_40.DEL.bed>erds_SmuFin_Simu_Normal_40.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_40.DUP.bed>erds_SmuFin_Simu_Normal_40.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_50.DEL.bed>erds_SmuFin_Simu_Normal_50.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Normal_50.DUP.bed>erds_SmuFin_Simu_Normal_50.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_10.DEL.bed>erds_SmuFin_Simu_Tumor_10.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_10.DUP.bed>erds_SmuFin_Simu_Tumor_10.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_20.DEL.bed>erds_SmuFin_Simu_Tumor_20.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_20.DUP.bed>erds_SmuFin_Simu_Tumor_20.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_30.DEL.bed>erds_SmuFin_Simu_Tumor_30.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_30.DUP.bed>erds_SmuFin_Simu_Tumor_30.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_40.DEL.bed>erds_SmuFin_Simu_Tumor_40.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_40.DUP.bed>erds_SmuFin_Simu_Tumor_40.DUP_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_50.DEL.bed>erds_SmuFin_Simu_Tumor_50.DEL_geno.bed
awk '{print $1,$2,$3,"het\t"}' erds_SmuFin_Simu_Tumor_50.DUP.bed>erds_SmuFin_Simu_Tumor_50.DUP_geno.bed

RM_Normal_From_Tumor(ppre+'erds/')
bed_size_control(ppre+'erds/')
make_alt_folder(ppre+'erds/',ppre+'alt_beds/')



#process SMuFin:
SMufin_process(ppre+'SMuFin/')
bed_size_control(ppre+'SMuFin/')
make_alt_folder(ppre+'SMuFin/',ppre+'alt_beds/')

#Process Ref:
bed_size_control(ppre+'ref_files/')
make_ref_folder(ppre+'ref_files/',ppre+'alt_beds/')


Produce.Pseudo.ROC.stats.py --path_ref ref_files --path_in alt_beds --appdix .min100.max1000000000.bed --RO_cff 0.5






plot(c(1,5),c(0,1),xlab='',main='Somatic',bty="n",ylab='Sensitivity',xaxt='n',type='n',)
data_b=data_new[data_new[,8]=='Somatic',]
cols=c('red','blue','darkblue','darkgreen','purple','orange','black')
algorithms=c('SVelter','Delly','DellyRP','Lumpy','Pindel','erds','SMuFin')
linetype <- c(1:7) 
plotchar <- seq(18,18+7,1)
for(k1 in algorithms){
  temp=data_b[data_b[,10]==k1,]
  temp=temp[order(temp[,9]),]
  color=cols[match(k1,algorithms)]
  ltype=linetype[match(k1,algorithms)]
  ptype=plotchar[match(k1,algorithms)]
  lines(c(1:5),temp[,7],col=color,pch=ptype,lty=ltype,type='b')
}
axis(1,paste('RD',seq(10,50,by=10),sep=''),cex=0.8,at=c(1:5))
legend(c(3,5),c(0,0.6),algorithms,col=cols,lty=linetype,pch=plotchar,cex=0.7)







