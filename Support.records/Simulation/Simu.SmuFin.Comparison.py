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
RM_Normal_From_Tumor(ppre+'SVelter/')
bed_size_control(ppre+'SVelter/')
make_alt_folder(ppre+'SVelter/',ppre+'alt_beds/')

#Process Delly:
vcf_to_bed(ppre+'Delly/')
order_bed(ppre+'Delly/')
combind_beds(ppre+'Delly/')
RM_Normal_From_Tumor(ppre+'Delly/')
bed_size_control(ppre+'Delly/')
make_alt_folder(ppre+'Delly/',ppre+'alt_beds/')

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
for k1 in os.listdir(ppre+'Pindel_Old/vcf/'):
	if k1.split('.')[-1]=='vcf' and k1.split('.')[-2]=='LargerThan5':
		if not k1.split('.')[1] in data_hash.keys():
			data_hash[k1.split('.')[1]]={}
		if not k1.split('.')[0].split('_')[-1] in data_hash[k1.split('.')[1]].keys():
			data_hash[k1.split('.')[1]][k1.split('.')[0].split('_')[-1]]=[]
		data_hash[k1.split('.')[1]][k1.split('.')[0].split('_')[-1]].append(k1)

for k1 in data_hash.keys():
	for k2 in data_hash[k1].keys():
		fout=ppre+'Pindel_Old/'+k1+'_'+k2+'.vcf'
		os.system(r'''cat %s > %s'''%(' '.join([ppre+'Pindel_Old/vcf/'+x for x in data_hash[k1][k2]]),fout))

vcf_to_bed(ppre+'Pindel_Old/')
order_bed(ppre+'Pindel_Old/')
combind_beds(ppre+'Pindel_Old/')
RM_Normal_From_Tumor(ppre+'Pindel_Old/')
bed_size_control(ppre+'Pindel_Old/')
make_alt_folder(ppre+'Pindel_Old/',ppre+'alt_beds/')

#process erds:
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


Produce.Pseudo.ROC.stats.py --path_ref ref_files --path_in alt_beds --appdix .min100.max1000000000.bed
Produce.Pseudo.ROC.stats.py --path_ref ref_files --path_in alt_beds --appdix .min5.max500.bed



