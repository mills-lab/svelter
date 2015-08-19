import os
#process ref files:
min_cffs=[50,100,200,300,400,500,600,700,800,900,1000,1500,2000]
path_ref='/mnt/EXT/Mills-data/xuefzhao/projects.axiom/1000Genome.Validation/Del.Bed/'
path_in='/mnt/EXT/Mills-data/xuefzhao/projects.axiom/1000Genome.Validation/vcf_files/'
script=path_in+'script/Bed.Size.Filter.py'
ref_in=[]
for k1 in os.listdir(path_ref):
    if k1.split('.')[-1]=='bed' and k1.split('.')[-2]=='Mappable':
        ref_in.append(path_ref+k1)

script=path_in+'script/Bed.Size.Filter.py'
for k1 in ref_in:
	for k2 in min_cffs:
		os.system(r'''python %s -i %s --min %d --max %s'''%(script,k1,k2,1000000000) )

#SVelter results:
path_in='/mnt/EXT/Mills-data/xuefzhao/projects.axiom/1000Genome.Validation/vcf_files/'
path1=path_ref
bed_in=[]
import os
for k1 in os.listdir(path_in):
     if k1.split('.')[-1]=='bed' and k1.split('.')[-2]=='Mappable':
             bed_in.append(path_in+k1)

script=path_in+'script/Bed.Size.Filter.py'
for k1 in bed_in:
	for k2 in min_cffs:
		os.system(r'''python %s -i %s --min %d --max %s'''%(script,k1,k2,1000000000) )

path2=path_in+'for.temp.use.20150116'
if not os.path.isdir(path2):
	os.system(r'''mkdir %s'''%(path2))

os.system(r'''mv %s %s'''%(path_in+'*.min*max*.bed',path2))
script=path_in+'script/Produce.Pseudo.ROC.stats.py'
for k1 in min_cffs:
	appdix='.Mappable.min'+str(k1)+'.max1000000000.bed'
	os.system(r'''python %s --path_ref %s --path_in %s --appdix %s'''%(script,path1,path2,appdix))

script=path_in+'script/Stat.Reorganize.py'
for k1 in os.listdir(path2):
	if k1.split('.')[-1]=='Stats':
		os.system(r'''python %s -i %s'''%(script,path2+'/'+k1))











