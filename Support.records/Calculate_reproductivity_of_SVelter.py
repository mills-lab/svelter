#Calculate reproductivity of SVelter:
#python
import os
ppre='/scratch/remills_flux/xuefzhao/NA12878.NGS/'
key_word='hg19_rep'
for k1 in os.listdir(ppre):
	if key_word in k1:
		path1=ppre+k1+'/'
		for k2 in os.listdir(path1):
			if k2.split('.')[-1]=='vcf':
				os.system(r'''SV.Simple.Output.Process.py vcf-to-bed --input %s'''%(path1+k2))

for k1 in os.listdir(ppre):
	if key_word in k1:
		path1=ppre+k1+'/'
		all_files=os.listdir(path1)
		for k2 in all_files:
			if k2.split('.')[-1]=='bed':
				os.system(r'''mv %s %s'''%(path1+k2,path1+k2.replace('NA12878_S1.','NA12878_S1_')))


#linux
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep1/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep2/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep3/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep4/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep5/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep6/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep7/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep8/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref hg19_rep0/ --path_in hg19_rep9/ --appdix .bed --RO_cff 0.1


Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep0/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep1/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep2/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep3/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep4/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep5/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep6/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep7/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep8/ --appdix .bed --RO_cff 0.1
Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in hg19_rep9/ --appdix .bed --RO_cff 0.1
cp Integrate.stat.R hg19_rep0/
cp Integrate.stat.R hg19_rep1/
cp Integrate.stat.R hg19_rep2/
cp Integrate.stat.R hg19_rep3/
cp Integrate.stat.R hg19_rep4/
cp Integrate.stat.R hg19_rep5/
cp Integrate.stat.R hg19_rep6/
cp Integrate.stat.R hg19_rep7/
cp Integrate.stat.R hg19_rep8/
cp Integrate.stat.R hg19_rep9/


cd hg19_rep0/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep1/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep2/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep3/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep4/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep5/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep6/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep7/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep8/
Rscript Integrate.stat.R 
cd ..
cd hg19_rep9/
Rscript Integrate.stat.R 
cd ..

cat hg19_rep0/Pseudo.ROC.Integrated.Stats hg19_rep1/Pseudo.ROC.Integrated.Stats hg19_rep2/Pseudo.ROC.Integrated.Stats hg19_rep3/Pseudo.ROC.Integrated.Stats hg19_rep4/Pseudo.ROC.Integrated.Stats hg19_rep5/Pseudo.ROC.Integrated.Stats hg19_rep6/Pseudo.ROC.Integrated.Stats hg19_rep7/Pseudo.ROC.Integrated.Stats hg19_rep8/Pseudo.ROC.Integrated.Stats hg19_rep9/Pseudo.ROC.Integrated.Stats >  Pseudo.ROC.Integrated.Stats


#!/usr/bin/env R

#!R
data=read.table('Pseudo.ROC.Stats',header=T)
data=data[data[,4]>0,]
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
write.table(data2,'Pseudo.ROC.Integrated.Stats',quote=F,col.names=F,row.names=F)


java_tool='/nfs/turbo/dcmb-brainsom/technical/application/jdk1.8.0_73/bin/java -jar'
picard_tool='/nfs/turbo/dcmb-brainsom/technical/application/picard-tools-2.1.0/picard.jar AddOrReplaceReadGroups'
#input='/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/alignment/het.RD10.sorted.bam' 
#output='/mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/alignment/het.RD10.RG.sorted.bam' 
ppre='/nfs/remills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/'
pbs_path=ppre+'pbs/pbs_add_ReadGroup/'
path_mkdir(pbs_path)
bam_path=ppre+'alignment/'
for k1 in os.listdir(bam_path):
        if k1.split('.')[-1]=='bam':
                input_file=bam_path+k1
                output_file=bam_path+k1.replace('.sorted.bam','.RG.sorted.bam')
                IDs='_'.join(k1.split('/')[-1].split('.')[:2])
                fpbs=pbs_path+k1+'.pbs'
                JobToDo=k1
                write_flux_new_pbs_header(fpbs,JobToDo)
                fo=open(fpbs,'a')
                print >>fo, ' '.join([java_tool,picard_tool,'I='+input_file,'O='+output_file,'RGID='+IDs,'RGLB='+IDs,'RGPU=unit1','RGPL=illumina','RGSM=20'])
                fo.close()





