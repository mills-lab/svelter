import os
def path_modify(path):
	if not path[-1]=='/':
		path+='/'
	return path

def Vali_files_readin(Vali_Folder):
	out=[]
	Vali_Folder=path_modify(Vali_Folder)
	for k1 in os.listdir(Vali_Folder):
		if k1.split('.')[-1]=='PacVal':
			out.append(Vali_Folder+k1)
	return out

def Vali_file_analysis(Vali_file,out,Score_Cff):
	fin=open(Vali_file)
	for line in fin:
		pin=line.strip().split()
		if float(pin[-1])>Score_Cff:
			if not pin[0] in out.keys():
				out[pin[0]]={}
			if not int(pin[1]) in out[pin[0]].keys():
				out[pin[0]][int(pin[1])]={}
			if not int(pin[2]) in out[pin[0]][int(pin[1])].keys():
				out[pin[0]][int(pin[1])][int(pin[2])]=float(pin[-1])
			else:
				out[pin[0]][int(pin[1])][int(pin[2])]=max([float(pin[-1]),out[pin[0]][int(pin[1])][int(pin[2])]])
	return out

def order_vali_hash(Vali_data):
	out={}
	for k1 in Vali_data.keys():
		out[k1]=[]
		for k2 in sorted(Vali_data[k1].keys()):
			for k3 in sorted(Vali_data[k1][k2].keys()):
				out[k1].append([k2,k3,Vali_data[k1][k2][k3]])
	return out

def unify_vali_hash(Vali_list):
	out={}
	for k1 in Vali_list.keys():
		out[k1]=[]
		for k2 in Vali_list[k1]:
			if out[k1]==[]:
				out[k1].append(k2)
			else:
				if k2[0]<out[k1][-1][-2]:
					out[k1][-1]+=k2
				else:
					out[k1].append(k2)
	return out

def compare_Reciprocal(list1,list2):
	if list2[0]>list1[1]: 
		return 0
	elif list2[1]<list1[0]:
		return 0
	else:
		return float(sorted(list1[:2]+list2[:2])[2]-sorted(list1[:2]+list2[:2])[1])/float(max([list2[1]-list2[0],list1[1]-list1[0]]))

def compare_list(list):
	rec=len(list)/3
	sub=[]
	for i in range(rec):
		sub.append(list[(i*3):(i*3+3)])
	sub2=[]
	test=[]
	test2=[]
	for a in range(len(sub)):
		for b in range(len(sub)):
			if b>a:
				i=sub[a]
				j=sub[b]
				if compare_Reciprocal(i,j)>0.5:
					if i[-1]>j[-1]:
						if not i in test2:
							test2.append(i)
							sub2+=i
					else:
						if not j in test2:
							test2.append(j)
							sub2+=j
					test.append(i)
					test.append(j)
	for i in sub:
		if not i in test:
			sub2+=i
	return sub2

def recursion_compare(list_old):
	while True:
		list_new=compare_list(list_old)
		if list_new==list_old: break
		else:
			list_old=list_new
	return list_new

def subsplit_list(compared_list,length):
	rec=len(compared_list)/length
	sub=[]
	for i in range(rec):
		sub.append(compared_list[(i*length):(i*length+length)])
	return sub

ppre='/scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/comparison_diff_algorithm/'
Score_Cff=0.1
Vali_Folder=ppre+'Extra_bed/'
Vali_files=Vali_files_readin(Vali_Folder)
Vali_data={}
for k1 in Vali_files:
	Vali_data=Vali_file_analysis(k1,Vali_data,Score_Cff)

Vali_list=order_vali_hash(Vali_data)
Vali_Uni=unify_vali_hash(Vali_list)
fo=open(ppre+'Validated.del.bed','w')
for k1 in Vali_Uni:
	for k2 in Vali_Uni[k1]:
		if len(k2)==3:
			print >>fo, ' '.join([str(i) for i in [k1]+k2])
		else:
			compared_list=recursion_compare(k2)
			out_list=subsplit_list(compared_list,3)
			for x in out_list:
				print >>fo, ' '.join([str(i) for i in [k1]+x])

fo.close()



#linux
cp ref_bed/Personalis_1000_Genomes_deduplicated_deletions.Mappable.min100.max10000000000.bed ref_NIST_and_VaLoR/
cp Validated.del.bed ref_NIST_and_VaLoR/
cat Personalis_1000_Genomes_deduplicated_deletions.Mappable.min100.max10000000000.bed Validated.del.bed >NIST_VaLoR_dels.bed
SV.Simple.Output.Process.py Mappable-Control --ref-prefix /scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/reference_SVelter/genome.Mappable.bed --input NIST_VaLoR_dels.bed

ln -s NIST_VaLoR_dels.Mappable.min100.max1000000000.bed   Delly.NA12878_S1_bam_DEL.DEL.Mappable.min100.max1000000000.bed
ln -s NIST_VaLoR_dels.Mappable.min100.max1000000000.bed   erds.NA12878_S1.DEL.Mappable.min100.max1000000000.bed
ln -s NIST_VaLoR_dels.Mappable.min100.max1000000000.bed   Lumpy.NA12878_S1.DEL.Mappable.min100.max1000000000.bed
ln -s NIST_VaLoR_dels.Mappable.min100.max1000000000.bed   Pindel.NA12878_S1.DEL.Mappable.min100.max1000000000.bed
ln -s NIST_VaLoR_dels.Mappable.min100.max1000000000.bed   Svelter.NA12878_S1.DEL.Mappable.min100.max1000000000.bed  
#python
prefix1='SV.Simple.Output.Process.py Mappable-Control --ref-prefix /scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/reference_SVelter/genome.Mappable.bed --input '
prefix2='SV.Simple.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input'
import os
ppre='/scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/comparison_diff_algorithm/alt_bed/'
for k1 in os.listdir(ppre):
	if k1.split('.')[-1]=='bed' and not 'Mappa' in k1:
		os.system(r'''%s %s'''%(prefix1,ppre+k1))

for k1 in os.listdir(ppre):
	if k1.split('.')[-1]=='bed' and 'Mappa' in k1:
		os.system(r'''%s %s'''%(prefix2,ppre+k1))

ln -s Personalis_1000_Genomes_deduplicated_deletions.Mappable.min100.max10000000000.bed	Delly.NA12878_S1_bam_DEL.DEL.Mappable.min100.max1000000000.bed
ln -s Personalis_1000_Genomes_deduplicated_deletions.Mappable.min100.max10000000000.bed	erds.NA12878_S1.DEL.Mappable.min100.max1000000000.bed
ln -s Personalis_1000_Genomes_deduplicated_deletions.Mappable.min100.max10000000000.bed	Lumpy.NA12878_S1.DEL.Mappable.min100.max1000000000.bed
ln -s Personalis_1000_Genomes_deduplicated_deletions.Mappable.min100.max10000000000.bed	Pindel.NA12878_S1.DEL.Mappable.min100.max1000000000.bed
ln -s Personalis_1000_Genomes_deduplicated_deletions.Mappable.min100.max10000000000.bed	Svelter.NA12878_S1.DEL.Mappable.min100.max1000000000.bed


Produce.Pseudo.ROC.stats.py --path_ref ref_bed/ --path_in alt_bed/ --appdix .Mappable.min100.max1000000000.bed --RO_cff 0.5
Produce.Pseudo.ROC.stats.py --path_ref ref_NIST_and_VaLoR/ --path_in alt_bed/ --appdix .Mappable.min100.max1000000000.bed --RO_cff 0.5

#R
data=read.table('NIST.Pseudo.ROC.Mappable.min100.max1000000000.Stats',header=T)
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
write.table(data2,'NIST.Integrated.Pseudo.ROC.Mappable.min100.max1000000000.Stats',quote=F,col.names=T,row.names=F)

data=read.table('NIST_and_VaLoR.Pseudo.ROC.Mappable.min100.max1000000000.Stats',header=T)
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
write.table(data2,'NIST_and_VaLoR.Integrated.Pseudo.ROC.Mappable.min100.max1000000000.Stats',quote=F,col.names=T,row.names=F)



