#!/usr/bin/python

#!python
#script usd to compare vcf files
#command='python VCF2bed.1kgp.mergedSV.py -i /mnt/EXT/Mills-data/xuefzhao/projects.axiom/NA12878.Validation/ALL.wgs.mergedSV.v3.20130502.svs.genotypes.NA12878.vcf'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:',)
dict_opts=dict(opts)
file_in_1=dict_opts['-i']

def end_cordi_calcu(pin):
	#info=pin
	chromo=pin[0]
	start=int(pin[1])
	end=0
	for x in pin[7].split(';'):
		if 'END' in x.split('='):
			end=int(x.split('=')[1])
	if not end==0:
		return [chromo,start,end]
	else:
		return 'Error'

def genotype_report(pin):
	gt_index=pin[8].split(':').index('GT')
	geno=pin[9].split(':')[gt_index]
	if '/' in geno:
		return 	[int(i) for i in geno.split('/')]
	elif '|' in geno:
		return 	[int(i) for i in geno.split('|')]

def genotype_calcu(pin):
	gt_index=pin[8].split(':').index('GT')
	geno=pin[9].split(':')[gt_index]
	out=0
	if '/' in geno:
		if geno.split('/')==['0','0']:
			out='Error'
		elif sorted(geno.split('/'))==['0','1']:
			out='het'
		elif sorted(geno.split('/'))==['1','1']:
			out='homo'
	elif '|' in geno:
		if geno.split('|')==['0','0']:
			out='Error'
		elif sorted(geno.split('|'))==['0','1']:
			out='het'
		elif sorted(geno.split('|'))==['1','1']:
			out='homo'
	return out

def SV_type_decide(pin):
	flag=0
	for x in pin[7].split(';'):
		if 'SVTYPE' in x:
			flag=1
			return x.split('=')[1]
	if flag==0:
		return 'Error'

fin=open(file_in_1)
SV_Hash={}
while True:
	pin=fin.readline().strip().split()
	if not pin: break
	if not pin[0][0]=='#': 
		SV_Type=SV_type_decide(pin)
		if not SV_Type=='Error': 
			if SV_Type=='CNV':
				allele_types=[1]
				for x in pin[4].split(','):
					allele_types.append(int(x.replace('<CN','').replace('>','')))
				gt_num=genotype_report(pin)
				gt_CN=[allele_types[x] for x in gt_num]
				if sorted(gt_CN)==[0,1]:
					SV_Type='DEL'
					gt_info='het'
				elif sorted(gt_CN)==[0,0]:
					SV_Type='DEL'
					gt_info='homo'
				else:
					if sorted(gt_CN)[0]==1 and sorted(gt_CN)[1]>1:
						SV_Type='DUP'
						gt_info='het'
					elif  sorted(gt_CN)[0]>1:
						SV_Type='DUP'
						gt_info='homo'
				pos_info=end_cordi_calcu(pin)
				if not pos_info=='Error':
						if not SV_Type in SV_Hash.keys():
							SV_Hash[SV_Type]={}
						if not pos_info[0] in SV_Hash[SV_Type].keys():
							SV_Hash[SV_Type][pos_info[0]]={}
						if not int(pos_info[1]) in SV_Hash[SV_Type][pos_info[0]].keys():
							SV_Hash[SV_Type][pos_info[0]][int(pos_info[1])]=[]
						SV_Hash[SV_Type][pos_info[0]][int(pos_info[1])].append(pos_info+[gt_info])
			else:
				pos_info=end_cordi_calcu(pin)
				gt_info=genotype_calcu(pin)
				if not pos_info=='Error':
					if not gt_info=='Error':
						if not SV_Type in SV_Hash.keys():
							SV_Hash[SV_Type]={}
						if not pos_info[0] in SV_Hash[SV_Type].keys():
							SV_Hash[SV_Type][pos_info[0]]={}
						if not int(pos_info[1]) in SV_Hash[SV_Type][pos_info[0]].keys():
							SV_Hash[SV_Type][pos_info[0]][int(pos_info[1])]=[]
						SV_Hash[SV_Type][pos_info[0]][int(pos_info[1])].append(pos_info+[gt_info])
				else:
					pos_info=[pin[0],int(pin[1]),int(pin[1])]
					gt_info=genotype_calcu(pin)
					if not gt_info=='Error':
						if not SV_Type in SV_Hash.keys():
							SV_Hash[SV_Type]={}
						if not pos_info[0] in SV_Hash[SV_Type].keys():
							SV_Hash[SV_Type][pos_info[0]]={}
						if not int(pos_info[1]) in SV_Hash[SV_Type][pos_info[0]].keys():
							SV_Hash[SV_Type][pos_info[0]][int(pos_info[1])]=[]
						SV_Hash[SV_Type][pos_info[0]][int(pos_info[1])].append(pos_info+[gt_info])

fin.close()
for k1 in SV_Hash.keys():
	fo=open(file_in_1.replace('.vcf','.'+k1+'.bed'),'a')
	for k2 in sorted(SV_Hash[k1].keys()):
		for k3 in sorted(SV_Hash[k1][k2].keys()):
			for k4 in SV_Hash[k1][k2][k3]:
				print >>fo, ' '.join([str(i) for i in k4])
	fo.close()





