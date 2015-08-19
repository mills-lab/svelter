#!/usr/bin/python

#!python
#command='python BPvcf.Mappable.filter.py -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/sv.rec/het.het.TRA.rec.BPLink.vcf -r /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.Mappable'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import time
import datetime
opts,args=getopt.getopt(sys.argv[1:],'r:i:o:h:',['ref=','prefix=','appdix='])
dict_opts=dict(opts)
ref_path='/'.join(dict_opts['-r'].split('/')[:-1]+['/'])
ref_hash={}
for k1 in os.listdir(ref_path):
	if dict_opts['-r'].split('/')[-1] in k1:
		fin=open(ref_path+k1)
		for line in fin:
			pin=line.strip().split()
			if not pin[0] in ref_hash.keys():
				ref_hash[pin[0]]=[]
			ref_hash[pin[0]].append([int(pin[1]),int(pin[2])])
		fin.close()


fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.vcf','.2.vcf'),'w')
for line in fin:
	pin=line.strip().split()
	if pin[0][0]=='#':
		print >>fo, '\t'.join(pin)
	else:
		chrom=pin[0]
		if chrom in ref_hash.keys():
			start=int(pin[2].split('_')[1])
			if '.rec.' in dict_opts['-i'].split('/')[-1]:
				end=int(pin[2].split('_')[-3])
			else:
				end=int(pin[2].split('_')[-4])
			flag=0
			for x in ref_hash[chrom]:
				if x[0]-1<start and x[1]+1>end:
					flag+=1
			if not flag==0:
				print >>fo, '\t'.join(pin)

fin.close()
fo.close()


