#!/usr/bin/python

#!python
#python /nfs/remills-scratch/datasets/Simulation.Xuefang/compare.scripts/bed2Mappbalebed.py -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/Delly/Delly.RD.5.RL.101.DEL.DEL.bed --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.Mappable

#sys.argv=command.split()[1:]
import os
import sys
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref='])
dict_opts=dict(opts)
ref_path='/'.join(dict_opts['--ref'].split('/')[:-1])
ref_prefix=dict_opts['--ref'].split('/')[-1]
ref_files=[]
for k1 in os.listdir(ref_path):
	if ref_prefix in k1 and k1.split('.')[-1]=='bed':
		ref_files.append(ref_path+'/'+k1)

ref_hash={}
for k1 in ref_files:
	chrom=k1.split('/')[-1].replace('.bed','').replace(ref_prefix+'.','')
	ref_hash[chrom]=[]
	fin=open(k1)
	for line in fin:
		pin=line.strip().split()
		ref_hash[chrom].append([int(pin[1]),int(pin[2])])
	fin.close()

file_in=dict_opts['-i']
file_out=file_in.replace('.bed','.Mappable.bed')
fin=open(file_in)
in_hash={}
chroms=[]
for line in fin:
	pin=line.strip().split()
	if not pin[0] in in_hash.keys():
		in_hash[pin[0]]=[]
		chroms.append(pin[0])
	in_hash[pin[0]].append([int(i) for i in pin[1:3]]+pin[3:])

fin.close()

fo=open(file_out,'w')
for k1 in chroms:
	if k1 in in_hash.keys() and k1 in ref_hash.keys():
		for k2 in in_hash[k1]:
			flag=0
			for k3 in ref_hash[k1]:
				if k3[0]-1<k2[0] and k3[1]+1>k2[1]:
					flag+=1
			if not flag==0:
				print >>fo, ' '.join([str(i) for i in [k1]+k2])

fo.close()


