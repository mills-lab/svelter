#!/usr/bin/env python

#!python
#command='python bedpe2bed.py -i /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.het/Lumpy/Simulate.het.RD.10.RL.101.sorted.bedpe --ref /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta'

#sys.argv=command.split()[1:]
import os
import sys
import numpy
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref='])
dict_opts=dict(opts)
out_hash={}
fin=open(dict_opts['-i'])
for line in fin:
	pin=line.strip().split()
	chrom=pin[0]
	x1=int(numpy.mean([int(pin[1]),int(pin[2])]))
	x2=int(numpy.mean([int(pin[4]),int(pin[5])]))
	sv_key=pin[10].split(':')[1][:3]
	if not sv_key in out_hash.keys():
		out_hash[sv_key]={}
	if not chrom in out_hash[sv_key].keys():
		out_hash[sv_key][chrom]={}
	if not x1 in out_hash[sv_key][chrom].keys():
		out_hash[sv_key][chrom][x1]=[]
	if not x2 in out_hash[sv_key][chrom][x1]:
		out_hash[sv_key][chrom][x1].append(x2)

fin.close()
chromos=[]
ref=dict_opts['--ref']
fin=open(ref+'.fai')
for line in fin:
	pin=line.strip().split()
	chromos.append(pin[0])

fin.close()
for k1 in out_hash.keys():
	fo=open(dict_opts['-i'].replace('bedpe',k1+'.bed'),'w')
	for k2 in chromos:
		if k2 in out_hash[k1].keys():
			for k3 in sorted(out_hash[k1][k2].keys()):
				for k4 in sorted(out_hash[k1][k2][k3]):
					print >>fo, '\t'.join([str(i) for i in [k2,k3,k4]])

fo.close()


