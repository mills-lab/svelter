#!/usr/bin/python

#!python
#python pindel.vcf2bed.py -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/Pindel/Pindel.het.RD.10.LargerThan10.vcf
#sys.argv=command.split()[1:]
def genotype_calcu(pin):
	gt_index=pin[8].split(':').index('GT')
	geno=pin[9].split(':')[gt_index]
	out=0
	if '/' in geno:
		if geno.split('/')==['0','0']:
			out='Error!'
		elif sorted(geno.split('/'))==['0','1']:
			out='het'
		elif sorted(geno.split('/'))==['1','1']:
			out='homo'
	elif '|' in geno:
		if geno.split('|')==['0','0']:
			out='Error!'
		elif sorted(geno.split('|'))==['0','1']:
			out='het'
		elif sorted(geno.split('|'))==['1','1']:
			out='homo'
	return out

import os
import sys
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
DELOut=dict_opts['-i'].replace('.vcf','.DEL.bed')
DUPOut=dict_opts['-i'].replace('.vcf','.DUP.bed')
INVOut=dict_opts['-i'].replace('.vcf','.INV.bed')
fdel=open(DELOut,'w')
fdup=open(DUPOut,'w')
finv=open(INVOut,'w')
for line in fin:
	pin=line.strip().split()
	if not pin[0][0]=='#':
		chrom=pin[0]
		start=pin[1]
		info=pin[7]
		flag=0
		gt=genotype_calcu(pin)
		if not gt=='Error!':
			for x in info.split(';'):
				if 'END' in x:
					end=x.split('=')[1]
					flag+=1
				if 'SVTYPE' in x:
					SV=x.split('=')[1]
					flag+=1
			if flag==2:
				if SV=='DEL':
					print >>fdel, ' '.join([chrom,start,end,gt])
				elif SV=='DUP' or SV=='DUP:TANDEM':
					print >>fdup, ' '.join([chrom,start,end,gt])
				elif SV=='INV':
					print >>finv, ' '.join([chrom,start,end,gt])

fin.close()
fdel.close()
fdup.close()
finv.close()



