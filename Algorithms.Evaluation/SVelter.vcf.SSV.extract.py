#command='python SVelter.py --ref /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -i /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.comp/comp.sim -o /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.comp/comp'
#command='python SVelter.vcf.SSV.extract.py -i /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.het/SVelter/Simulate.het.RD.20.RL.101.vcf'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import re
import pickle
import time
import datetime
import random
opts,args=getopt.getopt(sys.argv[1:],'i:o:',['ref='])
dict_opts=dict(opts)
filein=dict_opts['-i']
fileout=filein.replace('.vcf','.simple.vcf')
fin=open(filein)
fo=open(fileout,'w')
for line in fin:
	pin=line.strip().split()
	if pin[0][0]=='#':
		print >>fo, '\t'.join(pin)
	else:
		if pin[2][-1]=='S':
			print >>fo, '\t'.join(pin)

fin.close()
fo.close()



