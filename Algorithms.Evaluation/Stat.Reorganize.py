#command='python Stat.Reorganize.py -i /mnt/EXT/Mills-data/xuefzhao/projects.axiom/1000Genome.Validation/vcf_files/Pseudo.ROC.Mappable.min100.max1000000000.Stats'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import re
import pickle
import time
import datetime
import random
import numpy
import glob
import numpy as np
from scipy.stats import scoreatpercentile
opts,args=getopt.getopt(sys.argv[1:],'i:o:',['ref='])
dict_opts=dict(opts)
filein=dict_opts['-i']
fileout=filein+'.WG'
fin=open(filein)
fin.readline().strip().split()
in_hash={}
for line in fin:
	pin=line.strip().split()
	if not pin[0] in in_hash.keys():
		in_hash[pin[0]]={}
	if not pin[1] in in_hash[pin[0]].keys():
		in_hash[pin[0]][pin[1]]=[[],[],[]]
	in_hash[pin[0]][pin[1]][0].append(int(pin[3]))
	in_hash[pin[0]][pin[1]][1].append(int(pin[4]))
	in_hash[pin[0]][pin[1]][2].append(int(pin[5]))

fin.close()

for k1 in in_hash.keys():
	for k2 in in_hash[k1].keys():
		temp=[numpy.sum(i) for i in in_hash[k1][k2]]
		in_hash[k1][k2]=temp

fo=open(fileout,'w')
for k1 in in_hash.keys():
	for k2 in in_hash[k1].keys():
		print >>fo, ' '.join([k1,k2]+[str(i) for i in in_hash[k1][k2]])

fo.close()









