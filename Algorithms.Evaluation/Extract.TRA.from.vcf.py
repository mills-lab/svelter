#abstract all translocations from coverge file:
#command='python SVelter.py -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/SVelter/SVelter_het_RD40.vcf'
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
opts,args=getopt.getopt(sys.argv[1:],'i:o:',['ppre=','ref='])
dict_opts=dict(opts)
filein=dict_opts['-i']
fileout=dict_opts['-i'].replace('.vcf','.tra.vcf')
fin=open(filein)
fo=open(fileout,'w')
for line in fin:
	pin=line.strip().split()
	if pin[0][0]=='#':
		continue
	else:
		if pin[6]=='PASS':
			if ']' in pin[4] or '[' in pin[4]:
				print >>fo, '\t'.join(pin)

fin.close()
fo.close()








