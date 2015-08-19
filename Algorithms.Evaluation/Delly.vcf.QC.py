#!/usr/bin/env python

#!python
#command='Delly.vcf.QC.py -i Delly_het_RD10_DEL.vcf'
#sys.argv=command.split()
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
import scipy
from scipy import stats
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.vcf','.PassQC.vcf'),'w')
for line in fin:
	pin=line.strip().split()
	if pin[0][0]=='#':
		print >>fo, '\t'.join(pin)
	else: 
		if pin[6]=='PASS':
			print >>fo, '\t'.join(pin)

fin.close()
fo.close()


