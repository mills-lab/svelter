#!/usr/bin/env python

#!python
#command='sam2fasta -i temp.sam'
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
from scipy.stats import scoreatpercentile
opts,args=getopt.getopt(sys.argv[1:],'i:',)
dict_opts=dict(opts)
file_in_1=dict_opts['-i']
fin=open(file_in_1)
fo=open(file_in_1.replace('.sam','.fa'),'w')
for line in fin:
	pin=line.strip().split()
	if not pin[0][0]=='@':
		print >>fo, '>'+pin[0]
		print >>fo, pin[9]

fin.close()
fo.close()



