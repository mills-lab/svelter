#!/usr/bin/env python

#!python
#command='python /nfs/remills-data/xuefzhao/SVelter/SVelter.py -p /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/Delly'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import numpy
import re
opts,args=getopt.getopt(sys.argv[1:],'p:')
dict_opts=dict(opts)

if not dict_opts['-p'][-1]=='/':
	dict_opts['-p']+='/'

for x in os.listdir(dict_opts['-p']):
	if x.split('.')[-1]=='vcf':
		info1='_'.join(x.replace('.RL.101.','.').split('.')[:-1])+'.vcf'
		os.system(r'''cp %s %s'''%(dict_opts['-p']+x,dict_opts['-p']+info1))



