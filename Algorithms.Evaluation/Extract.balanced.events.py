#abstract balanced events from coverge file:
#command='python SVelter.py --ppre /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/bp_files/NA12878_S1/SPCff10.CluCff15.AlignCff0.0/ -o /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878.balanced.red --ref /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/reference/genome.fa'
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
pathin=dict_opts['--ppre']
if not pathin[-1]=='/':
	pathin+='/'

chromos=[]
ref=dict_opts['--ref']
fin=open(ref+'.fai')
for line in fin:
	pin=line.strip().split()
	chromos.append(pin[0])

fin.close()
fout=dict_opts['-o']
qual_cff=-50
def dist_flag(pin1):
	temp1=[]
	for x in pin1:
		if x in chromos:
			temp1.append([x])
		else:
			temp1[-1].append(int(x))
	flag=0
	for x in temp1:
		for x2 in range(len(x)-2):
			if x[x2+2]-x[x2+1]<100:
				flag+=1
	return flag

def bp_to_let(del_info_unit):
	flag=0
	for i in del_info_unit[0]:
	    if i in chromos or not i.isdigit():
	        flag+=1
	if not flag==0:
		letter=''.join([chr(i+97) for i in range(len(del_info_unit[0])-2*flag)])
		letters='/'.join([letter,letter])
		return letters
	else:
		return 0

def balance_SV_flag(let1,let2):
	flag=0
	if not let2==let1:
		for x in let1.split('/')[0]:
			if x in let2:
				if not x=='^':
					if not let2.count(x)==2:
						flag+=1
			else:
				flag+=1
	else:
		flag+=1
	return flag

fo=open(fout,'w')
for k1 in os.listdir(pathin):
	if k1.split('.')[-1]=='coverge':
		fin=open(pathin+k1)
		while True:
			pin1=fin.readline().strip().split()
			if not pin1: break
			pin2=fin.readline().strip().split()
			pin3=fin.readline().strip().split()
			pin4=fin.readline().strip().split()
			pin5=fin.readline().strip().split()
			if float(pin4[-1])-float(pin3[-1])>qual_cff:
				if dist_flag(pin1)==0:
					let1=bp_to_let([pin1])
					let2=pin2[0]
					if balance_SV_flag(let1,let2)==0:
						print >>fo, ' '.join(pin1+[let1,let2])
		fin.close()

fo.close()






