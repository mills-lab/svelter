#!/usr/bin/env python

#!python
#command='SVrec2CSVrec.py -i SVelter_SAMN02744161_flag1.SV.rec'
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

def simple_del_decide(k1,k2):
	flag=0
	if not '^' in k2:
		for x in range(len(k2.split('/')[0])-1):
			if ord(k2.split('/')[0][x+1])-ord(k2.split('/')[0][x])>0:
				continue
			else:
				flag+=1
		for x in range(len(k2.split('/')[1])-1):
			if ord(k2.split('/')[1][x+1])-ord(k2.split('/')[1][x])>0:
				continue
			else:
				flag+=1
	else:
		flag+=1
	return flag

def simple_del_region(k1,k2):
	if simple_del_decide(k1,k2)==0:
		out=[[],[]]
		for x in k1.split('/')[0]:
			if not x in k2.split('/')[0]:
				out[0].append(x)
		for x in k1.split('/')[1]:
			if not x in k2.split('/')[1]:
				out[1].append(x)
		if out[0]==out[1]:
			return 0
		elif out[0]==[] and not out[1]==[]:
			return 0
		elif out[1]==[] and not out[0]==[]:
			return 0
		else:
			return 1
	else:
		return 1

def simple_inv_decide(k1,k2):
	if k2.replace('^','')==k1:
		return 0
	else:
		return 1

def letter_separate(k2):
	out=[[],[]]
	for x1 in k2.split('/')[0]:
		if not x1=='^':
			out[0].append(x1)
		else:
			out[0][-1]+=x1
	for x1 in k2.split('/')[1]:
		if not x1=='^':
			out[1].append(x1)
		else:
			out[1][-1]+=x1
	return out

def letter_recombine(k2):
	k2_new=letter_separate(k2)
	out=[[],[]]
	if k2_new[0]==[]:
		out=[[],[k2_new[1][0]]]
	if k2_new[1]==[]:
		out=[[k2_new[0][0]],[]]
	if not k2_new[0]==[] and not k2_new[1]==[]:
		out=[[k2_new[0][0]],[k2_new[1][0]]]
	for x1 in k2_new[0]:
		if not '^' in x1 and not '^' in out[0][-1]:
			if ord(x1)-ord(out[0][-1][-1])==1:
				out[0][-1]+=x1
			else:
				out[0].append(x1)
		elif '^' in x1 and '^' in out[0][-1]:
			if ord(x1[0])-ord(out[0][-1][0])==-1:
				out[0][-1]=x1[0]+out[0][-1]
			else:
				out[0].append(x1)
		else:
			out[0].append(x1)
	for x1 in k2_new[1]:
		if not '^' in x1 and not '^' in out[1][-1]:
			if ord(x1)-ord(out[1][-1][-1])==1:
				out[1][-1]+=x1
			else:
				out[1].append(x1)
		elif '^' in x1 and '^' in out[1][-1]:
			if ord(x1[1])-ord(out[1][-1][1])==-1:
				out[1][-1]=x1[1]+out[1][-1]
			else:
				out[1].append(x1)
		else:
			out[1].append(x1)
	if not out[0]==[]:
		del out[0][0]
	if not out[1]==[]:
		del out[1][0]
	return out

def simple_tandemdup_decide(k1,k2):
	out=letter_recombine(k2)
	out2=[[],[]]
	for x in out[0]:
		if out2[0]==[]:
			out2[0].append(x)
		else:
			if not x==out2[0][-1][-len(x):]:
				out2[0].append(x)
	for x in out[1]:
		if out2[1]==[]:
			out2[1].append(x)
		else:
			if not x==out2[1][-1][-len(x):]:
				out2[1].append(x)
	if '/'.join([''.join(out2[0]),''.join(out2[1])])==k1:
		return 0
	else:
		return 1

def simple_dup_decide2(k1,k2):
	out=[[],[]]
	for x in k2.split('/')[0]:
		if out[0]==[]:
			out[0].append(x)
		else:
			if not x==out[0][-1]:
				out[0].append(x)
	for x in k2.split('/')[1]:
		if out[1]==[]:
			out[1].append(x)
		else:
			if not x==out[1][-1]:
				out[1].append(x)
	if '/'.join([''.join(out[0]),''.join(out[1])])==k1:
		return 0
	else:
		return 1

opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.SV.rec','.CSV.rec'),'w')
for line in fin:
        pin=line.strip().split()
        if not pin[0]=='a/a':
                if float(pin[-1])>0:
                        if simple_del_region(pin[0],pin[1])==0:
                                continue
                                #print pin[:2]
                        elif simple_inv_decide(pin[0],pin[1])==0:
                                continue
                                #print pin[:2]
                        elif simple_tandemdup_decide(pin[0],pin[1])==0:
                                print pin[:2]
                        elif simple_dup_decide2(pin[0],pin[1])==0:
                                print pin[:2]
                        else:
                                print >>fo, '\t'.join(pin[:-1])

fin.close()
fo.close()

