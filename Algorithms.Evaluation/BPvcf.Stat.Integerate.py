#!/usr/bin/python

#!python
#command='python SVelter.vcf2BPvcf.py -p /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/SVelter/'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import time
import datetime
import numpy
def allele_decide(pin):
	if '|' in pin[-1]:
		if pin[-1]=='1|0':
			return ['a']
		elif pin[-1]=='0|1':
			return ['b']
		elif pin[-1]=='0|0':
			return []
		elif pin[-1]=='1|1':
			return ['a','b']
	elif '/' in pin[-1]:
		if pin[-1]=='1/0':
			return ['a']
		elif pin[-1]=='0/1':
			return ['b']
		elif pin[-1]=='0/0':
			return []
		elif pin[-1]=='1/1':
			return ['a','b']

def hash_readin(hash_file):
	out={}
	fin=open(hash_file)
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='#':
			chromo=pin[0]
			start=int(pin[1])
			end=int(pin[2].split('_')[-3])
			key=pin[2]
			allele=allele_decide(pin)
			if not chromo in out.keys():
				out[chromo]={}
			if not key in out[chromo].keys():
				out[chromo][key]={}
			for x in allele:
				if not x in out[chromo][key].keys():
					out[chromo][key][x]=[]
				out[chromo][key][x].append([pin[1],pin[4]])
	fin.close()
	return out

def SVelter_hash_readin(hash_file):
	out={}
	fin=open(hash_file)
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='#':
			chromo=pin[0]
			start=int(pin[1])
			#end=int(pin[2].split('_')[-3])
			key='_'.join(pin[2].split('_')[:-3]+pin[2].split('_')[-2:])
			allele=allele_decide(pin)
			if not chromo in out.keys():
				out[chromo]={}
			if not key in out[chromo].keys():
				out[chromo][key]={}
			for x in allele:
				if not x in out[chromo][key].keys():
					out[chromo][key][x]=[]
				out[chromo][key][x].append([pin[1],pin[4]])
	fin.close()
	return out

def keys_dissect(k2):
	k2_pos={}
	for x in k2.split('_')[:-2]:
		if x in chromos:
			chr_cur=x
			k2_pos[chr_cur]=[]
		else:
			k2_pos[chr_cur].append(int(x))
	return k2_pos

def alt_keys_filter(k2,alt_keys):
	#eg k2: k2='1_185365626_185368946_185477605_ab/ab_ab/ba'
	#eg alt_keys: ['1_185365458_185477721_a/a_a/aa', '1_185369009_185477437_a/a_/a', '1_185365709_185368769_a/a_/a']
	out=[]
	k2_pos=keys_dissect(k2)
	for k3 in alt_keys:
		k3_pos=keys_dissect(k3)
		flag=0
		if k2_pos.keys()==k3_pos.keys():
			for x1 in k2_pos.keys():
				for x2 in k2_pos[x1]:
					for x3 in k3_pos[x1]:
						if abs(x3-x2)<100:
							flag+=1
		if not flag==0:
			out.append(k3)
	return out

def descript_compare(des1,des2):
	#eg: des1=']20:21417806]A', des2=']20:21417806]A'
	if des1==des2:
		return 1
	else:
		flag=0
		if ']' in des1 and not '[' in des1:
			if  ']' in des2 and not '[' in des2:
				if abs(int(des1.split(':')[1].split(']')[0].split('[')[0])-int(des2.split(':')[1].split(']')[0].split('[')[0]))<100:
					flag=1
					return 1 
		elif '[' in des1 and not ']' in des1:
			if  '[' in des2 and not ']' in des2:
				if abs(int(des1.split(':')[1].split(']')[0].split('[')[0])-int(des2.split(':')[1].split(']')[0].split('[')[0]))<100:
					flag=1
					return 1 
		elif '[' in des1 and ']' in des1:
			if  '[' in des2 and ']' in des2:
				if abs(int(des1.split(':')[1].split(']')[0].split('[')[0])-int(des2.split(':')[1].split(']')[0].split('[')[0]))<100:
					flag=1
					return 1 
		if flag==0: 
			return 0

opts,args=getopt.getopt(sys.argv[1:],'p:r:i:o:h:',['ref=','prefix=','appdix='])
dict_opts=dict(opts)
if not dict_opts['-p'][-1]=='/':
	dict_opts['-p']+='/'

out_stat={}
for k1 in os.listdir(dict_opts['-p']):
	if os.path.isfile(dict_opts['-p']+k1):
		if 'BPLink' in k1 and k1.split('.')[-1]=='stat' and not k1=='BPLink.Integrated.stat':
			out_stat[k1.replace('.stat','')]=[]
			fin=open(dict_opts['-p']+k1)
			TP=[]
			AllP=[]
			for line in fin:
				pin=line.strip().split()
				if pin[-1].split('_')[0] in [str(i) for i in range(23)]:
					TP.append(int(pin[0]))
					AllP.append(int(pin[1]))
			fin.close()
			TPNum=numpy.sum(TP)
			AllPNum=numpy.sum(AllP)
			fin=os.popen(r'''wc -l %s'''%(dict_opts['-p']+k1.replace('stat','vcf')))
			TotalP=int(fin.readline().strip().split()[0])
			fin.close()
			out_stat[k1.replace('.stat','')].append(float(TPNum)/float(AllPNum))
			out_stat[k1.replace('.stat','')].append(float(TPNum)/float(TotalP))

fo=open(dict_opts['-p']+'BPLink.Integrated.stat','w')
for k1 in sorted(out_stat.keys()):
	print>>fo, ' '.join([str(i) for i in [k1]+out_stat[k1]])
	print ' '.join([str(i) for i in [k1]+out_stat[k1]])

fo.close()


