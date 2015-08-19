#!/usr/bin/python

#!python
#command='python SVelter.vcf2BPvcf.py --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.fasta -r /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.homo/sv.rec/homo.homo.TRA.rec.BPLink.vcf -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.homo/SVelter/SVelter_homo_RD_10_RL_101_sorted.BPLink.vcf'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import time
import datetime
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

opts,args=getopt.getopt(sys.argv[1:],'r:i:o:h:',['ref=','prefix=','appdix='])
dict_opts=dict(opts)
ref=dict_opts['--ref']
chromos=[]
fref=open(ref+'.fai')
for line in fref:
	pin=line.strip().split()
	chromos.append(pin[0])

fref.close()
ref_file=dict_opts['-r']
in_file=dict_opts['-i']
ref_hash=hash_readin(ref_file)
in_hash=SVelter_hash_readin(in_file)

compare_hash={}
for k1 in ref_hash.keys():
	compare_hash[k1]={}
	for k2 in ref_hash[k1].keys():
		if not k2 in compare_hash[k1].keys():
			compare_hash[k1][k2]=[]
		ref_pos=[int(k2.split('_')[1]),int(k2.split('_')[-3])]
		alt_keys=[]
		if k1 in in_hash.keys():
			for k3 in in_hash[k1].keys():
				alt_pos=[int(k3.split('_')[1]),int(k3.split('_')[-3])]
				if alt_pos[0]>ref_pos[1]: continue
				elif alt_pos[1]<ref_pos[0]: continue
				else:
					alt_keys.append(k3)
		alt_k3=alt_keys_filter(k2,alt_keys)
		temp_compare=[]
		ref_recs=[]
		alt_recs=[]
		for a1 in ref_hash[k1][k2].keys():
			ref_recs+=ref_hash[k1][k2][a1]
		for a2 in alt_k3:
			for a3 in in_hash[k1][a2].keys():
				alt_recs+=in_hash[k1][a2][a3]
		temp_compare=[0,len(ref_recs),len(alt_recs)]
		alt_used_recs=alt_recs
		for ax in ref_recs:
			for ay in alt_recs:
				if abs(int(ax[0])-int(ay[0]))<100 and descript_compare(ax[1],ay[1])==1:
						alt_recs.remove(ay)
						temp_compare[0]+=1
						break
		compare_hash[k1][k2]=temp_compare+alt_k3
		#print temp_compare+alt_k3

fout=dict_opts['-i'].replace('vcf','stat')
fo=open(fout,'w')
for x in chromos:
	if x in compare_hash.keys():
		for y in compare_hash[x]:
			print >>fo, ' '.join([str(i) for i in compare_hash[x][y]+[y]])
			#print ' '.join([str(i) for i in compare_hash[x][y]+[y]])

fo.close()












