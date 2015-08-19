#!/usr/bin/env python

#!Python
#Usage:
#SVelter4.StructureResolvation.py --ppre workdir --ref ref.fa -s sample.bam -f file.txt --NullModel S
#this code is used to build the null model of insert length distribution (based on bimodal distribution), read depth distribution (based on NB dist, 100bp/bin, corrected by GC content), number of clipped reads, and number of read pairs with abnormai directions(FF). reads to build the null morday would come from cn2 regions larger than a customized size(1kb as default.)
#output data structures: folder named "InputDataSet", containing all input bam file;	folder "NullDist", containing four input folders: 'InsertLengthDist','ReadDepthDist','SplitReadsDist','AbnormalDirectionDist'
#Usage:
#For debug use only
#command='python code --NullModel S --ppre /nfs/remills-data/xuefzhao/projects.axiom/ -f /nfs/remills-data/xuefzhao/projects.axiom/bp_files/NA12878_S1/SPCff5.CluCff10.AlignCff0.0/chr2/NA12878_S1.chr2.7702.txt -s /nfs/remills-scratch/datasets/platinum/alignment/NA12878_S1.bam -n 100000 -o /nfs/remills-data/xuefzhao/projects.axiom/bp_files/NA12878_S1/SPCff5.CluCff10.AlignCff0.0/chr2/ --ref /nfs/remills-scratch/reference/hg19_platinum/genome.fa --MoveFlag 2'
#only invert, insert and delete are done, allow insert undeleted block
#sys.argv=command.split()[1:]
#eg on chr3:
#chr3	162507460	162512134	162525742	162545363	162547646	162547656	162626332
import os
import sys
import getopt
import random
import scipy
import math
import numpy
import pickle
from math import sqrt,pi,exp
import scipy
from scipy.stats import norm
import time
import datetime
time1=time.time()
import itertools

#CluNum=10
#Ln_CluLen=500

def pdf_calculate(x,alpha,mean1,mean2,std1,std2,upper_limit,lower_limit,Penalty_For_InsertLengthZero):
	Alpha=numpy.min([alpha,1-alpha])
	if mean1<mean2:
		Mean1=mean1
		Std1=std1
		Mean2=mean2
		Std2=std2
	elif mean1>mean2: 
		Mean1=mean2
		Std1=std2
		Mean2=mean1
		Std2=std1
	if not Alpha==0:
		if x<upper_limit and x>lower_limit:
			return math.log(Alpha/math.sqrt(2*math.pi*math.pow(Std1,2))*math.exp(-math.pow((x-Mean1),2)/(2*math.pow(Std1,2)))+(1-Alpha)/math.sqrt(2*math.pi*math.pow(Std2,2))*math.exp(-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))))
		elif x>=upper_limit:
			test1=math.pow((x-Mean1),2)/(2*math.pow(Std1,2))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
			if test1<200:
				return math.log(Alpha)-0.5*math.log(2*math.pi*math.pow(Std1,2))-math.pow((x-Mean1),2)/(2*math.pow(Std1,2))+math.log(1+(1-Alpha)*Std1/(Alpha*Std2)*math.exp(test1))
			elif test1>200:
				return math.log(1-Alpha)-0.5*math.log(2*math.pi*math.pow(Std2,2))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
		elif x<=lower_limit and not x==0:
			test2=-math.pow((x-Mean1),2)/(2*math.pow(Std1,2))+math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
			if test2<200:
				return math.log(1-Alpha)-0.5*math.log(2*math.pi*math.pow(Std2,2))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))+math.log(1+Alpha*Std2/((1-Alpha)*Std1)*math.exp(test2))
			if test2>200:
				return  math.log(Alpha)-0.5*math.log(2*math.pi*math.pow(Std1,2))-math.pow((x-Mean1),2)/(2*math.pow(Std1,2))
		elif x==0:
			return Penalty_For_InsertLengthZero
	else:
		if not x==0:
			return math.log(1-Alpha)-math.log(math.sqrt(2*math.pi*math.pow(Std2,2)))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
		elif x==0:
			return Penalty_For_InsertLengthZero

def bimodal_cdf_solver(cdf,alpha,mean1,mean2,std1,std2):
	if cdf>0.5:
		x=numpy.min([mean1,mean2])
  		while True:
			x+=0.1
			fc=alpha*(norm.cdf((x-mean1)/std1))+(1-alpha)*(norm.cdf((x-mean2)/std2))
			if abs(fc-cdf)<0.001:
				return x
 	elif cdf<0.5:
		x=numpy.max([mean1,mean2])
		while True:
			x-=0.1
			fc=alpha*(norm.cdf((x-mean1)/std1))+(1-alpha)*(norm.cdf((x-mean2)/std2))
			if abs(fc-cdf)<0.001:
				return x

def norm_cdf_solver(cdf,mean,std):
	x=mean
	if cdf>0.5:
  		while True:
			x+=0.1
			fc=(norm.cdf((x-mean)/std))
			if abs(fc-cdf)<0.001:
				return x
 	elif cdf<0.5:
		while True:
			x-=0.1
			fc=(norm.cdf((x-mean)/std))
			if abs(fc-cdf)<0.001:
				return x

def cdf_solver_application(Insert_Len_Stat,cdf):
	fstat=open(Insert_Len_Stat)
	temp=fstat.readline()
	temp=fstat.readline()
	temp=fstat.readline()
	if model_comp=='S':
		data1=fstat.readline().strip().split()
		fstat.close()
		flank_out=int(round(norm_cdf_solver(float(cdf),float(data1[1]),float(data1[2]))))
	elif model_comp=='C':
		data1=fstat.readline().strip().split()
		fstat.readline()
		data2=fstat.readline().strip().split()
		flank_out=int(round(bimodal_cdf_solver(float(cdf),float(data1[0]),float(data1[1]),float(data2[1]),float(data1[2]),float(data2[2]))))
		fstat.close()
	return flank_out

def GC_Stat_ReadIn(BamN,GC_Stat_Path,affix):
	GC_Stat_File=GC_Stat_Path+'/'+BamN+'.'+genome_name+affix
	f_GC_stat=open(GC_Stat_File)
	CN2_Region={}
	Chromos=f_GC_stat.readline().strip().split()
	GC_Content=f_GC_stat.readline().strip().split()
	for key_1 in Chromos:
			CN2_Region[key_1]={}
			for key_2 in GC_Content:
					CN2_Region[key_1][int(key_2)]=f_GC_stat.readline().strip().split()
	f_GC_stat.close()
	return [CN2_Region,Chromos,GC_Content]

def Reads_Direction_Detect(flag):
#flag is the number on the second position of each read in bam file
	#	if int(pi[1])&2>0: #two both mapped
	#	elif int(pi[1])&4>0: #the read itself mapped, mate not
	#	elif int(pi[1])&8>0: #mate mapped, the read itself not
	flag2=int(flag)
	if int(flag2)&16==0: 
			direct_1='+'
	elif int(flag2)&16>0:
			direct_1='-'
	if int(flag2)&32==0:
			direct_2='+'
	elif int(flag2)&32>0: 
			direct_2='-'
	return([direct_1,direct_2])

def cigar2reaadlength(cigar):
	import re
	pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
	cigars=[]
	for m in pcigar.finditer(cigar):
		cigars.append((m.groups()[0],m.groups()[1]))
	MapLen=0
	for n in cigars:
		if n[1]=='M' or n[1]=='D' or n[1]=='N':
			MapLen+=int(n[0])
	return MapLen

def Reads_block_assignment_2(bps,letters,block1,block2,flank):
	#In this function, we assign a read to the block where majority of the read falls in
	#'majority' is defined by the relative length. if the larger part of a read fells in block a, then the whole read are in block a 
	#to save time, we use this for now
	#may complex it later
	bps2=[int(i) for i in bps]
	relative_bps=[i-numpy.min(bps2) for i in bps2]
	if Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)==Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank):
		return Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)
	elif not Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)==Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank):
		length_left=Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)[-1]-numpy.min([block1,block2])
		length_right=-Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank)[-2]+numpy.max([block1,block2])
		if not length_left<length_right:
			return Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)
		elif length_left<length_right:
			return Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank)

def letters_bps_produce(letters,bps,flank):
	letters_bps={}
	letters_relative_bps={}
	letters_bps['left']=[bps[0]-flank,bps[0]]
	letters_relative_bps['left']=[-flank,0]
	for i in range(len(bps)-1):
		letters_relative_bps[letters[i]]=[bps[i]-bps[0],bps[i+1]-bps[0]]
		letters_bps[letters[i]]=[bps[i],bps[i+1]]
	letters_bps['right']=[bps[-1],bps[-1]+flank]
	letters_relative_bps['right']=[bps[-1]-bps[0],bps[-1]-bps[0]+flank]
	return [letters_bps,letters_relative_bps]

def Reads_block_assignment_3(bps,letters,letters_relative_bps,block1,block2,flank):
	bps2=[int(i) for i in bps]
	relative_bps=[-flank]+[i-numpy.min(bps2) for i in bps2]+[bps2[-1]-bps2[0]+flank]
	relative_letters=['left']+letters+['right']
	meanbl=int(numpy.mean([block1,block2]))
	if meanbl<relative_bps[0]:
		return ['leftError',relative_bps[0],relative_bps[0]]
	if not meanbl<relative_bps[-1]:
		return ['rightError', relative_bps[-1], relative_bps[-1]]
	else:
		resultbl=[relative_letters[i] for i in range(len(relative_bps)-1) if meanbl<relative_bps[i+1] and not meanbl<relative_bps[i]]
		return resultbl+letters_relative_bps[resultbl[0]]

def Reads_block_assignment_1(bps,letters,position):
	if position<bps[0]-2*flank or position>bps[-1]+2*flank:
		return '0'
	else:
		if position<bps[0]+1 and position>bps[0]-2*flank-1:
			return letters[0]
		elif position>bps[-1]-1 and position<bps[-1]+2*flank+1:
			return letters[-1]
		else:
			for i in range(len(bps)-1):
				if not position<bps[i] and not position>bps[i+1]:
					return letters[i]

def RD_Index_ReadIn(ppre_Path,BamN, chromo, region):
	if not ppre_Path[-1]=='/':
		ppre_Path+='/'
	path_in=ppre_Path+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/RD_Stat/'
	file_in=BamN+'.'+chromo+'.RD.index'
	fin=open(path_in+file_in)
	pos1=int(region[0])
	pos2=int(region[1])
	while True:
		pin1=fin.readline().strip().split()
		if not pin1: break
		pin2=fin.readline().strip().split()
		reg1=int(pin1[0].split(':')[1].split('-')[0])
		reg2=int(pin1[0].split(':')[1].split('-')[1])
		if not pos1<reg1 and not pos2>reg2:
			break

def Full_Info_of_Reads_Product_3(Initial_Bam,temp_bp,temp_let,bamChr,target_region,Chr_Link):
	Letter_Double={}
	Pair_ThroughBP=[]
	Double_Read_ThroughBP=[]
	Single_Read_ThroughBP=[]
	blackList=[]
	fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,bamChr,target_region[0]-flank,target_region[-1]+flank))
	num_of_reads=0
	while True:
		pbam=fbam.readline().strip().split()
		if not pbam: break
		if int(pbam[1])&4>0: continue
		if int(pbam[1])&1024>0:continue
		if not int(pbam[4])>QCAlign or int(pbam[1])&512>0:
			blackList.append(pbam[0])
			continue
		if pbam[0] in blackList: continue
		num_of_reads+=1
		if int(pbam[1])&8>0 or not pbam[6]=='=':
			pos1=int(pbam[3])+low_qual_edge
			pos2=int(pbam[3])+cigar2reaadlength(pbam[5])-low_qual_edge
			block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
			block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
			if block1==block2:
				BlockCov[block1]+=cigar2reaadlength(pbam[5])
			else:
				reg1a=temp_bp[temp_let.index(block1)]
				reg1b=temp_bp[temp_let.index(block1)+1]
				reg2a=temp_bp[temp_let.index(block2)]
				reg2b=temp_bp[temp_let.index(block2)+1]
				rela_1=pos1-low_qual_edge-temp_bp[temp_let.index(block1)]
				rela_2=pos2+low_qual_edge-temp_bp[temp_let.index(block2)]
				Single_Read_ThroughBP.append([block1,rela_1,block2,rela_2,pbam[5]])
			if not pbam[6]=='=':			   
				if not pbam[0] in Chr_Link:
					Chr_Link[pbam[0]]=[pbam[1:9]]
				else:
					Chr_Link[pbam[0]]+=[pbam[1:9]]
		elif int(pbam[1])&8==0:
			if pbam[6]=='=':
				if not pbam[0] in Letter_Double.keys():
					Letter_Double[pbam[0]]=[pbam[:9]]
				else:
					if not pbam[:9] in Letter_Double[pbam[0]]:
						Letter_Double[pbam[0]]+=[pbam[:9]]
						if int(Letter_Double[pbam[0]][0][3])<int(Letter_Double[pbam[0]][1][3]):
							pos1=int(Letter_Double[pbam[0]][0][3])+low_qual_edge
							pos2=int(Letter_Double[pbam[0]][1][3])+cigar2reaadlength(Letter_Double[pbam[0]][1][5])-low_qual_edge
						else:
							pos1=int(Letter_Double[pbam[0]][1][3])+low_qual_edge
							pos2=int(Letter_Double[pbam[0]][0][3])+cigar2reaadlength(Letter_Double[pbam[0]][0][5])-low_qual_edge
						block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
						block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
						if block1==block2:
							BlockCov[block1]+=cigar2reaadlength(Letter_Double[pbam[0]][0][5])
							del Letter_Double[pbam[0]]
							blackList.append(pbam[0])
	fbam.close()
	for key in Letter_Double.keys():
		if key in blackList:
			del Letter_Double[key]
			continue
		if len(Letter_Double[key])==2:
			pos1=int(Letter_Double[key][0][3])
			pos2=int(Letter_Double[key][1][3])
			if not pos1>pos2:
				pos1=int(Letter_Double[key][0][3])
				pos1b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
				pos2=int(Letter_Double[key][1][3])
				pos2b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
				direct_temp=Reads_Direction_Detect(Letter_Double[key][0][1])
			elif pos1>pos2:
				pos1=int(Letter_Double[key][1][3])
				pos1b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
				pos2=int(Letter_Double[key][0][3])
				pos2b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
				direct_temp=Reads_Direction_Detect(Letter_Double[key][1][1])
			block1=Reads_block_assignment_1(temp_bp,temp_let,pos1+low_qual_edge)
			block2=Reads_block_assignment_1(temp_bp,temp_let,pos2+low_qual_edge)
			block1b=Reads_block_assignment_1(temp_bp,temp_let,pos1b-low_qual_edge)
			block2b=Reads_block_assignment_1(temp_bp,temp_let,pos2b-low_qual_edge)
			rela_1=pos1-temp_bp[temp_let.index(block1)]
			rela_2=pos2-temp_bp[temp_let.index(block2)]
			rela_1b=pos1b-temp_bp[temp_let.index(block1b)]
			rela_2b=pos2b-temp_bp[temp_let.index(block2b)]
			if block1==block1b and block2==block2b:
				Pair_ThroughBP.append([block1,rela_1,rela_1b, block2,rela_2,rela_2b]+direct_temp)
			else:
				Double_Read_ThroughBP.append([block1,rela_1,block1b,rela_1b, block2,rela_2,block2b,rela_2b]+direct_temp)
			del Letter_Double[key]
		elif len(Letter_Double[key])==1:
			if Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))==Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3])+cigar2reaadlength(Letter_Double[key][0][5])):
				BlockCov[Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))]+=cigar2reaadlength(Letter_Double[key][0][5])
				del Letter_Double[key]
	Initial_DR_Penal=0
	for j in Pair_ThroughBP:
		if not j[-2:]==['+', '-']:
			Initial_DR_Penal+=1
	for j in Double_Read_ThroughBP:
		if not j[-2:]==['+', '-']:
			Initial_DR_Penal+=1
	for j in Pair_ThroughBP:
		Initial_Cov[j[0]]+=j[2]-j[1]
		Initial_Cov[j[3]]+=j[5]-j[4]
	for j in Single_Read_ThroughBP:
		Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
		Initial_Cov[j[2]]+=j[3]
	for j in Double_Read_ThroughBP:
		if j[0]==j[2]:
			Initial_Cov[j[0]]+=j[3]-j[1]
		else:
			Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
			Initial_Cov[j[2]]+=j[3]
		if j[4]==j[6]:
			Initial_Cov[j[4]]+=j[7]-j[5]
		else:
			Initial_Cov[j[4]]+=temp_bp[temp_let.index(j[4])+1]-temp_bp[temp_let.index(j[4])]-j[5]
			Initial_Cov[j[6]]+=j[7]
	for j in Pair_ThroughBP:
		Initial_IL.append(temp_bp[temp_let.index(j[3])]-temp_bp[temp_let.index(j[0])]-j[1]+j[5])
	for j in Double_Read_ThroughBP:
		Initial_IL.append(temp_bp[temp_let.index(j[6])]-temp_bp[temp_let.index(j[0])]-j[1]+j[7])
	return [Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,num_of_reads,Initial_DR_Penal]

def Full_Info_of_Reads_Product(Initial_Bam,bps,total_bps,total_letters,bamChr,flank,QCAlign,ReadLength,chr_link):
	#	letters=[chr(97+i) for i in range(len(bps)-1)]
	temp_bp=total_bps
	temp_let=total_letters
	BlockCov={}
	for j in temp_let:
		BlockCov[j]=0
	Letter_Double={}
	Pair_ThroughBP=[]
	Double_Read_ThroughBP=[]
	Single_Read_ThroughBP=[]
	blackList=[]
	fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,bamChr,bps[0]-flank,bps[-1]+flank))
	#	chr_link={}
	#	num_of_reads=0
	while True:
		pbam=fbam.readline().strip().split()
		if not pbam: break
		if int(pbam[1])&4>0: continue
		if int(pbam[1])&1024>0:continue
		if int(pbam[1])&512>0:
			blackList.append(pbam[0])
			continue
		if not int(pbam[4])>QCAlign:
			continue
		if pbam[0] in blackList: continue
	#		num_of_reads+=1
		if int(pbam[1])&8>0 or not pbam[6]=='=':
			pos1=int(pbam[3])+low_qual_edge
			pos2=int(pbam[3])+cigar2reaadlength(pbam[5])-low_qual_edge
			block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
			block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
			if block1==block2:
				BlockCov[block1]+=cigar2reaadlength(pbam[5])
			else:
				rela_1=pos1-low_qual_edge-temp_bp[temp_let.index(block1)]
				rela_2=pos2+low_qual_edge-temp_bp[temp_let.index(block2)]
				Single_Read_ThroughBP.append([block1,rela_1,block2,rela_2,pbam[5]])
			if not pbam[6]=='=':
				if not pbam[0] in chr_link.keys():
					chr_link[pbam[0]]=[pbam[1:9]]
				else:
					chr_link[pbam[0]]+=[pbam[1:9]]
		elif int(pbam[1])&8==0:
			if pbam[6]=='=':
				if not pbam[0] in Letter_Double.keys():
					Letter_Double[pbam[0]]=[pbam[:9]]
				else:
					if not pbam[:9] in Letter_Double[pbam[0]]:
						Letter_Double[pbam[0]]+=[pbam[:9]]
						if int(Letter_Double[pbam[0]][0][3])<int(Letter_Double[pbam[0]][1][3]):
							pos1=int(Letter_Double[pbam[0]][0][3])+low_qual_edge
							pos2=int(Letter_Double[pbam[0]][1][3])+cigar2reaadlength(Letter_Double[pbam[0]][1][5])-low_qual_edge
						else:
							pos1=int(Letter_Double[pbam[0]][1][3])+low_qual_edge
							pos2=int(Letter_Double[pbam[0]][0][3])+cigar2reaadlength(Letter_Double[pbam[0]][0][5])-low_qual_edge
						block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
						block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
						if block1==block2:
							BlockCov[block1]+=cigar2reaadlength(Letter_Double[pbam[0]][0][5])
							BlockCov[block1]+=cigar2reaadlength(Letter_Double[pbam[0]][1][5])
							del Letter_Double[pbam[0]]
							blackList.append(pbam[0])
	fbam.close()
	for key in Letter_Double.keys():
		if key in blackList:
			del Letter_Double[key]
			continue
		if len(Letter_Double[key])==2:
			pos1=int(Letter_Double[key][0][3])
			pos2=int(Letter_Double[key][1][3])
			if not pos1>pos2:
				pos1=int(Letter_Double[key][0][3])
				pos1b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
				pos2=int(Letter_Double[key][1][3])
				pos2b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
				direct_temp=Reads_Direction_Detect(Letter_Double[key][0][1])
			elif pos1>pos2:
				pos1=int(Letter_Double[key][1][3])
				pos1b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
				pos2=int(Letter_Double[key][0][3])
				pos2b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
				direct_temp=Reads_Direction_Detect(Letter_Double[key][1][1])
			block1=Reads_block_assignment_1(temp_bp,temp_let,pos1+low_qual_edge)
			block2=Reads_block_assignment_1(temp_bp,temp_let,pos2+low_qual_edge)
			block1b=Reads_block_assignment_1(temp_bp,temp_let,pos1b-low_qual_edge)
			block2b=Reads_block_assignment_1(temp_bp,temp_let,pos2b-low_qual_edge)
			rela_1=pos1-temp_bp[temp_let.index(block1)]
			rela_2=pos2-temp_bp[temp_let.index(block2)]
			rela_1b=pos1b-temp_bp[temp_let.index(block1b)]
			rela_2b=pos2b-temp_bp[temp_let.index(block2b)]
			if block1==block1b and block2==block2b:
				Pair_ThroughBP.append([block1,rela_1,rela_1b, block2,rela_2,rela_2b]+direct_temp)
			else:
				Double_Read_ThroughBP.append([block1,rela_1,block1b,rela_1b, block2,rela_2,block2b,rela_2b]+direct_temp)
			del Letter_Double[key]
		elif len(Letter_Double[key])==1:
			if Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][7]))==0:
				if Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))==Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3])+cigar2reaadlength(Letter_Double[key][0][5])):
					BlockCov[Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))]+=cigar2reaadlength(Letter_Double[key][0][5])
					del Letter_Double[key]
	Initial_DR_Penal=0
	for j in Pair_ThroughBP:
		if not j[-2:]==['+', '-']:
			Initial_DR_Penal+=1
	for j in Double_Read_ThroughBP:
		if not j[-2:]==['+', '-']:
			Initial_DR_Penal+=1
	Initial_Cov={}
	for j in temp_let:
		Initial_Cov[j]=0
	for j in Pair_ThroughBP:
		Initial_Cov[j[0]]+=j[2]-j[1]
		Initial_Cov[j[3]]+=j[5]-j[4]
	for j in Single_Read_ThroughBP:
		Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
		Initial_Cov[j[2]]+=j[3]
	for j in Double_Read_ThroughBP:
		if j[0]==j[2]:
			Initial_Cov[j[0]]+=j[3]-j[1]
		else:
			Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
			Initial_Cov[j[2]]+=j[3]
		if j[4]==j[6]:
			Initial_Cov[j[4]]+=j[7]-j[5]
		else:
			Initial_Cov[j[4]]+=temp_bp[temp_let.index(j[4])+1]-temp_bp[temp_let.index(j[4])]-j[5]
			Initial_Cov[j[6]]+=j[7]
	Initial_IL=[]
	for j in Pair_ThroughBP:
		Initial_IL.append(temp_bp[temp_let.index(j[3])]-temp_bp[temp_let.index(j[0])]-j[1]+j[5])
	for j in Double_Read_ThroughBP:
		Initial_IL.append(temp_bp[temp_let.index(j[6])]-temp_bp[temp_let.index(j[0])]-j[1]+j[7])
	Initial_ILPenal=[]
	for j in Initial_IL:
		Initial_ILPenal+=[pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)/len(Initial_IL)]
	return [Initial_DR_Penal,Initial_ILPenal,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,BlockCov,Initial_Cov,Letter_Double]

def Single_Rec_Read_Locate(Letter_Double_rec,temp_bp, temp_let):
	Pair_ThroughBP=[]
	Double_Read_ThroughBP=[]
	Single_Read_ThroughBP=[]
	Initial_IL=[]
	BlockCov={}
	Initial_Cov={}
	Initial_DR_Penal=0
	for j in temp_let:
		BlockCov[j]=0
	for key in Letter_Double_rec.keys():
		if len(Letter_Double_rec[key])==1:
			pos1=int(Letter_Double_rec[key][0][3]) 
			pos2=int(Letter_Double_rec[key][0][7])
			bamChr=Letter_Double_rec[key][0][2]
			fbamtemp=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,bamChr,pos2,pos2+ReadLength))
			while True:
				pbam=fbamtemp.readline().strip().split()
				if not pbam: break
				flag=0
				if pbam[0]==key:
					Letter_Double_rec[key]+=[pbam[:9]]
					flag+=1
				if flag==1:
					break
			fbamtemp.close()		 
	for key in Letter_Double_rec.keys():
		if len(Letter_Double_rec[key])==2:
			pos1=int(Letter_Double_rec[key][0][3])
			pos2=int(Letter_Double_rec[key][1][3])
			if not pos1>pos2:
				pos1=int(Letter_Double_rec[key][0][3])
				pos1b=pos1+cigar2reaadlength(Letter_Double_rec[key][0][5])
				pos2=int(Letter_Double_rec[key][1][3])
				pos2b=pos2+cigar2reaadlength(Letter_Double_rec[key][1][5])
				direct_temp=Reads_Direction_Detect(Letter_Double_rec[key][0][1])
			elif pos1>pos2:
				pos1=int(Letter_Double_rec[key][1][3])
				pos1b=pos2+cigar2reaadlength(Letter_Double_rec[key][1][5])
				pos2=int(Letter_Double_rec[key][0][3])
				pos2b=pos1+cigar2reaadlength(Letter_Double_rec[key][0][5])
				direct_temp=Reads_Direction_Detect(Letter_Double_rec[key][1][1])
			if not pos1<temp_bp[0]-flank+1 and not pos2b>temp_bp[-1]+flank-1:
				block1=Reads_block_assignment_1(temp_bp,temp_let,pos1+low_qual_edge)
				block2=Reads_block_assignment_1(temp_bp,temp_let,pos2+low_qual_edge)
				block1b=Reads_block_assignment_1(temp_bp,temp_let,pos1b-low_qual_edge)
				block2b=Reads_block_assignment_1(temp_bp,temp_let,pos2b-low_qual_edge)
				rela_1=pos1-temp_bp[temp_let.index(block1)]
				rela_2=pos2-temp_bp[temp_let.index(block2)]
				rela_1b=pos1b-temp_bp[temp_let.index(block1b)]
				rela_2b=pos2b-temp_bp[temp_let.index(block2b)]
				if block1==block1b==block2==block2:
					BlockCov[block1]+=cigar2reaadlength(Letter_Double_rec[key][0][5])
				else:						
					if block1==block1b and block2==block2b:
						Pair_ThroughBP.append([block1,rela_1,rela_1b, block2,rela_2,rela_2b]+direct_temp)
					else:
						Double_Read_ThroughBP.append([block1,rela_1,block1b,rela_1b, block2,rela_2,block2b,rela_2b]+direct_temp)
				del Letter_Double_rec[key]
	for j in Pair_ThroughBP:
		if not j[-2:]==['+', '-']:
			Initial_DR_Penal+=1
	for j in Double_Read_ThroughBP:
		if not j[-2:]==['+', '-']:
			Initial_DR_Penal+=1
	for j in temp_let:
		Initial_Cov[j]=0
	for j in Pair_ThroughBP:
		Initial_Cov[j[0]]+=j[2]-j[1]
		Initial_Cov[j[3]]+=j[5]-j[4]
	for j in Single_Read_ThroughBP:
		Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
		Initial_Cov[j[2]]+=j[3]
	for j in Double_Read_ThroughBP:
		if j[0]==j[2]:
			Initial_Cov[j[0]]+=j[3]-j[1]
		else:
			Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
			Initial_Cov[j[2]]+=j[3]
		if j[4]==j[6]:
			Initial_Cov[j[4]]+=j[7]-j[5]
		else:
			Initial_Cov[j[4]]+=temp_bp[temp_let.index(j[4])+1]-temp_bp[temp_let.index(j[4])]-j[5]
			Initial_Cov[j[6]]+=j[7]
	Initial_IL=[]
	for j in Pair_ThroughBP:
		Initial_IL.append(temp_bp[temp_let.index(j[3])]-temp_bp[temp_let.index(j[0])]-j[1]+j[5])
	for j in Double_Read_ThroughBP:
		Initial_IL.append(temp_bp[temp_let.index(j[6])]-temp_bp[temp_let.index(j[0])]-j[1]+j[7])
	Initial_ILPenal=[]
	for j in Initial_IL:
		Initial_ILPenal+=[pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)/len(Initial_IL)]
	return [Initial_DR_Penal,Initial_ILPenal,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,BlockCov,Initial_Cov]

def letter_rearrange(bps2):
	chr_letter_bp={}
	let_start=96
	for i in bps2:
		if not i[0] in chr_letter_bp.keys():
			chr_letter_bp[i[0]]={}
		for j in range(len(i))[1:-1]:
			chr_letter_bp[i[0]][chr(let_start+j)]=[]
			if int(i[j+1])-int(i[j])<10*flank:
				chr_letter_bp[i[0]][chr(let_start+j)]+=[int(i[j]),int(i[j+1])]
			else:
				chr_letter_bp[i[0]][chr(let_start+j)]+=[int(i[j]),int(i[j])+flank,int(i[j+1])-flank,int(i[j+1])]
		let_start+=len(i)-2
	return chr_letter_bp

def letter_GC_ReadIn(chr_letter_bp):
	block_GC_temp={}
	filein=ref_prefix+'.GC_Content'
	block_range={}
	GC_hash={}
	test_flag=0
	for i in chr_letter_bp.keys():
		if not os.path.isfile(filein+'.'+str(i)):
			test_flag+=1
	if test_flag==0:
		for i in chr_letter_bp.keys():
			GC_hash[i]={}
			block_range[i]=[]
			for j in chr_letter_bp[i].keys():
				block_range[i]+=chr_letter_bp[i][j]
			block_range[i]=[min(block_range[i]),max(block_range[i])]
			fin=open(filein+'.'+str(i))
			while True:
				pin=fin.readline().strip().split()
				if not pin: break
				pin2=fin.readline().strip().split()
				if 'chr' in block_range.keys()[0]:
					if not 'chr' in pin[0]:
						pin[0]='chr'+pin[0]
				elif not 'chr' in block_range.keys()[0]:
					if 'chr' in pin[0]:
						pin[0]=pin[0][3:]
				if pin[0] in block_range.keys():
					if not int(pin[2])<block_range[pin[0]][0] and not int(pin[1])>block_range[pin[0]][1]:
						GC_hash[pin[0]][pin[1]+'-'+pin[2]]=pin2
			fin.close()
		for k1 in chr_letter_bp.keys():
			block_GC_temp[k1]={}
			for k2 in GC_hash[k1].keys():
				bl2=[int(k2.split('-')[0]),int(k2.split('-')[1])]
				for k3 in chr_letter_bp[k1].keys():
					if min(chr_letter_bp[k1][k3])>bl2[0]-1 and max(chr_letter_bp[k1][k3])<bl2[1]+1:
						block_GC_temp[k1][k3]=GC_hash[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/100:(max(chr_letter_bp[k1][k3])-bl2[0])/100+1]
					elif min(chr_letter_bp[k1][k3])>bl2[0]-1 and max(chr_letter_bp[k1][k3])>bl2[1]:
						if not k3 in block_GC_temp[k1].keys():
							block_GC_temp[k1][k3]=GC_hash[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/100:]
						else:
							block_GC_temp[k1][k3]+=GC_hash[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/100:]						
					elif min(chr_letter_bp[k1][k3])<bl2[0] and max(chr_letter_bp[k1][k3])>bl2[0]-1:
						if not k3 in block_GC_temp[k1].keys():
							block_GC_temp[k1][k3]=GC_hash[k1][k2][:(max(chr_letter_bp[k1][k3])-bl2[0])/100+1]
						else:
							block_GC_temp[k1][k3]+=GC_hash[k1][k2][:(max(chr_letter_bp[k1][k3])-bl2[0])/100+1]						
					elif min(chr_letter_bp[k1][k3])<bl2[0]+1 and max(chr_letter_bp[k1][k3])>bl2[1]-1:
						if not k3 in block_GC_temp[k1].keys():
							block_GC_temp[k1][k3]=GC_hash[k1][k2]
						else:
							block_GC_temp[k1][k3]+=GC_hash[k1][k2]					   
		for k1 in block_GC_temp.keys():
			for k2 in block_GC_temp[k1].keys():
				if not block_GC_temp[k1][k2]==[]:
					block_GC_temp[k1][k2]=numpy.mean([float(k3) for k3 in block_GC_temp[k1][k2]])
				else:
					return 'error'
		return block_GC_temp
	else:
		return 'error'

def letter_RD_ReadIn(chr_letter_bp):
	test_flag=0
	for k1 in chr_letter_bp.keys():
		filein=dict_opts['--ppre']+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/RD_Stat/'+BamN+'.'+k1+'.RD.index'
		if not os.path.isfile(filein):
			test_flag+=1
	if test_flag==0:
		out={}
		RD_hash={}
		block_range={}
		for i in chr_letter_bp.keys():
			RD_hash[i]={}
			out[i]={}
			block_range[i]=[]
			for j in chr_letter_bp[i].keys():
				block_range[i]+=chr_letter_bp[i][j]
			block_range[i]=[min(block_range[i]),max(block_range[i])]
		for k1 in chr_letter_bp.keys():
			filein=dict_opts['--ppre']+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/RD_Stat/'+BamN+'.'+k1+'.RD.index'
			fin=open(filein)
			while True:
				pin=fin.readline().strip().split()
				if not pin: break
				pin2=fin.readline().strip().split()
				bl2=[int(pin[0].split(':')[1].split('-')[0]),int(pin[0].split(':')[1].split('-')[1])]
				if not bl2[1]<block_range[k1][0]+1 and not bl2[0]>block_range[k1][1]-1:
					RD_hash[k1][str(bl2[0])+'-'+str(bl2[1])]=pin2
			fin.close()
		for k1 in chr_letter_bp.keys():
			for k2 in RD_hash[k1].keys():
				bl2=[int(k2.split('-')[0]),int(k2.split('-')[1])]
				for j in sorted(chr_letter_bp[k1].keys()):
					if not j in out[k1].keys():
						out[k1][j]=[]
					if len(chr_letter_bp[k1][j])==4:
						bl1=chr_letter_bp[k1][j][1:-1]
						if bl1[0]>bl2[0]-1 and bl1[1]<bl2[1]+1:
							out[k1][j]+=RD_hash[k1][k2][(bl1[0]-bl2[0])/100:(bl1[1]-bl2[0])/100+1]
						elif bl1[0]>bl2[0]-1 and bl1[1]>bl2[1]:
							out[k1][j]+=RD_hash[k1][k2][(bl1[0]-bl2[0])/100:]
						elif bl1[0]<bl2[0] and bl1[1]<bl2[1]+1:
							out[k1][j]+=RD_hash[k1][k2][:(bl1[1]-bl2[0])/100+1]
						elif bl1[0]<bl2[0] and bl1[1]>bl2[1]:
							out[k1][j]+=RD_hash[k1][k2]
		for k1 in out.keys():
			for k2 in out[k1].keys():
				if out[k1][k2]==[]:
					out[k1][k2]=0
				else:
					out[k1][k2]=numpy.mean([float(k3) for k3 in out[k1][k2]])
		return out
	else:
		return 'error'

def letter_range_report(chr_letter_bp):
	ks_range={}
	for k1 in chr_letter_bp.keys():
		ks_range[k1]=[]
		for k2 in chr_letter_bp[k1].keys():
			ks_range[k1]+=chr_letter_bp[k1][k2]
		ks_range[k1]=[min(ks_range[k1]),max(ks_range[k1])]
	for k1 in chr_letter_bp.keys():
		chr_letter_bp[k1]['left']=[ks_range[k1][0]-flank,ks_range[k1][0]]
		chr_letter_bp[k1]['right']=[ks_range[k1][1],ks_range[k1][1]+flank]

def block_Read_From_Bam(chr_letter_bp):
	blocks_out={}
	for k1 in chr_letter_bp.keys():
		k_block=[]
		blocks_out[k1]=[]
		for k2 in chr_letter_bp[k1].keys():
			if not k2 in ['left','right']:
				k_block.append(k2)
		blocks_out[k1].append(chr_letter_bp[k1]['left'][:]+['left'])
		k_block.sort()
		for k2 in k_block+['right']:
			if len(chr_letter_bp[k1][k2])==2:
				blocks_out[k1][-1]+=chr_letter_bp[k1][k2]+[k2]
			elif len(chr_letter_bp[k1][k2])==4:
				blocks_out[k1][-1]+=chr_letter_bp[k1][k2][:2]+[k2]
				blocks_out[k1].append(chr_letter_bp[k1][k2][2:]+[k2])
	for k1 in blocks_out.keys():
		for k2 in blocks_out[k1]:
			k2.sort()
	return blocks_out

tolerance_bp=10
def pos_block_assign(block_bps_chr,read_pos):
	read_new=[]
	for i in read_pos:
		if type(i)==type(1):
			read_new.append(i)
	for j in range(len(read_new)/2):
		read_new[2*j]=read_new[2*j]+tolerance_bp
		read_new[2*j+1]=read_new[2*j+1]-tolerance_bp
	for i in read_new:
		if type(i)==type(1):
			for j in block_bps_chr.keys():
				if j=='left':
					if i<block_bps_chr[j][1]:
						read_pos.append(j)
				elif j=='right':
					if i>block_bps_chr[j][0]-1:
						read_pos.append(j)
				else:
					if block_bps_chr[j][0]-1<i and block_bps_chr[j][1]>i:
						read_pos.append(j)

def block_Info_ReadIn(chr_letter_bp,blocks_read_in):
	block_bps={}
	block_rds={}
	for k1 in chr_letter_bp.keys():
		block_bps[k1]={}
		block_rds[k1]={}
		for k2 in chr_letter_bp[k1].keys():
			block_bps[k1][k2]=[min(chr_letter_bp[k1][k2]),max(chr_letter_bp[k1][k2])]
			block_rds[k1][k2]=0
	Pair_ThroughBP={}
	Double_Read_ThroughBP={}
	Single_Read_ThroughBP={}
	total_rec={}
	rd_low_qual={}
	for k1 in chr_letter_bp.keys():
		Pair_ThroughBP[k1]=[]
		Double_Read_ThroughBP[k1]=[]
		Single_Read_ThroughBP[k1]=[]
		rd_low_qual[k1]={}
		for k2 in blocks_read_in[k1]:
			k2a=[]
			k2b=[]
			for k3 in k2:
				if type(k3)==type(1):
					k2a.append(k3)
				else:
					k2b.append(k3)
			fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,k1,min(k2a)-flank,max(k2a)+flank))
			blackList=[]
			temp_rec={}	
			temp_rec_LowQual={}		
			while True:
				pbam=fbam.readline().strip().split()
				if not pbam: break
				if int(pbam[1])&4>0: continue
				if int(pbam[1])&1024>0:continue
				if int(pbam[1])&512>0:
					blackList.append(pbam[0])
					continue
				#if not int(pbam[4])>QCAlign:continue
				if pbam[0] in blackList: continue
				if not int(pbam[4])>QCAlign:
					if not pbam[0] in temp_rec_LowQual.keys():
						temp_rec_LowQual[pbam[0]]=[]
					if not pbam[1:9] in temp_rec_LowQual[pbam[0]]:
						temp_rec_LowQual[pbam[0]]+=[pbam[1:9]]
				else:
					if not pbam[0] in temp_rec.keys():
						temp_rec[pbam[0]]=[]
					if not pbam[1:9] in temp_rec[pbam[0]]:
						temp_rec[pbam[0]]+=[pbam[1:9]]
			fbam.close()
			flank_region=[]
			for k3 in k2b:
				flank_region+=block_bps[k1][k3]
			flank_region=[min(flank_region),max(flank_region)]
			for k3 in temp_rec_LowQual.keys():
				for k4 in temp_rec_LowQual[k3]:
					read_pos=[int(k4[2]),int(k4[2])+cigar2reaadlength(k4[4])]
					pos_block_assign(block_bps[k1],read_pos)
					if read_pos[-1]==read_pos[-2]:
						if not read_pos[-1] in rd_low_qual[k1].keys():
							rd_low_qual[k1][read_pos[-1]]=0
						rd_low_qual[k1][read_pos[-1]]+=(read_pos[1]-read_pos[0])
					else:
						if not read_pos[-2] in rd_low_qual[k1].keys():
							rd_low_qual[k1][read_pos[-2]]=0
						if not read_pos[-1] in rd_low_qual[k1].keys():
							rd_low_qual[k1][read_pos[-1]]=0
						rd_low_qual[k1][read_pos[-2]]+=block_bps[k1][read_pos[-2]][1]-read_pos[0]
						rd_low_qual[k1][read_pos[-1]]+=-block_bps[k1][read_pos[-1]][0]+read_pos[1]
			for k3 in temp_rec.keys():
				if len(temp_rec[k3])>2:
					test_rec=[int(temp_rec[k3][0][7])]
					test_rec2=[temp_rec[k3][0]]
					test_let=0
					for k4 in temp_rec[k3][1:]:
						delflag=0
						for k5 in test_rec:
							if int(k4[7])+k5==0:
								test_let+=1
								k6=k3+chr(96+test_let)
								temp_rec[k6]=[test_rec2[test_rec.index(k5)],k4]
								del test_rec2[test_rec.index(k5)]
								del test_rec[test_rec.index(k5)]
								delflag+=1
						if delflag==0:
							test_rec.append(int(k4[7]))
							test_rec2.append(k4)
					temp_rec[k3]=test_rec2 
			for k3 in temp_rec.keys():
				if len(temp_rec[k3])==1:
					del_flag=0
					k4=temp_rec[k3][0]
					read_pos=[int(k4[2]),int(k4[2])+cigar2reaadlength(k4[4])]
					mate_pos=[int(k4[6]),int(k4[6])+ReadLength]
					if 'left' in k2b and mate_pos[1]<flank_region[0]:
						del_flag+=1
					elif 'right' in k2b and mate_pos[0]>flank_region[0]:
						del_flag+=1
					#elif not mate_pos[1]<flank_region[0] and not mate_pos[0]>flank_region[1]:
					#	del_flag+=1
					if del_flag>0:
						del temp_rec[k3]
						pos_block_assign(block_bps[k1],read_pos)
						if read_pos[-1]==read_pos[-2]:
							block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
						else:
							Single_Read_ThroughBP[k1].append(read_pos)
					else:
						if not k3 in total_rec.keys():
							total_rec[k3]=[k4]
						else:
							total_rec[k3]+=[k4]
				elif len(temp_rec[k3])==2: 
					if int(temp_rec[k3][0][7])==0 or int(temp_rec[k3][1][7])==0:
						continue
					if int(temp_rec[k3][0][7])+int(temp_rec[k3][1][7])==0 and int(temp_rec[k3][0][7])<0:
						temp_rec[k3]=[temp_rec[k3][1],temp_rec[k3][0]]
					read_pos=[int(temp_rec[k3][0][2]),int(temp_rec[k3][0][2])+cigar2reaadlength(temp_rec[k3][0][4]),int(temp_rec[k3][1][2]),int(temp_rec[k3][1][2])+cigar2reaadlength(temp_rec[k3][1][4])]+Reads_Direction_Detect(temp_rec[k3][0][0])
					#print temp_rec[k3]
					#if k3 in test2:
					#	print read_pos
					if read_pos[0]>read_pos[2]:
						read_pos=read_pos[2:4]+read_pos[:2]+[read_pos[-1],read_pos[-2]]
					pos_block_assign(block_bps[k1],read_pos)
					if read_pos[6]==read_pos[7]==read_pos[8]==read_pos[9]:
						block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
						block_rds[k1][read_pos[-1]]+=read_pos[3]-read_pos[2]
					elif read_pos[8]==read_pos[9] and read_pos[6]==read_pos[7]:
						Pair_ThroughBP[k1].append(read_pos[:6]+[read_pos[6],read_pos[8]])
					else:
						Double_Read_ThroughBP[k1].append(read_pos)
					del temp_rec[k3]
					#if k3 in test2:
					#	print read_pos
	for k3 in total_rec.keys():
		if len(total_rec[k3])==1: 
			del_flag=0
			k4=total_rec[k3][0]
			read_pos=[int(k4[2]),int(k4[2])+cigar2reaadlength(k4[4])]
			mate_pos=[int(k4[6]),int(k4[6])+ReadLength]
			if 'left' in k2b and mate_pos[1]<flank_region[0]:
				del_flag+=1
			elif 'right' in k2b and mate_pos[0]>flank_region[0]:
				del_flag+=1
			elif not mate_pos[1]<flank_region[0] and not mate_pos[0]>flank_region[1]:
				del_flag+=1
			if del_flag>0:
				del total_rec[k3]
				pos_block_assign(block_bps[k1],read_pos)
				if read_pos[-1]==read_pos[-2]:
					block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
				else:
					Single_Read_ThroughBP[k1].append(read_pos)
		elif len(total_rec[k3])==2:
			read_pos=[int(total_rec[k3][0][2]),int(total_rec[k3][0][2])+cigar2reaadlength(total_rec[k3][0][4]),int(total_rec[k3][1][2]),int(total_rec[k3][1][2])+cigar2reaadlength(total_rec[k3][1][4])]+Reads_Direction_Detect(total_rec[k3][0][0])
			#print read_pos
			if read_pos[0]>read_pos[2]:
				read_pos=read_pos[2:4]+read_pos[:2]+[read_pos[-1],read_pos[-2]]
			pos_block_assign(block_bps[k1],read_pos)
			if read_pos[6]==read_pos[7]==read_pos[8]==read_pos[9]:
				block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
				block_rds[k1][read_pos[-1]]+=read_pos[3]-read_pos[2]
			elif read_pos[8]==read_pos[9] and read_pos[6]==read_pos[7]:
				Pair_ThroughBP[k1].append(read_pos[:6]+[read_pos[6],read_pos[8]])
			else:
				Double_Read_ThroughBP[k1].append(read_pos)
			del total_rec[k3]
	#print total_rec
	direction_penal=0
	block_rd2={}
	block_bp2=block_bps
	for k1 in block_rds.keys():
		block_rd2[k1]={}
		for k2 in block_rds[k1].keys():
			block_rd2[k1][k2]=0
	for i2 in Pair_ThroughBP.keys():
		for i in Pair_ThroughBP[i2]:
			if not i[4:6]==['+','-']:
				direction_penal+=1
			block_rd2[i2][i[6]]+=i[1]-i[0]
			block_rd2[i2][i[7]]+=i[3]-i[2]
	for i2 in Double_Read_ThroughBP.keys():
		for i in Double_Read_ThroughBP[i2]:
			if i[6]==i[7]:
				block_rd2[i2][i[6]]+=i[1]-i[0]
				block_rd2[i2][i[8]]+=-i[2]+block_bp2[i2][i[8]][1]
				block_rd2[i2][i[9]]+=i[3]-block_bp2[i2][i[9]][0]
				#if -i[2]+block_bp2[i2][i[8]][1]>200 and i[8]=='a':
					#print i
				#if i[3]-block_bp2[i2][i[9]][0]>200 and i[9]=='a':
					#print i
			elif i[8]==i[9]:
				block_rd2[i2][i[8]]+=i[3]-i[2]
				block_rd2[i2][i[6]]+=-i[0]+block_bp2[i2][i[6]][1]
				block_rd2[i2][i[7]]+=i[1]-block_bp2[i2][i[7]][0]
				#if -i[0]+block_bp2[i2][i[6]][1]>101:
					#print i
				#if i[1]-block_bp2[i2][i[7]][0]>101:
					#print i
			else:
				block_rd2[i2][i[6]]+=-i[0]+block_bp2[i2][i[6]][1]
				block_rd2[i2][i[7]]+=i[1]-block_bp2[i2][i[7]][0]
				block_rd2[i2][i[8]]+=-i[2]+block_bp2[i2][i[8]][1]
				block_rd2[i2][i[9]]+=i[3]-block_bp2[i2][i[9]][0]
	for i2 in Single_Read_ThroughBP.keys():
		for i in Single_Read_ThroughBP[i2]:
			block_rd2[i2][i[2]]+=-i[0]+block_bp2[i2][i[2]][1]
			block_rd2[i2][i[3]]+=i[1]-block_bp2[i2][i[3]][0]
	for k1 in rd_low_qual.keys():
		for k2 in rd_low_qual[k1].keys():
			block_rds[k1][k2]+=rd_low_qual[k1][k2]
	return [block_rds,block_rd2,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP]

def total_rd_calcu(letter_RD2,letter_GC,chr_letter_bp,block_rd2):
	out={}
	for k1 in block_rd2.keys():
		for k2 in block_rd2[k1].keys():
			if not k2 in ['left','right']:
				out[k2]=[]
	for k2 in letter_RD2.keys():
		out[k2].append(letter_RD2[k2])
	for k1 in block_rd2.keys():
		for k2 in block_rd2[k1].keys():
			if not k2 in ['left','right']:
				if not chr_letter_bp[k1][k2][-1]==chr_letter_bp[k1][k2][0]:
					out[k2].append(float(block_rd2[k1][k2])/float(chr_letter_bp[k1][k2][-1]-chr_letter_bp[k1][k2][0]))
	for k2 in out.keys():
		out[k2]=numpy.sum(out[k2])
	GC2={}
	for k1 in letter_GC.keys():
		for k2 in letter_GC[k1].keys():
			if not k2 in ['left','right']:
				GC2[k2]=letter_GC[k1][k2]
	out2=GC_RD_Adj_hash(GC_Median_Num,GC_Overall_Median_Num,k1,GC2,out)
	return out2

def GC_RD_Adj_hash(GC_Median_Num,GC_Overall_Median_Num,Chromo,GC_Content,Coverage):
	Coverage_af_Adj={}
	Overall_Median_Coverage=float(GC_Overall_Median_Num)
	for key_1 in GC_Content.keys():
			Be_Adj_RD=Coverage[key_1]
			GC_Con=int(round(GC_Content[key_1]*100))
			if GC_Con in GC_Median_Num.keys():
					Median_Coverage=GC_Median_Num[GC_Con]
					Af_Adj_RD=Be_Adj_RD*Overall_Median_Coverage/Median_Coverage
			elif not GC_Con in GC_Median_Num.keys():
					Af_Adj_RD=Overall_Median_Coverage
			Coverage_af_Adj[key_1]=Af_Adj_RD
	return 	Coverage_af_Adj

def DR_Penal_Calcu(read_info):
	DR_Penal=0
	for i2 in read_info[2].keys():
		for i in read_info[2][i2]:
			if i[4:6]==['+', '+'] or i[4:6]==['-', '-'] or i[4:6]==['-', '+']:
				DR_Penal+=1
	for i2 in read_info[3].keys():
		for i in read_info[3][i2]:
			if i[4:6]==['+', '+'] or i[4:6]==['-', '-'] or i[4:6]==['-', '+']:
				DR_Penal+=1
	return DR_Penal

def IL_Penal_Calcu(read_info):
	Initial_IL=[]
	for j2 in read_info[2].keys():
		for j in read_info[2][j2]:
			Initial_IL.append(j[3]-j[1])
	for j2 in read_info[3].keys():
		for j in read_info[3][j2]:
			Initial_IL.append(j[3]-j[1])
	Initial_ILPenal=[]
	for j in Initial_IL:
		Initial_ILPenal+=[pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)/len(Initial_IL)]
	return Initial_ILPenal

def original_bp_let_produce(chr_letter_bp,bps2):
	all_bp=[int(i) for i in bps2[0][1:]]
	all_let=[chr(97+i) for i in range(len(all_bp)-1)]
	origin=int(bps2[0][1])
	for i in bps2[1:]:
		temp_ks=[]
		for k2 in chr_letter_bp[i[0]].keys():
			if not k2 in ['left','right']:
				temp_ks.append(k2)
		for k2 in sorted(temp_ks):
			all_let.append(k2)
			all_bp.append(all_bp[-1]+chr_letter_bp[i[0]][k2][-1]-chr_letter_bp[i[0]][k2][0])
	return [all_bp,all_let]	

def rela_Pair_ThroughBP(chr_letter_bp,Pair_ThroughBP):
	out=[]
	for k1 in Pair_ThroughBP.keys():
		for k2 in Pair_ThroughBP[k1]:
			rela=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
						k2[1]-chr_letter_bp[k1][k2[6]][0],
				  k2[7],k2[2]-chr_letter_bp[k1][k2[7]][0],
						k2[3]-chr_letter_bp[k1][k2[7]][0],k2[4],k2[5]]
			out.append(rela)
	return out 

def rela_Pair_Double_Read_ThroughBP(chr_letter_bp,Double_Read_ThroughBP):
	out=[]
	for k1 in Double_Read_ThroughBP.keys():
		 for k2 in Double_Read_ThroughBP[k1]:
			rela=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
				  k2[7],k2[1]-chr_letter_bp[k1][k2[7]][0],
				  k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
				  k2[9],k2[3]-chr_letter_bp[k1][k2[9]][0],k2[4],k2[5]]
			out.append(rela)
	return out

def read_Pair_Single_Read_ThroughBP(chr_letter_bp,Single_Read_ThroughBP):
	out=[]
	for k1 in Single_Read_ThroughBP.keys():
		for k2 in Single_Read_ThroughBP[k1]:
			rela=[k2[2],k2[0]-chr_letter_bp[k1][k2[2]][0],
				  k2[3],k2[1]-chr_letter_bp[k1][k2[3]][0]]
			out.append(rela)
	return out

def Full_Info_of_Reads_Integrate(bps2):
	chr_letter_bp=letter_rearrange(bps2)
	letter_GC=letter_GC_ReadIn(chr_letter_bp)
	#print 'GC Content Read In'
	letter_RD=letter_RD_ReadIn(chr_letter_bp)
	#print 'Read Depth Read In'
	letter_range_report(chr_letter_bp)
	blocks_read_in=block_Read_From_Bam(chr_letter_bp)
	read_info=block_Info_ReadIn(chr_letter_bp,blocks_read_in)
	#print 'Read Info Read In'
	block_rds=read_info[0]
	block_rd2=read_info[1]
	letter_RD2={}
	for k1 in letter_RD.keys():
		for k2 in letter_RD[k1].keys():
				if len(chr_letter_bp[k1][k2])==4:
						letter_RD2[k2]=letter_RD[k1][k2]*(chr_letter_bp[k1][k2][2]-chr_letter_bp[k1][k2][1])/(chr_letter_bp[k1][k2][3]-chr_letter_bp[k1][k2][0])
				else:
						letter_RD2[k2]=letter_RD[k1][k2]
	for k1 in block_rds.keys():
		for k2 in block_rds[k1].keys():
				if not k2 in ['left','right']:
					if not chr_letter_bp[k1][k2][-1]==chr_letter_bp[k1][k2][0]:
						letter_RD2[k2]+=float(block_rds[k1][k2])/float(chr_letter_bp[k1][k2][-1]-chr_letter_bp[k1][k2][0])
	Pair_ThroughBP=rela_Pair_ThroughBP(chr_letter_bp,read_info[2])
	Double_Read_ThroughBP=rela_Pair_Double_Read_ThroughBP(chr_letter_bp,read_info[3])
	Single_Read_ThroughBP=read_Pair_Single_Read_ThroughBP(chr_letter_bp,read_info[4])
	Initial_RD=total_rd_calcu(letter_RD2,letter_GC,chr_letter_bp,block_rd2)
	DR_Penal=DR_Penal_Calcu(read_info)
	IL_Penal=IL_Penal_Calcu(read_info)
	letter_GC_out={}
	for k1 in letter_GC.keys():
			for k2 in letter_GC[k1].keys():
					letter_GC_out[k2]=letter_GC[k1][k2]
	return [letter_RD2,Initial_RD,DR_Penal,numpy.mean(IL_Penal),Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,letter_GC_out]+original_bp_let_produce(chr_letter_bp,bps2)

def GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,Chromo,GC_Content,Coverage):
	Coverage_af_Adj=[]
	GC_Content=[float(i) for i in GC_Content]
	Coverage=[float(i) for i in Coverage]
	Overall_Median_Coverage=float(GC_Overall_Median_Num)
	for key_1 in range(len(GC_Content)):
			Be_Adj_RD=Coverage[key_1]
			GC_Con=GC_Content[key_1]
			GC_Con=int(round(GC_Con*100))
			if GC_Con in GC_Median_Num.keys():
					Median_Coverage=GC_Median_Num[GC_Con]
					Af_Adj_RD=Be_Adj_RD*Overall_Median_Coverage/Median_Coverage
			elif not GC_Con in GC_Median_Num.keys():
					Af_Adj_RD=Overall_Median_Coverage
			Coverage_af_Adj+=[Af_Adj_RD]
	return Coverage_af_Adj

def GC_RD_Correction(chrbam):
	ref_index=ref_prefix+'GC_Content'
	cov_index=dict_opts['--ppre']+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/RD_Stat/'+BamN+'.'+chrbam+'.RD.index'
	fref=open(ref_index)
	fcov=open(cov_index)
	GC=[]
	RD=[]
	while True:
		pref1=fref.readline().strip().split()
		if not pref1: break
		pref2=fref.readline().strip().split()
		if pref1[0]==chrbam: 
			GC.append([int(pref1[1]),int(pref1[2])]+pref2)
	while True:
		pcov1=fcov.readline().strip().split()
		if not pcov1: break
		pcov2=fcov.readline().strip().split()
		reg1=int(pcov1[0].split(':')[1].split('-')[0])
		reg2=int(pcov1[0].split(':')[1].split('-')[1])
		RD.append([reg1,reg2]+pcov2)
	fref.close()
	fcov.close()
	rd_rec=0
	RD2={}
	for i in GC:
		RD2[i[0]]=[]
		region1=i[:2]
		while True:
			if rd_rec==len(RD): break
			region2=RD[rd_rec]
			if region2[1]-1<region1[1]:
				RD2[i[0]].append(region2[:2])
				rd_rec+=1
			else:
				break
	RD3=[]
	rd_rec=-1
	for i in GC:
		rd_rec+=1
		RD3.append(i[:2])
		j=RD2[i[0]][0]
		j2=RD[rd_rec]
		RD3[-1]+=j2[2:]		
		for j in RD2[i[0]][1:]:
			rd_rec+=1
			j2=RD[rd_rec]
			RD3[-1]+=j2[12:]
	RD4=[]
	for i in range(len(RD3)):
		RD4.append(RD3[i][:2])
		RD4[-1]+=GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,chrbam,GC[i][2:],RD3[i][2:])
	return [GC,RD3,RD4]

def GC_index(chrbam):
	ref_index=ref_prefix+'GC_Content'
	fref=open(ref_index)
	Cov=[]
	RD=[]
	while True:
		pref1=fref.readline().strip().split()
		if not pref1: break
		pref2=fref.readline().strip().split()
		if pref1[0]==chrbam: 
			Cov.append([int(pref1[1]),int(pref1[2])]+pref2)
	return Cov

def Reads_Direction_Detect(flag):
#flag is the number on the second position of each read in bam file
	#	if int(pi[1])&2>0: #two both mapped
	#	elif int(pi[1])&4>0: #the read itself mapped, mate not
	#	elif int(pi[1])&8>0: #mate mapped, the read itself not
	flag2=int(flag)
	if int(flag2)&16==0: 
			direct_1='+'
	elif int(flag2)&16>0:
			direct_1='-'
	if int(flag2)&32==0:
			direct_2='+'
	elif int(flag2)&32>0: 
			direct_2='-'
	if flag2&8==0:
		return([direct_1,direct_2])
	else:
		return([direct_1,'0'])

def Reads_block_assignment(bps,letters,block):
	if numpy.mean(block[:2])<bps[0]:
		return 'left'
	elif not numpy.mean(block[:2])<bps[-1]:
		return 'right'
	else:
		for i in letters:
			if not numpy.mean(block[:2])<bps[letters.index(i)] and numpy.mean(block[:2])<bps[letters.index(i)+1]:
				return i

def clusterNums(data, ClusterLen, direction):
	data.sort()
	out2=[]
	out3=[]
	out=[]
	for i in data:
		if len(out)==0:
			out.append(i)
		else:
			if i-out[-1]< ClusterLen:
				out.append(i)
			else:
				if out[-1]-out[0]<ClusterLen:
					if direction=='f':
						out2.append(out[-1])
						out3.append(len(out))
					elif direction=='r':
						out2.append(out[0])
						out3.append(len(out))
				else:
					temp=[0]
					lenTime=int((out[-1]-out[0])/ClusterLen)
					lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
					lenInter2=sorted(lenInter)[::-1][:lenTime]
					for j in range(len(lenInter)):
						if lenInter[j] in  lenInter2:
							temp.append(j+1)
					temp.append(len(out))
					if direction=='f':
						for k in range(len(temp)-1):
							out2.append(out[temp[k]:temp[k+1]][-1])
							out3.append(temp[k+1]-temp[k])
					elif direction=='r':
						for k in range(len(temp)-1):
							out2.append(out[temp[k]:temp[k+1]][0])
							out3.append(temp[k+1]-temp[k])
				out=[i]
	if out[-1]-out[0]<ClusterLen:
		if direction=='f':
			out2.append(out[-1])
			out3.append(len(out))
		elif direction=='r':
			out2.append(out[0])
			out3.append(len(out))
	else:
		temp=[0]
		lenTime=int((out[-1]-out[0])/ClusterLen)
		lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
		lenInter2=sorted(lenInter)[::-1][:lenTime]
		for j in range(len(lenInter)):
			if lenInter[j] in  lenInter2:
				temp.append(j+1)
		temp.append(len(out))
		if direction=='f':
			for k in range(len(temp)-1):
				out2.append(out[temp[k]:temp[k+1]][-1])
				out3.append(temp[k+1]-temp[k])
		elif direction=='r':
			for k in range(len(temp)-1):
				out2.append(out[temp[k]:temp[k+1]][0])
				out3.append(temp[k+1]-temp[k])
	return [out2,out3]

def clusterQC(hash, CluNum):
	out=[]
	for i in range(len(hash[1])):
		if not hash[1][i]<CluNum:
			out.append(hash[0][i])
	return out

def Single_Read_Assort_For_insert(Full_Info,bp_list,flank):
	relative_bps=[i-bp_list[0] for i in bp_list]
	letter_list=[chr(97+i) for i in range(len(bp_list)-1)]
	Block_and_Reads={}
	Block_and_Reads['left']=[]
	Block_and_Reads['right']=[]
	SingleR_Through=Full_Info[6]
	Pair_Through=Full_Info[4]
	Read_Through=Full_Info[5]
	for block in letter_list:
			Block_and_Reads[block]=[]
	for j in Pair_Through:
		Block_and_Reads[j[0]]=[j[1:3],j[3:]]
		Block_and_Reads[j[3]]=[j[4:6],j[:3]+j[6:8]]
	for j in Read_Through:
		Block_and_Reads[j[0]]=[]
	for key in Full_Info_of_Reads.keys():
			read_left=[int(i) for i in Full_Info_of_Reads[key][:2]]+[Full_Info_of_Reads[key][-2]]
			read_right=[int(i) for i in Full_Info_of_Reads[key][2:4]]+[Full_Info_of_Reads[key][-1]]
			assign_left=Reads_block_assignment_2(relative_bps,letter_list,read_left[0],read_left[1],flank)
			assign_right=Reads_block_assignment_2(relative_bps,letter_list,read_right[0],read_right[1],flank)
			New_Info=['_'.join([assign_left[0],str(int(co)-assign_left[1])]) for co in Full_Info_of_Reads[key][:2]]+['_'.join([assign_right[0],str(int(co)-assign_right[1])]) for co in Full_Info_of_Reads[key][2:4]]+Full_Info_of_Reads[key][4:]
			Block_and_Reads[assign_left[0]][key]=New_Info
			Block_and_Reads[assign_right[0]][key]=New_Info
	return Block_and_Reads

def complementary(seq):
	seq2=[]
	for i in seq:
		if i in 'ATGCN':
			seq2.append('ATGCN'['TACGN'.index(i)])
		elif i in 'atgcn':
			seq2.append('atgcn'['tacgn'.index(i)])
	return ''.join(seq2)		

def Insert_Seq_Pool_Prod_2(original_bp_list,ori_1_Seq,flank):
	ini_letters=['left']+['I'+chr(97+i) for i in range(len(original_bp_list)-1)]+['right']+['I'+chr(97+i)+'^' for i in range(len(original_bp_list)-1)]
	relative_bps=[0]+[j-original_bp_list[0]+flank for j in original_bp_list]+[original_bp_list[-1]+flank-original_bp_list[0]+flank]
	Insert_Seq_Pool={}
	for k in range(len(original_bp_list)+1):
		Insert_Seq_Pool[ini_letters[k]]=ori_1_Seq[relative_bps[k]:relative_bps[k+1]]
	for k in range(len(original_bp_list)+1,len(ini_letters)):
		Insert_Seq_Pool[ini_letters[k]]=complementary(ori_1_Seq[relative_bps[k-len(original_bp_list)]:relative_bps[k+1-len(original_bp_list)]])
	return Insert_Seq_Pool

def Median_Pick(number_list):
	#number_list: a list of number
	#Output: median of all numbers in the list
	if 2*(len(number_list)/2)+1==len(number_list):
		return float(sorted(number_list)[len(number_list)/2])
	elif 2*(len(number_list)/2)==len(number_list):
		return float(sorted(number_list)[len(number_list)/2-1]+sorted(number_list)[len(number_list)/2])/float(2)

def All_Reads_IL_Score_2a(Insert_Len_Stat,Letter_Double,Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero):
	#This function is used to calculate the IL_Penalty for Letter_Double only
	IL_Statistics=IL_Stat(Insert_Len_Stat)
	Before_Length=[]
	for key in Letter_Double.keys():
		LETT=Letter_Double[key]
		for lett in LETT:
			Before_Length.append(numpy.max([int(i) for i in lett[:4]])-numpy.min([int(i) for i in lett[:4]]))
	score=0
	for j in Before_Length:
		score+=pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
	return score

def All_Reads_IL_Score_2b(Insert_Len_Stat,Be_Letter_Through,Af_Letter_Through,Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero):
	#This function is used to calculate the IL_Penalty for Letter_Through only
	IL_Statistics=IL_Stat(Insert_Len_Stat)
	Before_Length=[]
	After_Length=[]
	for key in Be_Letter_Through.keys():
		if not key in Af_Letter_Through.keys(): continue
		else:
			Before_Info=Be_Letter_Through[key]
			After_Info=Af_Letter_Through[key]
			Before_Length.append(numpy.max([int(i) for i in Before_Info[:4]])-numpy.min([int(i) for i in Before_Info[:4]]))
			After_Length.append(numpy.max([int(i) for i in After_Info[:4]])-numpy.min([int(i) for i in After_Info[:4]]))
	score=0
	for j in After_Length:
		score+=pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
	return score

def All_Reads_IL_Score_2d(Insert_Len_Stat,Letter_Through,Af_Letter_Through,Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero):
	#take the average of multi mapped reads
	IL_Statistics=IL_Stat(Insert_Len_Stat)
	After_Length=[]
	score=0
	keys={}
	for key1 in Af_Letter_Through.keys():
		if not '_'.join(key1.split('_')[:-1]) in Letter_Through.keys():
			keys[key1]=[key1]
		else:
			if not '_'.join(key1.split('_')[:-1]) in keys.keys():
				keys['_'.join(key1.split('_')[:-1])]=[key1]
			else:
				keys['_'.join(key1.split('_')[:-1])].append(key1)
	for key2 in keys.keys():
		if len(keys[key2])==1:
			key3=keys[key2][0]
			After_Info=Af_Letter_Through[key3]
			After_Length.append(numpy.max([int(i) for i in After_Info[:4]])-numpy.min([int(i) for i in After_Info[:4]]))
			score+=pdf_calculate(After_Length[-1],IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
		else:
			temp_pdf=[]
			temp_Length=[]
			for key3 in keys[key2]:
				After_Info=Af_Letter_Through[key3]
				temp_Length.append(numpy.max([int(i) for i in After_Info[:4]])-numpy.min([int(i) for i in After_Info[:4]]))
				temp_pdf.append(pdf_calculate(temp_Length[-1],IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero))
			score+=sum([i/len(temp_pdf)for i in temp_pdf])
	return score

def All_Reads_IL_Score_2e(Insert_Len_Stat,Letter_Through,Af_Letter_Through,Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero):
#take the best of multi mapped reads
	IL_Statistics=IL_Stat(Insert_Len_Stat)
	After_Length=[]
	score=0
	keys={}
	for key1 in Af_Letter_Through.keys():
		if not '_'.join(key1.split('_')[:-1]) in Letter_Through.keys():
			keys[key1]=[key1]
		else:
			if not '_'.join(key1.split('_')[:-1]) in keys.keys():
				keys['_'.join(key1.split('_')[:-1])]=[key1]
			else:
				keys['_'.join(key1.split('_')[:-1])].append(key1)
	for key2 in keys.keys():
		if len(keys[key2])==1:
			key3=keys[key2][0]
			After_Info=Af_Letter_Through[key3]
			if not numpy.max([int(i) for i in After_Info[:4]])-numpy.min([int(i) for i in After_Info[:4]])>1000:
				After_Length.append(numpy.max([int(i) for i in After_Info[:4]])-numpy.min([int(i) for i in After_Info[:4]]))
			else:
				After_Length.append(1000)
			score+=pdf_calculate(After_Length[-1],IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
		else:
			temp_pdf=[]
			temp_Length=[]
			for key3 in keys[key2]:
				After_Info=Af_Letter_Through[key3]
				temp_Length.append(numpy.max([int(i) for i in After_Info[:4]])-numpy.min([int(i) for i in After_Info[:4]]))
				temp_pdf.append(pdf_calculate(temp_Length[-1],IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero))
			score+=numpy.max(temp_pdf)
	return score

def Direction_Penalize_2a(Coordinate_Info):		
	#This function is used to calculate the Direction_Penalize for Letter_Double only
	#Eg of Coordinate_Info: hash containing many: 'HWI-ST177_136:2:47:18313:117413': ['1032', '1132', '1407', '1507', '+', '-'], 
	Direction_Penalty=0
	for key in Coordinate_Info.keys():
		LETT=Coordinate_Info[key]
		for lett in LETT:
			if lett[4]==lett[5]:
				Direction_Penalty+=1
			else: continue
	return Direction_Penalty

def Direction_Penalize_2b(Coordinate_Info):		
	#This function is used to calculate the Direction_Penalize for Letter_Double only
	#Eg of Coordinate_Info: hash containing many: 'HWI-ST177_136:2:47:18313:117413': ['1032', '1132', '1407', '1507', '+', '-'], 
	Direction_Penalty=0
	for key in Coordinate_Info.keys():
		lett=Coordinate_Info[key]
		if lett[4]==lett[5]:
			Direction_Penalty+=1
		else: continue
	return Direction_Penalty

def Prob_Possion(Number, Mean):
	#This function is used to calculate the probability of RD of each region
	#a possion model is applied
	#prob is produced in a log scale
	#float is transferred to the closest integer
	if not Mean==0:
		probs=[]
		from scipy.stats import poisson
		lamda=Mean
		K=Number
		K2=Float_2_Integer(K)
		pdf=K2*math.log(lamda)-lamda-math.log(math.factorial(K2))
		return pdf
	elif Mean==0:
		return -Number

def Prob_NB(Number, Mean, Variance):
	if not Mean==0:
		P=1-Mean/Variance
		R=Float_2_Integer(Mean*(1-P)/P)
		K2=Float_2_Integer(Number)
		if P>0 and R>1 and K2+R-1>0 :
			log_pdf=math.log(math.factorial(K2+R-1))-math.log(math.factorial(K2))-math.log(math.factorial(R-1))+R*math.log(1-P)+K2*math.log(P)
		else:
			log_pdf=Prob_Possion(Number, Mean)
		return log_pdf
	elif Mean==0:
		return -Number

def Prob_Norm(Number, Mean, Variance):
	return math.log(1/sqrt(2*numpy.pi*Variance)*exp(-(Number-Mean)**2/2/Variance))

def Block_Assign_To_Letters(bp_list,letter_list,flank):
	#Eg of bp_list:[184569179, 184569775, 184571064, 184572009, 184572016]
	#Eg of letter_list:['a', 'b', 'c', 'd']
	#Eg of flank:446
	number_of_blocks=(numpy.max(bp_list)-numpy.min(bp_list)+2*flank)/100+1
	blocks={}
	bp_list_new=[bp_list[0]-flank]+bp_list+[bp_list[-1]+flank]
	relative_bp_list=[i-numpy.min(bp_list_new) for i in bp_list_new]
	bp_length=[(bp_list_new[i+1]-bp_list_new[i]) for i in range(len(bp_list_new)-1)]
	letter_list_new=['left']+letter_list+['right']
	bp_blocks=[[letter_list_new[j]]+range(relative_bp_list[j]/100,relative_bp_list[j+1]/100+1) for j in range(len(relative_bp_list)-1)]
	blocks_bp={}
	for i in range(number_of_blocks):
		blocks_bp[i+1]=[bp_list_new[0]+i*100,bp_list_new[0]+i*100+99]
		for j in bp_blocks:
			if i in j:
				blocks_bp[i+1].append(j[0])
	blocks_bp[0]=[blocks_bp[1][0]-100,blocks_bp[1][0]-1,'0']
	blocks_bp[number_of_blocks+1]=[blocks_bp[number_of_blocks][1]+1,blocks_bp[number_of_blocks][1]+100,'0']
	return blocks_bp

def GC_Content_Calculate(seq2):
	NumAT=0
	NumGC=0
	for n in seq2:
		if n=='A' or n=='T' or n=='a' or n=='t':
			NumAT+=1
		elif n=='G' or n=='C' or n=='g' or n=='c':
			NumGC+=1
	return [NumGC,NumAT]			

def Float_2_Integer(Number):
	if Number-int(Number)<0.5:
		K2=int(Number)
	elif not Number-int(Number)<0.5:
		K2=int(Number)+1
	return K2

def c_Coverage_Calculate_InfoList(Full_Info,Chromo,bp_MP,letter_MP,original_bp_list,flank):
	bp_M=[i-original_bp_list[0] for i in bp_MP[0]]
	bp_P=[i-original_bp_list[0] for i in bp_MP[1]]
	M_New_bp=[bp_M[0]-flank]+bp_M+[bp_M[-1]+flank]
	P_New_bp=[bp_P[0]-flank]+bp_P+[bp_P[-1]+flank]
	M_coverage=Block_Assign_To_Letters(bp_MP[0],letter_MP[0],flank)
	P_coverage=Block_Assign_To_Letters(bp_MP[1],letter_MP[1],flank)
	for key in M_coverage.keys():
		M_coverage[key].append(0)
	for key in P_coverage.keys():
		P_coverage[key].append(0)		
	for key in Half_Info.keys():
		Half=Half_Info[key]
		if Half[0]<-flank-100: continue
		else:
			if Half[-1]=='M':
				M_coverage[(Half[0]-(M_New_bp[0]))/100+1][-1]+=1
			elif Half[-1]=='P':
				P_coverage[(Half[0]-(P_New_bp[0]))/100+1][-1]+=1
	return [M_coverage,P_coverage]

def c_GCContent_Calculate_InfoList(Ori_1_Seq,original_bp_list,flank):
	region_length=original_bp_list[-1]-original_bp_list[0]+2*flank
	region_length_new=(region_length/100+1)*100-2*flank
	Number_Of_Blocks=len(Ori_1_Seq)/100
	GC_Content={}
	for i in range(Number_Of_Blocks):
		GC_Content[i+1]=GC_Content_Calculate(Ori_1_Seq[i*100:(i+1)*100])[0]
	return GC_Content

def c_Coverage_Calculate_2a(Letter_Single,Letter_Double,Chromo,original_bp_list,original_letters,flank):
	letter_list=original_letters
	bp_list=[i-original_bp_list[0] for i in original_bp_list]
	bp_list_new=[bp_list[0]-flank]+bp_list+[bp_list[-1]+flank]
	coverage=Block_Assign_To_Letters(bp_list,letter_list,flank)
	for key in coverage.keys():
			coverage[key].append(0)
	for key in Letter_Single.keys():
			for i in Letter_Single[key]:
					keynumL=(i[0]+flank)/100+1
					keynumR=(i[1]+flank)/100+1
					lenL=coverage[keynumL][1]-i[0]
					lenR=i[1]-coverage[keynumR][0]+1
					if lenL>lenR:
							coverage[keynumL][-1]+=1
					else:
							coverage[keynumR][-1]+=1
	for key in Letter_Double.keys():
		for i in Letter_Double[key]:
			keynumL=(i[0]+flank)/100+1
			keynumR=(i[1]+flank)/100+1
			if keynumL in coverage.keys() and keynumR in coverage.keys():
				lenL=coverage[keynumL][1]-i[0]
				lenR=i[1]-coverage[keynumR][0]+1
				if lenL>lenR:
						coverage[keynumL][-1]+=1
				else:
						coverage[keynumR][-1]+=1
			keynumL=(i[2]+flank)/100+1
			keynumR=(i[3]+flank)/100+1
			if keynumL in coverage.keys() and keynumR in coverage.keys():
				lenL=coverage[keynumL][1]-i[0]
				lenR=i[1]-coverage[keynumR][0]+1
				if lenL>lenR:
						coverage[keynumL][-1]+=1
				else:
						coverage[keynumR][-1]+=1
	return coverage

def c_Coverage_Calculate_2b(Letter_Through,Chromo,original_bp_list,original_letters,flank):
	#Eg of RD_Full_Info_of_Reads (a hash list) elements: 'HWI-ST177_136:2:1:7920:85270': [1202, 1302, 1443, 1543, '+', '-']
	letter_list=original_letters
	bp_list=[i-original_bp_list[0] for i in bp_MP[0]]
	bp_list_new=[bp_list[0]-flank]+bp_list+[bp_list[-1]+flank]
	coverage=Block_Assign_To_Letters(bp_list,letter_list,flank)
	for key in coverage.keys():
		coverage[key].append(0)
	for key in Letter_Through.keys():
		i=Letter_Through[key]
		keynumL=(i[0]+flank)/100+1
		keynumR=(i[1]+flank)/100+1
		lenL=coverage[keynumL][1]-i[0]
		lenR=i[1]-coverage[keynumR][0]+1
		if lenL>lenR:
			coverage[keynumL][-1]+=1
		elif lenL<lenR:
			coverage[keynumR][-1]+=1
		elif lenL==lenR:
			coverage[keynumL][-1]+=0.5
			coverage[keynumR][-1]+=0.5
		keynumL=(i[2]+flank)/100+1
		keynumR=(i[3]+flank)/100+1
		lenL=coverage[keynumL][1]-i[0]
		lenR=i[1]-coverage[keynumR][0]+1
		if lenL>lenR:
			coverage[keynumL][-1]+=1
		elif lenL<lenR:
			coverage[keynumR][-1]+=1
		elif lenL==lenR:
			coverage[keynumL][-1]+=0.5
			coverage[keynumR][-1]+=0.5
	return coverage

def c_Coverage_Calculate_2d(Full_Info,Chromo,bp_MP,letter_MP,original_bp_list,flank):
	#Eg of RD_Full_Info_of_Reads (a hash list) elements: 'HWI-ST177_136:2:1:7920:85270': [1202, 1302, 1443, 1543, '+', '-']
	bp_M=[i-original_bp_list[0] for i in bp_MP[0]]
	bp_P=[i-original_bp_list[0] for i in bp_MP[1]]
	M_New_bp=[bp_M[0]-flank]+bp_M+[bp_M[-1]+flank]
	P_New_bp=[bp_P[0]-flank]+bp_P+[bp_P[-1]+flank]
	M_coverage=Block_Assign_To_Letters(bp_MP[0],letter_MP[0],flank)
	P_coverage=Block_Assign_To_Letters(bp_MP[1],letter_MP[1],flank)
	for key in M_coverage.keys():
		M_coverage[key].append(0)
	for key in P_coverage.keys():
		P_coverage[key].append(0)		
	for key in Full_Info.keys():
		if not len(Full_Info[key])==8:
			Halfa=Full_Info[key][:2]+[Full_Info[key][4]]+[Full_Info[key][6]]
			Halfb=Full_Info[key][2:4]+[Full_Info[key][5]]+[Full_Info[key][6]]
			for Half in [Halfa,Halfb]:
				if Half[0]<-flank-100: continue
				else:
					if Half[-1]=='M':
						M_coverage[(Half[0]-(M_New_bp[0]))/100+1][-1]+=1
					elif Half[-1]=='P':
						P_coverage[(Half[0]-(P_New_bp[0]))/100+1][-1]+=1
		elif  len(Full_Info[key])==8:
			Halfa=Full_Info[key][:2]+[Full_Info[key][4]]+[Full_Info[key][6]]
			Halfb=Full_Info[key][2:4]+[Full_Info[key][5]]+[Full_Info[key][6]]
			for Half in [Halfa,Halfb]:
				if Half[0]<-flank-100: continue
				else:
					if Half[-1]=='M':
						M_coverage[(Half[0]-(M_New_bp[0]))/100+1][-1]+=float(Full_Info[key][7])
					elif Half[-1]=='P':
						P_coverage[(Half[0]-(P_New_bp[0]))/100+1][-1]+=float(Full_Info[key][7])
	return [M_coverage,P_coverage]

def c_Coverage_Calculate_2e(Af_Info,Chromo,bp_MP,letter_MP,original_bp_list,flank):
	#Eg of RD_Full_Info_of_Reads (a hash list) elements: 'HWI-ST177_136:2:1:7920:85270': [1202, 1302, 1443, 1543, '+', '-']
	hashM={}
	for i in letter_MP[0]:
		if not i[0] in hashM.keys():
			hashM[i[0]]=[i[0]]
			if (letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))>1:
				hashM[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))[1:]]
	hashP={}
	for i in letter_MP[1]:
		if not i[0] in hashP.keys():
			hashP[i[0]]=[i[0]]
			if (letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))>1:
				hashP[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))[1:]]	
	hashMPLetterBP={}
	hashMPLetterBP['M']={}
	hashMPLetterBP['P']={}
	for j in range(len(letter_MP[0])):
		hashMPLetterBP['M'][hashM[letter_MP[0][j][0]][0]]=[bp_MP[0][j],bp_MP[0][j+1]]
		hashM[letter_MP[0][j][0]].remove(hashM[letter_MP[0][j][0]][0])
	for j in range(len(letter_MP[1])):
		hashMPLetterBP['P'][hashP[letter_MP[1][j][0]][0]]=[bp_MP[1][j],bp_MP[1][j+1]]
		hashP[letter_MP[1][j][0]].remove(hashP[letter_MP[1][j][0]][0])
	hashM={}
	hashM['left']=['left']
	hashM['right']=['right']
	for i in letter_MP[0]:
		if not i[0] in hashM.keys():
			hashM[i[0]]=[i[0]]
			if (letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))>1:
				hashM[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))[1:]]
	hashP={}
	hashP['left']=['left']
	hashP['right']=['right']
	for i in letter_MP[1]:
		if not i[0] in hashP.keys():
			hashP[i[0]]=[i[0]]
			if (letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))>1:
				hashP[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))[1:]]	
	M_Coverage={}
	M_Coverage['left']=0
	for key_1 in hashMPLetterBP['M'].keys():
		M_Coverage[key_1]=[0 for i in range((hashMPLetterBP['M'][key_1][1]-hashMPLetterBP['M'][key_1][0])/100)]
		if ((hashMPLetterBP['M'][key_1][1]-hashMPLetterBP['M'][key_1][0])-(hashMPLetterBP['M'][key_1][1]-hashMPLetterBP['M'][key_1][0])/100*100)>30:
			M_Coverage[key_1].append(0)
	P_Coverage={}
	P_Coverage['left']=0
	for key_1 in hashMPLetterBP['P'].keys():
		P_Coverage[key_1]=[0 for i in range((hashMPLetterBP['P'][key_1][1]-hashMPLetterBP['P'][key_1][0])/100)]
		if ((hashMPLetterBP['P'][key_1][1]-hashMPLetterBP['P'][key_1][0])-(hashMPLetterBP['P'][key_1][1]-hashMPLetterBP['P'][key_1][0])/100*100)>30:
			P_Coverage[key_1].append(0)
	for key in Af_Info.keys():
		if Af_Info[key][0]==Af_Info[key][1]==Af_Info[key][2]==Af_Info[key][3]==(-flank/2):
			M_Coverage['left']+=0.5
			P_Coverage['left']+=0.5
		else:
			if key in Letter_Through.keys():
				if Af_Info[key][6]=='M':
					lele=hashM[Letter_Through[key][6]]
					rile=hashM[Letter_Through[key][9]]
					lebl=Af_Info[key][:2]
					ribl=Af_Info[key][2:4]
					for lele1 in lele:
						if lele1=='left' or lele1=='right': continue
						block=[lele2-bps[0] for lele2 in hashMPLetterBP['M'][lele1]]
						if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
							lebl=[k-block[0] for k in lebl]
							M_Coverage[lele1][lebl[0]/100]+=float(lebl[0]/100*100+100-lebl[0])/float(lebl[1]-lebl[0])
							if lebl[1]/100<len(M_Coverage[lele1]):
								M_Coverage[lele1][lebl[1]/100]+=float(lebl[1]-lebl[1]/100*100)/float(lebl[1]-lebl[0])
					for rile1 in rile:
						if rile1=='left' or rile1=='right':continue
						block=[rile2-bps[0] for rile2 in hashMPLetterBP['M'][rile1]]
						if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
							ribl=[k-block[0] for k in ribl]
							M_Coverage[rile1][ribl[0]/100]+=float(ribl[0]/100*100+100-ribl[0])/float(ribl[1]-ribl[0])
							if ribl[1]/100<len(M_Coverage[rile1]):
								M_Coverage[rile1][ribl[1]/100]+=float(ribl[1]-ribl[1]/100*100)/float(ribl[1]-ribl[0])
				if Af_Info[key][6]=='P':
					lele=hashP[Letter_Through[key][6]]
					rile=hashP[Letter_Through[key][9]]
					lebl=Af_Info[key][:2]
					ribl=Af_Info[key][2:4]
					for lele1 in lele:
						if lele1=='left' or lele1=='right': continue
						block=[lele2-bps[0] for lele2 in hashMPLetterBP['P'][lele1]]
						if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
							lebl=[k-block[0] for k in lebl]
							P_Coverage[lele1][lebl[0]/100]+=float(lebl[0]/100*100+100-lebl[0])/float(lebl[1]-lebl[0])
							if lebl[1]/100<len(P_Coverage[lele1]):
								P_Coverage[lele1][lebl[1]/100]+=float(lebl[1]-lebl[1]/100*100)/float(lebl[1]-lebl[0])
					for rile1 in rile:
						if rile1=='left' or rile1=='right':continue
						block=[rile2-bps[0] for rile2 in hashMPLetterBP['P'][rile1]]
						if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
							ribl=[k-block[0] for k in ribl]
							P_Coverage[rile1][ribl[0]/100]+=float(ribl[0]/100*100+100-ribl[0])/float(ribl[1]-ribl[0])
							if ribl[1]/100<len(P_Coverage[rile1]):
								P_Coverage[rile1][ribl[1]/100]+=float(ribl[1]-ribl[1]/100*100)/float(ribl[1]-ribl[0])
			if not key in Letter_Through.keys():	
					key2='_'.join(key.split('_')[:-1])
					if Af_Info[key][6]=='M':
						lele=hashM[Letter_Through[key2][6]]
						rile=hashM[Letter_Through[key2][9]]
						lebl=Af_Info[key][:2]
						ribl=Af_Info[key][2:4]
						for lele1 in lele:
							if lele1=='left' or lele1=='right': continue
							block=[lele2-bps[0] for lele2 in hashMPLetterBP['M'][lele1]]
							if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
								lebl=[k-block[0] for k in lebl]
								M_Coverage[lele1][lebl[0]/100]+=float(lebl[0]/100*100+100-lebl[0])/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
								if lebl[1]/100<len(M_Coverage[lele1]):
									M_Coverage[lele1][lebl[1]/100]+=float(lebl[1]-lebl[1]/100*100)/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
						for rile1 in rile:
							if rile1=='left' or rile1=='right':continue
							block=[rile2-bps[0] for rile2 in hashMPLetterBP['M'][rile1]]
							if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
								ribl=[k-block[0] for k in ribl]
								M_Coverage[rile1][ribl[0]/100]+=float(ribl[0]/100*100+100-ribl[0])/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
								if ribl[1]/100<len(M_Coverage[rile1]):
									M_Coverage[rile1][ribl[1]/100]+=float(ribl[1]-ribl[1]/100*100)/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
					if Af_Info[key][6]=='P':
						lele=hashP[Letter_Through[key2][6]]
						rile=hashP[Letter_Through[key2][9]]
						lebl=Af_Info[key][:2]
						ribl=Af_Info[key][2:4]
						for lele1 in lele:
							if lele1=='left' or lele1=='right': continue
							block=[lele2-bps[0] for lele2 in hashMPLetterBP['P'][lele1]]
							if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
								lebl=[k-block[0] for k in lebl]
								P_Coverage[lele1][lebl[0]/100]+=float(lebl[0]/100*100+100-lebl[0])/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
								if lebl[1]/100<len(P_Coverage[lele1]):
									P_Coverage[lele1][lebl[1]/100]+=float(lebl[1]-lebl[1]/100*100)/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
						for rile1 in rile:
							if rile1=='left' or rile1=='right':continue
							block=[rile2-bps[0] for rile2 in hashMPLetterBP['P'][rile1]]
							if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
								ribl=[k-block[0] for k in ribl]
								P_Coverage[rile1][ribl[0]/100]+=float(ribl[0]/100*100+100-ribl[0])/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
								if ribl[1]/100<len(P_Coverage[rile1]):
									P_Coverage[rile1][ribl[1]/100]+=float(ribl[1]-ribl[1]/100*100)/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
	return [M_Coverage,P_Coverage]

def block_RD_Calculate_2a(Initial_GCRD_Adj,original_bp_list,flank):
	allele_BP=[0]+[flank+j-original_bp_list[0] for j in original_bp_list]+[2*flank+original_bp_list[-1]-original_bp_list[0]]
	allele_Letter=['left']+[chr(97+i) for i in range(len(original_bp_list)-1)]
	allele_RD=[]
	for k in range(len(allele_Letter)):
		length=allele_BP[k+1]-allele_BP[k]
		block=[allele_BP[k],allele_BP[k+1]]
		temp=[]
		if not block[0]==block[0]/100*100:
			blf=float((block[0]/100+1)*100-block[0])/100*Initial_GCRD_Adj[block[0]/100+1][3]
			temp.append(blf)
			for m in range(block[0]/100+2,block[1]/100+1):
				temp.append(Initial_GCRD_Adj[m][3])
			if not block[1]==block[1]/100*100:
				brf=float(block[1]-block[1]/100*100)/100*Initial_GCRD_Adj[block[1]/100+1][3]
				temp.append(brf)
			allele_RD.append(numpy.sum(temp)/length*100)		
		elif block[0]==block[0]/100*100:
			for m in range(block[0]/100+1,block[1]/100+1):
				temp.append(Initial_GCRD_Adj[m][3])
			if not block[1]==block[1]/100*100:
				brf=float(block[1]-block[1]/100*100)/100*Initial_GCRD_Adj[block[1]/100+1][3]
				temp.append(brf)
			allele_RD.append(numpy.sum(temp)/length*100)		
	return allele_RD

def left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj,flank):
	left_blocks=range(flank/100+1)[1:]
	left_Changes=[Af_GCRD_Adj[j][3]-Through_GCRD_Adj[j][3] for j in left_blocks]
	return numpy.mean(left_Changes)

def All_Block_RD(Initial_block_RD,Af_GCRD_Adj,Af_block_RD,Af_Letter,flank):
	All_Letters=['left']+[chr(97+i) for i in range(len(Initial_block_RD)-1)]
	CNm=[1]+[0 for j in range(len(Initial_block_RD)-1)]
	CNp=[1]+[0 for j in range(len(Initial_block_RD)-1)]
	k=Af_Letter[0]
	for m in k:
		CNm[ord(m[0])-96]+=1
	k=Af_Letter[1]
	for m in k:
		CNp[ord(m[0])-96]+=1
	RDm=[(Initial_block_RD[0]+left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj[0],flank))/2]+[0 for j in range(len(Initial_block_RD)-1)]
	RDp=[(Initial_block_RD[0]+left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj[1],flank))/2]+[0 for j in range(len(Initial_block_RD)-1)]
	RDs=[RDm,RDp]
	for p in range(len(Af_Letter)):
		for q in range(len(Af_Letter[p])):
			RDs[p][ord(Af_Letter[p][q][0])-96]+=Af_block_RD[p][q]
	for r in range(len(Initial_block_RD))[1:]:
		if CNm[r]==CNp[r]:
			RDs[0][r]+=Initial_block_RD[r]/2
			RDs[1][r]+=Initial_block_RD[r]/2
		elif CNm[r]==0 and not CNp[r]==0:
			RDs[1][r]+=Initial_block_RD[r]
		elif CNp[r]==0 and not CNm[r]==0:
			RDs[0][r]+=Initial_block_RD[r]
		else:
			RDs[0][r]+=Initial_block_RD[r]*CNm[r]/(CNp[r]+CNm[r])
			RDs[1][r]+=Initial_block_RD[r]*CNp[r]/(CNp[r]+CNm[r])
	CNs=[CNm,CNp]
	return [CNs,RDs]

def All_Block_RD_2(Initial_block_RD,Af_block_RD,Af_Letter,bps,flank):
	RDs=[[],[]]
	CNs=[[],[]]
	for let in [chr(97+i) for i in range(len(bps)-1)]:
		CNs[0].append(Af_Letter[0].count(let)+Af_Letter[0].count(let+'^'))
		CNs[1].append(Af_Letter[1].count(let)+Af_Letter[1].count(let+'^'))
		if not CNs[0][-1]+CNs[1][-1]==0:
			RDs[0].append(Initial_block_RD[ord(let)-96]*CNs[0][-1]/(CNs[0][-1]+CNs[1][-1]))
			RDs[1].append(Initial_block_RD[ord(let)-96]*CNs[1][-1]/(CNs[0][-1]+CNs[1][-1]))
		if CNs[0][-1]+CNs[1][-1]==0:
			RDs[0].append(0)
			RDs[1].append(0)
	for key in	Af_block_RD[0].keys():
		if not key=='left' and not key=='right':
			RDs[0][ord(key.split('_')[0])-97]+=float(Af_block_RD[0][key])/float(bps[ord(key.split('_')[0])-96]-bps[ord(key.split('_')[0])-97])*100
	for key in	Af_block_RD[1].keys():
		if not key=='left' and not key=='right':
			RDs[1][ord(key.split('_')[0])-97]+=float(Af_block_RD[1][key])/float(bps[ord(key.split('_')[0])-96]-bps[ord(key.split('_')[0])-97])*100
	CNs[0]=[1]+CNs[0]
	CNs[1]=[1]+CNs[1]
	RDs[0]=[Af_block_RD[0]['left']+Initial_block_RD[0]/2]+RDs[0]
	RDs[1]=[Af_block_RD[1]['left']+Initial_block_RD[0]/2]+RDs[1]
	return [CNs,RDs]

def basic_block_RD(Initial_block_RD,original_letters,Af_Letter):
	blockRD={}
	for i in range(len(original_letters)):
		blockRD[original_letters[i]]=Initial_block_RD[i]
	after_letter=[]
	for j in Af_Letter:
		for k in j:
			after_letter.append(k[0])
	afterletter=sorted(after_letter)
	blockCN={}
	for m in after_letter:
		if not m in blockCN.keys():
			blockCN[m]=1
		elif m in blockCN.keys():
			blockCN[m]+=1
	Af_rd=[]
	for n in Af_Letter:
		Af_rd.append([])
		for o in n:
			Af_rd[-1].append(blockRD[o[0]]/blockCN[o[0]])
	return Af_rd

def IL_Stat_Simp(StatFile):
	#for normal distribution used as IL null model 
	fstat=open(StatFile)
	temp=fstat.readline()
	temp=fstat.readline()
	temp=fstat.readline()
	p1=fstat.readline().strip().split()
	Mean1=float(p1[1])
	STD1=float(p1[2])
	fstat.close()
	return [[Mean1,0,STD1,1,1,0],[1,Mean1,STD1]]

def IL_Stat(StatFile):	
	#for bimodal distribution used as IL null model 
	fstat=open(StatFile)
	temp=fstat.readline()
	temp=fstat.readline()
	temp=fstat.readline()
	p1=fstat.readline().strip().split()
	Prop1=float(p1[0])
	Mean1=float(p1[1])
	STD1=float(p1[2])
	fstat.readline()
	p1=fstat.readline().strip().split()
	Prop2=float(p1[0])
	Mean2=float(p1[1])
	STD2=float(p1[2])
	p1=fstat.readline().strip().split()
	p1=fstat.readline().strip().split()
	Prop3=float(p1[0])
	Mean3=float(p1[1])
	STD3=float(p1[2])
	fstat.close()
	return [[Mean1,Mean2,STD1,STD2,Prop1,Prop2],[Prop3,Mean3,STD3]]

def BPs_Coverage(Af_Letter,original_bp_list,original_letters,Letter_Through,Af_Info,flank):
	blocklen={}
	for i in range(len(original_bp_list)-1):
		blocklen[original_letters[i]]=original_bp_list[i+1]-original_bp_list[i]
	blocklen['left']=flank
	blocklen['right']=flank
	tempM=[blocklen[j[0]] for j in Af_Letter[0]]
	tempP=[blocklen[j[0]] for j in Af_Letter[1]]
	Af_BPs=[[-flank,0]+[sum(tempM[:(k+1)]) for k in range(len(tempM))],[-flank,0,]+[sum(tempP[:(k+1)]) for k in range(len(tempP))]]
	Af_BPs=[Af_BPs[0]+[Af_BPs[0][-1]+flank],Af_BPs[1]+[Af_BPs[1][-1]+flank]]
	Af_BP_Through=[[0 for i in range(len(Af_BPs[0]))],[0 for i in range(len(Af_BPs[1]))]]
	for key in Af_Info.keys():
		if Af_Info[key][6]=='M':
			tempbps=Af_BPs[0]
			leftMost=numpy.min([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
			rightMost=numpy.max([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
			for m in range(len(tempbps)-1):
				if tempbps[m+1]>leftMost and tempbps[m]<leftMost:
					for n in range(m,len(tempbps)-1):
						if tempbps[n+1]>rightMost and tempbps[n]<rightMost:
							for p in range(m+1,n+1):
								if len(Af_Info[key])==7:
									Af_BP_Through[0][p]+=1
								elif len(Af_Info[key])==8:
									Af_BP_Through[0][p]+=float(Af_Info[key][7])
		if Af_Info[key][6]=='P':
			tempbps=Af_BPs[1]
			leftMost=numpy.min([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
			rightMost=numpy.max([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
			for m in range(len(tempbps)-1):
				if tempbps[m+1]>leftMost and tempbps[m]<leftMost:
					for n in range(m,len(tempbps)-1):
						if tempbps[n+1]>rightMost and tempbps[n]<rightMost:
							for p in range(m+1,n+1):
								if len(Af_Info[key])==7:
									Af_BP_Through[1][p]+=1
								elif len(Af_Info[key])==8:
									Af_BP_Through[1][p]+=float(Af_Info[key][7])
	return	[Af_BP_Through[0][1:-1],Af_BP_Through[1][1:-1]]

def RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chromo,GC_Hash,RD_List,Let_List):
	Coverage_af_Adj=[]
	Letters=[['left']+Let_List[0]+['right'],['left']+Let_List[1]+['right']]
	Overall_Median_Coverage=float(GC_Overall_Median_Num)
	Coverage_af_Adj=RD_List
	#Theo_RD=GC_Overall_Median_Coverage[str(Chromo)]/2
	#Theo_Var=GC_Var_Coverage[str(Chromo)]/2
	Theo_RD=GC_Overall_Median_Coverage[str(Chromo)]
	Theo_Var=GC_Var_Coverage[str(Chromo)]
	Prob_out=[]
	if Let_List==[[], []]:
		for i in Initial_GCRD_Adj.keys():
			if not i in ['left','right']:
				Prob_out.append(Prob_Norm(Initial_GCRD_Adj[i]*2,0,Theo_Var))
	else:
		for i in Coverage_af_Adj:
			for j in i:
				Prob_out.append(Prob_Norm(j*2,Theo_RD,Theo_Var))
	return numpy.mean(Prob_out)

def delete_block_produce(Letter_List):
	delete_block=[x.upper() for x in Letter_List]
	return delete_block

def delete_BPs_produce(Letter_List):
	delete_BPs=[(j+1,j+2) for j in range(len(Letter_List))]
	return delete_BPs

def insert_block_produce(Letter_List,Letter_List_origin):
	Insert_Pool=[]
	for i in Letter_List_origin:
		if not i in Letter_List and not i+'^' in Letter_List:
			Insert_Pool.append(i)
	return(Insert_Pool)

def insert_BPs_produce(Letter_List,Letter_List_origin):
	Insert_Pool=[]
	for i in Letter_List_origin:
		if not i in Letter_List and not i+'^' in Letter_List:
			Insert_Pool.append(i)	
	Insert_Pool2=[[(ord(j)-96),(ord(j)-95)] for j in Insert_Pool]
	return(Insert_Pool2)

def invert_block_produce(Letter_List):
	invert_block=[x.upper() for x in Letter_List]
	return invert_block

def invert_BPs_produce(Letter_List):
	invert_BPs=[(j+1,j+2) for j in range(len(Letter_List))]
	return invert_BPs

def CopyPaste_block_produce(Letter_List):
	CopyPaste_block=[x.upper() for x in Letter_List]
	return CopyPaste_block

def CopyPaste_BPs_produce(Letter_List):
	CopyPaste_BPs=[tuple(range(len(Letter_List)+2)[1:]) for j in range(len(Letter_List)+1)]
	return CopyPaste_BPs

def CutPaste_block_produce(Letter_List):
	CutPaste_block=[x.upper() for x in Letter_List]
	return CutPaste_block

def CutPaste_BPs_produce(Letter_List):
	BPs=[i for i in range(len(Letter_List)+2)[1:]]
	Letter_List=[chr(96+i) for i in BPs][:-1]
	CutPaste_BPs=[BPs[0:delete_BPs_produce(Letter_List)[ord(j[0])-97][0]-1]+BPs[delete_BPs_produce(Letter_List)[ord(j[0])-97][1]:len(BPs)] for j in Letter_List]
	return CutPaste_BPs

def Move_Choose(Move_Sample_Pool,MoveFlag,Move_probs):
	import random
	random_Pool=range(10**4)
	random_Prob=[0]+Move_probs
	random_Prob=[i*10**4 for i in random_Prob]
	for i in range(len(random_Prob[1:])):
		random_Prob[i+1]+=random_Prob[i]
	Move_numA=random.choice(random_Pool)
	for i in range(len(random_Prob)-1):
		if Move_numA<random_Prob[i+1] and not Move_numA<random_Prob[i]:
			Move_M=Move_Sample_Pool[i]
	if MoveFlag==0:
			Move_P='x'
	elif MoveFlag==1:
			Move_P='x'
	elif MoveFlag==2:
		Move_numA=random.choice(random_Pool)
		for i in range(len(random_Prob)-1):
			if Move_numA<random_Prob[i+1] and not Move_numA<random_Prob[i]:
				Move_P=Move_Sample_Pool[i]
	return [Move_M,Move_P]

def Move_Choice_procedure_2(Move,Letter_List,Letter_List_origin,ChrAllele):
	#ChrAllele Eg: 2m, 2p
	#Only delete, invert, insert
	from random import choice
	#Choice_Pool=['delete','invert','CopyPaste','CutPaste','insert']
	if len(Letter_List)==0 and not Move=='insert':
		return 'ERROR!'
	else:
		if Move=='delete':
			Move_Pool2=[]
			for k in range(len(delete_block_produce(Letter_List))):
				maternal_Move=[ChrAllele,'1',str(delete_BPs_produce(Letter_List)[k][0]),str(delete_BPs_produce(Letter_List)[k][1]),'del']
				Move_Pool2.append(maternal_Move)
		elif Move=='invert':
			Move_Pool2=[]
			for k in range(len(invert_block_produce(Letter_List))):
				maternal_Move=[ChrAllele,'1',str(invert_BPs_produce(Letter_List)[k][0]),str(invert_BPs_produce(Letter_List)[k][1]),'inv']
				Move_Pool2.append(maternal_Move)
		elif Move=='insert':
			Move_Pool2=[]
			insert_p=random.choice(range(len(Letter_List)+2)[1:])
			insert_b=Letter_List_origin
			insert_bl=[[ord(insert_b1)-96,ord(insert_b1)-95] for insert_b1 in insert_b]
			for ip in insert_bl:
				maternal_Move=[ChrAllele,insert_p,ip[0],ip[1],'ins']
				Move_Pool2.append(maternal_Move)
		elif Move=='x':
			Move_Pool2=[]
	return Move_Pool2

def BPList_Invert(BP_List,Command):
	return BP_List

def BPList_Invert_Letter(Letter_List,Command):
	letters_origin=Letter_List
	letters_before=letters_origin[0:(int(Command[2])-1)]
	letters_invert=letters_origin[(int(Command[2])-1):(int(Command[3])-1)]
	letters_after=letters_origin[(int(Command[3])-1):]
	letters_invert_after=[]
	for j in range(len(letters_invert)):
		letters_invert_after.append(letters_invert[len(letters_invert)-1-j])
	for k in range(len(letters_invert_after)):
		if not letters_invert_after[k][-1]=='^':
			letters_invert_after[k]+='^'
		elif letters_invert_after[k][-1]=='^':
			letters_invert_after[k]=letters_invert_after[k][:-1]
	return letters_before+letters_invert_after+letters_after		

def BPList_Delete(BP_List,Command):
	bp1=BP_List[int(Command[2])-1]
	bp2=BP_List[int(Command[3])-1]
	delete=BP_List[(int(Command[2])-1):int(Command[3])]
	delete_length=delete[-1]-delete[0]
	delete_after=BP_List[int(Command[3]):]
	for i in range(len(delete_after)):
		delete_after[i]=delete_after[i]-delete_length
	return BP_List[0:(int(Command[2])-1)]+[delete[0]]+delete_after	

def BPList_Delete_Letter(Letter_List,Command):
	letters_origin=Letter_List
	letters_delete=letters_origin[(int(Command[2])-1):(int(Command[3])-1)]
	return letters_origin[0:(int(Command[2])-1)]+letters_origin[(int(Command[3]))-1:]		

def BPList_Insert_Pool_Sequence(BP_List_origin,Ref_Sequence_Path_File,flank):
#Eg:Ref_Sequence_Path_File='/home/xuefzhao/References/chr2:184569179-184572016.fa'
#This function produce a hash list that contains the original sequences
	fi=open(Ref_Sequence_Path_File)
	temp=fi.readline().strip().split()
	seq1=[]
	while True:
		seq_temp=fi.readline().rstrip()
		if not seq_temp: break
		else:
			seq1.append(seq_temp)
	seq_ori=''.join(seq1)
	relative_List=map(lambda a:a-BP_List_origin[0]+flank,BP_List_origin)
	Letter_List_origin=[chr(x+97) for x in range(len(BP_List_origin)-1)]
	seq_letter={}
	for i in range(len(BP_List_origin)-1):
		seq_letter[Letter_List_origin[i]]=seq_ori[relative_List[i]:relative_List[i+1]]
		seq_letter[Letter_List_origin[i]+'^']=complementary(seq_ori[relative_List[i]:relative_List[i+1]])
	fi.close()
	return seq_letter

def BPList_Insert_Pool_Length(BP_List_origin):
#Eg:Ref_Sequence_Path_File='/home/xuefzhao/References/chr2:184569179-184572016.fa'
#This function produce a hash list that contains the original sequences
	Letter_List_origin=[chr(x+97) for x in range(len(BP_List_origin)-1)]
	seq_letter={}
	for i in range(len(BP_List_origin)-1):
		seq_letter[Letter_List_origin[i]]=int(BP_List_origin[i+1])-int(BP_List_origin[i])
		seq_letter[Letter_List_origin[i]+'^']=int(BP_List_origin[i+1])-int(BP_List_origin[i])
	return seq_letter

def BPList_Insert(BP_List,Command,BP_List_origin):
	Insert_block=chr(int(Command[2])+96)
	IL=BPList_Insert_Pool_Length(BP_List_origin)[Insert_block]
	Insert_Before=BP_List[0:int(Command[1])]
	insert_diff=[x+IL for x in BP_List[(int(Command[1])-1):]]
	return Insert_Before+insert_diff

def BPList_Insert_Letter(Letter_List,Command):
	Insert_block=chr(int(Command[2])+96)
	Letters_Before=Letter_List[0:(int(Command[1])-1)]
	Letters_After=Letter_List[(int(Command[1])-1):]
	return Letters_Before+[Insert_block]+Letters_After

def BPList_CopyPaste(BP_List,Command):
	if int(Command[1])<=int(Command[2]):
		Block_Before=BP_List[0:int(Command[1])]
		Block_ip_bp1=BP_List[(int(Command[1])-1):int(Command[2])]
		Block_bp1_bp2=BP_List[(int(Command[2])-1):int(Command[3])]
		Block_After=BP_List[(int(Command[3])-1):]
 		Block_bp1_bp2_cp=[]
		Dist1=Block_bp1_bp2[0]-Block_Before[-1]
		for i in range(len(Block_bp1_bp2)):
			Block_bp1_bp2_cp.append(Block_bp1_bp2[i]-Dist1)
		Dist2=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
		for j in range(len(Block_ip_bp1)):
			Block_ip_bp1[j]=Block_ip_bp1[j]+Dist2
		for j in range(len(Block_bp1_bp2)):
			Block_bp1_bp2[j]=Block_bp1_bp2[j]+Dist2
		for j in range(len(Block_After)):
			Block_After[j]=Block_After[j]+Dist2
		return Block_Before[:-1]+Block_bp1_bp2_cp[:-1]+Block_ip_bp1[:-1]+Block_bp1_bp2[:-1]+Block_After
	if int(Command[1])>=int(Command[3]):
		Block_Before=BP_List[0:int(Command[2])]
		Block_bp1_bp2=BP_List[(int(Command[2])-1):int(Command[3])]
		Block_bp2_ip=BP_List[(int(Command[3])-1):int(Command[1])]
		Block_After=BP_List[(int(Command[1])-1):]
 		Block_bp1_bp2_cp=[]
		Dist1=Block_After[0]-Block_bp1_bp2[0]
		for i in range(len(Block_bp1_bp2)):
			Block_bp1_bp2_cp.append(Block_bp1_bp2[i]+Dist1)
		Dist2=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
		for j in range(len(Block_After)):
			Block_After[j]=Block_After[j]+Dist2
		return Block_Before[:-1]+Block_bp1_bp2[:-1]+Block_bp2_ip[:-1]+Block_bp1_bp2_cp[:-1]+Block_After	

def BPList_CopyPaste_Letter(Letter_List,Command):
	letters_origin=Letter_List
	letters_before=letters_origin[0:(int(Command[1])-1)]
	letters_copy=letters_origin[(int(Command[2])-1):(int(Command[3])-1)]
	letters_after=letters_origin[(int(Command[1])-1):]
	return letters_before+letters_copy+letters_after	

def BPList_CutPaste(BP_List,Command):
	BP_Line=Command
	if int(BP_Line[1])<int(BP_Line[2]):
		Block_Before=BP_List[0:int(BP_Line[1])]
		Block_ip_bp1=BP_List[(int(BP_Line[1])-1):int(BP_Line[2])]
		Block_bp1_bp2=BP_List[(int(BP_Line[2])-1):int(BP_Line[3])]
		Block_After=BP_List[(int(BP_Line[3])-1):]
		Dist1=Block_bp1_bp2[0]-Block_Before[-1]
		for i in range(len(Block_bp1_bp2)):
			Block_bp1_bp2[i]=Block_bp1_bp2[i]-Dist1
		Dist2=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
		for j in range(len(Block_ip_bp1)):
			Block_ip_bp1[j]=Block_ip_bp1[j]+Dist2
		return Block_Before[:-1]+Block_bp1_bp2[:-1]+Block_ip_bp1[:-1]+Block_After
	if int(BP_Line[1])>int(BP_Line[3]):
		Block_Before=BP_List[0:int(BP_Line[2])]
		Block_bp1_bp2=BP_List[(int(BP_Line[2])-1):int(BP_Line[3])]
		Block_bp2_ip=BP_List[(int(BP_Line[3])-1):int(BP_Line[1])]
		Block_After=BP_List[(int(BP_Line[1])-1):]
		Dist1=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
		for i in range(len(Block_bp2_ip)):
			Block_bp2_ip[i]=Block_bp2_ip[i]-Dist1
		Dist2=Block_After[0]-Block_bp1_bp2[-1]
		for j in range(len(Block_bp1_bp2)):
			Block_bp1_bp2[j]=Block_bp1_bp2[j]+Dist2
		return Block_Before[:-1]+Block_bp2_ip[:-1]+Block_bp1_bp2[:-1]+Block_After		

def BPList_CutPaste_Letter(Letter_List,Command):
	Letter_Line=Command
	if int(Letter_Line[1])<int(Letter_Line[2]):	
		Block_Before=Letter_List[:(numpy.min(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1)]
		Block_After=Letter_List[(numpy.max(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1):]		
		Block_cut=Letter_List[(int(Letter_Line[2])-1):(int(Letter_Line[3])-1)]
		Block_ip_cut=Letter_List[(int(Letter_Line[1])-1):(int(Letter_Line[2])-1)]
		return Block_Before+Block_cut+Block_ip_cut+Block_After
	if  int(Letter_Line[1])>int(Letter_Line[3]):
		Block_Before=Letter_List[:(numpy.min(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1)]
		Block_After=Letter_List[(numpy.max(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1):]		
		Block_cut=Letter_List[(int(Letter_Line[2])-1):(int(Letter_Line[3])-1)]
		Block_ip_cut=Letter_List[(int(Letter_Line[3])-1):(int(Letter_Line[1])-1)]
		return Block_Before+Block_ip_cut+Block_cut+Block_After	

def BPList_X(BP_List,Command):
	return BP_List

def BPList_X_Letter(Letter_List,Command):
	return Letter_List

def BPList_Rearrange(BP_List,Command,BP_List_origin):
	if Command[-1]=='del' or Command[-1]=='delete':
		return BPList_Delete(BP_List,Command)
	elif Command[-1]=='inv' or Command[-1]=='invert':
		return BPList_Invert(BP_List,Command)	
	elif Command[-1]=='ins' or Command[-1]=='insert':
		return BPList_Insert(BP_List,Command,BP_List_origin)
	elif Command[-1]=='copy+paste' or Command[-1]=='CopyPaste':
		return BPList_CopyPaste(BP_List,Command)
	elif Command[-1]=='cut+paste' or Command[-1]=='CutPaste':
		return BPList_CutPaste(BP_List,Command)
	elif Command[-1]=='x' or Command[-1]=='X':
		return BPList_X(BP_List,Command)

def LetterList_Rearrange(Letter_List,Command,BP_List_origin):
	if Command[-1]=='del' or Command[-1]=='delete':
		return BPList_Delete_Letter(Letter_List,Command)
	elif Command[-1]=='inv' or Command[-1]=='invert':
		return BPList_Invert_Letter(Letter_List,Command)	
	elif Command[-1]=='ins' or Command[-1]=='insert':
		return BPList_Insert_Letter(Letter_List,Command)
	elif Command[-1]=='copy+paste' or Command[-1]=='CopyPaste':
		return BPList_CopyPaste_Letter(Letter_List,Command)
	elif Command[-1]=='cut+paste' or Command[-1]=='CutPaste':
		return BPList_CutPaste_Letter(Letter_List,Command)
	elif Command[-1]=='x' or Command[-1]=='X':
		return BPList_X_Letter(Letter_List,Command)

def Rearrange_RefSeq(Ref_Sequence_Path_File,BP_List_origin,Letter_List_After,flank):
	Seq_pool=BPList_Insert_Pool_Sequence(BP_List_origin,Ref_Sequence_Path_File,flank)
	Seq_Out=[Seq_pool[x] for x in Letter_List_After]
	return(''.join(Seq_Out))

def Sequence_Rearrange_3(Insert_Seq,Af_Letter,original_bp_list,flank):
	#Eg of Af_Letter:['a','b','c']
	seql=[]
	seql.append(Insert_Seq['left'])
	for let in Af_Letter:
		seql.append(Insert_Seq['I'+let])
	seql.append(Insert_Seq['right'])
	return ''.join(seql)

def MP_Letter_Through_Read_Rearrange_3(Full_Info_of_Read,Af_Letter,original_bp):
	Full_Info_of_Read=[int(i) for i in Full_Info_of_Read[:4]]+Full_Info_of_Read[4:]
	Left_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
	Right_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
	if numpy.max(Left_Read[:2])>numpy.max(Right_Read[:2]):
		Left_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
		Right_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
	Original_Letter=[chr(i+97) for i in range(len(original_bp)-1)]
	Letter_Len={}
	for i in range(len(Original_Letter)):
		Letter_Len[Original_Letter[i]]=original_bp[i+1]-original_bp[i]
	bp_MP=[[0],[0]]
	for j in range(len(Af_Letter)):
		for k in Af_Letter[j]:
			if k=='le' or k=='ri':
				bp_MP[j].append(bp_MP[j][-1]+Letter_Len[k])
			else:
				bp_MP[j].append(bp_MP[j][-1]+Letter_Len[k[0]])
	Af_M_let=['left']+Af_Letter[0]+['right']
	Af_M_bpp=[-flank]+bp_MP[0]+[bp_MP[0][-1]+flank]
	Af_P_let=['left']+Af_Letter[1]+['right']
	Af_P_bpp=[-flank]+bp_MP[1]+[bp_MP[1][-1]+flank]
	LeftM_New=[]
	RightM_New=[]
	for i in range(len(Af_M_let)):
		if Left_Read[3][0]==Af_M_let[i]:
			LeftM_New.append([Af_M_bpp[i]+Left_Read[3][1],Af_M_bpp[i]+Left_Read[3][2],Left_Read[2]])
		elif Left_Read[3][0]+'^'==Af_M_let[i]:
			LeftM_New.append([Af_M_bpp[i+1]-Left_Read[3][2],Af_M_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
		elif Left_Read[3][0]==Af_M_let[i]+'^':
			LeftM_New.append([Af_M_bpp[i+1]-Left_Read[3][2],Af_M_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
		if Right_Read[3][0]==Af_M_let[i]:
			RightM_New.append([Af_M_bpp[i]+Right_Read[3][1],Af_M_bpp[i]+Right_Read[3][2],Right_Read[2]])
		elif Right_Read[3][0]+'^'==Af_M_let[i]:
			RightM_New.append([Af_M_bpp[i+1]-Right_Read[3][2],Af_M_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
		elif Right_Read[3][0]==Af_M_let[i]+'^':
			RightM_New.append([Af_M_bpp[i+1]-Right_Read[3][2],Af_M_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
	LeftP_New=[]
	RightP_New=[]
	for i in range(len(Af_P_let)):
		if Left_Read[3][0]==Af_P_let[i]:
			LeftP_New.append([Af_P_bpp[i]+Left_Read[3][1],Af_P_bpp[i]+Left_Read[3][2],Left_Read[2]])
		elif Left_Read[3][0]+'^'==Af_P_let[i]:
			LeftP_New.append([Af_P_bpp[i+1]-Left_Read[3][2],Af_P_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
		elif Left_Read[3][0]==Af_P_let[i]+'^':
			LeftP_New.append([Af_P_bpp[i+1]-Left_Read[3][2],Af_P_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
		if Right_Read[3][0]==Af_P_let[i]:
			RightP_New.append([Af_P_bpp[i]+Right_Read[3][1],Af_P_bpp[i]+Right_Read[3][2],Right_Read[2]])
		elif Right_Read[3][0]+'^'==Af_P_let[i]:
			RightP_New.append([Af_P_bpp[i+1]-Right_Read[3][2],Af_P_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
		elif Right_Read[3][0]==Af_P_let[i]+'^':
			RightP_New.append([Af_P_bpp[i+1]-Right_Read[3][2],Af_P_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
	IN_Read_New=[]
	OUT_Read_New=[]
	for l in LeftM_New:
		for r in RightM_New:
			if not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])<Cut_Lower and not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])>Cut_Upper and not l[2]==r[2]:
				IN_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['M'])
			else:
				OUT_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['M'])
	for l in LeftP_New:
		for r in RightP_New:
			if not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])<Cut_Lower and not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])>Cut_Upper and not l[2]==r[2]:
				IN_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['P'])
			else:
				OUT_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['P'])
	if not len(IN_Read_New)==0:
		return ['IN']+IN_Read_New
	elif len(IN_Read_New)==0 and not len(OUT_Read_New)==0:
		return ['OUT']+OUT_Read_New
	else:
		return ['NON',[-flank/2,-flank/2,-flank/2,-flank/2,'+','-','M'],[-flank/2,-flank/2,-flank/2,-flank/2,'+','-','P']]

def MP_Letter_Through_Read_Rearrange_4(Full_Info_of_Read,Letter_Len,Af_BP_rela,Af_Letter_rela,original_bp):
	Full_Info_of_Read=[int(i) for i in Full_Info_of_Read[:4]]+Full_Info_of_Read[4:]
	Left_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
	Right_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
	if numpy.max(Left_Read[:2])>numpy.max(Right_Read[:2]):
		Left_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
		Right_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
	Left_Reads=[]
	Right_Reads=[]
	Left_Info=Left_Read[3]
	Right_Info=Right_Read[3]
	if Left_Read[3][0]=='left':
		Left_Reads.append(Left_Read[:3]+['M',Af_Letter_rela[0][0]])
		Left_Reads.append(Left_Read[:3]+['P',Af_Letter_rela[1][0]])
	if Left_Read[3][0]=='right':
		Left_Reads.append([Af_BP_rela[0][-2]+j for j in Left_Read[3][1:]]+[Left_Read[2]]+['M',Af_Letter_rela[0][-1]])
		Left_Reads.append([Af_BP_rela[1][-2]+j for j in Left_Read[3][1:]]+[Left_Read[2]]+['P',Af_Letter_rela[1][-1]])
	if not Left_Read[3][0] in ['left','right']:
		for i in range(len(Af_Letter_rela)):
			for j in range(len(Af_Letter_rela[i]))[1:-1]:
				if Af_Letter_rela[i][j].split('_')[0][0]==Left_Read[3][0]:
					if not Af_Letter_rela[i][j].split('_')[0][-1]=='^':
						Left_Reads.append([Af_BP_rela[i][j]+k for k in Left_Read[3][1:]]+[Left_Read[2],['M','P'][i],Af_Letter_rela[i][j]])
					if Af_Letter_rela[i][j].split('_')[0][-1]=='^':
						Left_Reads.append(sorted([Af_BP_rela[i][j+1]-k for k in Left_Read[3][1:]])+[oppo_direct(Left_Read[2]),['M','P'][i],Af_Letter_rela[i][j]])
	if Right_Read[3][0]=='left':
		Right_Reads.append(Right_Read[:3]+['M',Af_Letter_rela[0][0]])
		Right_Reads.append(Right_Read[:3]+['P',Af_Letter_rela[1][0]])
	if Right_Read[3][0]=='right':
		Right_Reads.append([Af_BP_rela[0][-2]+j for j in Right_Read[3][1:]]+[Right_Read[2]]+['M',Af_Letter_rela[0][-1]])
		Right_Reads.append([Af_BP_rela[1][-2]+j for j in Right_Read[3][1:]]+[Right_Read[2]]+['P',Af_Letter_rela[1][-1]])
	if not Right_Read[3][0] in ['left','right']:
		for i in range(len(Af_Letter_rela)):
			for j in range(len(Af_Letter_rela[i]))[1:-1]:
				if Af_Letter_rela[i][j].split('_')[0][0]==Right_Read[3][0]:
					if not Af_Letter_rela[i][j].split('_')[0][-1]=='^':
						Right_Reads.append([Af_BP_rela[i][j]+k for k in Right_Read[3][1:]]+[Right_Read[2],['M','P'][i],Af_Letter_rela[i][j]])
					if Af_Letter_rela[i][j].split('_')[0][-1]=='^':
						Right_Reads.append(sorted([Af_BP_rela[i][j+1]-k for k in Right_Read[3][1:]])+[oppo_direct(Right_Read[2]),['M','P'][i],Af_Letter_rela[i][j]])
	In_Pair_Reads=['IN']
	Out_Pair_Reads=['OUT']
	for m in Left_Reads:
		for n in Right_Reads:
			if m[3]==n[3]:
				Pair_Read=m[:2]+n[:2]+[m[2],n[2],m[3],m[4],n[4]]
				if (numpy.max(Pair_Read[:4])-numpy.min(Pair_Read[:4]))>Cut_Lower and (numpy.max(Pair_Read[:4])-numpy.min(Pair_Read[:4]))<Cut_Upper and not Pair_Read[4]==Pair_Read[5]:
					In_Pair_Reads.append(Pair_Read)	
				else:	
					Out_Pair_Reads.append(Pair_Read)	
	output=[]
	if not len(In_Pair_Reads)==1:
		for i in In_Pair_Reads[1:]:
			output.append(i+[Left_Info[0],Left_Info[1]/100,float(Left_Info[1]/100*100+100-Left_Info[1])/float(Left_Info[2]-Left_Info[1]),Left_Info[2]/100,float(Left_Info[2]-Left_Info[2]/100*100)/float(Left_Info[2]-Left_Info[1])]+[Right_Info[0],Right_Info[1]/100,float(Right_Info[1]/100*100+100-Right_Info[1])/float(Right_Info[2]-Right_Info[1]),Right_Info[2]/100,float(Right_Info[2]-Right_Info[2]/100*100)/float(Right_Info[2]-Right_Info[1])])
	else:
		if not len(Out_Pair_Reads)==1:
			output1=[j for j in Out_Pair_Reads[1:] if not j[4]==j[5]]
			output2=[j for j in Out_Pair_Reads[1:] if j[4]==j[5]]
			if not len(output1)==0:
				pdfs=[pdf_calculate(numpy.max(output1[i][:4])-numpy.min(output1[i][:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for i in range(len(output1))]
				for j in range(len(pdfs)):
					if pdfs[j]==numpy.max(pdfs):
						output.append(output1[j]+[Left_Info[0],Left_Info[1]/100,float(Left_Info[1]/100*100+100-Left_Info[1])/float(Left_Info[2]-Left_Info[1]),Left_Info[2]/100,float(Left_Info[2]-Left_Info[2]/100*100)/float(Left_Info[2]-Left_Info[1])]+[Right_Info[0],Right_Info[1]/100,float(Right_Info[1]/100*100+100-Right_Info[1])/float(Right_Info[2]-Right_Info[1]),Right_Info[2]/100,float(Right_Info[2]-Right_Info[2]/100*100)/float(Right_Info[2]-Right_Info[1])])
			else:
				pdfs=[pdf_calculate(numpy.max(output2[i][:4])-numpy.min(output2[i][:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for i in range(len(output2))]
				for j in range(len(pdfs)):
					if pdfs[j]==numpy.max(pdfs):
						output.append(output2[j]+[Left_Info[0],Left_Info[1]/100,float(Left_Info[1]/100*100+100-Left_Info[1])/float(Left_Info[2]-Left_Info[1]),Left_Info[2]/100,float(Left_Info[2]-Left_Info[2]/100*100)/float(Left_Info[2]-Left_Info[1])]+[Right_Info[0],Right_Info[1]/100,float(Right_Info[1]/100*100+100-Right_Info[1])/float(Right_Info[2]-Right_Info[1]),Right_Info[2]/100,float(Right_Info[2]-Right_Info[2]/100*100)/float(Right_Info[2]-Right_Info[1])])
		else:
			output.append([-flank/2,-flank/2,-flank/2,-flank/2,'+','-','M','left_0','left_0','left',1,1,0,0,'left',1,1,0,0])
			output.append([-flank/2,-flank/2,-flank/2,-flank/2,'+','-','P','left_0','left_0','left',1,1,0,0,'left',1,1,0,0])
	return 	output

def oppo_direct(x):
	if x=='+':
		return '-'
	elif x=='-':
		return '+'

def Letter_Through_Rearrange_3(IL_Statistics,Letter_Through,Af_Letter,original_bp_list):
	#Through reads that got mapped to several spots are divided absolutely equally to each spot
	New_Info_of_Reads={}
	for key in Letter_Through.keys():
		after_key=MP_Letter_Through_Read_Rearrange_3(Letter_Through[key],Af_Letter,original_bp_list)
		if len(after_key)==2:
			New_Info_of_Reads[key]=after_key[1]
		elif len(after_key)>2:
			if after_key[0]=='IN' or after_key[0]=='NON':
				for i in range(len(after_key))[1:]:
					New_Info_of_Reads[key+'_'+str(i)]=after_key[i]+[str(float(1)/float(len(after_key)-1))]
			elif after_key[0]=='OUT':
				pdfs=[pdf_calculate(numpy.max(after_key[i][:4])-numpy.min(after_key[i][:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for i in range(len(after_key))[1:]]
				best_Out=[after_key[j+1] for j in range(len(pdfs)) if pdfs[j]==numpy.max(pdfs)]
				if len(best_Out)==1:
					New_Info_of_Reads[key]=best_Out[0]
				elif len(best_Out)>1:
					for i in range(len(best_Out)+1)[1:]:
						New_Info_of_Reads[key+'_'+str(i)]=best_Out[i-1]+[str(float(1)/float(len(best_Out)))]
	return New_Info_of_Reads

def complement(direction):
	if direction=='+':
		return '-'
	elif direction=='-':
		return '+'
	else:
		return 'error'

min_resolution=70
def candidate_QC_Control(Read_List):
	if Read_List==[]:
		return []
	else:
		Qual_Filter_1=[]
		for j in Read_List:
			if not j[1]-j[0]>ReadLength+min_resolution and j[1]-j[0]>0 and not j[3]-j[2]>ReadLength+min_resolution and j[3]-j[2]>0:
				Qual_Filter_1.append(j)
		if not Qual_Filter_1==[]:
			if len(Qual_Filter_1)==1:
				Qual_Filter_1[0]+=[pdf_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_1]
				return Qual_Filter_1
			else:
				Qual_Filter_2=[]
				for j2 in Qual_Filter_1:
					if j2[-2:]==['+','-']:
						  Qual_Filter_2.append(j2)
				if not Qual_Filter_2==[]:
					if len(Qual_Filter_2)==1:
						Qual_Filter_2[0]+=[pdf_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
						return Qual_Filter_2
					else:
						Qual_Filter_3=[]
						Qual_IL=[pdf_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
						for jq in range(len(Qual_IL)):
							if Qual_IL[jq]==max(Qual_IL) and not Qual_Filter_1[jq] in Qual_Filter_3:
								Qual_Filter_3.append(Qual_Filter_1[jq]+[max(Qual_IL)])
						return Qual_Filter_3					
				else:
					Qual_Filter_2=Qual_Filter_1
					if len(Qual_Filter_2)==1:
						Qual_Filter_2[0]+=[pdf_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
						return Qual_Filter_2
					else:
						Qual_Filter_3=[]
						Qual_IL=[pdf_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
						for jq in range(len(Qual_IL)):
							if Qual_IL[jq]==max(Qual_IL) and not Qual_Filter_1[jq] in Qual_Filter_3:
								Qual_Filter_3.append(Qual_Filter_1[jq]+[max(Qual_IL)])
						return Qual_Filter_3 
		else:
			return []

def candidate_QC_Control2(M_Read_List,P_Read_List):
	Qual_Filter_1=[]
	for i in M_Read_List:
		Qual_Filter_1.append(i+['m'])
	for i in P_Read_List:
		Qual_Filter_1.append(i+['p'])
	Qual_Filter_2=[]
	for i in Qual_Filter_1:
		if i[-4:-2]==['+','-']:
			Qual_Filter_2.append(i)
	if not Qual_Filter_2==[]:
		Qual_Filter_3=[]
		IL_Qual=[pdf_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
		for j in range(len(IL_Qual)):
			if IL_Qual[j]==max(IL_Qual) and not Qual_Filter_2[j] in Qual_Filter_3:
				Qual_Filter_3.append(Qual_Filter_2[j])
	else:
		Qual_Filter_2=Qual_Filter_1
		Qual_Filter_3=[]
		IL_Qual=[pdf_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
		for j in range(len(IL_Qual)):
			if IL_Qual[j]==max(IL_Qual) and not Qual_Filter_2[j] in Qual_Filter_3:
				Qual_Filter_3.append(Qual_Filter_2[j])
	return Qual_Filter_3

def Cov_Cal_Block(pos,bp,cov,perc):
	 for j in range(len(bp)-2):
		 if not pos[0]<bp[j] and pos[0]<bp[j+1]:
			if not pos[1]<bp[j] and pos[1]<bp[j+1]:
				cov[j]+=(pos[1]-pos[0])*perc
			elif not pos[1]<bp[j+1] and pos[1]<bp[j+2]:
				cov[j]+=(bp[j+1]-pos[0])*perc
				cov[j+1]+=(pos[1]-bp[j+1])*perc
			elif not pos[1]<temp_bp[0][j+2] and pos[1]<temp_bp[0][j+3]:
				cov[j]+=(bp[j+1]-pos[0])*perc
				cov[j+1]+=(bp[j+2]-bp[j+1])*perc
				cov[j+2]+=(pos[1]-bp[j+2])*perc
	 j=len(bp)-2
	 if not pos[0]<bp[j] and pos[0]<bp[j+1]:
		if not pos[1]<bp[j] and pos[1]<bp[j+1]:
			cov[j]+=(pos[1]-pos[0])*perc
		else:
			cov[j]+=(bp[j+1]-pos[0])*perc

def penal_calculate(Map_All,temp_bp, Af_Letter,Af_BP,RD_within_B,letters_numbers,BlockGC,NoMapPenal):
	out_rd=[[0 for i in temp_bp[0][:-1]],[0 for i in temp_bp[1][:-1]]]
	IL_Rec={}
	DR_Penal=0
	out_tb=[[0 for i in temp_bp[0]],[0 for i in temp_bp[1]]]
	for i in Map_All:
		if len(i)>4:
			if not i[6] in IL_Rec.keys():
				IL_Rec[i[6]]=i[8]
			else:
				IL_Rec[i[6]]+=i[8]
			if not i[4:6]==['+','-']:
				DR_Penal+=1
			if i[7]=='m':
				i_block=[]
				for k in i[:4]:
					if k<temp_bp[0][1]:
						i_block.append(0)
					elif k>temp_bp[0][-2]-1:
						i_block.append(len(temp_bp[0])-2)
					else:
						for j in range(len(temp_bp[0])-1)[1:-1]:
							if temp_bp[0][j]-1<k and temp_bp[0][j+1]>k:
								i_block.append(j)
				if i_block[0]==i_block[1] and i_block[2]==i_block[3]:
					out_rd[0][i_block[0]]+=(i[1]-i[0])*i[-1]
					out_rd[0][i_block[2]]+=(i[3]-i[2])*i[-1]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[0][k2]+=i[8]/(i_block[2]-i_block[1])
				elif not i_block[0]==i_block[1] and i_block[2]==i_block[3]:
					out_rd[0][i_block[0]]+=(temp_bp[0][i_block[0]+1]-i[0])*i[-1]
					out_rd[0][i_block[1]]+=(i[1]-temp_bp[0][i_block[1]])*i[-1]
					out_rd[0][i_block[2]]+=(i[3]-i[2])*i[-1]
					out_tb[0][i_block[1]]+=i[8]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[0][k2]+=i[8]
				elif i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
					out_rd[0][i_block[0]]+=(i[1]-i[0])*i[-1]
					out_rd[0][i_block[2]]+=(temp_bp[0][i_block[2]+1]-i[2])*i[-1]
					out_rd[0][i_block[3]]+=(i[3]-temp_bp[0][i_block[3]])*i[-1]
					out_tb[0][i_block[3]]+=i[8]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[0][k2]+=i[8]
				elif not i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
					out_rd[0][i_block[0]]+=(temp_bp[0][i_block[0]+1]-i[0])*i[-1]
					out_rd[0][i_block[1]]+=(i[1]-temp_bp[0][i_block[1]])*i[-1]
					out_rd[0][i_block[2]]+=(temp_bp[0][i_block[2]+1]-i[2])*i[-1]
					out_rd[0][i_block[3]]+=(i[3]-temp_bp[0][i_block[3]])*i[-1]
					out_tb[0][i_block[1]]+=i[8]
					out_tb[0][i_block[3]]+=i[8]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[0][k2]+=i[8]
			if i[7]=='p':
				i_block=[]
				for k in i[:4]:
					if k<temp_bp[1][1]:
						i_block.append(0)
					elif k>temp_bp[1][-2]-1:
						i_block.append(len(temp_bp[1])-2)
					else:
						for j in range(len(temp_bp[1])-1)[1:-1]:
							if temp_bp[1][j]-1<k and temp_bp[1][j+1]>k:
								i_block.append(j)
				if i_block[0]==i_block[1] and i_block[2]==i_block[3]:
					out_rd[1][i_block[0]]+=(i[1]-i[0])*i[-1]
					out_rd[1][i_block[2]]+=(i[3]-i[2])*i[-1]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[1][k2]+=i[8]
				elif not i_block[0]==i_block[1] and i_block[2]==i_block[3]:
					out_rd[1][i_block[0]]+=(temp_bp[1][i_block[0]+1]-i[0])*i[-1]
					out_rd[1][i_block[1]]+=(i[1]-temp_bp[1][i_block[1]])*i[-1]
					out_rd[1][i_block[2]]+=(i[3]-i[2])*i[-1]
					out_tb[1][i_block[1]]+=i[8]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[1][k2]+=i[8]
				elif i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
					out_rd[1][i_block[0]]+=(i[1]-i[0])*i[-1]
					out_rd[1][i_block[2]]+=(temp_bp[1][i_block[2]+1]-i[2])*i[-1]
					out_rd[1][i_block[3]]+=(i[3]-temp_bp[1][i_block[3]])*i[-1]
					out_tb[1][i_block[3]]+=i[8]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[1][k2]+=i[8]
				elif not i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
					out_rd[1][i_block[0]]+=(temp_bp[1][i_block[0]+1]-i[0])*i[-1]
					out_rd[1][i_block[1]]+=(i[1]-temp_bp[1][i_block[1]])*i[-1]
					out_rd[1][i_block[2]]+=(temp_bp[1][i_block[2]+1]-i[2])*i[-1]
					out_rd[1][i_block[3]]+=(i[3]-temp_bp[1][i_block[3]])*i[-1]
					out_tb[1][i_block[1]]+=i[8]
					out_tb[1][i_block[3]]+=i[8]
					if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
						for k2 in range(i_block[1]+1,i_block[2]+1):
							out_tb[1][k2]+=i[8]
		else:
			if i[2]=='m':
				i_block=[]
				for k in i[:2]:
					if k<temp_bp[0][1]:
						i_block.append(0)
					elif k>temp_bp[0][-2]-1:
						i_block.append(len(temp_bp[0])-2)
					else:
						for j in range(len(temp_bp[0])-1)[1:-1]:
							if temp_bp[0][j]-1<k and temp_bp[0][j+1]>k:
								i_block.append(j)
				if i_block[0]==i_block[1]:
					out_rd[0][i_block[0]]+=(i[1]-i[0])*i[-1]
				elif not i_block[0]==i_block[1]:
					out_rd[0][i_block[0]]+=(temp_bp[0][i_block[0]+1]-i[0])*i[-1]
					out_rd[0][i_block[1]]+=(i[1]-temp_bp[0][i_block[1]])*i[-1]
			if i[2]=='p':
				i_block=[]
				for k in i[:2]:
					if k<temp_bp[1][1]:
						i_block.append(0)
					elif k>temp_bp[1][-2]-1:
						i_block.append(len(temp_bp[1])-2)
					else:
						for j in range(len(temp_bp[1])-1)[1:-1]:
							if temp_bp[1][j]-1<k and temp_bp[1][j+1]>k:
								i_block.append(j)
				if i_block[0]==i_block[1]:
					out_rd[1][i_block[0]]+=(i[1]-i[0])*i[-1]
				elif not i_block[0]==i_block[1]:
					out_rd[1][i_block[0]]+=(temp_bp[1][i_block[0]+1]-i[0])*i[-1]
					out_rd[1][i_block[1]]+=(i[1]-temp_bp[1][i_block[1]])*i[-1]
	block_bps_chr={}
	block_bps_chr['m']={}
	block_bps_chr['p']={}
	if not Penalty_For_InsertLengthZero in IL_Rec.keys():
		IL_Rec[Penalty_For_InsertLengthZero]=NoMapPenal
	else:
		IL_Rec[Penalty_For_InsertLengthZero]+=NoMapPenal
	IL_Penal=0
	IL_Weight=0
	for i in IL_Rec.keys():
		IL_Penal+=i*IL_Rec[i]
		IL_Weight+=IL_Rec[i]
	if not IL_Weight==0:
		IL_Output=IL_Penal/IL_Weight
	else:
		IL_Output=0
	Num_Read_TB=[out_tb[0][1:-1],out_tb[1][1:-1]]
	TB_Pena_2_out=0
	Num_total_TB=[]
	for x in Num_Read_TB:
		Num_total_TB+=x
	if numpy.sum(Num_total_TB)>0:
		pvalue=scipy.stats.chisquare(Num_total_TB)[1]
	else:
		pvalue=0.0
	if pvalue>0:
		TB_Pena_2_out=numpy.log(pvalue)
	else:
		TB_Pena_2_out=-100000000
	#for i in Num_Read_TB:
	#	for j in i:
	#		if j<Through_BP_Minimum:
	#			TB_Pena_2_out+=1
	Af_Block_Len=[[flank]+[Af_BP[0][i+1]-Af_BP[0][i] for i in range(len(Af_BP[0])-1)]+[flank],[flank]+[Af_BP[1][i+1]-Af_BP[1][i] for i in range(len(Af_BP[1])-1)]+[flank]]
	out_rd=[[out_rd[0][i]/Af_Block_Len[0][i] for i in range(len(out_rd[0]))],[out_rd[1][i]/Af_Block_Len[1][i] for i in range(len(out_rd[1]))]]
	out_rd_new=[[(RD_within_B['left']-out_rd[0][0]-out_rd[1][0])/2.0+out_rd[0][0],
	(RD_within_B['right']-out_rd[0][-1]-out_rd[1][-1])/2.0+out_rd[0][-1]],
	[(RD_within_B['left']-out_rd[0][0]-out_rd[1][0])/2.0+out_rd[1][0],
	(RD_within_B['right']-out_rd[0][-1]-out_rd[1][-1])/2.0+out_rd[1][-1]]]
	out_rd=[[out_rd_new[0][0]]+out_rd[0][1:-1]+[out_rd_new[0][-1]],[out_rd_new[1][0]]+out_rd[1][1:-1]+[out_rd_new[1][-1]]]
	out_rd_within=[[RD_within_B[Af_Letter[0][i]]/letters_numbers[0][i] for i in range(len(Af_Letter[0]))],[RD_within_B[Af_Letter[1][i]]/letters_numbers[1][i] for i in range(len(Af_Letter[1]))]]
	out_rd_within[0]=[0]+out_rd_within[0]+[0]
	out_rd_within[1]=[0]+out_rd_within[1]+[0]	
	cov_bp2=[[out_rd[0][i]+out_rd_within[0][i] for i in range(len(out_rd[0]))],[out_rd[1][i]+out_rd_within[1][i] for i in range(len(out_rd[1]))]]
	Cov_GC=[[BlockGC[k] for k in Af_Letter[0]],[BlockGC[k] for k in Af_Letter[1]]]
	adj_cov_bp=[GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,chrom_N,Cov_GC[0],cov_bp2[0][1:-1]),GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,chrom_N,Cov_GC[1],cov_bp2[1][1:-1])]
	return [IL_Output,adj_cov_bp,DR_Penal,TB_Pena_2_out,Num_total_TB]

def TB_Cal_BP(pos,bp,tblist,perc):
	for i in range(len(bp)-2)[1:]:
		if not pos[0]<bp[i] and pos[1]<bp[i+1]:
			break
		elif not pos[0]<bp[i] and pos[0]<bp[i+1] and not pos[1]<bp[i+1] and pos[1]<bp[i+2]:
			tblist[i+1]+=perc
		elif not pos[0]<bp[i] and pos[0]<bp[i+1] and not pos[1]<bp[i+2] and pos[1]<bp[i+3]:
			tblist[i+1]+=perc
			tblist[i+2]+=perc
	i=len(bp)-2
	if not pos[0]<bp[i] and pos[0]<bp[i+1] and not pos[1]<bp[i+1]:
		tblist[i+1]+=perc

def Be_Info_1_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal):
	be_info_1=Be_Info[0]
	for j in be_info_1:
		jMapPenam=0
		j_m_new=[]
		if j[0] in temp_letter[0] and j[3] in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]]:
				for kb in Let_BP_Info['m'][j[3]]:
					j_m_temp=[j[1]+ka[0],j[2]+ka[0],j[4]+kb[0],j[5]+kb[0]]
					if j_m_temp[0]>j_m_temp[2]:
						j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[j[-1],j[-2]]
					else:
						j_m_temp+=[j[-2],j[-1]]
					j_m_new.append(j_m_temp)
		if j[0]+'^' in temp_letter[0] and j[3] in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]+'^']:
				for kb in Let_BP_Info['m'][j[3]]:
					j_m_temp=[ka[1]-j[2],ka[1]-j[1],kb[0]+j[4],kb[0]+j[5]]
					if j_m_temp[0]>j_m_temp[2]:
						j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[j[-1],complement(j[-2])]
					else:
						j_m_temp+=[complement(j[-2]),j[-1]]
					j_m_new.append(j_m_temp)
		if j[0] in temp_letter[0] and j[3]+'^' in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]]:
				for kb in  Let_BP_Info['m'][j[3]+'^']:
					j_m_temp=[j[1]+ka[0],j[2]+ka[0],kb[1]-j[5],kb[1]-j[4]]
					if j_m_temp[0]>j_m_temp[2]:
						j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[complement(j[-1]),j[-2]]
					else:
						j_m_temp+=[j[-2],complement(j[-1])]
					j_m_new.append(j_m_temp)
		if j[0]+'^' in temp_letter[0] and j[3]+'^' in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]+'^']:
				for kb in  Let_BP_Info['m'][j[3]+'^']:
					j_m_temp=[ka[1]-j[2],ka[1]-j[1],kb[1]-j[5],kb[1]-j[4]]
					if j_m_temp[0]>j_m_temp[2]:
						j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[complement(j[-1]),complement(j[-2])]
					else:
						j_m_temp+=[complement(j[-2]),complement(j[-1])]
					j_m_new.append(j_m_temp)
		j_m_3a=candidate_QC_Control(j_m_new)
		if j_m_3a==[]:
			jMapPenam+=1
		j_p_new=[]
		if j[0] in temp_letter[1] and j[3] in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]]:
				for kb in Let_BP_Info['p'][j[3]]:
					j_p_temp=[j[1]+ka[0],j[2]+ka[0],j[4]+kb[0],j[5]+kb[0]]
					if j_p_temp[0]>j_p_temp[2]:
						j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[j[-1],j[-2]]
					else:
						j_p_temp+=[j[-2],j[-1]]
					j_p_new.append(j_p_temp)
		if j[0]+'^' in temp_letter[1] and j[3] in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]+'^']:
				for kb in Let_BP_Info['p'][j[3]]:
					j_p_temp=[ka[1]-j[2],ka[1]-j[1],kb[0]+j[4],kb[0]+j[5]]
					if j_p_temp[0]>j_p_temp[2]:
						j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[j[-1],complement(j[-2])]
					else:
						j_p_temp+=[complement(j[-2]),j[-1]]
					j_p_new.append(j_p_temp)
		if j[0] in temp_letter[1] and j[3]+'^' in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]]:
				for kb in Let_BP_Info['p'][j[3]+'^']:
					j_p_temp=[j[1]+ka[0],j[2]+ka[0],kb[1]-j[5],kb[1]-j[4]]
					if j_p_temp[0]>j_p_temp[2]:
						j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[complement(j[-1]),j[-2]]
					else:
						j_p_temp+=[j[-2],complement(j[-1])]
					j_p_new.append(j_p_temp)
		if j[0]+'^' in temp_letter[1] and j[3]+'^' in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]+'^']:
				for kb in  Let_BP_Info['p'][j[3]+'^']:
					j_p_temp=[ka[1]-j[2],ka[1]-j[1],kb[1]-j[5],kb[1]-j[4]]
					if j_p_temp[0]>j_p_temp[2]:
						j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[complement(j[-1]),complement(j[-2])]
					else:
						j_p_temp+=[complement(j[-2]),complement(j[-1])]
					j_p_new.append(j_p_temp)
		j_p_3a=candidate_QC_Control(j_p_new)
		if j_p_3a==[]:
			jMapPenam+=1
		if jMapPenam==2:
			Total_Cov_For_Pen[j[0]]+=j[2]-j[1]
			Total_Cov_For_Pen[j[3]]+=j[5]-j[4]
			NoMapPenal+=2
		elif jMapPenam==1:
			if j_m_3a==[]:
				Map_P+=[jp3+['p']+[float(1)/float(len(j_p_3a))] for jp3 in j_p_3a]
			elif j_p_3a==[]:
				Map_M+=[jp3+['m']+[float(1)/float(len(j_m_3a))] for jp3 in j_m_3a]
		else:
			j_mp_4a=candidate_QC_Control2(j_m_3a,j_p_3a)
			if not j_mp_4a==[]:
				Map_Both+=[j4+[float(1)/float(len(j_mp_4a))] for j4 in j_mp_4a]
			else:
				Total_Cov_For_Pen[j[0]]+=j[2]-j[1]
				Total_Cov_For_Pen[j[3]]+=j[5]-j[4]
				NoMapPenal+=2
	return NoMapPenal

def Be_Info_2_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal):
	be_info_2=Be_Info[1]
	for j in be_info_2:
		jMapPenam=0
		j_m_new=[]
		if j[0] in temp_letter[0] and j[2] in temp_letter[0] and j[4] in temp_letter[0] and j[6] in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]]:
				for kb in Let_BP_Info['m'][j[2]]:
					for kc in Let_BP_Info['m'][j[4]]:
						for kd in Let_BP_Info['m'][j[6]]:
							j_info_new=[ka[0]+j[1],kb[0]+j[3],kc[0]+j[5],kd[0]+j[7]]
							if j_info_new[0]>j_info_new[2]:
								j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],j[-2]])
							else:
								j_m_new.append(j_info_new+[j[-2],j[-1]])
		if j[0]+'^' in temp_letter[0] and j[2]+'^' in temp_letter[0] and j[4] in temp_letter[0] and j[6] in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]+'^']:
				for kb in Let_BP_Info['m'][j[2]+'^']:
					for kc in Let_BP_Info['m'][j[4]]:
						for kd in Let_BP_Info['m'][j[6]]:
							j_info_new=[kb[1]-j[3],ka[1]-j[1],kc[0]+j[5],kd[0]+j[7]]
							if j_info_new[0]>j_info_new[2]:
								j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],complement(j[-2])])
							else:
								j_m_new.append(j_info_new+[complement(j[-2]),j[-1]])
		if j[0] in temp_letter[0] and j[2] in temp_letter[0] and j[4]+'^' in temp_letter[0] and j[6]+'^' in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]]:
				for kb in Let_BP_Info['m'][j[2]]:
					for kc in Let_BP_Info['m'][j[4]+'^']:
						for kd in Let_BP_Info['m'][j[6]+'^']:
							j_info_new=[ka[0]+j[1],kb[0]+j[3],kd[1]-j[7],kc[1]-j[5]]
							if j_info_new[0]>j_info_new[2]:
								j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),j[-2]])
							else:
								j_m_new.append(j_info_new+[j[-2],complement(j[-1])])
		if j[0]+'^' in temp_letter[0] and j[2]+'^' in temp_letter[0] and j[4]+'^' in temp_letter[0] and j[6]+'^' in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]+'^']:
				for kb in Let_BP_Info['m'][j[2]+'^']:
					for kc in Let_BP_Info['m'][j[4]+'^']:
						for kd in Let_BP_Info['m'][j[6]+'^']:
							j_info_new=[kb[1]-j[3],ka[1]-j[1],kd[1]-j[7],kc[1]-j[5]]
							if j_info_new[0]>j_info_new[2]:
								j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),complement(j[-2])])
							else:
								j_m_new.append(j_info_new+[complement(j[-2]),complement(j[-1])])
		j_m_3a=candidate_QC_Control(j_m_new)
		if j_m_3a==[]:
			jMapPenam+=1
		j_p_new=[]
		if j[0] in temp_letter[1] and j[2] in temp_letter[1] and j[4] in temp_letter[1] and j[6] in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]]:
				for kb in Let_BP_Info['p'][j[2]]:
					for kc in Let_BP_Info['p'][j[4]]:
						for kd in Let_BP_Info['p'][j[6]]:
							j_info_new=[ka[0]+j[1],kb[0]+j[3],kc[0]+j[5],kd[0]+j[7]]
							if j_info_new[0]>j_info_new[2]:
								j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],j[-2]])
							else:
								j_p_new.append(j_info_new+[j[-2],j[-1]])
		if j[0]+'^' in temp_letter[1] and j[2]+'^' in temp_letter[1] and j[4] in temp_letter[1] and j[6] in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]+'^']:
				for kb in Let_BP_Info['p'][j[2]+'^']:
					for kc in Let_BP_Info['p'][j[4]]:
						for kd in Let_BP_Info['p'][j[6]]:
							j_info_new=[kb[1]-j[3],ka[1]-j[1],kc[0]+j[5],kd[0]+j[7]]
							if j_info_new[0]>j_info_new[2]:
								j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],complement(j[-2])])
							else:
								j_p_new.append(j_info_new+[complement(j[-2]),j[-1]])
		if j[0] in temp_letter[1] and j[2] in temp_letter[1] and j[4]+'^' in temp_letter[1] and j[6]+'^' in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]]:
				for kb in Let_BP_Info['p'][j[2]]:
					for kc in Let_BP_Info['p'][j[4]+'^']:
						for kd in Let_BP_Info['p'][j[6]+'^']:
							j_info_new=[ka[0]+j[1],kb[0]+j[3],kd[1]-j[7],kc[1]-j[5]]
							if j_info_new[0]>j_info_new[2]:
								j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),j[-2]])
							else:
								j_p_new.append(j_info_new+[j[-2],complement(j[-1])])
		if j[0]+'^' in temp_letter[1] and j[2]+'^' in temp_letter[1] and j[4]+'^' in temp_letter[1] and j[6]+'^' in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]+'^']:
				for kb in Let_BP_Info['p'][j[2]+'^']:
					for kc in Let_BP_Info['p'][j[4]+'^']:
						for kd in Let_BP_Info['p'][j[6]+'^']:
							j_info_new=[kb[1]-j[3],ka[1]-j[1],kd[1]-j[7],kc[1]-j[5]]
							if j_info_new[0]>j_info_new[2]:
								j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),complement(j[-2])])
							else:
								j_p_new.append(j_info_new+[complement(j[-2]),complement(j[-1])])
		j_p_3a=candidate_QC_Control(j_p_new)
		if j_p_3a==[]:
			jMapPenam+=1
		if jMapPenam==2:
			if j[0]==j[2]:
				Total_Cov_For_Pen[j[0]]+=j[3]-j[1]
			else:
				Total_Cov_For_Pen[j[0]]+=Be_BP_Letter[j[0]]-j[1]
				Total_Cov_For_Pen[j[2]]+=j[3]
			if j[4]==j[6]:
				Total_Cov_For_Pen[j[4]]+=j[7]-j[5]
			else:
				Total_Cov_For_Pen[j[4]]+=Be_BP_Letter[j[4]]-j[5]
				Total_Cov_For_Pen[j[6]]+=j[7]
			NoMapPenal+=2
		elif jMapPenam==1:
			if j_m_3a==[]:
				Map_P+=[jp3+['p']+[float(1)/float(len(j_p_3a))] for jp3 in j_p_3a]
			elif j_p_3a==[]:
				Map_M+=[jp3+['m']+[float(1)/float(len(j_m_3a))] for jp3 in j_m_3a]
		else:
			j_mp_4a=candidate_QC_Control2(j_m_3a,j_p_3a)
			if not j_mp_4a==[]:
				Map_Both+=[j4+[float(1)/float(len(j_mp_4a))] for j4 in j_mp_4a]
			else:
				if j[0]==j[2]:
					Total_Cov_For_Pen[j[0]]+=j[3]-j[1]
				else:
					Total_Cov_For_Pen[j[0]]+=Be_BP_Letter[j[0]]-j[1]
					Total_Cov_For_Pen[j[2]]+=j[3]
				if j[4]==j[6]:
					Total_Cov_For_Pen[j[4]]+=j[7]-j[5]
				else:
					Total_Cov_For_Pen[j[4]]+=Be_BP_Letter[j[4]]-j[5]
					Total_Cov_For_Pen[j[6]]+=j[7]
				NoMapPenal+=2
	return NoMapPenal

def Be_Info_3_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal):
	be_info_3=Be_Info[2]
	for j in be_info_3:
		j_m_new=[]
		if j[0] in temp_letter[0] and j[2] in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]]:
				for kb in Let_BP_Info['m'][j[2]]:
					temp_single=[ka[0]+j[1],kb[0]+j[3]]
					if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
						j_m_new.append(temp_single)
		if j[0]+'^' in temp_letter[0] and j[2]+'^' in temp_letter[0]:
			for ka in Let_BP_Info['m'][j[0]+'^']:
				for kb in Let_BP_Info['m'][j[2]+'^']:
					temp_single=[kb[1]-j[3],ka[1]-j[1]]
					if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
						j_m_new.append(temp_single)
		j_p_new=[]
		if j[0] in temp_letter[1] and j[2] in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]]:
				for kb in Let_BP_Info['p'][j[2]]:
					temp_single=[ka[0]+j[1],kb[0]+j[3]]
					if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
						j_p_new.append(temp_single)
		if j[0]+'^' in temp_letter[1] and j[2]+'^' in temp_letter[1]:
			for ka in Let_BP_Info['p'][j[0]+'^']:
				for kb in Let_BP_Info['p'][j[2]+'^']:
					temp_single=[kb[1]-j[3],ka[1]-j[1]]
					if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
						j_p_new.append(temp_single)
		if not j_m_new+j_p_new==[]:
			for j2 in j_m_new:
				Map_Both.append(j2+['m',float(1)/float(len(j_m_new+j_p_new))])
			for j2 in j_p_new:
				Map_Both.append(j2+['p',float(1)/float(len(j_m_new+j_p_new))])
		else:
			Total_Cov_For_Pen[j[0]]=Be_BP_Letter[j[0]]-j[1]
			Total_Cov_For_Pen[j[2]]=j[3]
			NoMapPenal+=1
	return NoMapPenal

def Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC,RD_within_B):
	Total_Cov_For_Pen={}
	for key in RD_within_B.keys():
		Total_Cov_For_Pen[key]=0
	Map_M=[]
	Map_P=[]
	Map_Both=[]
	Let_BP_Info={}
	Let_BP_Info['m']={}
	Let_BP_Info['p']={}
	temp_letter=[['left']+Af_Letter[0]+['right'],['left']+Af_Letter[1]+['right']]
	temp_bp=[[Af_BP[0][0]-flank]+Af_BP[0]+[Af_BP[0][-1]+flank],[Af_BP[1][0]-flank]+Af_BP[1]+[Af_BP[1][-1]+flank]]
	for j1 in range(len(temp_letter[0])):
		j=temp_letter[0][j1]
		if not j in Let_BP_Info['m'].keys():
				Let_BP_Info['m'][j]=[[temp_bp[0][j1],temp_bp[0][j1+1]]]
		else:
				Let_BP_Info['m'][j]+=[[temp_bp[0][j1],temp_bp[0][j1+1]]]
	for j1 in range(len(temp_letter[1])):
		j=temp_letter[1][j1]
		if not j in Let_BP_Info['p'].keys():
			Let_BP_Info['p'][j]=[[temp_bp[1][j1],temp_bp[1][j1+1]]]
		else:
			Let_BP_Info['p'][j]+=[[temp_bp[1][j1],temp_bp[1][j1+1]]]
	letters_numbers=[[Af_Letter[0].count(i[0])+Af_Letter[1].count(i[0])+Af_Letter[0].count(i[0]+'^')+Af_Letter[1].count(i[0]+'^') for i in Af_Letter[0]],[Af_Letter[0].count(i[0])+Af_Letter[1].count(i[0])+Af_Letter[0].count(i[0]+'^')+Af_Letter[1].count(i[0]+'^') for i in Af_Letter[1]]]
	NoMapPenal=0
	IL_Rec={}
	DR_Rec=0
	cov_bp=[[0 for i in range(len(temp_letter[0]))],[0 for i in range(len(temp_letter[1]))]]
	cov_bp2=[]
	NoMapPenal=Be_Info_1_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal)
	NoMapPenal=Be_Info_2_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal)
	NoMapPenal=Be_Info_3_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal)
	best_structure_sign_flag=0
	for key in Total_Cov_For_Pen.keys():
		if Total_Cov_For_Pen[key]==0:
			del Total_Cov_For_Pen[key]
		else:
			Total_Cov_For_Pen[key]/=float(Be_BP_Letter[key])
	for key in RD_within_B.keys():
		if not key[-1]=='^' and not key in ['left','right','left^', 'right^']:
			if not key in Af_Letter[0]+Af_Letter[1] and not key+'^' in Af_Letter[0]+Af_Letter[1]:
				if not key in Total_Cov_For_Pen.keys():
					Total_Cov_For_Pen[key]=0
				Total_Cov_For_Pen[key]+=RD_within_B[key]
	if NoMapPenal>0:
		best_structure_sign_flag+=1
	for key1 in Total_Cov_For_Pen.keys():
		if Total_Cov_For_Pen[key1]>2.58*GC_Std_Coverage[chrom_N]:
			best_structure_sign_flag+=1
	if not Map_M+Map_P+Map_Both==[]:
		penals=penal_calculate(Map_M+Map_P+Map_Both,temp_bp,Af_Letter,Af_BP,RD_within_B,letters_numbers,BlockGC,NoMapPenal)
		if penals[2]>0:
			best_structure_sign_flag+=1
		return penals[:-1]+[NoMapPenal,Total_Cov_For_Pen,best_structure_sign_flag]+[penals[-1]]
	else:
		return 0

def Move_Decide_2(IL_List,RD_List,GC_Overall_Median_Coverage,GC_Var_Coverage):
	#here, we add direction as an extra part to IL penalty
	IL_Weight=1
	RD_Weight=1
	regulator=-numpy.max([(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight) for j in range(len(IL_List))])/5
	T_Penal=[(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight)/regulator for j in range(len(IL_List))]
	T2_Penal=[math.exp(k) for k in T_Penal]
	Normalized_Penal=[l/numpy.sum(T2_Penal) for l in T2_Penal]
	indicator=float(random.choice(range(1000)))/1000
	for j in range(len(IL_List)):
		cdf=numpy.sum(Normalized_Penal[:j+1])
		if cdf>indicator or cdf==indicator:
			return [j,IL_List[j]*IL_Weight+RD_List[j]*RD_Weight]

def Move_Decide_3(IL_List,RD_List,GC_Overall_Median_Coverage,GC_Var_Coverage):
	#pick the highest scored structure as the final
	IL_Weight=1
	RD_Weight=1
	regulator=-numpy.max([(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight) for j in range(len(IL_List))])/5
	T_Penal=[(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight)/regulator for j in range(len(IL_List))]
	T2_Penal=[math.exp(k) for k in T_Penal]
	Normalized_Penal=[l/numpy.sum(T2_Penal) for l in T2_Penal]
	j=Normalized_Penal.index(max(Normalized_Penal))
	return [j,IL_List[j]*IL_Weight+RD_List[j]*RD_Weight]

def Af_Letter_QC(Af_Letter,Copy_num_estimate):
	Copy_Num_Real={}
	for i in Copy_num_estimate.keys():
		Copy_Num_Real[i]=0
	for i in Af_Letter:
		for j in i:
			Copy_Num_Real[j[0]]+=1
	FLAG=0
	for i in Copy_num_estimate.keys():
		if abs(Copy_num_estimate[i]-Copy_Num_Real[i])>1:
			FLAG+=1
	return FLAG

def getPwd(length):
	s=''.join([chr(97+i) for i in range(length)])
	sSet=set(s)
	if len(sSet)<length:return 0
	result=tmp=list(sSet)
	for i in range(length-1):
		result=createPwd(result,tmp)
	return result

def createPwd(sList,kList):
	assert len(sList)*len(kList)>0
	upper,sList=len(sList),sList*len(kList)
	for i in range(len(kList)):
		for j in range(upper):
			sList[i*upper+j]+=kList[i]
	return sList

def changeLet(num):
	if num>1:
		let1=getPwd(num)
		let2=[]
		for k1 in let1:
			let2.append([])
			for k2 in k1:
				let2[-1].append(chr(97+(ord(k2)-97)%2))
		let3=[]
		for k1 in let2:
			if not k1 in let3:
				let3.append(k1)
		let4=[]
		for k1 in let3:
			temp=[]
			for k2 in k1:
				if k2=='a':
					temp.append('a')
				else:
					temp.append('a^')
			let4.append(temp)
		return let4
	elif num==1:
		return[['a'],['a^']]
	elif num==0:
		return [[]]

def struc_propose_single_block(num):
	if num<1:
		return 'error'
	else:
		struc_rec=[]
		for i in range(num/2+1):
			rec_a=changeLet(i)
			rec_b=changeLet(num-i)
			for k1 in rec_a:
				for k2 in rec_b:
					if not sorted([k1,k2]) in struc_rec:
						struc_rec.append(sorted([k1,k2]))
		return struc_rec

def bp_to_let(del_info_unit,chromos):
	flag=0
	for i in del_info_unit:
		if i  in chromos:
			flag+=1
	letter=''.join([chr(i+97) for i in range(len(del_info_unit)-2*flag)])
	letters='/'.join([letter,letter])
	return letters

def write_best_letter(bps_all,Best_Letter_Rec,Best_Score_Rec,Score_rec_hash,original_letters):
	fo=open(output_Score_File,'a')
	time2=time.time()
	Best_Letter_2=[]
	if not Score_rec_hash=={}:
		temp1=Best_Let_modify(Best_Letter_Rec,Best_Score_Rec,Score_rec_hash)
		Best_Letter_Rec=temp1[0]
		Best_Score_Rec=temp1[1]
	for bestletter in Best_Letter_Rec:
		if not sorted(bestletter) in Best_Letter_2:
			Best_Letter_2.append(sorted(bestletter))
	bps3=[]
	for bps in bps_all:
		bps3+=bps
	for bestletter in Best_Letter_2:
		if not '/'.join([''.join(original_letters),''.join(original_letters)])=='/'.join([''.join(bestletter[0]),''.join(bestletter[1])]):
			print >>fo, ' '.join([str(bp_ele) for bp_ele in bps3])
			print >>fo, '/'.join([''.join(bestletter[0]),''.join(bestletter[1])])
			print >>fo, 'Theoretical Best Score: '+str(Best_IL_Score+Best_RD_Score)
			print >>fo, 'Current Best Scure: '+str(Best_Score_Rec)
			print>>fo, 'Time Consuming:'+str(datetime.timedelta(seconds=(time2-time1)))
	fo.close()

def zero_RD_Process(run_flag):
	Best_Letter_Rec=[[[], []]]
	Af_Letter=[[],[]]
	Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
	Best_Score_Rec=Best_IL_Score+Best_RD_Score
	run_flag+=1
	return([Best_Letter_Rec,Best_Score_Rec,run_flag])

def Af_TB_Penal_Less_Important_Caldu(Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num):
	Af_TB_Rec=0
	Af_TB_Rec_list=Af_Info_all[-1]
	Af_TB_Penal_a=Af_Info_all[4]
	for x in Af_TB_Rec_list:
		if x <(IL_Statistics[0]*IL_Statistics[4]+IL_Statistics[1]*IL_Statistics[5])/ReadLength*GC_Overall_Median_Num/float(10):
			Af_TB_Rec+=1
	return -float(Af_TB_Penal_a)/float(num_of_reads)-float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)

def Af_TB_Penal_More_Important_Caldu(Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num):
	Af_TB_Penal=Af_Info_all[3]
	Af_TB_Penal_a=Af_Info_all[4]
	return -float(Af_TB_Penal_a)/float(num_of_reads)-Af_TB_Penal

def one_RD_Process(run_flag,Score_rec_hash):
	if MoveFlag==2:
		Letter_Candidates=[[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']],[['a','a'], []],[['a','a^'], []],[['a^','a'], []],[['a^','a^'], []]]
	elif MoveFlag==1:
		Letter_Candidates=[i for i in [[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']],[['a','a'], []],[['a','a^'], []],[['a^','a'], []],[['a^','a^'], []]] if i[0]==i[1]]
	elif MoveFlag==0:
		Letter_Candidates=[i for i in [[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']],[['a','a'], []],[['a','a^'], []],[['a^','a'], []],[['a^','a^'], []]] if ['a'] in i]
	P_IL=[]
	P_RD=[]
	P_DR=[]
	P_TB=[]
	Letter_Rec=[]
	BP_Rec=[]
	for Af_Letter in Letter_Candidates:
		Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
		for i in Af_Letter[0]:
				Af_BP[0].append(Af_BP[0][-1]+Be_BP_Letter[i])
		for i in Af_Letter[1]:
				Af_BP[1].append(Af_BP[1][-1]+Be_BP_Letter[i])
		Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
		if not Af_Info_all==0:
			Letter_Rec.append(Af_Letter)
			BP_Rec.append(Af_BP)
			Af_IL_Penal=Af_Info_all[0]
			Af_RD_Rec=Af_Info_all[1]
			Af_DR_Penal=(Af_Info_all[2])**2
			Af_TB_Penal_a=Af_Info_all[4]
			Af_TB_Penal=Af_TB_Penal_Less_Important_Caldu(Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num)
			Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
			for key in Af_Info_all[5].keys():
					Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
			P_IL.append(Af_IL_Penal-IL_GS)
			P_RD.append((Af_RD_Penal-RD_GS)*(len(Af_Letter[0])+len(Af_Letter[1])))
			P_DR.append(Af_DR_Penal/num_of_read_pairs)
			P_TB.append(Af_TB_Penal)
	Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
	Regu_IL=[i*K_IL_new for i in Regu_IL]
	Regu_RD=[P_RD[i]*(1-P_TB[i]) for i in range(len(P_RD))]	
	Regulator=1
	ILTemp=[j/Regulator for j in Regu_IL]
	RDTemp=[i for i in Regu_RD]
	if not ILTemp==[]:
		DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
		Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
		Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
		run_flag+=1
		for x in range(len(Letter_Rec)):
			xy=ILTemp[x]+RDTemp[x]
			if not xy in Score_rec_hash.keys():
				Score_rec_hash[xy]=[]
			Score_rec_hash[xy].append(Letter_Rec[x])			   
	else:
		Best_Letter_Rec=[]
		Best_Score_Rec=100
		run_flag+=1
	Score_rec_hash2=score_rec_hash_Modify_for_short_del(Score_rec_hash)
	return([Best_Letter_Rec,Best_Score_Rec,run_flag,Score_rec_hash2])

def score_rec_hash_Modify_for_short_del(Score_rec_hash):
	Score_rec_hash_new={}
	for x in sorted(Score_rec_hash.keys())[::-1][:1]:
		Score_rec_hash_new[x]=Score_rec_hash[x]
	for x in sorted(Score_rec_hash.keys())[::-1][1:]:
		Score_rec_hash_new[x-1.1]=Score_rec_hash[x]
	return Score_rec_hash_new

def two_RD_Process(run_flag,Score_rec_hash):
	Letter_Candidates=struc_propose_single_block(2)+struc_propose_single_block(3)+struc_propose_single_block(4)+struc_propose_single_block(5)	
	if MoveFlag==2:
		Letter_Candidates=Letter_Candidates
	elif MoveFlag==1:
		Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
	elif MoveFlag==0:
		Letter_Candidates=[i for i in Letter_Candidates if ['a'] in i]
	P_IL=[]
	P_RD=[]
	P_DR=[]
	P_TB=[]
	Letter_Rec=[]
	BP_Rec=[]
	for Af_Letter in Letter_Candidates:
		#print Af_Letter
		Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
		for i in Af_Letter[0]:
				Af_BP[0].append(Af_BP[0][-1]+Be_BP_Letter[i])
		for i in Af_Letter[1]:
				Af_BP[1].append(Af_BP[1][-1]+Be_BP_Letter[i])
		Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
		if not Af_Info_all==0:
			Letter_Rec.append(Af_Letter)
			BP_Rec.append(Af_BP)
			Af_IL_Penal=Af_Info_all[0]
			Af_RD_Rec=Af_Info_all[1]
			Af_DR_Penal=(Af_Info_all[2])**2
			Af_TB_Penal_a=Af_Info_all[4]
			Af_TB_Rec=Af_Info_all[3]
			Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
			Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
			for key in Af_Info_all[5].keys():
					Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
			P_IL.append(Af_IL_Penal)
			P_RD.append(Af_RD_Penal)
			P_DR.append(Af_DR_Penal/num_of_read_pairs)
			P_TB.append(Af_TB_Penal)
	Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
	#Regu_RD=[P_RD[i]*(1+TB_Weight*P_TB[i]) for i in range(len(P_RD))]
	Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
	Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
	Regu_RD=[i-RD_GS for i in Regu_RD]
	#Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
	Regulator=1
	ILTemp=[j/Regulator for j in Regu_IL]
	RDTemp=[i for i in Regu_RD]
	if not ILTemp==[]:
		DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
		Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
		Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
		run_flag+=1
		for x in range(len(Letter_Rec)):
			xy=ILTemp[x]+RDTemp[x]
			if not xy in Score_rec_hash.keys():
				Score_rec_hash[xy]=[]
			Score_rec_hash[xy].append(Letter_Rec[x])			   
	else:
		Best_Letter_Rec=[]
		Best_Score_Rec=100
		run_flag+=1						
	return([Best_Letter_Rec,Best_Score_Rec,run_flag,Score_rec_hash])

def few_RD_Process(run_flag,Score_rec_hash):
	Letter_Candidates=struc_propose_single_block(copy_num_a)+struc_propose_single_block(copy_num_b)
	if MoveFlag==2:
		Letter_Candidates=Letter_Candidates
	elif MoveFlag==1:
		Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
	elif MoveFlag==0:
		Letter_Candidates=[i for i in Letter_Candidates if ['a'] in i]
	P_IL=[]
	P_RD=[]
	P_DR=[]
	P_TB=[]
	Letter_Rec=[]
	BP_Rec=[]
	for Af_Letter in Letter_Candidates:
		Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
		for i in Af_Letter[0]:
				Af_BP[0].append(Af_BP[0][-1]+Be_BP_Letter[i])
		for i in Af_Letter[1]:
				Af_BP[1].append(Af_BP[1][-1]+Be_BP_Letter[i])
		Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
		if not Af_Info_all==0:
			Letter_Rec.append(Af_Letter)
			BP_Rec.append(Af_BP)
			Af_IL_Penal=Af_Info_all[0]
			Af_RD_Rec=Af_Info_all[1]
			Af_DR_Penal=(Af_Info_all[2])**2
			Af_TB_Penal_a=Af_Info_all[4]
			Af_TB_Rec=Af_Info_all[3]
			Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
			Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
			for key in Af_Info_all[5].keys():
					Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
			P_IL.append(Af_IL_Penal)
			P_RD.append(Af_RD_Penal)
			P_DR.append(Af_DR_Penal/num_of_read_pairs)
			P_TB.append(Af_TB_Penal)
	Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
	Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
	Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
	Regu_RD=[i-RD_GS for i in Regu_RD]
	#Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
	Regulator=1
	ILTemp=[j/Regulator for j in Regu_IL]
	RDTemp=[i for i in Regu_RD]
	if not ILTemp==[]:
		DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
		Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
		Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
		run_flag+=1
		for x in range(len(Letter_Rec)):
			xy=ILTemp[x]+RDTemp[x]
			if not xy in Score_rec_hash.keys():
				Score_rec_hash[xy]=[]
			Score_rec_hash[xy].append(Letter_Rec[x])			   
	else:
		Best_Letter_Rec=[]
		Best_Score_Rec=100
		run_flag+=1						
	return([Best_Letter_Rec,Best_Score_Rec,run_flag,Score_rec_hash])

def many_RD_Process(run_flag):
	Best_Letter_Rec=[[['a' for i in range(copy_num_a/2)],['a' for i in range(copy_num_a/2)]]]
	Best_Score_Rec=100
	run_flag+=1		
	return([Best_Letter_Rec,Best_Score_Rec,run_flag])

def two_block_a(cpNum_a,cpNum_b):
	out=[]
	for k1 in cpNum_a:
		for k2 in cpNum_b:
			out.append(['a' for k in range(int(k1))]+['b' for k in range(int(k2))])
	out2=[]
	for k1 in out:
		t1=[]
		ta=[]
		tb=[]
		for k2 in range(k1.count('a')+1):
			ta.append([['a' for k in range(k2)],['a' for k in range(k1.count('a')-k2)]])
		for k2 in range(k1.count('b')+1):
			tb.append([['b' for k in range(k2)],['b' for k in range(k1.count('b')-k2)]])
		for k2 in ta:
			for k3 in tb:
				if not sorted([k2[0]+k3[0],k2[1]+k3[1]]) in t1:
					t1.append(sorted([k2[0]+k3[0],k2[1]+k3[1]]))
		out2+=t1
	return out2

def two_block_b(all_Strucs):
	out=[]
	for k1 in all_Strucs:
		k1a=[]
		k1b=[]
		rect=0
		if k1[0]==[]:
			k1a2=[[]]
		else:
			for k2 in range(len(k1[0])+1):
				all_lets=k1[0]+['^' for ki in range(k2)]
				for k3 in itertools.permutations(all_lets,len(all_lets)):
					if not k3[0]=='^' and not k3 in k1a:
						k3flag=0
						for k4 in range(len(k3[:-1])):
							if k3[k4+1]==k3[k4]=='^':
								k3flag+=1
						if k3flag==0:
							k1a.append(k3)
			k1a2=[]
			for k2 in k1a:
				tk2=[]
				for k3 in k2:
					if k3=='^':
						tk2[-1]+=k3
					else:
						tk2.append(k3)
				k1a2.append(tk2)
		if k1[1]==[]:
			k1b2=[[]]
		else:
			for k2 in range(len(k1[1])+1):
				all_lets=k1[1]+['^' for ki in range(k2)]
				for k3 in itertools.permutations(all_lets,len(all_lets)):
					if not k3[0]=='^' and not k3 in k1b:
						k3flag=0
						for k4 in range(len(k3[:-1])):
							if k3[k4+1]==k3[k4]=='^':
								k3flag+=1
						if k3flag==0:
							k1b.append(k3)
			k1b2=[]
			for k2 in k1b:
				tk2=[]
				for k3 in k2:
					if k3=='^':
						tk2[-1]+=k3
					else:
						tk2.append(k3)
				k1b2.append(tk2)
		for ka in k1a2:
			for kb in k1b2:
				out.append([ka,kb])
	return out

def struc_produce_two_block(Copy_num_estimate):
	cpNum_a=[j for j in [Copy_num_estimate['a']+i for i in [-1,0,1]] if j>-1]
	cpNum_b=[j for j in [Copy_num_estimate['b']+i for i in [-1,0,1]] if j>-1]
	all_Strucs=two_block_a(cpNum_a,cpNum_b)
	all_Str2=two_block_b(all_Strucs)
	return all_Str2

def two_block_RD_Process(run_flag):
	Letter_Candidates=struc_produce_two_block(Copy_num_estimate)
	if MoveFlag==2:
		Letter_Candidates=Letter_Candidates
	elif MoveFlag==1:
		Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
	elif MoveFlag==0:
		Letter_Candidates=[i for i in Letter_Candidates if ['a','b'] in i]
	P_IL=[]
	P_RD=[]
	P_DR=[]
	P_TB=[]
	Letter_Rec=[]
	BP_Rec=[]
	for Af_Letter in Letter_Candidates:
		Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
		for i in Af_Letter[0]:
				Af_BP[0].append(Af_BP[0][-1]+Be_BP_Letter[i])
		for i in Af_Letter[1]:
				Af_BP[1].append(Af_BP[1][-1]+Be_BP_Letter[i])
		Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
		if not Af_Info_all==0:
			Letter_Rec.append(Af_Letter)
			BP_Rec.append(Af_BP)
			Af_IL_Penal=Af_Info_all[0]
			Af_RD_Rec=Af_Info_all[1]
			Af_DR_Penal=(Af_Info_all[2])**2
			Af_TB_Penal_a=Af_Info_all[4]
			Af_TB_Rec=Af_Info_all[3]
			Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
			Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
			for key in Af_Info_all[5].keys():
					Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
			P_IL.append(Af_IL_Penal)
			P_RD.append(Af_RD_Penal)
			P_DR.append(Af_DR_Penal/num_of_read_pairs)
			P_TB.append(Af_TB_Penal)
	Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
	Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
	Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
	Regu_RD=[i-RD_GS for i in Regu_RD]
	#Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
	Regulator=1
	ILTemp=[j/Regulator for j in Regu_IL]
	RDTemp=[i for i in Regu_RD]
	if not ILTemp==[]:
		DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
		Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
		Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
		run_flag+=1
	else:
		Best_Letter_Rec=[]
		Best_Score_Rec=100
		run_flag+=1
	return([Best_Letter_Rec,Best_Score_Rec,run_flag])

def qual_check_bps2(bps2):
	flag=0
	if len(bps2)==1 and len(bps2[0])==3:
		for x1 in bps2:
			for x2 in range(len(x1)-2):
				if int(x1[x2+2])-int(x1[x2+1])<100: 
					flag=1
	if flag==0:
		return 'right'
	else:
		return 'error'

def let_reclust(vec_in):
	if vec_in==[]:
		return []
	else:
		k2e=[]
		k2e=[vec_in[0]]
		for k3 in range(len(vec_in)-1):
			if '^' in vec_in[k3+1]:
				if '^' in vec_in[k3] and ord(vec_in[k3][0])-ord(vec_in[k3+1][0])==1:
					k2e[-1]+=vec_in[k3+1]
				else:
					k2e.append(vec_in[k3+1])
			else:
				if ord(vec_in[k3+1][0])-ord(vec_in[k3][0])==1 and not '^' in vec_in[k3]:
					k2e[-1]+=vec_in[k3+1]
				else:
					k2e.append(vec_in[k3+1])
		k2f=[]
		for k3 in k2e:
			if '^' in k3:
				k5=''
				for k4 in range(len(k3)/2):
					k5+=k3[2*k4]
				k6=k5[::-1]+'^'
				if not k6 in k2f:
					k2f.append(k6)
			else:
				k2f.append(k3)
		return k2f 

def simple_flag_SA(k1,k2):
	#print [k1,k2]
	if k2=='':
		return[[i for i in k1],[],[],0]
	else:
		temp=[]
		break_flag=0
		for i in k2:
			if not i=='^':
				temp.append(i)
			else:
				temp[-1]+=i
		temp2=[temp[0]]
		for i in range(len(temp[1:])):
			if not '^' in temp[i] and not '^' in temp[i+1] and ord(temp[i+1])-ord(temp[i])==1:
				temp2[-1]+=temp[i+1]
			elif '^' in temp[i] and '^' in temp[i+1] and ord(temp[i+1][0])-ord(temp[i][0])==-1:
				temp2[-1]=temp[i+1][0]+temp2[-1]
			else:
				temp2.append(temp[i+1]) 
		outdel=[]
		outinv=[]
		outdup=[]
		outtra=0
		for i in range(len(temp2)):
			j=temp2[i]
			if '^' in j:
				if not j.replace('^','') in outinv:
					outinv.append(j.replace('^',''))
				temp2[i]=j.replace('^','')
		temp3=''.join(temp2)
		for i in range(len(temp3)-1):
			if not temp3[i+1]=='/' and not temp3[i]=='/':
				if ord(temp3[i+1])-ord(temp3[i])<0:
					outtra=1
		if not temp3==k1:
			temp4=[]
			for i in temp3:
				if temp3.count(i)>1:
					if not i in outdup:
						outdup.append(i)
				if not i in temp4:
					temp4.append(i)
			if not ''.join(temp4)==k1:
				for i in k1:
					if not i in temp4:
						outdel.append(i)
		if not outdup==[]:
			dupuni=unit_produce(outdup)
			outdup2=[]
			k3=k2
			for i in dupuni:
				ia=i
				ib=''.join([j+'^' for j in i[::-1]])
				if len(i)>1:
					if temp2.count(ia)+temp2.count(ib)>2:
						outdup2.append([i,temp2.count(ia)+temp2.count(ib)])
						k3=k3.replace(ia,'')
						k3=k3.replace(ib,'')
				elif len(i)==1:
					if k3.count(ia)+k3.count(ib)>1:
						outdup2.append([i,k3.count(ia)])
						k3=k3.replace(ia,'')
						k3=k3.replace(ib,'')
		else:
			outdup2=[]
		outdup3=[]
		for i in outdup2:
			flag=0
			for j in outdup2:
				if not j==i:
					if len(j[0])>len(i[0]) and i[0] in j[0]:
						flag+=1
			if flag==0:
				outdup3.append(i)		
		if len(outdup3)>1:
			outdup3.sort()
			outdup4=[outdup3[0]]
			for i in range(len(outdup3)-1):
				if outdup3[i+1][-1]==outdup3[i][-1] and ord(outdup3[i+1][0][0])-ord(outdup3[i][0][-1])==1:
					outdup4[-1][0]+=outdup3[i+1][0]
				else:
					outdup4.append(outdup3[i+1])
		else:
			outdup4=outdup3
		return [outdel,outinv,outdup4,outtra]

def unit_produce(list):
	temp1=[sorted(list)[0]]
	for k1 in sorted(list)[1:]:
		if ord(k1)-ord(temp1[-1][-1])==1:
			temp1[-1]+=k1
		else:
			temp1.append(k1)
	temp2=[]
	for k1 in temp1:
		for k2 in range(len(k1)+1)[1:]:
			for k3 in range(len(k1)-k2+1):
				temp2.append(k1[k3:(k3+k2)])
	return temp2[::-1]

def simpler_struc_pick(list_of_structures,original_letters):
	#eg of list_of_structures:[[['b^', 'a^', 'b'], ['a', 'b']],[['b^', 'a', 'b'], ['a', 'b']]]
	for x in list_of_structures:
		x.sort()
	out=[]
	for x in list_of_structures:
		out.append(0)
		for y in x:
			temp=simple_flag_SA(''.join(original_letters),''.join(y))
			temp2={}
			for x in temp[0]:
				if not x in temp2.keys():
					temp2[x]=1
				else:
					temp2[x]+=1
			for x in temp[1]:
				if not x in temp2.keys():
					temp2[x]=1
				else:
					temp2[x]+=1
			for x in temp[2]:
				if not x[0] in temp2.keys():
					temp2[x[0]]=1
				else:
					temp2[x[0]]+=1
			for x in temp2.keys():
				out[-1]+=float(temp2[x]-1)/float(2)+1
	out2=[list_of_structures[i] for i in range(len(out)) if out[i]==min(out)]
	return out2

def Best_Let_modify(Best_Letter_Rec,Best_Score_Rec,Score_rec_hash):
	if not Score_rec_hash=={}:
		temp1=[Best_Score_Rec]
		for x in sorted(Score_rec_hash.keys())[::-1][:2]:
			if not x in temp1 and x-Best_Score_Rec> -1:
				temp1.append(x)
		temp2=[]
		for x in temp1:
			temp2+=Score_rec_hash[x]
		temp3=simpler_struc_pick(temp2,original_letters)
		if len(temp3)==1:
			Best_Letter_Rec=temp3
			for x in temp1:
				if temp3[0] in Score_rec_hash[x]:
					new_best_score=x
			return [Best_Letter_Rec,new_best_score]
		else:
			thash=[]
			for x in temp3:
				for y in temp1:
					if x in Score_rec_hash[y]:
						thash.append(y)
			thash2=[temp3[i] for i in range(len(temp3)) if thash[i]==min(thash)]
			Best_Letter_Rec=thash2
			new_best_score=min(thash)
			return [Best_Letter_Rec,new_best_score] 
	else:
		return [Best_Letter_Rec,Best_Score_Rec] 

def ori_let_Modi(Be_Info, ori_let2):
	Ab_Dir_L1_a=[]
	Ab_Dir_L1_b=[]
	Ab_Dir_L2_a=[]
	Ab_Dir_L2_b=[]
	norm_dir={}
	for y in Be_Info[0]:
		if y[-2:]==['+','+']:
			Ab_Dir_L1_a.append(y)
		elif y[-2:]==['-','-']:
			Ab_Dir_L2_a.append(y)
		#else:
		#	if not y[0] in norm_dir.keys():
		#		norm_dir[y[0]]=1
		#	else:
		#		norm_dir[y[0]]+=1
		#	if not y[3] in norm_dir.keys():
		#		norm_dir[y[3]]=1
		#	else:
		#		norm_dir[y[3]]+=1
	for y in Be_Info[1]:
		if y[-2:]==['+','+']:
			Ab_Dir_L1_b.append(y)
		elif y[-2:]==['-','-']:
			Ab_Dir_L2_b.append(y)
		#else:
		#	if not y[0] in norm_dir.keys():
		#		norm_dir[y[0]]=1
		#	else:
		#		norm_dir[y[0]]+=1
		#	if not y[2] in norm_dir.keys():
		#		norm_dir[y[2]]=1
		#	else:
		#		norm_dir[y[2]]+=1
		#	if not y[4] in norm_dir.keys():
		#		norm_dir[y[4]]=1
		#	else:
		#		norm_dir[y[4]]+=1
		#	if not y[6] in norm_dir.keys():
		#		norm_dir[y[6]]=1
		#	else:
		#		norm_dir[y[6]]+=1
	vote_inv_hash={}
	for x in Ab_Dir_L1_a:
		if not x[3] in vote_inv_hash.keys():
			vote_inv_hash[x[3]]=1
		else:
			vote_inv_hash[x[3]]+=1
	for x in Ab_Dir_L1_b:
		if not x[4] in vote_inv_hash.keys():
			vote_inv_hash[x[4]]=0.5
		else:
			vote_inv_hash[x[4]]+=0.5
		if not x[6] in vote_inv_hash.keys():
			vote_inv_hash[x[6]]=0.5
		else:
			vote_inv_hash[x[6]]+=0.5
	for x in Ab_Dir_L2_a:
		if not x[0] in vote_inv_hash.keys():
			vote_inv_hash[x[0]]=1
		else:
			vote_inv_hash[x[0]]+=1
	for x in Ab_Dir_L2_b:
		if not x[0] in vote_inv_hash.keys():
			vote_inv_hash[x[0]]=0.5
		else:
			vote_inv_hash[x[0]]+=0.5
		if not x[2] in vote_inv_hash.keys():
			vote_inv_hash[x[2]]=0.5
		else:
			vote_inv_hash[x[2]]+=0.5
	new_let2=[]
	for x in ori_let2:
		if x in vote_inv_hash.keys() and vote_inv_hash[x]>2:
			new_let2.append(x+'^')
		else:
			new_let2.append(x)
	new_group=[]
	for x in new_let2:
		if new_group==[]:
			new_group.append([x])
		else:
			if x.count('^')==new_group[-1][-1].count('^'):
				new_group[-1].append(x)
			else:
				new_group.append([x])
	out=[]
	for x in new_group:
		if '^' in x[0] and len(x)>1:
			out+=x[::-1]
		else:
			out+=x
	return out

Penalty_For_InsertLengthZero=-20 #Toy example,decides later	
opts,args=getopt.getopt(sys.argv[1:],'o:f:n:s:',['NullModel=','NullGenomeName=','ReadLen=','NIteration=','MoveFlag=','ppre=','sample=','ref=','QCSplit=','QCAlign='])
dict_opts=dict(opts)
if not '/' in dict_opts['-f']:
	dict_opts['-f']='./'+dict_opts['-f']

if not '--NullModel' in dict_opts.keys():
	model_comp='S'
else:
	if dict_opts['--NullModel'] in ['S','Simple']:
		model_comp='S'
	else:
		model_comp='C'

if '--MoveFlag' in dict_opts.keys():
	MoveFlag=int(dict_opts['--MoveFlag'])
else:
	MoveFlag=2

if '--QCAlign' in dict_opts.keys():
	QCAlign=int(dict_opts['--QCAlign'])
else:
	QCAlign=20

if '--NullGenomeName' in dict_opts.keys():
	genome_name=dict_opts['--NullGenomeName']
else:
	genome_name='genome'

if '--NIteration' in dict_opts.keys():
	Trail_Number=int(dict_opts['--NIteration'])
else:
	Trail_Number=100000

if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
	print 'SVelter-0.1		  Last Update:2014-10-27'
	print 'Required Parameters:'
	print '--ppre, workind directory of SVelter, eg: .../SVelter/' 
	print '--Ref, absolute path of Reference genome. eg: .../SVelter/Reference/genome.fa'
	print '-s, absolute path of input bam file. eg: .../SVelter/BamFiles/Sample.bam'
	print '-f, input txt file containing clustered bps.'
	print ' '
	print 'Optional Parameters:'
	print '--QCAlign, minimum alignment quality required for mapped reads in bam file; default: 20'
	print '--QCSplit, minimum alighment of clipped parts of reads considered as a soft clip; default: 20'
else:
	if not '--ppre' in dict_opts.keys():
		print 'Error: please specify working directory using: --ppre'
	else:
		if not dict_opts['--ppre'][-1]=='/':
			dict_opts['--ppre']+='/'
		workdir=dict_opts['--ppre']
		if not '-f' in dict_opts.keys():
			print 'Error: please specify input txt file using : -f'
		else:
			if not '-o' in dict_opts.keys():
				dict_opts['-o']='/'.join(dict_opts['-f'].split('/')[:-1])
			if not dict_opts['-o'][-1]=='/':
				dict_opts['-o']+='/'
			if not '--ref' in dict_opts.keys():
				print 'Error: please specify refrence genome using --ref'
			else:
				ref_ppre='/'.join(dict_opts['--ref'].split('/')[:-1])+'/'
				ref_prefix='.'.join(dict_opts['--ref'].split('.')[:-1])
				if not '-s' in dict_opts.keys():
					print 'Error: please specify either input file using -s'
				else:
					BamN=dict_opts['-s'].split('/')[-1].replace('.bam','')
					Input_File=dict_opts['-f']
					IL_Weight=1
					RD_Weight=1
					DR_Weight=5
					TB_Weight=5
					#correct_letters_Info=[dict_opts['--CLI']]
					Insert_Len_Stat=dict_opts['--ppre']+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/IL_Null/ILNull.'+BamN+'.'+genome_name+'.Bimodal'
					if not os.path.isfile(Insert_Len_Stat):
						print 'wrong para for --ppre: please copy the null model under -ppre'
					ReadLenFin=dict_opts['--ppre']+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/All_Stats/'+BamN+'.'+genome_name+'.Stats'
					if not os.path.isfile(ReadLenFin):
						print 'wrong para for --ppre: please copy the null model under -ppre'
					else:
						fin=open(ReadLenFin)
						for line in fin:
							pin=line.strip().split()
						fin.close()
						ReadLength=int(pin[-1].split(':')[-1])
					Initial_Bam_Name=BamN+'.bam'
					Initial_Bam=dict_opts['-s']
					#BamName_Unique=BamN+'.cn2.1000'
					#BamFile_stat='_'.join(BamName_Unique.split('.'))
					flank=cdf_solver_application(Insert_Len_Stat,0.95)
					Cut_Lower=cdf_solver_application(Insert_Len_Stat,0.005)
					Cut_Upper=cdf_solver_application(Insert_Len_Stat,0.995)
					if model_comp=='C':
						IL_Stat_all=IL_Stat(Insert_Len_Stat)
						IL_Statistics=IL_Stat_all[0]
						IL_Normal_Stat=IL_Stat_all[1]
					elif model_comp=='S':
						IL_Stat_all=IL_Stat_Simp(Insert_Len_Stat)
						IL_Statistics=IL_Stat_all[0]
						IL_Normal_Stat=IL_Stat_all[1]
					IL_Estimate=IL_Statistics[0]*IL_Statistics[4]+IL_Statistics[1]*IL_Statistics[5]
					IL_SD=((IL_Statistics[2]*IL_Statistics[4])**2+(IL_Statistics[3]*IL_Statistics[5])**2)**(0.5)
					IL_Penal_Two_End_Limit=min([pdf_calculate(IL_Estimate-3*IL_SD,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
					,pdf_calculate(IL_Estimate+3*IL_SD,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)])
					low_qual_edge=5
					fi=open(Input_File)
					bps_hash={}
					bps_temp=[]
					break_flag=0
					for line in fi:
						pi=line.strip().split()
						if pi==[] or len(pi)<3:
							if bps_temp==[]:
								continue
							else:
								bp_key=0
								for l1 in bps_temp:
									bp_key+=len(l1)
								if not bp_key in bps_hash.keys():
									bps_hash[bp_key]=[]
								bps_hash[bp_key].append(bps_temp)
								bps_temp=[]
								if bp_key<3: 
									#print pi
									pi=pi
						else:
							bps_temp.append(pi)
					fi.close()
					bps_hash_inter={}
					for k1 in bps_hash.keys():
						bps_hash_inter[k1]=[]
						for k2 in bps_hash[k1]:
							if not k2 in bps_hash_inter[k1]:
								bps_hash_inter[k1].append(k2)
					bps_hash=bps_hash_inter
					output_Score_File=dict_opts['-o']+'_'.join(dict_opts['-f'].split('/')[-1].split('.')[:-1]+['flag'+str(MoveFlag)])+'_'+'MP'+str(QCAlign)+'.coverge'
					fo=open(output_Score_File,'w')
					fo.close()
					chromos_all=[]
					refFile=dict_opts['--ref']
					if not os.path.isfile(refFile):
						print 'Error: reference file not well defined'
					else:
						refIndFile=refFile+'.fai'
						fin=open(refIndFile)
						for line in fin:
							pin=line.strip().split()
							chromos_all.append(pin[0])
						fin.close()
					for bpsk1 in sorted(bps_hash.keys()):
						for bps2 in bps_hash[bpsk1]:
							for i in bps2:
								if len(i)<3:
									i.append(str(int(i[-1])+100))
					GC_Stat_Path=dict_opts['--ppre']+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/RD_Stat'
					Affix_GC_Stat='_MP'+str(QCAlign)+'_GC_Coverage_ReadLength'
					GC_Stat=GC_Stat_ReadIn(BamN,GC_Stat_Path,Affix_GC_Stat)
					GC_Content_Coverage=GC_Stat[0]
					Chromosome=GC_Stat[1]
					Coverage=[int(k) for k in GC_Stat[2][1:]]
					GC_Overall_Median_Coverage={}
					GC_Overall_Median_Num=[]
					GC_Median_Coverage={}
					GC_Median_Num={}
					GC_Mean_Coverage={}
					GC_Std_Coverage={}
					GC_Var_Coverage={}
					for a in Chromosome:
						GC_Overall_temp=[]
						GC_Median_Coverage[a]={}
						for b in Coverage:
							if not b in GC_Median_Num.keys():
								GC_Median_Num[b]=[]
							if len(GC_Content_Coverage[a][b][0])==2: continue
							elif len(GC_Content_Coverage[a][b][0])>2:
									num_list=[float(c) for c in GC_Content_Coverage[a][b][0][2:].split(',')]
									GC_Median_Num[b]+=num_list
									GC_Overall_Median_Num+=num_list
									GC_Overall_temp=GC_Overall_temp+num_list
									if not Median_Pick(num_list)==0.0:
											GC_Median_Coverage[a][b]=Median_Pick(num_list)
						if len(GC_Overall_temp)==0: continue
						elif len(GC_Overall_temp)>0: 
								GC_Overall_Median_Coverage[a]=Median_Pick(GC_Overall_temp)
								GC_Mean_Coverage[a]=numpy.mean(GC_Overall_temp)
								GC_Std_Coverage[a]=numpy.std(GC_Overall_temp)
								GC_Var_Coverage[a]=(GC_Std_Coverage[a])**2
					GC_Overall_Median_Num=Median_Pick(GC_Overall_Median_Num)
					for a in GC_Median_Num.keys():
						if GC_Median_Num[a]==[]:
							GC_Median_Num[a]=GC_Overall_Median_Num
						else:
							GC_Median_Num[a]=Median_Pick(GC_Median_Num[a])
					ChrN_Median_Coverage={}
					for i in GC_Median_Coverage.keys():
						for j in GC_Median_Coverage[i].keys():
							if not j in ChrN_Median_Coverage.keys():
								ChrN_Median_Coverage[j]=[GC_Median_Coverage[i][j]]
							else:
								ChrN_Median_Coverage[j]+=[GC_Median_Coverage[i][j]]
					refFlag=0
					refgenome=dict_opts['--ref']
					reftest=open(refgenome)
					preftest=reftest.readline().strip().split()
					if len(preftest[0])>3:
						refFlag=1
					reftest.close()
					if refFlag==0:
						chrom_N='N'
						chrom_X='X'
						chrom_Y='Y'
					else:
						chrom_N='chrN'
						chrom_X='chrX'
						chrom_Y='chrY'		
					GC_Median_Coverage[chrom_N]={}
					if not chrom_X in GC_Median_Coverage.keys():
						GC_Median_Coverage[chrom_X]={}
					if not chrom_Y in GC_Median_Coverage.keys():
						GC_Median_Coverage[chrom_Y]={}
					for i in ChrN_Median_Coverage.keys():
						GC_Median_Coverage[chrom_N][i]=numpy.mean(ChrN_Median_Coverage[i])
						GC_Median_Coverage[chrom_Y][i]=numpy.mean(ChrN_Median_Coverage[i])
						GC_Median_Coverage[chrom_X][i]=numpy.mean(ChrN_Median_Coverage[i])
					GC_Overall_Median_Coverage[chrom_N]=numpy.mean([GC_Overall_Median_Coverage[key] for key in GC_Overall_Median_Coverage.keys()])
					GC_Overall_Median_Coverage[chrom_X]=numpy.mean([GC_Overall_Median_Coverage[key] for key in GC_Overall_Median_Coverage.keys()])
					GC_Overall_Median_Coverage[chrom_Y]=numpy.mean([GC_Overall_Median_Coverage[key] for key in GC_Overall_Median_Coverage.keys()])
					GC_Var_Coverage[chrom_N]=numpy.mean([GC_Var_Coverage[key] for key in GC_Var_Coverage.keys()])
					GC_Var_Coverage[chrom_X]=numpy.mean([GC_Var_Coverage[key] for key in GC_Var_Coverage.keys()])
					GC_Var_Coverage[chrom_Y]=numpy.mean([GC_Var_Coverage[key] for key in GC_Var_Coverage.keys()])
					GC_Mean_Coverage[chrom_N]=numpy.mean([GC_Mean_Coverage[key] for key in GC_Mean_Coverage.keys()])
					GC_Mean_Coverage[chrom_X]=numpy.mean([GC_Mean_Coverage[key] for key in GC_Mean_Coverage.keys()])
					GC_Mean_Coverage[chrom_Y]=numpy.mean([GC_Mean_Coverage[key] for key in GC_Mean_Coverage.keys()])
					GC_Std_Coverage[chrom_N]=numpy.mean([GC_Std_Coverage[key] for key in GC_Std_Coverage.keys()])
					GC_Std_Coverage[chrom_X]=numpy.mean([GC_Std_Coverage[key] for key in GC_Std_Coverage.keys()])
					GC_Std_Coverage[chrom_Y]=numpy.mean([GC_Std_Coverage[key] for key in GC_Std_Coverage.keys()])
					for bpsk1 in sorted(bps_hash.keys()):
						for bps2 in bps_hash[bpsk1]:
							if qual_check_bps2(bps2)=='right':
								Chromo=bps2[0][0]
								K_RD=GC_Std_Coverage[str(Chromo)]/GC_Mean_Coverage[str(Chromo)]
								K_IL=IL_Normal_Stat[2]/IL_Normal_Stat[1]
								K_RD_new=1
								K_IL_new=(K_IL/K_RD)**2
								IL_GS=Prob_Norm(IL_Normal_Stat[1],IL_Normal_Stat[1],IL_Normal_Stat[2]**2)
								RD_GS=Prob_Norm(GC_Mean_Coverage[str(Chromo)],GC_Mean_Coverage[str(Chromo)],GC_Std_Coverage[str(Chromo)]**2)
								#RD_GS=Prob_Norm(GC_Mean_Coverage[str(Chromo)]/2,GC_Mean_Coverage[str(Chromo)]/2,GC_Var_Coverage[str(Chromo)]/2)
								#print bps2
								for i in bps2:
									temp2=[int(j) for j in i[1:]]
									k=[i[0]]+sorted(temp2)
									k2=k[:2]
									for k3 in temp2:
										if not k3 in k2 and k3-k2[-1]>10:
											k2.append(k3)
									if len(k2)>2:
										bps2[bps2.index(i)]=k2
									else:
										del bps2[bps2.index(i)]
								original_bps_all=[]
								for obas in bps2:
									original_bps_all+=obas
								original_structure=bp_to_let(original_bps_all,chromos_all)
								chr_letter_tbp=letter_rearrange(bps2)
								letter_tGC=letter_GC_ReadIn(chr_letter_tbp)
								if letter_tGC=='error': continue
								letter_tRD=letter_RD_ReadIn(chr_letter_tbp)
								if letter_tRD=='error': continue
								chr_letter_bp={}
								letter_GC={}
								letter_RD={}
								for k1 in chr_letter_tbp.keys():
									chr_letter_bp[k1]={}
									letter_GC[k1]={}
									letter_RD[k1]={}
									for k2 in chr_letter_tbp[k1].keys():
										if k2 in letter_tGC[k1].keys() and k2 in letter_tRD[k1].keys() and not math.isnan(letter_tRD[k1][k2]) and not math.isnan(letter_tGC[k1][k2]):
											chr_letter_bp[k1][k2]=chr_letter_tbp[k1][k2]
											letter_GC[k1][k2]=letter_tGC[k1][k2]
											letter_RD[k1][k2]=letter_tRD[k1][k2]
								left_keys=[]
								for k1 in chr_letter_bp.keys():
									for k2 in chr_letter_bp[k1].keys():
										left_keys.append(k2)
								if not left_keys==[]:
									bamFlag=0
									bamtest=os.popen(r'''samtools view -H %s'''%(dict_opts['-s']))
									for line in bamtest:
										pbamtest=line.strip().split()
										if pbamtest[0]=='@SQ' and pbamtest[1].split(':')[0]=='SN':
											if len(pbamtest[1].split(':')[1])>3:
												bamFlag+=1
											break
									bps3={}
									for k1 in chr_letter_bp.keys():
										bps3[k1]={}
										for k2 in chr_letter_bp[k1].keys():
											bps3[k1][chr_letter_bp[k1][k2][0]]=[chr_letter_bp[k1][k2][0],chr_letter_bp[k1][k2][-1]]
									bps4={}
									for k1 in bps3.keys():
										if not bps3[k1]=={}:
											bps4[k1]=[[k1]+bps3[k1][sorted(bps3[k1].keys())[0]]]
											for k2 in range(len(bps3[k1].keys())-1):
												if bps3[k1][sorted(bps3[k1].keys())[k2+1]][0]==bps3[k1][sorted(bps3[k1].keys())[k2]][-1]:
													bps4[k1][-1]+=[bps3[k1][sorted(bps3[k1].keys())[k2+1]][-1]]
												else:
													bps4[k1].append([bps3[k1][sorted(bps3[k1].keys())[k2+1]]])
									bps2=[]
									for k1 in bps4.keys():
										for k2 in bps4[k1]:
											bps2.append(k2)
									Chr=bps2[0][0]
									Through_BP_Minimum=GC_Overall_Median_Coverage[str(Chr)]/10
									if len(Chr)>3:
										refChr0=Chr[3:]
										refChr1=Chr
									else:
										refChr0=Chr
										refChr1='chr'+Chr
									if refFlag==0:
										refChr=refChr0
									else:
										refChr=refChr1
									if bamFlag==0:
										bamChr=refChr0
									else:
										bamChr=refChr1
									bamtest.close()
									local_Minimum=[]
									block2_bin=[]
									GC_blocks_index=[]
									Full_Info=Full_Info_of_Reads_Integrate(bps2)
									RD_within_B=Full_Info[0]
									RD_within_B['left']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
									RD_within_B['right']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
									for i in RD_within_B.keys():
										RD_within_B[i+'^']=RD_within_B[i]
									Initial_GCRD_Adj=Full_Info[1]
									Initial_GCRD_Adj['left']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
									Initial_GCRD_Adj['right']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
									Best_IL_Score=0
									for j in range(Cut_Lower,Cut_Upper+1):
										Single_ILScore=pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
										Best_IL_Score+=Single_ILScore*exp(Single_ILScore)
									Best_RD_Score=0
									let_chr_rec={}
									for i in chr_letter_bp.keys():
										for j in chr_letter_bp[i].keys():
											if j in left_keys:
												let_chr_rec[j]=i
									for i in let_chr_rec.keys():
										Theo_RD=GC_Overall_Median_Coverage[str(let_chr_rec[i])]
										Theo_Var=GC_Var_Coverage[str(let_chr_rec[i])]
										for j in range(int(Theo_RD/2),int(Theo_RD/2*3+1)):
											single_ProbNB=Prob_Norm(j,Theo_RD,Theo_Var)
											Best_RD_Score+=single_ProbNB*exp(single_ProbNB)
									Block_CN_Upper={}
									Copy_num_estimate={}
									for i in Initial_GCRD_Adj.keys():
										if not i in ['left','right']:
											Copy_num_estimate[i]=round(Initial_GCRD_Adj[i]*2/GC_Mean_Coverage[Chr])
											if Initial_GCRD_Adj[i]<float(GC_Mean_Coverage[Chr])/10.0:
												Copy_num_estimate[i]=-1
									Copy_num_Check=[]
									for CNE in Copy_num_estimate.keys():
										if Copy_num_estimate[CNE]>5:
											Copy_num_Check.append(CNE)
									if Copy_num_Check==[]:
										median_CN=GC_Overall_Median_Coverage[chrom_N]/2
										for key in Initial_GCRD_Adj.keys():
											if not key in ['left','right']:
												Block_CN_Upper[key]=Initial_GCRD_Adj[key]/median_CN+2
										Initial_DR=Full_Info[2]
										Initial_IL=Full_Info[3]
										#Initial_Info=Full_Info[4]+Full_Info[5]+Full_Info[6]
										BlockGC=Full_Info[7]
										BlockGC['left']=0.476
										BlockGC['right']=0.476
										BlockGC2={}
										for key_B_GC in BlockGC.keys():
											BlockGC2[key_B_GC]=BlockGC[key_B_GC]
											BlockGC2[key_B_GC+'^']=BlockGC[key_B_GC]	
										original_letters=Full_Info[9]
										original_bp_list=Full_Info[8]
										Be_BP_Letter={}
										for let_key in original_letters:
											Be_BP_Letter[let_key]=original_bp_list[original_letters.index(let_key)+1]-original_bp_list[original_letters.index(let_key)]
										ori_let2=[]
										for i in original_letters:
											ori_let2.append(i)
										for i in original_letters:
											if Copy_num_estimate[i]<1:
												ori_let2.remove(i)
											elif Copy_num_estimate[i]>3:
												letter_copy=int(Copy_num_estimate[i]/2)
												for j in range(letter_copy)[1:]:
													ori_let2.append(i)
										ori_bp2=[original_bp_list[0]]
										for i in ori_let2:
											ori_bp2.append(ori_bp2[-1]+Be_BP_Letter[i])
										Initial_TB=0
										Initial_Move_Prob=[1.0/3,1.0/3,1.0/3]
										Pair_Through=Full_Info[4]
										Read_Through=Full_Info[5]
										SingleR_Through=Full_Info[6]
										bp_MP=[original_bp_list,original_bp_list]
										letter_MP=[original_letters,original_letters]
										Be_BP=[ori_bp2,ori_bp2]
										Be_BP_Letter['left']=flank
										Be_BP_Letter['right']=flank
										for let_key in Be_BP_Letter.keys():
											Be_BP_Letter[let_key+'^']=Be_BP_Letter[let_key]
										num_of_read_pairs=1
										for k1 in Be_BP_Letter.keys():
												if not k1[-1]=='^' and not k1 in ['left','right']:
														#print k1
														num_of_read_pairs+=Be_BP_Letter[k1]*RD_within_B[k1]/2/ReadLength
										num_of_read_pairs+=len(Full_Info[4])+len(Full_Info[5])+len(Full_Info[6])
										Be_Info=[Pair_Through,Read_Through,SingleR_Through]
										ori_let2=ori_let_Modi(Be_Info,ori_let2)
										Be_Letter=[ori_let2,ori_let2]
										Move_Step=0
										best_iterations=0
										Best_Score=float("-inf")
										Best_Letter=[]
										Best_BPs=[]
										score_record=[]
										#best_score_rec=[]
										total_read_Pairs_Num=(original_bp_list[-1]-original_bp_list[0])*GC_Mean_Coverage[Chr]/2/ReadLength
										num_of_reads=total_read_Pairs_Num
										Best_Score_Rec=0
										Score_rec_hash={}
										break_Iteration_Flag=0
										run_flag=0
										if len(Full_Info[9])==1:
											if Full_Info[1]['a']<GC_Mean_Coverage[Chr]/4 and Full_Info[2]<3:
												Run_Result=zero_RD_Process(run_flag)
												Best_Letter_Rec=Run_Result[0]
												Best_Score_Rec=Run_Result[1]
												run_flag=Run_Result[2]
												Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
											else:
												if Full_Info[1]['a']<GC_Mean_Coverage[Chr]:
													Run_Result=one_RD_Process(run_flag,Score_rec_hash)
													Best_Letter_Rec=Run_Result[0]
													Best_Score_Rec=Run_Result[1]
													run_flag=Run_Result[2]
													Score_rec_hash=Run_Result[3]
													#Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
												else:
													if Full_Info[1]['a']<2*GC_Mean_Coverage[Chr]:
														Run_Result=two_RD_Process(run_flag,Score_rec_hash)
														Best_Letter_Rec=Run_Result[0]
														Best_Score_Rec=Run_Result[1]
														run_flag=Run_Result[2]
														Score_rec_hash=Run_Result[3]
														#Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
													else:
														copy_num_a=int(float(Full_Info[1]['a'])/(float(GC_Mean_Coverage[Chr])/2))
														copy_num_b=int(float(Full_Info[1]['a'])/(float(GC_Mean_Coverage[Chr])/2))+1
														if copy_num_b<6:
															Run_Result=few_RD_Process(run_flag,Score_rec_hash)
															Best_Letter_Rec=Run_Result[0]
															Best_Score_Rec=Run_Result[1]
															run_flag=Run_Result[2]
															Score_rec_hash=Run_Result[3]
															#Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
														else:
															Run_Result=many_RD_Process(run_flag)
															Best_Letter_Rec=Run_Result[0]
															Best_Score_Rec=Run_Result[1]
															run_flag=Run_Result[2]
															Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
										elif len(Full_Info[9])==2:
											bl2_flag=0
											for keyCNE in Copy_num_estimate.keys():
												if not Copy_num_estimate[keyCNE]<2:
													bl2_flag+=1
											if bl2_flag==0:
												Run_Result=two_block_RD_Process(run_flag)
												Best_Letter_Rec=Run_Result[0]
												Best_Score_Rec=Run_Result[1]
												run_flag=Run_Result[2]
												Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
										if run_flag==0:
											speed_test=10
											t1_sptest=time.time()
											while True:
												if Move_Step>speed_test: break
												Move_Step+=1
												M_Key_0='_'.join(['Step',str(Move_Step),'M'])
												P_Key_0='_'.join(['Step',str(Move_Step),'P'])
												Move_Sample_Pool=['delete','invert','insert']
												Move_M_P=Move_Choose(Move_Sample_Pool,MoveFlag,Initial_Move_Prob)
												#print Best_Letter
												#if not Move_M_P[1]=='x':
												#print Be_Letter+[Move_Step]
												M_Move_Choices=Move_Choice_procedure_2(Move_M_P[0],Be_Letter[0],original_letters,'2m')
												P_Move_Choices=Move_Choice_procedure_2(Move_M_P[1],Be_Letter[1],original_letters,'2p')
												if M_Move_Choices=='ERROR!' or P_Move_Choices=='ERROR!':
														Move_Step-=1
														continue
												if not M_Move_Choices=='ERROR!' and not P_Move_Choices=='ERROR!':
														#print Be_Letter
														P_IL=[]
														P_RD=[]
														P_DR=[]
														P_TB=[]
														Letter_Rec=[]
														BP_Rec=[]
														for m in [['2m','1','1','1','X']]+M_Move_Choices:
																p=[str(Chr)+'p','1','1','1','X']
																Move_MP=[m,p]
																Af_BP=[BPList_Rearrange(Be_BP[0],m,original_bp_list),BPList_Rearrange(Be_BP[1],p,original_bp_list)]
																Af_Letter=[LetterList_Rearrange(Be_Letter[0],m,original_bp_list),LetterList_Rearrange(Be_Letter[1],p,original_bp_list)]
																if MoveFlag==1:
																	Af_Letter[1]=Af_Letter[0]
																	Af_BP[1]=Af_BP[0]
																if not Af_Letter_QC(Af_Letter,Copy_num_estimate)==0:continue
																if not Best_Score_Rec==0 and Af_Letter in Best_Letter_Rec: continue
																letter_num_flag=0
																for key in Block_CN_Upper.keys():
																		if (Af_Letter[0]+Af_Letter[1]).count(key)>Block_CN_Upper[key]:
																				letter_num_flag+=1
																if not letter_num_flag==0: continue
																Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
																if Af_Info_all==0:continue
																Letter_Rec.append(Af_Letter)
																BP_Rec.append(Af_BP)
																Af_IL_Penal=Af_Info_all[0]
																Af_RD_Rec=Af_Info_all[1]
																Af_DR_Penal=(Af_Info_all[2])**2
																Af_TB_Penal_a=Af_Info_all[4]
																Af_TB_Rec=Af_Info_all[3]
																Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
																Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
																for key in Af_Info_all[5].keys():
																		Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
																P_IL.append(Af_IL_Penal)
																P_RD.append(Af_RD_Penal)
																P_DR.append(Af_DR_Penal/num_of_read_pairs)
																P_TB.append(Af_TB_Penal)
														if len(P_IL)==0: continue
														Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
														Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
														Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
														Regu_RD=[i-RD_GS for i in Regu_RD]
														Regulator=1
														ILTemp=[j/Regulator for j in Regu_IL]
														RDTemp=[i for i in Regu_RD]
														DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
														DECISION=DECISION_Score[0]
														S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
														Be_Letter=Letter_Rec[DECISION]
														Be_BP=BP_Rec[DECISION]
														if not S_DECISION in Score_rec_hash.keys():
															Score_rec_hash[S_DECISION]=[]
														Score_rec_hash[S_DECISION].append(Be_Letter)
														if S_DECISION>Best_Score:
																Best_Letter=[Be_Letter]
																Best_BPs=[Be_BP]
																Best_Score=S_DECISION
																best_iterations=0
														elif S_DECISION==Best_Score:
																if not Be_Letter in Best_Letter:
																		Best_Letter+=[Be_Letter]
																		Best_BPs+=[Be_BP]
																best_iterations+=1
														else:
																best_iterations+=1
														score_record.append(S_DECISION)
														#best_score_rec.append(Best_Score)
														P_IL=[]
														P_RD=[]
														P_DR=[]
														P_TB=[]
														Letter_Rec=[]
														BP_Rec=[]
														if not P_Move_Choices==[]:
															for p in [['2p','1','1','1','X']]+P_Move_Choices:
																m=[str(Chr)+'m','1','1','1','X']
																Move_MP=[m,p]
																Af_BP=[BPList_Rearrange(Be_BP[0],m,original_bp_list),BPList_Rearrange(Be_BP[1],p,original_bp_list)]
																Af_Letter=[LetterList_Rearrange(Be_Letter[0],m,original_bp_list),LetterList_Rearrange(Be_Letter[1],p,original_bp_list)]
																if MoveFlag==1:
																	Af_Letter[1]=Af_Letter[0]
																	Af_BP[1]=Af_BP[0]
																if not Af_Letter_QC(Af_Letter,Copy_num_estimate)==0:continue
																if not Best_Score_Rec==0 and Af_Letter in Best_Letter_Rec: continue
																letter_num_flag=0
																for key in Block_CN_Upper.keys():
																		if (Af_Letter[0]+Af_Letter[1]).count(key)>Block_CN_Upper[key]:
																				letter_num_flag+=1
																if not letter_num_flag==0: continue
																Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
																if Af_Info_all==0:
																		continue
																Letter_Rec.append(Af_Letter)
																BP_Rec.append(Af_BP)
																Af_IL_Penal=Af_Info_all[0]
																Af_RD_Rec=Af_Info_all[1]
																Af_DR_Penal=(Af_Info_all[2])**2
																Af_TB_Penal_a=Af_Info_all[4]
																Af_TB_Rec=Af_Info_all[3]
																Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
																Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
																for key in Af_Info_all[5].keys():
																		Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N])
																P_IL.append(Af_IL_Penal)
																P_RD.append(Af_RD_Penal)
																P_DR.append(Af_DR_Penal/num_of_read_pairs)
																P_TB.append(Af_TB_Penal)
														if len(P_IL)==0: continue
														Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
														Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
														Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
														Regu_RD=[i-RD_GS for i in Regu_RD]
														Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
														Regulator=1
														ILTemp=[j/Regulator for j in Regu_IL]
														RDTemp=[i for i in Regu_RD]
														DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
														DECISION=DECISION_Score[0]
														S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
														Be_Letter=Letter_Rec[DECISION]
														Be_BP=BP_Rec[DECISION]
														if not S_DECISION in Score_rec_hash.keys():
															Score_rec_hash[S_DECISION]=[]
														Score_rec_hash[S_DECISION].append(Be_Letter)
														if S_DECISION>Best_Score:
																Best_Letter=[Be_Letter]
																Best_BPs=[Be_BP]
																Best_Score=S_DECISION
																best_iterations=0
														elif S_DECISION==Best_Score:
																if not Be_Letter in Best_Letter:
																		Best_Letter+=[Be_Letter]
																		Best_BPs+=[Be_BP]
																best_iterations+=1
														else:
																best_iterations+=1
														score_record.append(S_DECISION)
														#best_score_rec.append(Best_Score)
											t2_sptest=time.time()
											if t2_sptest-t1_sptest<30 or bpsk1<4:
												while True:
														if Move_Step>Trail_Number: break
														if best_iterations>100: 
																if Best_Score_Rec==0:
																		best_iterations=0
																		Best_Score_Rec=Best_Score
																		Best_Letter_Rec=Best_Letter
																		Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
																		Best_BPs_Rec=Best_BPs
																		Be_Letter=Best_Letter[0]
																		Be_BP=Best_BPs[0]
																		Best_Score-=100
																		#print 'current best structure:'+str(Best_Score_Rec)
																else:
																		if Best_Score<Best_Score_Rec:
																				break_Iteration_Flag=1
																				#print 'final best structure:'+str(Best_Score_Rec)
																		elif Best_Score==Best_Score_Rec:
																				break_Iteration_Flag=1
																				for i in Best_Letter:
																						if not i in Best_Letter_Rec:
																								Best_Letter_Rec.append(i)
																		else:
																				best_iterations=0
																				Best_Score_Rec=Best_Score
																				Best_Letter_Rec=Best_Letter
																				Best_BPs_Rec=Best_BPs
																				Be_Letter=Best_Letter[0]
																				Be_BP=Best_BPs[0]
																				Best_Score-=100
																				#print 'current best structure:'+str(Best_Score_Rec)
														if break_Iteration_Flag>0:
																break
														Move_Step+=1
														M_Key_0='_'.join(['Step',str(Move_Step),'M'])
														P_Key_0='_'.join(['Step',str(Move_Step),'P'])
														Move_Sample_Pool=['delete','invert','insert']
														Move_M_P=Move_Choose(Move_Sample_Pool,MoveFlag,Initial_Move_Prob)
														M_Move_Choices=Move_Choice_procedure_2(Move_M_P[0],Be_Letter[0],original_letters,'2m')
														P_Move_Choices=Move_Choice_procedure_2(Move_M_P[1],Be_Letter[1],original_letters,'2p')
														if M_Move_Choices=='ERROR!' or P_Move_Choices=='ERROR!':
																Move_Step-=1
																continue
														if not M_Move_Choices=='ERROR!' and not P_Move_Choices=='ERROR!':
																#print Be_Letter
																P_IL=[]
																P_RD=[]
																P_DR=[]
																P_TB=[]
																Letter_Rec=[]
																BP_Rec=[]
																for m in [['2m','1','1','1','X']]+M_Move_Choices:
																		p=[str(Chr)+'p','1','1','1','X']
																		Move_MP=[m,p]
																		Af_BP=[BPList_Rearrange(Be_BP[0],m,original_bp_list),BPList_Rearrange(Be_BP[1],p,original_bp_list)]
																		Af_Letter=[LetterList_Rearrange(Be_Letter[0],m,original_bp_list),LetterList_Rearrange(Be_Letter[1],p,original_bp_list)]
																		if MoveFlag==1:
																			Af_Letter[1]=Af_Letter[0]
																			Af_BP[1]=Af_BP[0]
																		if not Af_Letter_QC(Af_Letter,Copy_num_estimate)==0:continue
																		if not Best_Score_Rec==0 and Af_Letter in Best_Letter_Rec: continue
																		letter_num_flag=0
																		for key in Block_CN_Upper.keys():
																				if (Af_Letter[0]+Af_Letter[1]).count(key)>Block_CN_Upper[key]:
																						letter_num_flag+=1
																		if not letter_num_flag==0: continue
																		Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
																		if Af_Info_all==0:continue
																		Letter_Rec.append(Af_Letter)
																		BP_Rec.append(Af_BP)
																		Af_IL_Penal=Af_Info_all[0]
																		Af_RD_Rec=Af_Info_all[1]
																		Af_DR_Penal=(Af_Info_all[2])**2
																		Af_TB_Penal_a=Af_Info_all[4]
																		Af_TB_Rec=Af_Info_all[3]
																		Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
																		Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
																		for key in Af_Info_all[5].keys():
																				Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
																		P_IL.append(Af_IL_Penal)
																		P_RD.append(Af_RD_Penal)
																		P_DR.append(Af_DR_Penal/num_of_read_pairs)
																		P_TB.append(Af_TB_Penal)
																if len(P_IL)==0: continue
																Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
																Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
																Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
																Regu_RD=[i-RD_GS for i in Regu_RD]
																Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
																Regulator=1
																ILTemp=[j/Regulator for j in Regu_IL]
																RDTemp=[i for i in Regu_RD]
																DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
																DECISION=DECISION_Score[0]
																S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
																Be_Letter=Letter_Rec[DECISION]
																Be_BP=BP_Rec[DECISION]
																if not S_DECISION in Score_rec_hash.keys():
																	Score_rec_hash[S_DECISION]=[Be_Letter]
																else: 
																	if not Be_Letter in Score_rec_hash[S_DECISION]:
																		Score_rec_hash[S_DECISION].append(Be_Letter)
																if S_DECISION>Best_Score:
																		Best_Letter=[Be_Letter]
																		Best_BPs=[Be_BP]
																		Best_Score=S_DECISION
																		best_iterations=0
																elif S_DECISION==Best_Score:
																		if not Be_Letter in Best_Letter:
																				Best_Letter+=[Be_Letter]
																				Best_BPs+=[Be_BP]
																		best_iterations+=1
																else:
																		best_iterations+=1
																score_record.append(S_DECISION)
																#best_score_rec.append(Best_Score)
																P_IL=[]
																P_RD=[]
																P_DR=[]
																P_TB=[]
																Letter_Rec=[]
																BP_Rec=[]
																if not P_Move_Choices==[]:
																	for p in [['2p','1','1','1','X']]+P_Move_Choices:
																		m=[str(Chr)+'m','1','1','1','X']
																		Move_MP=[m,p]
																		Af_BP=[BPList_Rearrange(Be_BP[0],m,original_bp_list),BPList_Rearrange(Be_BP[1],p,original_bp_list)]
																		Af_Letter=[LetterList_Rearrange(Be_Letter[0],m,original_bp_list),LetterList_Rearrange(Be_Letter[1],p,original_bp_list)]
																		if MoveFlag==1:
																			Af_Letter[1]=Af_Letter[0]
																			Af_BP[1]=Af_BP[0]
																		if not Af_Letter_QC(Af_Letter,Copy_num_estimate)==0:continue
																		if not Best_Score_Rec==0 and Af_Letter in Best_Letter_Rec: continue
																		letter_num_flag=0
																		for key in Block_CN_Upper.keys():
																				if (Af_Letter[0]+Af_Letter[1]).count(key)>Block_CN_Upper[key]:
																						letter_num_flag+=1
																		if not letter_num_flag==0: continue
																		Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
																		if Af_Info_all==0:
																				continue
																		Letter_Rec.append(Af_Letter)
																		BP_Rec.append(Af_BP)
																		Af_IL_Penal=Af_Info_all[0]
																		Af_RD_Rec=Af_Info_all[1]
																		Af_DR_Penal=(Af_Info_all[2])**2
																		Af_TB_Penal_a=Af_Info_all[4]
																		Af_TB_Rec=Af_Info_all[3]
																		Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
																		Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
																		for key in Af_Info_all[5].keys():
																				Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N])
																		P_IL.append(Af_IL_Penal)
																		P_RD.append(Af_RD_Penal)
																		P_DR.append(Af_DR_Penal/num_of_read_pairs)
																		P_TB.append(Af_TB_Penal)
																if len(P_IL)==0: continue
																Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
																Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
																Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
																Regu_RD=[i-RD_GS for i in Regu_RD]
																Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
																Regulator=1
																ILTemp=[j/Regulator for j in Regu_IL]
																RDTemp=[i for i in Regu_RD]
																DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
																DECISION=DECISION_Score[0]
																S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
																Be_Letter=Letter_Rec[DECISION]
																Be_BP=BP_Rec[DECISION]
																if not S_DECISION in Score_rec_hash.keys():
																	Score_rec_hash[S_DECISION]=[Be_Letter]
																else: 
																	if not Be_Letter in Score_rec_hash[S_DECISION]:
																		Score_rec_hash[S_DECISION].append(Be_Letter)
																if S_DECISION>Best_Score:
																		Best_Letter=[Be_Letter]
																		Best_BPs=[Be_BP]
																		Best_Score=S_DECISION
																		best_iterations=0
																elif S_DECISION==Best_Score:
																		if not Be_Letter in Best_Letter:
																				Best_Letter+=[Be_Letter]
																				Best_BPs+=[Be_BP]
																		best_iterations+=1
																else:
																		best_iterations+=1
																score_record.append(S_DECISION)
																#best_score_rec.append(Best_Score)
											else:
												gaps=[]
												bps2_new=[]
												for k1 in bps2:
													gaps.append([])
													for k2 in range(len(k1)-2):
														gaps[-1].append(int(k1[k2+2])-int(k1[k2+1]))
												for k1 in range(len(gaps)):
													bps2_new.append([])
													chr_rec=bps2[k1][0]
													rec1=1
													for k2 in range(len(gaps[k1])):
														if gaps[k1][k2]==max(gaps[k1]):
															bps2_new[-1].append([chr_rec]+bps2[k1][rec1:(k2+2)])
															bps2_new[-1].append([chr_rec]+bps2[k1][(k2+1):(k2+3)])
															rec1=k2+2
													bps2_new[-1].append([chr_rec]+bps2[k1][rec1:])
												for k1 in bps2_new:
													for k2 in k1:
														bps_hash[max(bps_hash.keys())].append([k2])
												Best_Letter_Rec=[]
												Best_Score_Rec=100
										struc_to_remove=[]
										for bestletter in Best_Letter_Rec:
											if '/'.join([''.join(bestletter[0]),''.join(bestletter[1])])==original_structure:
												struc_to_remove.append(bestletter)
										Best_Letter_Rec=[i for i in Best_Letter_Rec if not i in struc_to_remove]
										if Best_Letter_Rec==[] and Best_Score_Rec==100:
											continue
										else:
											write_best_letter(bps2,Best_Letter_Rec,Best_Score_Rec,Score_rec_hash,original_letters)
									else:
										Score_rec_hash={}
										bps_new={}
										Initial_DR=Full_Info[2]
										Initial_IL=Full_Info[3]
										#Initial_Info=Full_Info[4]+Full_Info[5]+Full_Info[6]
										BlockGC=Full_Info[7]
										BlockGC['left']=0.476
										BlockGC['right']=0.476
										BlockGC2={}
										for key_B_GC in BlockGC.keys():
											BlockGC2[key_B_GC]=BlockGC[key_B_GC]
											BlockGC2[key_B_GC+'^']=BlockGC[key_B_GC]	
										original_letters=Full_Info[9]
										original_bp_list=Full_Info[8]
										for bl in Copy_num_Check:
											for blk1 in chr_letter_bp.keys():
												for blk2 in sorted(chr_letter_bp[blk1].keys()):
													if blk2==bl:
														bps2_temp=[blk1]+[chr_letter_bp[blk1][blk2][0],chr_letter_bp[blk1][blk2][-1]]
														copy_num_a=int(Copy_num_estimate[bl])/2
														copy_num_b=Copy_num_estimate[bl]-copy_num_a
														Best_Letter_Rec=[[['a' for i in range(copy_num_a)],['a' for i in range(copy_num_a)]]]
														Best_Score_Rec=100
														write_best_letter([bps2_temp],Best_Letter_Rec,Best_Score_Rec,Score_rec_hash,original_letters)
										for blk1 in chr_letter_bp.keys():
											bps_new[blk1]=[]
											for blk2 in sorted(chr_letter_bp[blk1].keys()):
												if not blk2 in Copy_num_Check:
													bps_new[blk1].append([chr_letter_bp[blk1][blk2][0],chr_letter_bp[blk1][blk2][-1]])
										bps_new_2=[]
										for k1 in bps_new.keys():
											for k2 in bps_new[k1]:
												if bps_new_2==[]:
													bps_new_2.append([k1]+k2)
												else:
													if k1==bps_new_2[-1][0] and k2[0]==bps_new_2[-1][-1]:
														bps_new_2[-1]+=k2[1:]
													else:
														bps_new_2.append([k1]+k2)
										for k1 in bps_new_2:
											bps_hash[max(bps_hash.keys())].append([k1])


#functions not in use
def Block_GC_Production(ori_1_Seq, bps, flank):
	letters=[chr(97+i) for i in range(len(bps)-1)]
	temp_let=['left']+letters+['right']
	temp_bp=[bps[0]-flank]+bps+[bps[-1]+flank]
	temp_bp=[bp_temp-temp_bp[0] for bp_temp in temp_bp]
	GC_Out={}
	for key in temp_let:
		seq_temp=ori_1_Seq[temp_bp[temp_let.index(key)]:temp_bp[temp_let.index(key)+1]]
		GC_Out[key]=float(seq_temp.count('G')+seq_temp.count('g')+seq_temp.count('C')+seq_temp.count('c'))/float(len(seq_temp))
	return GC_Out

def MoveMP(m,p):
	Move_MP=[m,p]
	Af_BP=[BPList_Rearrange(Be_BP[0],m,original_bp_list),BPList_Rearrange(Be_BP[1],p,original_bp_list)]
	Af_Letter=[LetterList_Rearrange(Be_Letter[0],m,original_bp_list),LetterList_Rearrange(Be_Letter[1],p,original_bp_list)]
	if not Af_Letter_QC(Af_Letter,Copy_num_estimate)==0:
		return 'error'
	else:
		if not Best_Score_Rec==0 and Af_Letter in Best_Letter_Rec: 
			return 'error'
		else:
			letter_num_flag=0
			for key in Block_CN_Upper.keys():
					if (Af_Letter[0]+Af_Letter[1]).count(key)>Block_CN_Upper[key]:
							letter_num_flag+=1
			if not letter_num_flag==0: 
				return 'error'
			else:
				Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
				if Af_Info_all==0:
					return 'error'
				else:
					Letter_Rec.append(Af_Letter)
					BP_Rec.append(Af_BP)
					Af_IL_Penal=Af_Info_all[0]
					Af_RD_Rec=Af_Info_all[1]
					Af_DR_Penal=(Af_Info_all[2])**2
					Af_TB_Penal_a=Af_Info_all[4]
					Af_TB_Rec=Af_Info_all[3]
					Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
					Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
					for key in Af_Info_all[5].keys():
							Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
					P_IL.append(Af_IL_Penal)
					P_RD.append(Af_RD_Penal)
					P_DR.append(Af_DR_Penal/num_of_read_pairs)
					P_TB.append(Af_TB_Penal)
					return 'done'

def RD_Penalize_2c(All_CNRD,GC_Overall_Median_Coverage,GC_Var_Coverage,Chr):
	Theo_RD=GC_Overall_Median_Coverage[str(Chr)]
	Theo_Var=GC_Var_Coverage[str(Chr)]
	Sigma_Square=(Theo_Var-Theo_RD)/Theo_RD**2
	RDs=All_CNRD[1]
	CNs=All_CNRD[0]
	DProb=[]
	MProb=[]
	PProb=[]
	for m in range(len(RDs[0])):
		DProb.append(Prob_Norm(RDs[0][m]+RDs[1][m],Theo_RD*(CNs[0][m]+CNs[1][m]),Theo_RD*(CNs[0][m]+CNs[1][m])+(Theo_RD*(CNs[0][m]+CNs[1][m]))**2*Sigma_Square))
		MProb.append(Prob_Norm(RDs[0][m],Theo_RD*CNs[0][m],Theo_RD*CNs[0][m]+(Theo_RD*CNs[0][m])**2*Sigma_Square))
		PProb.append(Prob_Norm(RDs[1][m],Theo_RD*CNs[1][m],Theo_RD*CNs[1][m]+(Theo_RD*CNs[1][m])**2*Sigma_Square))
	if numpy.sum(DProb)==0 and not numpy.sum(DProb+MProb+PProb)==0:
		DProb[0]=1
	if numpy.sum(MProb)==0 and not numpy.sum(DProb+MProb+PProb)==0:
		MProb[0]=1
	if numpy.sum(PProb)==0 and not numpy.sum(DProb+MProb+PProb)==0:
		PProb[0]=1
	if not numpy.sum(DProb+MProb+PProb)==0:
		return [numpy.mean([n for n in DProb if not n==0]),numpy.mean([n for n in MProb if not n==0]),numpy.mean([n for n in PProb if not n==0])]
	else: 
		return [-100,-100,-100]

def RD_Penalize_2d(All_CNRD,GC_Overall_Median_Coverage,GC_Var_Coverage,Chr):
	Theo_RD=GC_Overall_Median_Coverage[str(Chr)]
	Theo_Var=GC_Var_Coverage[str(Chr)]
	Sigma_Square=(Theo_Var-Theo_RD)/Theo_RD**2
	RDs=All_CNRD[1]
	CNs=All_CNRD[0]
	DProb=[]
	MProb=[]
	PProb=[]
	for m in range(len(RDs[0])):
		DProb.append(Prob_Possion(RDs[0][m]+RDs[1][m],Theo_RD*(CNs[0][m]+CNs[1][m])))
		MProb.append(Prob_Possion(RDs[0][m],Theo_RD*CNs[0][m]))
		PProb.append(Prob_Possion(RDs[1][m],Theo_RD*CNs[1][m]))
	if numpy.sum(DProb)==0 and not numpy.sum(DProb+MProb+PProb)==0:
		DProb[0]=1
	if numpy.sum(MProb)==0 and not numpy.sum(DProb+MProb+PProb)==0:
		MProb[0]=1
	if numpy.sum(PProb)==0 and not numpy.sum(DProb+MProb+PProb)==0:
		PProb[0]=1
	if not numpy.sum(DProb+MProb+PProb)==0:
		return [numpy.mean([n for n in DProb if not n==0]),numpy.mean([n for n in MProb if not n==0]),numpy.mean([n for n in PProb if not n==0])]
	else: 
		return [-100,-100,-100]

def Full_Info_of_Reads_Integrate_old(bps2):
	letter_rd={}
	let_start=96
	for i in bps2:
		letter_rd[i[0]]={}
		for j in range(len(i)-2):
			let_start+=1
			#print chr(let_start)
			if int(i[j+2])-int(i[j+1])>10*flank:
				letter_rd[i[0]][chr(let_start)]=[int(i[j+1])+flank,int(i[j+2])-flank]		
	for i in letter_rd.keys():
		if not letter_rd[i]=={}:
			temp_blockGC=GC_RD_Correction(i)[2]
			for j in letter_rd[i].keys():
				j2=letter_rd[i][j]
				#print j2
				temp_reg1=[]
				for k in temp_blockGC:
					if k[0]<j2[0]+1 and k[1]>j2[1]-1:
						temp_reg1+=k[2:][(j2[0]-k[0])/100:((j2[1]-k[0])/100+1)]
					elif k[0]<j2[0]+1 and k[1]>j2[0]-1 and k[1]<j2[1]+1:
						temp_reg1+=k[2:][(j2[0]-k[0])/100:]
					elif k[0]>j2[0]-1 and k[0]<j2[1]+1 and k[1]>j2[1]-1:
						temp_reg1+=k[2:][:((j2[1]-k[0])/100+1)]
					elif k[0]>j2[0]-1 and k[1]<j2[1]+1:
						temp_reg1+=k[2:]
				j2+=[numpy.mean([float(k2) for k2 in temp_reg1])]
	block_pool={}
	let_start=97
	for i in bps2:
		block_pool[i[0]]=[]
		bp_temp1=[int(j) for j in i[1:]]
		bp_temp2=[[bp_temp1[0]]]
		let_temp3=[]
		for j in bp_temp1[1:]:
			if j-bp_temp2[-1][-1]>10*flank:
				bp_temp2.append([j])
			else:
				bp_temp2[-1].append(j)
		if len(bp_temp2)==1:
			bp_temp3=bp_temp2
			let_temp3.append([chr(let_start+k) for k in range(len(bp_temp2[0])-1)])
			let_start+=len(bp_temp2[0])-1
		else:
			j=bp_temp2[0]
			bp_temp3=[j+[j[-1]+flank]]
			for j in bp_temp2[1:-1]:
				bp_temp3.append([j[0]-flank]+j+[j[-1]+flank])
			j=bp_temp2[-1]
			bp_temp3.append([j[0]-flank]+j)
			for k in bp_temp3:
				let_temp3.append([])
				let_start-=1
				for k2 in range(len(k)-1):
					let_start+=1
					let_temp3[-1].append(chr(let_start))
		block_pool[i[0]].append(bp_temp3)
		block_pool[i[0]].append(let_temp3)
		let_start=ord(let_temp3[-1][-1])+1
	DR_Penal=0
	IL_Penal=[]
	Pair_ThroughBP=[]
	Double_Read_ThroughBP=[]
	Single_Read_ThroughBP=[]
	chr_link={}
	letter2=[]
	for i in bps2:
		for j in block_pool[i[0]][1]:
			letter2+=j
	letter_total=[]
	for i in letter2:
		if not i in letter_total:
			letter_total.append(i)
	letter2=[]
	for i in bps2:
		letter2.append([])
		for j in block_pool[i[0]][1]:
			for k in j:
				if not k in letter2[-1]:
					letter2[-1].append(k)
	bps2LR=[]
	letter2LR=[]
	for i in bps2:
		i_temp=[int(k) for k in i[1:]]
		bps2LR+=[[i_temp[0]-flank]+i_temp+[i_temp[-1]+flank]]
	for i in letter2:
		letter2LR.append(['left']+i+['right'])
	BlockCov={}#coverage of reads total fall within a block
	Initial_Cov={}
	Letter_Double_rec={}
	for i in letter2LR:
		for j in i:
			BlockCov[j]=0
			Initial_Cov[j]=0
	for i in block_pool.keys():
		for j in bps2:
			if i==j[0]:
				bps2_lr=bps2LR[bps2.index(j)]
				letter2_lr=letter2LR[bps2.index(j)]
		for j in range(len(block_pool[i][0])):
			temp_bp_s=block_pool[i][0][j]
			temp_let_s=block_pool[i][1][j]
			temp_info=Full_Info_of_Reads_Product(Initial_Bam,temp_bp_s,bps2_lr,letter2_lr,i,flank,QCAlign,ReadLength,chr_link)
			DR_Penal+=temp_info[0]
			IL_Penal+=temp_info[1]
			Pair_ThroughBP+=temp_info[2]
			Double_Read_ThroughBP+=temp_info[3]
			Single_Read_ThroughBP+=temp_info[4]
			for key in temp_info[5].keys():
				BlockCov[key]+=temp_info[5][key]
			for key in temp_info[6].keys():
				Initial_Cov[key]+=temp_info[6][key]
			for key in temp_info[7].keys():
				if not key in Letter_Double_rec.keys():
					Letter_Double_rec[key]=temp_info[7][key]
				else:
					Letter_Double_rec[key]+=temp_info[7][key]
	temp_info=Single_Rec_Read_Locate(Letter_Double_rec,bps2_lr, letter2_lr)
	DR_Penal+=temp_info[0]
	IL_Penal+=temp_info[1]
	Pair_ThroughBP+=temp_info[2]
	Double_Read_ThroughBP+=temp_info[3]
	Single_Read_ThroughBP+=temp_info[4]
	for key in temp_info[5].keys():
		BlockCov[key]+=temp_info[5][key]
	for key in temp_info[6].keys():
		Initial_Cov[key]+=temp_info[6][key]
	for i in chr_link.keys():
		if len(chr_link[i])==1:
			del chr_link[i]
			continue
		else:
			if chr_link[i][0][1] in block_pool.keys() and chr_link[i][0][5] in block_pool.keys():
				pos1=[chr_link[i][0][1],int(chr_link[i][0][2]),int(chr_link[i][0][2])+cigar2reaadlength(chr_link[i][0][4]),Reads_Direction_Detect(chr_link[i][0][0])[0]]
				pos2=[chr_link[i][1][1],int(chr_link[i][1][2]),int(chr_link[i][1][2])+cigar2reaadlength(chr_link[i][1][4]),Reads_Direction_Detect(chr_link[i][1][0])[0]]
				if pos2[0]==bps2[0][0]:
					pos2=[chr_link[i][0][1],int(chr_link[i][0][2]),int(chr_link[i][0][2])+cigar2reaadlength(chr_link[i][0][4]),Reads_Direction_Detect(chr_link[i][0][0])[0]]
					pos1=[chr_link[i][1][1],int(chr_link[i][1][2]),int(chr_link[i][1][2])+cigar2reaadlength(chr_link[i][1][4]),Reads_Direction_Detect(chr_link[i][1][0])[0]]
				reg1a=Reads_block_assignment_1(bps2LR[0],letter2LR[0],pos1[1]+low_qual_edge)
				reg1b=Reads_block_assignment_1(bps2LR[0],letter2LR[0],pos1[2]-low_qual_edge)
				reg2a=Reads_block_assignment_1(bps2LR[1],letter2LR[1],pos2[1]+low_qual_edge)
				reg2b=Reads_block_assignment_1(bps2LR[1],letter2LR[1],pos2[2]-low_qual_edge)
				if reg1a==reg1b and reg2a==reg2b:
					region1a=bps2LR[0][letter2LR[0].index(reg1a)]
					region2a=bps2LR[1][letter2LR[1].index(reg2a)]
					Pair_ThroughBP.append([reg1a,pos1[1]-region1a,pos1[2]-region1a,reg2a,pos2[1]-region2a,pos2[2]-region2a,pos1[-1],pos2[-1]])
					Initial_Cov[reg1a]+=cigar2reaadlength(chr_link[i][0][4])
					Initial_Cov[reg2a]+=cigar2reaadlength(chr_link[i][1][4])
				else:
					region1a=bps2LR[0][letter2LR[0].index(reg1a)]
					region1b=bps2LR[0][letter2LR[0].index(reg1b)]
					region2a=bps2LR[1][letter2LR[1].index(reg2a)]
					region2b=bps2LR[1][letter2LR[1].index(reg2b)]
					Initial_Cov[reg1a]+=cigar2reaadlength(chr_link[i][0][4])-(pos1[2]-region1b)
					Initial_Cov[reg1b]+=pos1[2]-region1b
					Initial_Cov[reg2a]+=cigar2reaadlength(chr_link[i][1][4])-(pos2[2]-region2b)
					Initial_Cov[reg2b]+=pos2[2]-region2b
					Double_Read_ThroughBP.append([reg1a,pos1[1]-region1a,reg1b,pos1[2]-region1b,reg2a,pos2[1]-region2a,reg2b,pos2[2]-region2b,pos1[-1],pos2[-1]])
	Block_GC={}
	for i in block_pool.keys():
		for j in range(len(block_pool[i][1])):
			for k in block_pool[i][1][j]:
				if not k in Block_GC.keys():
					Block_GC[k]=[]
				region=[block_pool[i][0][j][block_pool[i][1][j].index(k)],block_pool[i][0][j][block_pool[i][1][j].index(k)+1]]
				ref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(refgenome,i, region[0],region[1]))
				ref_seq=[]
				ref.readline().strip().split()
				for line in ref:
					ref_seq+=line.strip().split()
				ref2=''.join(ref_seq)
				if len(ref2)>0:
					Block_GC[k].append((float(ref2.count('c'))+float(ref2.count('g'))+float(ref2.count('C'))+float(ref2.count('G')))/float(len(ref2)))
				else:
					Block_GC[k].append(0.5)
	BlockGC={}
	for i in Block_GC.keys():
		BlockGC[i]=numpy.mean(Block_GC[i])
	for i2 in bps2:
		let_ttt=letter2[bps2.index(i2)]
		adj_ttt=GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,i2[0],[BlockGC[i] for i in let_ttt],[BlockCov[i] for i in let_ttt])
		for j in let_ttt:
			BlockCov[j]=adj_ttt[let_ttt.index(j)]
		adj_ttt=GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,i2[0],[BlockGC[i] for i in let_ttt],[Initial_Cov[i] for i in let_ttt])
		for j in let_ttt:
			Initial_Cov[j]=adj_ttt[let_ttt.index(j)]
	for i in letter_rd.keys():
		for j in letter_rd[i].keys():
			BlockCov[j]+=letter_rd[i][j][-1]*(letter_rd[i][j][1]-letter_rd[i][j][0])
	BlockCov2={}
	Initial_Cov2={}
	for i in letter2LR:
		for j in i[1:-1]:
			BlockCov2[j]=float(BlockCov[j])/float(bps2LR[letter2LR.index(i)][i.index(j)+1]-bps2LR[letter2LR.index(i)][i.index(j)])
			Initial_Cov2[j]=float(Initial_Cov[j])/float(bps2LR[letter2LR.index(i)][i.index(j)+1]-bps2LR[letter2LR.index(i)][i.index(j)])
	bps_total=[int(k) for k in bps2[0][1:]]
	for i in BlockCov2.keys():
		Initial_Cov2[i]+=BlockCov2[i]
	for i in bps2[1:]:
		for j in range(len(i)-2):
			bps_total.append(bps_total[-1]+int(i[j+2])-int(i[j+1]))
	return [BlockCov2,Initial_Cov2,DR_Penal,numpy.mean(IL_Penal),Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,BlockGC,letter_total,bps_total]

def RD_Check_Filein(x):
	RD_Check_Filein=workdir+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/RD_Stat/'+dict_opts['-s'].split('/')[-1].replace('.bam','')+'.'+x[0]+'.RD.index'
	fin=open(RD_Check_Filein)
	start=int(x[1])
	end=int(x[-1])
	while True:
		pin1=fin.readline().strip().split()
		pin2=fin.readline().strip().split()
		pstart=int(pin1[0].split(':')[1].split('-')[0])
		pend=int(pin1[0].split(':')[1].split('-')[1])
		if start>pstart-1 and end<pend+1:
			rd_check=[]
			for y in range(int(float(start-pstart)/100.0)-10,int(float(end-pstart)/100.0)+10):
				rd_check.append(float(pin2[y]))
			cn_check=[]
			for y in rd_check:
				cn_check+=[int(round(y/GC_Mean_Coverage[chrom]*2))]
			ave_mean=[]
			for y in range(len(cn_check)-2):
				ave_mean.append(numpy.mean([cn_check[y],cn_check[y+1],cn_check[y+2]]))
			dire_check=[]
			for y in range(3,(len(rd_check)-3)):
				if abs(cn_check[y]-ave_mean[y-3])<abs(cn_check[y]-ave_mean[y+1]):
					dire_check+=[-1]
				else:
					dire_check+=[1]
			rec=0
			bp_check_rec=[]
			for y in range(len(dire_check)-1):
				if [dire_check[y],dire_check[y+1]]==[-1,1]:
					mean_cn_check=int(round(numpy.mean(cn_check[rec:y+4])))
					if not mean_cn_check==2:
						bp_check_rec.append([rec,y+4])
					rec=y+4
			bps_check_found=[]
			for y in bp_check_rec:
				bps_check_found.append([start-1000+100*(y[0]-1),start-1000+100*(y[1]+1)])

def Chr_Ref_Seq_Produde(ln,Ref_Seq_File):
	ftest=open(Ref_Seq_File)
	ptest=ftest.readline().strip().split()	
	Seq2=[]
	if len(ln[0][0])>3:
		flag=1
	else:
		flag=0
	for i in ln:
		if flag==0:
			Chromo=ptest[0][1:-1]+i[0]
		else:
			Chromo=i[0]
		seq=os.popen(r'''samtools faidx %s %s:%d-%d'''%(Ref_Seq_File,refChr,i[1],i[2]))
		temp=seq.readline().strip().split()
		lines=[]
		while True:
			line=seq.readline().strip().split()
			if not line: break
			lines.append(line)
		seq.close()
		Seq1=lines[0][0]
		for j in range(len(lines))[1:]:
			Seq1=''.join([Seq1,lines[j][0]])
			Seq2.append(Seq1)
	ftest.close()
	return Seq2

def Ref_Seq_Produce(Chromo,bp_list,flank,Ref_Seq_File,num_of_bps):
	#Chromo=Chr, target chromosome
	#BamN: DG187, DG196... name of sample
	#eg of bp_list:[184569179, 184569775, 184571064, 184572009, 184572016]
	#Eg of flank:  flank : 446
	seq_length=numpy.max(bp_list)-numpy.min(bp_list)+2*flank
	num_of_block=seq_length/100+1
	region_left=numpy.min(bp_list)-flank+1
	region_right=numpy.max(bp_list)+flank
	seq=os.popen(r'''samtools faidx %s %s:%d-%d'''%(Ref_Seq_File,refChr,region_left,region_right))
	temp=seq.readline().strip().split()
	lines=[]
	while True:
		line=seq.readline().strip().split()
		if not line: break
		lines.append(line)
	Seq1=lines[0][0]
	for j in range(len(lines))[1:]:
		Seq1=''.join([Seq1,lines[j][0]])
	return Seq1

def Bam_Quality_Control(Initial_Bam,bps,Chr,flank,QCAlign):
	fbam=os.popen('samtools view %s %s:%d-%d'%(Initial_Bam,bamChr,numpy.min(bps)-flank+1,numpy.max(bps)+flank))
	Reads_Name_Unique=[]
	Neg_Reads_InsertLength={}
	Pos_Reads_InsertLength={}
	Neu_Reads_InsertLength={}
	Mate_Unmapped_IL0={}
	Read_Unmapped_IL0={}
	while True:
		pi=fbam.readline().strip().split('\t')[:10]
		if not pi[0]: break
		if not int(pi[4])>QCAlign: continue
		if int(pi[1])&4>0 and int(pi[1])&8>0: continue
		elif int(pi[1])&4>0 and int(pi[1])&8==0:
			Read_Unmapped_IL0[pi[0]]=pi
		elif int(pi[1])&4==0 and int(pi[1])&8>0: 
			Mate_Unmapped_IL0[pi[0]]=pi
		elif int(pi[1])&4==0 and int(pi[1])&8==0:
			if not pi[0] in Reads_Name_Unique:
				Reads_Name_Unique.append(pi[0])			
			if int(pi[8])<0:
				Neg_Reads_InsertLength[pi[0]]=pi
			elif int(pi[8])>0:
				Pos_Reads_InsertLength[pi[0]]=pi
			else:
				Neu_Reads_InsertLength[pi[0]]=pi
	return [Reads_Name_Unique,Neg_Reads_InsertLength,Pos_Reads_InsertLength,Neu_Reads_InsertLength,Mate_Unmapped_IL0,Read_Unmapped_IL0]


