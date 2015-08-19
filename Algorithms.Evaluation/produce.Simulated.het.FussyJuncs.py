#!/usr/bin/env python

#!python
#command='python SVelter.py --ref /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -i /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.het/het.sim -o /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.het/het'
#command='python SVelter.py --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.fasta -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/het.sim -o /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/het'
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
def sv_rec_2(sv_info):
	for k1ab in sorted(sv_info.keys()):
		for k2ab in sv_info[k1ab].keys():
			if not k2ab==k1ab:
				k1aba=k1ab.split('/')[0]
				k2aba=k2ab.split('/')[0]
				k2abb=k2ab.split('/')[1]
				flaga=[]
				flagb=[]
				test=[[],[]]
				if flaga==[] and not k1aba==k2aba:
					if k2aba=='':
						csv1=[[i for i in k1aba],[],[],0]
					else:
						csv1=simple_flag_SA(k1aba,k2aba)
					add_csv_info(csv1,1,k1ab,k2ab)
				if flagb==[] and not k1aba==k2abb:
					if k2abb=='':
						csv1=[[i for i in k2abb],[],[],0]
					else:
						csv1=simple_flag_SA(k1aba,k2abb)
					add_csv_info(csv1,2,k1ab,k2ab)

score_Cff=-20
def hash_reorder():
	for ka1 in del1.keys():
		if not ka1 in sv_out.keys():
			sv_out[ka1]={}
		for ka2 in del1[ka1]:
			#fref=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,ka1,str(ka2[0]+1),str(ka2[0]+1)))
			#tre=fref.readline().strip().split()
			#REF_AL=fref.readline().strip().split()[0]
			REF_AL='N'
			Pass_Sign='PASS'
			if ka2[3]<score_Cff:
				Pass_Sign='LowQual'
			if ka2[2]=='heta':
				GenoType='1|0'
			elif ka2[2]=='hetb':
				GenoType='0|1'
			elif ka2[2]=='homo':
				GenoType='1|1'
			ka_new=[ka1,ka2[0],ka2[-1],REF_AL,'<DEL>',ka2[3],Pass_Sign,'SVTYPE=DEL;END='+str(ka2[1]),'GT',GenoType]
			if not ka2[-1] in sv_out[ka1].keys():
				sv_out[ka1][ka2[-1]]=[]
			if not ka_new in sv_out[ka1][ka2[-1]]:
				sv_out[ka1][ka2[-1]].append(ka_new)
	for ka1 in inv1.keys():
		if not ka1 in sv_out.keys():
			sv_out[ka1]={}
		for ka2 in inv1[ka1]:
			#fref=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,ka1,str(ka2[0]+1),str(ka2[0]+1)))
			#tre=fref.readline().strip().split()
			#REF_AL=fref.readline().strip().split()[0]
			REF_AL='N'
			Pass_Sign='PASS'
			if ka2[3]<score_Cff:
				Pass_Sign='LowQual'
			if ka2[2]=='heta':
				GenoType='1|0'
			elif ka2[2]=='hetb':
				GenoType='0|1'
			elif ka2[2]=='homo':
				GenoType='1|1'
			ka_new=[ka1,ka2[0],ka2[-1],REF_AL,'<INV>',ka2[3],Pass_Sign,'SVTYPE=INV;END='+str(ka2[1]),'GT',GenoType]
			if not ka2[-1] in sv_out[ka1].keys():
				sv_out[ka1][ka2[-1]]=[]
			if not ka_new in sv_out[ka1][ka2[-1]]:
				sv_out[ka1][ka2[-1]].append(ka_new)
	for ka1 in dup1.keys():
		if not ka1 in sv_out.keys():
			sv_out[ka1]={}
		for ka2 in dup1[ka1]:
			#fref=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,ka1,str(ka2[0]+1),str(ka2[0]+1)))
			#tre=fref.readline().strip().split()
			#REF_AL=fref.readline().strip().split()[0]
			REF_AL='N'
			CopyNumber=str(ka2[-1])
			Pass_Sign='PASS'
			if ka2[3]<score_Cff:
				Pass_Sign='LowQual'
			if ka2[2]=='heta':
				GenoType='1|0'
			elif ka2[2]=='hetb':
				GenoType='0|1'
			elif ka2[2]=='homo':
				GenoType='1|1'
			ka_new=[ka1,ka2[0],ka2[-2],REF_AL,'<DUP>',ka2[3],Pass_Sign,'SVTYPE=DUP;END='+str(ka2[1]),'GT:CN',GenoType+':'+CopyNumber]
			if not ka2[-2] in sv_out[ka1].keys():
				sv_out[ka1][ka2[-2]]=[]
			if not ka_new in sv_out[ka1][ka2[-2]]:
				sv_out[ka1][ka2[-2]].append(ka_new)
	for ka1 in tra1.keys():
		ks1=ka1.split('_')[0]
		ks2='_'.join(ka1.split('_')[:-1])
		SV_Score=float(ka1.split('_')[-1])
		Pass_Sign='PASS'
		if SV_Score<score_Cff:
			Pass_Sign='LowQual'
		if not ks1 in sv_out.keys():
			sv_out[ks1]={}
		if not ks2 in sv_out[ks1].keys():
			sv_out[ks1][ks2]=[]
		for ka2 in tra1[ka1].keys():
			hetx='het'+ka2
			if ka2=='a':
				GenoType='1|0'
			elif ka2=='b':
				Genotype='0|1'
			for ka3 in tra1[ka1][ka2]:
				ka_new=ka3[:2]+[ks2,ka3[2]]+ka3[3:]+[SV_Score,Pass_Sign,'SVTYPE=TRA','GT',GenoType]
				if not ka_new in sv_out[ks1][ks2]:
					sv_out[ks1][ks2].append(ka_new)

def write_VCF_header(output_file):
	fo=open(output_file,'w')
	print output_file
	print>>fo, '##fileformat=VCFv4.1'
	print>>fo,'##fileDate='+time.strftime("%Y%m%d")
	print>>fo,'##reference=hg19'
	print>>fo,'##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">'
	print>>fo,'##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">'
	print>>fo,'##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">'
	print>>fo,'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
	print>>fo,'##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">'
	print>>fo,'##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">'
	print>>fo,'##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
	print>>fo,'##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">'
	print>>fo,'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
	print>>fo,'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
	print>>fo,'##FILTER=<ID=LowQual,Description="Score of final structural - Theoretical Score <-50">'
	print>>fo,'##ALT=<ID=DEL,Description="Deletion">'
	print>>fo,'##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">'
	print>>fo,'##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">'
	print>>fo,'##ALT=<ID=DUP,Description="Duplication">'
	print>>fo,'##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">'
	print>>fo,'##ALT=<ID=INS,Description="Insertion of novel sequence">'
	print>>fo,'##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">'
	print>>fo,'##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">'
	print>>fo,'##ALT=<ID=INV,Description="Inversion">'
	print>>fo,'##ALT=<ID=CNV,Description="Copy number variable region">'
	print>>fo,'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
	print>>fo,'##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">'
	print>>fo,'##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">'
	print>>fo,'##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">'
	print>>fo,'\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',output_file.split('/')[-1].replace('.vcf','')])
	fo.close()

def write_VCF_main(output_file):
	fo=open(output_file,'a')
	print output_file
	sv_reorganize={}
	for k1 in sv_out.keys():
		sv_reorganize[k1]={}
		for k2 in sv_out[k1].keys():
			start=int(k2.split('_')[1])
			if not start in sv_reorganize[k1].keys():
				sv_reorganize[k1][start]={}
			SVtemp_a=[]
			SVtemp_b=[]
			for k3 in sv_out[k1][k2]:
				if not k3[:-1] in SVtemp_a:
					SVtemp_a.append(k3[:-1])
					SVtemp_b.append([k3[-1]])
				else:
					SVtemp_b[SVtemp_a.index(k3[:-1])].append(k3[-1])
			SVtemp=[]
			sv_reorganize[k1][start][k2]=[]
			for k3 in range(len(SVtemp_a)):
				if len(SVtemp_b[k3])==2 and SVtemp_b[k3] in [['0|1', '1|0'],['1|0', '0|1']]:
					SVtemp_b[k3]=['1|1']
			for k3 in range(len(SVtemp_a)):
				for k4 in SVtemp_b[k3]:
					sv_reorganize[k1][start][k2].append(SVtemp_a[k3]+[k4])
	for k1 in chromos:
		if k1 in sv_reorganize.keys():
			for k2 in sorted(sv_reorganize[k1].keys()):
				for k3 in sorted(sv_reorganize[k1][k2].keys()):
					for k4 in sv_reorganize[k1][k2][k3]:
						if k4[3]=='N':
							k4[3]=ref_base_returnN(ref,k4[0],k4[1])
						print >>fo, '\t'.join([str(i) for i in k4])
	fo.close()

def simple_flag_SA(k1,k2):
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
				if temp2.count(ia)+temp2.count(ib)>1:
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
	return [outdel,outinv,outdup2,outtra]

def add_csv_info(csv1,flag_sex,k1,k2):
	#flag_sex=1: Maternal
	#flag_sex=2: Paternal
	if flag_sex==1:
		del_let=[csv1[0],[]]
		inv_let=[csv1[1],[]]
		dup_let=[csv1[2],[]]
	else:
		del_let=[[],csv1[0]]
		inv_let=[[],csv1[1]]
		dup_let=[[],csv1[2]]
	for k3 in sv_info[k1][k2]:
		del_info_add(k3,del_let)
		inv_info_add(k3,inv_let)
		dup_info_2_add(k3,dup_let)
	if csv1[3]==1:
		tra_info_add(k1,k2)

def del_info_add(k3,del_let):
	tempa=bp_to_hash(k3[:-1],del_let[0])
	tempb=bp_to_hash(k3[:-1],del_let[1])
	for k1 in tempa:
	    if k1 in tempb:
	        tempc='hom'
	        tempb.remove(k1)
	    else:
	        tempc='heta'
	    if not k1[0] in del1.keys():
	        del1[k1[0]]=[]
	    del1[k1[0]].append(k1[1:]+[tempc,k3[-1],'_'.join(k3[:-1])])
	for k1 in tempb:
	    if not k1[0] in del1.keys():
	        del1[k1[0]]=[]
	    del1[k1[0]].append(k1[1:]+['hetb',k3[-1],'_'.join(k3[:-1])])

def dup_info_add(k3,dup_let):
    #dup_let=[k2i,k2j]
    for k2x in dup_let:
        for k4 in k2x:
            temp=bp_to_hash(k3[:-1],[i for i in k4])
            for k5 in temp:
                if not k5[0] in dup1.keys():
                    dup1[k5[0]]=[]
                dup1[k5[0]].append(k5[1:]+[k3[-1],'_'.join(k3[:-1]),k2a.count(k4)])

def dup_info_2_add(k3,dup_let):
	temprec=-1
	for k2x in dup_let:
		temprec+=1
		hetx=['heta','hetb'][temprec]
		for k4 in k2x:
			temp=bp_to_hash(k3[:-1],[i for i in k4[0]])
			for k5 in temp:
			    if not k5[0] in dup1.keys():
			        dup1[k5[0]]=[]
			    if k4[1]>1:
				    dup1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]),k4[1]])

def inv_info_add(k3,inv_let):
    #inv_let=[k2m,k2n]
    temprec=-1
    for k2x in inv_let:
    	temprec+=1
    	hetx=['heta','hetb'][temprec]
        for k4 in k2x:
            temp=bp_to_hash(k3[:-1],[i for i in k4])
            for k5 in temp:
                if not k5[0] in inv1.keys():
                    inv1[k5[0]]=[]
                inv1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1])])

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

def dup_let_recombind(vec_in):
    if vec_in==[]:
        return []
    else:
        vec2=sorted(vec_in)
        vec=[[vec2[0]]]
        for ka in vec2[1:]:
            if ord(ka)-ord(vec[-1][-1])==1:
                vec[-1].append(ka)
            else:
                vec.append([ka])
        vec3=[]
        for ka in vec:
            if len(ka)==1:
                vec3.append(ka)
            else:
                for kb in range(2,len(ka)+1):
                    for kc in ka[:(1-kb)]:
                        vec3.append([])
                        for kd in range(kb):
                            vec3[-1].append(ka[ka.index(kc)+kd])                
        vec4=[''.join(i) for i in vec3]
        return vec4

def comp_info_reorganize(k1,k2):
    del_let=[[],[]]
    dup_let=[[],[]]
    inv_let=[[],[]]
    tra_let=[[],[]]
    k2a=k2.split('/')[0]
    k2b=k2.split('/')[1]
    k2c=[]
    k2d=[]
    for k3 in k2a:
        if not k3=='^':
            k2c.append(k3)
        else:
            k2c[-1]+=k3
    for k3 in k2b:
        if not k3=='^':
            k2d.append(k3)
        else:
            k2d[-1]+=k3
    for k3 in k1.split('/')[0]:
        if k2a.count(k3)==0:
            del_let[0].append(k3)
        if k2b.count(k3)==0:
            del_let[1].append(k3)
        if k2a.count(k3)>1:
            dup_let[0].append(k3)
        if k2b.count(k3)>1:
            dup_let[1].append(k3)
    k2e=let_reclust(k2c)
    k2f=let_reclust(k2d)
    k2g=dup_let_recombind(dup_let[0])
    k2h=dup_let_recombind(dup_let[1])
    k2i=[]
    #integreated dup sections
    k2j=[]
    #integreated dup sections
    for k3 in k2g:
        flag1=0
        for k4 in k2e:
            if k3 in k4:
                flag1+=1
        if flag1>1:
            k2i.append(k3)
    for k3 in dup_let[0]:
        if k2e.count(k3[0])+k2e.count(k3[0]+'^')>0:
            if not k3[0] in k2i:
                k2i.append(k3[0])
    for k3 in k2h:
        flag1=0
        for k4 in k2e:
            if k3 in k4:
                flag1+=1
        if flag1>1:
            k2j.append(k3)
    for k3 in dup_let[1]:
        if k2e.count(k3[0])+k2e.count(k3[0]+'^')>0:
            if not k3[0] in k2j:
                k2j.append(k3[0])
    k2m=[]
    for k3 in k2e:
        if k3[-1]=='^':
            k2m.append(k3)
    k2n=[]
    for k3 in k2f:
        if k3[-1]=='^':
            k2n.append(k3)
    for k3 in sv_info[k1][k2]:
        del_info_add(k3,del_let)
        dup_info_add(k3,[k2i,k2j])
        inv_info_add(k3,[k2m,k2n])

def bp_to_hash(bp_list,sv_let):
	bp_hash={}
	block_rec=0
	block_hash=[]
	sv_let=[i[0] for i in sv_let]
	for a3 in bp_list:
	    if a3 in chromos or not a3.isdigit():
	        block_hash.append([a3])
	    else:
	        block_hash[-1].append(a3)
	for a3 in block_hash:
	    for a4 in range(len(a3)-2):
	        bp_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
	        block_rec+=1
	out=[]
	if not sv_let==[]:
	    if len(sv_let)==1:
	        out=[bp_hash[sv_let[0]]]
	    else:
	        out.append(bp_hash[sv_let[0]])
	        for ka in range(len(sv_let)-1):
	            if ord(sv_let[ka+1])-ord(sv_let[ka])==1 and bp_hash[sv_let[ka+1]][0]==bp_hash[sv_let[ka]][0]:
	                out[-1]+=bp_hash[sv_let[ka+1]][1:]
	            else:
	                out.append(bp_hash[sv_let[ka+1]])
	out2=[]
	for ka in out:
	    out2.append([ka[0],int(ka[1]),int(ka[-1])])
	return out2

def tra_info_add(k1,k2):
	for k3 in sv_info[k1][k2]:
		SV_ID='_'.join([str(i) for i in k3])
		tra1[SV_ID]={}
		k2a=k2.split('/')[0]
		k2b=k2.split('/')[1]
		bp_hash={}
		block_rec=0
		block_hash=[]
		for a3 in k3[:-1]:
		    if a3 in chromos or not a3.isdigit():
		        block_hash.append([a3])
		    else:
		        block_hash[-1].append(a3)
		for a3 in block_hash:
		    for a4 in range(len(a3)-2):
		        bp_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
		        block_rec+=1
		for a3 in bp_hash.keys():
		    temp=[]
		    for a4 in bp_hash[a3][1:]:
		        temp.append(int(a4)-1)
		        temp.append(int(a4))
		    bp_hash[a3][1:]=temp
		#ref_allele['left']=[ref_allele[k1[0]][0]]
		#ref_allele['right']=[ref_allele[k1[-1]][1]]
		bp_hash['left']=[bp_hash[k1[0]][0],bp_hash[k1[0]][1],bp_hash[k1[0]][2]]
		bp_hash['right']=[bp_hash[k1[-1]][0],bp_hash[k1[-1]][3],bp_hash[k1[-1]][4]]
		ref_allele={}
		for a3 in bp_hash.keys():
			ref_allele[a3]=[bp_hash[a3][0]]
			for a4 in bp_hash[a3][1:]:
				ref_allele[a3].append(ref_base_returnN(ref,bp_hash[a3][0],a4))
		if not k2a==k1.split('/')[0] and del_flag_SA(k1.split('/')[0],k2a)==0:
		    flag1=0#flag1==0:w/o inversion in the alt structure
		    if '^' in k2a:
		        flag1+=1
		    flag2=0#flag2==0:w/o duplication in the alt structure
		    for j in k2a:
		        if k2a.count(j)>1:
		            flag2+=1
		    flag3=0 #flag3==0: w/o translocation
		    if len(k2a)>1:
		        for i in range(len(k2a)-1):
		            if not ord(k2a[i+1])>ord(k2a[i]):
		                flag3+=1
		    if flag1+flag2+flag3==0:
		        heta_Del_block=[]
		        for a1 in k1.split('/')[0]:
		            if not a1 in k2a:
		                heta_Del_block.append(a1)      
		        tra1[SV_ID]['a']=[]
		        block_hash=[]
		        del_hash={}
		        block_rec=0
		        for a3 in a2[0]:
		            if a3 in chromos:
		                block_hash.append([a3])
		            else:
		                block_hash[-1].append(a3)
		        for a3 in block_hash:
		            for a4 in range(len(a3)-2):
		                del_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
		                block_rec+=1                                
		        if not heta_Del_block==[]:
		            a_heta=0
		            heta_Del_new=[heta_Del_block[0]]
		            while True:
		                a_heta+=1
		                if a_heta==len(heta_Del_block):break
		                if ord(heta_Del_block[a_heta])-ord(heta_Del_block[a_heta-1])==1 and del_hash[heta_Del_block[a_heta]][0]==del_hash[heta_Del_block[a_heta-1]][0]:
		                    heta_Del_new[-1]+=heta_Del_block[a_heta]
		                else:
		                    heta_Del_new.append(heta_Del_block[a_heta])
		            for a3 in heta_Del_new:
		                a4=a3[0]
		                tra1[SV_ID]['a'].append(['DEL',del_hash[a4][0],int(del_hash[a4][1]),ref_allele[a4][2]])
		                a4=a3[-1]
		                tra1[SV_ID]['a'][-1].append(int(del_hash[a4][2])-1)
		    else:
		        tra1[SV_ID]['a']=[]
		        t1=[]
		        for a3 in k2a:
		            if not a3=='^':
		                t1.append(a3)
		            else:
		                t1[-1]+=a3
		        t2=[t1[0]]
		        for a3 in t1[1:]:
		            if not '^' in a3 and not '^' in t2[-1] and ord(a3)-ord(t2[-1][-1])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-1]][0]:
		                t2[-1]+=a3
		            elif '^' in a3 and '^' in t2[-1] and ord(t2[-1][-2])-ord(a3[0])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-2]][0]:
		                t2[-1]+=a3
		            else:
		                t2.append(a3)
		        a3='left'
		        a4=t2[0]
		        l_chr=bp_hash[a3][0]
		        r_chr=bp_hash[a4[0]][0]
		        if not '^' in a4:
		            if not a4[0]==k1[0]:
			            tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4[0]][2]])
			            tra1[SV_ID]['a'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4[0]][2])+'['])
		        elif '^' in a4:
		            tra1[SV_ID]['a'].append([r_chr, bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3][1])+']'])
		            tra1[SV_ID]['a'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
		        for t3 in range(len(t2)-1):
		            a3=t2[t3]
		            a4=t2[t3+1]
		            l_chr=bp_hash[a3[0]][0]
		            r_chr=bp_hash[a4[0]][0]
		            if not '^' in a3 and not '^' in a4:
		                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4[0]][2]])
		                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['])
		            elif '^' in a3 and not '^' in a4:
		                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4[0]][2]])
		                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['+ref_allele[a3[-2]][2]])
		            elif not '^' in a3 and '^' in a4:
		                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'])
		                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
		            elif '^' in a3 and '^' in a4:
		                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['])
		                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2], ']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'+ref_allele[a3[-2]][2]])
		        if len(t2)>1:
		            a3=t2[t3+1]
		        else:
		            a3=t2[0]
		        a4='right'
		        l_chr=bp_hash[a3[0]][0]
		        r_chr=bp_hash[a4][0]
		        if not '^' in a3:
		            if not a3[-1]==k1[-1]:
		                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4][2]])
		                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['])
		        if '^' in a3:
		            tra1[SV_ID]['a'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4][2]])
		            tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['+ref_allele[a3[-2]][2]])
		        #print [k1,k2]
		if not k2b==k1.split('/')[1] and del_flag_SA(k1.split('/')[1],k2b)==0:
		    flag1=0#flag1==0:w/o inversion in the alt structure
		    if '^' in k2b:
		        flag1+=1
		    flag2=0#flag2==0:w/o duplication in the alt structure
		    for j in k2b:
		        if k2b.count(j)>1:
		            flag2+=1
		    flag3=0 #flag3==0: w/o translocation
		    if len(k2b)>1:
		        for i in range(len(k2b)-1):
		            if not ord(k2b[i+1])>ord(k2b[i]):
		                flag3+=1
		    if flag1+flag2+flag3==0:
		        heta_Del_block=[]
		        for a1 in k1.split('/')[1]:
		            if not a1 in k2b:
		                heta_Del_block.append(a1)      
		        tra1[SV_ID]['b']=[]
		        block_hash=[]
		        del_hash={}
		        block_rec=0
		        for a3 in a2[0]:
		            if a3 in chromos:
		                block_hash.append([a3])
		            else:
		                block_hash[-1].append(a3)
		        for a3 in block_hash:
		            for a4 in range(len(a3)-2):
		                del_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
		                block_rec+=1                                
		        if not heta_Del_block==[]:
		            a_heta=0
		            heta_Del_new=[heta_Del_block[0]]
		            while True:
		                a_heta+=1
		                if a_heta==len(heta_Del_block):break
		                if ord(heta_Del_block[a_heta])-ord(heta_Del_block[a_heta-1])==1 and del_hash[heta_Del_block[a_heta]][0]==del_hash[heta_Del_block[a_heta-1]][0]:
		                    heta_Del_new[-1]+=heta_Del_block[a_heta]
		                else:
		                    heta_Del_new.append(heta_Del_block[a_heta])
		            for a3 in heta_Del_new:
		                a4=a3[0]
		                tra1[SV_ID]['b'].append(['DEL',del_hash[a4][0],int(del_hash[a4][1]),ref_allele[a4][2]])
		                a4=a3[-1]
		                tra1[SV_ID]['b'][-1].append(int(del_hash[a4][2])-1)
		    else:
		        tra1[SV_ID]['b']=[]
		        t1=[]
		        for a3 in k2b:
		            if not a3=='^':
		                t1.append(a3)
		            else:
		                t1[-1]+=a3
		        t2=[t1[0]]
		        for a3 in t1[1:]:
		            if not '^' in a3 and not '^' in t2[-1] and ord(a3)-ord(t2[-1][-1])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-1]][0]:
		                t2[-1]+=a3
		            elif '^' in a3 and '^' in t2[-1] and ord(t2[-1][-2])-ord(a3[0])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-2]][0]:
		                t2[-1]+=a3
		            else:
		                t2.append(a3)
		        a3='left'
		        a4=t2[0]
		        l_chr=bp_hash[a3][0]
		        r_chr=bp_hash[a4[0]][0]
		        if not '^' in a4:
		            if not a4[0]==k1[0]:
		                 tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4[0]][2]])
		                 tra1[SV_ID]['b'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4[0]][2])+'['])
		        elif '^' in a4:
		            tra1[SV_ID]['b'].append([r_chr, bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3][1])+']'])
		            tra1[SV_ID]['b'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
		        for t3 in range(len(t2)-1):
		            a3=t2[t3]
		            a4=t2[t3+1]
		            l_chr=bp_hash[a3[0]][0]
		            r_chr=bp_hash[a4[0]][0]
		            if not '^' in a3 and not '^' in a4:
		                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4[0]][2]])
		                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['])
		            elif '^' in a3 and not '^' in a4:
		                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4[0]][2]])
		                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['+ref_allele[a3[-2]][2]])
		            elif not '^' in a3 and '^' in a4:
		                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'])
		                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
		            elif '^' in a3 and '^' in a4:
		                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['])
		                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2], ']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'+ref_allele[a3[-2]][2]])
		        if len(t2)>1:
		            a3=t2[t3+1]
		        else:
		            a3=t2[0]                                
		        a4='right'
		        l_chr=bp_hash[a3[0]][0]
		        r_chr=bp_hash[a4][0]
		        if not '^' in a3:
		            if not a3[-1]==k1[-1]:
		                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4][2]])
		                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['])
		        if '^' in a3:
		            tra1[SV_ID]['b'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4][2]])
		            tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['+ref_allele[a3[-2]][2]])

def sv_homo_initial():
	sv_homo_info['DEL']=[]
	sv_homo_info['DUP']=[]
	sv_homo_info['INV']=[]
	sv_homo_info['TRA']=[]
	sv_homo_info['DUP:TANDOM']=[]

def produce_keys(key):
	if key=='DEL':
		ka='a/a'
		kb='/'
	elif key=='DUP:TANDOM':
		ka='a/a'
		dup_num=random.sample(range(2,20),1)
		kb='/'.join([''.join(['a' for i in range(dup_num[0])]),''.join(['a' for i in range(dup_num[0])])])
	elif key=='INV':
		ka='a/a'
		kb='a^/a^'
	elif key=='TRA':
		ka='ab/ab'
		kb='ba/ba'
	elif key=='DUP':
		ka='ab/ab'
		kb='aba/aba'
	return [ka,kb]

def sv_homo_produce():
	for k1 in SV_region:
		sv_len=k1[2]-k1[1]
		k2=k1[-1]
		sv_homo_info[k2].append(k1+produce_keys(k2))

def sv_het_produce():
	for k1 in sv_homo_info.keys():
		sv_het_info[k1]=[]
		for k2 in sv_homo_info[k1]:
			allele=random.choice(range(2))
			alle_poor=[k2[-2].split('/')[0],k2[-1].split('/')[0]]
			k2[-1]='/'.join([alle_poor[allele],alle_poor[1-allele]])
			sv_het_info[k1].append(k2)

def sv_rec_homo_produce():
	for k1 in sv_homo_info.keys():
		fo=open(dict_opts['-o']+'.homo.'+k1+'.rec','w')
		print dict_opts['-o']+'.homo.'+k1+'.rec'
		for k2 in sv_homo_info[k1]:
			print >>fo, ' '.join([str(i) for i in k2])
		fo.close()

def sv_rec_het_produce():
	for k1 in sv_het_info.keys():
		fo=open(dict_opts['-o']+'.het.'+k1+'.rec','w')
		print dict_opts['-o']+'.het.'+k1+'.rec'
		for k2 in sv_het_info[k1]:
			print >>fo, ' '.join([str(i) for i in k2])
		fo.close()

def sv_info_rewrite(sv_h_info):
	for k1 in sv_h_info.keys():
		for k2 in sv_h_info[k1]:
			if not k2[-2] in sv_info.keys():
				sv_info[k2[-2]]={}
			if not k2[-1] in sv_info[k2[-2]].keys():
				sv_info[k2[-2]][k2[-1]]=[]
			sv_info[k2[-2]][k2[-1]].append([str(i) for i in k2[:-3]]+[0.0])

def sv_stat_calcu(sv_hash,key):
	out=[]
	for k1 in sv_hash[key]:		
		sv_min=int(k1[1])
		sv_max=int(k1[2])
		sv_int=(int(k1[2])-int(k1[1]))/3
		out.append([k1[0],sv_min,sv_min+sv_int, sv_max-sv_int,sv_max])
	return out

def sv_size_pick(sv_stat):
	out=[]
	for k1 in sv_stat:
		out+=[random.choice(range(int(k1[1]),int(k1[2]))) for i in range(int(k1[0]/3))]
		out+=[random.choice(range(int(k1[2]),int(k1[3]))) for i in range(int(int(k1[0])-int(k1[0]/3))/2)]
		out+=[random.choice(range(int(k1[3]),int(k1[4]))) for i in range(int(k1[0])-int(k1[0]/3)-int(int(k1[0])-int(k1[0]/3))/2)]
	permute=random.sample(out,len(out))
	return out

def chromos_readin(refs):
	fin=open(refs+'.fai')
	chromos=[]
	chromo_length=[]
	genome_length=0
	for line in fin:
		pin=line.strip().split()
		chromos.append(pin[0])
		genome_length+=int(pin[1])
		chromo_length.append(int(pin[1]))
	fin.close()
	chromo_num_region=[]
	for k1 in chromo_length:
		chromo_num_region.append(int(round(float(k1)/float(genome_length)*sv_total_num)))
	chrom_to_remove=[]
	out_num_region=[]
	out_chromos=[]
	out_length={}
	for i in range(len(chromo_num_region)):
		if chromo_num_region[i]>1:
			out_chromos.append(chromos[i])
			out_num_region.append(chromo_num_region[i])
			out_length[chromos[i]]=chromo_length[i]
	return [genome_length]+[out_chromos]+[out_num_region]+[out_length]

def sv_hash_add(list_in,key):
	for i in list_in:
		if not i in sv_hash.keys():
			sv_hash[i]=[key]
		else:
			sv_hash[i]+=[key]

def sv_region_pick():
	#pick random regions across the genome
	SV_region=[]
	rec=-1
	sv_size=del_size+dup_size+inv_size+tra_size+dup2_size
	sv_size=random.sample(sv_size,len(sv_size))
	for k1 in range(len(chromos)):
		chromosome=chromos[k1]
		num_region=chromo_num_region[k1]
		range_region=chromo_length[chromosome]
		temp_start_region=sorted(random.sample(range(1000, range_region-1000),num_region+1))
		temp_end_region=[]
		for k2 in range(num_region):
			start=temp_start_region[k2]
			start2=temp_start_region[k2+1]
			if start2-start<1000: continue
			rec+=1
			temp_sv_size=sv_size[rec]
			sv_type=sv_hash[sv_size[rec]][0]
			del sv_hash[sv_size[rec]][0]
			end=start+temp_sv_size
			if not end<start2-300:
				end=random.choice(range(start,int(numpy.mean([start,start2]))))
			if sv_type=='TRA':
				end2=random.choice(range(end+100,start2-100))
			temp_end_region.append(end)
			if sv_type=='TRA':
				SV_region.append([chromos[k1],start,end,end2,sv_type])
			else:
				SV_region.append([chromos[k1],start,end,sv_type])
	return SV_region

def ref_base_returnN(ref,chromo,pos):
	return 'N'

def ref_base_readin(ref,chromo,pos):
	fref=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,chromo,str(pos),str(pos)))
	tre=fref.readline().strip().split()
	REF_AL=fref.readline().strip().split()
	if not REF_AL==[]:
		return REF_AL[0]
	else:
		return 'N'

def del_flag_SA(k1,k2):
	out=0
	if not '^' in k2:
		flagdup=0
		for i in k2:
			if k2.count(i)>1:
				flagdup+=1
		if flagdup==0:
			flagtra=0
			for i in range(len(k2)-1):
				if ord(k2[i+1])-ord(k2[i])<1:
					flagtra+=1
			if flagtra==0:
				if not k1==k2:
					out=1	
	return out

def order_SV_Homo_write(sv_info):
	for k1 in sv_info.keys():
		for k2 in sv_info[k1].keys():
			for k3 in sv_info[k1][k2]:
				if not k3[0] in order_SV_Pos.keys():
					order_SV_Pos[k3[0]]={}
				if not int(k3[1]) in order_SV_Pos[k3[0]].keys():
					order_SV_Pos[k3[0]][int(k3[1])]=[]
				order_SV_Pos[k3[0]][int(k3[1])].append([[k3[0]]+[int(i) for i in k3[1:-1]],[k2.split('/')[0]]])

def order_SV_Het_write(sv_info):
	for k1 in sv_info.keys():
		for k2 in sv_info[k1].keys():
			for k3 in sv_info[k1][k2]:
				if not k3[0] in order_SV_Pos.keys():
					order_SV_Pos[k3[0]]={}
				if not int(k3[1]) in order_SV_Pos[k3[0]].keys():
					order_SV_Pos[k3[0]][int(k3[1])]=[]
				order_SV_Pos[k3[0]][int(k3[1])].append([[k3[0]]+[int(i) for i in k3[1:-1]],[k2.split('/')[0],k2.split('/')[1],k1.split('/')[0]]])

def insert_read_decide(bp_list):
	#decide which class to simulate, ClassI~71%, ClassII~29%
	SV_class_decide=random.choice(range(100)) 
	if SV_class_decide>70:#if ClassII
		sub_class_decide=random.choice(range(100)) 
		if sub_class_decide<60:#2-20bp micro insertion of random seqs
			return produce_random_seqs(random.choice(range(2,20)))
		else:	#over 20bp insertion
			sub2_class_decide=random.choice(range(100)) 
			if sub2_class_decide<25: #25%, 20-50bp random seqs 
				return produce_random_seqs(random.choice(range(20,50)))
			elif sub2_class_decide<50: #25%, 20-50bp seqs from another chromosome
				temp=[]
				for x in seq_ins_pools.keys():
					if not x==bp_list[0]:
						temp.append(x)
				return random.choice(seq_ins_pools[random.choice(temp)])
			else: #50%, 20-50bp seqs from the same chromosome
				return random.choice(seq_ins_pools[bp_list[0]])
	else:#if ClassI
		return ''

def Ref_Alt_Produce(ChromoList,bp_list,letter_new,Ref_Seq_File):
	#Chromo=Chr, target chromosome
	#BamN: DG187, DG196... name of sample
	#eg of bp_list:[184569179, 184569775, 184571064, 184572009, 184572016]
	#Eg of flank:  flank : 446
	if letter_new=='':
		return insert_read_decide(bp_list)
	else:
		bp_hash={}
		bp_seq=[]
		for k1 in bp_list:
			if k1 in ChromoList:
				bp_seq.append([k1])
			else:
				bp_seq[-1].append(k1)
		rec=0
		for k1 in bp_seq:
			for k2 in range(len(k1)-2):
				rec+=1
				bp_hash[chr(96+rec)]=[k1[0],k1[k2+1],k1[k2+2]]
		letter_seq={}
		for k1 in bp_hash.keys():
			Chromo=bp_hash[k1][0]
			region_left=bp_hash[k1][1]
			region_right=bp_hash[k1][2]
			seq=os.popen(r'''samtools faidx %s %s:%d-%d'''%(Ref_Seq_File,Chromo,region_left,region_right))
			seq.readline().strip().split()
			lines=[]
			while True:
				line=seq.readline().strip().split()
				if not line: break
				lines.append(line)
			Seq1=lines[0][0]
			for j in range(len(lines))[1:]:
				Seq1=''.join([Seq1,lines[j][0]])
			letter_seq[k1]=Seq1
			letter_seq[k1+'^']=reverse(complementary(Seq1))
		new_Seq=''
		new_letter=[]
		for k1 in letter_new:
			if not k1=='^':
				new_letter.append(k1)
			else:
				new_letter[-1]+=k1
		for k1 in new_letter:
			new_Seq+=letter_seq[k1]
			new_Seq+=insert_read_decide(bp_list)
		return new_Seq

def Ref_Ref_Produce(Chromo,bp_list,Ref_Seq_File):
	start=int(bp_list[0])
	end=int(bp_list[-1])
	new1_ref=''
	fin=os.popen(r'''samtools faidx %s %s:%d-%d'''%(Ref_Seq_File, Chromo, start,end))
	fin.readline().strip().split()
	for line in fin:
		pin=line.strip().split()
		new1_ref+=pin[0]
	fin.close()
	return new1_ref

def reverse(seq):
	seq2=[]
	for i in seq[::-1]:
		seq2.append(i)
	return ''.join(seq2)		

def complementary(seq):
	seq2=[]
	for i in seq:
		if i in 'ATGCN':
			seq2.append('ATGCN'['TACGN'.index(i)])
		elif i in 'atgcn':
			seq2.append('atgcn'['tacgn'.index(i)])
	return ''.join(seq2)		

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

def fasta_homo_write(fasta_out):
	fo=open(fasta_out,'w')
	print fasta_out
	for k1 in chromos:
		print >>fo, '>'+k1
		new1_ref=''
		rec1_start=0
		for k2 in sorted(order_SV_Pos[k1].keys()):
			rec1_start+=1
			k3=order_SV_Pos[k1][k2]
			start=int(k3[0][0][1])
			end=int(k3[0][0][-1])
			new1_ref+=Ref_Ref_Produce(k1,[rec1_start,start-1],ref)
			new1_ref+=Ref_Alt_Produce(chromos,k3[0][0],k3[0][1][0],ref)
			rec1_start=end
		rec1_start+=1
		new1_ref+=Ref_Ref_Produce(k1,[rec1_start,chromo_length[k1]],ref)
		new1_seq=[]
		for k1 in range(len(new1_ref)/60):
			new1_seq.append(new1_ref[k1*60:(k1+1)*60])
		new1_seq.append(new1_ref[(k1+1)*60:])
		for k1 in new1_seq:
			if not k1=='':
				print >>fo, k1
	fo.close()

def fasta_het_write_a(fasta_out):
	fo1=open(fasta_out.replace('.het.fa','.het1.fa'),'w')
	#fo2=open(fasta_out.replace('.het.fa','.het2.fa'),'w')
	fo1.close()
	#fo2.close()
	print fasta_out.replace('.het.fa','.het1.fa')
	#print fasta_out.replace('.het.fa','.het2.fa')
	for k1 in chromos:
		fo1=open(fasta_out.replace('.het.fa','.het1.fa'),'a')
		#fo2=open(fasta_out.replace('.het.fa','.het2.fa'),'a')
		print >>fo1, '>'+k1
		#print >>fo2, '>'+k1
		new1_ref=''
		rec1_start=0
		#new2_ref=''
		#rec2_start=0
		for k2 in sorted(order_SV_Pos[k1].keys()):
			print [k1,k2]
			rec1_start+=1
			k3=order_SV_Pos[k1][k2]
			start=int(k3[0][0][1])
			end=int(k3[0][0][-1])
			new1_ref+=Ref_Ref_Produce(k1,[rec1_start,start-1],ref)
			if not k3[0][1][0]==k3[0][1][2]:
				new1_ref+=Ref_Alt_Produce(chromos,k3[0][0],k3[0][1][0],ref)
			else:
				new1_ref+=Ref_Ref_Produce(k1,[start,end],ref)
			rec1_start=end
			#rec2_start+=1
			#new2_ref+=Ref_Ref_Produce(k1,[rec2_start,start-1],ref)
			#if not k3[0][1][1]==k3[0][1][2]:
			#	new2_ref+=Ref_Alt_Produce(chromos,k3[0][0],k3[0][1][1],ref)
			#else:
			#	new2_ref+=Ref_Ref_Produce(k1,[start,end],ref)
			#rec2_start=end
		rec1_start+=1
		#rec2_start+=1
		new1_ref+=Ref_Ref_Produce(k1,[rec1_start,chromo_length[k1]],ref)
		new1_seq=[]
		for ka1 in range(len(new1_ref)/60):
			new1_seq.append(new1_ref[ka1*60:(ka1+1)*60])
		new1_seq.append(new1_ref[(ka1+1)*60:])
		for ka1 in new1_seq:
			if not ka1=='':
				print >>fo1, ka1
		#new2_ref+=Ref_Ref_Produce(k1,[rec2_start,chromo_length[k1]],ref)
		#new2_seq=[]
		#for ka1 in range(len(new2_ref)/60):
		#	new2_seq.append(new2_ref[ka1*60:(ka1+1)*60])
		#new2_seq.append(new2_ref[(ka1+1)*60:])
		#for ka1 in new2_seq:
		#	if not ka1=='':
		#		print >>fo2, ka1
		fo1.close()
		#fo2.close()

def fasta_het_write_b(fasta_out):
	#fo1=open(fasta_out.replace('.het.fa','.het1.fa'),'w')
	fo2=open(fasta_out.replace('.het.fa','.het2.fa'),'w')
	#fo1.close()
	fo2.close()
	#print fasta_out.replace('.het.fa','.het1.fa')
	print fasta_out.replace('.het.fa','.het2.fa')
	for k1 in chromos:
		#fo1=open(fasta_out.replace('.het.fa','.het1.fa'),'a')
		fo2=open(fasta_out.replace('.het.fa','.het2.fa'),'a')
		#print >>fo1, '>'+k1
		print >>fo2, '>'+k1
		#new1_ref=''
		#rec1_start=0
		new2_ref=''
		rec2_start=0
		for k2 in sorted(order_SV_Pos[k1].keys()):
			print [k1,k2]
			k3=order_SV_Pos[k1][k2]
			start=int(k3[0][0][1])
			end=int(k3[0][0][-1])
			#rec1_start+=1
			#new1_ref+=Ref_Ref_Produce(k1,[rec1_start,start-1],ref)
			#if not k3[0][1][0]==k3[0][1][2]:
			#	new1_ref+=Ref_Alt_Produce(chromos,k3[0][0],k3[0][1][0],ref)
			#else:
			#	new1_ref+=Ref_Ref_Produce(k1,[start,end],ref)
			#rec1_start=end
			rec2_start+=1
			new2_ref+=Ref_Ref_Produce(k1,[rec2_start,start-1],ref)
			if not k3[0][1][1]==k3[0][1][2]:
				new2_ref+=Ref_Alt_Produce(chromos,k3[0][0],k3[0][1][1],ref)
			else:
				new2_ref+=Ref_Ref_Produce(k1,[start,end],ref)
			rec2_start=end
		#rec1_start+=1
		rec2_start+=1
		#new1_ref+=Ref_Ref_Produce(k1,[rec1_start,chromo_length[k1]],ref)
		#new1_seq=[]
		#for ka1 in range(len(new1_ref)/60):
		#	new1_seq.append(new1_ref[ka1*60:(ka1+1)*60])
		#new1_seq.append(new1_ref[(ka1+1)*60:])
		#for ka1 in new1_seq:
		#	if not ka1=='':
		#		print >>fo1, ka1
		new2_ref+=Ref_Ref_Produce(k1,[rec2_start,chromo_length[k1]],ref)
		new2_seq=[]
		for ka1 in range(len(new2_ref)/60):
			new2_seq.append(new2_ref[ka1*60:(ka1+1)*60])
		new2_seq.append(new2_ref[(ka1+1)*60:])
		for ka1 in new2_seq:
			if not ka1=='':
				print >>fo2, ka1
		#fo1.close()
		fo2.close()

def Sample_info_ReadIn(Sam_File):
	fi=open(Sam_File)
	for line in fi:
		pin=line.strip().split()
		if not pin==[]:
			if not pin[0] in sv_hash.keys():
				sv_hash[pin[0]]=[]
				sv_hash[pin[0]].append([int(i) for i in pin[1:]])
				sv_hash[pin[0]][-1][0]=int(sv_hash[pin[0]][-1][0]*1.25)
			else:
				sv_hash[pin[0]].append([int(i) for i in pin[1:]])
				sv_hash[pin[0]][-1][0]=int(sv_hash[pin[0]][-1][0]*1.25)
	fi.close()

def sv_total_num_calcu():
	sv_total_num=0
	for k1 in del_stat:
		sv_total_num+=k1[0]
	for k1 in dup_stat:
		sv_total_num+=k1[0]
	for k1 in inv_stat:
		sv_total_num+=k1[0]
	for k1 in tra_stat:
		sv_total_num+=k1[0]
	for k1 in dup2_stat:
		sv_total_num+=k1[0]		
	return sv_total_num

def pick_random_seqs(ref,sv_total_num,chromo_length):
	#12% of all SVs have micro insrts at both /either ends
	#double number of seqs would be randomly picked from genome as long micro-insertions
	num_micro_ins_over20bp=float(sv_total_num)*0.12*2
	genome_length=0
	chromos_num_regions={}
	chrom_seqs={}
	for x in chromo_length.keys():
		if not 'GL' in x and not x in ['X','Y','MT']:
			genome_length+=chromo_length[x]
	for x in chromo_length.keys():
		if not 'GL' in x and not x in ['X','Y','MT']:
			chromos_num_regions[x]=float(chromo_length[x])/float(genome_length)*num_micro_ins_over20bp
	for x in chromos_num_regions.keys():
		chrom_seqs[x]=[]
		int_num=int(round(chromos_num_regions[x]))
		seq_pick=random.sample(range(10000,chromo_length[x]-10000),int_num)
		for y in sorted(seq_pick):
			length_pick=random.sample(range(20,50),1)[0]
			seqs=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,x,y,y+length_pick))
			seqs.readline()
			test=seqs.readline().strip()
			if not 'NNNNNNNN' in test:
				chrom_seqs[x].append(test)
			seqs.close()
	return chrom_seqs

def produce_random_seqs(length):
	out=[]
	for x in range(length):
		out.append(random.choice(['A','T','G','C']))
	return ''.join(out)

opts,args=getopt.getopt(sys.argv[1:],'i:o:',['ref='])
dict_opts=dict(opts)
Sam_File=dict_opts['-i']
sv_hash={}
Sample_info_ReadIn(Sam_File)
del_stat=sv_stat_calcu(sv_hash,'DEL')
dup_stat=sv_stat_calcu(sv_hash,'DUP:TANDOM')
dup2_stat=sv_stat_calcu(sv_hash,'DUP')
dup3_stat=[]
for i in dup2_stat:
	dup3_stat.append([i[0]]+[j+1000 for j in i[1:]])

dup2_stat=dup3_stat
inv_stat=sv_stat_calcu(sv_hash,'INV')
tra_stat=sv_stat_calcu(sv_hash,'TRA')
sv_total_num=sv_total_num_calcu()
del_size=sv_size_pick(del_stat)
dup_size=sv_size_pick(dup_stat)
dup2_size=sv_size_pick(dup2_stat)
inv_size=sv_size_pick(inv_stat)
tra_size=sv_size_pick(tra_stat)
refs=dict_opts['--ref']
ref=refs
if not os.path.isfile(refs):
	print 'Wrong reference genome !'

if not os.path.isfile(refs+'.fai'):
	print 'reference genome not indexed !'

chromos_TOTAL=chromos_readin(refs)
genome_length=chromos_TOTAL[0]
chromos=chromos_TOTAL[1]
chromo_num_region=chromos_TOTAL[2]
chromo_length=chromos_TOTAL[3]
sv_hash={}
sv_hash_add(del_size,'DEL')
sv_hash_add(dup_size,'DUP:TANDOM')
sv_hash_add(dup2_size,'DUP')
sv_hash_add(inv_size,'INV')
sv_hash_add(tra_size,'TRA')
SV_region=sv_region_pick()
SV_region_filter=[]
for x in SV_region:
	if x[-1]=='DUP' and x[2]-x[1]<1100: continue
	else:
		SV_region_filter.append(x)

SV_region=SV_region_filter
sv_homo_info={}
sv_homo_initial()
sv_homo_produce()
sv_het_info={}
sv_het_produce()
for y in range(len(sv_het_info['DUP'])):
	x=sv_het_info['DUP'][y]
	if x[2]-x[1]<2000:
		z=random.choice([x[1]+1000,x[2]-1000])
	else:
		z=random.choice(range(x[1]+800,x[1]+1200)+range(x[2]-1200,x[2]-800))
	sv_het_info['DUP'][y]=x[:2]+[z]+x[2:]

sv_rec_het_produce()
sv_info={}
sv_info_rewrite(sv_het_info)
dup1={}
inv1={}
del1={}
tra1={}
sv_rec_2(sv_info)
sv_out={}
hash_reorder()
vcf_out=dict_opts['-o']+'.vcf'
write_VCF_header(vcf_out)
write_VCF_main(vcf_out)
fasta_out=dict_opts['-o']+'.het.fa'
#produce fasta file containing all sv file for homo svs
order_SV_Pos={}
order_SV_Het_write(sv_info)
seq_ins_pools=pick_random_seqs(ref,sv_total_num,chromo_length)
fasta_het_write_a(fasta_out)
fasta_het_write_b(fasta_out)
os.system(r'''samtools faidx %s'''%(fasta_out.replace('.het.fa','.het1.fa')))
os.system(r'''samtools faidx %s'''%(fasta_out.replace('.het.fa','.het2.fa')))


