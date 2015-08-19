#!/usr/bin/python

#!python
#command='Comp.Compare.Delly.py --ref_genome /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.fasta -p /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/Delly/ --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/comp.SV.rec'
#sys.argv=command.split()
import os
import sys
import getopt
import numpy
import re
import random
import pickle
import time
import datetime
opts,args=getopt.getopt(sys.argv[1:],'p:',['ref=','ref_genome='])
dict_opts=dict(opts)
def end_cordi_calcu(pin):
	#info=pin
	chromo=pin[0]
	start=int(pin[1])
	end=0
	for x in pin[7].split(';'):
		if 'END' in x.split('='):
			end=int(x.split('=')[1])
	if not end==0:
		return [chromo,start,end]
	else:
		return 'Error!'

def genotype_calcu(pin):
	gt_index=pin[8].split(':').index('GT')
	geno=pin[9].split(':')[gt_index]
	out=0
	if '/' in geno:
		if geno.split('/')==['0','0']:
			out='Error!'
		elif sorted(geno.split('/'))==['0','1']:
			out='het'
		elif sorted(geno.split('/'))==['1','1']:
			out='homo'
	elif '|' in geno:
		if geno.split('|')==['0','0']:
			out='Error!'
		elif sorted(geno.split('|'))==['0','1']:
			out='het'
		elif sorted(geno.split('|'))==['1','1']:
			out='homo'
	return out

def geno_produce(pin):
	genotype=genotype_calcu(pin)
	if genotype=='homo':
		return key_homo
	elif genotype=='het':
		return key_het
	else:
		return 'Error'

def block_match(info_list):
	#eg of info_list=[['1', 26493213, 26493320, 26493631, '0.0', 'ab/ab', 'ab/b^ab'], ['1', 26493213, 26493320, 26493632, 9.9245387286030002, 'ab/ab', 'ab/b^ab']]
	#format: [ref_info_list,in_info_list]
	ref_bps=info_list[0][:-3]
	in_bps=info_list[1][:-3]
	ref_in_bps=[]
	rec1=1
	rec2=1
	while True:
		if rec1==len(ref_bps) and rec2==len(in_bps): break
		if rec1==len(ref_bps):
			while rec2<len(in_bps):
				ref_in_bps.append([0,in_bps[rec2]])
				rec2+=1
			continue
		if rec2==len(in_bps):
			while rec1<len(ref_bps):
				ref_in_bps.append([ref_bps[rec1],0])
				rec1+=1
			continue
		if abs(ref_bps[rec1]-in_bps[rec2])<51:
			ref_in_bps.append([ref_bps[rec1],in_bps[rec2]])
			rec1+=1
			rec2+=1
		elif ref_bps[rec1]-in_bps[rec2]>50:
			ref_in_bps.append([0,in_bps[rec2]])
			rec2+=1
		elif in_bps[rec2]-ref_bps[rec1]>50:
			ref_in_bps.append([ref_bps[rec1],0])
			rec1+=1
	ref_new=[ref_bps[0]]+[i[0] for i in ref_in_bps]
	in_new=[in_bps[0]]+[i[1] for i in ref_in_bps]
	ref_old_new=let_compare(ref_bps,ref_new)
	in_old_new=let_compare(in_bps,in_new)
	out=[[],[],[],[],[]]
	for x in ref_in_bps:
		if not x[0]==0:
			out[0].append(x[0])
		else:
			out[0].append(x[1])
	out[1]=''.join([''.join(ref_old_new[x]) for x in info_list[0][-2]])
	out[2]=''.join([''.join(ref_old_new[x]) for x in info_list[0][-1]])
	out[3]=''.join([''.join(in_old_new[x]) for x in info_list[1][-2]])
	out[4]=''.join([''.join(in_old_new[x]) for x in info_list[1][-1]])
	return out

def modify_comp_new(comp_new):
	for x in comp_new.keys():
		for y in comp_new[x].keys():
			for z in comp_new[x][y]:
				z[2]='/'.join(sorted(z[2].split('/')))
				z[4]='/'.join(sorted(z[4].split('/')))

def write_comp_new(comp_new,fout):
	fo=open(fout,'w')
	for x in comp_new.keys():
		for y in comp_new[x].keys():
			for z in comp_new[x][y]:
				print >>fo,  ' '.join([str(i) for i in [x]+z[0]+z[1:]])
	fo.close()

def chromos_name_readin(ref_file):
	ref_index=ref_file+'.fai'
	ref=ref_file
	chromos=[]
	fin=open(ref_index)
	for line in fin:
		pin=line.strip().split()
		chromos.append(pin[0])
	fin.close()
	return chromos

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

def hash_to_list(in_hash):
	ref_list={}
	for k1 in in_hash.keys():
		ref_list[k1]=[]
		for k2 in sorted(in_hash[k1].keys()):
			for k3 in in_hash[k1][k2]:
				ref_list[k1].append(k3)
	return ref_list

def bp_check(list1,list2):
	#eg: list1=['1', 3615657, 3638623]; list2=['1', 3615658, 3638624, 3751340];
	flag=0
	rec=[]
	for x in list1[1:]:
		for y in list2[1:]:
			if abs(x-y)<50:
				flag+=1
				rec.append(list2.index(y))
	return [flag,rec]

def compare_list(in_list,ref_list):
	out={}
	for k1 in in_list.keys():
		if k1 in ref_list.keys():
			print k1
			out[k1]={}
			for k2a in ref_list[k1]:
				for k2b in in_list[k1]:
					if k2b[-4]<k2a[1]: continue
					elif k2b[1]>k2a[-4]:continue
					else:
						if bp_check(k2a[:-3],k2b[:-3])>0:
							if not '_'.join([str(i) for i in k2a[:-3]]) in out[k1].keys():
								out[k1]['_'.join([str(i) for i in k2a[:-3]])]=[]
							out[k1]['_'.join([str(i) for i in k2a[:-3]])].append([k2a,k2b])
	return out	

def chrom_num_check(list):
	#eg of list: ['1', '770745', '770845', '6', '170920617', '170920717', 9.3618710041200011, 'ab/ab', 'aa/b']
	flag=0
	for x in list:
		if x in chromos:
			flag+=1
	return flag

def in_list_chrom_filter(in_list):
	out={}
	for x in in_list.keys():
		out[x]=[]
		for y in in_list[x]:
			if chrom_num_check(y)==1:
				out[x].append(y)
	return out

def list_transfor(in_list):
	for k1 in in_list.keys():
		for k2 in in_list[k1]:
			k2[1:-3]=[int(i) for i in k2[1:-3]]

def tra_csv_info_add(k1,k2,reorder_hash,tra1):
	for k3a in reorder_hash[k1][k2]:
		k3=k3a
		SV_ID='_'.join([str(i) for i in k3[:-1]])
		tra1[SV_ID]={}
		k2a=k2.split('/')[0]
		k2b=k2.split('/')[1]
		bp_hash={}
		block_rec=0
		block_hash=[]
		for a3 in k3[:-1]:
		    if a3 in chromos:
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
		bp_hash['left']=[bp_hash[k1[0]][0],bp_hash[k1[0]][1],bp_hash[k1[0]][2]]
		bp_hash['right']=[bp_hash[k1[-1]][0],bp_hash[k1[-1]][3],bp_hash[k1[-1]][4]]
		ref_allele={}
		for a3 in bp_hash.keys():
			ref_allele[a3]=[bp_hash[a3][0]]
			for a4 in bp_hash[a3][1:]:
				ref_allele[a3].append('N')
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

def simple_flag_SA(k1,k2):
	print [k1,k2]
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

def let_compare(old_bps,new_bps):
	in_old_new={}
	rec=0
	in_old_let=bp_to_let([[str(i) for i in old_bps]])
	in_new_let=bp_to_let([[str(i) for i in new_bps]])
	for x in in_new_let.split('/')[0]:
		temp=[new_bps[ord(x)-96],new_bps[ord(x)-95]]
		if temp==[0,0]: rec=rec
		elif temp[0]==0: rec=rec
		else: rec+=1
		if rec>0 and rec<len(old_bps)-1:
			in_old_new[x]=chr(96+rec)
		else:
			in_old_new[x]='0'
	out={}
	for x in in_old_new.keys():
		if not in_old_new[x] in out.keys():
			out[in_old_new[x]]=[]
		out[in_old_new[x]].append(x)
	for x in out.keys():
		out[x].sort()
	out['^']='^'
	out['/']='/'
	return out

def block_match(info_list):
	#eg of info_list=[['1', 26493213, 26493320, 26493631, '0.0', 'ab/ab', 'ab/b^ab'], ['1', 26493213, 26493320, 26493632, 9.9245387286030002, 'ab/ab', 'ab/b^ab']]
	#format: [ref_info_list,in_info_list]
	ref_bps=info_list[0][:-3]
	in_bps=info_list[1][:-3]
	ref_in_bps=[]
	rec1=1
	rec2=1
	while True:
		if rec1==len(ref_bps) and rec2==len(in_bps): break
		if rec1==len(ref_bps):
			while rec2<len(in_bps):
				ref_in_bps.append([0,in_bps[rec2]])
				rec2+=1
			continue
		if rec2==len(in_bps):
			while rec1<len(ref_bps):
				ref_in_bps.append([ref_bps[rec1],0])
				rec1+=1
			continue
		if abs(ref_bps[rec1]-in_bps[rec2])<51:
			ref_in_bps.append([ref_bps[rec1],in_bps[rec2]])
			rec1+=1
			rec2+=1
		elif ref_bps[rec1]-in_bps[rec2]>50:
			ref_in_bps.append([0,in_bps[rec2]])
			rec2+=1
		elif in_bps[rec2]-ref_bps[rec1]>50:
			ref_in_bps.append([ref_bps[rec1],0])
			rec1+=1
	ref_new=[ref_bps[0]]+[i[0] for i in ref_in_bps]
	in_new=[in_bps[0]]+[i[1] for i in ref_in_bps]
	ref_old_new=let_compare(ref_bps,ref_new)
	in_old_new=let_compare(in_bps,in_new)
	out=[[],[],[],[],[]]
	for x in ref_in_bps:
		if not x[0]==0:
			out[0].append(x[0])
		else:
			out[0].append(x[1])
	out[1]=''.join([''.join(ref_old_new[x]) for x in info_list[0][-2]])
	out[2]=''.join([''.join(ref_old_new[x]) for x in info_list[0][-1]])
	out[3]=''.join([''.join(in_old_new[x]) for x in info_list[1][-2]])
	out[4]=''.join([''.join(in_old_new[x]) for x in info_list[1][-1]])
	return out

def chromos_name_readin(ref_file):
	ref_index=ref_file+'.fai'
	ref=ref_file
	chromos=[]
	fin=open(ref_index)
	for line in fin:
		pin=line.strip().split()
		chromos.append(pin[0])
	fin.close()
	return chromos

path_in=dict_opts['-p']
if not path_in[-1]=='/':
	path_in+='/'

ref=dict_opts['--ref_genome']
chromos=chromos_name_readin(ref)
file_hash={}
for k1 in os.listdir(path_in):
	if k1.split('.')[-1]=='vcf':
		if 'DEL' in k1:
			x1=k1.replace('DEL.vcf','')
			if not x1 in file_hash.keys():
				file_hash[x1]=[]
			file_hash[x1].append(k1)
		if 'DUP' in k1:
			x1=k1.replace('DUP.vcf','')
			if not x1 in file_hash.keys():
				file_hash[x1]=[]
			file_hash[x1].append(k1)
		if 'INV' in k1:
			x1=k1.replace('INV.vcf','')
			if not x1 in file_hash.keys():
				file_hash[x1]=[]
			file_hash[x1].append(k1)
		if 'TRA' in k1:
			x1=k1.replace('TRA.vcf','')
			if not x1 in file_hash.keys():
				file_hash[x1]=[]
			file_hash[x1].append(k1)

ref_hash={}
fin=open(dict_opts['--ref'])
for line in fin:
	pin=line.strip().split()
	if not pin[0] in ref_hash.keys():
		ref_hash[pin[0]]={}
	if not int(pin[1]) in ref_hash[pin[0]].keys():
		ref_hash[pin[0]][int(pin[1])]=[]
	if not [pin[0]]+[int(i) for i in pin[1:4]] in ref_hash[pin[0]][int(pin[1])]:
		ref_hash[pin[0]][int(pin[1])].append(pin)

fin.close()
for k1 in file_hash.keys():
	delly_hash={}
	for k2 in file_hash[k1]:
		if not 'TRA' in k2:
			if 'DEL' in k2:
				keya='a/a'
				key_het='a/'
				key_homo='/'
			elif 'INV' in k2:
				keya='a/a'
				key_het='a/a^'
				key_homo='a^/a^'
			elif 'DUP' in k2:
				keya='a/a'			
				key_het='aa/aa'
				key_homo='a/aa'
			fin=open(path_in+k2)
			for line in fin:
				pin=line.strip().split()
				if not pin[0][0]=='#': 
					if pin[6]=='PASS':
						pos=end_cordi_calcu(pin)
						key_geno=geno_produce(pin)
						if not key_geno=='Error':
							if not pos[0] in delly_hash.keys():
								delly_hash[pos[0]]={}
							if not pos[1] in delly_hash[pos[0]].keys():
								delly_hash[pos[0]][pos[1]]=[]
							delly_hash[pos[0]][pos[1]].append(pos+[0,keya,key_geno])
			fin.close()
	in_list=hash_to_list(delly_hash)
	ref_list=hash_to_list(ref_hash)
	in_list=in_list_chrom_filter(in_list)
	list_transfor(in_list)
	list_transfor(ref_list)
	comp_all=compare_list(in_list,ref_list)
	comp_new={}
	for kx in comp_all.keys():
		comp_new[kx]={}
		for k2 in comp_all[kx].keys():
			comp_new[kx][k2]=[]
			for k3 in comp_all[kx][k2]:
				if bp_check(k3[0][:-3],k3[1][:-3])[0]>0:
					if not block_match(k3) in comp_new[kx][k2]:
						comp_new[kx][k2].append(block_match(k3))
	modify_comp_new(comp_new)
	write_comp_new(comp_new,dict_opts['-p']+k1+'comp.SV.compare.report')

