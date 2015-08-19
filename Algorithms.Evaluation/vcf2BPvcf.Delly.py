#!/usr/bin/python

#!python
#command='python SVelter.vcf2BPvcf.py --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.fasta -p /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/Delly/ --prefix Delly_comp_RD_10_'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import time
import datetime
opts,args=getopt.getopt(sys.argv[1:],'p:i:o:h:',['ref=','prefix=','appdix='])
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
		return 'Error'

def genotype_calcu(pin):
	gt_index=pin[8].split(':').index('GT')
	geno=pin[9].split(':')[gt_index]
	out=0
	if '/' in geno:
		if geno.split('/')==['0','0']:
			out='Error'
		elif sorted(geno.split('/'))==['0','1']:
			out='het'
		elif sorted(geno.split('/'))==['1','1']:
			out='homo'
	elif '|' in geno:
		if geno.split('|')==['0','0']:
			out='Error'
		elif sorted(geno.split('|'))==['0','1']:
			out='het'
		elif sorted(geno.split('|'))==['1','1']:
			out='homo'
	return out

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

def let_to_bp_list(letter_list,bp_hash):
    out=[]
    if letter_list==[]:
        return []
    else:
        for k1 in letter_list:
            out.append([])
            for k2 in k1:
                if not k2=='^':
                    if out[-1]==[]:
                        out[-1]+=[bp_hash[k2][0],bp_hash[k2][2],bp_hash[k2][4]]
                    else:
                        if bp_hash[k2][0]==out[-1][0]:
                            out[-1]+=[bp_hash[k2][2],bp_hash[k2][4]]
                        else:
                            out.append([bp_hash[k2][0],bp_hash[k2][2],bp_hash[k2][4]])
        out2=[]
        for k1 in out:
            out2.append([k1[0],k1[1],k1[-1]])
        return out2

def let_to_lets(letter_list):
    out=[]
    for k1 in letter_list:
        if k1==[]:
            out.append(k1)
        else:
            k2=[]
            if not '^' in k1[0]:
                for k3 in k1:
                    if k2==[]:
                        k2.append(k3)
                    else:
                        if ord(k3)-ord(k2[-1][-1])==1:
                            k2[-1]+=k3
                        else:
                            k2.append(k3)
            else:
                for k3 in k1:
                    if k2==[]:
                        k2.append(k3)
                    else:
                        if -ord(k3[0])+ord(k2[-1][0])==1:
                            k2[-1]=k3[0]+k2[-1]
                        else:
                            k2.append(k3)
            out.append(k2)
    return out

def read_in_structures(filein):
    fin=open(filein)
    while True:
        pin1=fin.readline().strip().split()
        if not pin1: break
        if pin1[0]=='Total': break
        pin2=fin.readline().strip().split()
        pin3=fin.readline().strip().split()
        pin4=fin.readline().strip().split()
        pin5=fin.readline().strip().split()
        if pin3[0]=='Theoretical' and pin4[0]=='Current' and pin5[0]=='Time':
            #if float(pin4[-1])-float(pin3[-1])>cutoff:
                let1=bp_to_let([pin1])
                if not let1==0:
	                let2='/'.join(sorted(pin2[0].split('/')))
	                if not let1 in sv_info.keys():
	                    sv_info[let1]={}
	                if not let2 in sv_info[let1].keys():
	                    sv_info[let1][let2]=[]
	                if not pin1 in sv_info[let1][let2]:
	                    sv_info[let1][let2].append(pin1+[float(pin4[-1])-float(pin3[-1])])
    fin.close()

def let_recombine(let_list):
	out=[]
	if not let_list==[]:
		out=[let_list[0]]
		for i in let_list[1:]:
			if ord(i)-ord(out[-1][-1])==1:
				out[-1]+=i
			else:
				out.append(i)
	return out

def inv_flag(k1,k2):
	if '^' in k2:
		return 1
	else:
		return 0

def inv_flag_SA(k1,k2):
	#SA:single allele
	out=0
	if '^' in k2:
		if k2.replace('^','')==k1:
			out+=1
	return out

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

def dup_flag(k1,k2):
	test=0
	for i in k2.split('/')[0]:
		if k2.split('/')[0].count(i)>1:
			test=1
	for i in k2.split('/')[1]:
		if k2.split('/')[1].count(i)>1:
			test=1
	return test

def dup_flag_SA(k1,k2):
	test=0
	for i in k2:
		if k2.count(i)>1:
			if k2.replace(''.join([i for j in range(k2.count(i))]),i)==k1:
				test=1
	return test

def tra_flag(k1,k2):
	test=0
	k2a=[i for i in k2.split('/')[0]]
	k2b=[i for i in k2.split('/')[1]]
	for i in range(len(k2a)-1):
		if ord(k2a[i+1])-ord(k2a[i])<0:
			test=1
	for i in range(len(k2b)-1):
		if ord(k2b[i+1])-ord(k2b[i])<0:
			test=1
	return test

def tra_flag_SA(k1,k2):
	test=0
	k2a=[i for i in k2]
	for i in range(len(k2a)-1):
		if ord(k2a[i+1])-ord(k2a[i])<0:
			test=1
	return test

def del_flag(k1,k2):
	if inv_flag(k1,k2)+dup_flag(k1,k2)+tra_flag(k1,k2)==0:
		return 1
	else:
		return 0

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

def allele_DUP_decide(x1):
	flag_Dup=0
	for x2 in x1:
		if x1.count(x2)>1:
			flag_Dup+=1
	if flag_Dup==0:
		return 'NO_DUP'
	else:
		return 'DUP'

def allele_TRA_decide(x1):
	flag_TRA=0
	if len(x1)<2:
		return 'NO_TRA'
	else:
		for x2 in range(len(x1)-1):
			if ord(x1[x2+1])-ord(x1[x2])<0:
				flag_TRA+=1
		if flag_TRA==0:
			return 'NO_TRA'
		else:
			return 'TRA'

def simple_DEL_decide(k1,k2):
	out='Not_Simple_DEL' 
	if not k2==k1:
		if not '^' in k2:
			flag_Dup=0
			for x1 in k2.split('/'):
				if allele_DUP_decide(x1)=='DUP':
					flag_Dup+=1
			if flag_Dup==0:
				flag_TRA=0
				for x1 in k2.split('/'):
					if allele_TRA_decide(x1)=='TRA':
						flag_TRA+=1
				if flag_TRA==0:
					out='Simple_DEL'
	return out

def simple_DEL_type(k1,k2):
	if simple_DEL_decide(k1,k2)=='Simple_DEL':
		del_part=[[],[]]
		for x1 in k1.split('/')[0]:
			if not x1 in k2.split('/')[0]:
				del_part[0].append(x1)
		for x1 in k1.split('/')[1]:
			if not x1 in k2.split('/')[1]:
				del_part[1].append(x1)
		return del_part
	else:
		return 'Error'

def simple_INV_decide(k1,k2):
	out='Not_Simple_INV' 
	if not k2==k1:
		if '^' in k2:
			k2_new=letter_recombine(k2)
			if '/'.join([''.join([i.replace('^','') for i in k2_new[0]]),''.join([i.replace('^','') for i in k2_new[1]])])==k1:
				out='Simple_INV'
	return out

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

def simple_INV_type(k1,k2):
	k2_new=letter_separate(k2)
	if simple_INV_decide(k1,k2)=='Simple_INV':
		out=[[],[]]
		for x1 in k2_new[0]:
			if '^' in x1:
				out[0].append(x1)
		for x1 in k2_new[1]:
			if '^' in x1:
				out[1].append(x1)
		out2=[[x1.replace('^','') for x1 in out[0]],[x1.replace('^','') for x1 in out[1]]]
		return out2
	else:
		return 'Error'

def simple_DUP_decide(k1,k2):
	out='Not_Simple_DUP' 
	if not k2==k1:
		if not '^' in k2:
			if simple_DEL_decide(k1,k2)=='Not_Simple_DEL':
				k2_new=letter_separate(k2)
				k3=[[],[]]
				for x1 in k2_new[0]:
					if k3[0]==[]:
						k3[0].append(x1)
					elif not x1==k3[0][-1]:
						k3[0].append(x1)
				for x1 in k2_new[1]:
					if k3[1]==[]:
						k3[1].append(x1)
					elif not x1==k3[1][-1]:
						k3[1].append(x1)
				if k3==letter_separate(k1):
					out='Simple_DUP'
	return out

def simple_DUP_type(k1,k2):
	if simple_DUP_decide(k1,k2)=='Simple_DUP':
		out=[[],[]]
		for x1 in k2.split('/')[0]:
			if k2.split('/')[0].count(x1)>1:
				if not x1 in out[0]:
					out[0].append(x1)
		for x1 in k2.split('/')[1]:
			if k2.split('/')[1].count(x1)>1:
				if not x1 in out[1]:
					out[1].append(x1)	
		return out	
	else:
		return 'Error'

def simple_TRA_decide(k1,k2):
	out='Not_Simple_TRA' 
	if '/'.join([''.join(sorted(k2.split('/')[0])),''.join(sorted(k2.split('/')[1]))])==k1:
		out='simple_TRA'
	return out

def simple_flag_SA(k1,k2):
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
			if ord(temp3[i+1])-ord(temp3[i])<0:
				outtra=1
		if not temp3==k1:
			temp4=[]
			for i in temp3:
				if temp3.count(i)>2:
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
					if k3.count(ia)+k3.count(ib)>2:
						outdup2.append([i,k3.count(ia)])
						k3=k3.replace(ia,'')
						k3=k3.replace(ib,'')
		else:
			outdup2=[]
		return [outdel,outinv,outdup2,outtra]

def sv_rec_3(sv_info):
	for k1ab in sv_info.keys():
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

def sv_rec(sv_info):
	for k1ab in sv_info.keys():
		for k2ab in sv_info[k1ab].keys():
			if not k2ab==k1ab:
				if del_flag(k1ab,k2ab)==1:
					#simple deletion
					delM=[]
					delP=[]
					for i in k1ab.split('/')[0]:
						if not i in k2ab.split('/')[0]:
							delM.append(i)
						if not i in k2ab.split('/')[1]:
							delP.append(i)
					for k3 in sv_info[k1ab][k2ab]:
						del_info_add(k3,[delM,delP])
					#dhom=[]
					#dhet=[[],[]]
					#for i in delM:
					#	if i in delP:
					#		dhom.append(i)
					#	else:
					#		dhet[0].append(i)
					#for i in delP:
					#	if i in delM:
					#		continue
					#	else:
					#		dhet[1].append(i)
					#dhom2=let_recombine(dhom)
					#dhet2=[let_recombine(dhet[0]),let_recombine(dhet[1])]
				else:
					if inv_flag(k1ab,k2ab)+dup_flag(k1ab,k2ab)==0:
						tra_info_add(k1ab,k2ab)
						if del_flag_SA(k1ab.split('/')[0],k2ab.split('/')[0])==1:
							delM=[]
							delP=[]
							for i in k1ab.split('/')[0]:
								if not i in k2ab.split('/')[0]:
									delM.append(i)
							for k3 in sv_info[k1ab][k2ab]:
								del_info_add(k3,[delM,delP])
						if del_flag_SA(k1ab.split('/')[1],k2ab.split('/')[1])==1:
							delM=[]
							delP=[]
							for i in k1ab.split('/')[0]:
								if not i in k2ab.split('/')[1]:
									delP.append(i)
							for k3 in sv_info[k1ab][k2ab]:
								del_info_add(k3,[delM,delP])
					else:
						k1aba=k1ab.split('/')[0]
						k2aba=k2ab.split('/')[0]
						k2abb=k2ab.split('/')[1]
						flaga=[]
						flagb=[]
						if del_flag_SA(k1aba,k2aba)==1:#simple del on one allele
							delM=[]
							delP=[]
							for i in k1ab.split('/')[0]:
								if not i in k2ab.split('/')[0]:
									delM.append(i)
							for k3 in sv_info[k1ab][k2ab]:
								del_info_add(k3,[delM,delP])
							flaga.append('del')
						if del_flag_SA(k1aba,k2abb)==1:#simple del on one allele
							delM=[]
							delP=[]
							for i in k1ab.split('/')[0]:
								if not i in k2ab.split('/')[1]:
									delP.append(i)
							for k3 in sv_info[k1ab][k2ab]:
								del_info_add(k3,[delM,delP])
							flagb.append('del')
						if dup_flag_SA(k1aba,k2aba)==1:#simple dup on one allele
							dupM=[]
							dupP=[]
							for i in k1aba:
								if k2aba.count(i)>1:
									dupM.append(i)
							for k3 in sv_info[k1ab][k2ab]:
								dup_info_add(k3,[dupM,dupP])
							flaga.append('dup')
						if dup_flag_SA(k1aba,k2abb)==1:#simple dup on one allele
							dupM=[]
							dupP=[]
							for i in k1aba:
								if k2abb.count(i)>1:
									dupP.append(i)
							for k3 in sv_info[k1ab][k2ab]:
								dup_info_add(k3,[dupM,dupP])
							flagb.append('dup')
						if inv_flag_SA(k1aba,k2aba)==1:#simple inv on one allele
							invM=[]
							invP=[]
							for i in range(len(k2aba)):
								if k2aba[i]=='^':
									invM.append(k2aba[i-1])
							for k3 in sv_info[k1ab][k2ab]:
								inv_info_add(k3,[invM,invP])
							flaga.append('inv')
						if inv_flag_SA(k1aba,k2abb)==1:#simple inv on one allele
							invM=[]
							invP=[]
							for i in range(len(k2abb)):
								if k2abb[i]=='^':
									invP.append(k2abb[i-1])
							for k3 in sv_info[k1ab][k2ab]:
								inv_info_add(k3,[invM,invP])
							flagb.append('inv')
						if flaga==[] and not k1aba==k2aba:
							csv1=simple_flag_SA(k1aba,k2aba)
							add_csv_info(csv1,1,k1ab,k2ab)
						if flagb==[] and not k1aba==k2abb:
							csv1=simple_flag_SA(k1aba,k2abb)
							add_csv_info(csv1,2,k1ab,k2ab)

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

def del_info_reorganize(k1,k2):
    #k1: original normal structure
    #k2: rearranged structure
    del_let=[[],[]]
    for k3 in k1.split('/')[0]:
        if not k3 in k2.split('/')[0]:
            del_let[0].append(k3)
    for k3 in k1.split('/')[1]:
        if not k3 in k2.split('/')[1]:
            del_let[1].append(k3)
    for k3 in sv_info[k1][k2]:
        del_bp=[]
        if not del_let[0]==[]:
            del_bp.append(bp_to_hash(k3,del_let[0]))
        else:
            del_bp.append([])
        if not del_let[1]==[]:
            del_bp.append(bp_to_hash(k3,del_let[1]))
        else:
            del_bp.append([])            
        if del_bp[0]==del_bp[1]:
            for k4 in del_bp[0]:
                if not k4[0] in del1.keys():
                    del1[k4[0]]=[]
                if not [int(k4[1]),int(k4[-1]),'hom'] in del1[k4[0]]:
                    del1[k4[0]].append([int(k4[1]),int(k4[-1]),'hom'])
        else:
            for k5 in del_bp:
                for k4 in k5:
                    if not k4[0] in del1.keys():
                        del1[k4[0]]=[]
                    if not [int(k4[1]),int(k4[-1]),'het'] in del1[k4[0]]:
                        del1[k4[0]].append([int(k4[1]),int(k4[-1]),'het'])

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
	    del1[k1[0]].append(k1[1:]+[tempc,k3[-1],'_'.join(k3[:-1]+['S'])])
	for k1 in tempb:
	    if not k1[0] in del1.keys():
	        del1[k1[0]]=[]
	    del1[k1[0]].append(k1[1:]+['hetb',k3[-1],'_'.join(k3[:-1]+['S'])])

def del_csv_info_add(k3,del_let):
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
	    del1[k1[0]].append(k1[1:]+[tempc,k3[-1],'_'.join(k3[:-1]+['C'])])
	for k1 in tempb:
	    if not k1[0] in del1.keys():
	        del1[k1[0]]=[]
	    del1[k1[0]].append(k1[1:]+['hetb',k3[-1],'_'.join(k3[:-1]+['C'])])

def dup_info_add(k3,dup_let):
    #dup_let=[k2i,k2j]
    for k2x in dup_let:
        for k4 in k2x:
            temp=bp_to_hash(k3[:-1],[i for i in k4])
            for k5 in temp:
                if not k5[0] in dup1.keys():
                    dup1[k5[0]]=[]
                dup1[k5[0]].append(k5[1:]+[k3[-1],'_'.join(k3[:-1])])

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
				    dup1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S']),k4[1]])

def dup_csv_info_2_add(k3,dup_let):
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
				    dup1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['C']),k4[1]])

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
                inv1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S'])])

def inv_csv_info_add(k3,inv_let):
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
                inv1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['C'])])

def tra_csv_info_add(k1,k2):
	for k3 in sv_info[k1][k2]:
		SV_ID='_'.join([str(i) for i in k3]+['C'])
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
				#ref_allele[a3].append(ref_base_readin(ref,bp_hash[a3][0],a4))
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

def tra_info_add_new(k1,k2):
	for k3 in sv_info[k1][k2]:
		k3=[str(i) for i in k3]
		SV_ID='_'.join(k3[:-1]+[k1,k2])
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
				#ref_allele[a3].append(ref_base_readin(ref,bp_hash[a3][0],a4))
				ref_allele[a3].append('N')
		if k2a=='':
		    a3='left'
		    a4='right'
		    l_chr=bp_hash[a3][0]
		    r_chr=bp_hash[a4][0]
		    tra1[SV_ID]['a']=[]
		    tra1[SV_ID]['a'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4][2]])
		    tra1[SV_ID]['a'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4][2])+'['])
		else:
			if not k2a==k1.split('/')[0]:
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
		if k2b=='':
		    a3='left'
		    a4='right'
		    l_chr=bp_hash[a3][0]
		    r_chr=bp_hash[a4][0]
		    tra1[SV_ID]['b']=[]
		    tra1[SV_ID]['b'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4][2]])
		    tra1[SV_ID]['b'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4][2])+'['])
		else:
			if not k2b==k1.split('/')[1]:
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

def hash_reorder_new():
	for ka1 in tra1.keys():
		ks1=ka1.split('_')[0]
		ks2='_'.join(ka1.split('_'))
		SV_Score=0.0
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
			if ka2=='b':
				GenoType='0|1'
			for ka3 in tra1[ka1][ka2]:
				ka_new=ka3[:2]+[ks2,ka3[2]]+ka3[3:]+[SV_Score,Pass_Sign,'SVTYPE=TRA','GT',GenoType]
				#print ka_new
				if not ka_new in sv_out[ks1][ks2]:
					sv_out[ks1][ks2].append(ka_new)

def write_VCF_header(output_file):
	fo=open(output_file,'w')
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
	print>>fo,'##ALT=<ID=DUP,Description="Duplication">'
	print>>fo,'##ALT=<ID=INV,Description="Inversion">'
	print>>fo,'##ALT=<ID=TRA,Description="Translocation">'
	print>>fo,'##ALT=<ID=INS,Description="Insertion">'
	print>>fo,'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
	print>>fo,'##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">'
	print>>fo,'##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">'
	print>>fo,'##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">'
	print>>fo,'\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',output_file.split('/')[-1].replace('.vcf','')])
	fo.close()

def ref_base_readin(ref,chromo,pos):
	fref=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,chromo,str(pos),str(pos)))
	tre=fref.readline().strip().split()
	REF_AL=fref.readline().strip().split()
	if not REF_AL==[]:
		return REF_AL[0]
	else:
		return 'N'

def write_VCF_main(output_file):
	fo=open(output_file,'a')
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
						#if k4[3]=='N':
						#	k4[3]=ref_base_readin(ref,k4[0],k4[1])
						if not k4[-1] in ['1|1','1/1','0|0','0/0']:
							print >>fo, '\t'.join([str(i) for i in k4])
						elif k4[-1] in ['1|1','1/1']:
							print >>fo, '\t'.join([str(i) for i in k4[:-1]+['1|0']])
							print >>fo, '\t'.join([str(i) for i in k4[:-1]+['0|1']])
						#print '\t'.join([str(i) for i in k4])
	fo.close()

def ROC_produce_files(ref_file,samp_file):
    ref_hash={}
    samp_hash={}
    out={}
    for i in chromos:
        ref_hash[i]=[]
        samp_hash[i]=[]
    fin=open(ref_file)
    for line in fin:
        pin=line.strip().split()
        ref_hash[pin[0]].append([int(i) for i in pin[1:3]])
    fin.close()
    fin=open(samp_file)
    for line in fin:
        pin=line.strip().split()
        samp_hash[pin[0]].append([int(i) for i in pin[1:3]])
    fin.close()
    for k1 in chromos:
        flag1=0
        if not ref_hash[k1]==[]:
            out[k1]=[]
            for k2 in ref_hash[k1]:
                flag2=0
                for k3 in samp_hash[k1]:
                    if k3[1]<k2[0]: continue
                    elif k3[0]>k2[1]: continue
                    else:
                        if float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/float(max(k2[1]-k2[0],k3[1]-k3[0]))>0.5:
                            flag2+=1
                if flag2>0:
                    flag1+=1
            out[k1]=[flag1,len(ref_hash[k1]),len(samp_hash[k1]),float(flag1)/float(len(ref_hash[k1]))]
    return out

def MissedSV_Produce_files(ref_file,samp_file):
    ref_hash={}
    samp_hash={}
    out={}
    for i in chromos:
        ref_hash[i]=[]
        samp_hash[i]=[]
    fin=open(ref_file)
    for line in fin:
        pin=line.strip().split()
        ref_hash[pin[0]].append([int(i) for i in pin[1:3]])
    fin.close()
    fin=open(samp_file)
    for line in fin:
        pin=line.strip().split()
        samp_hash[pin[0]].append([int(i) for i in pin[1:3]])
    fin.close()
    for k1 in chromos:
        flag1=0
        if not ref_hash[k1]==[]:
            out[k1]=[]
            for k2 in ref_hash[k1]:
                flag2=0
                for k3 in samp_hash[k1]:
                    if k3[1]<k2[0]: continue
                    elif k3[0]>k2[1]: continue
                    else:
                        if float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/float(max(k2[1]-k2[0],k3[1]-k3[0]))>0.5:
                            flag2+=1
                if flag2>0:
                    flag1+=1
                else:
                	out[k1].append(k2)
    return out

def bed_writing(filename,hash):
	fo=open(filename,'w')
	for k3 in chromos:
		if k3 in hash.keys():
			for k4 in hash[k3]:
				print >>fo, ' '.join([str(i) for i in k4])
	fo.close()

def ROC_writing(filename,hash):
	fo=open(filename,'w')
	for k1 in hash.keys():
		for k2 in hash[k1].keys():
			for k3 in chromos:
				if k3 in hash[k1][k2].keys():
					print >>fo, ' '.join([str(i) for i in [k1,k2,k3]+hash[k1][k2][k3]])
	fo.close()

def MissSV_writing(filename,hash):
	fo=open(filename,'w')
	for k1 in hash.keys():
		for k2 in hash[k1].keys():
			for k3 in chromos:
				if k3 in hash[k1][k2].keys():
					for k4 in hash[k1][k2][k3]:
						print >>fo, ' '.join([str(i) for i in [k3]+k4+[k1,k2]])
	fo.close()

def MissSV_Compare(File1,File2):
	#report missed SVs in File1 while not in File2
	hash1={}
	hash2={}
	for k1 in chromos:
		hash1[k1]={}
		hash2[k1]={}
	fin=open(File1)
	for line in fin:
		pin=line.strip().split()
		if not pin[3] in hash1[pin[0]].keys():
			hash1[pin[0]][pin[3]]={}
		if not pin[4].upper() in hash1[pin[0]][pin[3]].keys():
			hash1[pin[0]][pin[3]][pin[4].upper()]=[]
		hash1[pin[0]][pin[3]][pin[4].upper()].append([pin[1],pin[2]])
	fin.close()
	fin=open(File2)
	for line in fin:
		pin=line.strip().split()
		if not pin[3] in hash2[pin[0]].keys():
			hash2[pin[0]][pin[3]]={}
		if not pin[4].upper() in hash2[pin[0]][pin[3]].keys():
			hash2[pin[0]][pin[3]][pin[4].upper()]=[]
		hash2[pin[0]][pin[3]][pin[4].upper()].append([pin[1],pin[2]])
	fin.close()
	hash3={}
	for k1 in hash1.keys():
		hash3[k1]={}
		for k2 in hash1[k1].keys():
			hash3[k1][k2]={}
			if k2 in hash2[k1].keys():
				for k3 in hash1[k1][k2].keys():
					hash3[k1][k2][k3]=[]
					if k3 in hash2[k1][k2].keys():
						for k4 in hash1[k1][k2][k3]:
							if not k4 in hash2[k1][k2][k3]:
								hash3[k1][k2][k3].append(k4)
					else:
						hash3[k1][k2][k3]=hash1[k1][k2][k3]
			else:
				hash3[k1][k2]=hash1[k1][k2]
	fo=open(File1+'.vs.'+File2.split('/')[-1],'w')
	for k1 in chromos:
		if k1 in hash3.keys():
			for k2 in hash3[k1].keys():
				for k3 in hash3[k1][k2].keys():
					for k4 in hash3[k1][k2][k3]:
						print >>fo, ' '.join([str(i) for i in [k1]+k4+[k2,k3]])
	fo.close()

def hash_add():
	for k1 in sv_hash2.keys():
		chrom=k1.split('_')[0]
		if chrom in sv_out.keys():
			if k1 in sv_out[chrom].keys():
				sv_out[chrom][k1]+=sv_hash2[k1]
			else:
				sv_out[chrom][k1]=sv_hash2[k1]

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

if not dict_opts['-p'][-1]=='/':
	dict_opts['-p']+='/'

filesIn=[]
for k1 in os.listdir(dict_opts['-p']):
	if dict_opts['--prefix'] in k1 and k1.split('.')[-1]=='vcf' and not 'BPLink' in k1:
		filesIn.append(dict_opts['-p']+k1)

ref=dict_opts['--ref']
chromos=chromos_name_readin(ref)

sv_hash={}
sv_hash2={}
for filein in filesIn:
	fin=open(filein)
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='#':
			if pin[6]=='PASS':
				pos=end_cordi_calcu(pin)
				gt=genotype_calcu(pin)
				if pos[2]-pos[1]>99 and not pos=='Error' and not gt=='Error':
					if 'DEL' in filein.split('/')[-1]:
						k1='a/a'
						if gt=='het':
							k2='/a'
						elif gt=='homo':
							k2='/'
						if not k1 in sv_hash.keys():
							sv_hash[k1]={}
						if not k2 in sv_hash[k1].keys():
							sv_hash[k1][k2]=[]
						sv_hash[k1][k2].append(pos+[0.0]+[pin[2]])
					elif 'INV' in filein.split('/')[-1]:
						k1='a/a'
						if gt=='het':
							k2='a/a^'
						elif gt=='homo':
							k2='a^/a^'
						if not k1 in sv_hash.keys():
							sv_hash[k1]={}
						if not k2 in sv_hash[k1].keys():
							sv_hash[k1][k2]=[]
						sv_hash[k1][k2].append(pos+[0.0]+[pin[2]])
					elif 'DUP' in filein.split('/')[-1]:
						#CopyNum=int(pin[-1].split(':')[1])
						k1='a/a'
						if gt=='het':
							k2='a/aa'
						elif gt=='homo':
							k2='aa/aa'
						if not k1 in sv_hash.keys():
							sv_hash[k1]={}
						if not k2 in sv_hash[k1].keys():
							sv_hash[k1][k2]=[]
						sv_hash[k1][k2].append(pos+[0.0]+[pin[2]])
					else:
						if not pin[2] in sv_hash2.keys():
							sv_hash2[pin[2]]=[]
						sv_hash2[pin[2]].append(pin)

sv_info=sv_hash
del1={}
dup1={}
inv1={}
tra1={}
score_Cff=-50
for k1 in sv_info.keys():
	for k2 in sv_info[k1].keys():
		tra_info_add_new(k1,k2)

sv_out={}
hash_reorder_new()
hash_add()
output_file=dict_opts['-p']+dict_opts['--prefix']+'.BPLink.vcf'
write_VCF_header(output_file)
write_VCF_main(output_file)







