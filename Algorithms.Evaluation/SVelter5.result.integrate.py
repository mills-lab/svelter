#!/usr/bin/env python

#!Python
#Usage
#SVelter5.result.integrate.py --ppre workdir --ref ref.fa --RSPath path_of_.coverge_files -o out.vcf
#command='python /nfs/remills-data/xuefzhao/SVelter/SVelter5.result.integrate.py --ppre /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/ --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.fasta --RSPath /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/bp_files/Simulate.het.RD.20.RL.101.sorted/SPCff4.CluCff6.AlignCff0.0/ -o /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/SVelter/Simulate.het.RD.20.RL.101.vcf'
#--RSPath: input path under which all .coverge file will be used for analysis
#if not --RSPath specified, will take: workdir(as specified by '--ppre')/bp_files/NA19240/SPCff5.CluCff10.AlignCff0.2
#--OutFile: absolute path of output file, if not specified, output will be;workdir/bp_files//NA19240.SPCff5.CluCff10.AlignCff0.2.vcf
#QCCff: minimum score for a structure to pass quality control
#sys.argv=command.split()
import os
import sys
import getopt
import glob
import time
import datetime
opts,args=getopt.getopt(sys.argv[1:],'i:o:h:',['RSPath=','OutFile=','help=','ppre=','PathBam=','sample=','chr=','SamplePbs=', 'ref=','NullSplitLength=','NullILCI=','NullRDCI=','NullTBCI=','NullILCff=','NullSPCff=','NullDRCff=','NullBed=','NullBedLength=','NullSampleLength=','NullSampleNumber=','ExcludeBed=','IncludeBed=','ToolMappingQ=','FileMappingQ=','NullSamplePercentage=','SplitLength=','BPSPCff=','BPLNCff=','BPAlignQC=','BPAlignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)
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
			if k2.count(i)>2:
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
			k3=[[],[]]
			for x1 in k2.split('/')[0]:
				if k3[0]==[]:
					k3[0].append(x1)
				else:
					if not x1 in k3[0][0]:
						k3[0][0]+=x1
					else:
						if not x1==k3[0][0][-1]:
							index1=k3[0][0].index(x1)-k3[0][0].index(k3[0][0][-1])-(ord(x1)-ord(k3[0][0][-1]))
							index2=1-(ord(x1)-ord(k3[0][0][-1]))
							if abs(index1)>abs(index2): 
								k3[0][0]=k3[0][0].replace(x1,'')
								k3[0][0]+=x1
			for x1 in k2.split('/')[1]:
				if k3[1]==[]:
					k3[1].append(x1)
				else:
					if not x1 in k3[1][0]:
						k3[1][0]+=x1
					else:
						if not x1==k3[1][0][-1]:
							index1=k3[1][0].index(x1)-k3[1][0].index(k3[1][0][-1])-(ord(x1)-ord(k3[1][0][-1]))
							index2=1-(ord(x1)-ord(k3[1][0][-1]))
							if abs(index1)>abs(index2): 
								k3[1][0]=k3[1][0].replace(x1,'')
								k3[1][0]+=x1
				if k3==letter_recombine(k1):
					out='Simple_DUP'
	return out

def index_element(k1,ele):
	out=[]
	for x in range(len(k1)-len(ele)+1):
		if k1[x:(x+len(ele))]==ele:
			out.append(x)
	return out

def simple_DUP_type(k1,k2):
	if simple_DUP_decide(k1,k2)=='Simple_DUP':
		k3=letter_recombine(k2)
		out=[[],[]]
		for x1 in k3[0]:
			index_pos=index_element(k2.split('/')[0],x1)
			if len(index_pos)>1:
				tandem_flag=0
				for y1 in range(len(index_pos)-1):
					if index_pos[y1+1]-index_pos[y1]>len(x1):
						tandem_flag+=1
				if tandem_flag==0:
					if not x1+'_Tandem' in out[0]:
						out[0].append(x1+'_Tandem')
				else:
					if not x1+'_Disperse' in out[0]:
						out[0].append(x1+'_Disperse')
		for x1 in k3[1]:
			index_pos=index_element(k2.split('/')[1],x1)
			if len(index_pos)>1:
				tandem_flag=0
				for y1 in range(len(index_pos)-1):
					if index_pos[y1+1]-index_pos[y1]>len(x1):
						tandem_flag+=1
				if tandem_flag==0:
					if not x1+'_Tandem' in out[1]:
						out[1].append(x1+'_Tandem')
				else:
					if not x1+'_Disperse' in out[1]:
						out[1].append(x1+'_Disperse')
		return out	
	else:
		return 'Error'

def simple_TRA_decide(k1,k2):
	out='Not_Simple_TRA' 
	if '/'.join([''.join(sorted(k2.split('/')[0])),''.join(sorted(k2.split('/')[1]))])==k1:
		out='simple_TRA'
	return out

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

def sv_rec_2(sv_info):
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
	if simple_DEL_decide(k1,k2)=='Simple_DEL':
		for k3 in sv_info[k1][k2]:
			del_info_add(k3,del_let)
	elif simple_DUP_decide(k1,k2)=='Simple_DUP':
		dup_subtype=simple_DUP_type(k1,k2)
		dis_dup_let=[[],[]]
		tan_dup_let=[[],[]]
		for x in range(2):
			if not dup_subtype[x]==[]:
				for y in dup_subtype[x]:
					if y.split('_')[1]=='Disperse':
						dis_dup_let[x].append([y.split('_')[0],k2.split('/')[x].count(y.split('_')[0])])
					else:
						tan_dup_let[x].append([y.split('_')[0],k2.split('/')[x].count(y.split('_')[0])])
		if not tan_dup_let==[[],[]]:
			for k3 in sv_info[k1][k2]:
				dup_info_2_add(k3,tan_dup_let)
		if not dis_dup_let==[[],[]]:
			for k3 in sv_info[k1][k2]:
				disperse_dup_info_2_add(k3,dis_dup_let)
	elif simple_INV_decide(k1,k2)=='Simple_INV':
		for k3 in sv_info[k1][k2]:
			inv_info_add(k3,inv_let)
	elif simple_TRA_decide(k1,k2)=='simple_TRA':
		tra_info_add(k1,k2)
	else:
		for k3 in sv_info[k1][k2]:
			del_csv_info_add(k3,del_let)
			inv_csv_info_add(k3,inv_let)
			dup_csv_info_2_add(k3,dup_let)
		if csv1[3]==1:
			tra_csv_info_add(k1,k2)

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

def disperse_dup_info_2_add(k3,dup_let):
	temprec=-1
	for k2x in dup_let:
		temprec+=1
		hetx=['heta','hetb'][temprec]
		for k4 in k2x:
			temp=bp_to_hash(k3[:-1],[i for i in k4[0]])
			for k5 in temp:
			    if not k5[0] in disperse_dup.keys():
			        disperse_dup[k5[0]]=[]
			    if k4[1]>1:
				    disperse_dup[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S']),k4[1]])

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
		if not SV_ID in tra1.keys():
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
				ref_allele[a3].append(ref_base_readin(ref,bp_hash[a3][0],a4))
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
		        if not 'a' in tra1[SV_ID].keys():
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
		    	if not 'a' in tra1[SV_ID].keys():
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
		        if not 'b' in tra1[SV_ID].keys():  	    
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
		    	if not 'b' in tra1[SV_ID].keys():
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

def tra_info_add(k1,k2):
	for k3 in sv_info[k1][k2]:
		SV_ID='_'.join([str(i) for i in k3]+['S'])
		if not SV_ID in tra1.keys():
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
				ref_allele[a3].append(ref_base_readin(ref,bp_hash[a3][0],a4))
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
		        if not 'a' in tra1[SV_ID].keys(): 	   
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
		    	if not 'a' in tra1[SV_ID].keys():
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
		        if not 'b' in tra1[SV_ID].keys(): 	  
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
		    	if not 'b' in tra1[SV_ID].keys():
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

def hash_collaps():
	for k1 in sv_out.keys():
		for k2 in sv_out[k1].keys():
			if len(sv_out[k1][k2])>1:
				temp=[]
				temp2=[]
				for k3 in sv_out[k1][k2]:
					if not k3[:-1] in temp:
						temp.append(k3[:-1])
						temp2.append([k3[-1]])
					else:
						temp2[temp.index(k3[:-1])].append(k3[-1])
				for k3 in range(len(temp2)):
					if len(temp2[k3])>1:
						if sorted([temp2[k3][0].split(':')[0],temp2[k3][1].split(':')[0]])==['0|1', '1|0']:
							if not ':' in temp2[k3][0]:
								temp2[k3]=['1|1']
							else:	
								temp2[k3]=['1|1:'+str(int(temp2[k3][0].split(':')[1])+int(temp2[k3][1].split(':')[1]))]
				temp3=[]
				for k3 in range(len(temp2)):
					temp3.append(temp[k3]+temp2[k3])
				sv_out[k1][k2]=temp3

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

def hash_collaps2():
	temp={}
	for k1 in sv_out.keys():
		temp[k1]={}
		for k2 in sv_out[k1].keys():
			for k3 in sv_out[k1][k2]:
				pos=end_cordi_calcu(k3)
				if not pos[1] in temp[k1].keys():
					temp[k1][pos[1]]={}
				if not pos[2] in temp[k1][pos[1]].keys():
					temp[k1][pos[1]][pos[2]]=[]
				temp[k1][pos[1]][pos[2]].append([k1,k2,k3])
	out={}
	for k1 in temp.keys():
		out[k1]={}
		for k2 in temp[k1].keys():
			if len(temp[k1][k2])>1:
				flag=1
				for k3 in temp[k1][k2].keys():
					for k4 in temp[k1][k2][k3]:
							if not k4[2][4]=='<DUP>':
								flag=0
				if flag==1:
					for k4 in temp[k1][k2].keys():
						if not k4==max(temp[k1][k2].keys()):
							for k5 in temp[k1][k2][k4]:
								del sv_out[k1][k5[2][2]][sv_out[k1][k5[2][2]].index(k5[2])]
								if sv_out[k1][k5[2][2]]==[]:
									del sv_out[k1][k5[2][2]]

def hash_collaps3():
	for k1 in sv_out.keys():
		for k2 in sv_out[k1].keys():
			if len(sv_out[k1][k2])>1:
				temp1=[]
				temp2=[]
				for k3 in range(len(sv_out[k1][k2])):
					if not sv_out[k1][k2][k3][:5]+sv_out[k1][k2][k3][6:-1] in temp1:
						temp1.append(sv_out[k1][k2][k3][:5]+sv_out[k1][k2][k3][6:-1])
						temp2.append(sv_out[k1][k2][k3])
					else:
						continue
				sv_out[k1][k2]=temp2

def hash_reorder():
	for ka1 in del1.keys():
		if not ka1 in sv_out.keys():
			sv_out[ka1]={}
		for ka2 in del1[ka1]:
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
			ka_new=[ka1,ka2[0],ka2[-2],REF_AL,'<DUP:TANDEM>',ka2[3],Pass_Sign,'SVTYPE=DUP;END='+str(ka2[1]),'GT:CN',GenoType+':'+CopyNumber]
			if not ka2[-2] in sv_out[ka1].keys():
				sv_out[ka1][ka2[-2]]=[]
			if not ka_new in sv_out[ka1][ka2[-2]]:
				sv_out[ka1][ka2[-2]].append(ka_new)
	for ka1 in disperse_dup.keys():
		if not ka1 in sv_out.keys():
			sv_out[ka1]={}
		for ka2 in disperse_dup[ka1]:
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
		ks2='_'.join(ka1.split('_')[:-2]+[ka1.split('_')[-1]])
		SV_Score=float(ka1.split('_')[-2])
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
						if k4[3]=='N':
							k4[3]=ref_base_readin(ref,k4[0],k4[1])
						print >>fo, '\t'.join([str(i) for i in k4])
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

def dup_collaps(dup1):
	out={}
	for k1 in dup1.keys():
		out[k1]={}
		for k2 in dup1[k1]:
			if not k2[0] in out[k1].keys():
				out[k1][k2[0]]={}
			if not k2[1] in out[k1][k2[0]].keys():
				out[k1][k2[0]][k2[1]]=[]
			out[k1][k2[0]][k2[1]].append(k2)
	out=dup_redun_refine(out)
	out=dup_homo_refine(out)
	out=dup_name_refine(out)
	out2=dup_collaps_refine2(out)
	out_final={}
	for k1 in out2.keys():
		out_final[k1]=[]
		for k2 in sorted(out2[k1].keys()):
			for k3 in sorted(out2[k1][k2].keys()):
				for k4 in out2[k1][k2][k3]:
					out_final[k1].append(k4)
	return out_final

def dup_redun_refine(out):
	for k1 in out.keys():
		for k2 in sorted(out[k1].keys()):
			for k3 in sorted(out[k1][k2].keys()):
				if len(out[k1][k2][k3])>1:
					temp_dup=[]
					temp_dup2=[]
					for k4 in out[k1][k2][k3]:
						if not k4[:3]+k4[4:-1] in temp_dup:
							temp_dup.append(k4[:3]+k4[4:-1])
							temp_dup2.append(k4)
					out[k1][k2][k3]=temp_dup2
	return out

def dup_homo_refine(out):
	for k1 in out.keys():
		for k2 in sorted(out[k1].keys()):
			for k3 in sorted(out[k1][k2].keys()):
				if len(out[k1][k2][k3])>1:
					temp_dup=[]
					temp_dup2=[]
					for k4 in out[k1][k2][k3]:
						if not k4[:2]+k4[4:] in temp_dup:
							temp_dup.append(k4[:2]+k4[4:])
							temp_dup2.append(k4)
						else:
							temp_dup2[temp_dup.index(k4[:2]+k4[4:])][2]='homo'
					out[k1][k2][k3]=temp_dup2
	return out

def dup_name_refine(out):
	for k1 in out.keys():
		for k2 in sorted(out[k1].keys()):
			for k3 in sorted(out[k1][k2].keys()):
				if len(out[k1][k2][k3])>1:
					out[k1][k2][k3]=[out[k1][k2][k3][0]]
	return out

def dup_collaps_refine(out):
	outout={}
	for k1 in out.keys():
		outout[k1]={}
		out2=[]
		for k2 in sorted(out[k1].keys()):
			for k3 in sorted(out[k1][k2].keys()):
				out2.append([k2,k3])
		out3=[out2[0]]
		for k2 in out2[1:]:
			add_flag=0
			for k3 in range(len(out3)):
				if k2[0]==out3[k3][-1]:
					out3[k3]+=k2
					add_flag=1
			if add_flag==0:
				out3.append(k2)
		for k2 in out3:
			if len(k2)>2:
				rec_temp=[]
				for k3 in range(len(k2)/2):
					rec_temp+=out[k1][k2[k3*2]][k2[k3*2+1]]
				rec_t2=[rec_temp[0]]
				for k3 in range(len(rec_temp)-1):
					if abs(rec_temp[k3+1][-1]-rec_temp[k3][-1])<2:
						rec_t2[-1][1]=rec_temp[k3+1][1]
					else:
						rec_t2.append(rec_temp[k3+1])
			else:
				rec_t2=out[k1][k2[0]][k2[1]]
			for k3 in rec_t2:
				if not k3[0] in outout[k1].keys():
					outout[k1][k3[0]]={}
				if not k3[1] in outout[k1][k3[0]].keys():
					outout[k1][k3[0]][k3[1]]=[]
				outout[k1][k3[0]][k3[1]].append(k3)
	return outout

def dup_collaps_refine2(out):
	out2=dup_collaps_refine(out)
	for k1 in out2.keys():
		for k2 in out2[k1].keys():
			for k3 in out2[k1][k2].keys():
				if len(out2[k1][k2][k3])>1:
					test=[]
					for k4 in out2[k1][k2][k3]:
						test.append(k4[2])
					if 'homo' in test:
						out2[k1][k2][k3]=[out2[k1][k2][k3][test.index('homo')]]
					else:
						out2[k1][k2][k3]=[out2[k1][k2][k3][0]]
	return out2

def SV_Info_Write_svelter(sv_info):
	temp1={}
	for k1 in sv_info.keys():
		for k2 in sv_info[k1].keys():
			for k3 in sv_info[k1][k2]:
				if not k3[0] in temp1.keys():
					temp1[k3[0]]={}
				if not int(k3[1]) in temp1[k3[0]].keys():
					temp1[k3[0]][int(k3[1])]={}
				if not int(k3[-2]) in temp1[k3[0]][int(k3[1])].keys():
					temp1[k3[0]][int(k3[1])][int(k3[-2])]=[]
				temp1[k3[0]][int(k3[1])][int(k3[-2])].append(k3+[k1,k2])
	fo=open(output_file.replace('.vcf','.svelter'),'w')
	print >>fo, '\t'.join(['chr','start','end','bp_info','ref','alt','score'])
	for k1 in chromos:
		if k1 in temp1.keys():
			for k2 in sorted(temp1[k1].keys()):
				for k3 in sorted(temp1[k1][k2].keys()):
					for k4 in temp1[k1][k2][k3]:
						chrom_svelter=k1
						bp_start_svelter=k2
						bp_end_svelter=k3
						bps_info_svelter=':'.join(k4[:-3])
						struc_ref_svelter=k4[-2]
						struc_alt_svelter=k4[-1]
						score_svelter=k4[-3]
						print >>fo, '\t'.join([str(i) for i in [chrom_svelter,bp_start_svelter,bp_end_svelter,bps_info_svelter,struc_ref_svelter,struc_alt_svelter,score_svelter]])
	fo.close()



if not '--QCCff' in dict_opts:
	score_Cff=-20
else:
	score_Cff=int(dict_opts['--QCCff'])

if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
	print 'SVelter-0.1          Last Update:2014-08-20'
	print 'Required Parameters:'
	print '--ppre, workind directory of SVelter, eg: .../SVelter/' 
	print '--RSPath, path of .coverage files'
	print '-o, absolute path of output file'
	print '--ref, absolute path of reference genome'
	print 'Optional Parameters:'
else:
	if not '--ppre' in dict_opts.keys():
		print 'Error: please specify working directory using: --ppre'
	else:
		workdir=dict_opts['--ppre']
		if not workdir[-1]=='/':
		    workdir+='/'
		if not '--RSPath' in dict_opts.keys():
			print 'Error: please specify path of input .coverge files using --RSPath'
		else:
			if '--RSPath' in dict_opts.keys():
				#--RS: resolved structure, path of namely .coverge files
				if not dict_opts['--RSPath'][-1]=='/':
					dict_opts['--RSPath']+='/'
				InputPath=[dict_opts['--RSPath']]
			else:
				InputPath=[]
				for fi1 in os.listdir(workdir+'bp_files'):
					path1=workdir+'bp_files/'+fi1+'/'
					if os.path.isdir(path1):
						for fi2 in os.listdir(path1):
							path2=path1+fi2+'/'
							InputPath.append(path2)
			ref_file=0
			if '--ref' in dict_opts.keys():
				ref_file=dict_opts['--ref']
			else:
				ref_path=workdir+'reference/'
				for k1 in os.listdir(ref_path):
					if k1.split('.')[-1] in ['fa','fasta']:
						ref_file=ref_path+k1
			if ref_file==0:
				print 'Error: no valid reference genome available; pls specify ref genome using --ref'
			else:
				if not '-o' in dict_opts.keys():
					print 'Error: please specify output file using -o'
				else:
					output_file=dict_opts['-o']
					ref_index=ref_file+'.fai'
					ref=ref_file
					chromos=[]
					fin=open(ref_index)
					for line in fin:
						pin=line.strip().split()
						chromos.append(pin[0])
					fin.close()
					for path2 in InputPath:
						#output_file=path2+'SV.Output.vcf'
						sv_info={}
						for k3 in os.listdir(path2):
							if k3.split('.')[-1]=='coverge':
								read_in_structures(path2+k3)
						SV_Info_Write_svelter(sv_info)
						dup1={}
						disperse_dup={}
						inv1={}
						del1={}
						tra1={}
						sv_rec_2(sv_info)
						dup1=dup_collaps(dup1)
						sv_out={}
						hash_reorder()
						hash_collaps()
						hash_collaps2()
						hash_collaps3()
						write_VCF_header(output_file)
						write_VCF_main(output_file)
                        #temp_rec=workdir+dict_opts['--ref']+'.temp'
                        #fo=open(temp_rec,'w')
                        #print >>fo,' '.join(['SVelter','5','done'])
                        #fo.close() 


