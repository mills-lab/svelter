#!/usr/bin/env python

#!python
#SV.Complex.Output.Process.py SVelter --reference /mnt/EXT/Mills-scratch2/reference/GRCh37/human_g1k_v37.fasta --input-path /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/SVelter/bp_files.comp.het.RD10.sorted.bam/comp.het.RD10.sorted/SPCff3.CluCff3.AlignCff0.0 --ref-sv /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp_het.rec/comp_het.SV.rec
#sys.argv=command.split()
import os
import sys
script_name=sys.argv[0]
if len(sys.argv)<2:
	print 'SV.Complex.Output.Process.py'
	print ''
	print 'this script is used to process output of different algorithms'
	print ''
	print 'Usage:'
	print 'SV.Complex.Output.Process.py [options]  <parameters>'
	print ' '
	print 'Options:'
	print 'SVelter'
	print 'Delly'
	print 'Lumpy'
	print 'Pindel'
	print 'report2stat' 
	print 'comparison'   
	print ' '
	print 'Parameters:'
	print ' '
	print 'Parameters for SVelter/Delly/Lumpy/Pindel:'
	print '--reference: reference genome'
	print '--input-path: path where outputs were kept'
	print '--ref-sv: pre set complex SVs to compare to'
	print ' '
	print 'Parameters for report2stat:'
	print '--reference: reference genome'
	print '--report: .report files'
	print ''
	print 'Parameters for comparison:'
	print '--path: folder contains all .report files'
	print ''
	print 'Parameters for stat-integrate:'
	print '--stat: .stat file'
	print ''
else:
	function_name=sys.argv[1]
	if function_name=='SVelter':
		import getopt
		import numpy
		import re
		import random
		import pickle
		import time
		import datetime
		opts,args=getopt.getopt(sys.argv[2:],'p:',['ref-sv=','reference=','input-path='])
		dict_opts=dict(opts)
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
					#print k1
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
		def modify_comp_new(comp_new):
			for x in comp_new.keys():
				for y in comp_new[x].keys():
					for z in comp_new[x][y]:
						z[2]='/'.join(sorted(z[2].split('/')))
						z[4]='/'.join(sorted(z[4].split('/')))
		def write_comp_new(comp_new,fout):
			#fout=dict_opts['--input-path']+'SVelter.SV.compare.report'
			fo=open(fout,'w')
			for x in comp_new.keys():
				for y in comp_new[x].keys():
					for z in comp_new[x][y]:
						print >>fo,  ' '.join([str(i) for i in [x]+z[0]+z[1:]])
			fo.close()
		def SVelter_info_readin(input_path):
			out={}
			for x in os.listdir(input_path):
				if x.split('.')[-1]=='svelter':
					out[x]={}
					fin=open(input_path+x)
					pin=fin.readline().strip().split()
					for line in fin:
						pin=line.strip().split()
						if not pin[0] in out[x].keys():
							out[x][pin[0]]={}
						if not int(pin[1]) in out[x][pin[0]].keys():
							out[x][pin[0]][int(pin[1])]=[]
						out[x][pin[0]][int(pin[1])].append(pin[3].split(':')+[pin[-1],pin[-3],pin[-2]])
					fin.close()
			return out
		ref_hash={}
		fin=open(dict_opts['--ref-sv'])
		for line in fin:
			pin=line.strip().split()
			if not pin[0] in ref_hash.keys():
				ref_hash[pin[0]]={}
			if not int(pin[1]) in ref_hash[pin[0]].keys():
				ref_hash[pin[0]][int(pin[1])]=[]
			if not [pin[0]]+[int(i) for i in pin[1:4]] in ref_hash[pin[0]][int(pin[1])]:
				ref_hash[pin[0]][int(pin[1])].append(pin)
		fin.close()
		chromos=chromos_name_readin(dict_opts['--reference'])
		if not dict_opts['--input-path'][-1]=='/':
			dict_opts['--input-path']+='/'
		in_hash_all=SVelter_info_readin(dict_opts['--input-path'])
		for file_name in in_hash_all.keys():
			in_hash=in_hash_all[file_name]
			in_list=hash_to_list(in_hash)
			ref_list=hash_to_list(ref_hash)
			in_list=in_list_chrom_filter(in_list)
			list_transfor(in_list)
			list_transfor(ref_list)
			comp_all=compare_list(in_list,ref_list)
			comp_new={}
			x=file_name.split('.')[2]
			for k1 in comp_all.keys():
				comp_new[k1]={}
				for k2 in comp_all[k1].keys():
					comp_new[k1][k2]=[]
					for k3 in comp_all[k1][k2]:
						if bp_check(k3[0][:-3],k3[1][:-3])[0]>0:
							#print k3
							comp_new[k1][k2].append(block_match(k3))
			modify_comp_new(comp_new)
			write_comp_new(comp_new,dict_opts['--input-path']+'SVelter.'+file_name.replace('.vcf','')+'.comp.SV.compare.report')
	if function_name=='Delly':
		import getopt
		import numpy
		import re
		import random
		import pickle
		import time
		import datetime
		opts,args=getopt.getopt(sys.argv[2:],'p:',['ref-sv=','reference=','input-path='])
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
		path_in=dict_opts['--input-path']
		if not path_in[-1]=='/':
			path_in+='/'
		ref=dict_opts['--reference']
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
		fin=open(dict_opts['--ref-sv'])
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
			write_comp_new(comp_new,dict_opts['--input-path']+'Delly.'+k1+'.comp.SV.compare.report')
	if function_name=='Lumpy':
		import getopt
		import numpy
		import re
		import random
		import pickle
		import time
		import datetime
		opts,args=getopt.getopt(sys.argv[2:],'p:',['ref-sv=','reference=','input-path='])
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
		path_in=dict_opts['--input-path']
		if not path_in[-1]=='/':
		  path_in+='/'
		ref=dict_opts['--reference']
		chromos=chromos_name_readin(ref)
		file_hash={}
		for k1 in os.listdir(path_in):
		  if k1.split('.')[-1]=='bedpe' and os.path.isfile(path_in+k1):
		    x1=k1.replace('.bedpe','')
		    if not x1 in file_hash.keys():
		        file_hash[x1]=[]
		    file_hash[x1].append(k1)
		ref_hash={}
		fin=open(dict_opts['--ref-sv'])
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
		    fin=open(dict_opts['--input-path']+k2)
		    for line in fin:
		      pin=line.strip().split()
		      if 'DELETION' in pin[10]:
		        keya='a/a'
		        key_het='a/'
		        key_homo='/'
		      elif 'INVERSION' in pin[10]:      
		        keya='a/a'
		        key_het='a/a^'
		        key_homo='a^/a^'
		      elif 'DUPLICATION' in pin[10]:
		        keya='a/a'      
		        key_het='aa/aa'
		        key_homo='a/aa'
		      pos=[pin[0],int(numpy.median([int(pin[1]),int(pin[2])])),int(numpy.median([int(pin[4]),int(pin[5])]))]
		      key_geno=key_het
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
		  write_comp_new(comp_new,dict_opts['--input-path']+'Lumpy.'+k1+'.comp.SV.compare.report')
	if function_name=='Pindel':
		import getopt
		import numpy
		import re
		import random
		import pickle
		import time
		import datetime
		opts,args=getopt.getopt(sys.argv[2:],'p:',['ref-sv=','reference=','input-path='])
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
		def sv_type_detect(pin):
		  flag=0
		  for x in pin[7].split(';'):
		    if 'SVTYPE' in x:
		      flag=1
		      return x.split('=')[1]
		  if flag==0:
		    return 'Error'
		path_in=dict_opts['--input-path']
		if not path_in[-1]=='/':
		  path_in+='/'
		ref=dict_opts['--reference']
		chromos=chromos_name_readin(ref)
		file_hash={}
		for k1 in os.listdir(path_in):
		  if k1.split('.')[-1]=='vcf' and os.path.isfile(path_in+k1):
		    x1=k1.replace('.vcf','')
		    if not x1 in file_hash.keys():
		        file_hash[x1]=[]
		    file_hash[x1].append(k1)
		ref_hash={}
		fin=open(dict_opts['--ref-sv'])
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
		    fin=open(dict_opts['--input-path']+k2)
		    for line in fin:
		      pin=line.strip().split()
		      if not pin[0][0]=='#':
		        SV_TYPE=sv_type_detect(pin)
		        if not SV_TYPE=='Error':
		          if SV_TYPE=='DEL':
		            keya='a/a'
		            key_het='a/'
		            key_homo='/'
		          elif SV_TYPE=='INV':      
		            keya='a/a'
		            key_het='a/a^'
		            key_homo='a^/a^'
		          elif 'DUP' in SV_TYPE:
		            keya='a/a'      
		            key_het='aa/aa'
		            key_homo='a/aa'
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
		  write_comp_new(comp_new,dict_opts['--input-path']+'Pindel.'+k1+'.comp.SV.compare.report')
	if function_name=='report2stat':
		import getopt
		import glob
		import time
		import datetime
		opts,args=getopt.getopt(sys.argv[2:],'i:o:h:',['reference=','report=','appdix='])
		dict_opts=dict(opts)
		filesIn=dict_opts['--report']
		ref=dict_opts['--reference']
		chromos=[]
		fref=open(ref+'.fai')
		for line in fref:
			pin=line.strip().split()
			chromos.append(pin[0])
		fref.close()
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
		def tra_csv_info_add(k1,k2,tra1):
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
		def tra_info_add_new(k1,k2,tra1):
			for k3 in sv_info[k1][k2]:
				k3=[str(i) for i in k3]
				SV_ID=k3[-1]
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
				if k2a=='':
				    a3='left'
				    a4='right'
				    l_chr=bp_hash[a3][0]
				    r_chr=bp_hash[a4][0]
				    if not 'a' in tra1[SV_ID].keys():
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
				if k2b=='':
				    a3='left'
				    a4='right'
				    l_chr=bp_hash[a3][0]
				    r_chr=bp_hash[a4][0]
				    if not 'b' in tra1[SV_ID].keys():
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
		def hash_reorder_new():
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
		def hash_add():
			for k1 in sv_hash2.keys():
				chrom=k1.split('_')[0]
				if chrom in sv_out.keys():
					if k1 in sv_out[chrom].keys():
						sv_out[chrom][k1]+=sv_hash2[k1]
					else:
						sv_out[chrom][k1]=sv_hash2[k1]
		def let_add(k1,k2,k3):
			#eg: k1='bc/bc',k2='bcbc/bcbc',k3=['20', 24825005, 24825102, 24825872, 24931024]
			kx=k1.split('/')
			ky=k2.split('/')
			k_all=bp_to_let([[str(i) for i in k3]])
			kz=[]
			for x in k_all:
				if not x=='/':
					if not x in k1:
						kz.append(x)
				else:
					break
			#if not ord(k1[0])==97:
			#	for x in range(ord(k1[0]))[97:]:
			#		kz.append(chr(x))
			if not kz==[]:
				kz2=[kz[0]]
				if len(kz)>1:
					for x in range(len(kz)-1):
						if ord(kz[x+1])-ord(kz2[-1][-1])==1:
							kz2[-1]+=kz[x+1]
						else:
							kz2.append(kz[x+1])
				for x in kz2:
					if 'a' in x:
						kx[0]=''.join(x)+kx[0]
						kx[1]=''.join(x)+kx[1]
						ky[0]=''.join(x)+ky[0]
						ky[1]=''.join(x)+ky[1]
					else:
						kx[0]=kx[0]+''.join(x)
						kx[1]=kx[1]+''.join(x)
						ky[0]=ky[0]+''.join(x)
						ky[1]=ky[1]+''.join(x)
				return ['/'.join(kx),'/'.join(ky)]
			else:
				return [k1,k2]
		sv_hash_ref={}
		sv_hash_samp={}
		sv_hash2={}
		fin=open(filesIn)
		for line in fin:
			pin=line.strip().split()
			SV_ID='_'.join([str(i) for i in pin[:-4]])
			Score=0
			k_new=let_add(pin[-2],pin[-1],pin[:-4])
			if not k_new[0] in sv_hash_samp.keys():
				sv_hash_samp[k_new[0]]={}
			if not k_new[1] in sv_hash_samp[k_new[0]].keys():
				sv_hash_samp[k_new[0]][k_new[1]]=[]
			sv_hash_samp[k_new[0]][k_new[1]].append([pin[0]]+[int(i) for i in pin[1:-4]]+[Score,SV_ID])
			k_new=let_add(pin[-4],pin[-3],pin[:-4])	
			if not k_new[0] in sv_hash_ref.keys():
				sv_hash_ref[k_new[0]]={}
			if not k_new[1] in sv_hash_ref[k_new[0]].keys():
				sv_hash_ref[k_new[0]][k_new[1]]=[]
			sv_hash_ref[k_new[0]][k_new[1]].append([pin[0]]+[int(i) for i in pin[1:-4]]+[Score,SV_ID])
		hash_test={}
		fin=open(filesIn)
		for line in fin:
			pin=line.strip().split()
			key1='_'.join([str(i) for i in pin[:-4]])
			if not key1 in hash_test.keys():
				hash_test[key1]=[]
			hash_test[key1].append(pin[-4:])
		fin.close()
		del1={}
		dup1={}
		inv1={}
		tra1={}
		score_Cff=-50
		tra_ref={}
		tra_samp={}
		sv_info=sv_hash_samp
		for k1 in sv_info.keys():
			for k2 in sv_info[k1].keys():
				print [k1,k2]
				tra_info_add_new(k1,k2,tra_samp)
		sv_info=sv_hash_ref
		for k1 in sv_info.keys():
			for k2 in sv_info[k1].keys():
				print [k1,k2]
				tra_info_add_new(k1,k2,tra_ref)
		tra_ref2={}
		for k1 in tra_ref.keys():
			tra_ref2[k1]={}
			for k2 in tra_ref[k1].keys():
				tra_ref2[k1][k2]=[]
				for k3 in tra_ref[k1][k2]:
					if not k3 in tra_ref2[k1][k2]:
						tra_ref2[k1][k2].append(k3)
		def descript_compare(des1,des2):
			#eg: des1=']20:21417806]A', des2=']20:21417806]A'
			if des1==des2:
				return 1
			else:
				flag=0
				if ']' in des1 and not '[' in des1:
					if  ']' in des2 and not '[' in des2:
						if abs(int(des1.split(':')[1].split(']')[0].split('[')[0])-int(des2.split(':')[1].split(']')[0].split('[')[0]))<10:
							flag=1
							return 1 
				elif '[' in des1 and not ']' in des1:
					if  '[' in des2 and not ']' in des2:
						if abs(int(des1.split(':')[1].split(']')[0].split('[')[0])-int(des2.split(':')[1].split(']')[0].split('[')[0]))<10:
							flag=1
							return 1 
				elif '[' in des1 and ']' in des1:
					if  '[' in des2 and ']' in des2:
						if abs(int(des1.split(':')[1].split(']')[0].split('[')[0])-int(des2.split(':')[1].split(']')[0].split('[')[0]))<10:
							flag=1
							return 1 
				if flag==0: 
					return 0
		compare={}
		for k1 in tra_ref2.keys():
			compare[k1]=[]
			if k1 in tra_samp.keys():
				compScores=[]
				ref_list=[]
				samp_list=[]
				for kx in tra_ref2[k1].keys():
					ref_list+=tra_ref2[k1][kx]
				for ky in tra_samp[k1].keys():
					samp_list+=tra_samp[k1][ky]
				comp_score=0
				compScores.append([comp_score,len(ref_list),len(samp_list)])
				for km in ref_list:
					score_flag=0
					for kn in samp_list:
						if km[0]==kn[0] and abs(km[1]-kn[1])<100:
							if descript_compare(km[-1],kn[-1])==1:
								comp_score+=1
								score_flag+=1
								del samp_list[samp_list.index(kn)]
						if not score_flag==0: break
				compScores[-1][0]=comp_score
				print compScores+[k1]
			compare[k1]=compScores
		fout=open(filesIn.replace('.report','.stat'),'w')
		for k1 in compare.keys():
			for k2 in compare[k1]:
				print >>fout, ' '.join([str(i) for i in k2]+[k1]+['_'.join(i) for i in hash_test[k1]])
		fout.close()
	if function_name=='comparison':
		import getopt
		import numpy
		import re
		import random
		import pickle
		import time
		import datetime
		def key_modification(key_level_1):
			if not key_level_1[0]=='a':
				diff=ord(key_level_1[0])-ord('a')
				out=''
				for x in key_level_1:
					if not x in ['/','_','^']:
						out+=chr(ord(x)-diff)
					else:
						out+=x
			else:
				out=key_level_1
			return out
		opts,args=getopt.getopt(sys.argv[2:],'p:x:',['path=','ref-sv='])
		dict_opts=dict(opts)
		if not dict_opts['--path'][-1]=='/':
			dict_opts['--path']+='/'
		ref_sv_hash={}
		fin=open(dict_opts['--ref-sv'])
		for line in fin:
			pin=line.strip().split()
			if not pin[0] in ref_sv_hash.keys():
				ref_sv_hash[pin[0]]={}
			if not int(pin[1]) in ref_sv_hash[pin[0]].keys():
				ref_sv_hash[pin[0]][int(pin[1])]={}
			if not int(pin[-4]) in ref_sv_hash[pin[0]][int(pin[1])].keys():
				ref_sv_hash[pin[0]][int(pin[1])][int(pin[-4])]=[]
			ref_sv_hash[pin[0]][int(pin[1])][int(pin[-4])].append(pin[-2])
			ref_sv_hash[pin[0]][int(pin[1])][int(pin[-4])].append(pin[-1])
		fin.close()
		all_stat_hash={}
		for k1 in os.listdir(dict_opts['--path']):
			if k1.split('.')[-1]=='stat' and not 'integrated' in k1:
				alt_sv_hash={}
				fin=open(dict_opts['--path']+k1)
				for line in fin:
					pin=line.strip().split()
					if not pin[3].split('_')[0] in alt_sv_hash.keys():
						alt_sv_hash[pin[3].split('_')[0]]={}
					if not int(pin[3].split('_')[1]) in alt_sv_hash[pin[3].split('_')[0]].keys():
						alt_sv_hash[pin[3].split('_')[0]][int(pin[3].split('_')[1])]={}
					if not int(pin[3].split('_')[-1]) in alt_sv_hash[pin[3].split('_')[0]][int(pin[3].split('_')[1])].keys():
						alt_sv_hash[pin[3].split('_')[0]][int(pin[3].split('_')[1])][int(pin[3].split('_')[-1])]=[]
					alt_sv_hash[pin[3].split('_')[0]][int(pin[3].split('_')[1])][int(pin[3].split('_')[-1])].append(pin[-1])
				fin.close()
				alt_ref_sv_hash={}
				for x1 in alt_sv_hash.keys():
					if x1 in ref_sv_hash.keys():
						alt_ref_sv_hash[x1]={}
						for x2 in sorted(alt_sv_hash[x1].keys()):
							alt_ref_sv_hash[x1][x2]={}
							for x3 in sorted(alt_sv_hash[x1][x2].keys()):
								alt_ref_sv_hash[x1][x2][x3]=[alt_sv_hash[x1][x2][x3]]
								alt_info_temp=[x2,x3]
								for x4 in sorted(ref_sv_hash[x1].keys()):
									for x5 in sorted(ref_sv_hash[x1][x4].keys()):
										ref_info_temp=[x4,x5]
										if not x5<x2 and not x4>x3:
											alt_ref_sv_hash[x1][x2][x3]+=ref_sv_hash[x1][x4][x5]
				fin=open(dict_opts['--path']+k1)
				for line in fin:
					pin=line.strip().split()
					key_level_1=key_modification('_'.join(pin[-1].split('_')[:2]))
					key_level_2=alt_ref_sv_hash[pin[3].split('_')[0]][int(pin[3].split('_')[1])][int(pin[3].split('_')[-1])]
					if len(key_level_2)==3:
						key_level_1='_'.join(key_level_2[1:])
					if not key_level_1 in all_stat_hash.keys():
						all_stat_hash[key_level_1]={}
					if not k1 in all_stat_hash[key_level_1].keys():
						all_stat_hash[key_level_1][k1]=[[],[]]
					all_stat_hash[key_level_1][k1][0].append(float(pin[0])/float(pin[1]))
					if not pin[2]=='0':
						all_stat_hash[key_level_1][k1][1].append(float(pin[0])/float(pin[2]))
					else:
						all_stat_hash[key_level_1][k1][1].append(1)
				fin.close()
		fo=open(dict_opts['--path']+'BPLink.Integrated.Categorized'+'.stat','w')
		for k1 in sorted(all_stat_hash.keys()):
			for k2 in sorted(all_stat_hash[k1].keys()):
				print >>fo, ' '.join([str(i) for i in [k1,numpy.mean(all_stat_hash[k1][k2][0]),numpy.mean(all_stat_hash[k1][k2][1]),k2]])
		fo.close()
	if function_name=='stat-integrate':
		import getopt
		import numpy
		import re
		import random
		opts,args=getopt.getopt(sys.argv[2:],'i:x:',['stat='])
		dict_opts=dict(opts)
		fin=open(dict_opts['--stat'])
		out_hash={}
		for line in fin:
			pin=line.strip().split()
			if not '_'.join(pin[-1].split('_')[:2]) in out_hash.keys():
				out_hash['_'.join(pin[-1].split('_')[:2])]=[]
			out_hash['_'.join(pin[-1].split('_')[:2])].append(pin)
		fin.close()
		rec=0
		if not os.path.isdir(dict_opts['--stat']+'.sub/'):
			os.system(r'''mkdir %s'''%(dict_opts['--stat']+'.sub/'))
		for x in out_hash.keys():
			rec+=1
			fo=open(dict_opts['--stat']+'.sub/'+dict_opts['--stat'].split('/')[-1].replace('.stat','.'+str(rec)+'.stat'),'w')
			for y in out_hash[x]:
				print >>fo, ' '.join(y)
			fo.close()
		fo=open(dict_opts['--stat'].replace('.stat','.integrated.stat'),'w')
		for x in out_hash.keys():
			Sens=[]
			Spec=[]
			for y in out_hash[x]:
				if not y[2]=='0':
					Sens.append(float(y[0])/float(y[1]))
					Spec.append(float(y[0])/float(y[2]))
				else:
					print y
					Sens.append(0)
					Spec.append(1)
			print >>fo, ' '.join([str(i) for i in [numpy.mean(Sens),numpy.mean(Spec),x]])
		fo.close()

