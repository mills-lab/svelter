#command='python vcf.BPlink.compare.py --in /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/Delly/Delly_comp_RD10.BPLink.vcf --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/sv.rec/comp.vcf'
#command='python vcf.BPlink.compare.py --in /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/SVelter/SVelter_comp_RD40.BPlink.vcf --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/sv.rec/comp.vcf'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import datetime
opts,args=getopt.getopt(sys.argv[1:],'i:o:h:',['ref=','in='])
dict_opts=dict(opts)
fileref=dict_opts['--ref']
filein=dict_opts['--in']
def vcf_readin(vcf_file_name):
	ref_hash={}
	fin=open(vcf_file_name)
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='#':
			if pin[6]=='PASS':
				if not pin[2] in ref_hash.keys():
					ref_hash[pin[2]]=[]
				ref_hash[pin[2]].append([pin[0],pin[1],pin[4]])
	fin.close()
	return ref_hash

def hash_reorder_ref(in_hash):
	out_hash={}
	for k1 in in_hash.keys():
		if 'S' in k1 or 'C' in k1:
			chrom=k1.split('_')[0]
			start=int(k1.split('_')[1])
			end=int(k1.split('_')[-2])
		else:
			chrom=k1.split('_')[0]
			start=int(k1.split('_')[1])
			end=int(k1.split('_')[-1])			
		if not chrom in out_hash.keys():
			out_hash[chrom]={}
		if not start in out_hash[chrom].keys():
			out_hash[chrom][start]={}
		if not end in out_hash[chrom][start].keys():
			out_hash[chrom][start][end]=[]
		out_hash[chrom][start][end]+=[k1]+in_hash[k1]
	return out_hash

def hash_reorder_in(in_hash):
	out_hash={}
	for k1 in in_hash.keys():
		chrom=k1.split('_')[0]
		start=int(k1.split('_')[1])
		end=int(k1.split('_')[-2])
		if not chrom in out_hash.keys():
			out_hash[chrom]={}
		if not start in out_hash[chrom].keys():
			out_hash[chrom][start]={}
		if not end in out_hash[chrom][start].keys():
			out_hash[chrom][start][end]=[]
		out_hash[chrom][start][end]+=[k1]+in_hash[k1]
	return out_hash

def hash_reorder2(in_hash):
	ref1={}
	for k1 in in_hash.keys():
		ref1[k1]=[]
		for k2 in sorted(in_hash[k1].keys()):
			for k3 in sorted(in_hash[k1][k2].keys()):
				ref1[k1].append([k2,k3])
	return ref1

import numpy
def hash_link(ref_order,in_order):
	ref1=hash_reorder2(ref_order)
	in1=hash_reorder2(in_order)
	out1={}
	out_not_assigned=[]
	for k1 in ref1.keys():
		if k1 in in1.keys():
			out1[k1]={}
			for k2 in in1[k1]:
				for k3 in ref1[k1]:
					if k2[0]>k3[1] or k2[1]<k3[0]: continue
					else:
						if float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/min([float(k2[1]-k2[0]),float(k3[1]-k3[0])])>0.5:								
							if not '_'.join([str(i) for i in k3]) in out1[k1].keys():
								out1[k1]['_'.join([str(i) for i in k3])]=[]
							out1[k1]['_'.join([str(i) for i in k3])].append(k2)
			for k3 in ref1[k1]:
				if not '_'.join([str(i) for i in k3]) in out1[k1].keys():
					out1[k1]['_'.join([str(i) for i in k3])]=[]
	scores={}
	score2={}
	for k1 in out1.keys():
		scores[k1]=[]
		score2[k1]=[]
		for k2 in out1[k1].keys():
			ref_fac=ref_order[k1][int(k2.split('_')[0])][int(k2.split('_')[1])]
			temp_fac=[[],[]]
			for k3 in out1[k1][k2]:
				in_fac=[]
				in_fac+=in_order[k1][k3[0]][k3[1]]
				tscore=overlap_score(ref_fac,in_fac)
				if not tscore[0]==0:
					temp_fac[0].append(tscore[0])
					temp_fac[1].append(tscore[1])
			if not temp_fac==[[],[]]:
				scores[k1].append(numpy.mean(temp_fac[0]))
				score2[k1].append(numpy.mean(temp_fac[1]))
			else:
				scores[k1].append(0)
				score2[k1].append(0)
	for k1 in scores.keys():
		scores[k1]=numpy.mean(scores[k1])
		score2[k1]=numpy.sum(score2[k1])		
	return [scores,score2]

def hash_link2(ref_order,in_order):
	ref1=hash_reorder2(ref_order)
	in1=hash_reorder2(in_order)
	out1={}
	out_not_assigned=[]
	for k1 in ref1.keys():
		if k1 in in1.keys():
			out1[k1]={}
			for k2 in in1[k1]:
				for k3 in ref1[k1]:
					if k2[0]>k3[1] or k2[1]<k3[0]: continue
					else:
						if float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/min([float(k2[1]-k2[0]),float(k3[1]-k3[0])])>0.5:								
							if not '_'.join([str(i) for i in k3]) in out1[k1].keys():
								out1[k1]['_'.join([str(i) for i in k3])]=[]
							out1[k1]['_'.join([str(i) for i in k3])].append(k2)
			for k3 in ref1[k1]:
				if not '_'.join([str(i) for i in k3]) in out1[k1].keys():
					out1[k1]['_'.join([str(i) for i in k3])]=[]
	scores={}
	score2={}
	for k1 in out1.keys():
		scores[k1]=[]
		score2[k1]=[]
		for k2 in out1[k1].keys():
			ref_fac=ref_order[k1][int(k2.split('_')[0])][int(k2.split('_')[1])]
			temp_fac=[[],[]]
			in_fac=[]
			for k3 in out1[k1][k2]:
				in_fac+=in_order[k1][k3[0]][k3[1]]
			tscore=overlap_score(ref_fac,in_fac)
			if not out1[k1][k2]==[]:
				scores[k1].append(tscore[0])
				score2[k1].append(tscore[1])
			else:
				scores[k1].append(0)
				score2[k1].append(0)
	for k1 in scores.keys():
		scores[k1]=numpy.mean(scores[k1])
		score2[k1]=numpy.sum(score2[k1])		
	return [scores,score2]

def end_calcu(info):
	#eg of info:']18:67545136]N'
	out=info.split(':')[1]
	if ']' in out:
		return int(out.split(']')[0])
	if '[' in out:
		return int(out.split('[')[0])
	else:
		return 'error'

def sv_calcu(info):
	out=info.split(':')
	out2=[]
	if '[' in out[0]:
		out2.append('[')
	elif ']' in out[0]:
		out2.append(']')
	else:
		out2.append('error')
	if '[' in out[1]:
		out2.append('[')
	elif ']' in out[1]:
		out2.append(']')
	else:
		out2.append('error')
	return out2

comp_CI=[-51,51]
def overlap_score(ref_fac,in_fac):
	ref_new={}
	for x in ref_fac:
		if type(x)==type([]):
			if not int(x[1]) in ref_new.keys():
				ref_new[int(x[1])]={}
			if not end_calcu(x[2])=='error':
				if not end_calcu(x[2]) in ref_new[int(x[1])].keys():
					ref_new[int(x[1])][end_calcu(x[2])]=[]
			else:
				print x
			ref_new[int(x[1])][end_calcu(x[2])].append(x)
	in_new={}
	for x in in_fac:
		if type(x)==type([]):
			if not int(x[1]) in in_new.keys():
				in_new[int(x[1])]={}
			if not end_calcu(x[2])=='error':
				if not end_calcu(x[2]) in in_new[int(x[1])].keys():
					in_new[int(x[1])][end_calcu(x[2])]=[]
			else:
				print x
			in_new[int(x[1])][end_calcu(x[2])].append(x)
	ref_n2=[]
	for x in sorted(ref_new.keys()):
		for y in sorted(ref_new[x].keys()):
			ref_n2+=[[x,y]+i for i in ref_new[x][y]]
	in_n2=[]
	for x in sorted(in_new.keys()):
		for y in sorted(in_new[x].keys()):
			in_n2+=[[x,y]+i for i in in_new[x][y]]
	t3=0
	for x1 in ref_n2:
		temp=[]
		temp2=sv_calcu(x1[4])
		t2=0
		for y1 in in_n2:
			if abs(y1[0]-x1[0])<comp_CI[1]:
				if abs(y1[1]-x1[1])<comp_CI[1]:
					temp.append(sv_calcu(y1[4]))
		for y1 in temp:
			if y1==temp2:
				t2=1
		t3+=t2
	score=float(t3)/float(len(ref_n2))
	score2=t3
	return [score,score2]

def calcu_per_chrom(ref_hash):
	out={}
	for x in ref_hash.keys():
		if not x.split('_')[0] in out.keys():
			out[x.split('_')[0]]=len(ref_hash[x])
		else:
			out[x.split('_')[0]]+=len(ref_hash[x])
	return out

def hash_name_link(ref_order,in_order):
	ref_list={}
	in_list={}
	for k1 in ref_order.keys():
		ref_list[k1]=[]
		for k2 in sorted(ref_order[k1].keys()):
			for k3 in ref_order[k1][k2].keys():
				ref_list[k1].append([k2,k3])
	for k1 in in_order.keys():
		in_list[k1]=[]
		for k2 in sorted(in_order[k1].keys()):
			for k3 in in_order[k1][k2].keys():
				in_list[k1].append([k2,k3])
	link_hash={}
	for k1 in ref_list.keys():
		if k1 in in_list.keys():
			link_hash[k1]={}
			for k2 in ref_list[k1]:
				link_hash[k1]['_'.join([str(i) for i in k2])]=[]
				for k3 in in_list[k1]:
					if k3[0]>k2[1]: break
					elif k3[1]<k2[0]: continue
					else:
							link_hash[k1]['_'.join([str(i) for i in k2])].append(k3)
	score_sens={}
	score_allnum={}
	score_Not_missed={}
	for k1 in link_hash.keys():
		score_sens[k1]=[]
		score_allnum[k1]=[]
		score_Not_missed[k1]=0
		for k2 in link_hash[k1].keys():
			ref_detail=ref_order[k1][int(k2.split('_')[0])][int(k2.split('_')[1])]
			ref_d2=[]
			temp_score=[]
			for x in ref_detail:
				if type(x)==type([]):
					ref_d2.append(x)
			for k3 in link_hash[k1][k2]:
				in_d2=[]
				in_detail=in_order[k1][k3[0]][k3[1]]
				for x in in_detail:
					if type(x)==type([]):
						in_d2.append(x)
				temp_score.append(overlap_score(ref_d2,in_d2))
			if not link_hash[k1][k2]==[]:
				score_sens[k1].append(max([i[0] for i in temp_score]))
				score_allnum[k1].append(max([i[1] for i in temp_score]))
			else:
				score_Not_missed[k1]+=1
		score_Not_missed[k1]=1-float(score_Not_missed[k1])/float(len(link_hash[k1]))
	for k1 in score_sens.keys():
		score_sens[k1]=numpy.mean(score_sens[k1])
		score_allnum[k1]=numpy.sum(score_allnum[k1])
	return [score_sens,score_allnum,score_Not_missed]

def hash_name_link2(ref_order,in_order):
	ref_list={}
	in_list={}
	for k1 in ref_order.keys():
		ref_list[k1]=[]
		for k2 in sorted(ref_order[k1].keys()):
			for k3 in ref_order[k1][k2].keys():
				ref_list[k1].append([k2,k3])
	for k1 in in_order.keys():
		in_list[k1]=[]
		for k2 in sorted(in_order[k1].keys()):
			for k3 in in_order[k1][k2].keys():
				in_list[k1].append([k2,k3])
	link_hash={}
	for k1 in ref_list.keys():
		if k1 in in_list.keys():
			link_hash[k1]={}
			for k2 in ref_list[k1]:
				link_hash[k1]['_'.join([str(i) for i in k2])]=[]
				for k3 in in_list[k1]:
					if k3[0]>k2[1]: break
					elif k3[1]<k2[0]: continue
					else:
						if float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/float(max([k2[1]-k2[0],k3[1]-k3[0]]))>0.5:
							link_hash[k1]['_'.join([str(i) for i in k2])].append(k3)
	score_sens={}
	score_allnum={}
	for k1 in link_hash.keys():
		score_sens[k1]=[]
		score_allnum[k1]=[]
		for k2 in link_hash[k1].keys():
			ref_detail=ref_order[k1][int(k2.split('_')[0])][int(k2.split('_')[1])]
			ref_d2=[]
			temp_score=[]
			for x in ref_detail:
				if type(x)==type([]):
					ref_d2.append(x)
			for k3 in link_hash[k1][k2]:
				in_d2=[]
				in_detail=in_order[k1][k3[0]][k3[1]]
				for x in in_detail:
					if type(x)==type([]):
						in_d2.append(x)
				temp_score.append(overlap_score(ref_d2,in_d2))
			if not link_hash[k1][k2]==[]:
				score_sens[k1].append(sum([i[0] for i in temp_score]))
				score_allnum[k1].append(numpy.sum([i[1] for i in temp_score]))
			else:
				score_sens[k1].append(0)
				score_allnum[k1].append(0)
	for k1 in score_sens.keys():
		score_sens[k1]=numpy.mean(score_sens[k1])
		score_allnum[k1]=numpy.sum(score_allnum[k1])
	return [score_sens,score_allnum]


ref_hash=vcf_readin(fileref)
in_hash=vcf_readin(filein)

ref_order=hash_reorder_ref(ref_hash)
in_order=hash_reorder_in(in_hash)


out=hash_name_link(ref_order,in_order)

#hash_name_link2(ref_order,in_order)
#out2=hash_link2(ref_order,in_order)
out_ref=calcu_per_chrom(ref_hash)
out_in=calcu_per_chrom(in_hash)
fout=filein.replace('.vcf','.stats')
fo=open(fout,'w')
#print >>fo, ' '.join(['sample','chrom','sens','spec','#inref','#allFound'])
for k1 in out[0].keys():
	print >>fo, ' '.join([filein.split('/')[-1].split('.')[0],k1,str(out[0][k1]),str(float(out[1][k1])/float(out_in[k1])),str(out[2][k1]),str(out_ref[k1]),str(out_in[k1])])
	print 		' '.join([filein.split('/')[-1].split('.')[0],k1,str(out[0][k1]),str(float(out[1][k1])/float(out_in[k1])),str(out[2][k1]),str(out_ref[k1]),str(out_in[k1])])

fo.close()




