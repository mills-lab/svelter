#!/usr/bin/python

#!Python
#Script to analysis and realign pacbio reads aligned by blasr
#Usage
#command='Realign.Pacbio -i SRR1502785.sam'
#sys.argv=command.split()
import os
import sys
import getopt
import re
opts,args=getopt.getopt(sys.argv[1:],'i:')
dict_opts=dict(opts)

min_length=50
def cigar_modify(cigar):
	temp=[[]]
	for x in cigar:
		if ord(x)>64 and ord(x)<91:
			temp[-1].append(x)
			temp.append([])
		else:
			if temp[-1]==[]:
				temp[-1].append(x)
			else:
				temp[-1][-1]+=x
	temp.remove([])
	for x in temp:
		if x[1]=='S' and not int(x[0])>min_length:
			x[1]='I'
	return ''.join(''.join(y) for y in temp)

def cigar_integrate(pin):
	pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
	cigar=cigar_modify(pin[5])
	pos_start=int(pin[3])
	cigars=[]
	for m in pcigar.finditer(cigar):
		cigars.append((m.groups()[0],m.groups()[1]))
	read_start=0
	read_end=0
	Seq_pos=[]
	Read_Pos=[]
	Map_Pos=[]
	rec_seq=0
	rec_map=int(pin[3])
	for x in cigars:
		if x[1]=='I':
			rec_seq1=rec_seq
			rec_seq=rec_seq1+int(x[0])
			Seq_pos.append(pin[9][rec_seq1:rec_seq])
			Read_Pos.append([rec_seq1,rec_seq])
			Map_Pos.append([rec_map,rec_map])
		if x[1]=='D' or x[1]=='N':
			Seq_pos.append('')
			Read_Pos.append([rec_seq,rec_seq])
			rec_map1=rec_map
			rec_map=rec_map1+int(x[0])
			Map_Pos.append([rec_map1,rec_map])
		if x[1]=='M':
			rec_seq1=rec_seq
			rec_seq=rec_seq1+int(x[0])
			rec_map1=rec_map
			rec_map=rec_map1+int(x[0])
			Seq_pos.append(pin[9][rec_seq1:rec_seq])
			Read_Pos.append([rec_seq1,rec_seq])
			Map_Pos.append([rec_map1,rec_map])
		if x[1]=='S':
			rec_seq1=rec_seq
			rec_seq=rec_seq1+int(x[0])
			rec_map1=rec_map
			rec_map=rec_map1
			Seq_pos.append(pin[9][rec_seq1:rec_seq])
			Read_Pos.append([rec_seq1,rec_seq])
			Map_Pos.append([rec_map1,rec_map])
	Cigar_group=[]
	for x in cigars:
		if int(x[0])>min_length and not x[1]=='M':
			Cigar_group.append([x])
			Cigar_group.append([])
		else:
			if Cigar_group==[]:
				Cigar_group.append([x])
			else:
				Cigar_group[-1]+=[x]
	Seq2_pos=[]
	Read2_Pos=[]
	Map2_Pos=[]
	rec=-1
	for x in Cigar_group:
		Seq2_pos.append([])
		Read2_Pos.append([])
		Map2_Pos.append([])
		for y in x:
			rec+=1
			Seq2_pos[-1].append(Seq_pos[rec])
			Read2_Pos[-1].append(Read_Pos[rec])
			Map2_Pos[-1].append(Map_Pos[rec])
	chopped_reads=[''.join(x) for x in Seq2_pos if not x==[]]
	chopped_Read_Pos=[[x[0][0],x[-1][1]] for x in Read2_Pos if not x==[]]
	chopped_Map_Pos=[[x[0][0],x[-1][1]] for x in Map2_Pos if not x==[]]
	return [chopped_reads,chopped_Read_Pos,chopped_Map_Pos,[x for x in Cigar_group if not x==[]]]

def assign_readname(pin,out):
	if not pin[0]+'.0' in out:
		out.append(pin[0]+'.0')
		return pin[0]+'.0'
	else:
		rec=0
		while True:
			rec+=1
			if not pin[0]+'.'+str(rec) in out:
				out.append(pin[0]+'.'+str(rec))
				return pin[0]+'.'+str(rec)
			else:
				continue

def chop_x_for_fasta(x):
	if len(x)>60:
		lines=len(x)/60
		out=[]
		for k in range(lines):
			out.append(x[(k*60):((k+1)*60)])
		out.append(x[((k+1)*60):])
		return out
	else:
		return [x]

def reads_read_in(filein,fasta_out,rec_out):
	#fasta_out=filein.replace('.sam','.sr.extracted.fasta')
	#sam_out=filein.replace('.sam','.sr.extracted.rec')
	#if not os.path.isfile(fasta_out):
	#	fo=open(fasta_out,'w')
	#	fo.close()
	#if not os.path.isfile(sam_out):
	#	fo=open(sam_out,'w')
	#	fo.close()
	fin=open(filein)
	read_hash={}
	out=[]
	fo1=open(fasta_out,'a')
	#fo2=open(rec_out,'a')
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='@': 
			read_name=assign_readname(pin,out)
			chopped_parts=cigar_integrate(pin)
			#sam_info={}
			#fasta_reads={}
			if len(chopped_parts[0])==1:
				print >>fo_final, '\t'.join(pin)
				print pin
			else:
				rec2=0
				for x in range(len(chopped_parts[2])):
					sam_name=read_name+'.'+str(rec2)
					rec2+=1
					if chopped_parts[2][x][0]==chopped_parts[2][x][1]:
						if not chopped_parts[0][x]=='':
							print >>fo1,'>'+sam_name
							for y in chop_x_for_fasta(chopped_parts[0][x]):
								print >>fo1,y
						else:
							print >>rec_final,'\t'.join([str(i) for i in chopped_parts[2][x]+['DEL',sam_name]])
					else:
						new_reads=[sam_name]+pin[1:3]+[str(chopped_parts[2][x][0]),pin[4],''.join([''.join(i) for i in chopped_parts[3][1]])]+pin[6:9]+[pin[9][chopped_parts[1][x][0]:chopped_parts[1][x][1]],pin[10][chopped_parts[1][x][0]:chopped_parts[1][x][1]]]+pin[11:-1]+[pin[-1][chopped_parts[1][x][0]:chopped_parts[1][x][1]]]
						print >>fo_final, '\t'.join(new_reads)
	fin.close()
	fo1.close()

def sam_file_chop_up(filein,num_records):
	rec=0
	fin=open(filein)
	header=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='@': 
			break
		header.append(pin)
	fin.close()
	recs=[]
	fin=open(filein)
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='@': 
			if len(recs)<num_records:
				recs.append(pin)
			else:
				fo=open(filein.replace('.sam','.'+str(rec)+'.sam'),'w')
				for x in header:
					print >>fo, '\t'.join(x)
				for x in recs:
					print >>fo, '\t'.join(x)
				fo.close()
				recs=[]
				rec+=1
	fin.close()
	return rec

def file_start(file_out):
	fo=open(file_out,'w')
	fo.close()

def write_sam_header(file_in,file_out):
	fin=open(file_in)
	fo=open(file_out,'w')
	for line in fin:
		pin=line.strip().split('\t')
		if not pin[0][0]=='@':
			break
		else:
			print >>fo, '\t'.join(pin)
	fin.close()
	fo.close()

ref='/mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta'
filein=dict_opts['-i']
file_final_out=filein.replace('.sam','.sr.realigned.sam')
rec_final_out=filein.replace('.sam','.sr.realigned.rec')
write_sam_header(filein,file_final_out)
fo_final=open(file_final_out,'a')
rec_final=open(rec_final_out,'w')
#num_record=2000
#num_sam=sam_file_chop_up(filein,num_record)
sam_file_in=filein
fa_out=sam_file_in.replace('.sam','.sr.extracted.fasta')
sa_out=sam_file_in.replace('.sam','.sr.extracted.sam')
rec_out=sam_file_in.replace('.sam','.sr.extracted.rec')
file_start(fa_out)
reads_read_in(sam_file_in,fa_out,rec_out)
while True:
	f1=os.popen(r'''wc -l %s'''%(fa_out))
	p1=f1.readline().strip().split()
	f1.close()
	if int(p1[0])==0:
		break
	else:
		print int(p1[0])
	os.system(r'''blasr %s %s -sam -sa %s -out %s -bestn 2 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping soft'''%(fa_out,ref,ref+'.sa',sa_out))
	file_start(fa_out)
	reads_read_in(sa_out,fa_out,rec_out)

fo_final.close()
rec_final.close()
