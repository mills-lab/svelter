#!/usr/bin/env python

#!Python
#this code is used to build the null model of insert length distribution (based on bimodal distribution), read depth distribution (based on NB dist, 100bp/bin, corrected by GC content), number of clipped reads, and number of read pairs with abnormai directions(FF). reads to build the null morday would come from cn2 regions larger than a customized size(1kb as default.)
#output data structures: folder named "InputDataSet", containing all input bam file;    folder "NullModel", containing four input folders: 'InsertLengthDist','ReadDepthDist','SplitReadsDist','AbnormalDirectionDist', 'ThroughBreakPointDist'
#Usage:
#SVelter0.Ref.Index.py --ref ref --ppre workdir --ex exclude.bed
#For debug use only
#command='/mnt/EXT/Mills-data/xuefzhao/SVelter/SVelter0.Ref.Index.py --ppre /mnt/EXT/Mills-data/xuefzhao/SVelter/ --ref /mnt/EXT/Mills-data/xuefzhao/SVelter/reference/genome.fa --ex /mnt/EXT/Mills-data/xuefzhao/SVelter/tools/Exclude.bed -s /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/BamFiles/NA12878_S1.bam --chr chr1'
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
def File_Path_Modify(dict_opts):
	if not dict_opts.has_key('-p'):
		dict_opts['-p']='./'
	elif not dict_opts['-p'][-1]=='/':
		dict_opts['-p']+='/'
	return dict_opts

def Read_Block_From_Position(bps,letters,bounce_length,position,end):
	bps2=[int(i) for i in bps]
	relative_bps=[i-min(bps2) for i in bps2]
	if end=='right':
			if position-bounce_length<relative_bps[0]:
				return(['left']+[int(relative_bps[0])-flank,int(relative_bps[0])])
			elif position-bounce_length>relative_bps[-1]:
				return(['right']+[int(relative_bps[-1]),int(relative_bps[-1])+flank])
			else:
				for j in range(len(letters)):
					if position-bounce_length in range(relative_bps[j],relative_bps[j+1]+1):
	#					print [letters[j]]+[relative_bps[j]]+[relative_bps[j+1]]
						return [letters[j]]+[relative_bps[j]]+[relative_bps[j+1]]
	if end=='left':
			if position+bounce_length<relative_bps[0]:
				return(['left']+[int(relative_bps[0])-flank,int(relative_bps[0])])
			elif position+bounce_length>relative_bps[-1]:
				return(['right']+[int(relative_bps[-1]),int(relative_bps[-1])+flank])
			else:
				for j in range(len(letters)):
					if position+bounce_length in range(relative_bps[j],relative_bps[j+1]+1):
						return([letters[j]]+[relative_bps[j]]+[relative_bps[j+1]])

def Pick_CN2_Region(CN2_File,bamF,Length_Limit):
	fi=os.popen('samtools view -h %s'%(bamF))
	Chromosome=[]
	while True:
		pi=fi.readline().strip().split()
		if not pi[0][0]=='@': break
		if pi[0]=='@SQ':
			Chromosome.append(pi[1].split(':')[1])
	CN2_Region={}
	for chrom in Chromosome:
		CN2_Region[chrom]=[]
	fcn2=open(CN2_File)
	while True:
		pcn2=fcn2.readline().strip().split()
		if not pcn2: break
		if not int(pcn2[2])-int(pcn2[1])<Length_Limit:
			CN2_Region[pcn2[0]].append(pcn2)
			#print pcn2
		elif int(pcn2[2])-int(pcn2[1])<Length_Limit: continue 
	return [Chromosome,CN2_Region]

def Fasta_To_Sequence(fasta_file):
	ffa=open(fasta_file)
	ffa.readline().strip().split()
	sequence=[]
	while True:
		seq=ffa.readline().strip().split()
		if not seq: break
		sequence.append(seq)
	seq2=sequence[0][0]
	for seq3 in sequence[1:]:
		seq2=''.join([seq2,seq3[0]])
	if 'N' in seq2:
		return 'ERROR!'
	else:
		return seq2

def GC_Content_Calculate(seq2):
	NumAT=0
	NumGC=0
	for n in seq2:
		if n=='A' or n=='T':
			NumAT+=1
		elif n=='G' or n=='C':
			NumGC+=1
	return [NumGC,NumAT]			

def Region_Coverage_Calculate(sam_file,Number_Of_Windows,Region_Info):
    #sam file here has no header!
    fsam=open(sam_file)
    coverage={}
    num_of_wind=(int(Region_Info[2])-int(Region_Info[1]))/100
    for i in range(num_of_wind):
        coverage[i]=[int(Region_Info[1])+i*100,int(Region_Info[1])+i*100+99,0]
    coverage[-1]=[0,0,0]
    while True:
        psam=fsam.readline().strip().split()
        if not psam: break
        Read_Region=[int(psam[3]),int(psam[3])+len(psam[9])]
        left_block=(Read_Region[0]-int(Region_Info[1]))/100
        if left_block in coverage.keys():
            left_length=coverage[left_block][1]-Read_Region[0]+1
        else:
            left_length=0
        right_block=(Read_Region[1]-int(Region_Info[1]))/100
        if right_block in coverage.keys():
            right_length=Read_Region[1]-coverage[right_block][0]
        else:
            right_length=0
        if left_block<0:
            left_block=-1
        if right_block<0:
            right_block=-1
        if left_block>max(coverage.keys()):
            left_block=-1
        if right_block>max(coverage.keys()):
            right_block=-1
        if left_block==right_block:
            coverage[left_block][-1]+=Read_Region[1]-Read_Region[0]
        else:
            coverage[left_block][-1]+=left_length
            coverage[right_block][-1]+=right_length
            if right_block-left_block>1:
                for k in range(left_block+1,right_block):
                    coverage[k][-1]+=100
    del coverage[-1]
    for k1 in coverage.keys():
        coverage[k1][-1]=float(coverage[k1][-1])/100
    return coverage	

def cigar2reaadlength(cigar):
	pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
	cigars=[]
	for m in pcigar.finditer(cigar):
		cigars.append((m.groups()[0],m.groups()[1]))
	MapLen=0
	for n in cigars:
		if n[1]=='M' or n[1]=='D' or n[1]=='N':
			MapLen+=int(n[0])
	return MapLen

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

def cigar2split(cigar):
	pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
	cigars=[]
	for m in pcigar.finditer(cigar):
		cigars.append((m.groups()[0],m.groups()[1]))
	MapLen=[]
	maplen=0
	if cigars[0][1]=='S':
		MapLen.append(0)
	for n in cigars[1:]:
		if n[1]=='M' or n[1]=='D' or n[1]=='N':
			maplen+=int(n[0])
		if n[1]=='S': 
			MapLen.append(maplen-1)
	return MapLen

def cigar2splitlength(cigar):
	pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
	cigars=[]
	for m in pcigar.finditer(cigar):
		cigars.append((m.groups()[0],m.groups()[1]))
	splitlen=[]
	for i in cigars:
		if i[1]=='S':
			splitlen.append(int(i[0]))
	return splitlen

def path_modify(path):
    if not path[-1]=='/':
        path+='/'
    return path

def ex_file_read_in():
    if '--ex' in dict_opts.keys():
        ex_file=dict_opts['--ex']
    else:
        ex_file=workdir+'tools/Exclude.bed'
    return ex_file

def Gap_Hash_Ref1_read_in(Gap_Refs):
    Gap_Hash_Ref1={}
    for Gap_Ref1 in Gap_Refs:
        fgap=open(Gap_Ref1)
        for line in fgap:
            pgap=line.strip().split()
            if not pgap[0] in Gap_Hash_Ref1.keys():
                Gap_Hash_Ref1[pgap[0]]=[]
            Gap_Hash_Ref1[pgap[0]].append(pgap[1:4])
        fgap.close()
    return Gap_Hash_Ref1

def write_ExcludeBed(ExcludeBed):
    if not os.path.isfile(ExcludeBed):
        fo=open(ExcludeBed,'w')
        for chr_ex in chromos:
            print >>fo, ' '.join([chr_ex,'0','0'])
        fo.close()

def final_regions_decide(Gap_Hash_Ref1,hash_Cor,chrom):
    GapHash=Gap_Hash_Ref1[chrom]
    GapHash1={}
    for key1 in GapHash:
        GapHash1[int(key1[0])]=[int(key1[0]),int(key1[1])]
    GapHash=[]
    for key1 in sorted(GapHash1.keys()):
        GapHash.append(GapHash1[key1])
    for key1 in range(len(GapHash)-1):
        if GapHash[key1+1][0]<GapHash[key1][1]:
            GapHash[key1]=[min(GapHash[key1]+GapHash[key1+1]), max(GapHash[key1]+GapHash[key1+1])]
            GapHash[key1+1]=GapHash[key1]
    for key1 in GapHash:
        if GapHash.count(key1)>1:
            for key2 in range(GapHash.count(key1)-1):
                del GapHash[GapHash.index(key1)]
    flag1=0
    flag2=0
    Cor_Gap={}
    while True:
        if flag1==len(GapHash) or flag2==len(hash_Cor): break
        elif flag2==len(hash_Cor)-1 and not hash_Cor[flag2][1]>GapHash[flag1][0]: break
        else:
            if not  GapHash[flag1][0]>hash_Cor[flag2][1] and not GapHash[flag1][1]<hash_Cor[flag2][0]:
                if not flag2 in Cor_Gap.keys():
                    Cor_Gap[flag2]=[]
                Cor_Gap[flag2].append(flag1)
                flag1+=1
            elif GapHash[flag1][0]>hash_Cor[flag2][1]:
                flag2+=1
            elif GapHash[flag1][1]<hash_Cor[flag2][0]:
                flag1+=1
    Cor_Gap2=[]
    for key1 in range(len(hash_Cor)):
        if not key1 in Cor_Gap.keys():
            Cor_Gap2.append(hash_Cor[key1])
        else:
            CGapTemp=[hash_Cor[key1][0]]
            for key2 in Cor_Gap[key1]:
                CGapTemp+=GapHash[key2]
            CGapTemp.append(hash_Cor[key1][1])
            if CGapTemp[-2]>CGapTemp[-1]:
                del CGapTemp[-1]
                del CGapTemp[-1]
            for key3 in range(len(CGapTemp)/2):
                Cor_Gap2.append([CGapTemp[2*key3],CGapTemp[2*key3+1]])
    return Cor_Gap2

#Default Settings:
script_name=sys.argv[0]
opts,args=getopt.getopt(sys.argv[1:],'h:c:s:',['help=','ppre=','PathBam=','chr=','SamplePbs=', 'ref=','NullSplitLength=','NullILCff=','NullSPCff=','NullDRCff=','NullBed=','NullBedLength=','NullSampleLength=','NullSampleNumber=','ex=','IncludeBed=','ToolMappingQ=','FileMappingQ=','NullSamplePercentage=','SplitLength=','BPSPCff=','BPLNCff=','BPAlignQC=','BPAlignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)
Code_path='/'.join(sys.argv[0].split('/')[:-1])+'/'

time1=time.time()
if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
    print 'SVelter-0.1          Last Update:2014-08-20'
    print ''
    print 'Required Parameters:'
    print '--ppre, workind directory of SVelter, eg: .../SVelter/' 
    print '--ref, absolute path of reference genome. eg: .../SVelter/reference/genome.fa'
    print 'Optional Parameters'
    print '--ex, absolute path of bed file indicating regions to be excluded from analysis. If not provided, no mappable regions will be excluded'
    print '--chr, target chromosome could be specified if only partial genome is interested'
else:
    if not '--ppre' in dict_opts.keys():
        print 'Error: please specify working directory using --ppre'
    else:
        workdir = path_modify(dict_opts['--ppre'])
        ref_file=0
        if not '--ref' in dict_opts.keys():
            print 'Error: please specify refrence genome using --ref'
        else:
            ref_file=dict_opts['--ref']
            ref_path='/'.join(ref_file.split('/')[:-1])+'/'
            ref_index=ref_file+'.fai'
            if not os.path.isfile(ref_index):
                print 'Error: reference genome not indexed'
            else:
                if not workdir in ref_file:
                    ref_path=workdir+'reference/'
                    if not os.path.isdir(ref_path):
                        os.system(r'''mkdir %s'''%(ref_path))
                    os.system(r'''ln -s %s %s'''%(ref_file,ref_path))
                    os.system(r'''ln -s %s %s'''%(ref_index,ref_path))
                    ref_file=ref_path+ref_file.split('/')[-1]
                    ref_index=ref_file+'.fai'
                chromos=[]
                Chromo_Length={}
                fin=open(ref_index)
                for line in fin:
                    pin=line.strip().split()
                    chromos.append(pin[0])
                    Chromo_Length[pin[0]]=int(pin[1])
                fin.close()
                if '--chr' in dict_opts.keys():
                    chromos=[dict_opts['--chr']]
                ExcludeBed=ex_file_read_in()
                write_ExcludeBed(ExcludeBed)
                Gap_Refs=[ExcludeBed]
                Gap_Hash_Ref1=Gap_Hash_Ref1_read_in(Gap_Refs)
                Gap_Hash={}
                for chr_ex in chromos:
                    Gap_Hash[chr_ex]=[]
                for chrom in chromos:
                    fout_Name='.'.join(ref_file.split('.')[:-1])+'.Mappable.'+chrom+'.bed'
                    fout_N2='.'.join(ref_file.split('.')[:-1])+'.GC_Content.'+chrom
                    if not os.path.isfile(fout_Name):
                        fref=os.popen(r'''samtools faidx %s %s:'''%(ref_file,chrom))
                        fref.readline().strip().split()
                        while True:
                            pref=fref.readline().strip().split()
                            if not pref:break
                            Gap_Hash[chrom].append(pref[0])
                        fref.close()
                        #print fout_Name
                        fout=open(fout_Name,'w')
                        #print fout_N2
                        fout2=open(fout_N2,'w')
                        hash_key=chrom
                        if not Gap_Hash[hash_key]==[]:
                            hash_to_Seq=''.join(Gap_Hash[hash_key])
                            hash_Cor=[]
                            hash_cal=0
                            if not hash_to_Seq[0] in ['N','n']:
                                hash_Cor.append([0])
                            for hts in hash_to_Seq[1:]:
                                hash_cal+=1
                                if hts=='N' or hts=='n':
                                    if len(hash_Cor)==0: continue
                                    else:
                                        if len(hash_Cor[-1])==1:
                                            hash_Cor[-1].append(hash_cal)
                                else:
                                    if len(hash_Cor)==0:
                                        hash_Cor.append([hash_cal])
                                    elif len(hash_Cor[-1])==2:
                                        hash_Cor.append([hash_cal])
                            if len(hash_Cor[-1])==1:
                                hash_Cor[-1].append(hash_cal)
                            hakey2=hash_key
                            Cor_Gap2=final_regions_decide(Gap_Hash_Ref1,hash_Cor,chrom)
                            for key1 in Cor_Gap2:
                                if key1[1]-key1[0]>1000:
                                    print >>fout, ' '.join([hakey2, str(key1[0]),str(key1[1])])
                                    print >>fout2, ' '.join([hakey2, str(key1[0]),str(key1[1])])
                                    GC_out=[]
                                    for key2 in range((key1[1]-key1[0])/100):
                                        GC_region=hash_to_Seq[(key1[0]+key2*100):(key1[0]+(key2+1)*100)]
                                        GC_out.append(str(float(GC_region.count('g')+GC_region.count('G')+GC_region.count('c')+GC_region.count('C'))/100.0))
                                    GC_region=hash_to_Seq[(key1[0]+(key2+1)*100):key1[1]]
                                    if len(GC_region)>0:
                                        GC_out.append(str(float(GC_region.count('g')+GC_region.count('G')+GC_region.count('c')+GC_region.count('C'))/float(len(GC_region))))
                                    print >>fout2, ' '.join(GC_out)
                        fout.close()
                        fout2.close()
                    fout=open(fout_Name)
                    total_out={}
                    for line in fout:
                        pout=line.strip().split()
                        if not pout[0] in total_out:
                            total_out[pout[0]]=[]
                        total_out[pout[0]].append(pout[1:3])
                    fout.close()

time2=time.time()
print 'time consuming:'+str(time2-time1)
