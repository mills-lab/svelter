#!/usr/bin/python

#!Python
#this code is used to build the null model of insert length distribution (based on bimodal distribution), read depth distribution (based on NB dist, 100bp/bin, corrected by GC content), number of clipped reads, and number of read pairs with abnormai directions(FF). reads to build the null morday would come from cn2 regions larger than a customized size(1kb as default.)
#output data structures: folder named "InputDataSet", containing all input bam file;    folder "NullModel", containing four input folders: 'InsertLengthDist','ReadDepthDist','SplitReadsDist','AbnormalDirectionDist', 'ThroughBreakPointDist'
#Usage:
#SVelter1.NullModel.py --ppre workdir -s sample.bam --ref ref.fa --cn cn2.bed
#For debug use only
#command='python SVelter1.NullModel.py --NullModel Simple --ppre wirkdir -s sample.bam --ref genome.fa --cn2 cn2_files'
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
        seq2+=seq3[0]
        if ''.join(['N' for x in range(100)]) in seq2: break
    if ''.join(['N' for x in range(100)]) in seq2:
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

#Default Settings:
opts,args=getopt.getopt(sys.argv[1:],'o:h:s:',['NullModel=','cn=','ex=','help=','ppre=','PathBam=','sample=','chr=','SamplePbs=', 'ref=','NullSplitLength=','NullILCff=','NullSPCff=','NullDRCff=','NullBedLength=','NullSampleLength=','NullSampleNumber=','ExcludeBed=','IncludeBed=','ToolMappingQ=','FileMappingQ=','NullSamplePercentage=','SplitLength=','BPSPCff=','BPLNCff=','BPAlignQC=','BPAlignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)
if not '--NullModel' in dict_opts.keys():
    model_comp='S'
else:
    if dict_opts['--NullModel'] in ['S','Simple']:
        model_comp='S'
    else:
        model_comp='C'

if '--ReadLength' in dict_opts.keys():
    ReadLength=int(dict_opts['--ReadLength'])
    ReadLength_Flag=1
else:
    ReadLength=0
    ReadLength_Flag=0

if '--QCAlign' in dict_opts.keys():
    QCAlign=int(dict_opts['--QCAlign'])
else:
    QCAlign=20

if '--QCSplit' in dict_opts.keys():
    QCSplit=int(dict_opts['--QCSplit'])
else:
    QCSplit=20

if '--NullSplitLength' in dict_opts.keys():
    NullSplitLen_perc=int(dict_opts['--NullSplitLength'])
else:
    NullSplitLen_perc=0.9

if '--NullILCI' in dict_opts.keys():
    NullILCIs=dict_opts['--NullILCI']
else:
    NullILCIs=[0.025,0.05,0.95,0.975]

if '--NullRDCI' in dict_opts.keys():
    NullRDCIs=dict_opts['--NullRDCI']
else:
    NullRDCIs=[0.025,0.05,0.95,0.975]

if '--NullTBCI' in dict_opts.keys():
    NullTBCIs=dict_opts['--NullTBCI']
else:
    NullTBCIs=[0.0001,0.0005,0.9999,0.9995]

if '--NullILCff'in dict_opts.keys():
    NullILCff=dict_opts['--NullILCff']
else:
    NullILCff=0.999

if '--NullSPCff' in dict_opts.keys():
    NullSPCff=dict_opts['--NullSPCff']
else:
    NullSPCff=0.999

if '--NullDRCff' in dict_opts.keys():
    NullDRCff=dict_opts['--NullDRCff']
else:
    NullDRCff=0.999

if '--KeepFile' in dict_opts.keys():
    KeepFile=dict_opts['--KeepFile']
else:
    KeepFile='Yes'

if '--KeepFigure' in dict_opts.keys():
    KeepFigure=dict_opts['--KeepFigure']
else:
    KeepFigure='No'

if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
    print 'SVelter-0.1          Last Update:2014-08-20'
    print ''
    print 'Required Parameters:'
    print '--ppre, workind directory of SVelter, eg: .../SVelter/' 
    print '--ref, absolute path of reference genome. eg: .../SVelter/reference/genome.fa'
    print '-s, absolute path of input bam file. eg: .../SVelter/BamFiles/sample.bam'
    #print '--PathBam, absolute path of folder containing all input bam files. eg: .../SVelter/BamFiles/'
    print ''
    print 'Optional Parameters:'
    print '--ex, bed file containing genomic regions to exclude. eg: Excludable.bed'
    print '--cn, bed file containing genomic regions where no SVs have yet been reported, based on which null model were built'
    print '--NullModel, specify which stat model to be fitted on each parameter. if --NullModel==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
    print '--NullBedLength, minimum requirement for length of regions used to build null model; default: 2kb'
    print '--NullSamplePercentage, sampling percentage of regions from Nullbed'
    print '--NullSampleLength, if not --NullBed provided, SVelter will randomly sample regions on genome to build null model; --NullSampleLength specify the length of region picked; default: 5kb'
    print '--NullSampleNumber, if not --NullBed provided, SVelter will randomly sample regions on genome to build null model; --NullSampleNumber specify the number of region picked; default: 10k'
else:
    if not '--ppre' in dict_opts.keys():
        print 'Error: please specify working directory using: --ppre'
    else:
        workdir=dict_opts['--ppre']
        if not workdir[-1]=='/':
            workdir+='/'
        if not '-s' in dict_opts.keys():
            print 'Error: please specify either input file using -s'
        else:
            #if '-s' in dict_opts.keys():
            bam_path='/'.join(dict_opts['-s'].split('/')[:-1])+'/'
            bam_files=[dict_opts['-s']]
            #else:
            #    bam_path=dict_opts['--PathBam']
            #    if not bam_path[-1]=='/':
            #        bam_path+='/'
            #    bam_files=[]
            #    for file in os.listdir(bam_path):
            #        if file.split('.')[-1]=='bam':
            #            bam_files.append(bam_path+file)
            if not '--ref' in dict_opts.keys():
                print 'Error: please specify refrence genome using --ref'
            else:
                ref_path='/'.join(dict_opts['--ref'].split('/')[:-1])+'/'
                ref_file=dict_opts['--ref']
                ref_index=dict_opts['--ref']+'.fai'
                if not os.path.isfile(ref_index):
                    print 'Error: reference genome not indexed'
                else:
                    if '--NullSamplePercentage' in dict_opts.keys():
                        SamplingPercentage=float(dict_opts['--NullSamplePercentage'])
                    else:
                        SamplingPercentage=0.1
                    chromos=[]
                    fin=open(ref_index)
                    for line in fin:
                        pin=line.strip().split()
                        chromos.append(pin[0])
                    fin.close()
                    if '--cn' in dict_opts.keys():
                        cn2_file=dict_opts['--cn']
                        if '--NullBedLength' in dict_opts.keys():
                            cn2_length=int(dict_opts['--NullBedLength'])
                        else:
                            cn2_length=2000
                    else:
                        if '-x' in dict_opts.keys():
                            cn2_path=dict_opts['-x']
                            if not cn2_path[-1]=='/':
                                cn2_path+='/'
                        else:
                            cn2_path=workdir+'SVelter.tools/'
                        if not os.path.isdir(cn2_path):
                            os.system(r'''mkdir %s'''%(cn2_path))
                        if not 'chr' in chromos[0]:
                            cn2_file=cn2_path+'CN2.g1k.bed'
                        else:
                            cn2_file=cn2_path+'CN2.hg19.bed'
                        if os.path.isfile(cn2_file):
                            if '--NullBedLength' in dict_opts.keys():
                                cn2_length=int(dict_opts['--NullBedLength'])
                            else:
                                cn2_length=2000
                        else:
                            if not '--NullSampleLength' in dict_opts.keys():
                                dict_opts['--NullSampleLength']=5000
                            else:
                                dict_opts['--NullSampleLength']=int(dict_opts['--NullSampleLength'])
                            if not '--NullSampleNumber' in dict_opts.keys():
                                dict_opts['--NullSampleNumber']=10000
                            else:
                                dict_opts['--NullSampleNumber']=int(dict_opts['--NullSampleNumber'])
                            cn2_length=dict_opts['--NullSampleLength']-100
                            whole_genome={}
                            fref=open(ref_index)
                            for line in fref:
                                pref=line.strip().split()
                                whole_genome[pref[0]]=[int(pref[1])]
                            fref.close()
                            len_genome=0
                            for i in whole_genome.keys():
                                len_genome+=whole_genome[i][0]
                            fo=open(cn2_file,'w')
                            for i in whole_genome.keys():
                                num_i=int(float(whole_genome[i][0])/float(len_genome)*dict_opts['--NullSampleNumber'])
                                reg_i=[random.randint(1,whole_genome[i][0]-dict_opts['--NullSampleLength']) for j in range(num_i)]
                                for j in sorted(reg_i):
                                    print >>fo, ' '.join([i,str(j),str(j+dict_opts['--NullSampleLength']-1)])
                            fo.close()
                            SamplingPercentage=1
                    if not '--chr' in dict_opts.keys():
                        chr_flag=0
                    elif '--chr' in dict_opts.keys():
                        chr_flag=1
                        chrom_single=dict_opts['--chr']
                        if 'chr' in chrom_single and not 'chr' in chromos[0]:
                            if not chrom_single.replace('chr','') in chromos:
                                chromo2=[]
                            else:
                                chromo2=[chrom_single.replace('chr','')]
                        if not 'chr' in chrom_single and 'chr' in chromos[0]:
                            if not 'chr'+chrom_single in chromos:
                                chromo2=[]
                            else:
                                chromo2=['chr'+chrom_single]
                        else:
                            if not chrom_single in chromos:
                                chromo2=[]
                            else:
                                chromo2=[chrom_single]
                        chromos=chromo2
                    if chromos==[]:
                        print 'Error: --chr chromosome not in reference'
                    else:
                        if not '--chr' in dict_opts.keys():
                            genome_name='genome'
                        else:
                            genome_name=chromos[0]
                        if not '-o' in dict_opts.keys():
                            out_path=workdir
                        else:
                            out_path=dict_opts['-o']
                            if not out_path[-1]=='/':
                                out_path+='/'
                        NullPath=out_path+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/'
                        if not '--PathSVelter' in dict_opts.keys():
                            #print 'Error: please specify code directory using: --PathSVelter'
                            Script_Path='/'.join(sys.argv[0].split('/')[:-1])
                        else:
                            Script_Path=dict_opts['--PathSVelter']
                        if not Script_Path[-1]=='/':
                            Script_Path+='/'
                        if not os.path.isdir(NullPath):
                            os.system(r'''mkdir %s'''%(NullPath))
                        for bamF in bam_files:
                            if ReadLength_Flag==0:
                                ReadLengthHash={}
                            outputfile=NullPath+bamF.split('/')[-1].replace('.bam','')+'.'+genome_name+'.null'
                            fo=open(outputfile,'w')
                            print >>fo, ' '.join(['position','GCContent','ReadDepth','SplitReads','AbnormalDirection','ThroughBP'])
                            fo.close()
                            SplitLength={}
                            InsertLength={}
                            Ref_Flag=0
                            if 'chr' in chromos[0]:
                                Ref_Flag=1
                            fcn2=open(cn2_file)
                            chr_cn2=[]
                            while True:
                                pcn2=fcn2.readline().strip().split()
                                if not pcn2: break
                                if not len(pcn2)==3: break
                                if pcn2[0] in chromos:
                                    #print pcn2
                                    if not pcn2[0] in chr_cn2:
                                        chr_cn2.append(pcn2[0])
                                    if (int(pcn2[2])-int(pcn2[1]))>cn2_length and (int(pcn2[2])-int(pcn2[1]))<10**6:
                                        if not random.choice(range(100))>SamplingPercentage*100:
                                            freadpairs=os.popen('''samtools view %s %s:%d-%d'''%(bamF,pcn2[0],int(pcn2[1])+100,int(pcn2[2])-100))
                                            while True:
                                                preadpair=freadpairs.readline().strip().split()
                                                if not preadpair: break
                                                if not int(preadpair[4])>QCAlign: continue
                                                if preadpair[8]=='0': continue
                                                if not abs(int(preadpair[8])) in InsertLength.keys():
                                                    InsertLength[abs(int(preadpair[8]))]=1
                                                else: 
                                                    InsertLength[abs(int(preadpair[8]))]+=1
                                                if 'S' in preadpair[5]:
                                                    SplitLen=cigar2splitlength(preadpair[5])
                                                    for s in SplitLen:
                                                        if not s in SplitLength.keys():
                                                            SplitLength[s]=1
                                                        else:
                                                            SplitLength[s]+=1
                            #calculate the minimum requirement for the length of soft clip 
                            if not chr_cn2==[]:
                                SplitLengthPath=NullPath+'SplitLength'
                                if not os.path.isdir(SplitLengthPath):
                                        os.system(r'''mkdir %s'''%(SplitLengthPath))
                                SplitLengthOutput=SplitLengthPath+'/'+bamF.split('/')[-1].replace('.bam','')+'.'+genome_name+'.SplitLength'
                                #print SplitLengthOutput
                                fslo=open(SplitLengthOutput,'w')
                                print >> fslo, ' '.join(['Length_of_Split_Region', 'Time_of_Observation'])
                                Total_Split_Reads=0
                                for key in sorted(SplitLength.keys()):
                                        print >> fslo, ' '.join([str(key),str(SplitLength[key])])
                                        Total_Split_Reads+=SplitLength[key]
                                fslo.close()
                                Sub_Split_Reads=0
                                if NullSplitLen_perc>1:
                                    SplitLenPNum=int(NullSplitLen_perc)
                                elif NullSplitLen_perc<1:
                                    for key in sorted(SplitLength.keys()):
                                        Sub_Split_Reads+=SplitLength[key]
                                        if float(Sub_Split_Reads)/float(Total_Split_Reads) > NullSplitLen_perc:
                                            SplitLenPNum=key
                                            break
                                #>>> SplitLenPNum=8
                                #For Insert Length CI search: bimodal distribution
                                ILMedian=numpy.median(InsertLength.keys())
                                for key in InsertLength.keys():
                                        if key > ILMedian*10 or key==0:
                                                del InsertLength[key]
                                TotalILNum=0
                                for key in InsertLength.keys():
                                        TotalILNum+=InsertLength[key]
                                NullILCILeft=[]
                                NullILCIRight=[]
                                SubILNumleft=0
                                SubILNumright=0
                                NciLeft=0
                                NciRight=0
                                for keyn in range(len(sorted(InsertLength.keys()))):
                                        SubILNumleft+=InsertLength[sorted(InsertLength.keys())[keyn]]
                                        SubILNumright+=InsertLength[sorted(InsertLength.keys())[-(keyn+1)]]
                                        if NciLeft<len(NullILCIs)/2:
                                                if SubILNumleft<NullILCIs[NciLeft]*float(TotalILNum): continue
                                                if not SubILNumleft<NullILCIs[NciLeft]*float(TotalILNum): 
                                                        if len(NullILCILeft)==NciLeft:
                                                                NullILCILeft.append(sorted(InsertLength.keys())[keyn])
                                                                NciLeft+=1
                                                                #print NciLeft
                                        if NciRight<(len(NullILCIs)/2):
                                                if SubILNumright<NullILCIs[NciRight]*float(TotalILNum): continue
                                                if not SubILNumright<NullILCIs[NciRight]*float(TotalILNum):
                                                        if len(NullILCIRight)==NciRight:
                                                                NullILCIRight.append(sorted(InsertLength.keys())[-(keyn+1)])
                                                                NciRight+=1
                                                                #print NciRight
                                        if NciLeft==len(NullILCIs)/2 and NciRight==len(NullILCIs)/2: break
                                NullILCI=NullILCILeft+sorted(NullILCIRight)
                                ILNullDensity={}
                                RDNullDensity={}
                                DRNullDensity={}
                                TBNullDensity={}
                                SplitNullDensity={}
                                GC_Content={}
                                fcn2=open(cn2_file)
                                while True:
                                    pcn2=fcn2.readline().strip().split()
                                    if not pcn2: break
                                    if not len(pcn2)==3: break
                                    if pcn2[0] in chromos:
                                        if (int(pcn2[2])-int(pcn2[1]))>cn2_length and (int(pcn2[2])-int(pcn2[1]))<10**6:
                                            if random.choice(range(100))<SamplingPercentage*100:
                                                #print pcn2
                                                pos=range(int(pcn2[1])+100,(int(pcn2[2])+1)-100)
                                                pos0=[pcn2[0]+'_'+str(i) for i in pos]
                                                RDNull=[0 for i in pos]
                                                SplitNull=[0 for i in pos]
                                                DRNull=[0 for i in pos]
                                                TBNull=[0 for i in pos]
                                                ILNull=[0 for i in pos]
                                                readInf=[]
                                                freadpairs=os.popen('''samtools view %s %s:%d-%d'''%(bamF,pcn2[0],int(pcn2[1])+100,int(pcn2[2])-100))
                                                while True:
                                                        preadpair=freadpairs.readline().strip().split()
                                                        if not preadpair: break
                                                        if not int(preadpair[4])>QCAlign: continue
                                                        if ReadLength_Flag==0:
                                                            if not preadpair[9]=='*':
                                                                if not len(preadpair[9]) in ReadLengthHash.keys():
                                                                    ReadLengthHash[len(preadpair[9])]=1
                                                                else:
                                                                    ReadLengthHash[len(preadpair[9])]+=1
                                                        for i in range(max(int(preadpair[3])-pos[0],0),min(int(preadpair[3])+cigar2reaadlength(preadpair[5])-pos[0], pos[-1]-pos[0])):
                                                                RDNull[i]+=1
                                                        if int(preadpair[8])<NullILCI[0] or int(preadpair[8])>NullILCI[-1]:
                                                                if int(preadpair[3]) in pos: 
                                                                        ILNull[int(preadpair[3])-pos[0]]+=1
                                                        if not preadpair[5].find('S')==-1:
                                                                splitpos=[i+int(preadpair[3]) for i in cigar2split(preadpair[5])]
                                                                splitlent=cigar2splitlength(preadpair[5])
                                                                for j in range(len(splitpos)):
                                                                        if not splitlent[j]<SplitLenPNum and splitpos[j] in pos:
                                                                                SplitNull[splitpos[j]-pos[0]]+=1
                                                        if preadpair[6]=='=':
                                                                if not Reads_Direction_Detect(preadpair[1])==['+', '-']:
                                                                        abdrpos=min([int(preadpair[3]),int(preadpair[7])])-int(pcn2[1])-100
                                                                        if abdrpos>-1 and abdrpos<(int(pcn2[2])-int(pcn2[1])-200):
                                                                                DRNull[abdrpos]+=1
                                                                if not preadpair[0] in readInf:
                                                                        for j in range(max(int(preadpair[3])-pos[0],0), min(int(preadpair[3])+int(preadpair[8])-pos[0],pos[-1]-pos[0])):
                                                                                TBNull[j]+=1
                                                                        readInf.append(preadpair[0])
                                                fref=os.popen('''samtools faidx %s %s:%d-%d'''%(dict_opts['--ref'],pcn2[0],int(pcn2[1])+100,int(pcn2[1])+(int(pcn2[2])-int(pcn2[1]))/100*100+99-100))
                                                tref=fref.readline().strip().split()
                                                REFSEQUENCE=fref.readline().strip().split()
                                                while True:
                                                        pref=fref.readline().strip().split()
                                                        if not pref: break
                                                        REFSEQUENCE=[''.join(REFSEQUENCE+pref)]
                                                GCTempContent=[REFSEQUENCE[0][(100*i):(100*i+100)].count('G')+REFSEQUENCE[0][(100*i):(100*i+100)].count('C')+REFSEQUENCE[0][(100*i):(100*i+100)].count('g')+REFSEQUENCE[0][(100*i):(100*i+100)].count('c') for i in range(len(REFSEQUENCE[0])/100)]
                                                fo=open(outputfile,'a')
                                                GCNull=[GCTempContent[0] for i in range(100)]
                                                for j in GCTempContent[1:]:
                                                        GCNull+=[j for i in range(100)]
                                                for k in range(len(pos)):
                                                        if not RDNull[k]==0:
                                                                print >>fo, ' '.join([str(pos0[k]),str(GCNull[k]),str(RDNull[k]),str(SplitNull[k]),str(DRNull[k]),str(TBNull[k])])
                                                for i in range(len(RDNull)):
                                                        if not GCNull[i] in GC_Content.keys():
                                                                GC_Content[GCNull[i]]=[RDNull[i]]
                                                        if GCNull[i] in GC_Content.keys():
                                                                GC_Content[GCNull[i]].append(RDNull[i])
                                                fo.close()
                                                for k in range(len(pos)):
                                                        if not ILNull[k] in ILNullDensity.keys():
                                                                ILNullDensity[ILNull[k]]=1
                                                        elif ILNull[k] in ILNullDensity.keys():
                                                                ILNullDensity[ILNull[k]]+=1
                                                        if not RDNull[k] in RDNullDensity.keys():
                                                                RDNullDensity[RDNull[k]]=1
                                                        elif RDNull[k] in RDNullDensity.keys():
                                                                RDNullDensity[RDNull[k]]+=1
                                                        if not SplitNull[k] in SplitNullDensity.keys():
                                                                SplitNullDensity[SplitNull[k]]=1
                                                        elif SplitNull[k] in SplitNullDensity.keys():
                                                                SplitNullDensity[SplitNull[k]]+=1
                                                        if not DRNull[k] in DRNullDensity.keys():
                                                                DRNullDensity[DRNull[k]]=1
                                                        elif DRNull[k] in DRNullDensity.keys():
                                                                DRNullDensity[DRNull[k]]+=1
                                                for k in range(len(pos)/100)[1:]:
                                                        if not TBNull[k*100] in TBNullDensity.keys():
                                                                TBNullDensity[TBNull[k*100]]=1
                                                        elif TBNull[k*100] in TBNullDensity.keys():
                                                                TBNullDensity[TBNull[k*100]]+=1
                                if 0 in RDNullDensity.keys():
                                        del RDNullDensity[0]
                                if 0 in TBNullDensity.keys():
                                        del TBNullDensity[0]
                                OverallRDDenominator=0
                                OverallRDNumeritor=0
                                if not RDNullDensity=={}:
                                    for key in RDNullDensity.keys():
                                        if not key==0:
                                            OverallRDNumeritor+=key*RDNullDensity[key]
                                            OverallRDDenominator+=RDNullDensity[key]
                                    OverallRDNullDensity=float(OverallRDNumeritor)/float(OverallRDDenominator)
                                    fbRD=open(outputfile)
                                    pbRD=fbRD.readline().strip().split()
                                    RD_Af_Adj={}
                                    for key in GC_Content.keys():
                                            GC_Content[key]=numpy.mean(GC_Content[key])
                                    while True:
                                            pbRD=fbRD.readline().strip().split()
                                            if not pbRD: break
                                            if not len(pbRD)==6 : break
                                            RDAfAdj=int(pbRD[2])*OverallRDNullDensity/GC_Content[int(pbRD[1])]
                                            if not int(RDAfAdj) in RD_Af_Adj.keys():
                                                    RD_Af_Adj[int(RDAfAdj)]=1
                                            elif int(RDAfAdj) in RD_Af_Adj.keys():
                                                    RD_Af_Adj[int(RDAfAdj)]+=1
                                    #For Adjusted Read Depth CI search: negative binomial distribution
                                    RDMedian=numpy.median(RD_Af_Adj.keys())
                                    for key in RD_Af_Adj.keys():
                                            if key > RDMedian*10 or key==0:
                                                    del RD_Af_Adj[key]
                                    TotalRDNum=0
                                    for key in RD_Af_Adj.keys():
                                            TotalRDNum+=RD_Af_Adj[key]
                                    NullRDCILeft=[]
                                    NullRDCIRight=[]
                                    SubRDNumleft=0
                                    SubRDNumright=0
                                    NciLeft=0
                                    NciRight=0
                                    for keyn in range(len(sorted(RD_Af_Adj.keys()))):
                                            SubRDNumleft+=RD_Af_Adj[sorted(RD_Af_Adj.keys())[keyn]]
                                            SubRDNumright+=RD_Af_Adj[sorted(RD_Af_Adj.keys())[-(keyn+1)]]
                                            if NciLeft<len(NullRDCIs)/2:
                                                    if SubRDNumleft<NullRDCIs[NciLeft]*float(TotalRDNum): continue
                                                    if not SubRDNumleft<NullRDCIs[NciLeft]*float(TotalRDNum): 
                                                            if len(NullRDCILeft)==NciLeft:
                                                                    NullRDCILeft.append(sorted(RD_Af_Adj.keys())[keyn])
                                                                    NciLeft+=1
                                                                    #print NciLeft
                                            if NciRight<(len(NullRDCIs)/2):
                                                    if SubRDNumright<NullRDCIs[NciRight]*float(TotalRDNum): continue
                                                    if not SubRDNumright<NullRDCIs[NciRight]*float(TotalRDNum):
                                                            if len(NullRDCIRight)==NciRight:
                                                                    NullRDCIRight.append(sorted(RD_Af_Adj.keys())[-(keyn+1)])
                                                                    NciRight+=1
                                                                    #print NciRight
                                            if NciLeft==len(NullRDCIs)/2 and NciRight==len(NullRDCIs)/2: break
                                    NullRDCI=NullRDCILeft+sorted(NullRDCIRight)
                                    #For Through Break Points CI search: bimodal distribution
                                    TBMedian=numpy.median(TBNullDensity.keys())
                                    for key in TBNullDensity.keys():
                                            if key > TBMedian*10 or key==0:
                                                    del TBNullDensity[key]
                                    TotalTBNum=0
                                    for key in TBNullDensity.keys():
                                            TotalTBNum+=TBNullDensity[key]
                                    NullTBCILeft=[]
                                    NullTBCIRight=[]
                                    SubTBNumleft=0
                                    SubTBNumright=0
                                    NciLeft=0
                                    NciRight=0
                                    for keyn in range(len(sorted(TBNullDensity.keys()))):
                                            SubTBNumleft+=TBNullDensity[sorted(TBNullDensity.keys())[keyn]]
                                            SubTBNumright+=TBNullDensity[sorted(TBNullDensity.keys())[-(keyn+1)]]
                                            if NciLeft<len(NullTBCIs)/2:
                                                    if SubTBNumleft<NullTBCIs[NciLeft]*float(TotalTBNum): continue
                                                    if not SubTBNumleft<NullTBCIs[NciLeft]*float(TotalTBNum): 
                                                            if len(NullTBCILeft)==NciLeft:
                                                                    NullTBCILeft.append(sorted(TBNullDensity.keys())[keyn])
                                                                    NciLeft+=1
                                                                    #print NciLeft
                                            if NciRight<(len(NullTBCIs)/2):
                                                    if SubTBNumright<NullTBCIs[NciRight]*float(TotalTBNum): continue
                                                    if not SubTBNumright<NullTBCIs[NciRight]*float(TotalTBNum):
                                                            if len(NullTBCIRight)==NciRight:
                                                                    NullTBCIRight.append(sorted(TBNullDensity.keys())[-(keyn+1)])
                                                                    NciRight+=1
                                                                    #print NciRight
                                            if NciLeft==len(NullTBCIs)/2 and NciRight==len(NullTBCIs)/2: break
                                    NullTBCI=NullTBCILeft+sorted(NullTBCIRight)
                                    #For Aberrant Insert Length Read Pair:
                                    TotalILNum=0
                                    for key in ILNullDensity.keys():
                                            TotalILNum+=ILNullDensity[key]
                                    ILITotal=0
                                    for key in sorted(ILNullDensity.keys()):
                                            ILITotal+=ILNullDensity[key]
                                            if float(ILITotal)/float(TotalILNum)>NullILCff: break
                                    if sorted(ILNullDensity.keys()).index(key)>4 or len(ILNullDensity.keys())<4:
                                        ILIPoint=sorted(ILNullDensity.keys())[:sorted(ILNullDensity.keys()).index(key)+1]
                                        ILIPoint2=[float(ILNullDensity[i])/float(TotalILNum) for i in ILIPoint]
                                        ILIPoint3=[]
                                        for k in range(len(ILIPoint)):
                                            ILIPoint3.append(sum(ILIPoint2[:(k+1)]))
                                    else:
                                        ILIPoint=sorted(ILNullDensity.keys())[:4]
                                        ILIPoint2=[float(ILNullDensity[i])/float(TotalILNum) for i in ILIPoint]
                                        ILIPoint3=[]
                                        for k in range(len(ILIPoint)):
                                            ILIPoint3.append(sum(ILIPoint2[:(k+1)]))
                                    #For Split Reads point search: 
                                    TotalSplitNum=0
                                    for key in SplitNullDensity.keys():
                                            TotalSplitNum+=SplitNullDensity[key]
                                    SplitITotal=0
                                    for key in sorted(SplitNullDensity.keys()):
                                            SplitITotal+=SplitNullDensity[key]
                                            if float(SplitITotal)/float(TotalSplitNum)>NullSPCff: break
                                    if sorted(SplitNullDensity.keys()).index(key)>4 or len(SplitNullDensity.keys())<4:
                                        SplitIPoint=sorted(SplitNullDensity.keys())[:sorted(SplitNullDensity.keys()).index(key)+1]
                                        SplitIPoint2=[float(SplitNullDensity[i])/float(TotalSplitNum) for i in SplitIPoint]
                                        SplitIPoint3=[]
                                        for k in range(len(SplitIPoint)):
                                            SplitIPoint3.append(sum(SplitIPoint2[:(k+1)]))
                                    else:
                                        SplitIPoint=sorted(SplitNullDensity.keys())[:4]
                                        SplitIPoint2=[float(SplitNullDensity[i])/float(TotalSplitNum) for i in SplitIPoint]
                                        SplitIPoint3=[]
                                        for k in range(len(SplitIPoint)):
                                            SplitIPoint3.append(sum(SplitIPoint2[:(k+1)]))
                                    #For Aberrant Direction Read Pair search: 
                                    TotalDRNum=0
                                    for key in DRNullDensity.keys():
                                            TotalDRNum+=DRNullDensity[key]
                                    DRITotal=0
                                    for key in sorted(DRNullDensity.keys()):
                                            DRITotal+=DRNullDensity[key]
                                            if float(DRITotal)/float(TotalDRNum)>NullDRCff: break
                                    if sorted(DRNullDensity.keys()).index(key)>4 or len(DRNullDensity.keys())<4:
                                        DRIPoint=sorted(DRNullDensity.keys())[:(sorted(DRNullDensity.keys()).index(key)+1)]
                                        DRIPoint2=[float(DRNullDensity[i])/float(TotalDRNum) for i in DRIPoint]
                                        DRIPoint3=[]
                                        for k in DRIPoint:
                                            DRIPoint3.append(sum(DRIPoint2[:(k+1)]))
                                    else:
                                        DRIPoint=sorted(DRNullDensity.keys())[:4]
                                        DRIPoint2=[float(DRNullDensity[i])/float(TotalDRNum) for i in DRIPoint]
                                        DRIPoint3=[]
                                        for k in DRIPoint:
                                            DRIPoint3.append(sum(DRIPoint2[:(k+1)]))
                                    outputfileStat=NullPath+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.Stats'
                                    #print outputfileStat
                                    fos=open(outputfileStat,'w')
                                    print >>fos, 'Insert Lenth CIs'
                                    print >>fos, ' '.join([str(i) for i in NullILCIs])
                                    print >>fos, ' '.join([str(i) for i in NullILCI])
                                    print >>fos, 'Read Depth CIs'
                                    print >>fos, ' '.join([str(i) for i in NullRDCIs])
                                    print >>fos, ' '.join([str(i) for i in NullRDCI])
                                    print >>fos, 'Number of Reads Going Through a Break Points CIs'
                                    print >>fos, ' '.join([str(i) for i in NullTBCIs])
                                    print >>fos, ' '.join([str(i) for i in NullTBCI])
                                    print >>fos, 'Number of Read Pairs with Aberrant Insert Length'
                                    print >>fos, ' '.join([str(i) for i in ILIPoint])
                                    print >>fos, ' '.join([str(i) for i in ILIPoint3])
                                    print >>fos, 'Number of Split Reads'
                                    print >>fos, ' '.join([str(i) for i in SplitIPoint])
                                    print >>fos, ' '.join([str(i) for i in SplitIPoint3])
                                    print >>fos, 'Number of Read Pairs with Aberrant Direction'
                                    print >>fos, ' '.join([str(i) for i in DRIPoint])
                                    print >>fos, ' '.join([str(i) for i in DRIPoint3])
                                    if ReadLength_Flag==0:
                                        ReadLengthTag=0
                                        ReadLengthOut=0
                                        for RLK1 in ReadLengthHash.keys():
                                            if ReadLengthHash[RLK1]>ReadLengthTag:
                                                ReadLengthOut=RLK1
                                                ReadLengthTag=ReadLengthHash[RLK1]
                                        print >>fos, 'Read Length Of Reads'+':'+str(ReadLengthOut)
                                    fos.close()
                                    outputfileIL=NullPath+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.density.null'
                                    #print outputfileIL
                                    foIL=open(outputfileIL,'w')
                                    print >>foIL, ' '.join(['InsertLength','Percentage'])
                                    for s in ILNullDensity.keys():
                                            print >>foIL, ' '.join([str(s), str(ILNullDensity[s])])
                                    print >>foIL, ' '.join(['SplitReads','Percentage'])
                                    for s in SplitNullDensity.keys():
                                            print >>foIL, ' '.join([str(s), str(SplitNullDensity[s])])
                                    print >>foIL, ' '.join(['AbnormalDirection','Percentage'])
                                    for d in DRNullDensity.keys():
                                            print >>foIL, ' '.join([str(d),str(DRNullDensity[d])])
                                    print >>foIL, ' '.join(['InsertLength','Frequency'])
                                    for l in InsertLength.keys():
                                            print >>foIL, ' '.join([str(l),str(InsertLength[l])])
                                    print >>foIL, ' '.join(['ReadDepth','Frequency'])
                                    for r in RD_Af_Adj.keys():
                                            print >>foIL, ' '.join([str(r),str(RD_Af_Adj[r])])
                                    print >>foIL, ' '.join(['ThroughBreakPoint','Frequency'])
                                    for t in TBNullDensity.keys():
                                            print >>foIL, ' '.join([str(t),str(TBNullDensity[t])])
                                    print >>foIL,' '.join(['BinPosition','GC_Content'])
                                    for b in GC_Content.keys():
                                        if not b==0:
                                            print >>foIL, ' '.join([str(b), str(GC_Content[b])])
                                    foIL.close()
                                    RFigureDRSplit=Script_Path+'SVelter1.NullModel.Figure.a.r'	
                                    InsertLenNullTemp=NullPath+'InsertLenNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    #print InsertLenNullTemp
                                    fIL=open(InsertLenNullTemp,'w')
                                    for dr in ILNullDensity.keys():
                                            print >> fIL, ' '.join([str(dr),str(ILNullDensity[dr])])
                                    fIL.close()
                                    InsertLenNullfigure1='.'.join(InsertLenNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    InsertLenNullfigure2='.'.join(InsertLenNullTemp.split('.')[:-1])+'.2.jpg'
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,InsertLenNullTemp,InsertLenNullfigure1,BoxPlotColor,InsertLenNullfigure2))
                                    DRNullTemp=NullPath+'DirectionNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    #print DRNullTemp
                                    fDR=open(DRNullTemp,'w')
                                    for dr in DRNullDensity.keys():
                                            print >> fDR, ' '.join([str(dr),str(DRNullDensity[dr])])
                                    fDR.close()
                                    DRNullfigure1='.'.join(DRNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    DRNullfigure2='.'.join(DRNullTemp.split('.')[:-1])+'.2.jpg'
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,DRNullTemp,DRNullfigure1,BoxPlotColor,DRNullfigure2))
                                    SplitNullTemp=NullPath+'SplitNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    #print SplitNullTemp
                                    fSP=open(SplitNullTemp,'w')
                                    for sp in SplitNullDensity.keys():
                                            print >> fSP, ' '.join([str(sp),str(SplitNullDensity[sp])])
                                    fSP.close()
                                    SplitNullfigure1='.'.join(SplitNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    SplitNullfigure2='.'.join(SplitNullTemp.split('.')[:-1])+'.2.jpg'
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,SplitNullTemp,SplitNullfigure1,BoxPlotColor,SplitNullfigure2))
                                    if model_comp=='C':
                                        RFigureDRSplit2=Script_Path+'SVelter1.NullModel.Figure.d.r'	
                                    else:
                                        RFigureDRSplit2=Script_Path+'SVelter1.NullModel.Figure.d2.r'    
                                    RDNullTemp=NullPath+'RDNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fRD=open(RDNullTemp,'w')
                                    for rd in RD_Af_Adj.keys():
                                            print >> fRD, ' '.join([str(rd),str(RD_Af_Adj[rd])])
                                    fRD.close()
                                    RDNullfigure1='.'.join(RDNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    RDNullfigure2='.'.join(RDNullTemp.split('.')[:-1])+'.NegativeBinomial'
                                    os.system('''Rscript %s %s %s %s %s %s'''%(RFigureDRSplit2,RDNullTemp,RDNullfigure1,BoxPlotColor,lineColor,RDNullfigure2))
                                    ILNullTemp=NullPath+'ILNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fIL=open(ILNullTemp,'w')
                                    for il in InsertLength.keys():
                                            print >> fIL, ' '.join([str(il),str(InsertLength[il])])
                                    fIL.close()
                                    ILNullfigure1='.'.join(ILNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    ILNullfigure2='.'.join(ILNullTemp.split('.')[:-1])+'.Bimodal'
                                    os.system('''Rscript %s %s %s %s %s %s'''%(RFigureDRSplit2,ILNullTemp,ILNullfigure1,BoxPlotColor,lineColor,ILNullfigure2))
                                    TBNullTemp=NullPath+'TBNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fTB=open(TBNullTemp,'w')
                                    for tb in TBNullDensity.keys():
                                        print >> fTB, ' '.join([str(tb),str(TBNullDensity[tb])])
                                    fTB.close()
                                    TBNullfigure1='.'.join(TBNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    TBNullfigure2='.'.join(TBNullTemp.split('.')[:-1])+'.Bimodal'
                                    os.system('''Rscript %s %s %s %s %s %s'''%(RFigureDRSplit2,TBNullTemp,TBNullfigure1,BoxPlotColor,lineColor,TBNullfigure2))
                                    os.system('''rm  %s'''%(InsertLenNullTemp))
                                    os.system('''rm  %s'''%(DRNullTemp))
                                    os.system('''rm  %s'''%(SplitNullTemp))
                                    os.system('''rm  %s'''%(ILNullTemp))
                                    os.system('''rm  %s'''%(TBNullTemp))
                                    os.system('''rm  %s'''%(RDNullTemp))
                                    if not os.path.isdir(NullPath+'IL_Null'):
                                            os.system('''mkdir %s'''%(NullPath+'IL_Null'))
                                    os.system('''mv %s %s'''%(InsertLenNullfigure1,NullPath+'IL_Null'))
                                    os.system('''mv %s %s'''%(InsertLenNullfigure2,NullPath+'IL_Null'))
                                    os.system('''mv %s %s'''%(ILNullfigure1,NullPath+'IL_Null'))
                                    os.system('''mv %s %s'''%(ILNullfigure2,NullPath+'IL_Null'))
                                    if not os.path.isdir(NullPath+'RD_Null'):
                                            os.system('''mkdir %s'''%(NullPath+'RD_Null'))
                                    os.system('''mv %s %s'''%(RDNullfigure1,NullPath+'RD_Null'))
                                    os.system('''mv %s %s'''%(RDNullfigure2,NullPath+'RD_Null'))
                                    if not os.path.isdir(NullPath+'TB_Null'):
                                            os.system('''mkdir %s'''%(NullPath+'TB_Null'))
                                    os.system('''mv %s %s'''%(TBNullfigure1,NullPath+'TB_Null'))
                                    os.system('''mv %s %s'''%(TBNullfigure2,NullPath+'TB_Null'))
                                    if not os.path.isdir(NullPath+'DR_Null'):
                                            os.system('''mkdir %s'''%(NullPath+'DR_Null'))
                                    os.system('''mv %s %s'''%(DRNullfigure1,NullPath+'DR_Null'))
                                    os.system('''mv %s %s'''%(DRNullfigure2,NullPath+'DR_Null'))
                                    if not os.path.isdir(NullPath+'Split_Null'):
                                            os.system('''mkdir %s'''%(NullPath+'Split_Null'))
                                    os.system('''mv %s %s'''%(SplitNullfigure1,NullPath+'Split_Null'))
                                    os.system('''mv %s %s'''%(SplitNullfigure2,NullPath+'Split_Null'))
                                    if not os.path.isdir(NullPath+'All_Stats'):
                                            os.system('''mkdir %s'''%(NullPath+'All_Stats'))
                                    os.system('''mv %s %s'''%(outputfileIL,NullPath+'All_Stats'))
                                    os.system('''mv %s %s'''%(outputfileStat,NullPath+'All_Stats'))
                                    if KeepFile=='no' or KeepFile=='N' or KeepFile=='No' or KeepFile=='n':
                                            os.system(r'''rm %s'''%(outputfile))
                                    else:
                                            os.system('''mv %s %s'''%(outputfile,NullPath+'All_Stats'))
                                    if KeepFigure in ['no','N','No','n']:
                                        for root, dirnames, filenames in os.walk(NullPath):
                                            for f1 in dirnames:
                                                for f2 in os.listdir(NullPath+f1):
                                                    if f2.split('.')[-1] in ['jpg','pdf','png']:
                                                        os.system(r''' rm %s'''%(NullPath+f1+'/'+f2))
                                                        #print f2
                                    CN2_File=cn2_file
                                    Ref_Seq_File=dict_opts['--ref']
                                    #DGAP=bamF.split('/')[-1].split('.')[0]
                                    Window_Size=100
                                    Mini_CN2_Region=int(cn2_length)
                                    Length_Limit=int(cn2_length)
                                    CN2_Region={} #key of hash CN2_Region is the name of each chromosome
                                    for chrom in chromos:
                                            CN2_Region[chrom]={} #key of CN2_Region[chrom] is GC_content
                                            for con in range(101):
                                                    CN2_Region[chrom][con]=[]
                                    fcn2=open(CN2_File)
                                    temp_Name='temp.Null1.'+bamF.split('/')[-1]
                                    Ref_Flag=0
                                    ftest=open(Ref_Seq_File)
                                    ptest=ftest.readline().strip().split()
                                    if len(ptest[0])>3:
                                        Ref_Flag=1
                                    while True:
                                            pcn2=fcn2.readline().strip().split()
                                            if not len(pcn2)==3: break
                                            Chromosome=pcn2[0]
                                            if Chromosome in CN2_Region.keys():
                                                #if Bam_Flag==0 and 'chr' in Chromosome:
                                                #    Chromosome=Chromosome[3:]
                                                #elif Bam_Flag==1 and not 'chr' in Chromosome:
                                                #    Chromosome='chr'+Chromosome
                                                if int(pcn2[2])-int(pcn2[1])<Length_Limit: continue
                                                if not int(pcn2[2])-int(pcn2[1])<Length_Limit:
                                                    fasta_file=NullPath+temp_Name+'.fa'
                                                    os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(Ref_Seq_File,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),fasta_file))                    
                                                    Seq1=Fasta_To_Sequence(fasta_file)
                                                    if Seq1=='ERROR!':continue
                                                    if not Seq1=='ERROR!':
                                                        #print pcn2
                                                        sam_file=NullPath+temp_Name+'.sam'
                                                        os.system(r'''samtools view %s %s:%d-%d > %s'''%(bamF,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),sam_file))
                                                        Number_Of_Windows=len(Seq1)/100
                                                        GC_Content={}
                                                        for i in range(len(Seq1)/100+1)[1:]:				
                                                                Seq2=Seq1[(i-1)*100:i*100]
                                                                GC_Content[i]=GC_Content_Calculate(Seq2)
                                                        coverage=Region_Coverage_Calculate(sam_file,Number_Of_Windows,pcn2)
                                                        for j in GC_Content.keys():
                                                            if j in coverage.keys():
                                                                CN2_Region[Chromosome][GC_Content[j][0]].append(coverage[j][-1])
                                    os.system(r'''rm %s'''%(NullPath+temp_Name+'.fa'))
                                    os.system(r'''rm %s'''%(NullPath+temp_Name+'.sam'))
                                    Output_File=NullPath+'RD_Stat/'+bamF.split('/')[-1].replace('.bam','')+'.'+genome_name+'_MP'+str(QCAlign)+'_GC_Coverage_ReadLength'
                                    Output_Path=NullPath+'RD_Stat/'
                                    if not os.path.isdir(Output_Path):
                                        os.system(r'''mkdir %s'''%(Output_Path))
                                    fo=open(Output_File,'w')
                                    print >>fo, ' '.join(chromos)
                                    print >>fo, ' '.join([str(i) for i in range(101)])
                                    for key_1 in chromos:
                                            for key_2 in range(101):
                                                    print >>fo, ':'.join(['@',','.join(str(j) for j in CN2_Region[key_1][key_2])])
                                    fo.close()
                            else:
                                print 'no info on selected chromosomes'
                        #temp_rec=workdir+dict_opts['-s'].split('/')[-1].replace('.bam','.temp')
                        #fo=open(temp_rec,'w')
                        #print >>fo,' '.join(['SVelter','1','done'])
                        #fo.close() 



