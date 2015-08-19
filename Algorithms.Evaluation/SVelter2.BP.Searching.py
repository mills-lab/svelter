#!/usr/bin/env python

#!Python
#this code is used to build the Null model of insert length distribution (based on bimodal distribution), read depth distribution (based on NB dist, 100bp/bin, corrected by GC content), number of clipped reads, and number of read pairs with abnormai directions(FF). reads to build the Null morday would come from cn2 regions larger than a customized size(1kb as default.)
#Usage:
#For debug use only
#SVelter2.BP.Searching.py --ppre workdir --ref ref.fa -s sample.bam --chr chromosome/genome
#--ILModel: Bimodal/Normal
#--RDModel: NegativeBinomial/Normal
#--TBModel Bimodal/Normal
#command='python /nfs/remills-data/xuefzhao/SVelter/SVelter2.BP.Searching.py --ppre /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.null/ -s /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.null/BamFiles/Simulate.null.RD.10.RL.101.sorted.bam --chr 2 --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/reference.flux/human_g1k_v37.fasta --NullModel S'
#sys.argv=command.split()[1:]
#Suggestion: do not specify '-o' here. output files would be under: --ppre +'Breakpoints/'
import os
import sys
import getopt
import numpy
import re
import random
import time
import datetime
ILCutoff=0.95
RDCutoff=0.95
TBCutoff=0.9999
SplitCutoff=0.999
ABILCutoff=0.99
DRCutoff=0.99
SPLCutoff=0.85
#ToolMappingQ='/mnt/EXT/Mills-data/xuefzhao/software/bigWigSummary'
#FileMappingQ='/mnt/EXT/Mills-data/xuefzhao/projects/wgEncodeCrgMapabilityAlign36mer.bigWig'
opts,args=getopt.getopt(sys.argv[1:],'o:h:s:',['NullModel=','help=','ppre=','PathBam=','sample=','chr=','SamplePbs=', 'ref=','NullGenomeName=','NullSplitLength=','NullILCI=','NullRDCI=','NullTBCI=','NullILCff=','NullSPCff=','NullDRCff=','NullBed=','NullBedLength=','NullSampleLength=','NullSampleNumber=','ExcludeBed=','IncludeBed=','ToolMappingQ=','FileMappingQ=','NullSamplePercentage=','SplitLength=','BPSPCff=','BPLNCff=','BPAlignQC=','BPAlignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)
def IL_Prob_Bimodal(IL_Length,ILStats):
	log_y1=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['bimodal']['std1']))-(float(IL_Length)-float(ILStats['bimodal']['Mean1']))**2/2/float(ILStats['bimodal']['std1'])**2
	log_y2=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['bimodal']['std2']))-(float(IL_Length)-float(ILStats['bimodal']['Mean2']))**2/2/float(ILStats['bimodal']['std2'])**2
	return log_y1*float(ILStats['bimodal']['Bimodal1'])+log_y2*float(ILStats['bimodal']['Bimodal2'])

def IL_Prob_Normal(IL_Length,ILStats):
	log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['normal']['std']))-(float(IL_Length)-float(ILStats['normal']['Mean']))**2/2/float(ILStats['normal']['std'])**2
	return log_y

def log_factorial(K):
	k=int(K)
	k_list=[0]
	for i in range(k+1)[2:]:
		k_list.append(k_list[i-2]+numpy.log(i))
	return k_list[-1]

def RD_Prob_NegativeBinomial(ReadDepth,RDStats):
	P=1-float(RDStats['Mean'])/(float(RDStats['STD'])**2)
	R=float(RDStats['Mean'])*(1-P)/P
	log_y=log_factorial(ReadDepth+R-1)-log_factorial(ReadDepth)-log_factorial(R-1)+R*numpy.log(1-P)+ReadDepth*numpy.log(P)
	return log_y

def RD_Prob_Normal(ReadDepth,RDStats):
	log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(RDStats['STD']))-(float(ReadDepth)-float(RDStats['Mean']))**2/2/float(RDStats['STD'])**2
	return log_y

def TB_Prob_Bimodal(TB_Length,TBStats):
	log_y1=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['bimodal']['std1']))-(float(TB_Length)-float(TBStats['bimodal']['Mean1']))**2/2/float(TBStats['bimodal']['std1'])**2
	log_y2=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['bimodal']['std2']))-(float(TB_Length)-float(TBStats['bimodal']['Mean2']))**2/2/float(TBStats['bimodal']['std2'])**2
	return log_y1*float(TBStats['bimodal']['Bimodal1'])+log_y2*float(TBStats['bimodal']['Bimodal2'])

def TB_Prob_Normal(TB_Length,TBStats):
	log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['normal']['std']))-(float(TB_Length)-float(TBStats['normal']['Mean']))**2/2/float(TBStats['normal']['std'])**2
	return log_y

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

def cigar2split(cigar):
    MapLen=[]
    pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
    cigars=[]
    for m in pcigar.finditer(cigar):
        cigars.append((m.groups()[0],m.groups()[1]))
    if cigars[0][1]=='S' or cigars[0][1]=='H':
        MapLen.append(0)
        cigars=cigars[1:]
    maplen=0
    for n in cigars:
        if n[1]=='M' or n[1]=='D' or n[1]=='N':
            maplen+=int(n[0])
        if n[1]=='S' or n[1]=='H': 
            MapLen.append(maplen)
    return MapLen

def cigar2splitlen(cigar):
    MapLen=[]
    pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
    cigars=[]
    for m in pcigar.finditer(cigar):
        cigars.append((m.groups()[0],m.groups()[1]))
    for i in cigars:
        if i[1]=='S' or i[1]=='H':
            MapLen.append(int(i[0]))
    return MapLen

def cigar2splitqual(cigar,qual):
    import re
    pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
    cigars=[]
    qual_out=[]
    for m in pcigar.finditer(cigar):
        cigars.append((m.groups()[0],m.groups()[1]))
    pos1=[0]
    for n in cigars:
        if n[1]=='S' or n[1]=='H':
            pos1.append(pos1[-1]+int(n[0]))
            pos1.append(pos1[-1])
        elif n[1]=='M' or n[1]=='D' or n[1]=='N':
            pos1[-1]+=int(n[0])
    pos1.remove(pos1[-1])
    qual2=[]
    for i in range(len(pos1)/2):
        qual3=qual[pos1[2*i]:pos1[2*i+1]]
        for j in qual3:
            qual2.append(ord(j))
    if not qual2==[]:
        qual_out.append(numpy.mean(qual2))
    else:
        qual_out.append(30)
    return numpy.max(qual_out)

def clusterNums(data, ClusterLen, direction):
    if data==[]:
        return [[],[]]
    else:
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

def clusterNum2(data1,data2,direction,numCff):
    data1.sort()
    out={}
    datatemp=[]
    if direction=='f':
        if data2[1][0]>numCff:
            datatemp.append([0,data2[0][0]])
        for i in range(len(data2[1]))[1:]:
            if data2[1][i]>numCff:
                datatemp.append([data2[0][i-1],data2[0][i]])
    elif direction=='r':
        for i in range(len(data2[1])-1):
            if data2[1][i]>numCff:
                datatemp.append([data2[0][i],data2[0][i+1]])
        if data2[1][-1]>numCff:
            datatemp.append([data2[0][-1],data2[0][-1]+1000])
    if len(datatemp)>0:
        flag1=0
        flag2=0
        if direction=='f':
            while True:
                if not data1[flag2] <datatemp[flag1][0] and not data1[flag2] >datatemp[flag1][1]:
                    if not datatemp[flag1][1] in out.keys():
                        out[datatemp[flag1][1]]=[data1[flag2]]
                    else:
                        out[datatemp[flag1][1]]+=[data1[flag2]]
                    flag2+=1
                elif data1[flag2] >datatemp[flag1][1]:
                    flag1+=1
                elif data1[flag2]<datatemp[flag1][0]:
                    flag2+=1
                if flag1==len(datatemp) or flag2==len(data1): break
        if direction=='r':
            while True:
                if not data1[flag2] <datatemp[flag1][0] and  not data1[flag2] >datatemp[flag1][1]:
                    if not datatemp[flag1][0] in out.keys():
                        out[datatemp[flag1][0]]=[data1[flag2]]
                    else:
                        out[datatemp[flag1][0]]+=[data1[flag2]]
                    flag2+=1
                elif data1[flag2] >datatemp[flag1][1]:
                    flag1+=1
                elif data1[flag2]<datatemp[flag1][0]:
                    flag2+=1
                if flag1==len(datatemp) or flag2==len(data1): break
        return out
    else:
        return {}

def clusterNum3(hash1,hash2):
    hash3={}
    for key in hash2.keys():
        hash3[key]=[0,0,0]
        for key2 in hash2[key]:
            hash3[key][0]+=hash1[key2][0]
            hash3[key][1]+=hash1[key2][1]
            hash3[key][2]+=hash1[key2][2]
    return hash3

def clusterNums4(hash, ClusterLen,direction):
    data1=clusterNums(hash.keys(),ClusterLen,direction)[0]
    if direction=='r':
        data1=data1[::-1]
    temp=[[] for i in data1]
    if direction=='f':
        rec_temp=0
        for k1 in sorted(hash.keys()):
            if k1<data1[rec_temp]:
                temp[rec_temp].append(k1)
            elif k1==data1[rec_temp]:
                temp[rec_temp].append(k1)
                rec_temp+=1
    elif direction=='r':
        rec_temp=0
        for k1 in sorted(hash.keys())[::-1]:
            if k1>data1[rec_temp]:
                temp[rec_temp].append(k1)
            elif k1==data1[rec_temp]:
                temp[rec_temp].append(k1)
                rec_temp+=1
    data2=[0 for i in data1]
    rec_temp=-1
    for k1 in temp:
        rec_temp+=1
        for k2 in k1:
            data2[rec_temp]+=hash[k2]
    return [data1,data2]

def clusterNumSP(hash,ClusterLen,numCff):
    #eg: cluSP= clusterNumSP(SP,SPCluMin,SplitMin)
    data=hash.keys()
    data.sort()
    out1=clusterNums(data, ClusterLen, 'f')
    out2=clusterNums(data, ClusterLen, 'r')
    out={}
    flag1=0
    flag2=0
    while True:
        if data[flag2]< out2[0][flag1]:
            flag2+=1
        elif not data[flag2]< out2[0][flag1] and not data[flag2]> out1[0][flag1]:
            if not out1[0][flag1] in out.keys():
                out[out1[0][flag1]]=[out2[0][flag1],hash[data[flag2]]]
            else:
                out[out1[0][flag1]][1]+=hash[data[flag2]]
            flag2+=1
        elif data[flag2]> out1[0][flag1]:
            flag1+=1
        if flag1==len(out1[0]) or flag2==len(hash): break
    for i in out.keys():
        if not out[i][1]>numCff:
            del out[i]
    return out 

def clusterNumLN(hash,ClusterLen,numCff):
    #eg: hash= KeyF3
    out={}
    for key in hash.keys():
        out[key]={}
        key2={}
        for key3 in hash[key]:
            if not key3[1] in key2.keys():
                key2[key3[1]]=key3[0]
            else:
                key2[key3[1]]+=key3[0]
        for key4 in key2.keys():
            if len(key2[key4])>numCff:
                if key4=='ff' or key4=='rf':
                    temp1=clusterNums(key2[key4],ClusterLen,'f')
                    for temp2 in range(len(temp1[1])):
                        if temp1[1][temp2]>numCff:
                            if not key4 in out[key].keys():
                                out[key][key4]=temp1[0][temp2]
                elif key4=='rr' or key4=='fr':
                    temp1=clusterNums(key2[key4],ClusterLen,'r')
                    for temp2 in range(len(temp1[1])):
                        if temp1[1][temp2]>numCff:
                            if not key4 in out[key].keys():
                                out[key][key4]=temp1[0][temp2]
                elif key4==0:
                    if len(key2[key4])>numCff:
                        out[key][key4]=0
                else:
                    if key4[-1]=='+':
                        temp=clusterNums(key2[key4],ClusterLen,'f')
                    else:
                        temp=clusterNums(key2[key4],ClusterLen,'r')
                    for temp2 in range(len(temp[1])):
                        if temp[1][temp2]>numCff:
                            if not key4 in out[key].keys():
                                if key4[-1]=='+':
                                    out[key][key4]=[temp[0],'f']
                                elif key4[-1]=='-':
                                    out[key][key4]=[temp[0],'r']
    return out

def clusterQC(hash, CluNum):
    out=[]
    for i in range(len(hash[1])):
        if not hash[1][i]<CluNum:
            out.append(hash[0][i])
    return out

def clusterSupVis2(dataRound, dataLeft, dataRight,keypick):
    flag1=0
    flag2=0
    out={}
    if keypick=='left':
        while True:
            if flag1==len(dataRound): break
            else:
                if flag2==len(dataLeft)-1 and dataRound[flag1]>dataRight[flag2]: break
                else:
                    if dataRound[flag1]<dataLeft[flag2]:
                        flag1+=1
                    elif dataRound[flag1]>dataRight[flag2]:
                        flag2+=1
                    else:
                        if not dataLeft[flag2] in out.keys():
                            out[dataLeft[flag2]]=[dataRound[flag1]]
                        else:
                            out[dataLeft[flag2]]+=[dataRound[flag1]]
                        flag1+=1
    if keypick=='right':
        while True:
            if flag1==len(dataRound): break
            else:
                if flag2==len(dataLeft)-1 and dataRound[flag1]>dataRight[flag2]: break
                else:
                    if dataRound[flag1]<dataLeft[flag2]:
                        flag1+=1
                    elif dataRound[flag1]>dataRight[flag2]:
                        flag2+=1
                    else:
                        if not dataRight[flag2] in out.keys():
                            out[dataRight[flag2]]=[dataRound[flag1]]
                        else:
                            out[dataRight[flag2]]+=[dataRound[flag1]]
                        flag1+=1
    return out

def clusterSupVis(dataRou,dataCen,ClusterLen):
    flag1=0
    flag2=0
    dataRou.sort()
    dataCen.sort()
    out={}
    while True:
        if flag1==len(dataRou): break
        else:
            if flag2==len(dataCen)-1:
                if dataRou[flag1]>dataCen[flag2]+ClusterLen: break
                elif dataRou[flag1]<dataCen[flag2]-ClusterLen: 
                    flag1+=1
                else:
                    if not dataCen[flag2] in out.keys():
                        out[dataCen[flag2]]=[dataRou[flag1]]
                    else:
                        out[dataCen[flag2]]+=[dataRou[flag1]]
                    flag1+=1
            else:
                if dataRou[flag1]<dataCen[flag2]-ClusterLen:
                    flag1+=1
                elif dataRou[flag1]>dataCen[flag2]+ClusterLen:
                    flag2+=1
                else:
                    if not dataCen[flag2] in out.keys():
                        out[dataCen[flag2]]=[dataRou[flag1]]
                    else:
                        out[dataCen[flag2]]+=[dataRou[flag1]]
                    flag1+=1
    return out

def clusterLN(data, LnHash, numCff):
    out={}
    for i in data:
        if not i in LnHash.keys():
            continue
        else:
            j=LnHash[i]
            key=j[0]+j[2]
            if not key in out.keys():
                out[key]=[j[1]]
            else:
                out[key]+=[j[1]]
    for key in out.keys():
        if len(out[key])<numCff:
            del out[key]
    return out

def clusterSupVis3(dataRound,dataCen):
    #dataRound and dataCen has to be sorted
    out={}
    for key1 in dataCen:
        min_num=abs(dataRound[0]-key1)
        min_rec=dataRound[0]
        for key2 in dataRound[1:]:
            alt_num=abs(key2-key1)
            alt_rec=key2
            if alt_num<min_num:
                min_num=alt_num
                min_rec=alt_rec
        out[key1]=min_rec
    return out

def BPFilter_ReadNumThroughBP(BamInput,chrom,bpCoordinate,flank):
    #direction could be 'f','r' and 'bi'
    flank=int(flank)
    bps=[bpCoordinate-flank,bpCoordinate+flank]
    if Bam_Flag==1:
        fin=os.popen(r'''samtools view %s chr%s:%d-%d'''%(BamInput,chrom,int(bps[0]),int(bps[1])))
    else:
        fin=os.popen(r'''samtools view %s %s:%d-%d'''%(BamInput,chrom,int(bps[0]),int(bps[1])))
    RDs=[0*i for i in range(bps[1]-bps[0])]
    while True:
        pbam=fin.readline().strip().split()
        if not pbam:break
        if int(pbam[4])<int(QCAlign): continue 
        if int(pbam[1])&4>0: continue
        ReadLen=cigar2reaadlength(pbam[5])
        if ReadLen<90:
            pos_pbam=[int(pbam[3]),int(pbam[3])+ReadLen]
        else:
            diff_pbam=(ReadLen-90)/2
            pos_pbam=[int(pbam[3])+diff_pbam,int(pbam[3])+ReadLen-diff_pbam]
        for j in range(pos_pbam[0]-bps[0],pos_pbam[1]-bps[0]+1):
            if j in range(len(RDs)):
                RDs[j]+=1
    fin.close()
    min_through=float('inf')
    for i in range(len(RDs)):
        if RDs[i]<min_through:
            min_through=RDs[i]
            min_record=i
    if min_through<RD_Mean/2:
        return min_record+int(bps[0])
    else:
        return -1

def clusterSupVis4(dataRou,left,right):
	out=[]
	for i in sorted(dataRou):
		if i<left:
			continue
		if i>right:
			break
		else:
			out.append(i)
	return out

def overlap_calcu(list1,list2):
    list1.sort()
    list2.sort()
    list_len=min(len(list1),len(list2))
    rec=0
    reca=0
    recb=0
    while True:
        if reca==len(list1) or recb==len(list2): break
        if list1[reca]==list2[recb]:
            rec+=1
            reca+=1
            recb+=1
        elif list1[reca]<list2[recb]:
            reca+=1
        elif list1[reca]>list2[recb]:
            recb+=1
    return float(rec)/float(list_len)

def lower_ILCff_calcu(IL_Null_Stat,CICff):
    fin=open(IL_Null_Stat)
    flag=0
    outhash={}
    totalnum=0
    for line in fin:
        pin=line.strip().split()
        if pin[0]=='ReadDepth' and  pin[1]=='Frequency':
            flag=0
        if flag==1:
            outhash[int(pin[0])]=int(pin[1])
            totalnum+=int(pin[1])
        if pin[0]=='InsertLength' and pin[1]=='Frequency':
            flag=1
    fin.close()
    subnum=0
    for k1 in sorted(outhash.keys()):
        subnum+=outhash[k1]
        if float(subnum)/float(totalnum)>CICff/2:
            break
    return k1

def bp_subgroup(list, Min_Distinguish_Len):
    if list==[]:
        return []
    else:
        out=[[list[0]]]
        list.sort()
        for x in range(len(list)-1):
            if list[x+1]-list[x]>Min_Distinguish_Len:
                out.append([list[x+1]])
            else:
                out[-1].append(list[x+1])
        return out

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

def file_start(filename):
    fRDind=open(filename,'w')
    fRDind.close()

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

def path_mkdir(path):
        if not os.path.isdir(path):
                os.system(r'''mkdir %s'''%(path))

Length_Limit=2000
CN2_Region={}

#ToolMappingQ='/mnt/EXT/Mills-data/xuefzhao/software/bigWigSummary'
if not '--NullModel' in dict_opts.keys():
    model_comp='S'
else:
    if dict_opts['--NullModel'] in ['S','Simple']:
        model_comp='S'
    else:
        model_comp='C'

if '--ToolMappingQ' in dict_opts.keys() and '--FileMappingQ' in dict_opts.keys():
    ToolMappingQ=dict_opts['--ToolMappingQ']
    FileMappingQ=dict_opts['--FileMappingQ']
    align_QCflag=1
else:
    align_QCflag=0

if '--BPAlignQC' in dict_opts.keys():
    BPAlignQC=float(dict_opts['--BPAlignQC'])
else:
    BPAlignQC=0.0

if '--QCAlign' in dict_opts.keys():
    QCAlign=int(dict_opts['--QCAlign'])
else:
    QCAlign=20

if '--QCSplit' in dict_opts.keys():
    QCSplit=int(dict_opts['--QCSplit'])
else:
    QCSplit=20

if '--NullSplitLength' in dict_opts.keys():
    NullSplitLen_perc=float(dict_opts['--NullSplitLength'])
else:
    NullSplitLen_perc=0.9

if '--BPAlignQCFlank' in dict_opts.keys():
    BPAlignQCFlank=int(dict_opts['--BPAlignQCFlank'])
else:
    BPAlignQCFlank=500

if NullSplitLen_perc>1:
    NullSplitLen_perc=float(NullSplitLen_perc)/float(ReadLength)

if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
    print 'SVelter-0.1          Last Update:2014-08-20'
    print 'Required Parameters:'
    print '--ppre, workind directory of SVelter, eg: .../SVelter/' 
    print '--ref, absolute path of reference genome. eg: .../SVelter/reference/genome.fa'
    print '-s, absolute path of input bam file. eg: .../SVelter/BamFiles/sample.bam'
    #print '--PathBam, absolute path of folder containing all input bam files. eg: .../SVelter/BamFiles/'
    print 'Optional Parameters:'
    print '--NullBed, a describing bed files containing genomic regions that are used to build Null model; usually regions known to have few structural variants'
    print '--NullBedLength, minimum requirement for length of regions used to build Null model; default: 2kb'
    print '--NullSamplePercentage, sampling percentage of regions from Nullbed'
    print '--NullSampleLength, if not --NullBed provided, SVelter will randomly sample regions on genome to build Null model; --NullSampleLength specify the length of region picked; default: 5kb'
    print '--NullSampleNumber, if not --NullBed provided, SVelter will randomly sample regions on genome to build Null model; --NullSampleNumber specify the number of region picked; default: 10k'
    print ''
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
            if '-s' in dict_opts.keys():
                bam_path='/'.join(dict_opts['-s'].split('/')[:-1])+'/'
                bam_files=[dict_opts['-s']]
            else:
                bam_path=dict_opts['--PathBam']
                if not bam_path[-1]=='/':
                    bam_path+='/'
                bam_files=[]
                for file in os.listdir(bam_path):
                    if file.split('.')[-1]=='bam':
                        bam_files.append(bam_path+file)
            if not '--ref' in dict_opts.keys():
                print 'Error: please specify refrence genome using --ref'
            else:
                ref_path='/'.join(dict_opts['--ref'].split('/')[:-1])+'/'
                ref_file=dict_opts['--ref']
                ref_index=dict_opts['--ref']+'.fai'
                if not os.path.isfile(ref_index):
                    print 'Error: reference genome not indexed'
                else:
                    chromos=[]
                    fin=open(ref_index)
                    for line in fin:
                        pin=line.strip().split()
                        chromos.append(pin[0])
                    fin.close()
                    ref_flag=0
                    if 'chr' in chromos[0]:
                        ref_flag=1
                    if '--chr' in dict_opts.keys():
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
                        if '--NullGenomeName' in dict_opts.keys():
                            genome_name=dict_opts['--NullGenomeName']
                        else:
                            genome_name='genome'
                        if not '-o' in dict_opts.keys():
                            out_path=workdir
                        else:
                            out_path=dict_opts['-o']
                            if not out_path[-1]=='/':
                                out_path+='/'
                        BPPath=out_path+'BreakPoints.'+dict_opts['-s'].split('/')[-1]+'/'
                        if not os.path.isdir(BPPath):
                            os.system(r'''mkdir %s'''%(BPPath))
                        NullPath=out_path+'NullModel.'+dict_opts['-s'].split('/')[-1]+'/'
                        for bamF in bam_files:
                            fin=open(NullPath+'All_Stats/'+bamF.split('/')[-1].replace('.bam','.')+genome_name+'.Stats')
                            for line in fin:
                                pin=line.strip().split()
                            ReadLength=int(pin[-1].split(':')[1])
                            fin.close()
                            for chrF in chromos:
                                bamF_Name=bamF.split('/')[-1].replace('.bam','')
                                floc_Name=BPPath+bamF_Name+'.'+chrF
                                Refloc_name='.'.join(ref_file.split('.')[:-1])+'.Mappable.'+chrF+'.bed'
                                if os.path.isfile(Refloc_name):
                                    BamInput=bamF
                                    ILStat=NullPath+'IL_Null/ILNull.'+bamF_Name+'.'+genome_name+'.Bimodal'
                                    IL_Null_Stat=NullPath+'All_Stats/'+bamF_Name+'.'+genome_name+'.density.null'
                                    RDStat=NullPath+'RD_Null/RDNull.'+bamF_Name+'.'+genome_name+'.NegativeBinomial'
                                    TBStat=NullPath+'TB_Null/TBNull.'+bamF_Name+'.'+genome_name+'.Bimodal'
                                    SPLenStat=NullPath+'SplitLength/'+bamF_Name+'.'+genome_name+'.SplitLength'
                                    AllStat=NullPath+'All_Stats/'+bamF_Name+'.'+genome_name+'.Stats'
                                    ILStats={}
                                    fILS=open(ILStat)
                                    if model_comp=='C':
                                        pILS1=fILS.readline().strip().split()
                                        pILS2=fILS.readline().strip().split()
                                        ILStats['stat']={}
                                        for i in range(len(pILS1)):
                                                ILStats['stat'][pILS1[i]]=pILS2[i]
                                        pILS1=fILS.readline().strip().split()
                                        pILS2=fILS.readline().strip().split()
                                        ILStats['bimodal']={}
                                        for i in range(len(pILS1)):
                                                ILStats['bimodal'][pILS1[i]]=pILS2[i]
                                        pILS1=fILS.readline().strip().split()
                                        pILS2=fILS.readline().strip().split()
                                        for i in range(len(pILS1)):
                                                ILStats['bimodal'][pILS1[i]]=pILS2[i]
                                        pILS1=fILS.readline().strip().split()
                                        pILS2=fILS.readline().strip().split()
                                        ILStats['normal']={}
                                        for i in range(len(pILS1)):
                                                ILStats['normal'][pILS1[i]]=pILS2[i]
                                    elif model_comp=='S':
                                        pILS1=fILS.readline().strip().split()
                                        pILS2=fILS.readline().strip().split()
                                        ILStats['stat']={}
                                        for i in range(len(pILS1)):
                                                ILStats['stat'][pILS1[i]]=pILS2[i]
                                        pILS1=fILS.readline().strip().split()
                                        pILS2=fILS.readline().strip().split()
                                        ILStats['normal']={}
                                        for i in range(len(pILS1)):
                                                ILStats['normal'][pILS1[i]]=pILS2[i]
                                        ILStats['bimodal']={}
                                        for i in range(len(pILS1)):
                                                ILStats['bimodal'][pILS1[i]]=pILS2[i]
                                        for i in range(len(pILS1)):
                                                ILStats['bimodal'][pILS1[i]]=pILS2[i]
                                    fILS.close()    
                                    #Read In Statistics for RD distribution, including mean, std, etc. 
                                    RDStats={}
                                    fRDS=open(RDStat)
                                    pRDS1=fRDS.readline().strip().split()
                                    pRDS2=fRDS.readline().strip().split()
                                    for i in range(len(pRDS1)):
                                            RDStats[pRDS1[i]]=pRDS2[i]
                                    fRDS.close()
                                    #Read In Statistics for TBdistribution, including mean, std, etc. 
                                    TBStats={}
                                    fTBS=open(TBStat)
                                    if model_comp=='C':
                                        pTBS1=fTBS.readline().strip().split()
                                        pTBS2=fTBS.readline().strip().split()
                                        TBStats['stat']={}
                                        for i in range(len(pTBS1)):
                                                TBStats['stat'][pTBS1[i]]=pTBS2[i]
                                        pTBS1=fTBS.readline().strip().split()
                                        pTBS2=fTBS.readline().strip().split()
                                        TBStats['bimodal']={}
                                        for i in range(len(pTBS1)):
                                                TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
                                        pTBS1=fTBS.readline().strip().split()
                                        pTBS2=fTBS.readline().strip().split()
                                        for i in range(len(pTBS1)):
                                                TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
                                        pTBS1=fTBS.readline().strip().split()
                                        pTBS2=fTBS.readline().strip().split()
                                        TBStats['normal']={}
                                        for i in range(len(pTBS1)):
                                                TBStats['normal'][pTBS1[i]]=pTBS2[i]
                                    elif model_comp=='S':
                                        pTBS1=fTBS.readline().strip().split()
                                        pTBS2=fTBS.readline().strip().split()
                                        TBStats['stat']={}
                                        for i in range(len(pTBS1)):
                                                TBStats['stat'][pTBS1[i]]=pTBS2[i]
                                        pTBS1=fTBS.readline().strip().split()
                                        pTBS2=fTBS.readline().strip().split()
                                        TBStats['normal']={}
                                        for i in range(len(pTBS1)):
                                                TBStats['normal'][pTBS1[i]]=pTBS2[i]
                                        TBStats['bimodal']={}
                                        for i in range(len(pTBS1)):
                                                TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
                                        for i in range(len(pTBS1)):
                                                TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
                                    fTBS.close()
                                    #Read In Statistics for Clipped Length, finding the minimum length considered for a real split read
                                    fSPLStat=open(SPLenStat)
                                    fSPLStat.readline()
                                    SPLength={}
                                    totalSPlit=0
                                    subSPlit=0 
                                    if NullSplitLen_perc>1:
                                        SPLCff=NullSplitLen_perc
                                    else:        
                                        while True:
                                                pSPLStat=fSPLStat.readline().strip().split()
                                                if not pSPLStat: break
                                                SPLength[int(pSPLStat[0])]=int(pSPLStat[1])
                                                totalSPlit+=int(pSPLStat[1])
                                        if not SPLength.keys()==[]:
                                            for key in SPLength.keys():
                                                    subSPlit+=SPLength[key]
                                                    if float(subSPlit)/float(totalSPlit)>NullSplitLen_perc: 
                                                            break
                                            SPLCff=key-1
                                            if SPLCff>ReadLength/10:
                                                SPLCff=ReadLength/10
                                        else:
                                            SPLCff=ReadLength/10
                                    #Read In CutOff Insormation for each parametre, like IL, RD, TB, Split, DR	
                                    fAllS=open(AllStat)
                                    pAllS=fAllS.readline().rstrip()
                                    #print pAllS+':'
                                    ILCIs={}
                                    pAllS1=fAllS.readline().strip().split()
                                    pAllS2=fAllS.readline().strip().split()
                                    for i in range(len(pAllS1)):
                                            ILCIs[pAllS1[i]]=pAllS2[i]
                                    pAllS=fAllS.readline().rstrip()
                                    #print pAllS+':'
                                    RDCIs={}
                                    pAllS1=fAllS.readline().strip().split()
                                    pAllS2=fAllS.readline().strip().split()
                                    for i in range(len(pAllS1)):
                                        RDCIs[pAllS1[i]]=pAllS2[i]
                                    pAllS=fAllS.readline().rstrip()
                                    #print pAllS+':'
                                    TBCIs={}
                                    pAllS1=fAllS.readline().strip().split()
                                    pAllS2=fAllS.readline().strip().split()
                                    for i in range(len(pAllS1)):
                                        TBCIs[pAllS1[i]]=pAllS2[i]
                                    pAllS=fAllS.readline().rstrip()
                                    #print pAllS+':'
                                    InsertLen={}
                                    pAllS1=fAllS.readline().strip().split()
                                    pAllS2=fAllS.readline().strip().split()
                                    for i in range(len(pAllS1)):
                                            InsertLen[pAllS1[i]]=pAllS2[i]
                                    for i in range(len(pAllS1)):
                                            if float(pAllS2[i])>ABILCutoff:
                                                    InsertLenMin=int(pAllS1[i])-1
                                                    break
                                    if InsertLenMin<5:
                                            InsertLenMin=5
                                    pAllS=fAllS.readline().rstrip()
                                    #print pAllS+':'
                                    SplitReads={}
                                    pAllS1=fAllS.readline().strip().split()
                                    pAllS2=fAllS.readline().strip().split()
                                    for i in range(len(pAllS1)):
                                            SplitReads[pAllS1[i]]=pAllS2[i]
                                    for i in range(len(pAllS1)):
                                            if float(pAllS2[i])>SplitCutoff:
                                                    SplitMin=int(pAllS1[i])-1
                                                    break
                                    if SplitMin<5:
                                            SplitMin=5
                                    #print 'Maximum Number of Split Reads allowed:'+str(SplitMin+1)
                                    pAllS=fAllS.readline().rstrip()
                                    #print pAllS+':'
                                    AbDirection={}
                                    pAllS1=fAllS.readline().strip().split()
                                    pAllS2=fAllS.readline().strip().split()
                                    for i in range(len(pAllS1)):
                                            AbDirection[pAllS1[i]]=pAllS2[i]
                                    for i in range(len(pAllS1)):
                                        if float(pAllS2[i])>DRCutoff:
                                            DRMin=int(pAllS1[i])-1
                                            break
                                    if DRMin<5:
                                            DRMin=5
                                    #print 'Maximum Number of Abnormal Direction Read Pairs allowed:'+str(DRMin+1)
                                    fbam=os.popen('''samtools view -H %s'''%(BamInput))
                                    if '--BPSPCff' in dict_opts.keys():
                                        BPSPCff=int(float(dict_opts['--BPSPCff']))
                                    else:
                                        BPSPCff=int(round(2.0*float(RDStats['Median'])/float(10)))
                                    if BPSPCff<3:
                                        BPSPCff=3
                                    if '--BPLNCff' in dict_opts.keys():
                                        BPLNCff=int(float(dict_opts['--BPLNCff']))
                                    else:
                                        BPLNCff=int(round(3.0*float(RDStats['Median'])/float(10)))
                                    if BPLNCff<3:
                                        BPLNCff=3
                                    SPCluMin=BPSPCff
                                    LnCluMin=BPLNCff
                                    if '--SPCluLen' in dict_opts.keys():
                                        SPCluLen=int(dict_opts['--SPCluLen'])
                                    else:
                                        SPCluLen=5
                                    SubLnCluMin=max([LnCluMin,SPCluMin])
                                    LinkCluMin=min([LnCluMin,SPCluMin])
                                    #QCSP=2
                                    #SubLnCluMin=min(SPCluMin,LnCluMin)
                                    ClusterLen=float(ILStats['normal']['Mean'])+2*float(ILStats['normal']['STD'])-ReadLength
                                    Min_Distinguish_Len=100
                                    subLnClusterLen=ClusterLen/2
                                    CICff=0.003
                                    ILCffs=[lower_ILCff_calcu(IL_Null_Stat,CICff),int(float(ILStats['normal']['Mean'])+3*float(ILStats['normal']['STD']))]
                                    #if float(ILStats['stat']['Mean'])/float(ILStats['stat']['STD'])<2: 
                                    #    ILCffs=[int(float(ILStats['stat']['Mean'])-float(ILStats['stat']['STD'])),int(float(ILStats['stat']['Mean'])+float(ILStats['stat']['STD']))]
                                    #elif float(ILStats['stat']['Mean'])/float(ILStats['stat']['STD'])<4: 
                                    #    ILCffs=[int(float(ILStats['stat']['Mean'])-2*float(ILStats['stat']['STD'])),int(float(ILStats['stat']['Mean'])+2*float(ILStats['stat']['STD']))]
                                    #else: 
                                    #    ILCffs=[int(float(ILStats['stat']['Mean'])-3*float(ILStats['stat']['STD'])),int(float(ILStats['stat']['Mean'])+3*float(ILStats['stat']['STD']))]
                                    BPOutputa=floc_Name+'.'+'.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.SPs'
                                    fout=open(BPOutputa,'w')
                                    fout.close()
                                    BPOutputb=floc_Name+'.'+'.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.LNs'
                                    fout=open(BPOutputb,'w')
                                    fout.close()
                                    #BPOutputc=floc_Name+'_.'+'_.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.cluFR'
                                    #fout=open(BPOutputc,'w')
                                    #fout.close()
                                    BPOutputd=floc_Name+'.'+'.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.chromLNs'
                                    fout=open(BPOutputd,'w')
                                    fout.close()
                                    BPOutpute='/'.join(BPOutputd.split('/')[:-1])+'/'+'.'.join(BPOutputd.split('/')[-1].split('.')[:-6]+BPOutputd.split('/')[-1].split('.')[-5:])
                                    if not os.path.isfile(BPOutpute):
                                        fout=open(BPOutpute,'w')
                                        fout.close()
                                    #BPOutpute=floc_Name+'_.'+'_.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.cluLN'
                                    #fout=open(BPOutpute,'w')
                                    #fout.close()
                                    #print 'Searching for break points of sample: '+bamF_Name+', chromosome:'+chrF
                                    time1=time.time()
                                    abtLink={}
                                    floc=open(Refloc_name)
                                    loc_rec={}
                                    for line in floc:
                                        ploc=line.strip().split()
                                        if not ploc[0] in loc_rec.keys():
                                            loc_rec[ploc[0]]=[]
                                            loc_rec[ploc[0]].append([int(ploc2) for ploc2 in ploc[1:]])
                                        else:
                                            loc_rec[ploc[0]].append([int(ploc2) for ploc2 in ploc[1:]])
                                    floc.close()
                                    if not loc_rec=={}:
                                        test_mul_RP=[]
                                        test_mul_SP=[]
                                        chrom=chrF
                                        mini_fout_Name=BPPath+bamF_Name+'.mini.'+loc_rec.keys()[0]+'.sam'
                                        mini_fout_N2=BPPath+bamF_Name+'.mini.'+chrom+'.bam'
                                        mini_fout_N3=BPPath+bamF_Name+'.mini.'+chrom+'.sorted'
                                        mini_fout_N4=BPPath+bamF_Name+'.mini.'+chrom+'.sorted.bam'
                                        if not os.path.isfile(mini_fout_N4):
                                            os.system(r''' samtools view -H %s -o %s'''%(BamInput,mini_fout_Name)) 
                                            RD_index_Path=NullPath+'RD_Stat/'
                                            if not os.path.isdir(RD_index_Path):
                                                os.system(r'''mkdir %s'''%(RD_index_Path))
                                            RD_index_File=RD_index_Path+bamF_Name+'.'+chrF+'.RD.index'
                                            fRDind=open(RD_index_File,'w')
                                            fRDind.close()
                                            ClusterLen2=int(ClusterLen/100+1)*100
                                            for loc in loc_rec[chrom]:         
                                                if not loc[1]-loc[0]>2*10**6:
                                                    loc2=[loc]
                                                else:
                                                    sublocNum=(loc[1]-loc[0])/10**6
                                                    loc2=[[loc[0],loc[0]+10**6+int(ClusterLen2)]]
                                                    for slNum in range(sublocNum)[1:]:
                                                        loc2.append([loc[0]+slNum*10**6-int(ClusterLen2),loc[0]+(slNum+1)*10**6+int(ClusterLen2)])
                                                    if sublocNum>1:
                                                        loc2.append([loc[0]+(slNum+1)*10**6-int(ClusterLen2),loc[1]])
                                                for real_region in loc2:
                                                        fmini=open(mini_fout_Name,'a')
                                                        fRDind=open(RD_index_File,'a')
                                                        #flink=open(mini_link_Name,'a')
                                                        #if loc_flag==1:
                                                        #    print >>fRDind,'str'+chrom+':'+str(real_region[0])+'-'+str(real_region[1])
                                                        #else:
                                                        print >>fRDind,chrom+':'+str(real_region[0])+'-'+str(real_region[1])                
                                                        RD_RealRegion=[0 for i in range((real_region[1]-real_region[0])/100+1)]
                                                        fbam=os.popen('''samtools view %s %s:%d-%d'''%(BamInput,chrom,real_region[0],real_region[1]))
                                                        while True:
                                                            pbam1=fbam.readline().strip()
                                                            if not pbam1: break
                                                            pbam=pbam1.split()
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect(pbam[1])
                                                            ReadLen=cigar2reaadlength(pbam[5])
                                                            pos1=int(pbam[3])
                                                            pos2=int(pbam[3])+ReadLen
                                                            if pos2>real_region[0] and pos1<real_region[1]:
                                                                if pos1<real_region[0] and pos2>real_region[0]:
                                                                    pos1=real_region[0]
                                                                if pos1<real_region[1] and pos2>real_region[1]:
                                                                    pos2=real_region[1]
                                                                block1=(pos1-real_region[0])/100
                                                                block2=(pos2-real_region[0])/100
                                                                if block1==block2:
                                                                    RD_RealRegion[block1]+=ReadLen
                                                                else:
                                                                    blpos1=pos1-(block1*100+real_region[0])
                                                                    blpos2=pos2-(block2*100+real_region[0])
                                                                    RD_RealRegion[block1]+=(100-blpos1)
                                                                    RD_RealRegion[block2]+=blpos2
                                                            if int(pbam[4])>int(QCAlign):             #fail quality control, skip
                                                                absIL=abs(int(pbam[8]))
                                                                QCFlag=0
                                                                link_flag=0
                                                                if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                                                                    QCFlag+=1
                                                                if DRtemp==['+','-'] and int(pbam[8])<0:
                                                                    QCFlag+=1
                                                                if DRtemp==['+','+'] or DRtemp==['-','-']:
                                                                    QCFlag+=1
                                                                if int(pbam[8])==0:
                                                                   QCFlag+=1
                                                                if not pbam[5].find('S')==-1:
                                                                    QCFlag+=1
                                                                if not QCFlag==0:
                                                                    print>>fmini, pbam1
                                                        fbam.close()
                                                        fmini.close()
                                                        for rfrr in range(len(RD_RealRegion[:-1])):
                                                            RD_RealRegion[rfrr]=str(float(RD_RealRegion[rfrr])/100.0)
                                                        if real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100 ==0:
                                                            del RD_RealRegion[-1]
                                                        else:
                                                            RD_RealRegion[-1]=str(float(RD_RealRegion[-1])/float(real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100))
                                                        print >>fRDind, ' '.join(RD_RealRegion)
                                                        fRDind.close()
                                            os.system(r'''samtools view -h -Sb %s -o %s'''%(mini_fout_Name,mini_fout_N2))
                                            os.system(r'''samtools sort %s %s'''%(mini_fout_N2,mini_fout_N3))
                                            os.system(r'''samtools index %s '''%(mini_fout_N4))
                                            os.system(r'''rm %s'''%(mini_fout_N2))
                                            os.system(r'''rm %s'''%(mini_fout_Name))
                                        if not os.path.isfile(mini_fout_N4):
                                            os.system(r''' samtools view -H %s -o %s'''%(BamInput,mini_fout_Name)) 
                                            RD_index_Path=NullPath+'RD_Stat/'
                                            if not os.path.isdir(RD_index_Path):
                                                os.system(r'''mkdir %s'''%(RD_index_Path))
                                            RD_index_File=RD_index_Path+bamF_Name+'.'+chrF+'.RD.index'
                                            GC_index_File=RD_index_Path+bamF_Name+'.'+chrF+'.GC.index'
                                            file_start(RD_index_File)
                                            file_start(GC_index_File)
                                            ClusterLen2=int(ClusterLen/100+1)*100
                                            CN2_Region[chrom]={}
                                            for con in range(101):
                                                    CN2_Region[chrom][con]=[]
                                            for loc in loc_rec[chrom]:  
                                                #print loc
                                                Chromosome=chrom
                                                if loc[1]-loc[0]<Length_Limit:continue
                                                if not loc[1]-loc[0]>2*10**6:
                                                    loc2=[loc]
                                                else:
                                                    sublocNum=(loc[1]-loc[0])/10**6
                                                    loc2=[[loc[0],loc[0]+10**6+int(ClusterLen2)]]
                                                    for slNum in range(sublocNum)[1:]:
                                                        loc2.append([loc[0]+slNum*10**6-int(ClusterLen2),loc[0]+(slNum+1)*10**6+int(ClusterLen2)])
                                                    if sublocNum>1:
                                                        loc2.append([loc[0]+(slNum+1)*10**6-int(ClusterLen2),loc[1]])
                                                for real_region in loc2:
                                                    pcn2=[chrom]+real_region
                                                    temp_Name='temp.Null1.'+bamF.split('/')[-1]+'.'+chrom
                                                    fasta_file=NullPath+temp_Name+'.fa'
                                                    os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref_file,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),fasta_file))                    
                                                    Seq1=Fasta_To_Sequence(fasta_file)
                                                    if Seq1=='ERROR!':continue
                                                    if not Seq1=='ERROR!':
                                                        #print pcn2
                                                        sam_file=NullPath+temp_Name+'.sam'
                                                        os.system(r'''samtools view %s %s:%d-%d > %s'''%(bamF,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),sam_file))
                                                        Number_Of_Windows=len(Seq1)/100
                                                        GC_Content={}
                                                        for i in range(len(Seq1)/100+1):                
                                                                Seq2=Seq1[(i)*100:(i+1)*100]
                                                                GC_Content[i]=GC_Content_Calculate(Seq2)
                                                        #print real_region        
                                                        fmini=open(mini_fout_Name,'a')
                                                        fRDind=open(RD_index_File,'a')
                                                        #flink=open(mini_link_Name,'a')
                                                        #if loc_flag==1:
                                                        #    print >>fRDind,'str'+chrom+':'+str(real_region[0])+'-'+str(real_region[1])
                                                        #else:
                                                        print >>fRDind,chrom+':'+str(real_region[0])+'-'+str(real_region[1])                
                                                        RD_RealRegion=[0 for i in range((real_region[1]-real_region[0])/100+1)]
                                                        fbam=os.popen('''samtools view %s %s:%d-%d'''%(BamInput,chrom,real_region[0],real_region[1]))
                                                        while True:
                                                            pbam1=fbam.readline().strip()
                                                            if not pbam1: break
                                                            pbam=pbam1.split()
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect(pbam[1])
                                                            ReadLen=cigar2reaadlength(pbam[5])
                                                            pos1=int(pbam[3])
                                                            pos2=int(pbam[3])+ReadLen
                                                            if pos2>real_region[0] and pos1<real_region[1]:
                                                                if pos1<real_region[0] and pos2>real_region[0]:
                                                                    pos1=real_region[0]
                                                                if pos1<real_region[1] and pos2>real_region[1]:
                                                                    pos2=real_region[1]
                                                                block1=(pos1-real_region[0])/100
                                                                block2=(pos2-real_region[0])/100
                                                                if block1==block2:
                                                                    RD_RealRegion[block1]+=ReadLen
                                                                else:
                                                                    blpos1=pos1-(block1*100+real_region[0])
                                                                    blpos2=pos2-(block2*100+real_region[0])
                                                                    RD_RealRegion[block1]+=(100-blpos1)
                                                                    RD_RealRegion[block2]+=blpos2
                                                            if int(pbam[4])>int(QCAlign):             #fail quality control, skip
                                                                absIL=abs(int(pbam[8]))
                                                                QCFlag=0
                                                                link_flag=0
                                                                if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                                                                    QCFlag+=1
                                                                if DRtemp==['+','-'] and int(pbam[8])<0:
                                                                    QCFlag+=1
                                                                if DRtemp==['+','+'] or DRtemp==['-','-']:
                                                                    QCFlag+=1
                                                                if int(pbam[8])==0:
                                                                   QCFlag+=1
                                                                if not pbam[5].find('S')==-1:
                                                                    QCFlag+=1
                                                                if not QCFlag==0:
                                                                    print>>fmini, pbam1
                                                        fbam.close()
                                                        fmini.close()
                                                        for rfrr in range(len(RD_RealRegion[:-1])):
                                                            RD_RealRegion[rfrr]=str(float(RD_RealRegion[rfrr])/100.0)
                                                        if real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100 ==0:
                                                            del RD_RealRegion[-1]
                                                        else:
                                                            RD_RealRegion[-1]=str(float(RD_RealRegion[-1])/float(real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100))
                                                        print >>fRDind, ' '.join(RD_RealRegion)
                                                        fRDind.close()
                                                        num_of_wind=(int(pcn2[2])-int(pcn2[1]))/100
                                                        coverage={}
                                                        for i in range(num_of_wind+1):
                                                            coverage[i]=[int(pcn2[1])+i*100,int(pcn2[1])+i*100+99,0]
                                                        coverage[-1]=[0,0,0]
                                                        for xyzw in range(len(RD_RealRegion)):
                                                            if xyzw in coverage.keys():
                                                                coverage[xyzw][2]=float(RD_RealRegion[xyzw])
                                                        for j in GC_Content.keys():
                                                            if j in coverage.keys():
                                                                CN2_Region[Chromosome][GC_Content[j][0]].append(coverage[j][-1])
                                            if os.path.isfile(NullPath+temp_Name+'.fa'):
                                                os.system(r'''rm %s'''%(NullPath+temp_Name+'.fa'))
                                            if os.path.isfile(NullPath+temp_Name+'.sam'):
                                                os.system(r'''rm %s'''%(NullPath+temp_Name+'.sam'))
                                            fo=open(GC_index_File,'w')
                                            #print >>fo, ' '.join(chromos)
                                            #print >>fo, ' '.join([str(i) for i in range(101)])
                                            for key_1 in [chrom]:
                                                    for key_2 in range(101):
                                                            print >>fo, ':'.join(['@',','.join(str(j) for j in CN2_Region[key_1][key_2])])
                                            fo.close()
                                            os.system(r'''samtools view -h -Sb %s -o %s'''%(mini_fout_Name,mini_fout_N2))
                                            os.system(r'''samtools sort %s %s'''%(mini_fout_N2,mini_fout_N3))
                                            os.system(r'''samtools index %s '''%(mini_fout_N4))
                                            os.system(r'''rm %s'''%(mini_fout_N2))
                                            os.system(r'''rm %s'''%(mini_fout_Name))
                                        temp_IL_Rec={}
                                        Link_IL_Rec={}
                                        for loc in loc_rec[chrom]: 
                                            #print loc
                                            if not loc[1]-loc[0]>2*10**6:
                                                loc2=[loc]
                                            else:
                                                sublocNum=int(float(loc[1]-loc[0])/10**6)
                                                loc2=[[loc[0],loc[0]+10**6+int(ClusterLen)]]
                                                for slNum in range(sublocNum)[1:]:
                                                    loc2.append([loc[0]+slNum*10**6-int(ClusterLen),loc[0]+(slNum+1)*10**6+int(ClusterLen)])
                                                if sublocNum>1:
                                                    loc2.append([loc[0]+(slNum+1)*10**6-int(ClusterLen),loc[1]])
                                            for real_region in loc2:
                                                #if real_region[1]-real_region[0]<1000:continue
                                                #else:
                                                    #print real_region        
                                                    fbam=os.popen('''samtools view %s %s:%d-%d'''%(mini_fout_N4,chrom,real_region[0],real_region[1]))
                                                    abInfF={}
                                                    abInfR={}
                                                    abLink={}
                                                    abInf3={}
                                                    LinkFR={}
                                                    LinkFF={}
                                                    LinkRR={}
                                                    LinkRF={}
                                                    LinkNM={}
                                                    LinkNM['F']=[]
                                                    LinkNM['R']=[]
                                                    LinkSP={}
                                                    abLinkSP={}
                                                    test=[]
                                                    LinkSPTemp={}
                                                    while True:
                                                        pbam=fbam.readline().strip().split()
                                                        if not pbam: break
                                                        if int(pbam[4])<int(QCAlign): continue             #fail quality control, skip
                                                        if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                        DRtemp=Reads_Direction_Detect(pbam[1])
                                                        ReadLen=cigar2reaadlength(pbam[5])
                                                        absIL=abs(int(pbam[8]))
                                                        signIL=0
                                                        posF=[]
                                                        if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                                                                signIL+=1
                                                        if not pbam[5].find('S')==-1 or not pbam[5].find('H')==-1:
                                                            ClippedLen=cigar2splitlen(pbam[5])
                                                            ClippedPos=cigar2split(pbam[5])
                                                            ClippedQual=cigar2splitqual(pbam[5],pbam[10])
                                                            if ClippedQual>QCSplit:
                                                                    ClipAbsPos=[]
                                                                    for c in range(len(ClippedLen)):
                                                                        if ClippedLen[c]>SPLCff:
                                                                                pos1=int(pbam[3])+ClippedPos[c]
                                                                                posF.append(pos1)
                                                                                pos2=int(pbam[7])
                                                                                if DRtemp[0]=='+':
                                                                                    if not pos1 in abInfF.keys():
                                                                                        abInfF[pos1]=[0,0,0,1]
                                                                                    else:
                                                                                        abInfF[pos1][3]+=1
                                                                                elif DRtemp[0]=='-':
                                                                                    if not pos1 in abInfR.keys():
                                                                                        abInfR[pos1]=[0,0,0,1]
                                                                                    else:
                                                                                        abInfR[pos1][3]+=1
                                                            if not pbam[0] in LinkSPTemp.keys():
                                                                LinkSPTemp[pbam[0]]=[[]]
                                                            else:
                                                                LinkSPTemp[pbam[0]]+=[[]]
                                                            for c in ClippedPos:
                                                                pos1=int(pbam[3])+c
                                                                LinkSPTemp[pbam[0]][-1]+=[pos1,DRtemp[0]]
                                                            LinkSPTemp[pbam[0]].append('S')
                                                        #else:
                                                        if posF==[]:
                                                            if DRtemp[0]=='+':
                                                                pos1=int(pbam[3])+len(pbam[9])
                                                            elif DRtemp[0]=='-':
                                                                pos1=int(pbam[3])
                                                        else:
                                                            pos1=posF[0]
                                                        pos2=int(pbam[7])
                                                        if pbam[0] in LinkSPTemp.keys():
                                                            LinkSPTemp[pbam[0]]+=[[pos1,DRtemp[0]]]
                                                        else:
                                                            LinkSPTemp[pbam[0]]=[[pos1,DRtemp[0]]]
                                                        if DRtemp==['+','-'] and int(pbam[8])>0:
                                                                if signIL>0:
                                                                        #pos1=int(pbam[3])+ReadLen
                                                                        #pos2=int(pbam[7])
                                                                        if not pos1 in LinkFR.keys():
                                                                                LinkFR[pos1]=[pos2]
                                                                        else:
                                                                                LinkFR[pos1]+=[pos2]
                                                                        if not pos1 in abInfF.keys():
                                                                                abInfF[pos1]=[1,0,0,0]
                                                                        else:
                                                                                abInfF[pos1][0]+=1
                                                                        if not pos2 in abInfR.keys():
                                                                           abInfR[pos2]=[1,0,0,0]
                                                                        else:
                                                                                abInfR[pos2][0]+=1
                                                        elif DRtemp==['+','-'] and int(pbam[8])<0:
                                                                #pos1=int(pbam[3])+ReadLen
                                                                #pos2=int(pbam[7])
                                                                if not pos1 in LinkFR.keys():
                                                                        LinkFR[pos1]=[pos2]
                                                                else:
                                                                        LinkFR[pos1]+=[pos2]
                                                                if not pos1 in abInfF.keys():
                                                                        abInfF[pos1]=[0,1,0,0]
                                                                else:
                                                                        abInfF[pos1][1]+=1
                                                                if not pos2 in abInfR.keys():
                                                                   abInfR[pos2]=[0,1,0,0]
                                                                else:
                                                                        abInfR[pos2][1]+=1
                                                                if signIL>0:
                                                                        abInfF[pos1][0]+=1
                                                                        abInfR[pos2][0]+=1
                                                        elif DRtemp==['+','+'] and not int(pbam[8])==0:
                                                                #pos1=int(pbam[3])+ReadLen
                                                                if not pos1 in abInfF.keys():
                                                                        abInfF[pos1]=[0,1,0,0]
                                                                else:
                                                                        abInfF[pos1][1]+=1
                                                                if signIL>0:
                                                                        abInfF[pos1][0]+=1
                                                                if not pbam[0] in abtLink.keys():
                                                                        abtLink[pbam[0]]=[pos1]
                                                                elif pbam[0] in abtLink.keys():
                                                                        abtLink[pbam[0]].append(pos1)
                                                                        if not min(abtLink[pbam[0]]) in LinkFF.keys():
                                                                                LinkFF[min(abtLink[pbam[0]])]=[max(abtLink[pbam[0]])]
                                                                        else:
                                                                                LinkFF[min(abtLink[pbam[0]])]+=[max(abtLink[pbam[0]])]
                                                                        del abtLink[pbam[0]]
                                                        elif DRtemp==['-','-'] and int(pbam[8])>0:
                                                                #pos1=int(pbam[3])
                                                                #pos2=int(pbam[7])
                                                                if not min(pos1,pos2) in LinkRR.keys():
                                                                        LinkRR[min(pos1,pos2)]=[max(pos1,pos2)]
                                                                else:
                                                                        LinkRR[min(pos1,pos2)]+=[max(pos1,pos2)]
                                                                if not pos1 in abInfR.keys():
                                                                        abInfR[pos1]=[0,1,0,0]
                                                                else:
                                                                        abInfR[pos1][1]+=1                                    
                                                                if not pos2 in abInfR.keys():
                                                                        abInfR[pos2]=[0,1,0,0]
                                                                else:
                                                                        abInfR[pos2][1]+=1
                                                                if signIL>0:
                                                                        abInfR[pos1][0]+=1
                                                                        abInfR[pos2][0]+=1
                                                        elif int(pbam[8])==0:
                                                                if int(pbam[1])&8>0:
                                                                        if DRtemp[0]=='+':
                                                                                #pos1=int(pbam[3])+ReadLen
                                                                                LinkNM['F'].append(pos1)
                                                                                if pos1>real_region[0] and pos1<real_region[1]:
                                                                                        if not pos1 in abInfF.keys():
                                                                                                abInfF[pos1]=[0,0,1,0]
                                                                                        else:
                                                                                                abInfF[pos1][2]+=1
                                                                        elif DRtemp[0]=='-':
                                                                                #pos1=int(pbam[3])
                                                                                LinkNM['R'].append(pos1)
                                                                                if pos1>real_region[0] and pos1<real_region[1]:
                                                                                        if not pos1 in abInfR.keys():
                                                                                                abInfR[pos1]=[0,0,1,0]
                                                                                        else:
                                                                                                abInfR[pos1][2]+=1
                                                                if not pbam[6]=='=':
                                                                        if DRtemp[0]=='+':
                                                                                #pos1=int(pbam[3])+ReadLen
                                                                                if pos1>real_region[0] and pos1<real_region[1]:
                                                                                        if not pos1 in abLink.keys():
                                                                                                abLink[pos1]=['f',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                        else:
                                                                                                abLink[pos1]+=['f',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                        if not pos1 in abInfF.keys():
                                                                                                abInfF[pos1]=[0,0,1,0]
                                                                                        else:
                                                                                                abInfF[pos1][2]+=1
                                                                        elif DRtemp[0]=='-':
                                                                                #pos1=int(pbam[3])
                                                                                if pos1>real_region[0] and pos1<real_region[1]:
                                                                                        if not pos1 in abLink.keys():
                                                                                                abLink[pos1]=['r',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                        else:
                                                                                                abLink[pos1]+=['r',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                        if not pos1 in abInfR.keys():
                                                                                                abInfR[pos1]=[0,0,1,0]
                                                                                        else:
                                                                                                abInfR[pos1][2]+=1
                                                    for k1 in LinkSPTemp.keys():
                                                        if not 'S' in LinkSPTemp[k1]:
                                                            del LinkSPTemp[k1]
                                                        else:
                                                            if len(LinkSPTemp[k1])==2:
                                                                del LinkSPTemp[k1]
                                                    for k1 in LinkSPTemp.keys():
                                                        tempk1=[]
                                                        for k2 in LinkSPTemp[k1]:
                                                            if not k2=='S':
                                                                tempk1+=[k2]
                                                        for k2 in range(len(tempk1[0])/2):
                                                            for k3 in range(len(tempk1[1])/2):
                                                                if tempk1[0][2*k2]<tempk1[1][2*k3]:
                                                                    tempk2=[tempk1[0][2*k2],tempk1[0][2*k2+1],tempk1[1][2*k3],tempk1[1][2*k3+1]]
                                                                    if [tempk2[1],tempk2[3]]==['+','-']:
                                                                        if not tempk2[0] in LinkFR.keys():
                                                                            LinkFR[tempk2[0]]=[]
                                                                        LinkFR[tempk2[0]].append(tempk2[2])
                                                                    if [tempk2[1],tempk2[3]]==['+','+']:
                                                                        if not tempk2[0] in LinkFF.keys():
                                                                            LinkFF[tempk2[0]]=[]
                                                                        LinkFF[tempk2[0]].append(tempk2[2])
                                                                    if [tempk2[1],tempk2[3]]==['-','-']:
                                                                        if not tempk2[0] in LinkRR.keys():
                                                                            LinkRR[tempk2[0]]=[]
                                                                        LinkRR[tempk2[0]].append(tempk2[2])
                                                                    if [tempk2[1],tempk2[3]]==['-','+']:
                                                                        if not tempk2[0] in LinkRF.keys():
                                                                            LinkRF[tempk2[0]]=[]
                                                                        LinkRF[tempk2[0]].append(tempk2[2])
                                                    for k1 in abInfF.keys():
                                                            if not abInfF[k1][-1]==0:
                                                                    if not k1 in LinkSP.keys():
                                                                            LinkSP[k1]=0
                                                                    LinkSP[k1]+=abInfF[k1][-1]
                                                    for k1 in abInfR.keys():
                                                            if not abInfR[k1][-1]==0:
                                                                    if not k1 in LinkSP.keys():
                                                                            LinkSP[k1]=0
                                                                    LinkSP[k1]+=abInfR[k1][-1]                            
                                                    for k1 in LinkFR.keys():
                                                            for k2 in LinkFR[k1]:
                                                                    if not k2 in LinkRF.keys():
                                                                            LinkRF[k2]=[]
                                                                    LinkRF[k2].append(k1)
                                                    out_pair_bp=[]
                                                    out_single_bp=[]
                                                    SP4S=[]
                                                    if not LinkSP=={}:
                                                        SP2S=[]
                                                        linkSP_First=clusterNums(sorted(LinkSP.keys()), SPCluLen, 'f')
                                                        for k1 in range(len(linkSP_First[0])):
                                                                if linkSP_First[1][k1]==1:
                                                                        if linkSP_First[0][k1] in abInfF.keys():
                                                                                if abInfF[linkSP_First[0][k1]]==[0,0,0,1]:
                                                                                        del abInfF[linkSP_First[0][k1]]
                                                                        if linkSP_First[0][k1] in abInfR.keys():
                                                                                if abInfR[linkSP_First[0][k1]]==[0,0,0,1]:
                                                                                        del abInfR[linkSP_First[0][k1]]
                                                        linkSPF=clusterNums(sorted(LinkSP.keys()), SPCluLen, 'f')[0]
                                                        linkSPR=clusterNums(sorted(LinkSP.keys()), SPCluLen, 'r')[0]
                                                        linkSPFR=clusterSupVis2(sorted(LinkSP.keys()),sorted(linkSPR),sorted(linkSPF),'left')
                                                        for k1 in linkSPFR.keys():
                                                                qual_num=0
                                                                qual_rec=[]
                                                                for k2 in linkSPFR[k1]:
                                                                        qual_num+=LinkSP[k2]
                                                                        qual_rec.append(LinkSP[k2])
                                                                if qual_num<SPCluMin:
                                                                        del linkSPFR[k1]
                                                                else:
                                                                        SP2S.append(linkSPFR[k1][qual_rec.index(max(qual_rec))])
                                                        if not SP2S==[]:
                                                                SP3S=[]
                                                                key=[]
                                                                if align_QCflag==1:
                                                                    for key1 in SP2S:
                                                                        QCLNSP_flag=0
                                                                        for aqb in [key1]:
                                                                            if aqb>BPAlignQCFlank:
                                                                                LNSPFR_aqb=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,aqb-BPAlignQCFlank,aqb+BPAlignQCFlank))
                                                                                LNSPFR_score=float(LNSPFR_aqb.read().strip())  
                                                                                if LNSPFR_score<float(BPAlignQC):
                                                                                        QCLNSP_flag+=1
                                                                        if not QCLNSP_flag==0:
                                                                                SP3S.append(key1)
                                                                SP4S=[]
                                                                for k1 in SP2S:
                                                                        if not k1 in SP3S:
                                                                                SP4S.append(k1)
                                                                SP4S.sort()
                                                                for i in range(len(SP4S)-1):
                                                                        if SP4S[i+1]-SP4S[i]<10:
                                                                                SP4S[i]=SP4S[i+1]
                                                                SP5S=[]
                                                                for i in SP4S:
                                                                        if not i in SP5S:
                                                                                SP5S.append(i)
                                                                SP4S=SP5S
                                                                #SP5SF=clusterNums(sorted(SP4S), 100, 'f')[0]
                                                                #SP5SR=clusterNums(sorted(SP4S), 100, 'r')[0]
                                                                #SP5FR=clusterSupVis2(sorted(SP4S),sorted(SP5SR),sorted(SP5SF),'left')
                                                                #SP6S=[]
                                                                #for k1 in sorted(SP5FR.keys()):
                                                                #        temp=[]
                                                                #        for k2 in SP5FR[k1]:
                                                                #                temp.append(LinkSP[k2])
                                                                #        SP6S.append(SP5FR[k1][temp.index(max(temp))])
                                                                #SP4S=SP6S
                                                                if not SP4S==[]:
                                                                        LinkSPF2=clusterSupVis2(sorted(abInfF.keys()), [i-ClusterLen for i in sorted(SP4S)], [i+10 for i in sorted(SP4S)],'right')
                                                                        for k1x in LinkSPF2.keys():
                                                                            key1=k1x-10
                                                                            LinkSPF2[key1]=[i for i in LinkSPF2[k1x]]
                                                                            del LinkSPF2[k1x]
                                                                        for k1x in LinkSPF2.keys():
                                                                            test1=bp_subgroup(LinkSPF2[k1x],Min_Distinguish_Len)
                                                                            if len(test1)>1:
                                                                                for k2x in test1[:-1]:
                                                                                    new_core_ele=0
                                                                                    for k3x in k2x:
                                                                                        if k3x in abInfF.keys():
                                                                                            new_core_ele+=numpy.sum(abInfF[k3x])
                                                                                    if new_core_ele>BPLNCff:
                                                                                        #new_core.append(max(k2x))
                                                                                        if not max(k2x) in LinkSPF2.keys():
                                                                                            LinkSPF2[max(k2x)]=k2x
                                                                                        else:
                                                                                            LinkSPF2[max(k2x)]+=k2x
                                                                                        if not max(k2x) in SP4S:
                                                                                            SP4S.append(max(k2x))
                                                                        LinkSPR2=clusterSupVis2(sorted(abInfR.keys()), [i-10 for i in sorted(SP4S)], [i+ClusterLen for i in sorted(SP4S)],'left')
                                                                        for k1x in LinkSPR2.keys():
                                                                            key1=k1x+10
                                                                            LinkSPR2[key1]=[i for i in LinkSPR2[k1x]]
                                                                            del LinkSPR2[k1x]
                                                                        for k1x in LinkSPR2.keys():
                                                                            test1=bp_subgroup(LinkSPR2[k1x],Min_Distinguish_Len)
                                                                            if len(test1)>1:
                                                                                for k2x in test1[1:]:
                                                                                    new_core_ele=0
                                                                                    for k3x in k2x:
                                                                                        if k3x in abInfR.keys():
                                                                                            new_core_ele+=numpy.sum(abInfR[k3x])
                                                                                    if new_core_ele>BPLNCff:
                                                                                        #new_core.append(min(k2x))
                                                                                        if not min(k2x) in LinkSPR2.keys():
                                                                                            LinkSPR2[min(k2x)]=k2x
                                                                                        else:
                                                                                            LinkSPR2[min(k2x)]+=k2x
                                                                                        if not min(k2x) in SP4S:
                                                                                            SP4S.append(min(k2x))
                                                                        for k1x in LinkSPF2.keys():
                                                                            temp_rec=LinkSPF2[k1x]
                                                                            LinkSPF2[k1x]=[LinkSPF2[k1x],[]]
                                                                            for k2 in temp_rec:
                                                                                if k2 in LinkFR.keys():
                                                                                    LinkSPF2[k1x][1]+=LinkFR[k2]
                                                                                if k2 in LinkFF.keys():
                                                                                    LinkSPF2[k1x][1]+=LinkFF[k2]
                                                                        for k1x in LinkSPR2.keys():
                                                                            temp_rec=LinkSPR2[k1x]
                                                                            LinkSPR2[k1x]=[LinkSPR2[k1x],[]]
                                                                            for k2 in temp_rec:
                                                                                if k2 in LinkRR.keys():
                                                                                    LinkSPR2[k1x][1]+=LinkRR[k2]
                                                                                if k2 in LinkRF.keys():
                                                                                    LinkSPR2[k1x][1]+=LinkRF[k2]
                                                                        LinkSP_To_Link={}
                                                                        for k1x in SP4S:
                                                                            if k1x in LinkSPF2.keys():
                                                                                LinkSP_To_Link[k1x]=[[],[]]
                                                                                LinkSP_To_Link[k1x][0]+=LinkSPF2[k1x][0]
                                                                                LinkSP_To_Link[k1x][1]+=LinkSPF2[k1x][1]
                                                                            if k1x in LinkSPR2.keys():
                                                                                if not k1x in LinkSP_To_Link.keys():
                                                                                    LinkSP_To_Link[k1x]=[[],[]]
                                                                                LinkSP_To_Link[k1x][0]+=LinkSPR2[k1x][0]
                                                                                LinkSP_To_Link[k1x][1]+=LinkSPR2[k1x][1]
                                                                            if k1x in LinkSP_To_Link.keys():
                                                                                LinkSP_To_Link[k1x][0].sort()
                                                                                LinkSP_To_Link[k1x][1].sort()
                                                                        for k1x in sorted(LinkSP_To_Link.keys()):
                                                                            for k2 in LinkSP_To_Link[k1x][1]:
                                                                                if k2 in LinkSP_To_Link.keys():
                                                                                    if not sorted([k1x,k2]) in out_pair_bp and not k1x==k2:
                                                                                        out_pair_bp.append(sorted([k1x,k2]))
                                                                                    elif k1x==k2:
                                                                                        if not k1x in out_single_bp:
                                                                                            out_single_bp.append(k1x)
                                                                        for k1x in sorted(LinkSP_To_Link.keys()):
                                                                            if not LinkSP_To_Link[k1x][1]==[]:
                                                                                for k2 in sorted(LinkSP_To_Link.keys())[sorted(LinkSP_To_Link.keys()).index(k1x):]:
                                                                                    if not LinkSP_To_Link[k2][1]==[]:
                                                                                        if k2==k1x: 
                                                                                            if not k1x in out_single_bp:
                                                                                                out_single_bp.append(k1x)
                                                                                        else:
                                                                                            overlap_rec=[overlap_calcu(LinkSP_To_Link[k1x][1],LinkSP_To_Link[k2][0]),overlap_calcu(LinkSP_To_Link[k1x][0],LinkSP_To_Link[k2][1])]
                                                                                            if overlap_rec[0]+overlap_rec[1]>0:
                                                                                                if not sorted([k1x,k2]) in out_pair_bp and not k1x==k2:
                                                                                                    out_pair_bp.append(sorted([k1x,k2]))
                                                                                                elif k1x==k2:
                                                                                                    out_single_bp.append(k1x)
                                                                        for k1x in sorted(LinkSP_To_Link.keys()):
                                                                            if LinkSP_To_Link[k1x][1]==[]:
                                                                                out_single_bp.append(k1x)
                                                                                del LinkSP_To_Link[k1x]
                                                                        for k1x in out_pair_bp:
                                                                            for k2 in k1x:
                                                                                if k2 in out_single_bp:
                                                                                    del out_single_bp[out_single_bp.index(k2)]
                                                                                if k2 in LinkSP_To_Link.keys():
                                                                                    del LinkSP_To_Link[k2]
                                                                        for k1x in out_single_bp:
                                                                            if k1x in LinkSP_To_Link.keys():
                                                                                del LinkSP_To_Link[k1x]
                                                                SP4S.sort()
                                                    LinkSP_To_Link={}
                                                    tempIL={}
                                                    if not abInfF=={}:
                                                            clu_a_F=clusterNums(abInfF.keys(), ClusterLen, 'f')[0]
                                                            if not clu_a_F==[]:
                                                                    clu_b_F=clusterNums(abInfF.keys(), ClusterLen, 'r')[0]
                                                                    clu_c_F=clusterSupVis2(sorted(abInfF.keys()), clu_b_F, [caf+10 for caf in clu_a_F],'right')
                                                                    if not clu_c_F=={}:
                                                                            for key2 in clu_c_F.keys():
                                                                                    key2b=key2-10
                                                                                    record=0
                                                                                    record2=0
                                                                                    for key3 in clu_c_F[key2]:
                                                                                            if not key3 >key2b:
                                                                                                    record+=sum(abInfF[key3][:3])
                                                                                            if abs(key2b-key3)<SPCluLen:
                                                                                                    record2+=abInfF[key3][3]
                                                                                    if not record+record2<LinkCluMin:
                                                                                            tempIL[key2b]=[[record,record2,'f'],clu_c_F[key2]]
                                                                                    del clu_c_F[key2]
                                                    if not abInfR=={}:
                                                            clu_a_R=clusterNums(abInfR.keys(), ClusterLen, 'r')[0]
                                                            if not clu_a_R==[]:
                                                                    clu_b_R=clusterNums(abInfR.keys(), ClusterLen, 'f')[0]
                                                                    clu_c_R=clusterSupVis2(sorted(abInfR.keys()), [car-10 for car in clu_a_R], clu_b_R,'left')
                                                                    if not clu_c_R=={}:
                                                                            for key2 in clu_c_R.keys():  
                                                                                    key2b=key2+10
                                                                                    record=0
                                                                                    record2=0
                                                                                    for key3 in clu_c_R[key2]:
                                                                                            if not key3<key2b:
                                                                                                    record+= sum(abInfR[key3][:3])   
                                                                                            if abs(key3-key2b)<SPCluLen:
                                                                                                    record2+=abInfR[key3][3]
                                                                                    if not record+record2<LinkCluMin:
                                                                                            tempIL[key2b]=[[record,record2,'r'],clu_c_R[key2]]
                                                                                    del clu_c_R[key2]
                                                    if not tempIL=={}:
                                                        for aqb in tempIL.keys():
                                                            if aqb<BPAlignQCFlank:
                                                                del tempIL[aqb]
                                                                continue
                                                            if align_QCflag==1:
                                                                LNSPFR_aqb=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,aqb-BPAlignQCFlank,aqb+BPAlignQCFlank))
                                                                tPairF_b=float(LNSPFR_aqb.read().strip()) 
                                                                if tPairF_b<float(BPAlignQC):
                                                                    del tempIL[aqb]
                                                    LinkIL={}
                                                    for k1 in tempIL.keys():
                                                        temp_mate_F={}
                                                        temp_mate_R={}
                                                        info_mate=0
                                                        if tempIL[k1][0][2]=='f':
                                                            for k2 in tempIL[k1][1]:
                                                                if k2 in LinkFF.keys():
                                                                    for k3 in LinkFF[k2]:
                                                                        if k3 in abInfF.keys():
                                                                            temp_mate_F[k3]=sum(abInfF[k3])
                                                                            #info_mate+=sum(abInfF[k3])
                                                                if k2 in LinkFR.keys():
                                                                    for k3 in LinkFR[k2]:
                                                                        if k3 in abInfR.keys():
                                                                            temp_mate_R[k3]=sum(abInfR[k3])
                                                                            #info_mate+=sum(abInfR[k3])
                                                        elif tempIL[k1][0][2]=='r':
                                                            for k2 in tempIL[k1][1]:
                                                                if k2 in LinkRF.keys():
                                                                    for k3 in LinkRF[k2]:
                                                                        if k3 in abInfF.keys():
                                                                            temp_mate_F[k3]=sum(abInfF[k3])
                                                                            #info_mate+=sum(abInfF[k3])
                                                                if k2 in LinkRR.keys():
                                                                    for k3 in LinkRR[k2]:
                                                                        if k3 in abInfR.keys():
                                                                            temp_mate_R[k3]=sum(abInfR[k3])
                                                                            #info_mate+=sum(abInfR[k3])
                                                        for k1x in temp_mate_F.keys():
                                                            info_mate+=temp_mate_F[k1x]
                                                        for k1x in temp_mate_R.keys():
                                                            info_mate+=temp_mate_R[k1x]
                                                        if not info_mate<LinkCluMin:
                                                            LinkIL[k1]=[[],[]]
                                                            if not temp_mate_F=={}:
                                                                LinkIL[k1][0]=clusterQC(clusterNums4(temp_mate_F, ClusterLen, 'f'),LinkCluMin)
                                                            if not temp_mate_R=={}:
                                                                LinkIL[k1][1]=clusterQC(clusterNums4(temp_mate_R, ClusterLen, 'r'),LinkCluMin)
                                                            else:continue
                                                    for k1 in LinkIL.keys():
                                                            for k2 in LinkIL[k1]:
                                                                    for k3 in k2:
                                                                        if k3>BPAlignQCFlank:
                                                                            if align_QCflag==1:
                                                                                tPairF_QC=0
                                                                                tPairF_a=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,k3-BPAlignQCFlank,k3+BPAlignQCFlank))
                                                                                tPairF_b=float(tPairF_a.read().strip())  
                                                                                if not tPairF_b<float(BPAlignQC):
                                                                                        if not [min([k3,k1]),max([k3,k1])] in out_pair_bp and not k1==k3:
                                                                                            out_pair_bp.append([min([k3,k1]),max([k3,k1])])
                                                                                        elif k1==k3:
                                                                                            out_single_bp.append(k1)
                                                                            else:
                                                                                        if not [min([k3,k1]),max([k3,k1])] in out_pair_bp and not k1==k3:
                                                                                                out_pair_bp.append([min([k3,k1]),max([k3,k1])])        
                                                                                        elif  k1==k3:     
                                                                                            out_single_bp.append(k1)                                  
                                                    if not out_pair_bp==[]:
                                                        temp_out_pair_bp=[]
                                                        out_BPmodify={}
                                                        for k1 in out_pair_bp:
                                                            for k2 in k1:
                                                                out_BPmodify[k2]=[]
                                                        if not SP4S==[]:
                                                            LBSP_tempIL=clusterSupVis3(sorted(SP4S),sorted(out_BPmodify.keys()))
                                                            for k1 in out_pair_bp:
                                                                temp_k1=[]
                                                                for k2 in k1:
                                                                    k3=LBSP_tempIL[k2]
                                                                    if abs(k2-k3)<ClusterLen:
                                                                        temp_k1.append(k3)
                                                                    else:
                                                                        temp_k1.append(k2)
                                                                temp_out_pair_bp.append(temp_k1)
                                                        out_pair_bp=temp_out_pair_bp
                                                        for k1 in out_pair_bp:
                                                            for k2 in k1:
                                                                if k2 in SP4S:
                                                                    del SP4S[SP4S.index(k2)]
                                                        out_pair_modify={}
                                                        for i in out_pair_bp:
                                                            if not i[0] in out_pair_modify.keys():
                                                                out_pair_modify[i[0]]=[]
                                                            if not i[1] in out_pair_modify[i[0]]:
                                                                out_pair_modify[i[0]].append(i[1])
                                                        if len(out_pair_modify)>1:
                                                            while True: 
                                                                if len(out_pair_modify)==1: break
                                                                out_pair_qc=[]
                                                                for i in range(len(sorted(out_pair_modify.keys()))-1):
                                                                    out_pair_qc.append(sorted(out_pair_modify.keys())[i+1]-sorted(out_pair_modify.keys())[i])
                                                                if min(out_pair_qc)>50:break
                                                                else:
                                                                    out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))+1]]+=out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                                                    del out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                                        for k1 in out_pair_modify.keys():
                                                            while True:
                                                                if len(out_pair_modify[k1])==1: break
                                                                out_pair_modify[k1].sort()
                                                                out_pair_qc=[]
                                                                for i in range(len(out_pair_modify[k1])-1):
                                                                    out_pair_qc.append(out_pair_modify[k1][i+1]-out_pair_modify[k1][i])
                                                                if min(out_pair_qc)>50:break
                                                                else:
                                                                    out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))+1]=out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))]
                                                                for i in out_pair_modify[k1]:
                                                                    if out_pair_modify[k1].count(i)>1:
                                                                        for j in range(out_pair_modify[k1].count(i)-1):
                                                                            del out_pair_modify[k1][out_pair_modify[k1].index(i)]
                                                        out_pair_numrec={}
                                                        for k1 in out_pair_modify.keys():
                                                            out_pair_numrec[k1]=[]
                                                            for k2 in [k1]+out_pair_modify[k1]:
                                                                if k2 in tempIL.keys() and not k2 in LinkSP.keys():
                                                                    out_pair_numrec[k1].append(len(tempIL[k2][1]))
                                                                elif k2 in LinkSP.keys() and not k2 in tempIL.keys():
                                                                    out_pair_numrec[k1].append(LinkSP[k2])
                                                                elif k2 in LinkSP.keys() and k2 in tempIL.keys():
                                                                    out_pair_numrec[k1].append(LinkSP[k2]+len(tempIL[k2][1]))
                                                                else:
                                                                    out_pair_numrec[k1].append(0)
                                                        fout=open(BPOutputb,'a')
                                                        for i in sorted(out_pair_modify.keys()):
                                                            for j in out_pair_modify[i]:                            
                                                                print >>fout, ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                                #print ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                                if i in tempIL.keys():
                                                                    del tempIL[i]
                                                                if j in tempIL.keys():
                                                                    del tempIL[j]
                                                        fout.close()
                                                        fout=open(BPOutputa,'a')
                                                        for i in sorted(out_single_bp+LinkSP_To_Link.keys()):
                                                            num=0
                                                            if i in abInfF.keys():
                                                                num+=sum(abInfF[i])
                                                            if i in abInfR.keys():
                                                                num+=sum(abInfR[i])
                                                            print >>fout, ' '.join([str(j) for j in [i,num]])
                                                            #print ' '.join([str(j) for j in [i,num]])
                                                        fout.close()
                                                    for i in tempIL.keys():
                                                        temp_IL_Rec[i]=tempIL[i]
                                                        temp_mate_F=[]
                                                        temp_mate_R=[]
                                                        for k2 in tempIL[i][1]:
                                                            if k2 in LinkFF.keys():
                                                                    temp_mate_F+=LinkFF[k2]
                                                            if k2 in LinkFR.keys():
                                                                    temp_mate_R+=LinkFR[k2]
                                                        temp_IL_Rec[i].append(temp_mate_F) 
                                                        temp_IL_Rec[i].append(temp_mate_R) 
                                                    tempLNF=[]
                                                    tempLNR=[]
                                                    for key in abLink.keys():
                                                            for key2 in range(len(abLink[key])/3):
                                                                    if abLink[key][3*key2]=='f':
                                                                            tempLNF.append(key)
                                                                    else:
                                                                            tempLNR.append(key)
                                                    for key in abLinkSP.keys():
                                                            for key2 in range(len(abLinkSP[key])/3):
                                                                    if abLinkSP[key][3*key2]=='f':
                                                                            tempLNF.append(key)
                                                                    else:
                                                                            tempLNR.append(key)
                                                    FtLNFR=clusterNums(tempLNF+tempLNR,ClusterLen,'f')[0]
                                                    RtLNFR=clusterNums(tempLNF+tempLNR,ClusterLen,'r')[0]
                                                    FRtLNFR=clusterSupVis2(sorted(tempLNF+tempLNR),RtLNFR,FtLNFR,'left')
                                                    for key1 in FRtLNFR.keys():
                                                            t_LNFR=[]
                                                            for key2 in FRtLNFR[key1]:
                                                                    if key2 in abLink.keys():
                                                                            t_LNFR+=abLink[key2]
                                                                    else:
                                                                            if abs(key2-key1)<SPCluLen or abs(key2-max(FRtLNFR[key1]))<SPCluLen:
                                                                                    t_LNFR+=abLinkSP[key2]
                                                            if len(t_LNFR)/3<LinkCluMin:
                                                                    del FRtLNFR[key1]
                                                            else:
                                                                if not t_LNFR in FRtLNFR[key1]:
                                                                    FRtLNFR[key1].append(t_LNFR)
                                                    FRtLNFRb=[]
                                                    for key1 in FRtLNFR.keys():
                                                        t1_LN=[]
                                                        t2_LN=[] 
                                                        t3_LN={}
                                                        t4out=[]
                                                        for key2 in FRtLNFR[key1][:-1]:
                                                                if key2 in abLink.keys():
                                                                    if not key2 in t1_LN:
                                                                        for key3 in range(len(abLink[key2])/3):
                                                                            if not abLink[key2][3*key3:3*(key3+1)] in t2_LN:
                                                                                t2_LN.append(abLink[key2][3*key3:3*(key3+1)])
                                                                                t1_LN.append(key2)
                                                        for key2 in t2_LN:
                                                                if not key2[-1].split('_')[0] in t3_LN.keys():
                                                                        t3_LN[key2[-1].split('_')[0]]={}
                                                                        t3_LN[key2[-1].split('_')[0]]['a']=[]
                                                                        t3_LN[key2[-1].split('_')[0]]['b']=[]
                                                                        t3_LN[key2[-1].split('_')[0]]['c']=[]
                                                                        t3_LN[key2[-1].split('_')[0]]['d']=[]
                                                                t3_LN[key2[-1].split('_')[0]]['a'].append(key2[0])
                                                                t3_LN[key2[-1].split('_')[0]]['b'].append(key2[1])
                                                                t3_LN[key2[-1].split('_')[0]]['c'].append(key2[-1].split('_')[1])
                                                                t3_LN[key2[-1].split('_')[0]]['d'].append(key2)                    
                                                        for key2 in t3_LN.keys():
                                                            if clusterQC(clusterNums(t3_LN[key2]['b'], ClusterLen, 'f'), LinkCluMin)==[]:
                                                                    del t3_LN[key2]
                                                            else:
                                                                    t4LN=clusterSupVis2(t3_LN[key2]['b'],clusterQC(clusterNums(t3_LN[key2]['b'], ClusterLen, 'r'), LinkCluMin),clusterQC(clusterNums(t3_LN[key2]['b'], ClusterLen, 'f'), LinkCluMin), 'left')                               
                                                                    for key5 in t4LN.keys():
                                                                            t4LNa=[]
                                                                            t4LNb=[]
                                                                            t4LNc=[]
                                                                            t4LNd=[]
                                                                            t4out=[]
                                                                            for key6 in t4LN[key5]:
                                                                                    t4LNa.append(t3_LN[key2]['a'][t3_LN[key2]['b'].index(key6)])
                                                                                    t4LNc.append(t3_LN[key2]['c'][t3_LN[key2]['b'].index(key6)])
                                                                                    t4LNd.append(key6)
                                                                                    t4LNb.append(t1_LN[t2_LN.index(t3_LN[key2]['d'][t3_LN[key2]['b'].index(key6)])])
                                                                            if not 'f' in t4LNa or float(t4LNa.count('r'))/float(t4LNa.count('f'))>5:
                                                                                    t4out+=[chrom,'r',min(t4LNb)]
                                                                            elif not 'r' in t4LNa or float(t4LNa.count('f'))/float(t4LNa.count('r'))>5:
                                                                                    t4out+=[chrom,'f',max(t4LNb)]
                                                                            else:
                                                                                    t4out+=[chrom,min(t4LNb),max(t4LNb)]
                                                                            if not '+' in t4LNc or float(t4LNc.count('-'))/float(t4LNc.count('+'))>5:
                                                                                    t4out+=[key2,'-',min(t4LNd)]
                                                                            elif not '-' in t4LNc or float(t4LNc.count('+'))/float(t4LNc.count('-'))>5:
                                                                                    t4out+=[key2,'+',max(t4LNd)]
                                                                            else:
                                                                                    t4out+=[key2,min(t4LNd),max(t4LNd)]
                                                        if not t4out==[] and not t4out in FRtLNFRb:
                                                                FRtLNFRb.append(t4out)             
                                                    if not FRtLNFRb==[]:
                                                            fout=open(BPOutputd,'a')
                                                            for keyfrt in FRtLNFRb:
                                                                    print >>fout, ' '.join([str(keyfr2) for keyfr2 in keyfrt])
                                                            fout.close()
                                        Link_IL_Rec={}
                                        for k1 in temp_IL_Rec.keys():
                                            temp_mate_F=temp_IL_Rec[k1][2]
                                            temp_mate_R=temp_IL_Rec[k1][3]
                                            if not len(temp_IL_Rec[k1][1])+len(temp_IL_Rec[k1][2])+len(temp_IL_Rec[k1][3])<LnCluMin:
                                                Link_IL_Rec[k1]=[clusterQC(clusterNums(temp_mate_F, ClusterLen, 'f'),LinkCluMin),
                                                                clusterQC(clusterNums(temp_mate_R, ClusterLen, 'r'),LinkCluMin)]
                                                del temp_IL_Rec[k1]
                                            else:continue
                                        for k1 in Link_IL_Rec.keys():
                                            for k2 in Link_IL_Rec[k1]:
                                                    for k3 in k2:
                                                        if k3>BPAlignQCFlank:
                                                            if align_QCflag==1:
                                                                tPairF_QC=0
                                                                tPairF_a=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,k3-BPAlignQCFlank,k3+BPAlignQCFlank))
                                                                tPairF_b=float(tPairF_a.read().strip())  
                                                                if not tPairF_b<float(BPAlignQC):
                                                                        if not [min([k3,k1]),max([k3,k1])] in out_pair_bp:
                                                                                out_pair_bp.append([min([k3,k1]),max([k3,k1])])
                                                            else:
                                                                        if not [min([k3,k1]),max([k3,k1])] in out_pair_bp:
                                                                                out_pair_bp.append([min([k3,k1]),max([k3,k1])])
                                        if not out_pair_bp==[]:
                                            temp_out_pair_bp=[]
                                            out_BPmodify={}
                                            for k1 in out_pair_bp:
                                                for k2 in k1:
                                                    out_BPmodify[k2]=[]
                                            if not SP4S==[]:
                                                LBSP_tempIL=clusterSupVis3(sorted(SP4S),sorted(out_BPmodify.keys()))
                                                for k1 in out_pair_bp:
                                                    temp_k1=[]
                                                    for k2 in k1:
                                                        k3=LBSP_tempIL[k2]
                                                        if abs(k2-k3)<ClusterLen:
                                                            temp_k1.append(k3)
                                                        else:
                                                            temp_k1.append(k2)
                                                    temp_out_pair_bp.append(temp_k1)
                                            out_pair_bp=temp_out_pair_bp
                                            for k1 in out_pair_bp:
                                                for k2 in k1:
                                                    if k2 in SP4S:
                                                        del SP4S[SP4S.index(k2)]
                                            out_pair_modify={}
                                            for i in out_pair_bp:
                                                if not i[0] in out_pair_modify.keys():
                                                    out_pair_modify[i[0]]=[]
                                                if not i[1] in out_pair_modify[i[0]]:
                                                    out_pair_modify[i[0]].append(i[1])
                                            if len(out_pair_modify)>1:
                                                while True: 
                                                    if len(out_pair_modify)==1: break
                                                    out_pair_qc=[]
                                                    for i in range(len(sorted(out_pair_modify.keys()))-1):
                                                        out_pair_qc.append(sorted(out_pair_modify.keys())[i+1]-sorted(out_pair_modify.keys())[i])
                                                    if min(out_pair_qc)>50:break
                                                    else:
                                                        out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))+1]]+=out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                                        del out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                            for k1 in out_pair_modify.keys():
                                                while True:
                                                    if len(out_pair_modify[k1])==1: break
                                                    out_pair_modify[k1].sort()
                                                    out_pair_qc=[]
                                                    for i in range(len(out_pair_modify[k1])-1):
                                                        out_pair_qc.append(out_pair_modify[k1][i+1]-out_pair_modify[k1][i])
                                                    if min(out_pair_qc)>50:break
                                                    else:
                                                        out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))+1]=out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))]
                                                    for i in out_pair_modify[k1]:
                                                        if out_pair_modify[k1].count(i)>1:
                                                            for j in range(out_pair_modify[k1].count(i)-1):
                                                                del out_pair_modify[k1][out_pair_modify[k1].index(i)]
                                            out_pair_numrec={}
                                            for k1 in out_pair_modify.keys():
                                                out_pair_numrec[k1]=[]
                                                for k2 in [k1]+out_pair_modify[k1]:
                                                    if k2 in tempIL.keys() and not k2 in LinkSP.keys():
                                                        out_pair_numrec[k1].append(len(tempIL[k2][1]))
                                                    elif k2 in LinkSP.keys() and not k2 in tempIL.keys():
                                                        out_pair_numrec[k1].append(LinkSP[k2])
                                                    elif k2 in LinkSP.keys() and k2 in tempIL.keys():
                                                        out_pair_numrec[k1].append(LinkSP[k2]+len(tempIL[k2][1]))
                                                    else:
                                                        out_pair_numrec[k1].append(0)
                                            fout=open(BPOutputb,'a')
                                            for i in sorted(out_pair_modify.keys()):
                                                for j in out_pair_modify[i]:                            
                                                    print >>fout, ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                    #print ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                    if i in tempIL.keys():
                                                        del tempIL[i]
                                                    if j in tempIL.keys():
                                                        del tempIL[j]
                                            fout.close()
                                        fout=open(BPOutputa,'a')
                                        for i in sorted(temp_IL_Rec.keys()):
                                            print >>fout, ' '.join([str(i),str(sum(temp_IL_Rec[i][0][:2]))])	
                                        fout.close()
                                        time2=time.time()
                                        #print 'Break point searching done for sample: '+bamF_Name+', chromosome:'+chrF
                                        #print 'Time Consuming:'+str(datetime.timedelta(seconds=(time2-time1)))
                                        os.system(r'''cat %s >> %s'''%(BPOutputd,BPOutpute))
                                        os.system(r'''rm %s'''%(BPOutputd))




def not_in_use2():
                                            for loc in loc_rec[chrom]:  
                                                print loc
                                                Chromosome=chrom
                                                if loc[1]-loc[0]<Length_Limit:continue
                                                else:
                                                    pcn2=[chrom]+loc
                                                    temp_Name='temp.Null1.'+bamF.split('/')[-1]+'.'+chrom
                                                    fasta_file=NullPath+temp_Name+'.fa'
                                                    os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref_file,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),fasta_file))                    
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
                                                fo=open(GC_index_File,'w')
                                                #print >>fo, ' '.join(chromos)
                                                #print >>fo, ' '.join([str(i) for i in range(101)])
                                                for key_1 in [chrom]:
                                                        for key_2 in range(101):
                                                                print >>fo, ':'.join(['@',','.join(str(j) for j in CN2_Region[key_1][key_2])])
                                                fo.close()
                                                if not loc[1]-loc[0]>2*10**6:
                                                    loc2=[loc]
                                                else:
                                                    sublocNum=(loc[1]-loc[0])/10**6
                                                    loc2=[[loc[0],loc[0]+10**6+int(ClusterLen2)]]
                                                    for slNum in range(sublocNum)[1:]:
                                                        loc2.append([loc[0]+slNum*10**6-int(ClusterLen2),loc[0]+(slNum+1)*10**6+int(ClusterLen2)])
                                                    if sublocNum>1:
                                                        loc2.append([loc[0]+(slNum+1)*10**6-int(ClusterLen2),loc[1]])
                                                for real_region in loc2:
                                                    #if real_region[1]-real_region[0]<1000:continue
                                                    #else:
                                                        print real_region        
                                                        fmini=open(mini_fout_Name,'a')
                                                        fRDind=open(RD_index_File,'a')
                                                        #flink=open(mini_link_Name,'a')
                                                        #if loc_flag==1:
                                                        #    print >>fRDind,'str'+chrom+':'+str(real_region[0])+'-'+str(real_region[1])
                                                        #else:
                                                        print >>fRDind,chrom+':'+str(real_region[0])+'-'+str(real_region[1])                
                                                        RD_RealRegion=[0 for i in range((real_region[1]-real_region[0])/100+1)]
                                                        fbam=os.popen('''samtools view %s %s:%d-%d'''%(BamInput,chrom,real_region[0],real_region[1]))
                                                        while True:
                                                            pbam1=fbam.readline().strip()
                                                            if not pbam1: break
                                                            pbam=pbam1.split()
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect(pbam[1])
                                                            ReadLen=cigar2reaadlength(pbam[5])
                                                            pos1=int(pbam[3])
                                                            pos2=int(pbam[3])+ReadLen
                                                            if pos2>real_region[0] and pos1<real_region[1]:
                                                                if pos1<real_region[0] and pos2>real_region[0]:
                                                                    pos1=real_region[0]
                                                                if pos1<real_region[1] and pos2>real_region[1]:
                                                                    pos2=real_region[1]
                                                                block1=(pos1-real_region[0])/100
                                                                block2=(pos2-real_region[0])/100
                                                                if block1==block2:
                                                                    RD_RealRegion[block1]+=ReadLen
                                                                else:
                                                                    blpos1=pos1-(block1*100+real_region[0])
                                                                    blpos2=pos2-(block2*100+real_region[0])
                                                                    RD_RealRegion[block1]+=(100-blpos1)
                                                                    RD_RealRegion[block2]+=blpos2
                                                            if int(pbam[4])>int(QCAlign):             #fail quality control, skip
                                                                absIL=abs(int(pbam[8]))
                                                                QCFlag=0
                                                                link_flag=0
                                                                if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                                                                    QCFlag+=1
                                                                if DRtemp==['+','-'] and int(pbam[8])<0:
                                                                    QCFlag+=1
                                                                if DRtemp==['+','+'] or DRtemp==['-','-']:
                                                                    QCFlag+=1
                                                                if int(pbam[8])==0:
                                                                   QCFlag+=1
                                                                if not pbam[5].find('S')==-1:
                                                                    QCFlag+=1
                                                                if not QCFlag==0:
                                                                    print>>fmini, pbam1
                                                        fbam.close()
                                                        fmini.close()
                                                        for rfrr in range(len(RD_RealRegion[:-1])):
                                                            RD_RealRegion[rfrr]=str(float(RD_RealRegion[rfrr])/100.0)
                                                        if real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100 ==0:
                                                            del RD_RealRegion[-1]
                                                        else:
                                                            RD_RealRegion[-1]=str(float(RD_RealRegion[-1])/float(real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100))
                                                        print >>fRDind, ' '.join(RD_RealRegion)
                                                        fRDind.close()


def not_in_use():
                                        if not os.path.isfile(mini_fout_N4):
                                            os.system(r''' samtools view -H %s -o %s'''%(BamInput,mini_fout_Name)) 
                                            RD_index_Path=NullPath+'RD_Stat/'
                                            if not os.path.isdir(RD_index_Path):
                                                os.system(r'''mkdir %s'''%(RD_index_Path))
                                            RD_index_File=RD_index_Path+bamF_Name+'.'+chrF+'.RD.index'
                                            GC_index_File=RD_index_Path+bamF_Name+'.'+chrF+'.GC.index'
                                            file_start(RD_index_File)
                                            file_start(GC_index_File)
                                            ClusterLen2=int(ClusterLen/100+1)*100
                                            CN2_Region[chrom]={}
                                            for con in range(101):
                                                    CN2_Region[chrom][con]=[]
                                            for loc in loc_rec[chrom]:  
                                                #print loc
                                                Chromosome=chrom
                                                if loc[1]-loc[0]<Length_Limit:continue
                                                if not loc[1]-loc[0]>2*10**6:
                                                    loc2=[loc]
                                                else:
                                                    sublocNum=(loc[1]-loc[0])/10**6
                                                    loc2=[[loc[0],loc[0]+10**6+int(ClusterLen2)]]
                                                    for slNum in range(sublocNum)[1:]:
                                                        loc2.append([loc[0]+slNum*10**6-int(ClusterLen2),loc[0]+(slNum+1)*10**6+int(ClusterLen2)])
                                                    if sublocNum>1:
                                                        loc2.append([loc[0]+(slNum+1)*10**6-int(ClusterLen2),loc[1]])
                                                for real_region in loc2:
                                                    pcn2=[chrom]+real_region
                                                    temp_Name='temp.Null1.'+bamF.split('/')[-1]+'.'+chrom
                                                    fasta_file=NullPath+temp_Name+'.fa'
                                                    os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref_file,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),fasta_file))                    
                                                    Seq1=Fasta_To_Sequence(fasta_file)
                                                    if Seq1=='ERROR!':continue
                                                    if not Seq1=='ERROR!':
                                                        #print pcn2
                                                        sam_file=NullPath+temp_Name+'.sam'
                                                        os.system(r'''samtools view %s %s:%d-%d > %s'''%(bamF,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),sam_file))
                                                        Number_Of_Windows=len(Seq1)/100
                                                        GC_Content={}
                                                        for i in range(len(Seq1)/100+1):                
                                                                Seq2=Seq1[(i)*100:(i+1)*100]
                                                                GC_Content[i]=GC_Content_Calculate(Seq2)
                                                        print real_region        
                                                        fmini=open(mini_fout_Name,'a')
                                                        fRDind=open(RD_index_File,'a')
                                                        #flink=open(mini_link_Name,'a')
                                                        #if loc_flag==1:
                                                        #    print >>fRDind,'str'+chrom+':'+str(real_region[0])+'-'+str(real_region[1])
                                                        #else:
                                                        print >>fRDind,chrom+':'+str(real_region[0])+'-'+str(real_region[1])                
                                                        RD_RealRegion=[0 for i in range((real_region[1]-real_region[0])/100+1)]
                                                        fbam=os.popen('''samtools view %s %s:%d-%d'''%(BamInput,chrom,real_region[0],real_region[1]))
                                                        while True:
                                                            pbam1=fbam.readline().strip()
                                                            if not pbam1: break
                                                            pbam=pbam1.split()
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect(pbam[1])
                                                            ReadLen=cigar2reaadlength(pbam[5])
                                                            pos1=int(pbam[3])
                                                            pos2=int(pbam[3])+ReadLen
                                                            if pos2>real_region[0] and pos1<real_region[1]:
                                                                if pos1<real_region[0] and pos2>real_region[0]:
                                                                    pos1=real_region[0]
                                                                if pos1<real_region[1] and pos2>real_region[1]:
                                                                    pos2=real_region[1]
                                                                block1=(pos1-real_region[0])/100
                                                                block2=(pos2-real_region[0])/100
                                                                if block1==block2:
                                                                    RD_RealRegion[block1]+=ReadLen
                                                                else:
                                                                    blpos1=pos1-(block1*100+real_region[0])
                                                                    blpos2=pos2-(block2*100+real_region[0])
                                                                    RD_RealRegion[block1]+=(100-blpos1)
                                                                    RD_RealRegion[block2]+=blpos2
                                                            if int(pbam[4])>int(QCAlign):             #fail quality control, skip
                                                                absIL=abs(int(pbam[8]))
                                                                QCFlag=0
                                                                link_flag=0
                                                                if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                                                                    QCFlag+=1
                                                                if DRtemp==['+','-'] and int(pbam[8])<0:
                                                                    QCFlag+=1
                                                                if DRtemp==['+','+'] or DRtemp==['-','-']:
                                                                    QCFlag+=1
                                                                if int(pbam[8])==0:
                                                                   QCFlag+=1
                                                                if not pbam[5].find('S')==-1:
                                                                    QCFlag+=1
                                                                if not QCFlag==0:
                                                                    print>>fmini, pbam1
                                                        fbam.close()
                                                        fmini.close()
                                                        for rfrr in range(len(RD_RealRegion[:-1])):
                                                            RD_RealRegion[rfrr]=str(float(RD_RealRegion[rfrr])/100.0)
                                                        if real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100 ==0:
                                                            del RD_RealRegion[-1]
                                                        else:
                                                            RD_RealRegion[-1]=str(float(RD_RealRegion[-1])/float(real_region[1]-real_region[0]-(real_region[1]-real_region[0])/100*100))
                                                        print >>fRDind, ' '.join(RD_RealRegion)
                                                        fRDind.close()
                                                        num_of_wind=(int(pcn2[2])-int(pcn2[1]))/100
                                                        coverage={}
                                                        for i in range(num_of_wind+1):
                                                            coverage[i]=[int(pcn2[1])+i*100,int(pcn2[1])+i*100+99,0]
                                                        coverage[-1]=[0,0,0]
                                                        for xyzw in range(len(RD_RealRegion)):
                                                            if xyzw in coverage.keys():
                                                                coverage[xyzw][2]=float(RD_RealRegion[xyzw])
                                                        for j in GC_Content.keys():
                                                            if j in coverage.keys():
                                                                CN2_Region[Chromosome][GC_Content[j][0]].append(coverage[j][-1])
                                            if os.path.isfile(NullPath+temp_Name+'.fa'):
                                                os.system(r'''rm %s'''%(NullPath+temp_Name+'.fa'))
                                            if os.path.isfile(NullPath+temp_Name+'.sam'):
                                                os.system(r'''rm %s'''%(NullPath+temp_Name+'.sam'))
                                            fo=open(GC_index_File,'w')
                                            #print >>fo, ' '.join(chromos)
                                            #print >>fo, ' '.join([str(i) for i in range(101)])
                                            for key_1 in [chrom]:
                                                    for key_2 in range(101):
                                                            print >>fo, ':'.join(['@',','.join(str(j) for j in CN2_Region[key_1][key_2])])
                                            fo.close()
                                            os.system(r'''samtools view -h -Sb %s -o %s'''%(mini_fout_Name,mini_fout_N2))
                                            os.system(r'''samtools sort %s %s'''%(mini_fout_N2,mini_fout_N3))
                                            os.system(r'''samtools index %s '''%(mini_fout_N4))
                                            os.system(r'''rm %s'''%(mini_fout_N2))
                                            os.system(r'''rm %s'''%(mini_fout_Name))
