#!/usr/bin/env python

#!Python
#Usage
#SVelter3.BPIntegrate.py --ppre workdir --ref ref.fa -s sample.bam 
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import numpy
import re
def clusterSupVis2(dataRound, dataLeft, dataRight):
    flag1=0
    flag2=0
    out={}
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
    return out

def SearchClosest(dataRound, hashcentre):
    data1=sorted(dataRound)
    data2=sorted(hashcentre.keys())
    for data in data2:
        if data1[0]-data>0:
            disa=data1[0]-data
            hashcentre[data]=[data1[0],disa]
        else:
            for flag2 in range(len(data1)-1):
                disa=data1[flag2]-data
                disb=data1[flag2+1]-data
                if not disa*disb>0:
                    if numpy.abs(disa)<numpy.abs(disb):
                        hashcentre[data]=[data1[flag2],disa]
                    else:
                        hashcentre[data]=[data1[flag2+1],disb]

def SearchClosest2(DataRou1, DataCen1,cutoff):
    dataRou=sorted(DataRou1)
    dataCen=sorted(DataCen1)
    out={}
    out2=[]
    for k1 in dataCen:
        for k2 in dataRou:
            if k2>k1 and dataRou[dataRou.index(k2)-1]<k1:
                k3=dataRou[dataRou.index(k2)-1]
                if abs(k2-k1)>abs(k3-k1):
                    out[k1]=k3
                else:
                    out[k1]=k2
                if abs(out[k1]-k1)>cutoff:
                    out2.append(k1)
                    del out[k1]
                break
    return [out,out2]

def clusterNums(Data1, ClusterLen, direction):
    if not Data1==[]:
        data=sorted(Data1)
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
    else:
        return [[],[]]

def clusterQC(hash, CluNum):
    out=[]
    for i in range(len(hash[1])):
        if not hash[1][i]<CluNum:
            out.append(hash[0][i])
    return out

def clusterSupVis3(dataRound,dataCen):
    #dataRound and dataCen has to be sorted
    out={}
    if not dataCen==[] and not dataRound==[]:
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

def clusterSupVis(DataRou1,DataCen1,ClusterLen):
    dataRou=sorted(DataRou1)
    dataCen=sorted(DataCen1)
    flag1=0
    flag2=0
    out={}
    if not dataRou==[] and not dataCen==[]:
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

def LNd_hash_readin(chrom_LN):
    LNd_hash={}
    fin=open(bps_in_path+chrom_LN)
    for line in fin:
        pin=line.strip().split()
        temp=sorted([pin[0],pin[3]])
        pin2=[temp[0],pin[pin.index(temp[0])+1],pin[pin.index(temp[0])+2],temp[1],pin[pin.index(temp[1])+1],pin[pin.index(temp[1])+2]]
        if not '_'.join([pin2[0],pin2[3]]) in LNd_hash.keys():
            tempkey='_'.join([pin2[0],pin2[3]])
            LNd_hash[tempkey]=[{},{}]
        if not pin2[1] in LNd_hash[tempkey][0].keys():
            LNd_hash[tempkey][0][pin2[1]]=[]
        LNd_hash[tempkey][0][pin2[1]]+=pin2[4:6]
        if not pin2[2] in LNd_hash[tempkey][0].keys():
            LNd_hash[tempkey][0][pin2[2]]=[]
        LNd_hash[tempkey][0][pin2[2]]+=pin2[4:6]
        if not pin2[4] in LNd_hash[tempkey][1].keys():
            LNd_hash[tempkey][1][pin2[4]]=[]
        LNd_hash[tempkey][1][pin2[4]]+=pin2[1:3]
        if not pin2[5] in LNd_hash[tempkey][1].keys():
            LNd_hash[tempkey][1][pin2[5]]=[]
        LNd_hash[tempkey][1][pin2[5]]+=pin2[1:3]
    fin.close()
    for k1 in LNd_hash.keys():
        for k2 in LNd_hash[k1][0].keys():
            if k2 in ['f','r','+','-']:
                del LNd_hash[k1][0][k2]
            else:
                temp=[]
                for k3 in LNd_hash[k1][0][k2]:
                    if not k3 in ['f','r','+','-'] and not k3 in temp:
                        temp.append(k3)
                LNd_hash[k1][0][k2]=temp
        for k2 in LNd_hash[k1][1].keys():
            if k2 in ['f','r','+','-']:
                del LNd_hash[k1][1][k2]
            else:
                temp=[]
                for k3 in LNd_hash[k1][1][k2]:
                    if not k3 in ['f','r','+','-'] and not k3 in temp:
                        temp.append(k3)
                LNd_hash[k1][1][k2]=temp
    tempall={}
    temp1=[]
    for k1 in LNd_hash.keys():
        temp1=[]
        for k2 in LNd_hash[k1][0].keys():
            temp2=[k2]
            for k3 in LNd_hash[k1][0][k2]:
                for k4 in LNd_hash[k1][0].keys():
                    if not k4==k2 and k3 in LNd_hash[k1][0][k4]:
                        if not k4 in temp2:
                            temp2.append(k4)
            if not temp2==[k2]:
                temp2.sort()
                temp3=[]
                for k3 in temp2:
                    for k4 in LNd_hash[k1][0][k3]:
                        if not k4 in temp3:
                            temp3.append(k4)
                temp3.sort()
                if not [k1.split('_')[0]]+temp2+[k1.split('_')[1]]+temp3 in temp1:
                    temp1.append([k1.split('_')[0]]+temp2+[k1.split('_')[1]]+temp3)
        for k2 in LNd_hash[k1][1].keys():
            temp2=[k2]
            for k3 in LNd_hash[k1][1][k2]:
                for k4 in LNd_hash[k1][1].keys():
                    if not k4==k2 and k3 in LNd_hash[k1][1][k4]:
                        if not k4 in temp2:
                            temp2.append(k4)
            if not temp2==[k2]:
                temp2.sort()
                temp3=[]
                for k3 in temp2:
                    for k4 in LNd_hash[k1][1][k3]:
                        if not k4 in temp3:
                            temp3.append(k4)
                temp3.sort()
                if not [k1.split('_')[1]]+temp2+[k1.split('_')[0]]+temp3 in temp1:
                    temp1.append([k1.split('_')[1]]+temp2+[k1.split('_')[0]]+temp3)
        tempall[k1]=temp1
        for k2 in tempall[k1]:
            if k2[0]==k1.split('_')[0]:
                for k3 in k2[1:]:
                    if k3 in k1.split('_'): break
                    if k3 in LNd_hash[k1][0].keys():
                        del LNd_hash[k1][0][k3]
            if k2[0]==k1.split('_')[1]:
                for k3 in k2[1:]:
                    if k3 in k1.split('_'): break
                    if k3 in LNd_hash[k1][1].keys():
                        del LNd_hash[k1][1][k3]
        temp2=[]
        clu1=clusterNums(sorted([int(i) for i in LNd_hash[k1][0].keys()]), 1000, 'f')
        cen1=clusterQC(clu1, 2)
        if not cen1==[]:
            sup1=clusterSupVis(sorted([int(i) for i in LNd_hash[k1][0].keys()]),cen1,1000)
            for k3 in sup1.keys():
                temp3=[k1.split('_')[0]]+[str(i) for i in sorted(sup1[k3])]+[k1.split('_')[1]]
                temp4=[]
                for k4 in sup1[k3]:
                    temp4+=LNd_hash[k1][0][str(k4)]
                temp4.sort()
                temp3+=temp4
                if not temp3 in temp2:
                    temp2.append(temp3)
        clu2=clusterNums(sorted([int(i) for i in LNd_hash[k1][1].keys()]), 1000, 'f')
        cen2=clusterQC(clu2, 2)
        if not cen2==[]:
            sup2=clusterSupVis(sorted([int(i) for i in LNd_hash[k1][1].keys()]),cen2,1000)
            for k3 in sup2.keys():
                temp3=[k1.split('_')[1]]+[str(i) for i in sorted(sup2[k3])]+[k1.split('_')[0]]
                temp4=[]
                for k4 in sup2[k3]:
                    temp4+=LNd_hash[k1][1][str(k4)]
                temp4.sort()
                temp3+=temp4
                if not temp3 in temp2:
                    temp2.append(temp3)
        tempall[k1]+=temp2
        for k2 in temp2:
            if k2[0]==k1.split('_')[0]:
                for k3 in k2[1:]:
                    if k3 in k1.split('_'): break
                    if k3 in LNd_hash[k1][0].keys():
                        del LNd_hash[k1][0][k3]
            if k2[0]==k1.split('_')[1]:
                for k3 in k2[1:]:
                    if k3 in k1.split('_'): break
                    if k3 in LNd_hash[k1][1].keys():
                        del LNd_hash[k1][1][k3]
    return [LNd_hash,tempall]

def LNe_hash_Filter(LNe_hash):
    out={}
    for k1 in LNe_hash.keys():
        tempall=[]
        for k2 in LNe_hash[k1]:
            khash={}
            kchr=k2[0]
            khash[kchr]=[]
            for k3 in k2[1:]:
                if k3 in allchromos or not type(k3) == type(1):
                    khash[k3]=[]
                    kchr=k3
                elif type(k3) == type(1) and int(k3)<23:
                    khash[k3]=[]
                    kchr=k3                    
                else:
                    if type(k3) == type(1) and int(k3)>22:
                        khash[kchr].append(int(k3))
            temp=[]
            if len(khash.keys())>1:
                for k3 in sorted(khash.keys()):
                    temp+=[k3]+sorted(khash[k3])
            if not temp in tempall and not temp==[]:
                tempall.append(temp)
        if not tempall==[]:
            out[k1]=tempall
    return out

def LNd_hash_pair(LNd_hash):
    out={}
    for k1 in LNd_hash.keys():
        if k1.split('_')[0] in SP_Info.keys() and k1.split('_')[1] in SP_Info.keys():
            out[k1]=[]
            cen1=[int(i) for i in LNd_hash[k1][0].keys()]
            if not cen1==[]:
                clu1=clusterSupVis(sorted(SP_Info[k1.split('_')[0]].keys()),cen1,10000)
                for k2 in cen1:
                    if k2 in clu1.keys():
                        temp=[]
                        for k3 in clu1[k2]:
                            temp.append(abs(k3-k2))
                        clu1[k2]=clu1[k2][temp.index(min(temp))]
                    else:
                        clu1[k2]=k2+100
            for k2 in LNd_hash[k1][0].keys():
                k3=sorted([int(i) for i in LNd_hash[k1][0][k2]])
                if len(k3)>1:
                    continue
                else:
                    if not k3==[]:
                        clu2=clusterSupVis(sorted(SP_Info[k1.split('_')[0]].keys()),k3,10000)
                        for k4 in k3:
                            if k4 in clu2.keys():
                                temp=[]
                                for k3 in clu2[k4]:
                                    temp.append(abs(k3-k4))
                                clu2[k4]=clu2[k4][temp.index(min(temp))]
                            else:
                                clu2[k4]=k4+100
                        LNd_hash[k1][0][k2]=[clu2.keys()[0],clu2[clu2.keys()[0]]]
            if not cen1==[]:
                for k2 in clu1.keys():
                    temp=[k1.split('_')[0]]+sorted([k2,clu1[k2]])+[k1.split('_')[1]]+sorted(LNd_hash[k1][0][str(k2)])
                    if not temp in out[k1]:
                        out[k1].append(temp)
    return out

def SP_Info_Cluster(SP_Info):
    for S_Chr in SP_Info.keys():
        bpCluster=clusterNums(sorted(SP_Info[S_Chr].keys()), 11, 'f')
        bpMilti=[]
        for k1 in range(len(bpCluster[1])):
            if bpCluster[1][k1]>1:
                bpMilti.append(bpCluster[0][k1])
        if not bpMilti==[]:
            bpClu2=clusterSupVis(sorted(SP_Info[S_Chr].keys()),sorted(bpMilti),11)
            for k1 in bpClu2.keys():
                temp=[]
                for k2 in bpClu2[k1]:
                    temp.append(SP_Info[S_Chr][k2])
                newk=bpClu2[k1][temp.index(max(temp))]
                for k2 in bpClu2[k1]:
                    if not k2==newk:
                        SP_Info[S_Chr][newk]+=SP_Info[S_Chr][k2]
                        del SP_Info[S_Chr][k2]
        else:
            bpClu2=[]
    del bpCluster
    del bpClu2

def SP_info_ReadIn(bps_hash,bps_in_path):
    single_chromos=bps_hash[S_Sample][S_Para].keys()
    SP_links={}
    SP_Link1={}
    SP_Info={}
    SP_Link2={}
    SP_Link3={}
    for S_Chr in single_chromos:
        SP_links[S_Chr]={}
        SP_Link1[S_Chr]=[]
        SP_Info[S_Chr]={}
        filein1=bps_hash[S_Sample][S_Para][S_Chr][0]
        fin=open(bps_in_path+filein1)
        for line in fin:
            pin=line.strip().split()
            if not int(pin[0]) in SP_Info[S_Chr].keys():
                SP_Info[S_Chr][int(pin[0])]=int(pin[1])
            else:
                SP_Info[S_Chr][int(pin[0])]=max(int(pin[1]),SP_Info[S_Chr][int(pin[0])])
            if not int(pin[2]) in SP_Info[S_Chr].keys():
                SP_Info[S_Chr][int(pin[2])]=int(pin[3])
            else:
                SP_Info[S_Chr][int(pin[2])]=max(int(pin[3]),SP_Info[S_Chr][int(pin[2])])
            if not pin[0]==pin[2]:
                SP_Link1[S_Chr]+=[pin[0],pin[2]]
                #SP_links[S_Chr][min([int(l) for l in [pin[0],pin[2]]])]=[]
            #SP_links[S_Chr][min([int(l) for l in [pin[0],pin[2]]])]+=[max([int(l) for l in [pin[0],pin[2]]])]
        fin.close()
        filein1=bps_hash[S_Sample][S_Para][S_Chr][1]
        fin=open(bps_in_path+filein1)
        for line in fin:
                pin=line.strip().split()
                if not int(pin[0]) in SP_Info[S_Chr].keys():
                    SP_Info[S_Chr][int(pin[0])]=int(pin[1])
                else:
                    SP_Info[S_Chr][int(pin[0])]=max(int(pin[1]),SP_Info[S_Chr][int(pin[0])])
        fin.close()
    SP_Info_Cluster(SP_Info)
    for S_Chr in SP_Link1.keys():
        SP_Link2[S_Chr]=[int(i) for i in SP_Link1[S_Chr]]
    for S_Chr in single_chromos:
        bpCluster=clusterSupVis(SP_Info[S_Chr].keys(),SP_Link2[S_Chr],10)
        SP_Link3[S_Chr]=[]
        for k1 in bpCluster.keys():
            if len(bpCluster[k1])==1 and bpCluster[k1][0]==k1:
                del bpCluster[k1]
        for k1 in range(len(SP_Link2[S_Chr])):
            k2=SP_Link2[S_Chr][k1]
            if k2 in bpCluster.keys():
                if len(bpCluster[k2])==1 and bpCluster[k2][0]==k2:
                    SP_Link3[S_Chr]+=[k2]
                else:
                    if len(bpCluster[k2])==1:
                        SP_Link3[S_Chr]+=[bpCluster[k2][0]]
                    else:
                        temp=[]
                        for k3 in bpCluster[k2]:
                            temp.append(SP_Info[S_Chr][k3])
                        k1_replace=bpCluster[k2][temp.index(max(temp))]
                        SP_Link3[S_Chr]+=[k1_replace]
            else:
                SP_Link3[S_Chr]+=[k2]
        SP_links[S_Chr]={}
        for k1 in range(len(SP_Link3[S_Chr])/2):
            minNum=min([SP_Link3[S_Chr][2*k1],SP_Link3[S_Chr][2*k1+1]])
            maxNum=max([SP_Link3[S_Chr][2*k1],SP_Link3[S_Chr][2*k1+1]])
            if not minNum in SP_links[S_Chr].keys():
                SP_links[S_Chr][minNum]=[]
            SP_links[S_Chr][minNum]+=[maxNum]
    return [SP_links,SP_Info,SP_Link3]

def SP_Info_Sort(SP_links):
    out={}
    for k1 in SP_links.keys():
        out[k1]=[]
        for k2 in sorted(SP_links[k1].keys()):
            if not sorted([k2]+SP_links[k1][k2]) in out[k1]:
                out[k1].append(sorted([k2]+SP_links[k1][k2]))
    return out

def BP_Pair_Collaps(SP_list_Chr,SP_Info_Chr):
    #eg:SP_list_Chr=SP_list['22'];SP_Info_Chr=SP_Info['22']
    out1={}
    out2=[]
    temp=[]
    out=[]
    Cff1=10**6
    for k1 in SP_list_Chr:
        k2=[]
        for k3 in k1:
            if not k3 in k2:
                k2.append(k3)
        k2.sort()
        if len(k2)>1:
            if len(k2)>4:
                out1[min(k2)]=sorted(k2)
            else:
                if k2[-1]-k2[0]<Cff1:
                    temp.append(k2)
                else:
                    out2.append(k2)
    if not temp==[]:
        temp2=temp[0]
        flag1=0
        while True:
            flag1+=1
            if flag1==len(temp): break
            if temp[flag1][0]-temp[flag1-1][-1]<500:
                temp2+=temp[flag1]
            else:
                out.append(temp2)
                temp2=temp[flag1]
    Cff2=5*(10**6)
    out3=[]
    if not out2==[]:
        for k1 in out2:
            k1.sort()
            k2=[[k1[0]]]
            for k3 in k1[1:]:
                if k3-k2[-1][-1]<Cff2:
                    k2[-1].append(k3)
                else:
                    k2.append([k3])
            k5=[]
            for k3 in k2:
                if len(k3)>1: 
                    k5+=k3
                else: 
                    k4=SearchClosest2(SP_Info_Chr.keys(),k3,Cff2)[0]
                    if not k4=={}:
                        k5+=sorted([k4.keys()[0],k4[k4.keys()[0]]])
                    else:
                        k5+=k3
            out3.append(k5)
    if not out+out3==[]:
        for k1 in out+out3:
            if not min(k1) in out1.keys():
                out1[min(k1)]=[]
            out1[min(k1)]+=k1
    reorder=[]
    for k1 in sorted(out1.keys()):
        reorder.append(out1[k1])
    return reorder

def Left_SP_Cluster(Left_SP_1,Cff):
    Left_SP=sorted(Left_SP_1)
    out=[[Left_SP[0]]]
    for k1 in Left_SP[1:]:
        if k1-out[-1][-1]<Cff:
            out[-1].append(k1)
        else:
            out.append([k1])
    out1=[]
    out2=[]
    for k1 in out:
        if len(k1)>1:
            out1.append(k1)
        else:
            out2+=k1
    return [out1,out2]

def pair_link(SP_LI_in):
    singletons={}
    SP_links=SP_LI_in[0]
    SP_Info=SP_LI_in[1]
    SP_Link1=SP_LI_in[2]
    SP_Single={}
    SP_link2={}
    SP_link4={}
    for k1 in SP_Info.keys():
        SP_Single[k1]=[]
        for k2 in SP_Info[k1].keys():
            if SP_Info[k1][k2]<2:
                del SP_Info[k1][k2]
            else:
                if not k2 in SP_Link1[k1]:
                    SP_Single[k1].append(k2)
        SP_assign=SearchClosest2(SP_Link1[k1], SP_Single[k1],5000)
        Left_SP=SP_assign[1]
        Paired_SP=SP_assign[0]
        if not Left_SP==[]:
            SingletonPair=Left_SP_Cluster(Left_SP,5000)
            singletons[k1]=SingletonPair[1]
            for k2 in Paired_SP.keys():
                if SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2] in SP_links[k1].keys():
                    tempk=SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2]
                elif SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2+1] in SP_links[k1].keys():
                    tempk=SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2+1]
                else: continue
                SP_links[k1][tempk].append(k2)
            for k2 in SingletonPair[0]:
                if not min(k2) in SP_links[k1].keys():
                    SP_links[k1][min(k2)]=[]
                SP_links[k1][min(k2)]+=k2
    SP_list=SP_Info_Sort(SP_links)
    collapsed_SP_list={}
    for k1 in SP_list.keys():
        collapsed_SP_list[k1]=BP_Pair_Collaps(SP_list[k1],SP_Info[k1])   
    return [collapsed_SP_list,singletons]

def single_split(singleton_chr,out_single):
    if len(singleton_chr)<5:
        out_single.append(singleton_chr)
    else:
        single_inter=[]
        for kx in range(len(singleton_chr[1:])):
            single_inter.append(singleton_chr[kx+1]-singleton_chr[kx])
        max_index=1
        max_rec=single_inter[1]
        for kx in range(2,len(single_inter)-1):
            if single_inter[kx]>max_rec:
                max_rec=single_inter[kx]
                max_index=kx
        temp=[singleton_chr[:(max_index+1)],singleton_chr[(max_index+1):]]
        for temp_list in temp:
            single_split(temp_list,out_single)

def singleton_split(singleton_chr):
    #singleton_chr=singletons['22']
    out_single=[]
    SingletonChr=sorted(singleton_chr)
    single_split(SingletonChr,out_single)
    return out_single

def write_bp_1(SP_link3,missed_pairs):
    fout_rec=[]
    max_num_rec=[[]]
    max_num=int(dict_opts['--batch'])
    for k1 in SP_link3.keys():
        for k2 in SP_link3[k1]:
            k2.sort()
            k3=[k2[0]]
            for k4 in k2:
                if not k4 in k3 and k4-k3[-1]>10:
                    k3.append(k4)
            k5=k3
            if len(k5)>1 and not len(k5)>8:
                if len(max_num_rec[-1])<max_num:
                    max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])
                elif len(k5)==1:
                    max_num_rec.append([])
                    max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])                                                
            else:
                singletons[k1]+=k5
        for k5 in missed_pairs[k1]:
            if len(max_num_rec[-1])<max_num:
                    max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])
            else:
                max_num_rec.append([])
                max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])                                                
    for k1 in singletons.keys():      
        if len(singletons[k1])>1:
            out_single=singleton_split(singletons[k1])
            for k2 in out_single:
                if len(max_num_rec[-1])<max_num:
                    max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k2])
                else:
                    max_num_rec.append([])
                    max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k2])
    rec=0
    for k1 in max_num_rec:
        rec+=1
        fout=para_sample_bps_fold+S_Sample+'.'+str(rec)+'.'+'txt'
        if not fout in fout_rec:
            fout_rec.append(fout)
        if not os.path.isfile(fout):
            fo=open(fout,'w')
        else:
            fo=open(fout,'a')
        for k2 in k1:
            print >>fo, ' '.join(k2)
            print >>fo, ' '
        fo.close()
    for fout in fout_rec:
        fin=open(fout)
        test=[]
        for line in fin:
            pin=line.strip().split()
            if not pin==[]:
                if not pin in test:
                    test.append(pin)
        fin.close()
        fo=open(fout,'w')
        for k1 in test:
            print >>fo, ' '.join(k1)
            print >>fo, ' '
        fo.close()

def bp_lists_length_decide(bps):
    #eg:bps=[[41577777, 41585823, 41601593, 41602950, 41604425, 41605770, 41613170, 41637479, 41656018, 41657463, 41691589]]
    flag=0
    for x in bps:
        if len(x)>8:
            flag+=1
    if flag==0:
        return 'TRUE'
    else:
        return 'FALSE'

def bp_list_separate(bps):
    #eg:bps=[[41577777, 41585823, 41601593, 41602950, 41604425, 41605770, 41613170, 41637479, 41656018, 41657463, 41691589]]
    while True:
        if bp_lists_length_decide(bps)=='TRUE':
            break
        elif bp_lists_length_decide(bps)=='FALSE':
            for bps2 in bps:
                if len(bps2)>8:
                    inter=[bps2[x+1]-bps2[x] for x in range(len(bps2)-1)]
                    inter1=max(inter[1:-1])
                    bps+=[bps2[:inter.index(inter1)+1],bps2[inter.index(inter1)+1:]]
                    del bps[bps.index(bps2)]

def bps_list_chop(SP_link4):
    for x in SP_link4.keys():
        for y in SP_link4[x]:
            if len(y)>8:
                #print y
                bps_y=[y]
                bp_list_separate(bps_y)
                SP_link4[x]+=bps_y

def write_bp_2(SP_link3,missed_pairs):
    #if not '--batch' in dict_opts.keys():
    fout_rec=[]
    for k1 in SP_link3.keys():
        fout=para_sample_bps_fold+S_Sample+'.'+str(k1)+'.'+'txt'
        if not fout in fout_rec:
            fout_rec.append(fout)
        if not os.path.isfile(fout):
            fo=open(fout,'w')
        else:
            fo=open(fout,'a')
        for k2 in SP_link3[k1]:
            k2.sort()
            k3=[k2[0]]
            for k4 in k2:
                if not k4 in k3 and k4-k3[-1]>10:
                    k3.append(k4)
            k5=k3
            if len(k5)>1 and not len(k5)>8:
                print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k5])
                print >>fo, ' '
            elif len(k5)==1:
                if k1 in singletons.keys():
                    singletons[k1]+=k5
                else:
                    singletons[k1]=k5
        for k2 in missed_pairs[k1]:
            k2.sort()
            print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
            print >>fo, ' '
        fo.close()
    for k1 in singletons.keys():
        fout=para_sample_bps_fold+S_Sample+'.'+str(k1)+'.'+'txt'
        if not fout in fout_rec:
            fout_rec.append(fout)
        if not os.path.isfile(fout):
            fo=open(fout,'w')
        else:
            fo=open(fout,'a')
        if len(singletons[k1])>1:
            out_single=singleton_split(singletons[k1])
            for k2 in out_single:
                print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
                print >>fo, ' '
        fo.close()
    for fout in fout_rec:
        fin=open(fout)
        test=[]
        for line in fin:
            pin=line.strip().split()
            if not pin==[]:
                if not pin in test:
                    test.append(pin)
        fin.close()
        fo=open(fout,'w')
        for k1 in test:
            print >>fo, ' '.join(k1)
            print >>fo, ' '
        fo.close()

def write_bp_3(SP_link3,missed_pairs):
    fout_rec=[]
    fout=para_sample_bps_fold+S_Sample+'.txt'
    if not fout in fout_rec:
        fout_rec.append(fout)        
    if not os.path.isfile(fout):
        fo=open(fout,'w')
    else:
        fo=open(fout,'a')
    fo.close()
    for k1 in SP_link3.keys():
        fo=open(fout,'a')
        for k2 in SP_link3[k1]:
            k2.sort()
            k3=[k2[0]]
            for k4 in k2:
                if not k4 in k3 and k4-k3[-1]>10:
                    k3.append(k4)
            k5=k3
            if len(k5)>1 and not len(k5)>8:
                print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k5])
                print >>fo, ' '
            elif len(k5)==1:
                if k1 in singletons.keys():
                    singletons[k1]+=k5
                else:
                    singletons[k1]=k5
        for k2 in missed_pairs[k1]:
            k2.sort()
            print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
            print >>fo, ' '
        fo.close()
    for k1 in singletons.keys():
        fo=open(fout,'a')
        if len(singletons[k1])>1:
            out_single=singleton_split(singletons[k1])
            for k2 in out_single:
                print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
                print >>fo, ' '
        fo.close()
    for fout in fout_rec:
        fin=open(fout)
        test=[]
        for line in fin:
            pin=line.strip().split()
            if not pin==[]:
                if not pin in test:
                    test.append(pin)
        fin.close()
        fo=open(fout,'w')
        for k1 in test:
            print >>fo, ' '.join(k1)
            print >>fo, ' '
        fo.close()

def missed_pair_check(SP_Check,SP_link3):
    out={}
    out2={}
    for k1 in SP_Check.keys():
        rec1=0
        rec2=0
        temp={}
        temp1=SP_Check[k1]
        temp1.sort()
        temp2=SP_link3[k1]
        for x in temp2:
            if not x[0] in temp.keys():
                temp[x[0]]=[]
            temp[x[0]].append(x)
        temp3=[]
        for x in sorted(temp.keys()):
            temp3+=temp[x]
        temp4={}
        while True:
            if rec1==len(temp1) or rec2==len(temp3): break
            if temp3[rec2][-1]-temp3[rec2][0]<10**7:
                if temp1[rec1]<temp3[rec2][0]-10**4:
                    rec1+=1
                elif temp1[rec1]>temp3[rec2][1]+10**4:
                    rec2+=1
                else:
                    if not '_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])]) in temp4.keys():
                        temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])]=[]
                    if not temp3[rec2] in temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])]:
                        temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])].append(temp3[rec2])
                    if not temp1[rec1] in temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])][0]:
                        temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])].append(temp1[rec1])
                    rec1+=1
            else:
                rec2+=1
        out[k1]=[]
        temp5=[]
        for x in temp4.keys():
            if len(temp4[x])==1:
                out[k1].append(sorted(temp4[x][0]))
                temp5+=sorted(temp4[x][0])
            else:
                out[k1].append(sorted(temp4[x][0]+temp4[x][1:]))
                temp5+=sorted(temp4[x][0]+temp4[x][1:])
        temp6=[]
        for x in temp1:
            if not x in temp5:
                temp6.append(x)
        temp6.sort()
        temp6b=[]
        for x in temp3:
            flagb=0
            for y in x:
                if not y in temp5:
                    flagb+=1
            if not flagb==0:
                temp6b.append(x)
        out[k1]+=temp6b
        data_cen=clusterNums(temp6, 10**6, 'r')[0]
        out2[k1]=[]
        if not data_cen==[]:
            temp7=[[data_cen[0]]]
            rec1=0
            for x in temp6:
                if x<temp7[-1][0]:
                    temp7[-1].append(x)
                else:
                    rec1+=1
                    if rec1<len(data_cen):
                        temp7.append([data_cen[rec1]])  
            for x in temp7:
                if len(x)>1:
                    out2[k1].append(x)
    return [out,out2]

min_length=100
def link_SP_link3(SP_link3):
    out={}
    for k1 in SP_link3.keys():
        temp2=[]
        temp1={}
        for k2 in SP_link3[k1]:
            if not k2[0] in temp1.keys():
                temp1[k2[0]]=[]
            temp1[k2[0]].append(k2)
        for k2 in sorted(temp1.keys()):
            for k3 in temp1[k2]:
                temp2.append(k3)
        out[k1]=[]
        for k2 in temp2:
            if out[k1]==[]:
                out[k1].append(k2)
            else:
                if k2[0]<out[k1][-1][-1]+10 and out[k1][-1][-1]-out[k1][-1][0]<2*10**6:
                    out[k1][-1]+=k2
                else:
                    out[k1].append(k2)
    out2={}
    for k1 in out.keys():
        out2[k1]=[]
        for k2 in out[k1]:
            k2.sort()
            temp=[k2[0]]
            for k3 in k2[1:]:
                if not k3 in temp and k3-temp[-1]>100:
                    temp.append(k3)
            if len(temp)>1:
                out2[k1].append(temp)
    for k1 in out2.keys():
        temp=[]
        for k3 in out2[k1]:
            temp2=[k3[0]]
            for k4 in k3[1:]:
                if k4-temp2[-1]>min_length:
                    temp2.append(k4)
                else:
                    if SP_Info[k1][k4]>SP_Info[k1][temp2[-1]]:
                        temp2[-1]=k4
                    else:
                        continue
            if len(temp2)>1:
                temp.append(temp2)
        out2[k1]=temp
    return out2

opts,args=getopt.getopt(sys.argv[1:],'o:h:s:',['NullModel=','batch=','help=','ppre=','PathBam=','sample=','chr=','SamplePbs=','ref=','NullSplitLength=','NullILCI=','NullRDCI=','NullTBCI=','NullILCff=','NullSPCff=','NullDRCff=','NullBed=','NullBedLength=','NullSampleLength=','NullSampleNumber=','NullSamplePercentage=','ExcludeBed=','IncludeBed=','ToolMappingQ=','FileMappingQ=','SplitLength=','BPSPCff=','BPLNCff=','BPalignQC=','BPalignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)
if not '--ReadLen' in dict_opts.keys():
    ReadLength=101
else:
    ReadLength=int(dict_opts['--ReadLen'])

if '--ToolMappingQ' in dict_opts.keys() and '--FileMappingQ' in dict_opts.keys():
    ToolMappingQ=dict_opts['--ToolMappingQ']
    FileMappingQ=dict_opts['--FileMappingQ']
    align_QCflag=1
else:
    align_QCflag=0

if '--BPalignQC' in dict_opts.keys():
    BPalignQC=float(dict_opts['--BPalignQC'])
else:
    BPalignQC=0.2

if '--QCAlign' in dict_opts.keys():
    QCAlign=int(dict_opts['--QCAlign'])
else:
    QCAlign=20

if '--QCSplit' in dict_opts.keys():
    QCSplit=int(dict_opts['--QCSplit'])
else:
    QCSplit=20

if '--NullSplitLength' in dict_opts.keys():
    Null_SplitLen_perc=float(dict_opts['--NullSplitLength'])
else:
    Null_SplitLen_perc=0.1

if '--BPalignQCFlank' in dict_opts.keys():
    BPalignQCFlank=int(dict_opts['--BPalignQCFlank'])
else:
    BPalignQCFlank=500

para_filter=[]
if '--BPSPCff' in dict_opts.keys() and '--BPLNCff' in dict_opts.keys() and '--BPalignQC' in dict_opts.keys():
    para_filter=['SPCff'+dict_opts['--BPSPCff']+'.CluCff'+dict_opts['--BPLNCff']+'.AlignCff'+dict_opts['--BPalignQC']]

#affix1='.'.join(['SPCff'+dict_opts['--SPCff'],'CluCff'+dict_opts['--CluCff'],'AlignCff'+dict_opts['--AlignQC'],'DUCff'+dict_opts['--DUCff'],'LNs'])
#affix2='.'.join(['SPCff'+dict_opts['--SPCff'],'CluCff'+dict_opts['--CluCff'],'AlignCff'+dict_opts['--AlignQC'],'DUCff'+dict_opts['--DUCff'],'SPs'])
#affix3='.'.join(['SPCff'+dict_opts['--SPCff'],'CluCff'+dict_opts['--CluCff'],'AlignCff'+dict_opts['--AlignQC'],'DUCff'+dict_opts['--DUCff'],'chromLNs'])
if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
    print 'SVelter-0.1          Last Update:2014-08-20'
    print ' '
    print 'Required Parameters:'
    print '--ppre, workind directory of SVelter, eg: .../SVelter/' 
    print '--ref, absolute path of reference genome. eg: .../SVelter/reference/genome.fa'
    print '-s, absolute path of input bam file. eg: .../SVelter/BamFiles/sample.bam'
    #print '--PathBam, absolute path of folder containing all input bam files. eg: .../SVelter/BamFiles/'
    print ' '
    print 'Optional Parameters:'
    print '--batch, specify number of bp clusters in each file; if 0 prvided,output files will be classified by chromosomes; if not --batch provided, all bp clusters will be writen into one txt file;'
else:
    if not '--ppre' in dict_opts.keys():
        print 'Error: please specify working directory using: --ppre'
    else:
        workdir=dict_opts['--ppre']
        if not workdir[-1]=='/':
            workdir+='/'
        if not '-o' in dict_opts.keys():
            out_path=workdir
        else:
            out_path=dict_opts['-o']
        if not out_path[-1]=='/':
            out_path+='/'
        bps_in_path=workdir+'BreakPoints.'+dict_opts['-s'].split('/')[-1]+'/'
        if not '-s' in dict_opts.keys() and not '--PathBam' in dict_opts.keys():
            print 'Error: please specify either input file using -s or input path using --PathBam'
        else:
            if '-s' in dict_opts.keys():
                bam_path='/'.join(dict_opts['-s'].split('/')[:-1])+'/'
                bam_files=[dict_opts['-s']]
                bam_names=[dict_opts['-s'].split('/')[-1].replace('.bam','')]
            else:
                print 'Error: please specify input bam file using -s'
                bam_path=dict_opts['--PathBam']
                if not bam_path[-1]=='/':
                    bam_path+='/'
                bam_files=[]
                bam_names=[]
                for file in os.listdir(bam_path):
                    if file.split('.')[-1]=='bam':
                        bam_files.append(bam_path+file)
                        bam_names.append(file.replace('.bam',''))
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
                    fref=open(ref_index)
                    for line in fref:
                        pref=line.strip().split()
                        chromos.append(pref[0])
                    fref.close()
                    allchromos=chromos
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
                    bps_hash={}
                    for i in bam_names:
                        bps_hash[i]={}
                    bps_folder=out_path+'bp_files.'+dict_opts['-s'].split('/')[-1]+'/'
                    if not os.path.isdir(bps_folder):
                        os.system(r'''mkdir %s'''%(bps_folder))
                    for i in bam_names:
                        bps_hash[i]={}
                        for file1 in os.listdir(bps_in_path):
                            if file1.split('.')[-1]=='LNs':
                                if i in file1: 
                                    keyj='.'.join(file1.split('.')[-5:-1])
                                    single_chr_name=file1.split('.')[-6]
                                    if not keyj in bps_hash[i].keys():
                                        bps_hash[i][keyj]={}
                                    if single_chr_name in chromos:
                                        bps_hash[i][keyj][single_chr_name]=[]
                                        bps_hash[i][keyj][single_chr_name].append(file1)
                                        bps_hash[i][keyj][single_chr_name].append(file1.replace('LNs','SPs'))
                    for i in bps_hash.keys():
                        if not os.path.isdir(bps_folder+i):
                            os.system(r'''mkdir %s'''%(bps_folder+i))
                        for j in bps_hash[i].keys():
                            if not os.path.isdir(bps_folder+i+'/'+j):
                                os.system(r'''mkdir %s'''%(bps_folder+i+'/'+j))
                    for S_Sample in bps_hash.keys():
                        if para_filter==[]:
                            para_filter=sorted(bps_hash[S_Sample].keys())[::-1]
                        for S_Para in para_filter:
                            rec=0
                            para_sample_bps_fold=bps_folder+S_Sample+'/'+S_Para+'/'
                            SP_LI_in=SP_info_ReadIn(bps_hash,bps_in_path)
                            SP_links=SP_LI_in[0]
                            SP_Info=SP_LI_in[1]
                            SP_Check=SP_LI_in[2]
                            SP_link3a=pair_link(SP_LI_in)
                            singletons=SP_link3a[1]
                            SP_link3=SP_link3a[0]
                            missed_pairs=missed_pair_check(SP_Check,SP_link3)
                            SP_link4=link_SP_link3(missed_pairs[0])
                            bps_list_chop(SP_link4)
                            missed_pairs=link_SP_link3(missed_pairs[1])
                            #pout=bps_in_path+'BPfiles/'
                            #if not os.path.isdir(pout):
                            #   os.system(r'''mkdir %s'''%(pout))
                            devide_start=10**6
                            if not '--batch' in dict_opts.keys():
                                write_bp_3(SP_link4,missed_pairs)
                            else:
                                if dict_opts['--batch']=='0':
                                    write_bp_2(SP_link4,missed_pairs)
                                else:
                                    write_bp_1(SP_link4,missed_pairs)
                            S_Chr=bps_hash[S_Sample][S_Para].keys()[0] 
                            chrom_LN=bps_hash[S_Sample][S_Para][S_Chr][1].replace('.'+S_Chr+'.SPCff','.SPCff').replace('SPs','chromLNs')
                            LNall_hash=LNd_hash_readin(chrom_LN)
                            LNd_hash=LNall_hash[0]
                            LNe_hash=LNe_hash_Filter(LNall_hash[1])
                            LNd_hash=LNd_hash_pair(LNd_hash)
                            LN_out=[]
                            for k1 in LNe_hash.keys():
                                for k2 in LNe_hash[k1]:
                                    temp=[]
                                    for k3 in k2:
                                        if k3 in chromos:
                                            temp.append([k3])
                                        else:
                                            if not temp==[]:
                                                temp[-1].append(str(k3))
                                    LN_out.append(temp)
                            for k1 in LNd_hash.keys():
                                for k2 in LNd_hash[k1]:
                                    temp=[]
                                    for k3 in k2:
                                        if k3 in chromos:
                                            temp.append([k3])
                                        else:
                                            if not temp==[]:
                                                temp[-1].append(str(k3))
                                    LN_out.append(temp)
                            if not '--batch' in dict_opts.keys():
                                fout=para_sample_bps_fold+S_Sample+'.txt'
                                if not os.path.isfile(fout):
                                    fo=open(fout,'w')
                                else:
                                    fo=open(fout,'a')
                                for k1 in LN_out:
                                    for k2 in k1:
                                        print >>fo, ' '.join(k2)
                                    print >>fo, ' '
                                fo.close()
                            else:
                                if dict_opts['--batch']=='0':
                                    fout=para_sample_bps_fold+S_Sample+'.'+'LN'+'.'+'txt'
                                    if not os.path.isfile(fout):
                                        fo=open(fout,'w')
                                    else:
                                        fo=open(fout,'a')
                                    for k1 in LN_out:
                                        for k2 in k1:
                                            print >>fo, ' '.join(k2)
                                        print >>fo, ' '
                                    fo.close()
                                else:
                                    max_num=int(dict_opts['--batch'])
                                    LN_out2=[[]]
                                    for k1 in LN_out:
                                        if len(LN_out2[-1])<max_num:
                                            LN_out2[-1].append(k1)
                                        else:
                                            LN_out2.append([])
                                            LN_out2[-1].append(k1)
                                    for k1 in LN_out2:
                                        rec+=1
                                        fout=para_sample_bps_fold+S_Sample+'.'+str(rec)+'.'+'txt'
                                        if not os.path.isfile(fout):
                                            fo=open(fout,'w')
                                        else:
                                            fo=open(fout,'a')
                                        for k2 in k1:
                                            for k3 in k2:
                                                print >>fo, ' '.join(k3)
                                            print >>fo, ' '
                                        fo.close()
