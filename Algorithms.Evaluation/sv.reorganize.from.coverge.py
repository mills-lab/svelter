#command='python sv.reorganize.from.coverge.py --ppre /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/bp_files/NA12878_S1/SPCff10.CluCff15.AlignCff0.0/ --ref /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/reference/genome.fa'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import time
import datetime
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:o:h:',['ppre=','ref='])
dict_opts=dict(opts)
ref=dict_opts['--ref']
path2=dict_opts['--ppre']
if not path2[-1]=='/':
	path2+='/'

chromos=[]
ref_index=ref+'.fai'
fin=open(ref_index)
for line in fin:
	pin=line.strip().split()
	chromos.append(pin[0])

fin.close()

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


sv_info={}
for k3 in os.listdir(path2):
	if k3.split('.')[-1]=='coverge':
		read_in_structures(path2+k3)

fout=path2+'sv.reorganized.rec'
fo=open(fout,'w')
for k1 in sv_info.keys():
	for k2 in sv_info[k1].keys():
		if not k2==k1:
			for k3 in sv_info[k1][k2]:
				print >>fo, ' '.join([k1,k2]+[str(i) for i in k3])

fo.close()







