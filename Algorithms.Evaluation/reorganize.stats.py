#command='python reorganize.stats.py --ppre /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.comp/Delly/'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import datetime
opts,args=getopt.getopt(sys.argv[1:],'i:o:h:',['ppre=','in='])
dict_opts=dict(opts)
ppre=dict_opts['--ppre']
if not ppre[-1]=='/':
	ppre+='/'

fo=open(ppre+'Pseudo.ROC.allSV.stats','w')
for k1 in os.listdir(ppre):
	if k1.split('.')[-1]=='stats' and k1.split('.')[-2]=='BPlink' and not k1=='Produce.Pseudo.ROC.stats':
		sample_name=k1.replace('.stats','')
		print sample_name
		fin=open(ppre+k1)
		temp=[0,0,0,0]
		rec=0
		for line in fin:
			pin=line.strip().split()
			if not pin[1]=='Y':
				temp[0]+=float(pin[2])
				temp[1]+=float(pin[3])
				temp[2]+=float(pin[4])
				temp[3]+=float(pin[5])
				rec+=1
		temp=[i/float(rec) for i in temp]
		print >>fo,' '.join([sample_name]+[str(i) for i in temp])

fo.close()


