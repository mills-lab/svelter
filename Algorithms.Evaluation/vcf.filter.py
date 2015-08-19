#vcf filter
#command='python vcf.filter.py --in /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/SVelter/SVelter_comp_RD60.RL.101.BPLink.vcf --cff -10'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import datetime
opts,args=getopt.getopt(sys.argv[1:],'i:o:h:',['cff=','in='])
dict_opts=dict(opts)
cutoff=float(dict_opts['--cff'])
filein=dict_opts['--in']
fin=open(filein)
fo=open(filein.replace('.vcf','.cff'+str(cutoff)+'.vcf'),'w')
for line in fin:
	pin=line.strip().split()
	if pin[0][0]=='#':
		print >>fo,'\t'.join(pin)
	else:
		if float(pin[5])>cutoff:
			print >>fo, '\t'.join(pin)

fin.close()
fo.close()


