#command='python vcf.BPlink.compare.py -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/Delly/Delly_comp_RD10.BPLink.vcf'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import glob
import datetime
opts,args=getopt.getopt(sys.argv[1:],'i:o:h:',['ref=','in='])
dict_opts=dict(opts)

vcf_in=dict_opts['-i']
fin=open(vcf_in)
fo=open(vcf_in.replace('.vcf','.PEcontrol.vcf'),'w')
cff=10
for line in fin:
	pin=line.strip().split()
	if not pin[0][0]=='#':
		qual=0
		for x in pin[7].split(';'):
			if x.split('=')[0]=='PE':
				qual=int(x.split('=')[1])
		if qual>cff:
			print >>fo, '\t'.join(pin)
	else:
			print >>fo, '\t'.join(pin)

fin.close()
fo.close()		



