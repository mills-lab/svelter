import os
def path_modify(path):
	if not path[-1]=='/':
		path+='/'
	return path

def Vali_files_readin(Vali_Folder):
	out=[]
	Vali_Folder=path_modify(Vali_Folder)
	for k1 in os.listdir(Vali_Folder):
		if k1.split('.')[-1]=='PacVal':
			out.append(Vali_Folder+k1)
	return out

def Vali_file_read(filein,Vali_hash,score_cff):
	fin=open(filein)
	pin=fin.readline().strip().split()
	for line in fin:
		pin=line.strip().split()
		if not pin[4] in Vali_hash.keys():
			Vali_hash[pin[4]]={}
		if not pin[5] in Vali_hash[pin[4]].keys():
			Vali_hash[pin[4]][pin[5]]=[0,0]
		if float(pin[-1])>score_cff:
			Vali_hash[pin[4]][pin[5]][0]+=1
		Vali_hash[pin[4]][pin[5]][1]+=1
	fin.close()
	return Vali_hash

ppre='/scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/VaLoR_Vali/'
fin=open(ppre+'NA12878.Subtype.SV')
SV_hash={}
for line in fin:
	pin=line.strip().split()
	if not pin==[]:
		if not pin[0] in SV_hash.keys():
			SV_hash[pin[0]]={}
		if not pin[1] in SV_hash[pin[0]].keys():
			SV_hash[pin[0]][pin[1]]={}
		if not pin[2] in SV_hash[pin[0]][pin[1]].keys():
			SV_hash[pin[0]][pin[1]][pin[2]]=[]

score_cff=0.1
Vali_Folder=ppre
Vali_files=Vali_files_readin(Vali_Folder)
Vali_hash={}
for k1 in Vali_files:
	Vali_hash=Vali_file_read(k1,Vali_hash,score_cff)

for k1 in SV_hash.keys():
	for k2 in SV_hash[k1].keys():
		for k3 in SV_hash[k1][k2].keys():
			SV_hash[k1][k2][k3]=Vali_hash[k2][k3]
			print ' '.join([str(i) for i in [k1,k2,k3]+SV_hash[k1][k2][k3]])

SVs=['INV_DUP','INV_DEL','DEL_DUP','DEL_DUP_INV']
Figure_hash={}
for k1 in SVs:
	Figure_hash[k1]=[0,0]
	for k2 in SV_hash[k1].keys():
		for k3 in SV_hash[k1][k2].keys():
			print [k1,k2,k3]
			Figure_hash[k1][0]+=SV_hash[k1][k2][k3][0]
			Figure_hash[k1][1]+=SV_hash[k1][k2][k3][1]

#print validated and all SV numbers for the four SVs listed above
for k1 in Figure_hash.keys():
	print ' '.join([str(i) for i in [k1]+Figure_hash[k1]])

Figure_hash={}
for k1 in SV_hash.keys():
	Figure_hash[k1]=[0,0]
	for k2 in SV_hash[k1].keys():
		for k3 in SV_hash[k1][k2].keys():
			Figure_hash[k1][0]+=SV_hash[k1][k2][k3][0]
			Figure_hash[k1][1]+=SV_hash[k1][k2][k3][1]

listed_SVs=[]
for k1 in SV_hash.keys():
	for k2 in SV_hash[k1].keys():
		for k3 in SV_hash[k1][k2].keys():
			listed_SVs.append([k2,k3])

for k1 in Vali_hash:
	for k2 in Vali_hash[k1].keys():
		if not [k1,k2] in listed_SVs:
			print ' '.join([str(i) for i in [k1,k2]+Vali_hash[k1][k2]])


for k1 in Figure_hash.keys():
	print ' '.join([str(i) for i in [k1]+Figure_hash[k1]+[float(Figure_hash[k1][0])/float(Figure_hash[k1][1])]])





ppre='/scratch/remills_flux/xuefzhao/CHM1/IL500/hg19/SVelter.rec4/'
fin=open(ppre+'CHM1.Subtype.SV')
SV_hash={}
for line in fin:
	pin=line.strip().split()
	if not pin==[]:
		if not pin[0] in SV_hash.keys():
			SV_hash[pin[0]]={}
		if not pin[1] in SV_hash[pin[0]].keys():
			SV_hash[pin[0]][pin[1]]={}
		if not pin[2] in SV_hash[pin[0]][pin[1]].keys():
			SV_hash[pin[0]][pin[1]][pin[2]]=[]

fin.close()

score_cff=0.2
Vali_Folder=ppre
Vali_files=Vali_files_readin(Vali_Folder)
Vali_hash={}
for k1 in Vali_files:
	Vali_hash=Vali_file_read(k1,Vali_hash,score_cff)

for k1 in SV_hash.keys():
	for k2 in SV_hash[k1].keys():
		for k3 in SV_hash[k1][k2].keys():
			if k3 in Vali_hash[k2].keys():
				SV_hash[k1][k2][k3]=Vali_hash[k2][k3]
				print ' '.join([str(i) for i in [k1,k2,k3]+SV_hash[k1][k2][k3]])

SVs=['DUP_INV','DEL_INV','DEL_DUP','DEL_DUP_INV']
Figure_hash={}
for k1 in SVs:
	if k1 in SV_hash.keys():
		Figure_hash[k1]=[0,0]
		for k2 in SV_hash[k1].keys():
			for k3 in SV_hash[k1][k2].keys():
				print [k1,k2,k3]
				if len(SV_hash[k1][k2][k3])>0:
					Figure_hash[k1][0]+=SV_hash[k1][k2][k3][0]
					Figure_hash[k1][1]+=SV_hash[k1][k2][k3][1]

for k1 in Figure_hash.keys():
	print ' '.join([str(i) for i in [k1]+Figure_hash[k1]])


Figure_hash={}
for k1 in SV_hash.keys():
	if k1 in SV_hash.keys():
		Figure_hash[k1]=[0,0]
		for k2 in SV_hash[k1].keys():
			for k3 in SV_hash[k1][k2].keys():
				print [k1,k2,k3]
				if len(SV_hash[k1][k2][k3])>0:
					Figure_hash[k1][0]+=SV_hash[k1][k2][k3][0]
					Figure_hash[k1][1]+=SV_hash[k1][k2][k3][1]

for k1 in Figure_hash.keys():
	print ' '.join([str(i) for i in [k1]+Figure_hash[k1]])




listed_SVs=[]
for k1 in SV_hash.keys():
	for k2 in SV_hash[k1].keys():
		for k3 in SV_hash[k1][k2].keys():
			listed_SVs.append([k2,k3])

for k1 in Vali_hash:
	for k2 in Vali_hash[k1].keys():
		if not [k1,k2] in listed_SVs:
			print ' '.join([str(i) for i in [k1,k2]+Vali_hash[k1][k2]])


fo=open('/scratch/remills_flux/xuefzhao/CHM1/IL500/hg19/SVelter/VaLoR_Vali/Integrated.Validated.Result.Cff0.2','w')
for k1 in sorted(Vali_hash.keys()):
	for k2 in sorted(Vali_hash[k1].keys()):
		print >>fo,' '.join([str(i) for i in [k1,k2]+Vali_hash[k1][k2]+[float(Vali_hash[k1][k2][0])/float(Vali_hash[k1][k2][1])]])

fo.close()





