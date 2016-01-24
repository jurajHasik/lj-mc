#!/usr/bin/python3
import sys
import math

print('Opening'+sys.argv[1])
f = open(sys.argv[1],'r') #First arg is filename
eql = 0
prl = 0

eAvg = 0
pAvg = 0
Elst = []
Plst = []
for line in f:
	fields = line.split()
	if fields[0]=='rT':
		print('rT: '+fields[1])
	if fields[0]=='rDens':
		print('rDens: '+fields[1])
	if fields[0]=='E_lrc:':
		print(line)
	if fields[0]=='EQ':
		eql=eql+1
	if fields[0]=='PR':
		prl = prl+1
		eAvg = eAvg + float(fields[4])
		Elst.append(float(fields[4]))
		pAvg = pAvg + float(fields[6]) 		
		Plst.append(float(fields[6]))
print('Num of samples: '+str(prl))
eAvg = eAvg / prl
pAvg = pAvg / prl
eVar = 0
pVar = 0
for i in range(prl):
	eVar = eVar + (Elst[i]-eAvg)*(Elst[i]-eAvg)
	pVar = pVar + (Plst[i]-pAvg)*(Plst[i]-pAvg)
eVar = eVar/prl
pVar = pVar/prl
print('E: {0:.5f} +/- {1:.5f}'.format(eAvg,math.sqrt(eVar)))
print('P: {0:.5f} +/- {1:.5f}'.format(pAvg,math.sqrt(pVar)))