'''
Created by Patrick Kilecdi
To sort and examine the run through logs
'''

f=open("restart_relaxations.dat","r")
lines=[]
st=f.readline()
while st!="":
	lines.append(st)
	st=f.readline()
f.close()
lines.sort(key=lambda line: line.split()[0])
f=open("restart_relaxations_sorted.dat","w")
for line in lines:
	f.write(line+"\n")
f.close()

