'''
Created by farren to re-relax structures lost
'''
#Find repeating replicas in relaxation file
filename = "restart_relaxations.in"
lines = [line.rstrip('\n') for line in open(filename)]
seen = set()
uniq_reps = []
for line in lines:
    if line.split()[0] not in seen:
        uniq_reps.append(line.split()[0])
        seen.add(line.split()[0])

test = [line.rstrip('\n') for line in open(filename)]
print str(test)

#Make array of text file
arr =[]
with open (filename) as f:
  lines = f.readlines()
  for line in lines:
      words = line.split()
      arr.append(words)	

#Make arrays of data from each replica
start_count = 0
end_count = 0
end_iter_count = 0
rep_ars =[[] for i in range(len(uniq_reps))]
iter =[[] for i in range(len(uniq_reps))]
final = [[] for i in range(len(uniq_reps))]
for i in range(len(uniq_reps)):
	for line in arr:
		if line[0] == uniq_reps[i]:
			rep_ars[i].append((line[0],line[1],line[2]))
for i in range(len(rep_ars)):
	for line in rep_ars[i]:
		iter[i].append(line[1])
	maxi = max(iter[i])
	print "Maxi", maxi
	for line in rep_ars[i]:
		if line[1] == maxi:
			final[i].append((line[0],line[1],line[2]))
for i in range(len(rep_ars)):
	for line in final[i]:
		print line[2]
		if line[2] == "started_relaxing:":
			start_count = end_count + 1
		if line[2] == "finished_relaxing:":
			end_count = end_count +1
		if line[2] == "finished_iteration:":
			end_iter_count = end_iter_count +1
        if start_count == end_count:
		if end_count == end_iter_count:
			print "done"
		# struct has been added to collection
		else:
			print "finished relaxing but not added to collection"
	                #struct finished relaxing but hasn't been added to collection
                        # Make new Strcutre() 
	elif start_count != end_count:
		print "ho!"
		#make Structure() from geometry in folder
		#send through FHI-aims which returns new_struct

	#pass struct throught cell checks
	#add to collection (don't bother about convergence)

	 	

