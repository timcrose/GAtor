import numpy as np
import os
import shutil
from math import sin, cos

num = 2
delta = 1
geometry_in='/Users/FarrenSCurtis/GA/GA/src/test_integrate/cube'
pool_dir='/Users/FarrenSCurtis/GA/GA/src/test_integrate/initpool'

def main(stoic,seed):
#def pool(geometry_in, num, delta):
	#Import space-delimited geometry.in into 5 columns
	name, x1, y1, z1, species = np.loadtxt(geometry_in, delimiter=' ', unpack=True, dtype='str')

	#Make tables of the strings "lattice_vector", "atom", and "atom_type"	
	name_latvec = name[0:3]	
	name_atom = name[3:]
	atom_type = species[3:]
	name_sig = ['a', 'b', 'c', 'a+b', 'a+c', 'b+c','a+b+c','a_d1','b_d1', 
		   'c_d1','ab_d1', 'ac_d1', 'bc_d1','abc_d1','a_d2','b_d2','c_d2',
		   'ab_d2','ac_d2','bc_d2','abc_d2']
	
	#Make table of x,y,z coordinates of initial geometry
	x_str=x1[3:]
	y_str=y1[3:]
	z_str=z1[3:]
        x_in=[float(i) for i in x_str]
        y_in=[float(i) for i in y_str]
        z_in=[float(i) for i in z_str]
	
	#Shift to Center of geometry coordinates
        x = x_in-np.average(x_in)
        y = y_in-np.average(y_in)
        z = z_in-np.average(z_in)
	r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

	#Find indices of two largest distances (used for signature)
	ind = np.argpartition(r, -2)[-2:]
	max1_ind = ind[1]
	max2_ind = ind[0]
	
	#Separate input lattice vectors
	latvecin_a=[float(x1[0]),float(y1[0]), float(z1[0])]
	latvecin_b=[0.0         ,float(y1[1]), float(z1[1])]
	latvecin_c=[0.0         ,0.0         , float(z1[2])]

	#Generate the new, random lattice vectors and compute distances for signature
	ax = [i for i in range(num)]
	ay = [i for i in range(num)]
        az = [i for i in range(num)]
        by = [i for i in range(num)]
        bz = [i for i in range(num)]
        cz = [i for i in range(num)]

	latcol_1 = [i for i in range(num)]
	latcol_2 = [i for i in range(num)]
	latcol_3 = [i for i in range(num)]
	
	a_dist = [i for i in range(num)]
	b_dist = [i for i in range(num)]
	c_dist = [i for i in range(num)]
	ab_dist = [i for i in range(num)]
	ac_dist = [i for i in range(num)]
	bc_dist = [i for i in range(num)]
	abc_dist = [i for i in range(num)]

	for i in range(num):
		ax[i]=latvecin_a[0]+np.random.uniform(-delta,delta)
		ay[i]=latvecin_a[1]+np.random.uniform(-delta,delta)	
		az[i]=latvecin_a[2]+np.random.uniform(-delta,delta)
		by[i]=latvecin_b[1]+np.random.uniform(-delta,delta)
		bz[i]=latvecin_b[2]+np.random.uniform(-delta,delta)
		cz[i]=latvecin_c[2]+np.random.uniform(-delta,delta)
		latcol_1[i]=[ax[i],0.0,0.0]
		latcol_2[i]=[ay[i],by[i],0.0]
		latcol_3[i]=[az[i],bz[i],cz[i]]
	
		a_dist[i]=np.sqrt(ax[i]**2+ay[i]**2+az[i]**2)
		b_dist[i]=np.sqrt(by[i]**2+bz[i]**2)	
		c_dist[i]=cz[i]
		ab_dist[i]=np.sqrt(ax[i]**2+(ay[i]+by[i])**2+(az[i]+bz[i])**2)
		ac_dist[i]=np.sqrt(ax[i]**2+ay[i]**2+(az[i]+cz[i])**2)
		bc_dist[i]=np.sqrt(by[i]**2+(bz[i]+cz[i])**2)
		abc_dist[i]=np.sqrt(ax[i]**2+(ay[i]+by[i])**2+(az[i]+bz[i]+cz[i])**2)

	#Generate random rotations of input coordinates
	theta = [i for i in range(num)]
	psi = [i for i in range(num)]
	phi = [i for i in range(num)]
	Rxyz = [i for i in range(num)]

	for i in range(num):
		theta[i] = np.random.rand(1) * np.pi * 2
        	psi[i] = np.random.rand(1) * np.pi * 2
        	phi[i] = np.random.rand(1) * np.pi * 2
		Rxyz[i] = np.array([ ((cos(theta[i]) * cos(psi[i])),
                             (-cos(phi[i]) * sin(psi[i])) + (sin(phi[i]) * sin(theta[i]) * cos(psi[i])),
                             (sin(phi[i]) * sin(psi[i])) + (cos(phi[i]) * sin(theta[i]) * cos(psi[i]))),

                             ((cos(theta[i]) * sin(psi[i])),
                             (cos(phi[i]) * cos(psi[i])) + (sin(phi[i]) * sin(theta[i]) * sin(psi[i])),
                             (-sin(phi[i]) * cos(psi[i])) + (cos(phi[i]) * sin(theta[i]) * sin(psi[i]))),

                             ((-sin(theta[i])),
                             (sin(phi[i]) * cos(theta[i])),
                             (cos(phi[i]) * cos(theta[i])))])

	#Rotate input coordinates with random rotation matrices
        rot_xyz = np.asarray(zip(x,y,z))
	xyz=np.zeros((num,len(rot_xyz),3))
	for i in range(len(Rxyz)):
		xyz[i]=rot_xyz
		for j in range(len(rot_xyz)):
			xyz[i][j] = np.dot(Rxyz[i],rot_xyz[j])
	
	#Initialize and compute signatures from 2 furthest atoms for each generated geometry
	a_d1 = [i for i in range(num)]
        b_d1 = [i for i in range(num)]
        c_d1 = [i for i in range(num)]
        ab_d1 = [i for i in range(num)]
        ac_d1 = [i for i in range(num)]
        bc_d1 = [i for i in range(num)]
	abc_d1 = [i for i in range(num)]

        a_d2 = [i for i in range(num)]
        b_d2 = [i for i in range(num)]
        c_d2 = [i for i in range(num)]
        ab_d2 = [i for i in range(num)]
        ac_d2 = [i for i in range(num)]
        bc_d2 = [i for i in range(num)]
        abc_d2 = [i for i in range(num)]
	
	sig = [i for i in range (num)]
	sort_sig = [i for i in range(num)]

	for i in range(num):
		a_d1[i] = np.sqrt((ax[i]-xyz[i][max1_ind][0])**2 + (ay[i]-xyz[i][max1_ind][1])**2 +
			 (az[i]-xyz[i][max1_ind][2])**2)
		a_d2[i] = np.sqrt((ax[i]-xyz[i][max2_ind][0])**2 + (ay[i]-xyz[i][max2_ind][1])**2 + 
			 (az[i]-xyz[i][max2_ind][2])**2)	

		b_d1[i] = np.sqrt((xyz[i][max1_ind][0])**2 + (by[i]-xyz[i][max1_ind][1])**2 + 
			 (bz[i]-xyz[i][max1_ind][2])**2)
		b_d2[i] = np.sqrt((xyz[i][max2_ind][0])**2 + (by[i]-xyz[i][max2_ind][1])**2 +
			 (bz[i]-xyz[i][max2_ind][2])**2)

		c_d1[i] = np.sqrt((xyz[i][max1_ind][0])**2 + (xyz[i][max1_ind][1])**2 +
			 (cz[i]-xyz[i][max1_ind][2])**2)
		c_d2[i] = np.sqrt((xyz[i][max2_ind][0])**2 + (xyz[i][max2_ind][1])**2 + 
			 (cz[i]-xyz[i][max2_ind][2])**2)

		ab_d1[i] = np.sqrt((ax[i]-2*xyz[i][max1_ind][0])**2+(ay[i]+by[i]-2*xyz[i][max1_ind][1])**2+
		       	(az[i]+bz[i]-2*xyz[i][max1_ind][2])**2)
		ab_d2[i] = np.sqrt((ax[i]-2*xyz[i][max2_ind][0])**2+(ay[i]+by[i]-2*xyz[i][max2_ind][1])**2+
			(az[i]+bz[i]-2*xyz[i][max2_ind][2])**2)

		ac_d1[i] = np.sqrt((ax[i]-2*xyz[i][max1_ind][0])**2+(ay[i]-2*xyz[i][max1_ind][1])**2+
			(az[i]+cz[i]-2*xyz[i][max1_ind][2])**2) 
      		ac_d2[i] = np.sqrt((ax[i]-2*xyz[i][max2_ind][0])**2+(ay[i]-2*xyz[i][max2_ind][1])**2+
			(az[i]+cz[i]-2*xyz[i][max2_ind][2])**2)
		
		bc_d1[i] = np.sqrt((-2*xyz[i][max1_ind][0])**2+(by[i]-2*xyz[i][max1_ind][1])**2+
			(bz[i]+cz[i]-2*xyz[i][max1_ind][2])**2)
		bc_d2[i] = np.sqrt((-2*xyz[i][max1_ind][0])**2+(by[i]-2*xyz[i][max2_ind][1])**2+
			(bz[i]+cz[i]-2*xyz[i][max2_ind][2])**2)
		
		abc_d1[i] = np.sqrt((ax[i]-3*xyz[i][max1_ind][0])**2+(ay[i]+by[i]-3*xyz[i][max1_ind][1])**2+
			(az[i]+bz[i]+cz[i]-3*xyz[i][max1_ind][2])**2)
		abc_d2[i] = np.sqrt((ax[i]-3*xyz[i][max2_ind][0])**2+(ay[i]+by[i]-3*xyz[i][max2_ind][1])**2+
			(az[i]+bz[i]+cz[i]-3*xyz[i][max2_ind][2])**2)
		
		sig[i] =[a_dist[i], b_dist[i], c_dist[i], ab_dist[i], ac_dist[i], bc_dist[i], abc_dist[i],
			a_d1[i], b_d1[i], c_d1[i], ab_d1[i], ac_d1[i], bc_d1[i], abc_d1[i],
			a_d2[i], b_d2[i], c_d2[i], ab_d2[i], ac_d2[i], bc_d2[i], abc_d2[i]]
		sort_sig[i] = np.sort(sig[i])
	
	#OUTPUT geometry.in and signature files
#	pool_dir='/home/curtis/Codes/GA/test_init_pool/initial_pool'
#	shutil.rmtree(pool_dir)

#	pool_dir='pwd'
	os.mkdir(pool_dir)
	output_path = pool_dir
	
	for index in range(num):
		output_lattice = zip(name_latvec, latcol_1[index], latcol_2[index], latcol_3[index])
		output_atoms = zip(name_atom, xyz[index], atom_type)
		output_sig = zip(sort_sig[index])
		#output_sig = zip(name_sig,sig[index],sort_sig[index]) 
		#used if you want the name of the distances and the unsorted list to 
		#also be printed in addition to the sorted ones
		folder_name = os.path.join(output_path, "geo{0}".format(index))
		
		try:
            		os.mkdir(folder_name)
        	except OSError, e: 
            		print("Error creating output dir: {0}".format(e))
		with open(os.path.join(folder_name, "geo{0}.in".format(index)), "w+") as f:
			f.write("# This is the rotated output geometry\n")
                        for line in output_lattice:
                                f.write(" ".join(str(x) for x in line) + "\n")
			for line in output_atoms:
                                f.write(" ".join(str(x) for x in line).replace('[','').replace(']','') + "\n")
		with open(os.path.join(folder_name, "geo{0}_sig".format(index)), "w+") as g:
			for line in output_sig:
				g.write(" ".join(str(x) for x in line) + "\n")

	print "Geometries and signatures created"
