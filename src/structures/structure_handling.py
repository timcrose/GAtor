# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:09:38 2015

@author: Patrick Kilecdi

Basic modification of crystal structures
Funtions that read in a struct() and necessary arguments, and return a struct()
"""
import sys # for debuging!
sys.path.append("/lustre/project/nmarom/gator_shared_duplicate/gator/src/") #for debuging!

import numpy
from structures import structure
from core import output
import copy
import exceptions
from core import file_handler
from core import user_input

lat_interp = {0:'lattice_vector_a',1:'lattice_vector_b',2:'lattice_vector_c'}
ui = user_input.get_config()
verbose = ui.verbose()
all_geo = ui.all_geo()
nmpc = ui.get_eval('unit_cell_settings','num_molecules')
olm = output.local_message

def cell_lower_triangular(struct,create_duplicate=True):
	'''
	Sets the cell back to lower triangular form
	Returns a boolean to indicate whether or not the reset was required
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)

	if (abs(struct.properties["lattice_vector_a"][1])<0.001 and 
            abs(struct.properties["lattice_vector_a"][2])<0.001 and 
	    abs(struct.properties["lattice_vector_b"][2])<0.001):
		return struct

	struct.properties.update(lattice_parameters(struct)) 
	#Add in case not calculated

	new_lattice=lattice_lower_triangular(struct)
	old_lattice=struct.get_lattice_vectors()
	rots = numpy.dot(numpy.transpose(new_lattice),numpy.linalg.inv(numpy.transpose(old_lattice)))
	struct=cell_transform_mat(struct,rots,create_duplicate=False)
	return struct

def cell_translation(struct,trans_vec,create_duplicate=True):
	'''
	Translate the entire structure by trans_vec
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	for i in range (len(struct.geometry)):
		for j in range (3):
			struct.geometry[i][j]+=trans_vec[j]
	return struct

def cell_reflection_z(struct,create_duplicate=True):
	'''
	Flips the cell's z axis 
	Should be sufficient to explore the enantiomers
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	for i in range(len(struct.geometry)):
		struct.geometry[i][2]=-struct.geometry[i][2]
	struct.properties["lattice_vector_a"][2]=-struct.properties["lattice_vector_a"][2]
	struct.properties["lattice_vector_b"][2]=-struct.properties["lattice_vector_b"][2]
	struct.properties["lattice_vector_c"][2]=-struct.properties["lattice_vector_c"][2]
	return struct
def cell_reflection_x(struct,create_duplicate=True):
	'''
	Flips the cell's x axis 
	Should be sufficient to explore the enantiomers
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	for i in range(len(struct.geometry)):
		struct.geometry[i][0]=-struct.geometry[i][0]
	struct.properties["lattice_vector_a"][0]=-struct.properties["lattice_vector_a"][0]
	struct.properties["lattice_vector_b"][0]=-struct.properties["lattice_vector_b"][0]
	struct.properties["lattice_vector_c"][0]=-struct.properties["lattice_vector_c"][0]
	return struct

def cell_rotation(struct,vec=None,theta_deg=None,theta_rad=None,phi_deg=None,phi_rad=None,origin=[0,0,0],deg=None,rad=None,create_duplicate=True):
	if create_duplicate:
		struct=copy.deepcopy(struct)    
	if (deg==None) and (rad==None):
		return False
	if (vec==None) and (((theta_deg==None) and (theta_rad==None)) or ((phi_deg==None) and (phi_rad==None))):
		return False
	if rad==None:
		rad=numpy.deg2rad(deg)
	if (theta_rad==None) and (theta_deg!=None):
		theta_rad=numpy.deg2rad(theta_deg)
	if (phi_rad==None) and (phi_deg!=None):
		phi_rad=numpy.deg2rad(phi_deg)
	if vec==None:
		vec=[numpy.sin(phi_rad)*numpy.cos(theta_rad),numpy.sin(phi_rad)*numpy.sin(theta_rad),numpy.cos(phi_rad)]
	else:
		l=(vec[0]**2+vec[1]**2+vec[2]**2)**0.5
		for j in range (3):
			vec[j]/=l
	c=numpy.cos(rad); s=numpy.sin(rad)
	x,y,z=vec
	mat=[[x*x*(1-c)+c,x*y*(1-c)-z*s,x*z*(1-c)+y*s],
		 [x*y*(1-c)+z*s,y*y*(1-c)+c,y*z*(1-c)-x*s],
		 [x*z*(1-c)-y*s,y*z*(1-c)+x*s,z*z*(1-c)+c]]
	if origin!=[0,0,0]:
		cell_translation(struct,[-j for j in origin],create_duplicate=False)
	for i in range (len(struct.geometry)):
		oldpos=[0,0,0]
		for j in range (3):
			oldpos[j]=struct.geometry[i][j]
		newpos=numpy.dot(mat,oldpos)
		for j in range(3):
			struct.geometry[i][j]=newpos[j]
	struct.properties["lattice_vector_a"]=numpy.dot(mat,struct.properties["lattice_vector_a"])
	struct.properties["lattice_vector_b"]=numpy.dot(mat,struct.properties["lattice_vector_b"])
	struct.properties["lattice_vector_c"]=numpy.dot(mat,struct.properties["lattice_vector_c"])
	if origin!=[0,0,0]:
		cell_translation(struct,origin,create_duplicate=False)
	return struct

def cell_transform_mat(struct,mat,origin=[0,0,0],create_duplicate=True):
	'''
	Transform a structure through a matrix form.
	Allows the input of an origin
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	if origin!=[0,0,0]:
		cell_translation(struct,[-j for j in origin],create_duplicate=False)
	for i in range (len(struct.geometry)):
		oldpos=[0,0,0]
		for j in range (3):
			oldpos[j]=struct.geometry[i][j]
		newpos=numpy.dot(mat,oldpos)
		for j in range(3):
			 struct.geometry[i][j]=newpos[j]
	struct.properties["lattice_vector_a"]=numpy.dot(mat,struct.properties["lattice_vector_a"])
	struct.properties["lattice_vector_b"]=numpy.dot(mat,struct.properties["lattice_vector_b"])
	struct.properties["lattice_vector_c"]=numpy.dot(mat,struct.properties["lattice_vector_c"])
	if origin!=[0,0,0]:
		cell_translation(struct,origin,create_duplicate=False)
	return struct


def cell_extension(struct,create_duplicate=True):
	if create_duplicate:
		struct=copy.deepcopy(struct)
	napm=int(len(struct.geometry)/nmpc)
	struct.geometry=numpy.concatenate((struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry))
	extension=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
	for i in range (1,8):
		for j in range (nmpc):
			mole_translation(struct,i*nmpc+j,napm,frac=extension[i],create_duplicate=False)
	return struct

def cell_populate_molecules(struct,mole_info_list):
	'''
	Reads in a list of molecule information and fills the structure with the molecules defined in ui.conf; 
        requires the mole_info_list to be comprehensive and matching with the number of molecules specified in unit_cell_settings
	mole_info_list should be a list of tuples defined as:
	[COM[0],COM[1],COM[2],is_mirror_reflection,vec[0],vec[1],vec[2],angle]
	'''
	ui=user_input.get_config()
	if len(mole_info_list)!=ui.get_eval("unit_cell_settings","num_molecules"):
		raise RuntimeError("In structure_handling.cell_populate_molecules, the number of molecules to be populated does not match")
		return False
	count=-1
	mole_list=ui.get_eval("unit_cell_settings","molecule_list")
	llist=[]
	for (molename,napm,occurance) in mole_list:
		geo=file_handler.get_molecule_geo(molename)
		for t in range (occurance):
			count+=1
			llist+=list_transformation(geo,mole_info_list[count])
	struct.geometry=structure.convert_array(llist)
	return True	

def lattice_lower_triangular(struct):
	'''
	Returns a list of lattice vectors that corresponds to the a, b, c, alpha, beta, gamma as specified by struct
	! In row vector form !
	'''
	lattice=[[0 for i in range (3)] for j in range (3)]
	a=struct.properties["a"]; b=struct.properties["b"]; c=struct.properties["c"]
	alpha=numpy.deg2rad(struct.properties["alpha"])
	beta=numpy.deg2rad(struct.properties["beta"])
	gamma=numpy.deg2rad(struct.properties["gamma"])
	lattice[0][0] = a
	lattice[0][1] = 0; lattice[0][2] = 0
	lattice[1][0] = numpy.cos(gamma)*b
	lattice[1][1] = numpy.sin(gamma)*b
	lattice[1][2] = 0
	lattice[2][0] = numpy.cos(beta)*c
	lattice[2][1] = (b*c*numpy.cos(alpha) - lattice[1][0]*lattice[2][0])/lattice[1][1]
	lattice[2][2] = (c**2 - lattice[2][0]**2 - lattice[2][1]**2)**0.5
	return lattice

def lattice_parameters(struct):
	'''
	Returns a dictionary of the lattice parameters a, b, c, alpha, beta, gamma
	'''
	parameters={}
	parameters["a"] = numpy.linalg.norm(struct.properties["lattice_vector_a"])
	parameters["b"] = numpy.linalg.norm(struct.properties["lattice_vector_b"])
	parameters["c"] = numpy.linalg.norm(struct.properties["lattice_vector_c"])
	parameters["alpha"] = numpy.rad2deg(numpy.arccos(numpy.dot(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])/parameters["b"]/parameters["c"]))
	parameters["beta"] = numpy.rad2deg(numpy.arccos(numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])/parameters["a"]/parameters["c"]))
	parameters["gamma"] = numpy.rad2deg(numpy.arccos(numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])/parameters["a"]/parameters["b"]))
	return parameters
			

def mole_translation(struct,mn,napm,vec=None,frac=None,create_duplicate=True):
	if (vec==None) and (frac==None):
		raise exceptions.RuntimeError("Please input at least one type of translation vector into structure_handling.mole_translation")
	if create_duplicate:
		struct=copy.deepcopy(struct)
	if vec==None:
		vec=[0,0,0]
		for i in range (3):
			for j in range (3):
				vec[j]+=struct.properties[lat_interp[i]][j]*frac[i]
	for i in range (mn*napm,mn*napm+napm):
		for j in range (3):
			struct.geometry[i][j]+=vec[j]
	return struct

def mole_recognize(struct):
	'''
	Interpret the geometry of struct in terms of the molecules given in ui.conf
	Returns a tuple of [COM[0],COM[1],COM[2],is_mirror_reflection,vec[0],vec[1],vec[2],angle,residual]
	where vec is the axis of rotation, and angle is the angle being rotated
	the orientation is determined according to the standard geometry
	found in run_calcs/molecules
	'''
	molecules=ui.get_eval('unit_cell_settings','molecule_list')
	sum=0;count=-1
	result=[]
	struct=copy.deepcopy(struct)
	for (molename,napm,occurance) in molecules:
		geo=file_handler.get_molecule_geo(molename)
		for j in range (occurance):
			com=cm_calculation(struct,range(sum,sum+napm))
			orient=mole_get_orientation(struct,geo,range(sum,sum+napm),com,create_duplicate=False)
			if orient==False:
				return False
			result.append(com+orient)
				#if a molecule is torn apart by relaxation, mole_get_orientation will return False, thus causing error for connecting the list
	return result

def mole_get_orientation(struct,atom_list,geo,com=None,tol=0.1,create_duplicate=True):
	'''
	Check if the list of atoms provided in struct fit the molecule specified by mole_name
	if yes, then the orientation of the molecule is given in a list
	[is_mirror_reflection,vec[0],vec[1],vec[2],angle,residual]
	Make sure molecule defined by geo has its COM at the origin
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
		geo = copy.deepcopy(geo)
	if com==None:
		com=cm_calculation(struct,atom_list)
	if len(atom_list)!=len(geo):
		raise RuntimeError("In structure_handling.mole_get_orientation, the length of atom_list does not match that of geo")
		return False
	for atom in atom_list:
		for j in range (3):
			struct.geometry[atom][j]-=com[j]
	lll=0
	match_molecule_length_requirement=0.0001
	match_molecule_cross_tolerance=0.001
	match_molecule_tolerance=tol
	result=[False,0,0,0,0,0]
	while lll<2:
		lll+=1
		vec1=vec2=None
		i=0
		rotation_axis=[0,0,0]
		chosen=None
		while (vec2==None) and (i<len(geo)):
			diff=[struct.geometry[atom_list[i]][j]-geo[i][j] for j in range (3)]
			leng=numpy.linalg.norm(diff)
			if leng<0.0001:
				rotation_axis=struct.geometry[atom_list[i]][:3]
				if numpy.linalg.norm(rotation_axis)>match_molecule_cross_tolerance:
					break
			if leng>match_molecule_length_requirement:
				if vec1==None:
					vec1=diff
					chosen=i
					i+=1
					continue
				if numpy.linalg.norm(numpy.cross(vec1,diff))>match_molecule_cross_tolerance:
					vec2=diff
					break
			i+=1
		if numpy.linalg.norm(rotation_axis)<match_molecule_cross_tolerance:
			if vec2==None:
				vec2=[0,0,0]
			rotation_axis=numpy.cross(vec1,vec2)
		rl=numpy.linalg.norm(rotation_axis)
		if rl>match_molecule_cross_tolerance:
			for j in range (3):
				rotation_axis[j]/=rl
			if chosen==None:
				chosen=0
				while (chosen<len(geo)) and (numpy.linalg.norm(numpy.cross(rotation_axis,struct.geometry[atom_list[chosen]]))<match_molecule_cross_tolerance):
					chosen+=1
			v1=[struct.geometry[atom_list[chosen]][x] for x in range(3)]
			v2=[geo[chosen][x] for x in range(3)]
			m1=numpy.dot(v1,rotation_axis)
			v1=[v1[j]-m1*rotation_axis[j] for j in range (3)]
			m2=numpy.dot(v2,rotation_axis)
			v2=[v2[j]-m2*rotation_axis[j] for j in range (3)]
			l1=numpy.linalg.norm(v1); l2=numpy.linalg.norm(v2)
			direction=numpy.cross(v1,v2)
			rad=numpy.arcsin(numpy.linalg.norm(direction)/(l1*l2))
			if numpy.sign(numpy.dot(v1,v2))==-1:
				rad=numpy.pi-rad
			if numpy.linalg.norm(direction)>match_molecule_cross_tolerance:
				rad*=numpy.sign(numpy.dot(direction,rotation_axis))
			geo_1=list_rotation(geo,vec=rotation_axis,rad=-rad,create_duplicate=True)
			resi=0
			for i in range (len(geo)):
				resi+=numpy.linalg.norm([struct.geometry[atom_list[i]][j]-geo_1[i][j] for j in range (3)])**2
#			print struct.get_geometry_atom_format()
			for j in range (3):
				result[j+1]=rotation_axis[j]
			result[4]=numpy.rad2deg(rad)
			result[5]=resi		

			if resi<match_molecule_tolerance:
				for j in range (3):
					result[j+1]=rotation_axis[j]
				result[4]=numpy.rad2deg(rad)
				result[5]=resi
				return result
		list_reflection_z(geo,False)
		result[0]=True
	return False

def cm_calculation (struct,atom_list):
	'''
	Reads in a list of atom
	Find the center of mass
	'''
#    print "this is atom_list",atom_list
	cm=[0,0,0]; tm=0;
	ui=user_input.get_config()
	for i in range(len(atom_list)):
		tm+=ui.get_eval("molar_mass",struct.geometry[atom_list[i]][3])
		for j in range (3):
			cm[j]+=ui.get_eval("molar_mass",struct.geometry[atom_list[i]][3])*struct.geometry[atom_list[i]][j]
	for j in range (3):
		cm[j]/=tm
#    print 'this is cm', cm
	return cm
	
def move_molecule_in (struct,nmpc=nmpc, create_duplicate=True):
	'''
	Translate the molecules by the cell vector such that their center of mass lies within the cell
	'''
	if create_duplicate:
		struct=copy.deepcopy(struct)
	napm=int(len(struct.geometry)/nmpc)
	lattice=[struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"]]
	lattice=numpy.transpose(lattice)
	latinv=numpy.linalg.inv(lattice)    
	for i in range (nmpc):
		cm=cm_calculation(struct,range(i*napm,i*napm+napm))
		frac=numpy.dot(latinv,cm)
		for j in range (0,3):
			lat=lat_interp[j]
			vec=struct.properties[lat]
			if (frac[j]<-0.0001):
				kk=int(-frac[j]+1)
				for k in range(i*napm,i*napm+napm):
					for l in range (3):
						struct.geometry[k][l]+=kk*vec[l]
			elif (frac[j]>0.99999):
				kk=int(frac[j]+0.00001)
				for k in range(i*napm,i*napm+napm):
					for l in range (3):
						struct.geometry[k][l]-=kk*vec[l]
	return struct

def angle(l1,l2):
	return (numpy.rad2deg(numpy.arccos(numpy.dot(l1,l2)/(numpy.linalg.norm(l1)*numpy.linalg.norm(l2)))))

def cell_modification(struct,napm=None,create_duplicate=True):
	'''
	Cell modification using Niggli reduction
	'''
	if create_duplicate:
		struct = copy.deepcopy(struct)
	lats = struct.get_lattice_vectors()
	from spglib import niggli_reduce
	reduced_lats = 	niggli_reduce(lats)
	del(struct.properties["lattice_vector_a"])
	del(struct.properties["lattice_vector_b"])
	del(struct.properties["lattice_vector_c"])
	struct.set_lattice_vectors(reduced_lats)
	nmpc = len(struct.geometry)/napm
	cell_lower_triangular(struct,False)
	move_molecule_in(struct,nmpc,False)
	return struct

def cell_modification_old (struct,replica=ui.get_replica_name(),create_duplicate=True):#Replica name is passed in for verbose output
	'''
	Method found in the 2011 Lonie paper
	Make the skewed angle correct
	'''
	napm=int(len(struct.geometry)/nmpc)
	if create_duplicate:
		struct=copy.deepcopy(struct)
	test=True
	run=0
	count=0
	try: #Updates the missing struct information
		gamma=struct.properties['gamma']
		a=struct.properties['a']
	except:
		struct.properties.update(lattice_parameters(struct))

	while test and run<10:
		test=False
		run+=1
		if (struct.properties['gamma']>120) or (struct.properties['gamma']<60):
			count+=1
		if count==1:
#			before=print_aims(struct)
			test=True
			if struct.properties['a']>=struct.properties['b']:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])/(numpy.linalg.norm(struct.properties["lattice_vector_b"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_a"][j]-=c*struct.properties["lattice_vector_b"][j]
				struct.properties['a']=numpy.linalg.norm(struct.properties['lattice_vector_a'])
			else:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])/(numpy.linalg.norm(struct.properties["lattice_vector_a"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_b"][j]-=c*struct.properties["lattice_vector_a"][j]
				struct.properties['b']=numpy.linalg.norm(struct.properties['lattice_vector_b'])
			struct.properties['alpha']=angle(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])
			struct.properties['beta']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])
			struct.properties['gamma']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])

		if (struct.properties['beta']>120) or (struct.properties['beta']<60):
			count+=1

		if count==1:
#			before=print_aims(struct)
			test=True
			if struct.properties['a']>=struct.properties['c']:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_c"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_a"][j]-=c*struct.properties["lattice_vector_c"][j]
				struct.properties['a']=numpy.linalg.norm(struct.properties['lattice_vector_a'])
			else:
				frac=numpy.dot(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_a"])**2)
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				c=numpy.round(frac)
				for j in range (3):
					struct.properties["lattice_vector_c"][j]-=c*struct.properties["lattice_vector_a"][j]
				struct.properties['c']=numpy.linalg.norm(struct.properties['lattice_vector_c'])
			struct.properties['alpha']=angle(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])
			struct.properties['beta']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])
			struct.properties['gamma']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])

		if (struct.properties['alpha']>120) or (struct.properties['alpha']<60):
			count+=1

		if count==1:
#			before=print_aims(struct)
			test=True
			#d=struct.properties
			if struct.properties['b']>=struct.properties['c']:
				frac=numpy.dot(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_c"])**2)
				c=numpy.round(frac)                
				#c=numpy.ceil(abs(frac))*numpy.sign(frac)
				for j in range (3):
					struct.properties["lattice_vector_b"][j]-=c*struct.properties["lattice_vector_c"][j]
				struct.properties['b']=numpy.linalg.norm(struct.properties['lattice_vector_b'])
			else:
				frac=numpy.dot(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])/(numpy.linalg.norm(struct.properties["lattice_vector_b"])**2)
				#c=numpy.round(frac)                
				c=numpy.ceil(abs(frac))*numpy.sign(frac)
				for j in range (3):
					struct.properties["lattice_vector_c"][j]-=c*struct.properties["lattice_vector_b"][j]
				struct.properties['c']=numpy.linalg.norm(struct.properties['lattice_vector_c'])
			struct.properties['alpha']=angle(struct.properties["lattice_vector_b"],struct.properties["lattice_vector_c"])
			struct.properties['beta']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_c"])
			struct.properties['gamma']=angle(struct.properties["lattice_vector_a"],struct.properties["lattice_vector_b"])

	struct=move_molecule_in(struct,create_duplicate=False)
	struct=cell_lower_triangular(struct,False)

	if count>0 and verbose:
		st="-- Cell modification required for structure " 
		if struct.struct_id!=None:
			st += str(struct.struct_id)
		else:
			st += "\n"
		
#		if all_geo: 
#			st += "Geometry after modification:\n"
#			st += print_aims(struct)
		output.local_message(st,replica)
	if count==0 and verbose:
		st = "-- Cell modification uncessary for structure "
		if struct.struct_id!=None:
			st += str(struct.struct_id)
		output.local_message(st,replica)
	return struct

def cell_check(struct,replica):
	'''
	Checks if a structure has reasonable cell volume, lattice vector lengths, and distance between molecules and atoms
	Should only be performed after cell_modification is done
	'''
	sname = "cell_check_settings"
#	output.local_message("--------Begin cell check--------",replica)
	struct = copy.deepcopy(struct)
	struct = cell_modification(struct,
				   int(struct.get_n_atoms()/nmpc),
				   create_duplicate=False)

	#Volume check
	standard_volume=ui.get_eval(sname,"target_volume")
	ucv = struct.get_unit_cell_volume()
	if standard_volume!=None:
		if verbose: output.local_message("-- Volume check:", replica)
		if verbose: output.local_message("-- New structure's volume: %f" % ucv,replica)
		upper_volume=ui.get_eval(sname,"volume_upper_ratio")*standard_volume
		lower_volume=ui.get_eval(sname,"volume_lower_ratio")*standard_volume
		if verbose: output.local_message("-- Allowed range: %f - %f" % (lower_volume,upper_volume),replica)

		if ucv > upper_volume or ucv < lower_volume:
			if verbose: 
				output.local_message("-- Failed volume check\n",replica)
			return False
		elif verbose:
			output.local_message("-- Passed volume check\n",replica)
	elif verbose:
		output.local_message("##No volume check is called##",replica)


	#Lattice vector principal component check
	#Note: this check assumes that the lattice vectors are put into lower triangular form
	#This is done in run_GA.py
	if ui.get_boolean(sname,"lattice_vector_check"):
		if verbose: output.local_message("Lattice vector principal component check:",replica)
		lv = ["lattice_vector_a","lattice_vector_b","lattice_vector_c"]
		if verbose:
			for i in range(3):
				output.local_message("-- Principal component of " + 
				lv[i]+": "+str(struct.properties[lv[i]][i]),replica)

		lower_bound = ui.get_eval(sname,"lattice_vector_lower_ratio")*ucv**(1/3.0)
		upper_bound = ui.get_eval(sname,"lattice_vector_upper_ratio")*ucv**(1/3.0)
		if ui.has_option(sname,"lattice_vector_lower_bound"):
			lower_bound = max(lower_bound,
			ui.get_eval(sname,"lattice_vector_lower_bound"))

		if ui.has_option(sname,"lattice_vector_upper_bound"):
			upper_bound = min(upper_bound,
			ui.get_eval(sname,"lattice_vector_upper_bound"))

#		cell_lower_triangular(struct,False)
		if verbose:
			output.local_message("-- Acceptable range: %f - %f"
					% (lower_bound, upper_bound),replica)

		for i in range (3):
			if struct.properties[lv[i]][i] > upper_bound\
			or struct.properties[lv[i]][i] < lower_bound:
				if verbose:
					output.local_message("-- %s failed principal component check\n" % lv[i],replica)
					return False
		if verbose:
			output.local_message("-- Passed principal component check\n", replica)
	elif verbose:
		output.local_message("##No lattice principal component check is called##", replica)


	#Now enters a series of closeness checks
	#These all require the cells to be orthogonalized
	if ui.has_option(sname,"COM_distance_check"):
		if verbose: output.local_message("COM distance check:",replica)
		lowerbound = ui.get_eval(sname,"COM_distance_check")
		if verbose:
			output.local_message("-- Minimum distance: %f " 
						% lowerbound,replica)
		if COM_distance_check(struct,nmpc,lowerbound):
			output.local_message("-- Passed COM distance check\n",replica)
		else:
			output.local_message("-- Failed COM distance check\n",replica)
			return False
	elif verbose:
		output.local_message("##No COM distance check is called##")

	if verbose:
		output.local_message("Combined interatomic distance check:",replica)
	if combined_distance_check(struct,replica):
		if verbose:
			output.local_message("-- Passed all atomic distance check(s)\n",replica)
	else:
		if verbose:
			output.local_message("-- Failed atomic distance check\n",replica)
		return False


	return True

def COM_distance_check(struct,nmpc=None,lowerbound=None):
	'''
	Calculates and returns the minimal distance between the center of mass of the molecules in the unit cell
	if lowerbound is specified, then returns a True system if the minimal distance is greater than the lowerbound
	'''
	if nmpc==None:
		nmpc=struct.properties["nmpc"]
	napm=int(len(struct.geometry)/nmpc)
	
	cmlist=[cm_calculation(struct,range(napm*i,napm*i+napm)) for i in range (nmpc)] #Calculates the center of mass of all molecules
	tr=([[0,0,0],[0,0,1],[0,0,-1],[0,1,0],[0,-1,0],[0,1,1],
             	[0,1,-1],[0,-1,1],[0,-1,-1],[1,0,0],[-1,0,0],[1,0,1],
		[1,0,-1],[-1,0,1],[-1,0,-1],[1,1,0],[1,-1,0],[-1,1,0],
		[-1,-1,0],[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],
		[-1,1,-1],[-1,-1,1],[-1,-1,-1]]) 
	#Each molecule will have to be compared with the 27 equivalent of the other one in order to conduct the distance check
	min_dist = min([numpy.linalg.norm(struct.properties["lattice_vector_a"]),numpy.linalg.norm(struct.properties["lattice_vector_b"]),numpy.linalg.norm(struct.properties["lattice_vector_c"])])
	if min_dist<lowerbound:
		return False


	for m1 in range (nmpc-1):
		for m2 in range (m1+1,nmpc):
			for tr_choice in range (27):
				new_cm=[cmlist[m2][j]+struct.properties["lattice_vector_a"][j]*tr[tr_choice][0]+struct.properties["lattice_vector_b"][j]*tr[tr_choice][1]+struct.properties["lattice_vector_c"][j]*tr[tr_choice][2] for j in range (3)] #Move the molecule by the fractional vector specified by tr[tr_choice]
				diff=[cmlist[m1][j]-new_cm[j] for j in range (3)]
				if min_dist==None or numpy.linalg.norm(diff)<min_dist:
					min_dist=numpy.linalg.norm(diff)
				if lowerbound!=None and min_dist<lowerbound:
					return False
	if lowerbound==None:
		return min_dist
	elif min_dist>=lowerbound:
		return True
	else:
		return False

def combined_distance_check(struct,replica):
	'''
	Combines full atomic distance check, interatomic distance check, 
	and specific radius check in one function
	'''
	struct = copy.deepcopy(struct)
	sname = "cell_check_settings"
	min_dist = ui.get_eval(sname,"full_atomic_distance_check")
	if verbose:
		olm("-- Full atomic distance check minimum distance between all pairs of atoms: %f" % min_dist, replica)
	
	if ui.has_option(sname,"interatomic_distance_check"):
		min_dist_2 = ui.get_eval(sname,"interatomic_distance_check")
		if verbose:
			olm("-- Interatomic distance check enforcing minimum distance interatomic pairs of atoms: %f" % min_dist_2, replica)
	else:
		min_dist_2 = None

	if ui.has_option(sname,"specific_radius_proportion"):
		sr = ui.get_eval(sname,"specific_radius_proportion")
		if verbose:
			olm("-- Specific radius check enforcing minimum proportion: %f" % sr, replica)
	else:
		sr = None

	
	#First conducts atomic distance check between atom pairs from the same molecule
	napm = struct.get_n_atoms()/nmpc
	for i in range (nmpc):
		for j in range(i*napm,(i+1)*napm-1):
			for k in range (j+1,(i+1)*napm):
				if struct.get_atom_distance(j,k) < min_dist:
					if verbose:
						olm("-- Atoms %i and %i too close"
						% (j,k) ,replica)
					return False

	if min_dist_2!=None:
		min_dist = max(min_dist,min_dist_2)

	total_atoms=len(struct.geometry)

	#Next make sure molecules won't run into its own replicas

	tr = [[x,y,z] for x in range (0,2) for y in range (0,2) for z in range (0,2) if x!=0 or y!=0 or z!=0]
	radius = ui.get_section_as_dict("specific_radius",eval=True)
	for a1 in range (total_atoms):
		for tr_choice in tr:
			new_apos = [struct.geometry[a1][j]+struct.properties["lattice_vector_a"][j]*tr_choice[0]+struct.properties["lattice_vector_b"][j]*tr_choice[1]+struct.properties["lattice_vector_c"][j]*tr_choice[2] for j in range (3)]
			for a2 in range (a1-a1%napm,a1-a1%napm+napm):
				diff = [struct.geometry[a2][j]-new_apos[j] for j in range (3)]
				dist = numpy.linalg.norm(diff)
				if dist < min_dist and verbose:
					olm("-- Atoms %i and %i too close"
					% (a1,a2), replica)
					olm("-- Distance: "+str(dist),replica)
	
				if dist<min_dist:
					return False

				a1_s = struct.geometry[a1][3].lower()
				a2_s = struct.geometry[a2][3].lower()
				if (sr!=None and a1_s in radius and a2_s in radius 
				and dist<(radius[a1_s]+radius[a2_s])*sr):
					if verbose:
						olm("-- Atoms %i and %i too close",\
						replica)
						olm("-- Specific radius proportion: %f" % (dist/radius[a1_s]/radius[a2_s]),replica)
					return False
	

	#Then compare each atom pair from different molecules
	tr=([[0,0,0],[0,0,1],[0,0,-1],[0,1,0],[0,-1,0],[0,1,1],
           [0,1,-1],[0,-1,1],[0,-1,-1],[1,0,0],[-1,0,0],[1,0,1],
           [1,0,-1],[-1,0,1],[-1,0,-1],[1,1,0],[1,-1,0],[-1,1,0],
           [-1,-1,0],[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],
           [-1,1,-1],[-1,-1,1],[-1,-1,-1]])
	struct = move_molecule_in(struct,nmpc=total_atoms) #Move in to prevent too far away from cell
	for a1 in range (total_atoms-napm):
		for tr_choice in tr:
			new_apos = [struct.geometry[a1][j]+struct.properties["lattice_vector_a"][j]*tr_choice[0]+struct.properties["lattice_vector_b"][j]*tr_choice[1]+struct.properties["lattice_vector_c"][j]*tr_choice[2] for j in range (3)]
			start = (a1/napm+1)*napm
			for a2 in range (start,total_atoms): 
			#Atoms should not be compared to those from the same molecule
				diff = [struct.geometry[a2][j]-new_apos[j] for j in range (3)]
				dist = numpy.linalg.norm(diff)

				if dist < min_dist and verbose:
					olm("-- Atoms %i and %i too close"
					% (a1,a2), replica)
					olm("-- Distance: "+str(dist),replica)

				if dist < min_dist:
					return False

				a1_s = struct.geometry[a1][3].lower()
				a2_s = struct.geometry[a2][3].lower()
				if (sr!=None and a1_s in radius and a2_s in radius 
				and dist<(radius[a1_s]+radius[a2_s])*sr):
					if verbose:
						olm("-- Atoms %i and %i too close",\
						replica)
						olm("-- Specific radius proportion: %f" % (dist/radius[a1_s]/radius[a2_s]),replica)
	
					return False
	
	return True
			

def list_rotation(llist,vec=None,theta_deg=None,theta_rad=None,phi_deg=None,phi_rad=None,origin=[0,0,0],deg=None,rad=None,create_duplicate=True):
	'''
	Same as cell_rotation, but now acts upon a list of coordinates
	'''
	if create_duplicate:
		llist=copy.deepcopy(llist)
	if (deg==None) and (rad==None):
		return False
	if (vec==None) and (((theta_deg==None) and (theta_rad==None)) or ((phi_deg==None) and (phi_rad==None))):
		return False
	if rad==None:
		rad=numpy.deg2rad(deg)
	if (theta_rad==None) and (theta_deg!=None):
		theta_rad=numpy.deg2rad(theta_deg)
	if (phi_rad==None) and (phi_deg!=None):
		phi_rad=numpy.deg2rad(phi_deg)
	if vec==None:
		vec=[numpy.sin(phi_rad)*numpy.cos(theta_rad),numpy.sin(phi_rad)*numpy.sin(theta_rad),numpy.cos(phi_rad)]
	else:
		l=(vec[0]**2+vec[1]**2+vec[2]**2)**0.5
		for j in range (3):
			vec[j]/=l
	c=numpy.cos(rad); s=numpy.sin(rad)
	x,y,z=vec
	mat=[[x*x*(1-c)+c,x*y*(1-c)-z*s,x*z*(1-c)+y*s],
		 [x*y*(1-c)+z*s,y*y*(1-c)+c,y*z*(1-c)-x*s],
		 [x*z*(1-c)-y*s,y*z*(1-c)+x*s,z*z*(1-c)+c]]
	if origin!=[0,0,0]:
		list_translation(llist,[-j for j in origin],create_duplicate=False)
	for i in range (len(llist)):
		oldpos=[0,0,0]
		for j in range (3):
			oldpos[j]=llist[i][j]
		newpos=numpy.dot(mat,oldpos)
		for j in range(3):
			llist[i][j]=newpos[j]
	if origin!=[0,0,0]:
		list_translation(llist,origin,create_duplicate=False)
	return llist

def list_translation(llist,trans_vec,create_duplicate=True):
	'''
	Similar to struct_translation, but acts upon a list of coordinates
	'''
	if create_duplicate:
		llist=copy.deepcopy(llist)
	for i in range (len(llist)):
		for j in range (3):
			llist[i][j]+=trans_vec[j]
	return llist

def list_reflection_z(llist,create_duplicate=True):
	'''
	Flip the sign of the z component of the coordinates specified by llist
	'''
	if create_duplicate:
		llist=copy.deepcopy(llist)
	for i in range(len(llist)):
		llist[i][2]=-llist[i][2]
	return llist

def list_transformation(llist,info,create_duplicate=True):
	'''
	Performs transformation as specified by info on llist
	info format: [COM[0],COM[1],COM[2],is_mirror_reflection,vec[0],vec[1],vec[2],angle]
	'''
	if create_duplicate:
		llist=copy.deepcopy(llist)
	if info[3]:
		list_reflection_z(llist,False)
	list_rotation(llist,vec=info[4:7],deg=angle,create_duplicate=False)
	list_translation(llist,info[0:3],False)
	return llist
	

 
def print_aims(struct):
	os=''        
	st="lattice_vector"
	for j in range(0,3):
		st=st+"{0:18.10f}".format(struct.properties["lattice_vector_a"][j])
	os=os+st+"\n"
	st="lattice_vector"
	for j in range(0,3):
		st=st+"{0:18.10f}".format(struct.properties["lattice_vector_b"][j])
	os=os+st+"\n"
	st="lattice_vector"
	for j in range(0,3):
		st=st+"{0:18.10f}".format(struct.properties["lattice_vector_c"][j])
	os=os+st+"\n"
	for i in range (len(struct.geometry)):
		st='atom'
		for j in range(0,3):
			st=st+"%18.10f" % (struct.geometry[i][j])
		os=os+st+' '+struct.geometry[i][3]+"\n"
	os+='\n'
	return os

def main():
	struct=structure.Structure()
	struct.set_property("lattice_vector_a",[10,0,0])
	struct.set_property("lattice_vector_b",[0,10,0])
	struct.set_property("lattice_vectory_c",[0,0,10])
	mole_info=[[1,2,3,False,0,0,1,72],[5,6,7,True,0,1,0,0]]
	cell_populate_mole(struct,mole_info)
	print print_aims(struct)


if __name__ == '__main__':
	main2()
	

