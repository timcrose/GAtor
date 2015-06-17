# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:09:38 2015

@author: Patrick Kilecdi

Basic modification of crystal structures
Funtions that read in a struct() and necessary arguments, and return a struct()
"""
import numpy
from structures import structure
from core import output
import copy
import exceptions
from core import user_input
molar_mass={'H':1,'C':12,'N':14,'O':16,'S':32}
lat_interp={0:'lattice_vector_a',1:'lattice_vector_b',2:'lattice_vector_c'}
ui=user_input.get_config()
verbose=ui.get_eval('run_settings','verbose')
nmpc=ui.get_eval('unit_cell_settings','num_molecules')
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

def cell_extension(struct,create_duplicate=True):
    if create_duplicate:
        struct=copy.deepcopy(struct)
    napm=int(len(struct.geometry)/nmpc)
    struct.geometry=numpy.concatenate((struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry,struct.geometry))
#    for i in range (3):
#        struct.geometry+=struct.geometry
    extension=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
    for i in range (1,8):
        for j in range (nmpc):
            mole_translation(struct,i*nmpc+j,napm,frac=extension[i],create_duplicate=False)
    return struct

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

def cm_calculation (struct,atom_list):
    '''
    Reads in a list of atom
    Find the center of mass
    '''
#    print "this is atom_list",atom_list
    cm=[0,0,0]; tm=0;
    for i in range(len(atom_list)):
        tm+=molar_mass[struct.geometry[atom_list[i]][3]]
        for j in range (3):
            cm[j]+=molar_mass[struct.geometry[atom_list[i]][3]]*struct.geometry[atom_list[i]][j]
    for j in range (3):
        cm[j]/=tm
#    print 'this is cm', cm
    return cm
    
def move_molecule_in (struct,create_duplicate=True):
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

def cell_modification (struct,replica,create_duplicate=True):#Replica name is passed in for verbose output
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
    while test and run<10:
        test=False
        run+=1
        if (struct.properties['gamma']>120) or (struct.properties['gamma']<60):
            count+=1
	    if count==1:
		before=print_aims(struct)
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
		before=print_aims(struct)
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
		before=print_aims(struct)
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
    struct=move_molecule_in(struct,False)
    if count>0 and verbose:
	str="------Cell modification required------\n"
	str+="Geometry before modification:\n"
	str+=before
	str+="Geometry after modification:\n"
	str+=print_aims(struct)
	output.local_message(str,replica)
    if count==0 and verbose:
	output.local_message("Cell checked! Reasonable angles observed!\n",replica)
    return struct
    
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
#        print "in structure_handling, this is struct.geometry[i]", struct.geometry[i]
        for j in range(0,3):
#            print "in structure_handling, this is struct.geometry[i][j]", struct.geometry[i][j]
#            st=st+"{0:18.10f}".format(struct.geometry[i][j])
            st=st+"%18.10f" % (struct.geometry[i][j])
        os=os+st+' '+struct.geometry[i][3]+"\n"
    os+='\n'
    return os

