import os
import sys
import pickle
import time
import string
import commands,subprocess
import re
import numpy as np
	
path_sep = os.sep
workdir = os.getcwd()

K_BOTZMAN = 0.001987191
ATOM_MASS = 1.6605389

"""
try assigning same weights to vol and hvap in the next run
"""
w_vol  = 250.0
w_hvap = 120.0
w_sfe  =  1.0

mass = 0.0
T_Sim = 0.0
density = 0.0
volume = 0.0
hvap = 0.0
Wrong_SFE = 10000.0 
dG_sfe = Wrong_SFE


List_Mol = []
List_volume = []
List_hvap = []
List_sfe = []
List_T_Sim = []
List_sfe_A = []

Job_List = []
gasJobList=[]

n_vdw_type = 0
vdw_list = []
vdw_list_map = []
df_dEmin_List = []
df_dRmin_List = []

def query_vdw_type(ChemName):
    for i in range(0,n_vdw_type):
        if ChemName in vdw_list_map[i]: #if(vdw_list[i] == ChemName):
            return i
    return -1
    
def read_vdw_type():
    n_vdw_type = 0
#    f=open(".." + path_sep + ".." + path_sep + "vdw-to-fit.txt","r")
    f=open(".." + path_sep + "vdw-to-fit.txt","r")
    line = f.readline()
    
    while 1:
        line = f.readline()
        if not line: break
        szBuff = line.split()
        if( (len(szBuff) >= 7) and (szBuff[0][0]) != '!') :
            vdw_list.append(szBuff[0])
            n_vdw_type = n_vdw_type + 1
            
            df_dEmin_List.append(0.0)
            df_dRmin_List.append(0.0)
            
    f.close()

    fin=open(".." + path_sep + "equivalent.txt","r")
    for line in fin:
        token=line.rstrip().split()
        if len(set(token).intersection(set(vdw_list)))>0:
                vdw_list_map.append(token)
    fin.close()
    return n_vdw_type


def Accumulate_Grad_dE(nRec,vdw_Names_list,df_Emin_list,df_Rmin_list,Scale):
    global df_dEmin_List,df_dRmin_List
    for i in range(0,nRec):
        Idx = query_vdw_type(vdw_Names_list[i])
        if(Idx >= 0):
            df_dEmin_List[Idx] = df_dEmin_List[Idx] + df_Emin_list[i]*Scale
            df_dRmin_List[Idx] = df_dRmin_List[Idx] + df_Rmin_list[i]*Scale


def read_vdw_grad(szName):
    nRec = 0
    vdw_Names_list = []
    df_Emin_list = []
    df_Rmin_list = []
    fIn = open(szName,"r")
    while(1):
        line = fIn.readline()
        if not line: break
        szBuff = line.split()
        if(len(szBuff) == 3):
            vdw_Names_list.append(szBuff[0])
            df_Emin_list.append(string.atof(szBuff[1]))
            df_Rmin_list.append(string.atof(szBuff[2]))
            nRec = nRec + 1
  
    fIn.close()

    return (nRec,vdw_Names_list,df_Emin_list,df_Rmin_list)

def get_experimenal_values(molecule):
    global mass,T_Sim,density,hvap,dG_sfe
    #pattern=re.compile(molecule)
    fin = open('/lcrc/project/Drude/chetan/workspace/allparameters-correct.dat','r')
    for line in fin:
	token = line.rstrip().split()
	if token[2]==molecule:
    		mass = float(token[3])
    		density = float(token[4])
    		T_Sim = float(token[5])
    		hvap = float(token[6])
    		dG_sfe = float(token[7])
    fin.close()

def isAtomTypeMatch(typeList,atomTypeToMatch,optAtomType):
    matchTypes=[]
    for ts in typeList:
	if optAtomType in ts:
	    matchTypes = ts
    return atomTypeToMatch in matchTypes
		

def getTemp(idMol):
    fin = open('/lcrc/project/Drude/chetan/workspace/allparameters-correct.dat','r')
    for line in fin:
        token = line.rstrip().split()
        if token[2]==idMol:
                T_Sim = token[5]
		break
    fin.close()
    return T_Sim	


n_vdw_type = read_vdw_type()

nMol = 0
f=open("../exp_set.txt","r")

atomType="OG311"
run = sys.argv[1]

#submit all liquid jobs before submitting gas jobs
allMols=[]
for line in f:
    mol_id = line.rstrip().split()[0]
    allMols.append(mol_id)
    newdir = run + path_sep + mol_id
    os.chdir(newdir)
   
    get_experimenal_values(mol_id)

    volume = ATOM_MASS * mass / density
    str_T_Sim = '%f' %(T_Sim)

    List_T_Sim.append(T_Sim)
    List_hvap.append(hvap)
    List_volume.append(volume)
    List_sfe.append(dG_sfe)
    List_Mol.append(mol_id)

    if(Done == 0):
        molnumber = mol_id[2:]

        if(os.path.isfile("./liquid/dE_dvdw_liquid.txt")):
	    mols_to_process.append(mol_id) 	
    nMol = nMol + 1

    os.chdir(workdir)

f.close()

os.chdir(workdir)

# To process the data
fLog = open("mol_obj_grad.txt","r")

obj_func = 0.0
for id in range(0,nMol):
    mol_id = List_Mol[id]
    T_Sim = List_T_Sim[id]
    RT = T_Sim * K_BOTZMAN
    newdir = run + path_sep + mol_id
    os.chdir(newdir)

    boxName = "liquid/drude-box.pdb"
    command1 = "tail -n 3 "+boxName+" | grep ""ATOM"" | awk '{print $5}'"
    process = subprocess.Popen(command1,stdout=subprocess.PIPE, shell=True)
    N_Mol_l = int(process.communicate()[0])

    if os.path.isfile("enthalpy-vol-summary.txt"):
	fin = open("enthalpy-vol-summary.txt","r")	
	for line in fin:
		token = line.rstrip().split()
		if token[0]=="(DeltaH)average:":H_vap_Sim = float(token[1])
		elif token[0]=="(Volume)average,std:":MolVol_l = float(token[1])/N_Mol_l

    ### Gradients for the gas phase are computed from Gas job script
   	
    dHere = os.getcwd()
    #print dHere
    df_Emin_L,df_Rmin_L = read_vdw_grad(dHere+'/liquid/dE_dvdw_liquid.txt',atomType)
    ScaleL = -2.0*w_hvap*(H_vap_Sim/List_hvap[id] - 1.0)/List_hvap[id]

    df_Emin_G,df_Rmin_G = read_vdw_grad(dHere+'/gas/dE_dvdw_gas.txt',atomType)
    ScaleG =  2.0*w_hvap*(H_vap_Sim/List_hvap[id] - 1.0)/List_hvap[id]

    egrad = (df_Emin_G*ScaleG) + (df_Emin_L*ScaleL) ### TODO fix this  

    df_Emin_LV,df_Rmin_LV = read_vdw_grad(dHere+'/liquid/dV_dvdw_liquid.txt')
    ScaleLV =  2.0*w_vol*(MolVol_l/List_volume[id] - 1.0)/(List_volume[id]*RT)
   
    egrad = (df_Emin_G*ScaleG) + (df_Emin_L*ScaleL) + (df_Emin_LV*ScaleLV)
    rgrad = (df_Rmin_G*ScaleG) + (df_Rmin_L*ScaleL) + (df_Rmin_LV*ScaleLV)  

    d_hvap = H_vap_Sim/List_hvap[id] - 1.0
    d_v = MolVol_l/List_volume[id] - 1.0
    obj_func = w_vol*d_v*d_v + w_hvap*d_hvap*d_hvap
        
    print mol_id,obj_func,egrad,rgrad
    print >> fLog, mol_id,obj_func,egrad,rgrad
 
    os.chdir(workdir)

fLog.close()

