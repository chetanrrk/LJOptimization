import os,sys,pickle
import numpy as np
import pandas as pd

dat = pd.read_excel("/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/scripts/clusters_drude_for_analysis.xlsx")

def readAtomMap(atomMap):	
    fin=open(atomMap,"r")
    cgen=[]
    drude=[]
    for line in fin:
	token=line.rstrip().split()
	if token[0]!="CGenFF":
	    cgen.append(token[0])
	    drude.append(token[1])
    fin.close()
    return [cgen,drude]

def readGAAMPFF(ff):
    fin=open(ff,"r")
    gpTypes=[]
    for line in fin:
	token=line.rstrip().split()
	try:
	    if token[0]=="ATOM":gpTypes.append(token[2]) 
	except IndexError: continue
    fin.close()	
    return gpTypes

def getClusterAssignment(drudeType):
    cluster=-1
    clusters=[]
    for c in dat.columns:
	for d in dat[c]:
	    if not pd.isnull(d):
		if drudeType ==d: return c
    return -1

def getDrudeToGaampMap(mols,drudeTypes):
    gpTypes=[]	
    for m in mols:
        atomMap = "/lcrc/project/Drude/chetan/OLD_GAAMP/reparameterize/"+m+"/cgenff/cgenff_drude_new2_fix.txt"
        gaampFF = "/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/"+m+"/drude-result/box/ff_avg.str"
        if os.path.isfile(atomMap):
           if os.path.isfile(gaampFF):
                gaampTypes= readGAAMPFF(gaampFF)
                cgen,drude = readAtomMap(atomMap)
                #for a in drude:
                for i in range(len(drude)):
                    a = drude[i]
                    gp = gaampTypes[i]
                    if a in drudeTypes:
                        gpTypes.append(gp)
    return gpTypes

"""
TODO:
1) Read GAAMP FF get GAAMP atom types
2) figure out cluster assignment for each atom type
3) fix atom type based on cluster assignment
4) write the initial and bounds for each LJ parameter
"""
if __name__=="__main__":
    mols = pickle.load(open("hcAromats.p","r"))
    drudeTypes=[]
    for m in mols: 
	atomMap = "/lcrc/project/Drude/chetan/OLD_GAAMP/reparameterize/"+m+"/cgenff/cgenff_drude_new2_fix.txt"
	gaampFF = "/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/"+m+"/drude-result/box/ff_avg.str"
	if os.path.isfile(atomMap):
	   if os.path.isfile(gaampFF):
		gaampTypes= readGAAMPFF(gaampFF)
	   	cgen,drude = readAtomMap(atomMap)
		#for a in drude:
		for i in range(len(drude)):
		    a = drude[i]
		    gp = gaampTypes[i]
		    if a not in drudeTypes: 
			drudeTypes.append(a)
	   else: print m,"=========gaamp prob========="	
	else: print m,"=========cgenff prob========="

    	  
    #print drudeTypes
    clusters={}
    for a in drudeTypes:
        cluster=getClusterAssignment(a)
        if cluster==-1: print "=========cluster prob========="
	else: 
	    if not clusters.has_key(cluster): clusters[cluster]=[a]
	    else: clusters[cluster].append(a)	

    """
    fout=open("atomtypes_to_optimize.txt","w")	
    for c in clusters:
	print>>fout,c,
	for a in clusters[c]:
	    print>>fout,a,
	fout.write("\n")
    fout.close()   
    """
   
    	 
    fout=open("gaamp_atomtypes_to_optimize.txt","w") 	 
    for c in clusters:
	drudeTypes= clusters[c]
	#print drudeTypes
    	identicalAtomTypes = getDrudeToGaampMap(mols,drudeTypes)
	uniqTypes=[]
	for ats in identicalAtomTypes:
	    if not ats in uniqTypes:uniqTypes.append(ats)
	print uniqTypes
	for a in uniqTypes: #identicalAtomTypes:
	     print>>fout,a,
	fout.write("\n")
    fout.close()		


