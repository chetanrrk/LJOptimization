import os
import sys
import subprocess
import string
import shutil
import numpy as np
import multiprocessing as mp

def partitionMols(nParts):
    mols=[m.rstrip().split()[0] for m in open("./exp_set.txt","r")]	
    partition=[]	  
    nMolsInAPart = int(np.ceil(len(mols)*1.0/nParts))
    print nMolsInAPart
    start=0
    end=nMolsInAPart 
    for i in range(nParts):
	partition.append(mols[start:end])
	start=end
	end = end+nMolsInAPart
    return partition

def writeMinInp(T_Sim, L_BOX):
    line="""
# INPUT
coordinates             drude-box.pdb 
structure               mol-drude-liq.xpsf
parameters              ff.str 
paraTypeCharmm          on

# INITIAL CONDITIONS

temperature     $TEMP

# output params
outputname      namd-npt-min.out
binaryoutput    no

outputEnergies  100


# PME

PME                     yes
PMETolerance            10e-6
PMEInterpOrder          6
PMEGridSpacing          1.0


# Periodic Boundary Conditions
cellBasisVector1    $L_BOX 0 0
cellBasisVector2    0 $L_BOX 0
cellBasisVector3    0 0 $L_BOX
cellOrigin          0.0 0.0 0.0

wrapAll                 on

# drude
drude           on
drudeTemp       1.0
drudeBondLen    0.2
drudeHardWall   yes

# langevin thermostat
langevinPiston        on
langevinPistonTarget  1.01325      ;# pressure in bar -> 1 atm
langevinPistonPeriod  200.         ;# oscillation period around 200 fs
langevinPistonDecay   100.         ;# oscillation decay time of 100 fs
langevinPistonTemp    $TEMP        ;# coupled to heat bath

langevin        on
langevinTemp    $TEMP
langevinDamping 2.0 


StrainRate              0.0 0.0 0.0
useGroupPressure        yes        ;# This is needed when using SHAKE
useflexiblecell         no
splitpatch              hydrogen
hgroupcutoff            2.8
stepspercycle           20
margin                  0.0

switching               on
switchdist              10.0
cutoff                  12.0
pairlistdist            16.0
LJCorrection            on

exclude                 scaled1-4
1-4scaling              1.0

# RESPA PROPAGATOR

timestep                1.0 #fs
#fullElectFrequency      2
#nonbondedFreq           2


# SHAKE
rigidbonds              all

# COM
commotion               no

minimize 10000

"""
    line=line.replace("$TEMP", T_Sim)
    line=line.replace("$L_BOX", L_BOX)
    fOut=open('namd_NPT_min','w')
    fOut.write(line)
    fOut.close()
    

def writeMDInp(T_Sim, L_BOX):
    line="""
# INPUT
coordinates             namd-npt-min.out.coor
structure               mol-drude-liq.xpsf
parameters              ff.str 
paraTypeCharmm          on

# INITIAL CONDITIONS

temperature     $TEMP

# output params
outputname      namd-npt.out
binaryoutput    no

outputEnergies  100

DCDfreq        100
DCDfile        npt-box.dcd

# PME

PME                     yes
PMETolerance            10e-6
PMEInterpOrder          6
PMEGridSpacing          1.0


# Periodic Boundary Conditions
cellBasisVector1    $L_BOX 0 0
cellBasisVector2    0 $L_BOX 0
cellBasisVector3    0 0 $L_BOX
cellOrigin          0.0 0.0 0.0

wrapAll                 on

# drude
drude           on
drudeTemp       1.0
drudeBondLen    0.2
drudeHardWall   yes


# langevin thermostat
langevinPiston        on
langevinPistonTarget  1.01325      ;# pressure in bar -> 1 atm
langevinPistonPeriod  200.         ;# oscillation period around 200 fs
langevinPistonDecay   100.         ;# oscillation decay time of 100 fs
langevinPistonTemp    $TEMP        ;# coupled to heat bath

langevin        on
langevinTemp    $TEMP
langevinDamping 2.0 


StrainRate              0.0 0.0 0.0
useGroupPressure        yes        ;# This is needed when using SHAKE
useflexiblecell         no
splitpatch              hydrogen
hgroupcutoff            2.8
stepspercycle           20
margin                  0.0

switching               on
switchdist              10.0
cutoff                  12.0
pairlistdist            16.0
LJCorrection            on

exclude                 scaled1-4
1-4scaling              1.0

# RESPA PROPAGATOR

timestep                0.5 # 1.0 #fs it's 1fs by default..change it to 1fs
#fullElectFrequency      2
#nonbondedFreq           2


# SHAKE
rigidbonds              all

# COM
commotion               no

run 1000000 ; #500000

"""
    line=line.replace("$TEMP", T_Sim)
    line=line.replace("$L_BOX", L_BOX)
    fOut=open('namd_NPT','w')
    fOut.write(line)
    fOut.close()


def writeMinJobs(mol_ids):
    fout=open("min_jobs.sh","w")
    print>>fout,"export cDir=/lcrc/project/Drude/chetan/software/NAMD_2.13_Linux-x86_64-TCP/\n"
    i=0
    for m in mol_ids:
    	print>>fout,"$cDir/charmrun $cDir/namd2 +p 16 ++nodelist ../namd2.nodelist"+str(i)+" ++mpiexec ++remote-shell mpiexec "+m+"/namd_NPT_min > "+m+"/log_min.txt &"
	i+=1
    print>>fout,"wait"	
    fout.close()
    os.system("chmod +x min_jobs.sh")	  
	
def runMin():
    os.system("./min_jobs.sh")

def writeMDJobs(mol_ids):
    fout=open("md_jobs.sh","w")
    print>>fout,"export cDir=/lcrc/project/Drude/chetan/software/NAMD_2.13_Linux-x86_64-TCP/\n"
    i=0
    for m in mol_ids:
        print>>fout,"$cDir/charmrun $cDir/namd2 +p 16 ++nodelist ../namd2.nodelist"+str(i)+" ++mpiexec ++remote-shell mpiexec "+m+"/namd_NPT > "+m+"/log.txt &"
	i+=1
    print>>fout,"wait"    
    fout.close()        
    os.system("chmod +x md_jobs.sh")     

def runMD():
    os.system("./md_jobs.sh")

def extractEnergies(mol_ids):
    """
    Number of energy/volume steps and dcd print out should exacty match!!  
    """	
    fout = open("getEnergies.sh","w")
    for m in mol_ids:
    	print>>fout,'grep "^ENERGY:" '+m+'/log.txt | awk '+"'{print $14}' | sed -e '1,1d' > "+m+"/E_namd_dcd.txt &"	
	print>>fout,'grep "^ENERGY:" '+m+'/log.txt | awk '+"'{print $19}' | sed -e '1,1d' > "+m+"/V_namd_dcd.txt &"
    fout.close()
    os.system("chmod +x getEnergies.sh; ./getEnergies.sh")

def get_energy_and_vol(mol_id):

    fin = open(mol_id[0]+"/log.txt","r")
    fout_ene = open(mol_id[0]+"/E_namd_dcd.txt", "w")
    fout_vol = open(mol_id[0]+"/V_namd_dcd.txt", "w")
    flag = 0 ### flag for "ETITLE:"	
    for line in fin:
	token = line.rstrip().split()
	try:
		if token[0] == "ETITLE:":
			flag = 1
		elif flag:
			if token[0] == "ENERGY:" and int(token[1]) > 0 :
				print >> fout_ene, token[13]
				print >> fout_vol, token[18] 		
	except IndexError:
		continue
    fin.close()
    fout_ene.close()
    fout_vol.close()

def extractEnergiesMP(mol_ids):
    """
    Number of energy/volume steps and dcd print out should exacty match!!  
    Applying python's MP to avoid problems when making subprocesses	
    """
    procs = len(mol_ids)
    pool = mp.Pool(processes=procs) # creates as many process as the mols

    mol_list = [ [m] for m in mol_ids]
    result = pool.map(get_energy_and_vol, mol_list)
    pool.close()
    pool.join()
 

"""""
def writeGrad(mol_ids,T_Sims,L_BOXES):
    fout=open("liq_grads.sh","w")
    print>>fout,"export liqG=/lcrc/project/Drude/chetan/workspace/ljopt-independant-code/code_opt_lj_grad_liquid/\n"
    for i in range(len(mol_ids)):
	m=mol_ids[i]
	print>>fout,"mpiexec -np 16 $liqG/grad_liquid_e_v_new "+m+"/ff.str "+m+"/mol-drude-liq.xpsf "+m+"/npt-box.dcd 10.0 12.0 100 "+str(T_Sims[i])+" "+str(L_BOXES[i])+"> "+m+"/log-grad.txt &"
    print>>fout,"wait"
    fout.close()
    os.system("chmod +x liq_grads.sh")
"""""

def computeGrad(mol_ids,T_Sims,L_BOXES):
    for i in range(len(mol_ids)):
	m=mol_ids[i]
	os.chdir(m)
	liqG = "/lcrc/project/Drude/chetan/workspace/ljopt-independant-code/code_opt_lj_grad_liquid/"
	command="mpiexec -np 16 "+liqG+"grad_liquid_e_v_new ff.str mol-drude-liq.xpsf npt-box.dcd 10.0 12.0 100 "+str(T_Sims[i])+" "+str(L_BOXES[i])+" > log-grad.txt"
	retcode = subprocess.call(command,shell=True)
	if retcode>0:
		raise Exception("Problem computing gradients for",m)
	os.chdir("../")

"""
def computeGrad():
    os.system("./liq_grads.sh")
"""

def getTemp(idMol):
    fin = open('/lcrc/project/Drude/chetan/workspace/allparameters-correct.dat','r')
    for line in fin:
        token = line.rstrip().split()
        if token[2]==idMol:
                T_Sim = token[5]
                break
    fin.close()
    return T_Sim

def query_box_size(szName):
    f=open(szName,"r")
    line = f.readline()
    line = f.readline()
    line = f.readline()
    szBuff = line.split()
    L_Box = string.atof(szBuff[1])
    return L_Box


###LiqMaster is called by the main python script
if __name__=="__main__":
	###reading this partition of the molecules
	jobNum = sys.argv[1]
	cwd = sys.argv[2] 
	boxPath1 = '/homes/huanglei/restored/huanglei/ff/LJ/compound/'
	boxPath2 = '/ini_ff/ini_box/namd-npt.out.xsc'
        os.system("cp "+cwd+"/mols"+jobNum+".txt .")
	mol_ids = [m.rstrip().split()[0] for m in open("mols"+jobNum+".txt","r")]
	T_Sims = [getTemp(m) for m in mol_ids]
	box_lens= [str(query_box_size(boxPath1+m+boxPath2)) for m in mol_ids]
	for i in range(len(mol_ids)):
		m = mol_ids[i]
		if not os.path.isdir(m):
			os.mkdir(m)
			os.chdir(m)
			os.system("cp "+cwd+'/'+m+'/liquid/* .')	
			writeMinInp(T_Sims[i], box_lens[i])
			writeMDInp(T_Sims[i], box_lens[i])
			os.chdir("../")
	writeMinJobs(mol_ids)
	writeMDJobs(mol_ids)
	runMin()
	runMD()

        extractEnergiesMP(mol_ids)

	computeGrad(mol_ids,T_Sims,box_lens)
	
	for m in mol_ids:
		os.system("cp "+m+"/dE* "+cwd+"/"+m+"/liquid/")
		os.system("cp "+m+"/dV* "+cwd+"/"+m+"/liquid/")
		os.system("cp "+m+"/*_namd_dcd.txt "+cwd+"/"+m+"/liquid/")
		os.system("cp "+m+"/log-grad.txt "+cwd+"/"+m+"/liquid/")


