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

ferr=open(workdir+"/errors_dump.txt","w")

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
vdw_list_map = {}
df_dEmin_List = []
df_dRmin_List = []


"""
Maintains a list of halogens to use their specific grad calculation code
"""
halos = [m.rstrip().split()[0] for m in open("../halos.txt", "r")]

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

    eqList=[]
    fin=open(".." + path_sep + "equivalent.txt","r")
    for line in fin:
        token=line.rstrip().split()
        #if len(set(token).intersection(set(vdw_list)))>0:
        #        vdw_list_map.append(token)
	eqList.append(token)
    fin.close()

    for i in range(len(vdw_list)):
        for v2 in eqList:
                if vdw_list[i] in v2:
                        vdw_list_map[i] = v2       

    return n_vdw_type


def Accumulate_Grad_dE(nRec,vdw_Names_list,df_Emin_list,df_Rmin_list,Scale):
    global df_dEmin_List,df_dRmin_List
    for i in range(0,nRec):
        Idx = query_vdw_type(vdw_Names_list[i])
        if(Idx >= 0):
            df_dEmin_List[Idx] = df_dEmin_List[Idx] + df_Emin_list[i]*Scale
            df_dRmin_List[Idx] = df_dRmin_List[Idx] + df_Rmin_list[i]*Scale

def my_mkdir(newdir):
    if(not os.path.exists(newdir)) :
        os.mkdir(newdir)


def get_output_shell(cmd):
    retcode, result = commands.getstatusoutput(cmd)
    if retcode != 0:
        print 'Error in running: ' + cmd
        exit(1)
    else :
#        print "get_output_shell: " + result
        return result

def Wait_All_Job_Done():
    while(1):
        Done = 1
        for jobid in Job_List :
            retcode, result = commands.getstatusoutput('squeue --job ' + jobid)
            if retcode == 0: 
		try:
			if result.split()[12]=="PD" or result.split()[12]=="R":		
                		Done = 0
                		break
		except IndexError:continue
	    else: print "problem checking job status of",jobid

        if(Done):
            break
        else:
            time.sleep(15)

def query_box_size(szName):
    f=open(szName,"r")
    line = f.readline()
    line = f.readline()
    line = f.readline()
    szBuff = line.split()
    L_Box = string.atof(szBuff[1])
    return L_Box

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
    fin = open('allparameters-correct.dat','r')
    for line in fin:
	token = line.rstrip().split()
	if token[2]==molecule:
    		mass = float(token[3])
    		density = float(token[4])
    		T_Sim = float(token[5])
    		hvap = float(token[6])
    		dG_sfe = float(token[7])
    fin.close()

def writeGasInp(mol_id,temp):
    fout = open("gas-prep.in","w")
    print>>fout,"molName = "+mol_id+"  ### residue name"
    print>>fout,"liq_box = drude-box.pdb   ### name of the liq box pdb"
    print>>fout,"xpsf = mol-drude-gas.xpsf     ### name of the xpsf"
    print>>fout,"ffStr = ff.str    ### name of the ff.str\n"
    print>>fout,"temp = "+temp+"     ### simulation temp"	
    print>>fout,"mdSteps = 1000000 ### num. of MD Steps"
    print>>fout,"timeStep = 1.0 ### in fs\n"
    print>>fout,"nDraws = 250   ### sets of mols randomly drawn"
    print>>fout,"molPicks = 2  ### num. of mols in a set"
    print>>fout,"precision = 0.01 ### precision*KT is the stopping criteria"
    print>>fout,"ensembleAvg = False ### boltz. weighted ensemble average used for stdErr\n"
    print>>fout,"nodes = 1 ## number of nodes used to run"
    print>>fout,"tasks_per_node = 8 ### number of processor per node"
    print>>fout,"runTime = 7:00:00 ## time allowed for single gas MD to run"
    fout.close()


def createinput_liquid_job(T_Sim,mol_id,L_BOX):
    line="""#!/bin/bash
#SBATCH --job-name="""+mol_id+"""_liq
#SBATCH --account=Drude
#SBATCH --partition=sball
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=7:00:00
#SBATCH -o liq.out
#SBATCH -e liq.err

export SCRATCH=/scratch/$SLURM_JOB_ID/liquid
export CWD=`pwd`
export cDir=/lcrc/project/Drude/chetan/software/NAMD_2.13_Linux-x86_64-TCP/
#export OMP_NUM_THREADS=1

mkdir -p $SCRATCH
cp * $SCRATCH
cd $SCRATCH

###minimizing first for some steps
$cDir/charmrun $cDir/namd2 +p 16 ++mpiexec ++remote-shell mpiexec namd_NPT_min > log_min.txt


###running md reading minimized frame
$cDir/charmrun $cDir/namd2 +p 16 ++mpiexec ++remote-shell mpiexec namd_NPT > log.txt


# Number of energy/volume steps and dcd print out should exacty match!!  
grep "^ENERGY:" log.txt | awk '{print $14}' | sed -e '1,1d' > E_namd_dcd.txt
grep "^ENERGY:" log.txt | awk '{print $19}' | sed -e '1,1d' > V_namd_dcd.txt

###skipping first 100 md frames as equilibration
mpiexec -np 16 /lcrc/project/Drude/chetan/workspace/ljopt-independant-code/code_opt_lj_grad_liquid/grad_liquid_e_v_new ff.str mol-liq.xpsf npt-box.dcd 10.0 12.0 100 $TEMP $LBOX> log-grad.txt

cp E_namd_dcd.txt $CWD
cp V_namd_dcd.txt $CWD
cp dE_dvdw_liquid.txt $CWD
cp dV_dvdw_liquid.txt $CWD
cp log-grad.txt $CWD

#tar -czvf npt-box.dcd.tgz npt-box.dcd 
#rm npt-box.dcd
#tar -czvf log.txt.tgz log.txt
#rm log.txt
"""
    line=line.replace("$TEMP", T_Sim)
    line=line.replace("$LBOX",L_BOX)	 
    return(line)


def createinput_namd_liquid(T_Sim, L_BOX,mol_id):
    line="""
# INPUT
coordinates             namd-npt-min.out.coor
structure               mol-liq.xpsf
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

timestep                2.0 #fs
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
    return(line)


def createinput_namd_liquid_min(T_Sim, L_BOX):
    line="""
# INPUT
coordinates             ini-box.pdb
structure               mol-liq.xpsf
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
    return(line)

def createinput_sfe_job(org_dir):
    line="""
#!/bin/bash
#PBS -l nodes=1:ppn=16
##PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
####PBS -q roux
#PBS -j oe
#PBS -N sfe

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodes.txt

cp * $TMPDIR
cd $TMPDIR

/homes/huanglei/prog/NAMD/NAMD_2.8_Linux-x86_64-multicore/namd2 +p16 namd_NPT > log-namd.txt
#~/prog/NAMD/NAMD_2.8_Linux-x86_64-multicore/namd2 +p8 namd_NPT > log-namd.txt

/homes/huanglei/prog/grad_vdw_sfe ff.str mol.xpsf npt-box.dcd 100 > log-grad-vdw-sfe.txt

/homes/huanglei/prog/de_vdw_sfe_pert ff.str mol.xpsf $org_dir/npt-box.dcd > E_int_0_B.txt
/homes/huanglei/prog/de_vdw_sfe_pert $org_dir/ff.str mol.xpsf npt-box.dcd > E_int_1_A.txt
/homes/huanglei/prog/de_vdw_sfe_pert ff.str mol.xpsf npt-box.dcd > E_int_1_B.txt

/homes/huanglei/prog/minus E_int_0_B.txt $org_dir/../E_int_0_A.txt 0.0 > state_0.txt
/homes/huanglei/prog/minus E_int_1_B.txt E_int_1_A.txt 1.0 > state_1.txt

cat state_0.txt > wham.in
cat state_1.txt >> wham.in
/homes/huanglei/prog/wham_fep wham.in wham.out
cat wham.out | awk '{print $2}' > ddG_sfe.txt

cp * $PBS_O_WORKDIR/

"""
    line=line.replace("$org_dir", org_dir)
    return(line)

def createinput_namd_sfe(L_BOX,MolSize):
    line="""
# INPUT
coordinates             ini-box.pdb 
structure               mol.xpsf
parameters              ff.str 
paraTypeCharmm          on

temperature     298.150000

# output params
outputname      namd-npt.out
binaryoutput    no

outputEnergies  1000

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
cellOrigin          0.0 0.0 0.1

wrapAll                 on


# langevin thermostat
langevinPiston        on
langevinPistonTarget  1.01325      ;# pressure in bar -> 1 atm
langevinPistonPeriod  200.         ;# oscillation period around 200 fs
langevinPistonDecay   100.         ;# oscillation decay time of 100 fs
langevinPistonTemp    298.150000   ;# coupled to heat bath

langevin        on
langevinTemp    298.150000
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
pairlistdist            15.0
LJCorrection            on

exclude                 scaled1-4
1-4scaling              1.0

timestep                2.0 #fs

# SHAKE
rigidbonds              all

# COM
commotion               no

tclForces           on

tclForcesScript {

  global kforce

  set natom   $MolSize
  set kforce  10.0 ;  # kcal/mol/AA Force constant for the restraint

# Output parameters

   set selecom {}
   for {set i 1 } {$i <= $natom } {incr i 1} {
        set id $i
        lappend selecom $id
        addatom $id
   }

## group indices with addgroup

  set grp [addgroup $selecom]

# Set reference value of com coordinates

  set xc0   0.000
  set yc0   0.000
  set zc0   0.000

   proc calcforces {} {
    global grp xc0 yc0 zc0
    global kforce
    global freq fout ofile

    loadcoords coords

    # DO COM RESTRAINT

    # get com coords
    set rc $coords($grp)
    set xc [lindex $rc 0]
    set yc [lindex $rc 1]
    set zc [lindex $rc 2]

    # calculate forces on com
    set diffx [expr {$xc0 - $xc}]
    set diffy [expr {$yc0 - $yc}]
    set diffz [expr {$zc0 - $zc}]
    set f {}
    set fx [expr {$kforce*$diffx}]
    set fy [expr {$kforce*$diffy}]
    set fz [expr {$kforce*$diffz}]
    lappend f $fx $fy $fz

    # distribute force proportionally to atoms
    # in group grp (NAMD does this automatically)
    addforce $grp $f

    # add energy
    addenergy [expr {0.5*$kforce*$diffx*$diffx}]
    addenergy [expr {0.5*$kforce*$diffy*$diffy}]
    addenergy [expr {0.5*$kforce*$diffz*$diffz}]

    return
  }
}


minimize 100
#run 2000000
run 1500000
"""
    line=line.replace("$L_BOX", L_BOX)
    line=line.replace("$MolSize", MolSize)
    return(line)


def getVDWParam(vdwToFit):
    fin = open(vdwToFit)
    params={}
    for line in fin:
        token=line.rstrip().split()
        if token[0]!="vdw_type":
           params[token[0]] = [t for t in token[1:]]
    return params

def isAtomTypeMatch(typeList,atomTypeToMatch,optAtomType):
    matchTypes=[]
    for ts in typeList:
	if optAtomType in ts:
	    matchTypes = ts
    return atomTypeToMatch in matchTypes
		

def updateVdw(ffStr,vdwParams):
    fin=open("../../equivalent.txt","r") ### clusters gaamp types for all molecules in a class used 
    atomTypes=[]
    for line in fin:
	token=line.rstrip().split()
	atomTypes.append(token)	
    fin.close() 

    fin = open(ffStr,"r")
    fout = open("ff_new.str","w")
    flag=0
    for line in fin:
        token = line.rstrip().split()
        try:
           #if token[1]=="EMIN":
           if token[0]=="NONBONDED":
              flag=1
           if flag==1:
              for k in vdwParams:
		   if isAtomTypeMatch(atomTypes,token[0],k):
                        line = line.replace(token[2],vdwParams[k][0])
                        line = line.replace(token[3],vdwParams[k][1])
                        #line = line.replace(token[5],'{0:.4f}'.format(float(vdwParams[k][0])/2.0))
                        #line = line.replace(token[6],vdwParams[k][1])
        
           if token[0]!="ioformat":
              print>>fout, line.rstrip()
        except IndexError:continue
    fin.close()
    fout.close()
    os.system("mv ff_new.str ff.str")


def submitNow(r=32,q=100):
	"""
	r= #running jobs
	q= #queued jobs
	"""
	time.sleep(15) ### short delay to let other jobs populate the queue
        i=0
        rJobs=0 ###running jobs
        pdJobs=0 ###pending jobs
        retcode, result = commands.getstatusoutput('squeue -u ac.rupakhetic')
        if not retcode:
                result= result.split()
                lim= len(result)/8
                for j in range(lim):
                        status = result[i:(i+8)][4]
                        if status=="R":rJobs+=1
                        elif status=="PD":pdJobs+=1
                        i = i+8

        #if rJobs+pdJobs<q and rJobs<r: return True
        if rJobs+pdJobs<r: return True
        else: return False


def getTemp(idMol):
    fin = open('allparameters-correct.dat','r')
    for line in fin:
        token = line.rstrip().split()
        if token[2]==idMol:
                T_Sim = token[5]
		break
    fin.close()
    return T_Sim	


def writeGasMaster(mol_id,T_Sim):
    fout = open("gas_master.sh","w")	
    print>>fout,"#!/bin/bash"
    print>>fout,"#SBATCH --job-name="+mol_id+"_gas_master"
    print>>fout,"#SBATCH --account=Drude"
    print>>fout,"#SBATCH --partition=sball"
    print>>fout,"#SBATCH --nodes=1"
    print>>fout,"#SBATCH --ntasks-per-node=16"
    print>>fout,"#SBATCH --time=7:00:00"
    print>>fout,"#SBATCH -o gas_master.out"	
    print>>fout,"#SBATCH -e gas_master.err"	
    print>>fout,"export SCRATCH=/lcrc/globalscratch/rupakhetic/$SLURM_JOB_ID/gas"
    if mol_id in halos:
       print>>fout,"export codeDir=../../halo_code/"
    else:
       print>>fout,"export codeDir=../../non_halo_code/"
    print>>fout,"export CWD=`pwd`"
    print>>fout,"mkdir -p $SCRATCH"
    print>>fout,"cp $CWD/* $SCRATCH"
    print>>fout,"cp $codeDir/GasJob-MPI.py $SCRATCH"
    print>>fout,"cd $SCRATCH\n"
    print>>fout,"export NODES=$SLURM_JOB_NODELIST"
    print>>fout,"echo $NODES"
    print>>fout,"echo $SLURM_JOB_ID"
    print>>fout,"echo $SLURM_JOB_NODELIST"
    print>>fout,"python $codeDir/writeNodeList.py $NODES\n"
    print>>fout,"python $codeDir/GasJob-MPI.py gas-prep.in\n"
    print>>fout,"cp gas_done.txt $CWD"
    print>>fout,"cp gas_sub_error.txt $CWD"
    print>>fout,"cp submitted_jobs.txt $CWD"
    print>>fout,"cp stdErr.txt $CWD"
    print>>fout,"cp dE_dvdw_gas.txt $CWD"
    print>>fout,"cp data_E_gas.txt $CWD\n"
    print>>fout,"cp stats_E_gas.txt $CWD\n"
    #print>>fout,"nohup python ../../../GasJob.py gas-prep.in &"	
    fout.close()    


def partitionMols(nParts,mols_to_process):
    partition=[]
    try:
        nMolsInAPart = int(np.ceil(len(mols_to_process)*1.0/nParts))
    except ZeroDivisionError:
        return partition
    print nMolsInAPart
    start=0
    end=nMolsInAPart
    for i in range(nParts):
        partition.append(mols_to_process[start:end])
        start=end
        end = end+nMolsInAPart
    return partition


"""
mol_ids: list of ids being processed
T_sims: simulation temperature for a list of mols
nodes: number of nodes to use
"""
def writeLiquidMaster(jobNum,nodes=20):
    fout = open("liq_master"+str(jobNum)+".sh","w")	
    print>>fout,"#!/bin/bash"
    print>>fout,"#SBATCH --job-name="+str(jobNum)+"_liq_master"
    print>>fout,"#SBATCH --account=Drude"
    print>>fout,"#SBATCH --partition=sball"
    print>>fout,"#SBATCH --nodes="+str(nodes)
    print>>fout,"#SBATCH --ntasks-per-node=16"
    print>>fout,"#SBATCH --time=7:00:00"
    print>>fout,"#SBATCH -o liq_master"+str(jobNum)+".out"	
    print>>fout,"#SBATCH -e liq_master"+str(jobNum)+".err"	
    print>>fout,"export SCRATCH=/lcrc/globalscratch/rupakhetic/$SLURM_JOB_ID/liquid"
    print>>fout,"export CWD=`pwd`"
    print>>fout,"export codeDir=/lcrc/project/Drude/chetan/workspace/ljopt-tests/drude_lj_prod/optCode/"
    print>>fout,"mkdir -p $SCRATCH"
    print>>fout,"cp $codeDir/LiqJobMPI.py $SCRATCH"
    print>>fout,"cd $SCRATCH\n"
    print>>fout,"export NODES=$SLURM_JOB_NODELIST"
    print>>fout,"echo $NODES"
    print>>fout,"echo $SLURM_JOB_NODELIST"
    print>>fout,"echo $SLURM_JOB_ID"
    print>>fout,"python $codeDir/writeNodeList.py $NODES\n"
    print>>fout,"python LiqJobMPI.py "+str(jobNum)+" $CWD\n"
    fout.close()    


def submitGasMD(exclude):
    fin=open("../exp_set.txt","r")
    mols=[]
    for line in fin:
	mol_id = line.rstrip().split()[0]
	if mol_id not in exclude:
		mols.append(mol_id)
    fin.close()
    
    for mol_id in mols:
    	newdir = workdir + path_sep + mol_id
    	os.chdir(newdir)
	str_T_Sim = getTemp(mol_id)
	gasHelper(mol_id,str_T_Sim)

def getNodeInfo(jobid):
	retcode, result = commands.getstatusoutput('squeue -j '+jobid)
	if retcode==0: return result.split()[-1]
	else: 
		print "Error checking job id",jobid
		return -999
   
def gasHelper(mol_id,str_T_Sim):
    	os.chdir('gas')
    	#if(not os.path.isfile("dE_dvdw_gas.txt")):
	if not os.path.isfile("gas_done.txt"):
    		print "submitting gas job for "+mol_id
        	writeGasInp(mol_id,str_T_Sim)
		writeGasMaster(mol_id,T_Sim)
		#os.system("bash gas_master.sh")
	
		while (True):
		   if submitNow():
			retcode, result = commands.getstatusoutput('sbatch gas_master.sh')
			if retcode != 0:
		   	    print >>ferr,'Error in submitting job ' + mol_id + ' for gas simulation.'	
		   	    exit(1)
			else: 
			    
			    jobid = result.split()[3]
		    	    gasJobList.append(jobid)
			    node = getNodeInfo(jobid)	
			    fout=open("nodeInfo.txt","w")
			    print>>fout,node
			    fout.close()
			    break ### get out of the while loop	
		   else:time.sleep(15)
  		

def waitAllGasJobs():
    while(1):
        Done = 1
	"""
	In case of restart, this flag file will be handy
	"""
	for m in os.listdir("."):
	    if m.startswith("m_"):
		if not os.path.isfile(m+"/gas/gas_done.txt"):
		   Done = 0 
		   break		
	
        for jobid in gasJobList:
            retcode, result = commands.getstatusoutput('squeue --job ' + jobid)
            if retcode == 0: 
		try:
			if result[12]=="PD" or result[12]=="R": ### pending or running status		
                		Done = 0
                		break
		except IndexError:continue
	    else:
		print "Error in processing Gas job",jobid	
		continue ## job is done if here 
	
        if(Done):
            break
        else:
            time.sleep(15)


if (os.path.isfile('done.txt')):
    Done = 1
else:
    Done = 0


n_vdw_type = read_vdw_type()

nMol = 0
f=open("../exp_set.txt","r")

#submit all liquid jobs before submitting gas jobs
mols_to_process=[] ## stores only mols whose derivatives are to be computed
allMols=[]
for line in f:
    mol_id = line.rstrip().split()[0]
    allMols.append(mol_id)
    newdir = workdir + path_sep + mol_id
    my_mkdir(newdir)
    os.chdir(newdir)
   
    get_experimenal_values(mol_id)

    #print mol_id
    #print mass
    #print density 
    #print hvap

    volume = ATOM_MASS * mass / density
    str_T_Sim = '%f' %(T_Sim)

    List_T_Sim.append(T_Sim)
    List_hvap.append(hvap)
    List_volume.append(volume)
    List_sfe.append(dG_sfe)
    List_Mol.append(mol_id)

    if(Done == 0):
        my_mkdir('gas')   
        my_mkdir('liquid')
        my_mkdir('sfe')

        
        """	
        # preparing the files for simulations ... added LPH to the halogens pdbs in a new folder and using those as source dir
        src_dir_bulk1 = '/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/'+mol_id+'/drude-result/box' #'/lcrc/project/Drude/chetan/FinalDrudeBoxes/FinalDrudeBoxes/' + mol_id
        src_dir_mol1 = '/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/'+mol_id+'/drude-result/box'  #'/lcrc/project/Drude/chetan/FinalDrudeMol/' + mol_id

        src_dir_bulk2 = "/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/halos_LPH_code/halos/"+mol_id+"/drude_results_new/box" #'/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/'+mol_id+'/drude-result/box' 
        src_dir_mol2 = "/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/halos_LPH_code/halos/"+mol_id+"/drude_results_new/box" #'/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/'+mol_id+'/drude-result/box'  
        """

	psfs = '../../psfs'
        molnumber = mol_id[2:]

        os.system('cp ../../ff_to_opt/ff_new_'+mol_id+'.str ff.str')  ##updated ff from previous optimization


	vdwParams = getVDWParam("../vdw-param.txt")

	updateVdw('ff.str',vdwParams) # update vdw parameters with vdw to be fitted 

        """
	#  check halos first b/c of LPH stuff
	if os.path.isfile(src_dir_mol2 + '/mol-1.pdb'):
		src_dir_mol = src_dir_mol2
		src_dir_bulk = src_dir_bulk2
	else:
		src_dir_mol = src_dir_mol1
		src_dir_bulk = src_dir_bulk1
        """
        
        #os.system('cp ' + src_dir_mol + '/mol-1.pdb ./gas/mol.pdb')
        src_dir_bulk = '../../train_all_pdbs/'+mol_id+'-drude-box.pdb'

        os.system('cp ff.str gas/')
        os.system('cp ff.str liquid/')
        os.system('cp ff.str sfe/')

        os.system('cp ' + src_dir_bulk+' ./liquid/drude-box.pdb')
        os.system('cp ' + src_dir_bulk+' ./gas/drude-box.pdb')
        #os.system('cp ' + src_dir_bulk+'/mol-drude-liq.xpsf ./liquid/')
        #os.system('cp ' + src_dir_mol+'/mol-drude-gas.xpsf ./gas/')
        os.system('cp ' + psfs+'/'+mol_id+'-drude-liq.xpsf ./liquid/mol-drude-liq.xpsf') ### atom types are updated
        os.system('cp ' + psfs+'/'+mol_id+'-drude-gas.xpsf ./gas/mol-drude-gas.xpsf')	### atom types are updated
	
        if(not os.path.isfile("./liquid/dE_dvdw_liquid.txt")):
	    mols_to_process.append(mol_id) 	
    nMol = nMol + 1

f.close()

os.chdir(workdir)

"""
lets submit liq jobs using new schedule
"""
liqJobNum = 1 ### counter of created liquid jobs
nParts = len(mols_to_process) ### number of liq jobs to be created
partMols = partitionMols(nParts,mols_to_process)	 
boxPath1 = '/homes/huanglei/restored/huanglei/ff/LJ/compound/'
boxPath2 = '/ini_ff/ini_box/namd-npt.out.xsc'

if len(partMols) > 0:
    nodes = len(partMols[0]) ### can optimize this further
else:
    nodes = 0
#print "partMols"
#print partMols;

for mols in partMols:
	if len(mols)<1: continue ### mols are done ... no need to resubmit
	fout=open("mols"+str(liqJobNum)+".txt","w")
	for m in mols: print>>fout,m
	fout.close()
	while (True):
	    writeLiquidMaster(str(liqJobNum),nodes)
	    if submitNow():
	       print "submitting liquid master for "+str(liqJobNum)
	       retcode, result = commands.getstatusoutput('sbatch liq_master'+str(liqJobNum)+'.sh')
	       if retcode != 0:
	            print >>ferr,'Error in submitting job ' + str(liqJobNum) + ' for liquid simulation.'
	            exit(1)
	       else :
	            jobid = result.split()[3]
	            Job_List.append(jobid)
	            node = getNodeInfo(jobid)       
	            fout=open("liq_master_nodeInfo"+str(liqJobNum)+".txt","w")
	            print>>fout,node
	            fout.close()
	            break   ### break out of this while loop
	    else: time.sleep(15)
	liqJobNum+=1

"""
submit gas jobs in parallel for all molecules after submitting liquid jobs
"""
exclude=[] ### to exclude mols that have derivatives computed
for m in allMols:
    if os.path.isfile(m+"/gas/dE_dvdw_gas.txt"):
	exclude.append(m)
if Done==0: submitGasMD(exclude)
print "submitted all gas jobs"

os.chdir(workdir)

if(Done == 0):
    waitAllGasJobs() ### makes main thread wait until all gas jobs are done
    Wait_All_Job_Done()

# To process the data

os.chdir(workdir)
fLog = open("log.txt","w")

obj_func = 0.0
for id in range(0,nMol):
    mol_id = List_Mol[id]
    T_Sim = List_T_Sim[id]
    RT = T_Sim * K_BOTZMAN
    newdir = workdir + path_sep + mol_id
    os.chdir(newdir)

    boxName = "liquid/drude-box.pdb"
    command1 = "tail -n 3 "+boxName+" | grep ""ATOM"" | awk '{print $5}'"
    process = subprocess.Popen(command1,stdout=subprocess.PIPE, shell=True)
    try:
        N_Mol_l = int(process.communicate()[0])
    except ValueError:
        process = subprocess.Popen(command1,stdout=subprocess.PIPE, shell=True)
        N_Mol_l = int(process.communicate()[0].split('\n')[0])
    #print os.getcwd() 
    if not os.path.isfile("enthalpy-vol-summary.txt"):
    	os.system("python ../../ComputeDeltaH_NAMD.py "+str(T_Sim))
    ###read the H and Vol txt file here but wait until done
    while True:	
	if os.path.isfile("enthalpy-vol-summary.txt"):
    	   fin = open("enthalpy-vol-summary.txt","r")	
	   for line in fin:
		token = line.rstrip().split()
		if token[0]=="(DeltaH)average:":H_vap_Sim = float(token[1])
		elif token[0]=="(Volume)average,std:":MolVol_l = float(token[1])/N_Mol_l
	   break		   				
	else:
	   time.sleep(180) 

    ### Gradients for the gas phase are computed from Gas job script
   	
    dHere = os.getcwd()
    #print dHere
    nRec,vdw_Names_list,df_Emin_list,df_Rmin_list = read_vdw_grad(dHere+'/liquid/dE_dvdw_liquid.txt')

    #print id, List_hvap[id]	
    Scale = -2.0*w_hvap*(H_vap_Sim/List_hvap[id] - 1.0)/List_hvap[id]
    Accumulate_Grad_dE(nRec,vdw_Names_list,df_Emin_list,df_Rmin_list,Scale)

    nRec,vdw_Names_list,df_Emin_list,df_Rmin_list = read_vdw_grad(dHere+'/gas/dE_dvdw_gas.txt')
    Scale =  2.0*w_hvap*(H_vap_Sim/List_hvap[id] - 1.0)/List_hvap[id]
    Accumulate_Grad_dE(nRec,vdw_Names_list,df_Emin_list,df_Rmin_list,Scale)

    nRec,vdw_Names_list,df_Emin_list,df_Rmin_list = read_vdw_grad(dHere+'/liquid/dV_dvdw_liquid.txt')
    Scale =  2.0*w_vol*(MolVol_l/List_volume[id] - 1.0)/(List_volume[id]*RT)
    Accumulate_Grad_dE(nRec,vdw_Names_list,df_Emin_list,df_Rmin_list,Scale)
    
    d_hvap = H_vap_Sim/List_hvap[id] - 1.0
    d_v = MolVol_l/List_volume[id] - 1.0
    obj_func = obj_func + w_vol*d_v*d_v + w_hvap*d_hvap*d_hvap
        
    #print mol_id,obj_func
	
    szBuff = '%-8s V0 = %8.3f  V = %8.3f   Hvap_0 = %6.3f Hvap = %6.3f\n' % (mol_id, List_volume[id], MolVol_l, List_hvap[id],H_vap_Sim)
    fLog.write(szBuff)

fLog.close()


os.chdir(workdir)

fOut = open("obj_func.txt","w")
szBuff = '%.8e\n' % (obj_func)
fOut.write(szBuff)
fOut.close()


fOut = open("grad_list.txt","w")
for i in range(0,n_vdw_type):
    szBuff = '%-8s  %20.14e %20.14e\n' % (vdw_list[i],df_dEmin_List[i],df_dRmin_List[i])
    fOut.write(szBuff)
fOut.close()

os.system('echo "Done" > done.txt')


