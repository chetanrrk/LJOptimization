import os
import sys
import pickle
import time
import string
import commands,subprocess
import re
	
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
T_Sim = '0.0'
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

    fin=open(".." + path_sep + "gaamp_atomtypes_to_optimize.txt","r")
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
            retcode, result = commands.getstatusoutput('qstat ' + jobid)
            if retcode == 0: # job is still running
                Done = 0
                break
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

def writeGasInp(mol_id,temp):
    fout = open("gas-prep.in","w")
    print>>fout,"molName = "+mol_id+"  ### residue name"
    print>>fout,"liq_box = drude-box.pdb   ### name of the liq box pdb"
    print>>fout,"xpsf = mol-drude-gas.xpsf     ### name of the xpsf"
    print>>fout,"ffStr = ff.str    ### name of the ff.str\n"
    print>>fout,"temp = "+temp+"     ### simulation temp"	
    print>>fout,"mdSteps = 1000000 ### num. of MD Steps"
    print>>fout,"timeStep = 1.0 ### in fs\n"
    print>>fout,"nDraws = 25   ### sets of mols randomly drawn"
    print>>fout,"molPicks = 20  ### num. of mols in a set"
    print>>fout,"precision = 0.01 ### precision*KT is the stopping criteria"
    print>>fout,"ensembleAvg = False ### boltz. weighted ensemble average used for stdErr\n"
    print>>fout,"nodes = 1 ## number of nodes used to run"
    print>>fout,"tasks_per_node = 16 ### number of processor per node"
    print>>fout,"runTime = 00:10:00 ## time allowed for single gas MD to run"
    fout.close()


def createinput_liquid_job(T_Sim,mol_id,L_BOX):
    line="""
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=05:00:00
#PBS -j oe
#PBS -N """+mol_id+"""_liq
#PBS -A Drude

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodes.txt

source /blues/gpfs/home/software/spack-0.10.1/opt/spack/linux-centos7-x86_64/gcc-4.8.5/intel-parallel-studio-cluster.2018.1-egcacagjokllneqafnpkfnp4njzklpsk/parallel_studio_xe_2018/psxevars.sh


###minimizing first for some steps

mpirun -n 16 namd2 namd_NPT_min > log_min.txt

###running md reading minimized frame
mpirun -n 16 namd2 namd_NPT > log.txt

# Number of energy/volume steps and dcd print out should exacty match!!  
grep "^ENERGY:" log.txt | awk '{print $14}' | sed -e '1,1d' > E_namd_dcd.txt
grep "^ENERGY:" log.txt | awk '{print $19}' | sed -e '1,1d' > V_namd_dcd.txt

cp * $PBS_O_WORKDIR/
cd $PBS_O_WORKDIR/

###skipping first 100 md frames as equilibration
mpiexec -np 16 /lcrc/project/Drude/chetan/workspace/ljopt-independant-code/code_opt_lj_grad_liquid/grad_liquid_e_v_new ff.str mol-drude-liq.xpsf npt-box.dcd 10.0 12.0 100 $TEMP $LBOX> log-grad.txt

tar -czvf npt-box.dcd.tgz npt-box.dcd 
rm npt-box.dcd
tar -czvf log.txt.tgz log.txt
rm log.txt
"""
    line=line.replace("$TEMP", T_Sim)
    line=line.replace("$LBOX",L_BOX)	 
    return(line)


def createinput_namd_liquid(T_Sim, L_BOX,mol_id):
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
drudeTemp       0.1
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
# drude
#drudeDamping    5.0


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

run 2000000 ; #500000

"""
    line=line.replace("$TEMP", T_Sim)
    line=line.replace("$L_BOX", L_BOX)
    if mol_id=="m_211":
	line=line.replace("timestep                1.0 #fs", "timestep                0.5 #fs")	
    return(line)


def createinput_namd_liquid_min(T_Sim, L_BOX):
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
drudeTemp       0.1
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
# drude
#drudeDamping    5.0


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
cellOrigin          0.0 0.0 0.0

wrapAll                 on

# drude
#drude           on
#drudeTemp       1
#drudeBondLen    0.2
#drudeBondConst  40000

# langevin thermostat
langevinPiston        on
langevinPistonTarget  1.01325      ;# pressure in bar -> 1 atm
langevinPistonPeriod  200.         ;# oscillation period around 200 fs
langevinPistonDecay   100.         ;# oscillation decay time of 100 fs
langevinPistonTemp    298.150000   ;# coupled to heat bath

langevin        on
langevinTemp    298.150000
langevinDamping 2.0 
# drude
#drudeDamping    5.0


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
    fin=open("../../gaamp_atomtypes_to_optimize.txt","r") ### clusters gaamp-drude types for all molecules in a class used 
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

def submitNow():
    numJobs=0
    retcode, result = commands.getstatusoutput("qstat -u rupakhetic")
    result = result.strip().split()
    for j in result:
	if j.endswith("bmgt1."):numJobs+=1
    if numJobs < 120:return True  


def getTemp(idMol):
    fin = open('/lcrc/project/Drude/chetan/workspace/allparameters-correct.dat','r')
    for line in fin:
        token = line.rstrip().split()
        if token[2]==idMol:
                T_Sim = token[5]
		break
    fin.close()
    return T_Sim	


def writeGasMaster(mol_id,T_Sim):
    fout = open("gas_master.sh","w")
    """	
    print>>fout,"#!/bin/bash"
    print>>fout,"#PBS -l nodes=1:ppn=1"
    print>>fout,"#PBS -l walltime=05:00:00"
    print>>fout,"#PBS -j oe"
    print>>fout,"#PBS -N "+mol_id+"_gas_master"
    print>>fout,"#PBS -A Drude\n"
    print>>fout,"cd $PBS_O_WORKDIR"
    print>>fout,"cat $PBS_NODEFILE > nodes.txt\n"
    print>>fout,"python ../../../GasJob.py gas-prep.in"
    """
    print>>fout,"nohup python ../../../GasJob.py gas-prep.in &"	
    fout.close()    

def submitGasMD():
    fin=open("../exp_set.txt","r")
    mols=[]
    for line in fin:
	mol_id = line.rstrip().split()[0]
	mols.append(mol_id)
    fin.close()
    
    for mol_id in mols:
    	newdir = workdir + path_sep + mol_id
    	os.chdir(newdir)
	str_T_Sim = getTemp(mol_id)
	gasHelper(mol_id,str_T_Sim)
   
def gasHelper(mol_id,str_T_Sim):
    	os.chdir('gas')
    	#if(not os.path.isfile("dE_dvdw_gas.txt")):
	if not os.path.isfile("gas_done.txt"):
    		print "submitting gas job for "+mol_id
        	writeGasInp(mol_id,str_T_Sim)
		writeGasMaster(mol_id,T_Sim)
		os.system("bash gas_master.sh")
		os.system("ps -f >> gas_master_job_ids.txt")
		time.sleep(15) ### allows time to check the queue for each gas job
		"""
		while (True):
		   if submitNow():
			retcode, result = commands.getstatusoutput('qsub gas_master.sh')
			if retcode != 0:
		   	    print >>ferr,'Error in submitting job ' + mol_id + ' for gas simulation.'	
		   	    exit(1)
			else: 
		    	    jobid = string.split(result,'.')
		    	    jobid = jobid[0]
		    	    gasJobList.append(jobid)
			    break ### get out of the while loop	
		   else:time.sleep(15)
  		"""

def waitAllGasJobs():
    while(1):
        Done = 1
	for m in os.listdir("."):
	    if m.startswith("m_"):
		if not os.path.isfile(m+"/gas/gas_done.txt"):
		   Done = 0 
		   break		
	"""
        for jobid in gasJobList:
            retcode, result = commands.getstatusoutput('qstat ' + jobid)
            if retcode == 0: # job is still running
                Done = 0
                break
	"""
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

for line in f:
    #mol_id = f.readline()
    #mol_id = mol_id.replace('\n','')
    #if not mol_id: break
    mol_id = line.rstrip().split()[0]

    newdir = workdir + path_sep + mol_id
    my_mkdir(newdir)
    os.chdir(newdir)
   
    #if os.path.isfile("liquid/dE_dvdw_liquid.txt") and os.path.isfile("liquid/dV_dvdw_liquid.txt"):
    # 	continue	

 
    get_experimenal_values(mol_id)

    print mol_id
    print mass
    print density 
    print hvap

    volume = ATOM_MASS * mass / density
    str_T_Sim = '%f' %(T_Sim)

    List_T_Sim.append(T_Sim)
    List_hvap.append(hvap)
    List_volume.append(volume)
    List_sfe.append(dG_sfe)
    List_Mol.append(mol_id)

    if(Done == 0):
	###open up these comments
        my_mkdir('gas')   
        my_mkdir('liquid')
        my_mkdir('sfe')
        
        # preparing the files for simulations
        src_dir_bulk = '/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/'+mol_id+'/drude-result/box' #'/lcrc/project/Drude/chetan/FinalDrudeBoxes/FinalDrudeBoxes/' + mol_id
        src_dir_mol = '/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/'+mol_id+'/drude-result/box'  #'/lcrc/project/Drude/chetan/FinalDrudeMol/' + mol_id
        molnumber = mol_id[2:]

        os.system('cp ../../ff_to_opt/ff_new_'+mol_id+'.str ff.str')  ##updated ff from previous optimization
	#os.system('cp '+src_dir_bulk+'/ff_avg.str ff.str')


	vdwParams = getVDWParam("../vdw-param.txt")

	updateVdw('ff.str',vdwParams) # update vdw parameters with vdw to be fitted 

        os.system('cp ' + src_dir_mol + '/mol-1.pdb ./gas/mol.pdb')
        os.system('cp ff.str gas/')
        os.system('cp ff.str liquid/')
        os.system('cp ff.str sfe/')

        L_Box = query_box_size('/homes/huanglei/restored/huanglei/ff/LJ/compound/'+mol_id+'/ini_ff/ini_box/namd-npt.out.xsc')
        str_L_Box = '%f' %(L_Box)

	###liq minimization	
	inpline=createinput_namd_liquid_min(str_T_Sim,str_L_Box)
        fOut=open('./liquid/namd_NPT_min','w')
        fOut.write(inpline)
        fOut.close()	
        
	###liq md
        inpline=createinput_namd_liquid(str_T_Sim,str_L_Box,mol_id)
        fOut=open('./liquid/namd_NPT','w')
        fOut.write(inpline)
        fOut.close()
        
        os.system('cp ' + src_dir_bulk+'/drude-box.pdb ./liquid/')
        os.system('cp ' + src_dir_bulk+'/drude-box.pdb ./gas/')
        os.system('cp ' + src_dir_bulk+'/mol-drude-liq.xpsf ./liquid/')
        os.system('cp ' + src_dir_mol+'/mol-drude-gas.xpsf ./gas/')
	
        inpline=createinput_liquid_job(str_T_Sim,mol_id,str_L_Box)
        fOut=open('./liquid/job.sh','w')
        fOut.write(inpline)
        fOut.close()

        os.chdir('./liquid/')
	
        if(not os.path.isfile("dE_dvdw_liquid.txt")):
            while (True):
                if submitNow():
	    	   print "submitting liquid job for "+mol_id
            	   retcode, result = commands.getstatusoutput('qsub job.sh')
            	   if retcode != 0:
                	print >>ferr,'Error in submitting job ' + mol_id + ' for liquid simulation.'
                	exit(1)
            	   else :
                   	jobid = string.split(result,'.')
                   	jobid = jobid[0]
                   	Job_List.append(jobid)
		   	break	### break out of this while loop
		else: time.sleep(15)

    nMol = nMol + 1

f.close()

os.chdir(workdir)

"""
submit gas jobs in parallel for all molecules after submitting liquid jobs
"""
if Done==0: submitGasMD()
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
    N_Mol_l = int(process.communicate()[0])

    print os.getcwd() 

    if not os.path.isfile("gas/dE_dvdw_gas.txt"): 
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

    ### Lets compute gradients for the gas jobs
    if not os.path.isfile("gas/dE_dvdw_gas.txt"): 
	os.chdir("gas/coords-vels")	 
    	L_Box = query_box_size('/homes/huanglei/restored/huanglei/ff/LJ/compound/'+mol_id+'/ini_ff/ini_box/namd-npt.out.xsc')
    	str_L_Box = '%f' %(L_Box)
	
    	#os.system("/lcrc/project/Drude/chetan/workspace/ljopt-independant-code/code_opt_lj_grad_gas/new/grad_gas_e ../ff.str ../mol-drude-gas.xpsf "+str(T_Sim)+" 2000000 0.0005 "+str_L_Box+" 100 > log.txt")
    	if not os.path.isfile("dE_dvdw_gas.txt"):
    		os.system("../../../../../code_grad_gas/grad_gas ../ff.str ../mol-drude-gas.xpsf "+str(T_Sim)+" 1000000 0.001 "+str_L_Box+" 0 > log.txt")
    	os.system("cp dE_dvdw_gas.txt ../")
    	os.chdir("../")
    	os.system("tar -zcf coords-vels.tar.gz coords-vels")
    	os.system("tar -zcf logs.tar.gz logs")
    	os.system("rm -r coords-vels logs")	
    	os.chdir("../")	

   	
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


