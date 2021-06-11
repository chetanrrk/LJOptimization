import os,sys,math,time
import string
import numpy as np
import commands
import pickle

ferr=open("gas_sub_error.txt","w")
fsubmit=open("submitted_jobs.txt","w")

baseDir = "." #os.getcwd()
configPath = baseDir+"/configPath"
logPath = baseDir+"/logs"
minLogPath = baseDir+"/min_logs"
jobsDir = baseDir+"/jobs"
coordsVelPath = baseDir+"/coords-vels"
pdbPath = baseDir+"/pdbs"


jobIDList = [] ### stores job IDs of all gas mols

molName = ""   ### residue name
liq_box = ""   ### name of the liq box pdb
nMols = 0      ### mols in the liq box will be assigned after reading input	
xpsf = ""      ### name of the xpsf
ffStr = ""     ### name of the ff.str

temp = 0.0     ### simulation temp
mdSteps = 0    ### Total MD Time Steps
timeStep = 0   ### time step in fs	    

nDraws = 0     ### sets of mols randomly drawn 
molPicks = 0   ### num. of mols in a set
ensembleAvg = False ### to apply bltz averaging
precision = 0.0 ### precision*KT is the stopping criteria  
#partition = "" ### partition to submit the job
nodes = 0      ### number of nodes used to run
tasks_per_node = 0   ### number of processor per node
#runTime = "00:00:00" ### time allowed for single gas MD to run in hh:mm:ss 

K_BOTZMAN = 0.001987191

def readInp(gasIn):
    fin = open(gasIn,"r")

    for line in fin:
	token = line.rstrip().split("=")
	try:
	   if token[0].strip()=="molName":
	   	global molName
	   	molName = token[1].split("#")[0].strip()
	   elif token[0].strip()=="liq_box":
	   	global liq_box 
	   	liq_box = token[1].split("#")[0].strip()
	   elif token[0].strip()=="xpsf":
	   	global xpsf 
	   	xpsf = token[1].split("#")[0].strip()
	   elif token[0].strip()=="ffStr":
	   	global ffStr 
	   	ffStr = token[1].split("#")[0].strip()	
           elif token[0].strip()=="temp":
           	global temp
           	temp = float(token[1].split("#")[0].strip())
           elif token[0].strip()=="mdSteps":
           	global mdSteps
           	mdSteps = int(token[1].split("#")[0].strip())
           elif token[0].strip()=="timeStep":
           	global timeStep
           	timeStep = float(token[1].split("#")[0].strip())
	   elif token[0].strip()=="nDraws":
	   	global nDraws 
	   	nDraws = int(token[1].split("#")[0].strip())
	   elif token[0].strip()=="molPicks":
	   	global molPicks 
	   	molPicks = int(token[1].split("#")[0].strip())
	   elif token[0].strip()=="ensembleAvg":
		global ensembleAvg
		ensembleAvg = token[1].split("#")[0].strip()
		if ensembleAvg == "True":ensembleAvg = True
		else: ensembleAvg = False
	   elif token[0].strip()=="precision":
		global precision
		precision = float(token[1].split("#")[0].strip())
	   #elif token[0].strip()=="partition":
	   #	global partition
	   #	partition = token[1].split("#")[0].strip()
	   elif token[0].strip()=="nodes":
	   	global nodes
	   	nodes = token[1].split("#")[0].strip()
	   elif token[0].strip()=="tasks_per_node":
	   	global tasks_per_node
	   	tasks_per_node = token[1].split("#")[0].strip()
	   #elif token[0].strip()=="runTime":
	   #	global runTime
	   #	runTime = token[1].split("#")[0].strip()
	except IndexError:continue 

    fin.close()		
   	
    global nMols 
    nMols = int( [l.rstrip().split()[4] for l in open(liq_box,"r") if len(l.rstrip().split())>=4 and l.rstrip().split()[0]=="ATOM"][-1])		
    #if nDraws*molPicks>nMols:sys.exit("Error: Total mols drawn is larger than total mols. Exiting!")

def writeNAMDConf(molName,molNum,xpsf,ffStr,temp):
    fout = open("gas_"+molName+"_"+molNum+".inp","w")
    print >> fout,"# INPUT"
    print >> fout,"coordinates             ../coords-vels/"+molNum+"_min.out.coor"	
    print >> fout,"structure               ../"+xpsf
    print >> fout,"parameters              ../"+ffStr
    print >> fout,"paraTypeCharmm          on\n"
    print >> fout,"temperature "+str(temp)+"\n"	
    print >> fout,"# output params"
    print >> fout,"outputname      ../coords-vels/"+molNum+".out"
    print >> fout,"binaryoutput    no"
    print >> fout,"outputEnergies  100\n"
    print >> fout,"DCDfreq   100"
    print >> fout,"DCDfile   ../coords-vels/"+molNum+"_gas.dcd\n"
    print >> fout,"#drude"
    print >> fout,"drude           on"
    print >> fout,"drudeTemp       1.0"
    print >> fout,"drudeBondLen    0.2"
    print >> fout,"drudeHardWall   yes\n"       #"drudeBondConst  500\n"    
    print >> fout,"langevin        on"
    print >> fout,"langevinTemp    "+str(temp)
    print >> fout,"langevinDamping 5.0"
    print >> fout,"langevinHydrogen on\n"
    print >> fout,"exclude                 scaled1-4"
    print >> fout,"1-4scaling              1.0"
    print >> fout,"switching               on"
    print >> fout,"switchdist              1000.0"
    print >> fout,"cutoff                  1200.0"
    print >> fout,"pairlistdist            1500.0\n"
    print >> fout,"#SHAKE"
    print >> fout,"rigidbonds              all"
    print >> fout,"timestep                "+str(timeStep)+" #fs"
    print >> fout,"nonbondedFreq           1"
    print >> fout,"fullElectFrequency      1"
    print >> fout,"stepspercycle           10\n"
    print >> fout,"useGroupPressure        yes ;# needed for rigidBonds"
    print >> fout,"useFlexibleCell         no"
    print >> fout,"useConstantArea         no"
    print >> fout,"exec hostname"
    print >> fout,"run "+str(mdSteps)+"\n"
    fout.close()	


def writeNAMDMinConf(molName,molNum,xpsf,ffStr,temp):
    fout = open("gas_"+molName+"_"+molNum+"_min.inp","w")
    print >> fout,"# INPUT"
    print >> fout,"coordinates             ../pdbs/mol_"+molNum+".pdb"
    print >> fout,"structure               ../"+xpsf
    print >> fout,"parameters              ../"+ffStr
    print >> fout,"paraTypeCharmm          on\n"
    print >> fout,"temperature "+str(temp)+"\n"
    print >> fout,"# output params"
    print >> fout,"outputname      ../coords-vels/"+molNum+"_min.out"
    print >> fout,"binaryoutput    no"
    print >> fout,"outputEnergies  100\n"
    print >> fout,"#drude"
    print >> fout,"drude           on"
    print >> fout,"drudeTemp       1.0"
    print >> fout,"drudeBondLen    0.2"
    print >> fout,"drudeHardWall   yes\n"       #"drudeBondConst  500\n"   
    print >> fout,"langevin        on"
    print >> fout,"langevinTemp    "+str(temp)
    print >> fout,"langevinDamping 5.0"
    print >> fout,"langevinHydrogen on\n"
    print >> fout,"exclude                 scaled1-4"
    print >> fout,"1-4scaling              1.0"
    print >> fout,"switching               on"
    print >> fout,"switchdist              1000.0"
    print >> fout,"cutoff                  1200.0"
    print >> fout,"pairlistdist            1500.0\n"
    print >> fout,"#SHAKE"
    print >> fout,"rigidbonds              all"
    print >> fout,"timestep                "+str(timeStep)+" #fs"
    print >> fout,"nonbondedFreq           1"
    print >> fout,"fullElectFrequency      1"
    print >> fout,"stepspercycle           10\n"
    print >> fout,"useGroupPressure        yes ;# needed for rigidBonds"
    print >> fout,"useFlexibleCell         no"
    print >> fout,"useConstantArea         no"
    print >> fout,"exec hostname"
    print >> fout,"minimize 1000 \n"
    fout.close()


def getPDB(molNum):
    fin = open(liq_box,"r")
    fout = open("mol_"+molNum+".pdb","w")
    for line in fin:
        token = line.rstrip().split()
        if len(token)>5:
            if token[4]==molNum:print>>fout,line.rstrip()
    print>>fout,"END"
    fin.close()
    fout.close()


def writeQMinJob(jobNum,molName,molNum,copyIdx):
    procs = int(tasks_per_node)
    if not os.path.isfile("min_job_"+jobNum+".sh"):
    	fout = open("min_job_"+jobNum+".sh","w")
    	print>>fout,"export cDir=/lcrc/project/Drude/chetan/software/NAMD_2.13_Linux-x86_64-TCP/\n"
    else: fout = open("min_job_"+jobNum+".sh","a")	

    print>>fout,"$cDir/charmrun $cDir/namd2 +p "+str(procs)+" ++nodelist namd2.nodelist"+str(copyIdx)+" ++mpiexec ++remote-shell mpiexec ../"+configPath+"/gas_"+molName+"_"+molNum+"_min.inp > ../"+minLogPath+"/"+molName+"_"+molNum+"_min.log &\n"

    fout.close()
	

def writeQMDJob(jobNum,molName,molNum,copyIdx):
    procs = int(tasks_per_node)
    if not os.path.isfile("md_job_"+jobNum+".sh"):
        fout = open("md_job_"+jobNum+".sh","w")
    	print>>fout,"export cDir=/lcrc/project/Drude/chetan/software/NAMD_2.13_Linux-x86_64-TCP/\n"
    else: fout = open("md_job_"+jobNum+".sh","a")

    print>>fout,"$cDir/charmrun $cDir/namd2 +p "+str(procs)+" ++nodelist namd2.nodelist"+str(copyIdx)+" ++mpiexec ++remote-shell mpiexec ../"+configPath+"/gas_"+molName+"_"+molNum+".inp > ../"+logPath+"/"+molName+"_"+molNum+".log &\n"

    fout.close()


def runMin(jobNum):
    if os.path.isdir(jobsDir):
    	os.system("cd "+jobsDir+"; bash min_job_"+jobNum+".sh; cd ../")    
    else: 
	print "problem with changing dir. Exiting calculations!",jobsDir
	sys.exit()
	

def runMD(jobNum):
    if os.path.isdir(jobsDir):
    	os.system("cd "+jobsDir+"; bash md_job_"+jobNum+".sh; cd ../")    
    else:
	print "problem with changing dir. Exiting calculations!",jobsDir
	sys.exit()

def writePotEne():
    os.chdir("./logs")
    for f in os.listdir("."):
        if f.endswith("log"):
            m = f.split(".")[0].split("_")[2]
            if os.path.isfile("E_namd_dcd_"+m+"_gas.txt"):os.system("rm E_namd_dcd_"+m+"_gas.txt")
            os.system("grep 'ENERGY:  ' "+f+" | awk '{print $14}' >> E_namd_dcd_"+m+"_gas.txt")	
    os.chdir("../") 

def query_box_size(szName):
    f=open(szName,"r")
    line = f.readline()
    line = f.readline()
    line = f.readline()
    szBuff = line.split()
    L_Box = string.atof(szBuff[1])
    return L_Box

def computeGrads():
    if not os.path.isfile("dE_dvdw_gas.txt"):
	os.chdir("./coords-vels")
	L_Box = query_box_size('/homes/huanglei/restored/huanglei/ff/LJ/compound/'+molName+'/ini_ff/ini_box/namd-npt.out.xsc')
	str_L_Box = '%f' %(L_Box)
	if not os.path.isfile("dE_dvdw_gas.txt"):
	    os.system("/lcrc/project/Drude/chetan/workspace/ljopt-tests/drude_lj_prod/code_grad_gas/grad_gas ../ff.str ../mol-drude-gas.xpsf "+str(temp)+" 1000000 0.001 "+str_L_Box+" 0 > log.txt")	
	os.system("cp dE_dvdw_gas.txt ../")
	os.chdir("../")

def run(nDraws,molPicks):
    idx = [i for i in range(1,nMols+1)]
    np.random.shuffle(idx)	

    stdErrs= []
    countErrs = 0 ## if stdErr for 3 consecutive draws are same we can safely stop drawing
    jobNum=1 ###number of gas jobs written	

    ###setting up directories to write outputs	
    if not os.path.isdir(configPath):os.mkdir(configPath)
    if not os.path.isdir(minLogPath):os.mkdir(minLogPath)
    if not os.path.isdir(logPath):os.mkdir(logPath)
    if not os.path.isdir(jobsDir):os.mkdir(jobsDir)
    if not os.path.isdir(coordsVelPath):os.mkdir(coordsVelPath)
    if not os.path.isdir(pdbPath):os.mkdir(pdbPath)

    os.system("mv namd2.nodelist* jobs")	

    fout = open("stdErr.txt","w")
    print >>fout, '{:10}'.format('#Draws'),'{:20}'.format('Std.Err'),'{:20}'.format('Threshold')	
    start = 0; end = molPicks
    #print "molName,liq_box",molName,liq_box
    #print "draws,picks",nDraws,molPicks
	
    for i in range(nDraws):
	if start > len(idx):
	   print "all mol copies tried...can't reduce error more"
	   break 	
	print "###########Draw "+str(i)+"###########"
	copyIdx=1 ### copy index in a run
	for j in idx[start:end]:
	    writeNAMDMinConf(molName,str(j),xpsf,ffStr,temp)				
	    writeNAMDConf(molName,str(j),xpsf,ffStr,temp)				
	    getPDB(str(j))	    		
    	    os.system("mv mol*.pdb pdbs")
    	    os.system("mv gas*inp "+configPath)
    	    writeQMinJob(str(jobNum),molName,str(j),copyIdx)
    	    writeQMDJob(str(jobNum),molName,str(j),copyIdx)
	    copyIdx+=1
	outF = open("min_job_"+str(jobNum)+".sh","a")
	print>>outF,"wait"
	outF.close()
	outF = open("md_job_"+str(jobNum)+".sh","a")
	print>>outF,"wait"
	outF.close()
    	os.system("mv *job*sh "+jobsDir)
	runMin(str(jobNum)) ## runs minimization on all picked copies
	runMD(str(jobNum))  ## runs md on all minimized copies		    
	
	start = end
	end = end+molPicks	
	jobNum+=1

	"""
	Checks if submitted gas jobs are done and waits for completion
        Check varience in Potential Energy and decide if next draw is required
	Break if Avg_EPot_Gas <= precision*KT else keep adding mols
	"""
	print "checking jobs"
	#checkGasJobs()			### checks and wait for gas jobs if needed
	allEne = collectAllPotEnergy()  ### all energies seen over all mols all MD steps
	stdErr = np.std(allEne)/np.sqrt(len(allEne))

	stdErrs.append(stdErr)


	threshold = precision*K_BOTZMAN*temp
	#print >>fout, str(i+1),stdErr, threshold
	print >>fout,'{:10}'.format(str(i+1)),'{:20}'.format(str(stdErr)),'{:20}'.format(str(threshold))
	if stdErr <= threshold:
	   print "termination criteria met...compute enthalpy now"		
	   break

	elif len(stdErrs)>=3:
	   if stdErr == stdErrs[-1] and stdErr == stdErrs[-2] : ### checking past two draws std errors
	   	print "saturated std error ... compute enthalpy now"
	   	break	  

	else:
	   print "std. error: "+str(stdErr)+" threshold: "+str(threshold) 	   
    
    fout.close() ##stderr file
    fout=open("stats_E_gas.txt","w")
    print>>fout,"mean,std: ",np.average(allEne),np.std(allEne)   
    fout.close()
    fout=open("data_E_gas.txt","w")
    print>>fout,np.average(allEne)
    fout.close()
    os.system("cp stats_E_gas.txt ../")
    writePotEne()	
    computeGrads()

def submitNow(r=32,q=100): #r= no.running jobs; q= no.queued jobs
        i=0
        rJobs=0 ###running jobs
        pdJobs=0 ###pending jobs
        retcode, result = commands.getstatusoutput('squeue -u rupakhetic')
        if not retcode:
                result= result.split()
                lim= len(result)/8
                for j in range(lim):
                        status = result[i:(i+8)][4]
                        if status=="R":rJobs+=1
                        elif status=="PD":pdJobs+=1
                        i = i+8

        #if rJobs+pdJobs<q and rJobs<32: return True
        if rJobs<32: return True
        else: return False

"""
NOTE: code needs to wait if it cannot submit job right away!! Else leads to all kinds of problems !!
"""
def submitJobs(molName,molNum):
    #os.chdir(jobsDir)
    os.system("cd "+jobsDir)
    while(True):
        if submitNow():
                retcode, result = commands.getstatusoutput("sbatch job_"+molName+"_"+molNum+".sh")
                if retcode != 0:
                   print >> ferr,"Error in submitting job_" + molName +"_"+ molNum +" for gas simulation."
                   exit(1)
                else:
                    try:
			jobid = result.split()[3]
			print>>fsubmit,jobid,molName +"_"+ molNum
                        jobIDList.append(jobid)
                        break   ### job is submitted and we can break out now
                    except IndexError:pass
        else: time.sleep(10)

    #os.chdir(baseDir)
    os.system("cd "+baseDir)


"""
Checks if all submitted gas calculation are done. Waits for all gas jobs to be done.
"""
def checkGasJobs():
    while(1):
        Done = 1
        for jobid in jobIDList:
            retcode, result = commands.getstatusoutput('squeue --job ' + jobid)	
            if retcode == 0: 
		if result.split()[12]=="PD" or result.split()[12]=="R":
               		Done = 0
                	break
	    else: print "problem checking gas job", jobid	
	    
        if(Done):
            break
        else:
            time.sleep(15)


def collectAllPotEnergy():
    ### Jobs are done if we are here ... get all pot. energies from log files
    os.chdir(logPath)		
    allEne = []
    for f in os.listdir("."):
	a_mol_energy = getNAMDEnergy(f)
	allEne = np.concatenate((allEne,a_mol_energy))		
    os.chdir("../")
    return allEne


def getNAMDEnergy(logF):
    energy=[]
    try:
    	fin = open(logF,"r")
    	for line in fin:
            token=line.rstrip().split()
            if len(token)>=14:
           	if token[0]=="ENERGY:":
              	   energy.append(float(token[13]))
    except IOError:print "problem with "+logF	
    fin.close()
    return energy


if __name__=="__main__":
   gasIn = sys.argv[1]
   readInp(gasIn) ##reads in inputs for gas phase calculation	 
   run(nDraws,molPicks) 
   fout=open("gas_done.txt","w")
   fout.close()
   ferr.close()
   fsubmit.close()

