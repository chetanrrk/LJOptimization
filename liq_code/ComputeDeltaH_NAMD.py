import os,sys,math,time
import commands,subprocess

base = os.getcwd()
logs = base +"/gas/logs/"
boxName = "liquid/drude-box.pdb" ###name of the liq box
T_Sim = float(sys.argv[1]) ###simulation temperature
K_BOTZMAN = 0.001987191
RT = T_Sim * K_BOTZMAN

command1 = "tail -n 3 "+boxName+" | grep ""ATOM"" | awk '{print $5}'"
process = subprocess.Popen(command1,stdout=subprocess.PIPE, shell=True)
nMols = int(process.communicate()[0])


def getEnergy(logF):

    command = "grep ERROR: "+logF+" | tail -n 1"
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    try:	
    	if process.communicate()[0].split()[0] == "ERROR:": ## some error in gas simulation
	    return -1
    except IndexError:
	pass
	
    energy=[]    
    fin = open(logF,"r")
    for line in fin:
        token=line.rstrip().split()
        if len(token)>=14:
           if token[0]=="ENERGY:":
              #if int(token[1])>=200:energy.append(float(token[13]))
              energy.append(float(token[13]))      
    fin.close()
    return energy        

def stats(arr):
    mean =  reduce(lambda x,y:x+y,arr)/len(arr)
    std = math.sqrt(reduce(lambda i,j: i+j, map(lambda x: (mean-x)*(mean-x), arr))/ len(arr))
    return mean,std 

def getVolume(logF):
    fin=open(logF,"r")
    vol=[]
    for line in fin:
        token=line.rstrip().split()
        if len(token)>=18:
            if token[0]=="ENERGY:":
                vol.append(float(token[18]))
    fin.close()
    return vol   


if __name__=="__main__":
    os.chdir(logs)
    avgEne = []
    for f in os.listdir("."):
        energy = getEnergy(f)
	if energy==-1:continue
        try:
            avgEne.append(stats(energy)[0])
        except TypeError:
	    continue

    ###write the pot energies in a separate file for each snapshot
    for f in os.listdir("."):
	if f.endswith("log"):
            m = f.split(".")[0].split("_")[2]
	    if os.path.isfile("E_namd_dcd_"+m+"_gas.txt"):os.system("rm E_namd_dcd_"+m+"_gas.txt")	
	    os.system("grep 'ENERGY:  ' "+f+" | awk '{print $14}' >> E_namd_dcd_"+m+"_gas.txt")

    ### average gas energies across all snapshots and all time step	
    os.chdir(base)
    gas_mean,gas_std = stats(avgEne)[0],stats(avgEne)[1]
    fout = open("enthalpy-vol-summary-tmp.txt","w")
    print >>fout,"(Gas)average,std:",gas_mean,gas_std
    fout2 = open("data_E_gas.txt","w")
    print >>fout2,gas_mean
    fout2.close()	

    liq_energy = getEnergy(base+"/liquid/log.txt")
    liq_mean,liq_std = stats(liq_energy)[0],stats(liq_energy)[1] ### deviation within the box 
    print >>fout,"(Liquid)average,std:",liq_mean,liq_std
    print >>fout,"(DeltaH)average:",gas_mean-(liq_mean/int(nMols))+RT
    
    vol = getVolume(base+"/liquid/log.txt")  
    vol_mean,vol_std = stats(vol)[0],stats(vol)[1]
    print >>fout,"(Volume)average,std:",vol_mean,vol_std
    fout.close()

    ### moving "-tmp" works as flag to signal completion of this calculation
    os.system("mv enthalpy-vol-summary-tmp.txt enthalpy-vol-summary.txt")


