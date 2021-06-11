import os,sys,pickle

def readOptLJ(optFFPath):
    ljs = {}
    opt_ljs = open(optFFPath+"/vdw-param.txt","r")
    for line in opt_ljs:
	token=line.rstrip().split()
	ljs[token[0]]=[token[1],token[2]]
    opt_ljs.close()	
    return ljs


def isAtomTypeMatch(typeList,atomTypeToMatch,optAtomType):
    matchTypes=[]
    for ts in typeList:
        if optAtomType in ts:
            matchTypes = ts
    return atomTypeToMatch in matchTypes


def updateVdw(ffStr,vdwParams,m,outF):
    fin=open("../gaamp_atomtypes_to_optimize.txt","r") ### clusters gaamp-drude types for all molecules in a class used 
    atomTypes=[]
    for line in fin:
        token=line.rstrip().split()
        atomTypes.append(token)
    fin.close()

    fin = open(ffStr,"r")
    fout = open(outF,"w")
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


mols = pickle.load(open("../alihatic_alchohols.p","r"))

for m in mols:
    optFFPath = "/lcrc/project/Drude/chetan/workspace/ljopt-tests/drude_lj_prod/alkanes/fullMultiTypes/run_17"
    ljs=readOptLJ(optFFPath)
    print "processing",m
    src_dir_bulk = '/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/'+m+'/drude-result/box'
    if not os.path.isfile(src_dir_bulk+'/ff_avg.str'): continue	
    os.system('cp '+src_dir_bulk+'/ff_avg.str ff_'+m+'.str')	
    updateVdw('ff_'+m+'.str',ljs,m,'temp_'+m+'_1.str')
    #break
    
    optFFPath = "/lcrc/project/Drude/chetan/workspace/ljopt-tests/drude_lj_prod/alkenes/run_23"
    ljs=readOptLJ(optFFPath)
    updateVdw('temp_'+m+'_1.str',ljs,m,'temp_'+m+'_2.str')    
    #break
 
    optFFPath = "/lcrc/project/Drude/chetan/workspace/ljopt-tests/drude_lj_prod/alkynes/run_10"
    ljs=readOptLJ(optFFPath)
    updateVdw('temp_'+m+'_2.str',ljs,m,'temp_'+m+'_3.str')    

    os.system("mv temp_"+m+"_3.str ff_new_"+m+".str")
    os.system("rm temp_*")


