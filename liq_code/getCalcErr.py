import os,sys
import numpy as np

minDir = sys.argv[1]
if __name__=="__main__":
    fin = open(minDir+"/log.txt","r")
    volErr=[];hVapErr=[]
    for line in fin:
	token=line.rstrip().split()
	volErr.append(np.abs(float(token[3])-float(token[6])))
	hVapErr.append(np.abs(float(token[9])-float(token[12])))
    fin.close()
    print "volErr:",np.average(volErr)
    print "hVapErr:",np.average(hVapErr)	

