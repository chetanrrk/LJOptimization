import os,sys


nodeList=[]
nodes=sys.argv[1]

def getNodeRange(prefix,start,end):
    for i in range(int(start),int(end)+1):
	nodeList.append(prefix+str(i))

def addNodes():
    if nodes.startswith("b["):
        prefix,tmpNodes=nodes.split("[")[0],nodes.split("[")[1].split("]")[0].split(",")
        for n in tmpNodes:	
            tmp = n.split("-")
            if len(tmp)==1:
		nodeList.append(prefix+tmp[0])
            else: 
        	getNodeRange(prefix,tmp[0],tmp[1])
    else:
         nodeList.append(nodes)

if __name__=="__main__":
    addNodes()
    for i in range(1,len(nodeList)+1):
    	fout = open("namd2.nodelist"+str(i),"w")
    	print>>fout,"group main"
    	print>>fout,"host",nodeList[i-1]
    	fout.close()

