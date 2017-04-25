from collections import Counter
import sys
#=======================================================================================
###Global variables
fastaFile = ""
#Usage
usage = "python readvidjilcloneout.py -f fastaFile\n"
#=======================================================================================
### Read parameters
def readParameters(args):
	global fastaFile
	for i in range(1,len(args)):
		if (args[i] == "-f"):
			fastaFile = args[i+1]
		elif (args[i] == "-h"):
			print usage
#=======================================================================================
### Check parameters
def checkParameters():
	if (fastaFile == ""):
		print "ERROR::Parameter -f fastaFile is required"
		sys.exit(1);
#=======================================================================================
def read(nomfichier):
	listresult=[]
	cloneNum="C"+(nomfichier.split("-")[1])
	print cloneNum
	file=open(nomfichier,"r")
	lines=file.readlines()
	file.close()
	for l in lines:
		if l[0:2]==">S":
			m=l.split("|")
			M=m[1]+"/"+m[2]+"/"+m[3]
			listresult.append(M)
	return dict(Counter(listresult)),cloneNum
#=======================================================================================
def writetofile(dico,cloneNum):
	filetowrite=open("outvidjil.txt","a")
	for key in dico.keys():
		gene=key.split("/")
		line=cloneNum+"\t"+gene[0]+"\t"+gene[1]+"\t"+gene[2]+"\t"+str(dico[key])+"\n"
		filetowrite.write(str(line))
	return "0"
#=======================================================================================
#										Main
#=======================================================================================
readParameters(sys.argv)
checkParameters()
dico,cloneNum=read(fastaFile)
writetofile(dico,cloneNum)

