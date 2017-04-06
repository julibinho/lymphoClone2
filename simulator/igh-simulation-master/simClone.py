from Bio import SeqIO
from random import *
import sys
import os


verbose = False
inputFile = ""
outputFile = ""
nbM = 5 #max number of mutations
n = 100 #number of sequences, default 1000
maxClones = 20 # max number of clones
cloneDist = [0.65, 0.23, 0.02]
randNbSeed = 1
cloneFile = ""
SEP_CLONE_FILE = "\t"
PERC_CLONE = 5


usage = "python [options] simClone.py -i <fasta sequences> -o <output file>\nwhere:\n"
usage +="\nBasic options are:\n-n <x> : number of sequences to be generated, default 1000\n"
usage +="-cloneDist <x>, clone distribution, default 0.65,0.23, 0.02\n"
usage +="-cloneFile <x>, file having clone distribution\n"
usage +="-nbM <x>, max number of mutations per sequence, default 5\n"
usage +="-maxClones <x>, max number of clones, default 20\n"
usage +="--v print log messages\n"



########################################################################
### Read parameters
def readParameters(args):
    global inputFile, outputFile, n, cloneDist, nbM, maxClones, cloneFile, verbose;
    
    if (len(args)==1):
        print (usage)
        sys.exit(1);
    
    for i in range(1,len(args)):
        if (args[i] == "-i"):
            inputFile = args[i+1]
        elif (args[i] == "-o"):
            outputFile = args[i+1]
        elif (args[i] == "-cloneFile"):
            cloneFile = args[i+1]
        elif (args[i] == "-cloneDist"):
            cloneDist = args[i+1].split(",")  
        elif (args[i] == "-n"):
            n = int(args[i+1])
        elif (args[i] == "-nbM"):
            nbM = int(args[i+1])
        elif (args[i] == "--v"):
            verbose = True
        elif (args[i] == "-maxClones"):
            maxClones = int(args[i+1])
        elif (args[i] == "-h" or args[i] == "-help"):
            print (usage)
            sys.exit(1);
########################################################################
### Check parameters
def checkParameters():
    if (inputFile == ""):
        print ("ERROR::Parameter -i inputFile is required, this is file produced by ighsim.py script\n")
        sys.exit(1);
    elif (outputFile == ""):
        print ("ERROR::Parameter -o outputFile is required, this file will contain the clone population\n")
        sys.exit(1);

    if (verbose):
        print ("Parameters ok ")
        print ("-i ", inputFile)
        print ("-o ", outputFile)
        print ("-Number of sequences ", n)
        print ("-cloneDistribution ", cloneDist)
        print ("-cloneFile ", cloneFile)
        print ("-nbM ", nbM)
        print ("-maxClones ", maxClones)
        print ("-verbose ", verbose)
        
########################################################################
### Read sequences to list
def readFastaList(inputFile):
    seqList = []
    for record in SeqIO.parse(inputFile, "fasta"):
        seqList.append(record)
    return seqList

########################################################################
### Save a file
def saveFile(fileName, info):
    f = open(fileName, 'w')
    f.write(info)
    f.close()
########################################################################
###  generateClone
def generateClone(seqList, n, cloneDist):
    sequences = ""
    total = sum(cloneDist)
    hashIds = {}
    #generating main clones, those specified into cloneDist list
    for i in cloneDist:
        totalSeq = int(i*n)
        if (totalSeq > 0):
            stop = False
            while not stop:
                randN = randint(0,len(seqList)-1) #take one sequence randomly 
                if randN not in hashIds.keys():
                    print ("generating clone " + str(i*100) + "% with " + str(totalSeq) + " seqs")
                    generateOneClone(seqList[randN], totalSeq, randN + 1)
                    hashIds[randN] = True
                    stop = True
                if  len (hashIds.keys()) == len(seqList):
                    stop = True
                    print ("Warning :: All sequences were used in clone generation")
        else:
            print ("Warning :: It is impossible to create the clone " + str(i*100) + "% with 0 sequences")
    #generating remaing clones
    remain = 1 - total
  
########################################################################
###  saveMajorCloneSeq
def saveMajorCloneSeq(majorSeq, nbseq, outputFile):
    info = "";  count = 1; 
    #save the sequence in a temporary file
    id = majorSeq.id.split('|')
    id = id[0]
    print (id)
    for i in range(nbseq):
        ann = majorSeq.id
        info += ">" + ann.replace(id, id + "-" + str(count)) + "\n" + majorSeq.seq + "\n" 
        count +=1
    saveFile(outputFile, str(info))
########################################################################
###  generateOneClone
def generateOneClone(majorSeq, nbseq, id):
    
    cmd = "";     
    outputFileTmp = outputFile + ".tmp"
    saveMajorCloneSeq(majorSeq,  nbseq, outputFileTmp)
    cmd = "python ighmut.py "+ str(nbM) + " " + outputFileTmp +  " " + outputFile + "." + str(id) +  " " + str(randNbSeed)
    os.system(cmd)
    os.system("rm " + outputFileTmp)

########################################################################
###  readCloneFile
def readCloneFile(cloneFile):   
    clnList = []
    file = open(cloneFile, "r")
    for line in file:
        if line != "\n" and line[0:12] != "Patient_name":
            arrayLine = line.split(SEP_CLONE_FILE)
            perc = arrayLine[PERC_CLONE].replace("%", "")
            perc = float(perc)/100
            clnList.append(perc)
    #print (clnList)
    return clnList
    file.close() 
########################################################################
## Main
########################################################################
readParameters(sys.argv)
checkParameters()

seqList = readFastaList(inputFile)
if (cloneFile != ""):
    cloneDist = readCloneFile(cloneFile)
generateClone(seqList, n, cloneDist)



