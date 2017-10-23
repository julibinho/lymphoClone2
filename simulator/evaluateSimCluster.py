import sys

inputFile = ""
totalSeq = 0
verbose = False

usage = "python  evaluateSimCluster.py -i <cluster output> -totalSeq <total of sequences in the priginal file>\n"



########################################################################
### Read parameters
def readParameters(args):
    global inputFile, totalSeq, verbose;
    
    if (len(args)< 1):
        print (usage)
        sys.exit(1);
    
    for i in range(1,len(args)):
        if (args[i] == "-i"):
            inputFile = args[i+1]
        elif (args[i] == "-totalSeq"):
            totalSeq = int(args[i+1])
        elif (args[i] == "--v"):
            verbose = True
        elif (args[i] == "-h" or args[i] == "-help"):
            print (usage)
            sys.exit(1);
########################################################################
### Check parameters
def checkParameters():
    if (inputFile == ""):
        print ("ERROR::Parameter -i inputFile is required, this is the formated output of clustering tool\n")
        sys.exit(1);


    if (verbose):
        print ("Parameters ok ")
        print ("-i ", inputFile)
        if totalSeq > 0:
        	print ("-totalSeq ", totalSeq)
        print ("-verbose ", verbose)
        


# =============================================================================
#               evaluate Measures
# =============================================================================		
def evaluateMeasures(cluster, hashCluster, totalSeq):	
	count = 0; tp = 0; sumTp = 0.0; sumFp = 0.0; sumFn = 0.0; 
	countSeq = 0
		
	for i in cluster:
		count +=1
		if (count % 5000 == 0): print "Processed ", count
		
		arraySeqIds = cluster[i]
		
		#find the majority label
		hashAux = {}; tp = 0;
		
		for j in arraySeqIds:			
			IDmember = j.split('-')[0]
			countSeq +=1
			if IDmember in hashAux.keys():
				hashAux[IDmember] += 1
			else:
				hashAux[IDmember] = 1
		#take the max label
		maxValue = 0; maxLabel = ""
		for d in hashAux:
			if hashAux[d] > maxValue:
				maxValue = hashAux[d]
				maxLabel = d

		#count correct label
		for j in arraySeqIds:
			IDmember = j.split('-')[0]
			if IDmember !="" and IDmember == maxLabel:
				tp +=1
		
		#compute measures
		totalCluster = hashCluster[maxLabel]
		fn = totalCluster - tp
		fp = len(arraySeqIds) - tp
		sumTp += tp; sumFp += fp; sumFn += fn; 
		#print (tp, fp, fn)
		
	#print (sumTp, sumFp, sumFn)
	pre = sumTp/(sumTp+ sumFp);
	if totalSeq > 0:
		rec = sumTp/(totalSeq);
	else:
		rec = sumTp/(sumTp+ sumFn);
	print ("Pre = ", pre, "Rec = ", rec , "F-score = ", 2*pre*rec/(pre + rec), "NumCluster = ", count, "NumSeqs = ", countSeq)

	
# =============================================================================
#               Read Result file
# =============================================================================	
def readResultFile(filename):
	hashCluster = {}; cluster = {}; count = 0;
	file = open(filename, "r")
	for line in file.readlines(): 
		count +=1
		if (count % 5000 == 0): print "Processed ", count
		IDcluster = line.split('\t')[0]
		members = line.split('\t')[1]
		if (members == "\n" or members == ""):
			print "Warnning:: Cluster ", IDcluster, " has no members"
		else: 
			arraySeqIds = members.split()
			cluster[IDcluster] = arraySeqIds
			for m in arraySeqIds:
				headSeq = m.split('-')[0]
				if (headSeq not in hashCluster.keys()):
					hashCluster[headSeq] = 1
				else:
					hashCluster[headSeq] += 1
	file.close()
	return 	hashCluster, cluster	
# =============================================================================
#               Main
# =============================================================================	
readParameters(sys.argv)
checkParameters()
hashCluster, cluster = readResultFile(inputFile)
evaluateMeasures(cluster, hashCluster, totalSeq)