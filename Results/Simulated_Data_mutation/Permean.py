#####
import sys

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###Global variables
comparisonfile= ""
Tool_name = "" 
Fasta_file=""
Output_file= ""
MutNumber=""


#Usage
usage = "python Permean.py [options] -c comparisonfile -t Tool_name -m MutNumber -f Fasta_file -o Output_file\n "
usage +="-t <x>\t: Name of tool, 1 for MixCR; 2 for Partis; 3 for Vidjil; 4 for High-Vquest\n"
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Read parameters
def readParameters(args):
	global comparisonfile
	global Tool_name
	global Fasta_file
	global Output_file
	global MutNumber
	for i in range(1,len(args)):
		if (args[i] == "-c"):
			comparisonfile = args[i+1]
		elif (args[i] == "-t"):
			Tool_name = args[i+1]
		elif (args[i] == "-f"):
			Fasta_file = args[i+1]
		elif (args[i] == "-o"):
			Output_file = args[i+1]
		elif (args[i] == "-m"):
			MutNumber = args[i+1]
		elif (args[i] == "-h"):
			print usage


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Check parameters
def checkParameters():
	print comparisonfile
	if (comparisonfile == ""):
		print "ERROR::Parameter -f comparisonfile is required."
		print "Use -h option for more information. \n"
		sys.exit(1);
	elif (Tool_name == ""):
		print "ERROR::Parameter -t Tool_name is required\n"
		sys.exit(1);
	elif (Fasta_file == ""):
		print "ERROR:: Parameter -f Fasta_File is required\n"
		sys.exit(1);
	elif(Output_file == ""):
		print "ERROR:: Parameter -o Output_File is required\n"
	elif(MutNumber == ""):
		print "ERROR:: Parameter -m MutNumber is required\n"


#*********************************************************************************
#********************************************************************************* 
# 				   Read files
#*********************************************************************************
#********************************************************************************* 



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
def write_to_file(tool_list,file_to_write):
	f=open(file_to_write,"w")
	f.write("\n")
	for pair in tool_list:
		for i in range(0,2):
			if pair[i] != []:
				l=pair[i][0]+" "+pair[i][1]+" "+pair[i][2]+" "+pair[i][3]+'\n'
			else:
				l="___"+"___"+"___"+'\n'
			f.write(l)
			if i==1:
				f.write("\n")
	f.close
	return 0
	
#=================================================================================
#				      Vidjil
#=================================================================================
def read_vidjil(vidjil_output):
	lines=read_output_file(vidjil_output)
	vidjil_list=[]
	for l in lines:
		if l[:2]==">S":
			splited=l.split("\t")
			Response,rearrangement_type=split_vidjil_seq_name(splited[0],MutNumber)
			Result=split_vidjil_result(rearrangement_type,splited[-1])
			vidjil_list.append([Response,Result])
	return vidjil_list

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def split_vidjil_seq_name(seqName,MutNumber):
	splited_by_esp=seqName.split(" ")	
	if int(MutNumber)==0:
		SeqName,V,D,J=splited_by_esp[0].split("|")
	else:
		SeqName,V,D,J,Mut=splited_by_esp[0].split("|")
	return ["#",V,D,J],splited_by_esp[-1]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def split_vidjil_result(rearrangement_type,SeqResult):
	list_of_result=[]
	Res= SeqResult.split(" ")
	if rearrangement_type =="VDJ":
		list_of_result=[">",Res[0],Res[2],Res[4]]
	elif rearrangement_type =="VJ":
		list_of_result=[">",Res[0],"___",Res[2]]
	elif rearrangement_type =="VD":
		list_of_result=[">",Res[0],Res[2],"___"]
	else:
		print "error in vidjil result :: see split_vidjil_result()"
	return list_of_result

#=================================================================================
#				      Mixcr
#=================================================================================

# mixcrformat ( need seq fasta file)
def readFastaMul(nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	seq=""
	nom=""
	lesSeq={}
	for l in lines:
		if l[0] == '>':
			if seq != "":
				lesSeq[seq] = nom
			nom=l[1:-1]
			seq=""
		else:
			seq=seq+l[:-1]
	if seq != "":
		lesSeq[seq] = nom
	return lesSeq
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def split_mixcr_seq_name(seqName,MutNumber):
	print seqName.split("|")
	if MutNumber == 0:
		SeqName,V,D,J=seqName.split("|")
	else:
		SeqName,V,D,J,mut=seqName.split("|")
	return ["#",V,D,J]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def split_mixcr_result(SeqResult):
	list_of_result=[">"]
	for gene in SeqResult:
		if gene != "":
			list_of_result.append(gene.split("(")[0])
		else:
			list_of_result.append("___")
	return  list_of_result

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_Mixcr(fasta_file,mixcr_output):
	print mixcr_output
	Dico_fasta_file=readFastaMul(fasta_file)
	lines=read_output_file(mixcr_output)
	Mixcr_list=[]
	for l in range(1,len(lines)):
		Response=split_mixcr_seq_name(Dico_fasta_file[lines[l].split("\t")[0]],int(MutNumber)) #Dico[DNA sequence]=name of sequence
		Result=split_mixcr_result(lines[l].split("\t")[2:5])
		Mixcr_list.append([Response,Result])
	return Mixcr_list

#=================================================================================
#				      Partis
#=================================================================================
def cvs2txt(cvsfile,txtfile):
	try:
   		my_input_file = open(cvsfile, "r")
	except IOError as e:
    		print("I/O error({0}): {1}".format(e.errno, e.strerror))

	if not my_input_file.closed:
    		text_list = [];
    		for line in my_input_file.readlines():
        		line = line.split(",", 44)
        		text_list.append(" ".join(line))
    		my_input_file.close()

	try:
    		my_output_file = open(txtfile, "w")
	except IOError as e:
   		print("I/O error({0}): {1}".format(e.errno, e.strerror))
	if not my_output_file.closed:
    		for line in text_list:
        		my_output_file.write("  " + line)
   		print 'File Successfully written.'
    		my_output_file.close()
	return 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def readFastaMulforpartis(lines,MutNumber):
	Diconom={}
	for l in lines:
		if l[0] == '>':
			info=l.split("|")
			if MutNumber==0:
				Diconom[info[0][1:]]=["#",info[1],info[2],info[3][:-1]]
			else:
				Diconom[info[0][1:]]=["#",info[1],info[2],info[3]]
	return Diconom

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def readpartis(lines):
	DicoRes={}
	for l in range(1,len(lines)):
		split=lines[l].split(" ")
		DicoRes[split[2]]=[">",split[3],split[4],split[5]]
	return DicoRes

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_partis(fastafile,partis_output): 
	partis_list=[]
	outputfile="output_partis_"+str(fastafile.split("/")[-1].split(".")[0])+".txt"
	cvs2txt(partis_output,outputfile)
	lines_fasta,lines_partis=read_output_file(fastafile),read_output_file(outputfile)
	Diconom,DicoRes=readFastaMulforpartis(lines_fasta,int(MutNumber)),readpartis(lines_partis)
	for key in Diconom.keys():
		partis_list.append([Diconom[key],DicoRes[key]])
	return partis_list



#=================================================================================
#			           High-Vquest
#=================================================================================
def split_name(seqName,MutNumber):
	SeqName,V,D,J,mut="","","","",""
	if "|"  in seqName:
		if MutNumber==0:
			SeqName,V,D,J=seqName.split("|")
		else:
			SeqName,V,D,J,mut=seqName.split("|")
	return ["#",V,D,J]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_highVquest(inputfile):
	indice=[3,13,9]
	highVquest_list=[]
	lines=read_output_file(inputfile)
	for l in range(1,len(lines)):
		splited_line=lines[l].split("\t")
		Response=split_name(splited_line[1],int(MutNumber))
		Result=[">"]
		for i in indice:
			gene=splited_line[i].split(" ")
			if gene[0]=="":
				Result.append("___")
			else:
				Result.append(gene[1])
		highVquest_list.append([Response,Result])
		#print len(highVquest_list)
	return highVquest_list	


#*********************************************************************************
#********************************************************************************* 
#				    Compare
#*********************************************************************************
#********************************************************************************* 

def read_comparison(filename):
	list_plus_allele,list_minus_allele=[],[]
	lines=read_output_file(filename)
	for l in lines:
		if l == "\n":
			rep,res="","" 
		elif l[0]=="#":
			rep=l
		elif l[0]=='>':
			res=l
		if rep != "" and res != "":
			list_plus_allele.append(counter(rep,res,"plus"))
			list_minus_allele.append(counter(rep,res,"minus"))
	#print list_plus_allele
	finallist_plus_allele= [sum(i) for i in zip(*list_plus_allele)]
	finallist_minus_allele= [sum(i) for i in zip(*list_minus_allele)]
	return finallist_plus_allele,finallist_minus_allele,len(list_minus_allele)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def counter(Rep,Res,type_of_analyse):			#type_of_analyse : +/- allele
	listcomp=[]
	rep,res=Rep[:-1].split(" "),Res[:-1].split(" ")
	if type_of_analyse == "plus":
		listcomp=plus_allele(rep,res)
	else:
		listcomp=minus_allele(rep,res)
	if listcomp[0]==1 and listcomp[2]==1: 	#VJ  comparison
		listcomp.append(1)
	else:
		listcomp.append(0)
	if listcomp[0]==1 and listcomp[1]==1 and listcomp[2]==1:	#VDJ comparison 
		listcomp.append(1)
	else:
		listcomp.append(0)
	return listcomp

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plus_allele(rep,res):
	listcomp=[]
	for i in range(1,4):
		if rep[i]==res[i]:
			listcomp.append(1)
		else :
			listcomp.append(0)
	return listcomp

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def minus_allele(rep,res):
	listcomp=[]
	for i in range(1,4):
		if rep[i].split("*")[0]==res[i].split("*")[0]:
			listcomp.append(1)
		else :
			listcomp.append(0)
	return listcomp

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def whichtooltorun(comparisonfile,Tool_name,Fasta_file):
	list_of_comparison=[]
	if Tool_name ==1:# Mixcr
		list_of_comparison=read_Mixcr(Fasta_file,comparisonfile)
	elif Tool_name ==2: #Partis
		list_of_comparison=read_partis(Fasta_file,comparisonfile)
	elif Tool_name ==3: #Vidjil
		list_of_comparison=read_vidjil(comparisonfile)
	elif Tool_name ==4: #High-Vquest
		list_of_comparison=read_highVquest(comparisonfile)
	else:
		print usage
	return list_of_comparison

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def name(Tool_name):

	ToolName=""
	if Tool_name ==1:# Mixcr
		ToolName="Mixcr"
	elif Tool_name ==2: #Partis
		ToolName="Partis"
	elif Tool_name ==3: #Vidjil
		ToolName="Vidjil"
	elif Tool_name ==4: #High-Vquest
		ToolName="High-Vquest"
	else:
		print usage
	return ToolName

#*********************************************************************************
#********************************************************************************* 
#				     Main
#*********************************************************************************
#********************************************************************************* 

readParameters(sys.argv)
checkParameters()
list_of_gene=['V','D','J','VJ','VDJ']
list_of_comparison = whichtooltorun(comparisonfile,int(Tool_name),Fasta_file)
write_to_file(list_of_comparison,Output_file)
finallist_plus_allele,finallist_minus_allele,tested=read_comparison(Output_file)
ToolName=name(int(Tool_name))


print 
print "=== Summary, performance test, %4s (With allele) ==="  % (ToolName)
print 
print "            tested     passed     failed     //tested     //all		"
for l in range(len((list_of_gene))):
	print "    %-5s     %4d       %4d       %4d       %4f       %4f   " % (list_of_gene[l],tested, finallist_plus_allele[l], int(tested)-int(finallist_plus_allele[l]),(float(finallist_plus_allele[l])/int(tested)),(int(finallist_plus_allele[l])/float(10000)))
print
print "===================================================================================="

print 
print "=== Summary, performance test, %4s (Without allele) ==="  % (ToolName)
print 
print "            tested     passed     failed     //tested     //all		"
for l in range(len((list_of_gene))):
	print "    %-5s     %4d       %4d       %4d       %4f       %4f   " % (list_of_gene[l],tested, finallist_minus_allele[l], int(tested)-int(finallist_minus_allele[l]),(float(finallist_minus_allele[l])/int(tested)),(int(finallist_minus_allele[l])/float(10000)))
print 
print "===================================================================================="


