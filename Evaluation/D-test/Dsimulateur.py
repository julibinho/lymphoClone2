# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:46:46 2017

@author: davi
"""

def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines
 
def readseq(lines):
    lesSeq={}
    seq=""
    nom=""
    for l in lines:
        if l[0] == '>':
            if seq != "":
                lesSeq[nom] = seq
            nom=l[1:-1].split("|")[0]         
            seq=""      
        else: 
            seq=seq+l[:-1]
    return lesSeq

def locD(lines):
    locD={}
    for l in range(1,len(lines)):
        split=lines[l].split("\t")
        print split
        start=int( split[4])+int( split[5])
        stop=start+int(split[6])
        locD[split[0]]=[ split[1], split[2], split[3],start,stop]
    return locD

def readD(lesSeq,locD):
    f=open("file_to_write_seq_plus10.txt","w")
    for key in lesSeq.keys():
        title=">"+str(locD[key][0])+"|"+str(locD[key][1])+"|"+str(locD[key][2])+"\n"
        seq=lesSeq[key][locD[key][3]:locD[key][4]]
        if seq !="" and len(seq)>10:
            f.write(title)
            f.write(seq)
            f.write("\n")
#sequenes

sequenes=read_output_file("simple_plus_indels.fas")
lesSeq=readseq(sequenes)

#localistion


localistion=read_output_file("simple_plus_indels.txt")
locD=locD(localistion)
readD(lesSeq,locD)