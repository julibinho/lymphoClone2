# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 11:45:18 2017

@author: davi
"""
import sys
import operator
from collections import Counter
import collections
import glob
import numpy as np
import random


#####################################################################

def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines

#####################################################################

def UniqSeqType(lines):
    SeqType_list=[]
    for l in lines:
        split=l.split('\t')
        SeqType_list.append([split[1],split[2],split[3]])
    uniq_SeqType_list = [list(t) for t in set(map(tuple, SeqType_list))]
    return uniq_SeqType_list

#####################################################################

def Clusters(lines):
    cluster={}
    for l in lines:
        split=l.split('\t')
        if split[0] in cluster.keys():
            cluster[split[0]].append([[split[1],split[2],split[3]],split[4][:-1]])
        else:
            cluster[split[0]]=[[[split[1],split[2],split[3]],split[4][:-1]]]
    return cluster
    
#####################################################################

# input : uniq_SeqType_list(list of lists)  : list of  rearrengements of sequences 
# output : dictance matrix  

def DistanceMatrix(uniq_SeqType_list):
    matriceD=np.zeros((len(uniq_SeqType_list),len(uniq_SeqType_list)))
    for i0 in range (0,len(uniq_SeqType_list),1 ):
        for i1 in range (0,len(uniq_SeqType_list),1 ):
            if i0!=i1:
                matriceD[i0][i1]=calculeDist(uniq_SeqType_list[i0],uniq_SeqType_list[i1])
    return matriceD
                
        
##################################################################### 

# input : 2 rearrangements of sequences in list format [V1,D1,J1] and [V2,D2,J2]
# output : distance between 2 sequences 
def calculeDist(list1,list2):
    inter=set(list1).intersection(list2)   
    return 1-(len(list(inter))/float(3))
    
#####################################################################
    
def Ai(i,cluster,key,matriceD,uniq_SeqType_list):
    compt=1
    dist=0
    for l in range(len(cluster[key])):
        if i!=cluster[key][l][0]:
            dist+=matriceD[uniq_SeqType_list.index(i)][uniq_SeqType_list.index(cluster[key][l][0])]
            compt+=1
    return dist/float(compt)

#####################################################################

def Bi(i,cluster,key,matriceD,uniq_SeqType_list):
    listb=[]
    for clu in cluster.keys():
        if key!=clu:
                compt=1
                dist=0
                for l in range(len(cluster[clu])):
                     dist+=matriceD[uniq_SeqType_list.index(i)][uniq_SeqType_list.index(cluster[clu][l][0])]
                     compt+=1
                listb.append(dist/float(compt))
    return min(listb)
            
#####################################################################

# SILHOUETTE

def silhouette(clusters,matriceD,uniq_SeqType_list):
    DicoSil={}
    num=1
    lengths=len(clusters)
    for key in clusters.keys():
        print num, "th from",lengths,"clusters"
        Sum=0
        for l in range(len(clusters[key])):
            i=clusters[key][l][0]
            ai=Ai(i,clusters,key,matriceD,uniq_SeqType_list) 
            bi=Bi(i,clusters,key,matriceD,uniq_SeqType_list)
            Sum+=calculeSil(ai,bi)
        DicoSil[key]=Sum/float(len(clusters[key]))
        num+=1
    return DicoSil

#####################################################################

def calculeSil(ai,bi):
    si=0
    if ai!=bi:
        if ai<bi:
            si=1-(ai/float(bi))
        else:
            si=(float(bi)/ai)-1
    return si

#####################################################################

# Calcule avrage silhouette width

def AvrSilWid(dicoS):
    Sum=0
    for key in dicoS.keys():
        Sum+=dicoS[key]
    return Sum/float(len(dicoS))
#####################################################################

lines=read_output_file("outPartis_0030_S22.txt")
uniq_SeqType_list=UniqSeqType(lines)
cluester=Clusters(lines)
matriceD=DistanceMatrix(uniq_SeqType_list)
dicoS=silhouette(cluester,matriceD,uniq_SeqType_list)
print AvrSilWid(dicoS)


