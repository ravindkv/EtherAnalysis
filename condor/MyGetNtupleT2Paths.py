
#//////////////////////////////////////////////////
#                                                 #
# Copy Ntuple root files from T2_IN_TIFR to T3.   #
# Merge them into a single root file at T3.       #
# Send the merged file back to T2_IN_TIFR.        #
# Store the full path of merged file for Analysis #
#                                                 #
#//////////////////////////////////////////////////


import os
import sys
import datetime


#USERS INPUTS
#-------------------------------
isMC = True
isData = True
#-------------------------------


#Sample paths at T2
MC_T2Paths = []
Data_T2Paths = []
def getSampName(line, sampName, sampPaths):
    if line.find(sampName)!=-1:
        s1 = line.replace("]","").split(",")
        for n in range(1, len(s1)-1):
            sampPaths.append(s1[n].replace(" ",""))

for line in open("ntupleT2Paths.txt"):
    line = line.strip()
    print line
    if len(line)!=0:
        getSampName(line, "MC", MC_T2Paths)
        getSampName(line, "DATA", Data_T2Paths)

def execme(command):
    print "\033[01;32m"+ "Excecuting: "+ "\033[00m",  command
    os.system(command)

#T2 path of one sample
#/cms/store/user/rverma/multicrab_09april17/MuMC_20170409/TTJets_MuMC_20170409/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_MuMC_20170409

#read T2 Paths of ntuples, merge them, send merged file back to T2
def mergeNtuples(T2Paths, range_):
    for m in range(range_):
        pathT2 = T2Paths[m].replace("'","")
        #print pathT2
        #local dir at T3 for the samples
        #dirT3 = pathT2.split("/")[-4]
        sampT3dir = pathT2.split("/")[-1]
        execme("rfdir  /dpm/indiacms.res.in/home/"+pathT2+ " >"+sampT3dir+".log")
	for line in open(sampT3dir+".log"):
	    line = line.strip()
	    sp1 = line.split(" ")
        execme("rfdir  /dpm/indiacms.res.in/home/"+pathT2+"/"+sp1[-1]+ " >"+sampT3dir+"_1.log")
	    
	for line in open(sampT3dir+"_1.log"):
	    line = line.strip()
	    sp2 = line.split(" ")
        execme("rfdir  /dpm/indiacms.res.in/home/"+pathT2+"/"+sp1[-1]+"/"+sp2[-1]+ " >"+sampT3dir+"_2.log")
	
	f = open(sampT3dir+"_Ntuples_Paths.log", "w")
	for line in open(sampT3dir+"_2.log"):
	    line = line.strip()
	    sp3 = line.split(" ")
	    print sp3[-1]
	    f.write(pathT2+"/"+sp1[-1]+"/"+sp2[-1]+"/"+sp3[-1]+"\n")
 	execme("rm "+sampT3dir+".log")	
 	execme("rm "+sampT3dir+"_1.log")	
 	execme("rm "+sampT3dir+"_2.log")	
#muon channel
range_MC = len(MC_T2Paths)
range_Data = len(Data_T2Paths)
print range_Data
if isMC:
    mergeNtuples(MC_T2Paths, range_MC)
if isData:
    mergeNtuples(Data_T2Paths, range_Data)

