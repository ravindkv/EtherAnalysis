
#//////////////////////////////////////////////////
#                                                 #
# Merge many into a single root fil.              #
#                                                 #
#//////////////////////////////////////////////////

import os
import sys
import datetime

def execme(command):
    print ""
    print "\033[01;32m"+ "Excecuting: "+ "\033[00m",  command
    print ""
    os.system(command)

#USERS INPUTS
is2016Data = True
is2017Data = False

if(is2016Data):
    #USERS INPUTS
    isMC = True
    isData = False
    mc = []
    mc.append("DYJetsToLL_M50")
    mc.append("DYJetsToLL_Pt100To250")
    mc.append("DYJetsToLL_Pt250To400")
    mc.append("DYJetsToLL_Pt400To650")
    mc.append("DYJetsToLL_Pt50To100") 
    mc.append("DYJetsToLL_Pt650ToInf")
    mc.append("TT")
    mc.append("WJetsToLNu")
    mc.append("WW")
    mc.append("WZ")
    mc.append("ZZ")
    mc.append("ExLepMuMuZ_M1000")
    mc.append("ExLepMuMuZ_M1250")
    mc.append("ExLepMuMuZ_M1500")
    mc.append("ExLepMuMuZ_M1750")
    mc.append("ExLepMuMuZ_M2000")
    mc.append("ExLepMuMuZ_M250") 
    mc.append("ExLepMuMuZ_M2500")
    mc.append("ExLepMuMuZ_M3000")
    mc.append("ExLepMuMuZ_M3500")
    mc.append("ExLepMuMuZ_M4000")
    mc.append("ExLepMuMuZ_M4500")
    mc.append("ExLepMuMuZ_M500") 
    mc.append("ExLepMuMuZ_M5000")
    mc.append("ExLepMuMuZ_M750") 
    mc.append("ExLepEEZ_M1000")
    mc.append("ExLepEEZ_M1250")
    mc.append("ExLepEEZ_M1500")
    mc.append("ExLepEEZ_M1750")
    mc.append("ExLepEEZ_M2000")
    mc.append("ExLepEEZ_M250")
    mc.append("ExLepEEZ_M2500")
    mc.append("ExLepEEZ_M3000")
    mc.append("ExLepEEZ_M3500")
    mc.append("ExLepEEZ_M4000")
    mc.append("ExLepEEZ_M4500")
    mc.append("ExLepEEZ_M500") 
    mc.append("ExLepEEZ_M5000")
    mc.append("ExLepEEZ_M750") 

    data = ["MuRunB2v2", 
            "MuRunCv1", 
            "MuRunDv1", 
            "MuRunEv1", 
            "MuRunFv1", 
            "MuRunGv1", 
            "MuRunH2v1", 
            "MuRunH3v1",
            "EleRunBver2v2", 
            "EleRunCv1", 
            "EleRunDv1", 
            "EleRunEv1", 
            "EleRunFv1", 
            "EleRunGv1", 
            "EleRunHver2v1",
            "EleRunHver3v1"]
    
    if isMC:
        for samp in range(len(mc)):
	        execme("hadd -k "+str(mc[samp])+"_Merged.root "+str(mc[samp])+"_MC*")
    
    if(isData):
        for samp in range(len(data)):
            execme("hadd -k "+str(data[samp])+"_Merged.root "+str(data[samp])+"*")
        execme("hadd -k all_muData.root MuRun*_Merged.root")
        execme("hadd -k all_eleData.root EleRun*_Merged.root")
    
    if(isMC):
    	execme("hadd -k all_DY_Pt.root DY*_Pt*Merged.root")
    	execme("hadd -k all_DY_M50.root DYJetsToLL_M50_Merged.root")
    	execme("hadd -k all_DY.root all_DY_Pt.root all_DY_M50.root")
    	execme("hadd -k all_TT.root TT_Merged.root")
    	execme("hadd -k all_WJets.root  WJ*_Merged.root")
    	execme("hadd -k all_VV.root WW_Merged.root WZ_Merged.root ZZ_Merged.root")
    	execme("hadd -k all_Bkg.root all_DY_M50.root all_TT.root all_WJets.root all_VV.root")
        execme("hadd -k all_ExLepMuMuZ_M250.root ExLepMuMuZ_M250_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M500.root ExLepMuMuZ_M500_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M750.root ExLepMuMuZ_M750_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M1000.root ExLepMuMuZ_M1000_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M1250.root ExLepMuMuZ_M1250_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M1500.root ExLepMuMuZ_M1500_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M1750.root ExLepMuMuZ_M1750_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M2000.root ExLepMuMuZ_M2000_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M2500.root ExLepMuMuZ_M2500_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M3000.root ExLepMuMuZ_M3000_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M3500.root ExLepMuMuZ_M3500_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M4000.root ExLepMuMuZ_M4000_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M4500.root ExLepMuMuZ_M4500_Merged.root")
        execme("hadd -k all_ExLepMuMuZ_M5000.root ExLepMuMuZ_M5000_Merged.root")
        execme("hadd -k all_ExLepEEZ_M250.root ExLepEEZ_M250_Merged.root")
        execme("hadd -k all_ExLepEEZ_M500.root ExLepEEZ_M500_Merged.root")
        execme("hadd -k all_ExLepEEZ_M750.root ExLepEEZ_M750_Merged.root")
        execme("hadd -k all_ExLepEEZ_M1000.root ExLepEEZ_M1000_Merged.root")
        execme("hadd -k all_ExLepEEZ_M1250.root ExLepEEZ_M1250_Merged.root")
        execme("hadd -k all_ExLepEEZ_M1500.root ExLepEEZ_M1500_Merged.root")
        execme("hadd -k all_ExLepEEZ_M1750.root ExLepEEZ_M1750_Merged.root")
        execme("hadd -k all_ExLepEEZ_M2000.root ExLepEEZ_M2000_Merged.root")
        execme("hadd -k all_ExLepEEZ_M2500.root ExLepEEZ_M2500_Merged.root")
        execme("hadd -k all_ExLepEEZ_M3000.root ExLepEEZ_M3000_Merged.root")
        execme("hadd -k all_ExLepEEZ_M3500.root ExLepEEZ_M3500_Merged.root")
        execme("hadd -k all_ExLepEEZ_M4000.root ExLepEEZ_M4000_Merged.root")
        execme("hadd -k all_ExLepEEZ_M4500.root ExLepEEZ_M4500_Merged.root")
        execme("hadd -k all_ExLepEEZ_M5000.root ExLepEEZ_M5000_Merged.root")

if(is2017Data):
    #USERS INPUTS
    isMuData = True
    isMuMC = True
    isEleData = False
    isEleMC = False
    
    mc = []
    mc.append("TT")
    mc.append("DYJetsToLL_M50")
    mc.append("WJetsToLNu")
    mc.append("WW")
    mc.append("WZ")
    mc.append("ZZ")
    
    if(isMuMC):
        mc.append("ExLepMuMuZ_M250")  
        mc.append("ExLepMuMuZ_M500")   
        mc.append("ExLepMuMuZ_M750")      
        mc.append("ExLepMuMuZ_M1000")   
        mc.append("ExLepMuMuZ_M1250")     
        mc.append("ExLepMuMuZ_M1500")   
        mc.append("ExLepMuMuZ_M1750")     
        mc.append("ExLepMuMuZ_M2000")   
        mc.append("ExLepMuMuZ_M2500")     
        mc.append("ExLepMuMuZ_M3000")   
        mc.append("ExLepMuMuZ_M3500")  
        mc.append("ExLepMuMuZ_M4000")
	mc.append("ExLepMuMuZ_M4500")
	mc.append("ExLepMuMuZ_M5000")
    else:
        mc.append("ExLepEEZ_M250")
        mc.append("ExLepEEZ_M500")
        mc.append("ExLepEEZ_M750")
        mc.append("ExLepEEZ_M1000")
        mc.append("ExLepEEZ_M1250")
        mc.append("ExLepEEZ_M1500")
        mc.append("ExLepEEZ_M1750")
        mc.append("ExLepEEZ_M2000")
	mc.append("ExLepEEZ_M2500")
	mc.append("ExLepEEZ_M3000")
	mc.append("ExLepEEZ_M3500")
	mc.append("ExLepEEZ_M4000")
	mc.append("ExLepEEZ_M4500")
	mc.append("ExLepEEZ_M5000")
        
    if(isMuData):
    	data = ["MuRunBv1", "MuRunCv1", "MuRunDv1", "MuRunEv1", "MuRunFv1", "MuRunGv1", "MuRunHv1"]
    
    if(isEleData): 
    	data = ["EleRunBver2v2", "EleRunCv1", "EleRunDv1", "EleRunEv1", "EleRunFv1", "EleRunGv1", "EleRunHver2v1", "EleRunHver3v1"]
    
    if isMuMC:
        for samp in range(len(mc)):
    	    execme("hadd -k "+str(mc[samp])+"_Merged.root "+str(mc[samp])+"_MuMC*")
    
    if isEleMC:
        for samp in range(len(mc)):
    	    execme("hadd -k "+str(mc[samp])+"_Merged.root "+str(mc[samp])+"_EleMC*")
    
    if(isMuData):
        for samp in range(len(data)):
            execme("hadd -k "+str(data[samp])+"_Merged.root "+str(data[samp])+"*")
        execme("hadd -k all_Data.root MuRun*_Merged.root")
    
    if(isEleData):
        for samp in range(len(data)):
            execme("hadd -k "+str(data[samp])+"_Merged.root "+str(data[samp])+"*")
        execme("hadd -k all_Data.root EleRun*_Merged.root")
    
    execme("hadd -k all_TT.root TT_Merged.root")
    execme("hadd -k all_VV.root WW_Merged.root WZ_Merged.root ZZ_Merged.root")
    execme("hadd -k all_WJets.root  WJ*_Merged.root")
    #execme("hadd -k all_WJets.root W1*_Merged.root W2*_Merged.root W3*_Merged.root W4*_Merged.root WJ*_Merged.root")
    execme("hadd -k all_DY_M50.root DYJ*_Merged.root")
    execme("hadd -k all_Bkg.root all_DY_M50.root all_TT.root all_WJets.root all_VV.root")
    if(isMuMC):
        execme("hadd -k all_ExLep_M250.root ExLepMuMuZ_M250_Merged.root")
        execme("hadd -k all_ExLep_M500.root ExLepMuMuZ_M500_Merged.root")
        execme("hadd -k all_ExLep_M750.root ExLepMuMuZ_M750_Merged.root")
        execme("hadd -k all_ExLep_M1000.root ExLepMuMuZ_M1000_Merged.root")
        execme("hadd -k all_ExLep_M1250.root ExLepMuMuZ_M1250_Merged.root")
        execme("hadd -k all_ExLep_M1500.root ExLepMuMuZ_M1500_Merged.root")
        execme("hadd -k all_ExLep_M1750.root ExLepMuMuZ_M1750_Merged.root")
        execme("hadd -k all_ExLep_M2000.root ExLepMuMuZ_M2000_Merged.root")
        execme("hadd -k all_ExLep_M2500.root ExLepMuMuZ_M2500_Merged.root")
        execme("hadd -k all_ExLep_M3000.root ExLepMuMuZ_M3000_Merged.root")
        execme("hadd -k all_ExLep_M3500.root ExLepMuMuZ_M3500_Merged.root")
        execme("hadd -k all_ExLep_M4000.root ExLepMuMuZ_M4000_Merged.root")
        execme("hadd -k all_ExLep_M4500.root ExLepMuMuZ_M4500_Merged.root")
        execme("hadd -k all_ExLep_M5000.root ExLepMuMuZ_M5000_Merged.root")
    else:
        execme("hadd -k all_ExLep_M250.root ExLepEEZ_M250_Merged.root")
        execme("hadd -k all_ExLep_M500.root ExLepEEZ_M500_Merged.root")
        execme("hadd -k all_ExLep_M750.root ExLepEEZ_M750_Merged.root")
        execme("hadd -k all_ExLep_M1000.root ExLepEEZ_M1000_Merged.root")
        execme("hadd -k all_ExLep_M1250.root ExLepEEZ_M1250_Merged.root")
        execme("hadd -k all_ExLep_M1500.root ExLepEEZ_M1500_Merged.root")
        execme("hadd -k all_ExLep_M1750.root ExLepEEZ_M1750_Merged.root")
        execme("hadd -k all_ExLep_M2000.root ExLepEEZ_M2000_Merged.root")
        execme("hadd -k all_ExLep_M2500.root ExLepEEZ_M2500_Merged.root")
        execme("hadd -k all_ExLep_M3000.root ExLepEEZ_M3000_Merged.root")
        execme("hadd -k all_ExLep_M3500.root ExLepEEZ_M3500_Merged.root")
        execme("hadd -k all_ExLep_M4000.root ExLepEEZ_M4000_Merged.root")
        execme("hadd -k all_ExLep_M4500.root ExLepEEZ_M4500_Merged.root")
        execme("hadd -k all_ExLep_M5000.root ExLepEEZ_M5000_Merged.root")


