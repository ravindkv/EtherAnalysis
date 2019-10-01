
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
    isMuData = False
    isMuMC = False
    isEleData = True
    isEleMC = True
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

    if(isMuMC):
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
    else:
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

    if(isMuData):
    	data = ["MuRunB2v2", "MuRunCv1", "MuRunDv1", "MuRunEv1", "MuRunFv1", "MuRunGv1", "MuRunH2v1", "MuRunH3v1"]
    
    if(isEleData): 
    	data = ["EleRunBver2v2", "EleRunCv1", "EleRunDv1", "EleRunEv1", "EleRunFv1", "EleRunGv1", "EleRunHver2v1"]
    
    
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
    
    execme("hadd -k all_DY_Pt.root DY*_Pt*Merged.root")
    execme("hadd -k all_DY_M50.root DYJetsToLL_M50_Merged.root")
    execme("hadd -k all_TT.root TT_Merged.root")
    execme("hadd -k all_WJets.root  WJ*_Merged.root")
    execme("hadd -k all_VV.root WW_Merged.root WZ_Merged.root ZZ_Merged.root")
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

if(is2017Data):
    #USERS INPUTS
    isMuData = True
    isMuMC = True
    isEleData = False
    isEleMC = False
    
    mc = []
    mc.append("TTJetsP")
    mc.append("DYJetsToLL")
    mc.append("WW")
    mc.append("WZ")
    mc.append("ZZ")
    
    if(isMuMC):
        mc.append("QCD_Pt-1000toInf_Mu")  
        mc.append("QCD_Pt-120to170_Mu")   
        mc.append("QCD_Pt-15to20_Mu")      
        mc.append("QCD_Pt-170to300_Mu")   
        mc.append("QCD_Pt-20to30_Mu")     
        mc.append("QCD_Pt-300to470_Mu")   
        mc.append("QCD_Pt-30to50_Mu")     
        mc.append("QCD_Pt-470to600_Mu")   
        mc.append("QCD_Pt-50to80_Mu")     
        mc.append("QCD_Pt-600to800_Mu")   
        mc.append("QCD_Pt-800to1000_Mu")  
        mc.append("QCD_Pt-80to120_Mu")    
    else:
        mc.append("QCD_Pt-15to20_EM")
        mc.append("QCD_Pt-20to30_EM")
        mc.append("QCD_Pt-30to50_EM")
        mc.append("QCD_Pt-50to80_EM")
        mc.append("QCD_Pt-80to120_EM")
        mc.append("QCD_Pt-120to170_EM")
        mc.append("QCD_Pt-170to300_EM")
        mc.append("QCD_Pt-300toInf_EM")
    
    
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
    
    execme("hadd -k all_TTJetsP.root TTJetsP_Merged.root")
    execme("hadd -k all_QCD.root QCD*_Merged.root")
    execme("hadd -k all_VV.root WW_Merged.root WZ_Merged.root ZZ_Merged.root")
    #execme("hadd -k all_WJets.root W1*_Merged.root W2*_Merged.root W3*_Merged.root W4*_Merged.root WJ*_Merged.root")
    execme("hadd -k all_DY.root DYJ*_Merged.root")
    execme("hadd -k all_ExLepMuMuZ_M250.root ExLepMuMuZ_M250_Merged.root")
    execme("hadd -k all_ExLepMuMuZ_M1500.root ExLepMuMuZ_M1500_Merged.root")
    execme("hadd -k all_ExLepMuMuZ_M2000.root ExLepMuMuZ_M2000_Merged.root")
    execme("hadd -k all_ExLepMuMuZ_M2500.root ExLepMuMuZ_M2500_Merged.root")
    execme("hadd -k all_ExLepMuMuZ_M4000.root ExLepMuMuZ_M4000_Merged.root")

