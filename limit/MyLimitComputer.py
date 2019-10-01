#//////////////////////////////////////////////////
#                                                 #
# 	Limit computation at 13 TeV 		  #
#                                                 #
#//////////////////////////////////////////////////
# python MyLimitComputer.py --ch mu  --mass 250

import os
import sys
import argparse

def execme(command):
    #print "\033[01;32m"+ "Excecuting: "+ "\033[00m",  command
    os.system(command)

#---------------------------------------------
#function to prepare data cards
#---------------------------------------------
def makeDataCards(IN_FILE_DIR, HIST_DIR, HIST_NAME, CHANNEL_NAME, SIG_MASS, SIG_LABEL, SIG_FILE):
    execme('root -l -q -b \"MyExLepDataCardMaker.C(\\\"'+IN_FILE_DIR+'\\\", \\\"'+HIST_DIR+'\\\",\\\"'+HIST_NAME+'\\\",\\\"'+CHANNEL_NAME+'\\\",'+str(SIG_MASS)+', \\\"'+SIG_LABEL+'\\\", \\\"'+SIG_FILE+'\\\")\"')

def moveDataCards(CHANNEL_NAME, HIST_ARRAY, MASS, LIMIT_DIR):
    execme('mkdir -p '+LIMIT_DIR)
    execme('mv *.root '+LIMIT_DIR)
    execme('mv *datacard_* '+LIMIT_DIR)

#---------------------------------------------
#function to calculate limits
#---------------------------------------------
def calcLimits(CHANNEL_NAME, COMB_DATACARD_NAME, MASS, isGOF):
    t2wDataCardName = "t2w_"+COMB_DATACARD_NAME.replace(".txt",".root")
    execme('text2workspace.py '+COMB_DATACARD_NAME+' -o '+t2wDataCardName )
    execme('combine --rAbsAcc 0.000001 '+t2wDataCardName+' -M AsymptoticLimits --mass '+str(MASS)+' --name _llstar_llZ_llq_13TeV_'+CHANNEL_NAME)
    if(isGOF):
        #execme('combine '+t2wDataCardName+' -M GoodnessOfFit --algo saturated -t 1000 -s -1 --mass '+str(MASS)+' --name _llstar_llZ_llq_13TeV_toy_'+CHANNEL_NAME)
        execme('combine '+t2wDataCardName+' -M GoodnessOfFit --algo saturated --mass '+str(MASS)+' --name _llstar_llZ_llq_13TeV_data_'+CHANNEL_NAME)

#---------------------------------------------------
#function to get datacards to be combined
#---------------------------------------------------
def getCardsToBeCombined(CHANNEL_ARRAY, IN_FILE_DIR_ARRAY, HIST_ARRAY, MASS):
    #make separate cards first
    for CH in range(len(CHANNEL_ARRAY)):
        for HIST in range(len(HIST_ARRAY)):
            SIG_LABEL = ""
            #SIG_LABEL = "ExLep"+str(MASS)
            SIG_FILE = "all_ExLep_M"+str(MASS)+".root"
            #makeDataCards(IN_FILE_DIR_ARRAY[CH], HIST_ARRAY[HIST][0],HIST_ARRAY[HIST][1], CHANNEL_ARRAY[CH], MASS, SIG_LABEL, SIG_FILE)
            makeDataCards(IN_FILE_DIR_ARRAY[CH], HIST_ARRAY[HIST][0],HIST_ARRAY[HIST][1]+str(MASS), CHANNEL_ARRAY[CH], MASS, SIG_LABEL, SIG_FILE)
    #store separate cards in an array
    COMB_CARD_CHANNEL_HIST_MASS = []
    for CH in range(len(CHANNEL_ARRAY)):
        COMB_CARD_HIST_MASS = []
        for HIST in range(len(HIST_ARRAY)):
            COMB_CARD_MASS = []
            COMB_CARD_MASS.append('datacard_llstar_llZ_llq_13TeV_'+CHANNEL_ARRAY[CH]+'_'+HIST_ARRAY[HIST][0]+'_'+HIST_ARRAY[HIST][1]+str(MASS)+'.txt')
            #COMB_CARD_MASS.append('datacard_llstar_llZ_llq_13TeV_'+CHANNEL_ARRAY[CH]+'_'+HIST_ARRAY[HIST][0]+'_'+HIST_ARRAY[HIST][1]+'_ExLep'+str(MASS)+'.txt')
            COMB_CARD_HIST_MASS.append(COMB_CARD_MASS)
        COMB_CARD_CHANNEL_HIST_MASS.append(COMB_CARD_HIST_MASS)    
    return COMB_CARD_CHANNEL_HIST_MASS

#---------------------------------------------------
#function to arrange datacards for combined limits
#---------------------------------------------------
def sortCardsForCombine(COMB_CARD_CHANNEL_HIST_MASS_ARRAY, CHANNEL_ARRAY, HIST_ARRAY, MASS):
    SORT_CARD = ' '
    COMB_CARD_MASS = []
    for CH in range(len(CHANNEL_ARRAY)):
        for HIST in range(len(HIST_ARRAY)):
            COMB_CARD_MASS.append(COMB_CARD_CHANNEL_HIST_MASS_ARRAY[CH][HIST][0])
    for STR in COMB_CARD_MASS:
        SORT_CARD = SORT_CARD+STR+' '
    return SORT_CARD

#---------------------------------------------------
#function to calculate combined limits
#---------------------------------------------------
def calcCombinedLimit(CHANNEL_ARRAY, IN_FILE_DIR_ARRAY, HIST_ARRAY,  MASS_ARRAY, isCondSub, isGOF):
    #print IN_FILE_DIR_ARRAY, HIST_ARRAY, MASS_ARRAY
    COMB_CHANNEL_NAME = '_'.join(CHANNEL_ARRAY)
    HIST_ARRAY_ = []
    for HIST in range(len(HIST_ARRAY)):
        HIST_ARRAY_.append(HIST_ARRAY[HIST][1])
    COMB_HIST_NAME = '_'.join(HIST_ARRAY_)
    for MASS in range(len(MASS_ARRAY)):
        LIMIT_DIR = "output/"+COMB_CHANNEL_NAME+"/"+HIST_ARRAY[0][0]+"/Mass"+str(MASS_ARRAY[MASS])
        if not isCondSub: execme('mkdir -p '+LIMIT_DIR)
        getCardsToBeCombined_ = getCardsToBeCombined(CHANNEL_ARRAY, IN_FILE_DIR_ARRAY, HIST_ARRAY, MASS_ARRAY[MASS])
        sortCardsForCombine_ = sortCardsForCombine(getCardsToBeCombined_, CHANNEL_ARRAY, HIST_ARRAY, MASS)
        print sortCardsForCombine_
        COMB_DATACARD_NAME = 'combine_datacard_llstar_llZ_llq_13TeV_'+COMB_CHANNEL_NAME+"_ExLep"+str(MASS_ARRAY[MASS])+".txt"
        if len(CHANNEL_ARRAY)>1 or len(HIST_ARRAY)>1:
            execme('combineCards.py '+sortCardsForCombine_+' > '+COMB_DATACARD_NAME)
        else: execme('cp '+sortCardsForCombine_+' '+COMB_DATACARD_NAME)
        calcLimits(COMB_CHANNEL_NAME, COMB_DATACARD_NAME, MASS_ARRAY[MASS], isGOF)
        if not isCondSub: moveDataCards(COMB_CHANNEL_NAME, HIST_ARRAY, MASS_ARRAY[MASS], LIMIT_DIR)


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--ch", default="mu", help="The channel name e.g. mu or ele or mu_ele")
    parser.add_argument("--cat", default=1, help="Type of event category")
    parser.add_argument("--mass", default=250, help="Mass of the excited lepton")
    parser.add_argument("--allMass", default=False, help="Mass of the excited lepton")
    parser.add_argument("--histDir",  default="ZTag1", help="histogram directory")
    parser.add_argument("--batch", default=False, help="Want to submit condor jobs")
    parser.add_argument("--isGOF", default=False, help="Want to determine the goodness-of-fit")
    args = parser.parse_args()
    
    #---------------------------------------------
    #USERS INPUTS
    #---------------------------------------------
    path_file_dir="/home/sthakur/t3store2/xRavindra/CMSSW_8_1_0/src/HiggsAnalysis/LimitCompExLepton/for2016Data/Limits/"
    muon_file_dir=path_file_dir+"stack_tmp_mu_new"
    ele_file_dir=path_file_dir+"stack_tmp_ele_new"
    #muon_file_dir=path_file_dir+"stack_for2016Data_20190117_Mu_Sys"
    #ele_file_dir=path_file_dir+"stack_for2016Data_20190117_Ele_Sys"

    hist_array = []
    #hist_array.append([args.histDir, "mlZ_max_sig"])
    hist_array.append([args.histDir, "A_mlZ_max_sig"])

    #mass_array = [250, 500]
    mass_array = [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
    if(args.ch=="mu"):
        in_channel= ["mu"]
        in_file = [muon_file_dir]
    if(args.ch=="ele"):
        in_channel= ["ele"]
        in_file = [ele_file_dir]
    if(args.ch=="mu_ele"):
        in_channel= ["mu", "ele"]
        in_file = [muon_file_dir, ele_file_dir]
 
    if(args.allMass): in_mass = mass_array
    else: in_mass = [args.mass]
    calcCombinedLimit(in_channel, in_file, hist_array, in_mass, args.batch, args.isGOF)

