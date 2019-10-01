#!/bin/bash
#REFERENCE
#https://github.com/florez/CONDOR

#//////////////////////////////////////
#                                     #
# Cond.sub runs Analyzer_TIFR.sh #
# Analyzer_TIFR.sh runs runMe.sh      #
# runMe.sh runs Analyzer.C       #
#                                     #
#//////////////////////////////////////

#------------------------------------------------
#create a directory where all the outputs will be
#stored, for different merged ntuple input files
#------------------------------------------------
ntupleT2Paths=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
#voms-proxy-init --voms cms

log="log_"
logDir=$log${ntupleT2Paths/.txt/""}
baseDir="/home/rverma/t3store3/AN-18-126/Analyze2016Data/CondorLog/"
mkdir -p $baseDir$logDir
outcond="$baseDir$logDir"

cp Cond.sub $outcond
cp Analyzer_TIFR.sh $outcond
cp $ntupleT2Paths $outcond
cd $outcond

#------------------------------------------------
#read the file, where paths of ntuples are stored
#do not put empty lines in mergedNtupleT2.txt
#------------------------------------------------
count=0
cat $ntupleT2Paths | while read ntupleT2Path
do
  #----------------------------------------------
  #print T2Paths of ntuple, on terminal
  #split the T2Paths of ntuples, into an array
  #get the  second last entry of the array
  #remove .root, from the input ntuple
  #----------------------------------------------
  ((count++))
  echo -e "\033[01;32m input ntuple=\033[00m" $count": " $ntupleT2Path
  IFS='/ ' read -r -a array <<< "$ntupleT2Path"
  len=${#array[@]}
  fifth_last=`expr $len - 4`
  sec_last=`expr $len - 1`
  ntuple=${array[$sec_last]}
  ##ntuple=${array[$fifth_last]}${array[$sec_last]}
  #echo $ntuple
  iFile=${ntuple/.root/""}
  #----------------------------------------------
  #copy condor scripts to each input ntuple dir
  #replace Cond.sub arguments, as per input
  #submit the condor jobs, for each ntuple
  #----------------------------------------------
  
  mkdir -p $iFile
  cp Cond.sub $iFile
  cp Analyzer_TIFR.sh $iFile
  cd $iFile 
  sed -i "s:FNAME:$ntupleT2Path:g" Cond.sub
  sed -i "s:OUTPUTFILE:$iFile:g" Cond.sub
  sed -i "s:OUTPUTDIR:$iFile:g" Cond.sub
  condor_submit Cond.sub
  cd ../
done
