#!/bin/bash

#REFERENCE
#https://github.com/florez/CONDOR

#------------------------------------------------
#pass the arguments to Analysis/runMe.sh script 
#these arguments will go to the Analyzer.C
#------------------------------------------------

inNtupleFile=$1
outAnalFile=$2
outAnalDir=$3
date

#------------------------------------------------
#this script runs on some remote condor machine.
#link lxplus to this remote machine, using scram.
#copy the compiled lxplus package to this machine

#//////////// T3 /////////////////////////////////
echo "CONDOR DIR: $_CONDOR_SCRATCH_DIR"
cd ${_CONDOR_SCRATCH_DIR}
cp -r /home/rverma/t3store3/AN-18-126/Analyze2016Data/CMSSW_10_4_0/ .

#------------------------------------------------
#copy the lxplus package to the remote machine
#and run the codes at remote machine
#------------------------------------------------
cd CMSSW_10_4_0/src/ExLep2016Analysis
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
./runMe.sh $inNtupleFile $outAnalFile $outAnalDir

#---------------------------------------------
#copy the output from remote machine to the lxplus
#or to any other place e.g. Tier-2
#Remove the package, after copying the output
#------------------------------------------------
echo "OUTPUT: "
ls ${_CONDOR_SCRATCH_DIR}/CMSSW_10_4_0/src/ExLep2016Analysis/13TeV/$outAnalDir

cp -rf ${_CONDOR_SCRATCH_DIR}/CMSSW_10_4_0/src/ExLep2016Analysis/13TeV/$outAnalDir/* /home/rverma/t3store3/AN-18-126/Analyze2016Data/CondorOut/

cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_4_0
echo "DONE"
date

