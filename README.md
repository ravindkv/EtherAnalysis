# Ether Analysis
  
### Use it in Non-CMSSSW area. Only ROOT6 is needed.
* git clone https://github.com/ravindkv/EtherAnalysis.git
* cd EtherAnalysis/src
Change complier to c++11 in line 11 of the Makefile
* make clean 
* make
* cd .. 

### Load shared libraries and run, in one go ### 
* root -l 'runMe.C("Analyzer")'

 
### Use it inside the CMSSW  ###  
* source /cvmfs/cms.cern.ch/cmsset_default.sh
* cmsrel CMSSW_10_4_0
* cd CMSSW_10_4_0/src
* cmsenv
* git clone https://github.com/ravindkv/EtherAnalysis.git
* cd EtherAnalysis/src
* make clean 
* make
* cd .. 

### Submit condor batch jobs  ###

* cd condor
* ./RunCond_TIFR.sh ntupleT2Paths.txt
