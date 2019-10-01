#ifndef __SAMPLEINFO_H__
#define __SAMPLEINFO_H__

#include "TROOT.h"
#include <string>
#include <vector>

class SampleInfo
{
 public:
  SampleInfo();
  ~SampleInfo();
  
  std::string sampleName;
  std::vector<double>pileup;
  std::vector<double>truepileup;
  ///std::vector<double>puWeights;
  ///std::vector<double>truepuWeights;
  //Number of particles in an event
  int hepNUP;
  //particle pdg id
  std::vector<int>hepIDUP;
  //particle status
  std::vector<int>hepISTUP;
  std::vector<double>topPtWeights;
  double gen_weight;
 private :

};
#endif 
