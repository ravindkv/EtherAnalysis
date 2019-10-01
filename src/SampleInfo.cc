#include "interface/SampleInfo.h"

SampleInfo::SampleInfo()
{
  sampleName="";
  pileup.clear();
  truepileup.clear();
  ///puWeights.clear();
  ///truepuWeights.clear();
  hepNUP = -99;
  hepIDUP.clear();
  hepISTUP.clear();
  topPtWeights.clear();
  gen_weight = 1.0;
}

SampleInfo::~SampleInfo()
{
}
