#ifndef __MYEVENT_H__
#define __MYEVENT_H__

/*

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aruna Kumar Nayak
//         Created:  Tue Aug 23 16:11:19 WEST 2011
// $Id: MyEvent.h,v 1.1 2012/01/05 11:12:21 anayak Exp $
//
//


// system include files
#include "TROOT.h"
#include <memory>
#include <string>
#include <vector>
#include <map>


// user include files
#include "TString.h"
#include "interface/MyJet.h"
#include "interface/MyMET.h"
#include "interface/MyVertex.h"
#include "interface/MyElectron.h"
#include "interface/MyMuon.h"
//#include "interface/MyTau.h"
#include "interface/SampleInfo.h"

class MyEvent
{
 public:
  MyEvent();
  ~MyEvent();
  
  enum SampleType {DATA = 0, 
      DY1JetsToLL=11, 
      DY2JetsToLL=12, 
      DY3JetsToLL=13, 
      DY4JetsToLL=14, 
      DYJetsToLL=10, 
      HplusM100=1,
      HplusM120=2,
      HplusM140=3,
      HplusM150=4,
      HplusM160=5,
      HplusM80=6,
      HplusM90=7,
      QCD_Pt_15to20=20,
      QCD_Pt_20to30=21,
      QCD_Pt_30to50=22,
      QCD_Pt_50to80=23,
      QCD_Pt_80to120=24,
      QCD_Pt_120to170=25,
      QCD_Pt_170to300=26,
      QCD_Pt_300to470=27,
      ST_s=31,
      ST_t=32,
      ST_tW=33,
      TTJets=30,
      W1JetsToLNu=41,
      W2JetsToLNu=42,
      W3JetsToLNu=43,
      W4JetsToLNu=44,
      WJetsToLNu=40,
      WW=50, 
      WZ=51, 
      ZZ=52, 
      OTHER=80};
  
  enum TTChannel {TTUNKNOWN=81, 
      TTALLJETS, 
      TTEJETS, 
      TTMUJETS, 
      TTTAUJETS, 
      TTEE, 
      TTEMU, 
      TTETAU, 
      TTMUMU, 
      TTMUTAU, 
      TTTAUTAU };
  
  enum ZChannel {ZUNKNOWN=61, 
      ZEE=62, 
      ZMUMU=63, 
      ZTAUTAU=64, 
      ZQQ=65};
  
  enum WChannel {WUNKNOWN=71, 
      WENU=72, 
      WMUNU=73, 
      WTAUNU=74, 
      WQQ=75};

  void Reset();
  
  // ---- General event information.
  unsigned int               runNb;
  unsigned int               eventNb;
  int                        lumiblock;
  
  std::vector<MyElectron> Electrons;
  std::vector<MyMuon> Muons;
  std::vector<MyJet> Jets;
//  std::vector<MyTau> Taus;
  std::vector<MyMET> mets;
  std::vector<MyVertex> PrimaryVtxs;
  SampleInfo sampleInfo;
  
  // ---- Monte Carlo information
  bool isData;
  // ---- Trigger info
  std::vector<std::string> hlt;

  int eventQuality;

 private :
  
};
#endif
