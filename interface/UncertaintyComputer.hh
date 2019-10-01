#ifndef _uncertaintycomputer_h_
#define _uncertaintycomputer_h_

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TTimeStamp.h"
#include <exception>

#ifdef _STANDALONE
#include "Reader.h"
#else
#include "interface/Reader.h"
#endif
#include "interface/BTagSF.hh"
#include "BTagCalibrationStandalone.h"

#endif

class UncertaintyComputer{

public :
  UncertaintyComputer()
  {
    btsf = new BTagSF(12345);
  }

   virtual ~UncertaintyComputer(){
   ///~UncertaintyComputer(){
     delete btsf;
  }
  double getJERSF(double eta, int jer=0);
  double jetPtWithJESJER(MyJet jet, int jes=0, int jer=0); 
  void  openCSVfile(const std::string &filename); 
  double DeltaR(MyLorentzVector aV, MyLorentzVector bV);
  
  // bTag SF, by event reweighting
  double getBTagPmcSys(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet);
  double getBTagPdataSys(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, int scale);
  
private :
  BTagSF* btsf;
  ClassDef(UncertaintyComputer, 1)
};
#endif
