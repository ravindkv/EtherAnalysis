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

#endif

class UncertaintyComputer{

public :
  UncertaintyComputer()
  {
  }

   virtual ~UncertaintyComputer(){
   ///~UncertaintyComputer(){
  }
  double getJERSF(double eta, int jer=0);
  double jetPtWithJESJER(MyJet jet, int jes=0, int jer=0); 
  void  openCSVfile(const std::string &filename); 
  double DeltaR(MyLorentzVector aV, MyLorentzVector bV);
  
private :
  ClassDef(UncertaintyComputer, 1)
};
#endif
