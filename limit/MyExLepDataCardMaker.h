#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> 

using namespace std;
double totLumi = 35.9;
//TString inFileDir="stack_for2016Data_20190117_Ele_Sys";
TString inFileDir="stack_for2016Data_20190117_Mu_Sys";
bool isMuChannel = true;
bool isEleChannel = false;
TFile* fTT      = TFile::Open(inFileDir+"/all_TT.root");

class MyExLepDataCardMaker{
  public:
  TH1F* getHisto(TFile *inRootFile, TString histPath, TString histName, TFile* fTT);
  TH1F* readHisto(TFile *inFile, TString histDir, TString histSubDir, TString inHistName, double sf, TString process);
  double getSysUnc(TH1F *hCentral, TH1F* hUp, TH1F* hDown);
  double getStatUnc(TH1F* hCentral, double sError = 0.0);
  private:
  double dont_use;  
};

//----------------------------------------//
//Variuos functions
//----------------------------------------//
TH1F*  MyExLepDataCardMaker:: getHisto(TFile *inRootFile, TString histPath, TString histName, TFile *fTT){
  TH1F* hist;
  TString fullPath = histPath+histName;
  string exception_msg (inRootFile->GetName()+TString("/")+fullPath+", does not exist");
  try{
    if(!(inRootFile->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    ///cout<<"WARNING:"<<e<<endl;
  }
  try{
    if(!(fTT->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    cout<<"\033[01;31mERROR: \033[0m"<<e<< endl;
    exit(0);
  }
  if(!(inRootFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(inRootFile->Get(fullPath))->Clone(histName);
  return hist;
}

//Read histos from input file. Return empty hist if the hist does not exist. 
TH1F* MyExLepDataCardMaker:: readHisto(TFile *inFile, TString histDir, TString histSubDir, TString inHistName, double sf, TString process){
  TH1F* hist;
  TString histPath = histDir+histSubDir;
  if(!(inFile->Get(histPath+inHistName))){
    hist = (TH1F*)(fTT->Get("base/Iso/ControlP/mll"));//initialise an empty hist
    hist->Reset();
  }else hist = (TH1F*)(inFile->Get(histPath+inHistName));
  hist->Scale(sf);
  cout<<setw(10)<<process<<setw(10)<<""<<histPath<<setw(15)<<inHistName<<setw(15)<<hist->Integral()<<endl;
  return hist;
}  

//get normalised uncertainity
double MyExLepDataCardMaker:: getSysUnc(TH1F *hCentral, TH1F* hUp, TH1F* hDown){
  return 1 + max(fabs(hUp->Integral() - hCentral->Integral()), fabs(hCentral->Integral() - hDown->Integral()))/hCentral->Integral();
}

//get statistical uncertainity
double MyExLepDataCardMaker:: getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00; 
  return statUnc;
}
