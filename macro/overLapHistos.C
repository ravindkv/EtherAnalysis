#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <cstring>

#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"
#include "TTree.h"

using namespace std;


//-------------------------------------------//
///INPUT FILES
//-------------------------------------------//
//data
TFile* fData    = TFile::Open("all_muData.root");
//bkg
TFile* fVV      = TFile::Open("all_VV.root");
TFile* fDY      = TFile::Open("all_DY.root");
TFile* fWJ      = TFile::Open("all_WJets.root");
TFile* fQCD     = TFile::Open("all_QCD.root");
TFile* fST      = TFile::Open("all_ST.root");
TFile* fTT      = TFile::Open("all_TTJetsP.root");
//signal
TFile *fWH80      = TFile::Open("all_Hplus80.root");
TFile *fWH90      = TFile::Open("all_Hplus90.root");
TFile *fWH100     = TFile::Open("all_Hplus100.root");
TFile *fWH120     = TFile::Open("all_Hplus120.root");
TFile *fWH140     = TFile::Open("all_Hplus140.root");
TFile *fWH150     = TFile::Open("all_Hplus150.root");
TFile *fWH155     = TFile::Open("all_Hplus155.root");
TFile *fWH160     = TFile::Open("all_Hplus160.root");
//data driven qcd
TFile* fQCD_dd = TFile::Open("all_QCD_dd.root"); 

//histogram name and range
TString histName = "pt_bjet" ;
Float_t xMin_ = 0.0 ;
Float_t xMax_ = 1000.0 ;


//-------------------------------------------//
//function to get histogram from root file
//-------------------------------------------//
TH1F* getHisto(TString file, TString histName, int histColor, Float_t xMin_ = 0.0, Float_t xMax_ = 500){

  TFile f(file,"READ");
  if(f.IsZombie() || f.IsZombie()){
    cout << "Cannot open file "<< endl;
    continue;
  }
  TH1F* h;
  h = (TH1F*)f.Get("base/Iso/BTag/"+histName);
  h->SetMarkerColor(kRed);
  h->SetMarkerStyle(1);
  h->SetMarkerSize(1.5);
  h->SetLineColor(histColor);
  h->SetLineWidth(3);

  h->GetXaxis()->SetLimits(xMin_, xMax_);
  h->GetYaxis()->SetTitleOffset(1.33);
  //h->SetMinimum(yMin_);
  //h->SetMaximum(yMax_);
  h->GetXaxis()->SetTitle(histName);
  h->GetYaxis()->SetTitle("Events");
  return h;
}

//-------------------------------------------//
//overlap histograms
//-------------------------------------------//
void overLapHistos( ){
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  
  TLegend* leg = new TLegend(0.19,0.60,0.46,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  
  TH1F* hTT = getHisto(fTT, histName, 1, xMin_, xMax_);
  TH1F* hWH80 = getHisto(fWH80, histName, 2, xMin_, xMax_);
  
  g11->Draw("ALP");
  g21->Draw("LPsame");
  g31->Draw("LPsame");
  g41->Draw("LPsame");
  
  c1->cd();
  gPad->Modified();
  //leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=13 TeV}{35.5 fb^{-1}}"));
  
  leg->AddEntry(g11,"cutSet1_CTagInc","PL");
  leg->AddEntry(g21,"cutSet2_CTagInc","PL");
  leg->AddEntry(g31,"cutSet3_CTagInc","PL");
  leg->AddEntry(g41,"cutSet4_CTagInc","PL");
  
  leg->Draw();
  gPad->RedrawAxis();
}
