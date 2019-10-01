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

bool isMuChannel = true;
bool isEleChannel = false;
bool isSaveHisto  = true;
//-------------------------------------------//
///INPUT FILES
//-------------------------------------------//
//data
/*
TFile* fData    = TFile::Open("all_muData.root");
//bkg
TFile* fVV      = TFile::Open("all_VV.root");
TFile* fDY      = TFile::Open("all_DY.root");
TFile* fWJ      = TFile::Open("all_WJets.root");
TFile* fQCD     = TFile::Open("all_QCD.root");
TFile* fST      = TFile::Open("all_ST.root");
*/
TFile* fTT      = TFile::Open("all_ExLep_M250.root");

//signal
TFile *fM250      = TFile::Open("all_ExLep_M250.root");
TFile *fM500      = TFile::Open("all_ExLep_M500.root");
TFile *fM750      = TFile::Open("all_ExLep_M750.root");
TFile *fM1000     = TFile::Open("all_ExLep_M1000.root");
TFile *fM1250     = TFile::Open("all_ExLep_M1250.root");
TFile *fM1500     = TFile::Open("all_ExLep_M1500.root");
TFile *fM1750     = TFile::Open("all_ExLep_M1750.root");
TFile *fM2000     = TFile::Open("all_ExLep_M2000.root");
TFile *fM2500     = TFile::Open("all_ExLep_M2500.root");
TFile *fM3000     = TFile::Open("all_ExLep_M3000.root");
TFile *fM3500     = TFile::Open("all_ExLep_M3500.root");
TFile *fM4000     = TFile::Open("all_ExLep_M4000.root");
TFile *fM4500     = TFile::Open("all_ExLep_M4500.root");
TFile *fM5000     = TFile::Open("all_ExLep_M5000.root");
bool isSM250sto = true;

//-------------------------------------------//
//function to get histogram from root file
//-------------------------------------------//
TH1F* getHisto(TFile *file, TString histName, int histColor, Float_t xMin_ = 0.0, Float_t xMax_ = 500){

  if(file->IsZombie() || file->IsZombie()){
    cout << "Cannot open file "<< endl;
  }
  TH1F* h;
  h = (TH1F*)file->Get("base/Iso/ZTag/"+histName);
  //h->SetMarkerColor(kRed);
  cout<<h->GetNbinsX()<<endl;
  //h->Rebin();
  cout<<h->GetNbinsX()<<endl;
  h->SetMarkerStyle(1);
  h->SetMarkerSize(1.5);
  h->SetLineColor(histColor);
  h->SetLineWidth(3);
  cout<<histName<<": Int = "<<h->Integral()<<", Mean = "<<h->GetMean()<<", RMS = "<<h->GetRMS()<<endl;
  h->GetXaxis()->SetRangeUser(xMin_, xMax_);
  h->GetYaxis()->SetRangeUser(0, 10);
  //h->GetXaxis()->SetTitle(histName);
  h->GetXaxis()->SetTitle("mass of fat-jet (GeV)");
  h->GetYaxis()->SetTitle("Events");
  h->GetYaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleSize(0.05); 
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05); 
  h->GetXaxis()->LabelsOption("u"); // extra
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetNdivisions(5);
  h->GetYaxis()->SetNdivisions(5);
  h->SetTitle("");
  return h;
}
TPaveText * paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font ){
  TPaveText *pt = new TPaveText(minX, minY, maxX, maxY, "brNDC"); // good_v1
  pt->SetBorderSize(size);
  pt->SetFillColor(fillColor);
  pt->SetFillStyle(style);
  pt->SetLineColor(lineColor);
  pt->SetTextFont(font);
  return pt;
}

void MyOverlayHist(){
  TString histName = "mass_fatjet";
  gStyle->SetFrameLineWidth(2);
  //pave text CMS box
  TPaveText *pt = new TPaveText(0.60,0.9054,0.90,0.9362, "brNDC"); 
  pt->SetBorderSize(1); pt->SetFillColor(19);
  pt->SetFillStyle(0); pt->SetTextSize(0.05);
  pt->SetLineColor(0); pt->SetTextFont(132);
  TText *text = pt->AddText("35.9 fb^{-1} (13 TeV)");
  text->SetTextAlign(11);
  
  //channel
  TPaveText *ch = paveText(0.25,0.80,0.30,0.85, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.05);
  if(isMuChannel) ch->AddText("#splitline{#mu + jets}{""}");
  if(isEleChannel) ch->AddText("#splitline{e + jets}{""}");
     
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas(histName, histName, 800, 800);
  TLegend* leg = new TLegend(0.65,0.40,0.80,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);
  
  TGaxis::SetMaxDigits(3);
  //overlay
  c1->cd();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13); 
  //gPad->SetTopMargin(0.11);
  double xMin_ = 0.0;
  double xMax_ = 400;
 
  //hWH80->GetYaxis()->SetRangeUser(1, 1.2* hWH80->GetMaximum());
  TH1F* hM250   = getHisto(fM250 , histName, 2, xMin_, xMax_);
  TH1F* hM500   = getHisto(fM500 , histName, 3, xMin_, xMax_);
  TH1F* hM750   = getHisto(fM750 , histName, 4, xMin_, xMax_);
  TH1F* hM1000  = getHisto(fM1000, histName, 52, xMin_, xMax_);
  TH1F* hM1250  = getHisto(fM1250, histName, 6, xMin_, xMax_);
  TH1F* hM1500  = getHisto(fM1500, histName, 7, xMin_, xMax_);
  TH1F* hM1750  = getHisto(fM1750, histName, 8, xMin_, xMax_);
  TH1F* hM2000  = getHisto(fM2000, histName, 11, xMin_, xMax_);
  /*
  TH1F* hM2500  = getHisto(fM2500, histName, 11, xMin_, xMax_);
  TH1F* hM3000  = getHisto(fM3000, histName, 11, xMin_, xMax_);
  TH1F* hM3500  = getHisto(fM3500, histName, 11, xMin_, xMax_);
  TH1F* hM4000  = getHisto(fM4000, histName, 11, xMin_, xMax_);
  TH1F* hM4500  = getHisto(fM4500, histName, 11, xMin_, xMax_);
  TH1F* hM5000  = getHisto(fM5000, histName, 11, xMin_, xMax_);
  */
  //hTT->Draw("HIST");
  hM250 ->Draw("HIST");
  hM500 ->Draw("HISTSAME");
  hM750 ->Draw("HISTSAME");
  hM1000->Draw("HISTSAME");
  hM1250->Draw("HISTSAME");
  hM1500->Draw("HISTSAME");
  hM1750->Draw("HISTSAME");
  hM2000->Draw("HISTSAME");
  /*
  hM2500->Draw("HISTSAME");
  hM3000->Draw("HISTSAME");
  hM3500->Draw("HISTSAME");
  hM4000->Draw("HISTSAME");
  hM4500->Draw("HISTSAME");
  hM5000->Draw("HISTSAME");
  */
  //leg->AddEntry(hTT, "ttbar","PL");
  leg->AddEntry(hM250 , "m_{l*} = 250 GeV","PL");
  leg->AddEntry(hM500 , "m_{l*} = 500 GeV","PL");
  leg->AddEntry(hM750 , "m_{l*} = 750  GeV","PL");
  leg->AddEntry(hM1000, "m_{l*} = 1000 GeV","PL");
  leg->AddEntry(hM1250, "m_{l*} = 1250 GeV","PL");
  leg->AddEntry(hM1500, "m_{l*} = 1500 GeV","PL");
  leg->AddEntry(hM1750, "m_{l*} = 1750 GeV","PL");
  leg->AddEntry(hM2000, "m_{l*} = 2000 GeV","PL");
  /*
  leg->AddEntry(hM2500, "m_{l*} = 2500 GeV","PL");
  leg->AddEntry(hM3000, "m_{l*} = 3000 GeV","PL");
  leg->AddEntry(hM3500, "m_{l*} = 3500 GeV","PL");
  leg->AddEntry(hM4000, "m_{l*} = 4000 GeV","PL");
  leg->AddEntry(hM4500, "m_{l*} = 4500 GeV","PL");
  leg->AddEntry(hM5000, "m_{l*} = 5000 GeV","PL");
  */
  leg->Draw();
  ch->Draw();
  pt->Draw();
  if(isSaveHisto){
    TString outFile("$PWD/");
    if(isMuChannel) outFile += "sig_"+histName+"_mu.png";
    if(isEleChannel) outFile += "sig_"+histName+"_ele.pdf";
    c1->SaveAs(outFile);
    //c1->Close();
  }
}
