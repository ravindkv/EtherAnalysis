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

void LimitPlotter(TString CHANNEL="mu", TString zTagDir="ZTag3",
         bool obs= false, bool isOut= true )
  {
  gStyle->SetFrameLineWidth(2);
  TCanvas *c1 = new TCanvas();
  gPad->SetLogy();
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  float X[] = {250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
  float obsY[]      = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY[]      = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
		                                              
  float expX1sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expX1sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sL_[]  = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sH_[]  = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
                                                              
		                                              
  float expX2sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expX2sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY2sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY2sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};

  int nMassPoints = 14;
  TString ch_hist = CHANNEL;
  TString massFiles [14] = {
  "Mass250/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH250.root",
  "Mass500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH500.root",
  "Mass750/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH750.root",
  "Mass1000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1000.root",
  "Mass1250/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1250.root",
  "Mass1500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1500.root",
  "Mass1750/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1750.root",
  "Mass2000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH2000.root",
  "Mass2500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH2500.root",
  "Mass3000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH3000.root",
  "Mass3500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH3500.root",
  "Mass4000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH4000.root",
  "Mass4500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH4500.root",
  "Mass5000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH5000.root"};
  
  double maxY = 1.0;
  for(unsigned int i = 0 ; i < nMassPoints; i++){
    TFile f("output/"+CHANNEL+"/"+zTagDir+"/"+massFiles[i],"READ"); 
    if(f.IsZombie()){
      cout << "Cannot open file for " << string(CHANNEL.Data()) << " and mass " << X[i] << endl;
      continue;
    }
    vector<double> sigMass;
    sigMass.push_back(250); 
    sigMass.push_back(500); 
    sigMass.push_back(750); 
    sigMass.push_back(1000); 
    sigMass.push_back(1250); 
    sigMass.push_back(1500); 
    sigMass.push_back(1750); 
    sigMass.push_back(2000); 
    sigMass.push_back(2500); 
    sigMass.push_back(3000); 
    sigMass.push_back(3500); 
    sigMass.push_back(4000); 
    sigMass.push_back(4500); 
    sigMass.push_back(5000);
    vector<double>sigXss;
    if(CHANNEL=="mu"){
      sigXss.push_back(0.00427)    ;
      sigXss.push_back(0.00291)    ;
      sigXss.push_back(0.001761)   ;
      sigXss.push_back(0.001177)   ;
      sigXss.push_back(0.0007263)  ;
      sigXss.push_back(0.0004267)  ;
      sigXss.push_back(0.0002654)  ;
      sigXss.push_back(0.0002021)  ;
      sigXss.push_back(0.00006755) ;
      sigXss.push_back(0.00002223) ;
      sigXss.push_back(0.000008166);
      sigXss.push_back(0.000003209);
      sigXss.push_back(0.000001133);
      sigXss.push_back(0.000000452);
    }
    if(CHANNEL=="ele"){
      sigXss.push_back(0.004407)  ;
      sigXss.push_back(0.002916)  ;
      sigXss.push_back(0.001607)  ;
      sigXss.push_back(0.001079)  ;
      sigXss.push_back(0.0007171) ;
      sigXss.push_back(0.0004856) ;
      sigXss.push_back(0.0002907) ;
      sigXss.push_back(0.0001908) ;
      sigXss.push_back(0.00005731);
      sigXss.push_back(0.00002497);
      sigXss.push_back(0.000009023);
      sigXss.push_back(0.000003897);
      sigXss.push_back(0.000001243);
      sigXss.push_back(0.000000461);
    }
    if(CHANNEL=="mu_ele"){
      sigXss.push_back(0.00427     + 0.004407)  ;
      sigXss.push_back(0.00291     + 0.002916)  ;
      sigXss.push_back(0.001761    + 0.001607)  ;
      sigXss.push_back(0.001177    + 0.001079)  ;
      sigXss.push_back(0.0007263   + 0.0007171) ;
      sigXss.push_back(0.0004267   + 0.0004856) ;
      sigXss.push_back(0.0002654   + 0.0002907) ;
      sigXss.push_back(0.0002021   + 0.0001908) ;
      sigXss.push_back(0.00006755  + 0.00005731);
      sigXss.push_back(0.00002223  + 0.00002497);
      sigXss.push_back(0.000008166 + 0.000009023);
      sigXss.push_back(0.000003209 + 0.000003897);
      sigXss.push_back(0.000001133 + 0.000001243);
      sigXss.push_back(0.000000452 + 0.000000461);
    }

    Double_t r;
    TTree* limit = (TTree*)f.Get("limit");
    limit->SetBranchAddress("limit",&r);
    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);
      //multiply by xss
      r = 1000*r*sigXss[i]; //1000 to convert xss into fb
      if(k==0) expY2sL[i] = r;
      if(k==1) expY1sL[i] = r;
      if(k==1) expY1sL_[i] = r;
      if(k==2) expY[i]    = r;
      if(k==3) expY1sH[i] = r;
      if(k==3) expY1sH_[i] = r;
      if(k==4) expY2sH[i] = r;
      if(k==5) obsY[i]    = r;
    }
    maxY = expY2sH[0];
  }
  cout<<std::setprecision(4)<<endl;
  cout<<"Mass:"<<setw(15)<<"base value"<<setw(15)<<"-2 #sigma"<<setw(15)<<"-1 #sigma"<<setw(15)<<"+1 #sigma"<<setw(15)<<"+2 #sigma"<<endl; 
  for(int i1 = 0 ; i1 < nMassPoints ; i1++){
  cout<<X[i1]<<setw(15)<<expY[i1]<<setw(15)<<expY2sL[i1]<<setw(15)<< expY1sL[i1]<<setw(15)<<expY1sH[i1]<<setw(15)<<expY2sH[i1]<<endl; 
    expY1sH[i1] = TMath::Abs(expY1sH[i1]-expY[i1]);
    expY1sL[i1] = TMath::Abs(expY1sL[i1]-expY[i1]);
    expY2sH[i1] = TMath::Abs(expY2sH[i1]-expY[i1]);
    expY2sL[i1] = TMath::Abs(expY2sL[i1]-expY[i1]);
  //cout<<setw(10)<<X[i1]<<setw(15)<<expY[i1]<<setw(15)<<expY2sL[i1]<<setw(15)<< expY1sL[i1]<<setw(15)<<expY1sH[i1]<<setw(15)<<expY2sH[i1]<<endl; 
  }
  for(int i1 = 0 ; i1 < nMassPoints ; i1++){
    expY1sH_[i1] = 100*TMath::Abs(expY1sH_[i1]-expY[i1]);
    expY1sL_[i1] = 100*TMath::Abs(expY1sL_[i1]-expY[i1]);
    //cout<<"$"<<std::setprecision(2)<<100*expY[i1]<<"^{+"<<expY1sH_[i1]<<"}"<<"_"<<"{-"<< expY1sL_[i1]<<"}"<<"$"<<endl;
    //cout<<endl;
  }
  
  TMultiGraph *mg = new TMultiGraph();
  TString ch_name = "ch";
  cout<<CHANNEL<<endl;
  if(CHANNEL=="mu") ch_name = "#mu";
  if(CHANNEL=="ele") ch_name = "e";
  if(CHANNEL=="mu_ele") ch_name = "lep";
  mg->SetMaximum(1.5*maxY);

  TGraphAsymmErrors* expected = new TGraphAsymmErrors(nMassPoints, X, expY, expX1sL ,expX1sL , expX1sL, expX1sL);
  TGraphAsymmErrors* oneSigma = new TGraphAsymmErrors(nMassPoints, X, expY, expX1sL, expX1sL,  expY1sL, expY1sH);
  TGraphAsymmErrors* twoSigma = new TGraphAsymmErrors(nMassPoints, X, expY, expX2sL, expX2sL,  expY2sL, expY2sH);
  TGraphAsymmErrors* observed = new TGraphAsymmErrors(nMassPoints, X, obsY, expX1sL ,expX1sL , expX1sL, expX1sL);

 
  //oneSigma->SetMarkerColor(kBlack);
  //oneSigma->SetMarkerStyle(kFullCircle);
  oneSigma->SetFillColor(kGreen+1);
  oneSigma->SetFillStyle(1001);

  //twoSigma->SetMarkerColor(kBlack);
  //twoSigma->SetMarkerStyle(kFullCircle);
  twoSigma->SetFillColor(kOrange);
  twoSigma->SetFillStyle(1001);

  expected->SetMarkerColor(kBlack);
  expected->SetMarkerStyle(23);
  //expected->SetMarkerSize(0.0);
  expected->SetLineColor(kBlack);
  expected->SetLineWidth(2);

  observed->SetMarkerColor(kBlue);
  observed->SetMarkerStyle(20);
  observed->SetLineColor(kBlue);
  observed->SetLineStyle(2);
  observed->SetLineWidth(4);

  mg->Add(twoSigma);
  mg->Add(oneSigma);
  mg->Add(expected);
  mg->Draw("ALP3");
  if(obs) mg->Add(observed);

  gPad->Modified();
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  //gPad->SetGridy();
  //gPad->SetGridx();
  gPad->SetRightMargin(0.05);
  gStyle->SetFrameLineWidth(3);
  //mg->SetMinimum(1.0);
  //mg->SetMaximum(yMax);
  TString fullProcess = "";
  if(CHANNEL=="mu"){
    mg->GetXaxis()->SetTitle("M_{#mu*} (GeV)");
    fullProcess = "#mu#mu*#rightarrow #mu#muZ #rightarrow 2#mu2q";   
  }
  if(CHANNEL=="ele"){
    mg->GetXaxis()->SetTitle("M_{e*} (GeV)");
    fullProcess = "ee*#rightarrow eeZ #rightarrow 2e2q";   
  }
  if(CHANNEL=="mu_ele"){
    mg->GetXaxis()->SetTitle("M_{l*} (GeV)");
    fullProcess = "ll*#rightarrow llZ #rightarrow 2l2q";   
  }
  
  //mg->GetYaxis()->CenterTitle();
  mg->GetXaxis()->SetLimits(0,5500);
  mg->GetYaxis()->SetTitleOffset(1.20);
  mg->GetYaxis()->SetNdivisions(6);
  mg->GetXaxis()->SetNdivisions(11);
  mg->GetYaxis()->SetRangeUser(0.01, 100);
  mg->GetYaxis()->SetMoreLogLabels(true);
  mg->GetXaxis()->SetTitleOffset(1.15);
  mg->GetYaxis()->SetTitle("#sigma (fb) #times BR("+fullProcess+")");
  mg->GetYaxis()->SetTitleSize(0.06);   
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetLabelSize(0.04);   
  mg->GetYaxis()->SetLabelSize(0.05);   
  mg->GetXaxis()->SetTickLength(0.04);
  mg->GetYaxis()->SetTickLength(0.04);

  TLegend* leg = new TLegend(0.65,0.60,0.85,0.87,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);
  leg->SetHeader("95% CL upper limits");
  leg->AddEntry(expected,"Median expected","LP");
  if(obs) leg->AddEntry(observed,"Observed","LP");
  leg->AddEntry(oneSigma, "68% expected","F");
  leg->AddEntry(twoSigma, "95% expected","F");
  leg->Draw();

 TPaveText *pl2 = new TPaveText(0.64,0.67,0.75,0.87, "brNDC");
  pl2->SetTextSize(0.052);
  pl2->SetFillColor(0);
  pl2->SetTextFont(132);
  pl2->SetBorderSize(0);
  pl2->SetTextAlign(11);
  pl2->AddText("t #rightarrow H^{#pm}b,");
  pl2->AddText("H^{+} #rightarrow c#bar{s}");
  pl2->AddText("BR(H^{+} #rightarrow c#bar{s}) = 1");
 
  //pave text CMS box
  TPaveText *pt = new TPaveText(0.15,0.9354,0.65,0.9362, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.06);
  pt->SetLineColor(0);
  pt->SetTextFont(132);
  //TText *text = pt->AddText("#sqrt{s}=13 TeV, 160 fb^{-1} ");
  TText *text = pt->AddText("CMS Preliminary, 13 TeV (35.9 fb^{-1}) ");
  //TText *text = pt->AddText(dir+":  CMS Preliminary,    #sqrt{s} = 13 TeV,    35.45 fb^{-1}; ");
  text->SetTextAlign(11);
  
  //pave text channel box
  TPaveText *ch = new TPaveText(0.85,0.9154898,0.9010067,0.9762187,"brNDC");
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetTextSize(0.06);
  ch->SetBorderSize(1);
  if(CHANNEL=="mu")ch->AddText("#mu + jets");
  if(CHANNEL=="ele") ch->AddText("e + jets");
  if(CHANNEL=="mu_ele") ch->AddText("l + jets");
  //pl2->Draw("SAME");
  pt->Draw("SAME");
  ch->Draw("SAME");
  leg->Draw("SAME");
  
  //-------------------------------
  //Draw theoretical limit
  //-------------------------------
  TFile fTh("theoryXss.root","READ");
  TGraph *graphLamb500 = (TGraph*)(fTh.Get("lambda_500"));
  TGraph *graphLamb2000 = (TGraph*)(fTh.Get("lambda_2000"));
  TGraph *graphLamb4000 = (TGraph*)(fTh.Get("lambda_4000"));
  TGraph *graphLamb10000 = (TGraph*)(fTh.Get("lambda_10000"));
  TGraph *graphLambMlstar = (TGraph*)(fTh.Get("lambda_Mlstar"));
  graphLamb500->SetLineColor(kBlue);
  graphLamb2000->SetLineColor(kCyan);
  graphLamb4000->SetLineColor(kGray);
  graphLamb10000->SetLineColor(kViolet);
  graphLambMlstar->SetLineColor(kRed);
  graphLamb500->SetLineWidth(3);
  graphLamb2000->SetLineWidth(3);
  graphLamb4000->SetLineWidth(3);
  graphLamb10000->SetLineWidth(3);
  graphLambMlstar->SetLineWidth(4);
  //TMultiGraph *mg2 = new TMultiGraph();
  graphLamb500->Draw("SAME");
  graphLamb2000->Draw("SAME");
  graphLamb4000->Draw("SAME");
  graphLamb10000->Draw("SAME");
  graphLambMlstar->Draw("SAME");
  //graphLambMlstar->SetLineColor(2);
  //graphLambMlstar->SetLineWidth(1504);
  //graphLambMlstar->SetFillStyle(3005);
  //graphLambMlstar->SetFillColor(45);
  //graphLambMlstar->SetFillStyle(3004);
  //mg2->Add(graphLamb500);
  //mg2->Add(graphLamb2000);
  //mg2->Add(graphLamb4000);
  //mg2->Add(graphLamb10000);
  //mg2->Add(graphLambMlstar);
  //mg2->Draw("BC");

  TLegend* grapLeg = new TLegend(0.30,0.60,0.40,0.80,NULL,"brNDC");
  //grapLeg->SetNColumns(2);
  grapLeg->SetBorderSize(0);
  grapLeg->SetTextSize(0.03);
  grapLeg->SetFillColor(0);
  grapLeg->AddEntry(graphLamb500, "#Lambda = 0.5 TeV", "L");
  grapLeg->AddEntry(graphLamb2000, "#Lambda = 2 TeV", "L");
  grapLeg->AddEntry(graphLamb4000, "#Lambda = 4 TeV", "L");
  grapLeg->AddEntry(graphLamb10000, "#Lambda = 10 TeV", "L");
  grapLeg->AddEntry(graphLambMlstar, "#Lambda = M_{l^{*}}", "L");
  grapLeg->Draw();

  gPad->RedrawAxis();
  TString outFile = "output_"+CHANNEL;
  TString outDir = "output/"+CHANNEL;
  gPad->SaveAs(outDir+"/"+outFile+".png");
  if(isOut){
    TFile *fout = new TFile(outDir+"/"+outFile+".root", "RECREATE");
    expected->Write("expected");
    observed->Write("observed");
    oneSigma->Write("oneSigma");
    twoSigma->Write("twoSigma");
    fout->Close();
  }
}

void MyLimitPlotter(){
  LimitPlotter("mu", "ZTag3",  false, false);
  LimitPlotter("ele","ZTag1", false, false);
  //LimitPlotter("mu_ele", false, true );
}
