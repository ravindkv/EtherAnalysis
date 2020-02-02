#include "MyEtherStackHisto.h"

void plotStackedHisto(TString chDir, TString baseDir, TString histDir, TString histName, TString xTitle,bool isData=false, bool isSig=false, double xmin=0, double xmax=10, double unc=false, bool isDYdd=false){
  MyEtherStackHisto MyESHist;
  string hist_name (histName);
  //Pad
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(3);
  TCanvas *c1 = new TCanvas(histName, histName, 700, 700);
  //TCanvas *c1 = new TCanvas("c1", "Data_MC", 400, 600);
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.30,0.30,0.98};
  if(isData){
    //c1->Divide(1, 2); c1->cd(1);
    c1->Divide(1, 2, 0, 0); c1->cd(1);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    if(isData) gPad->SetLogy(true);
    if(hist_name.find("mjj") != string::npos) gPad->SetLogy(false);
  }
  
  bool isMuChannel = false;
  bool isEleChannel = false;
  if(chDir.Contains("Muon")) isMuChannel = true;
  if(chDir.Contains("Electron")) isEleChannel = true;
  //-------------------------------
  // stack MC Bkg histo
  //-------------------------------
  THStack* hStack = new THStack("hStack","");
  TH1F* hVV = MyESHist.getHisto(fVV, chDir, baseDir, histDir, histName);
  TH1F* hMC = (TH1F*)hVV->Clone("hMC");
  int col_depth =0;
  hVV->SetFillColor(kOrange + col_depth);
  hStack->Add(hVV);
  TH1F* hWJ = MyESHist.stackHisto(fWJ, chDir, baseDir, histDir, histName, kViolet +col_depth , 1,   hStack, hMC);
  TH1F* hTT = MyESHist.stackHisto(fTT, chDir, baseDir, histDir, histName, kCyan+col_depth, 1,   hStack, hMC);

  // trim the histDir string
  std::string histDir_str;
  std::string histDir_orig(histDir);
  std::remove_copy(histDir_orig.begin(), histDir_orig.end(), std::back_inserter(histDir_str), '/');
  TString histDir_(histDir_str);
  //
  //-------------------------------
  // DY from Data
  //-------------------------------
  //qcd scale factors for data-driven DY
  double dySF = 1.0;
  double dyErr = 0.0;
  if(isDYdd){
    vector<double> sfAndErr;
    if(isMuChannel) sfAndErr = MyESHist.getTransFactDY2(fMuData, fNonDYBkg, chDir, baseDir, histDir, histName, xTitle, xmin, xmax);
    if(isEleChannel) sfAndErr = MyESHist.getTransFactDY2(fEleData, fNonDYBkg, chDir, baseDir, histDir, histName, xTitle, xmin, xmax);
    //vector<double> sfAndErr = MyESHist.getTransFactDY(fData, fTT, fWJ, fVV, chDir, baseDir, histDir, histName, xTitle, xmin, xmax);
    dySF = sfAndErr[0];
    dyErr = sfAndErr[1];
  }
  TH1F * hDYdd = MyESHist.getHisto(fDY, chDir, baseDir, histDir, histName);
  hDYdd->Reset(); // initialize empty hist
  if(isDYdd){
    hDYdd = MyESHist.getDataDrivenDY(chDir, baseDir, histDir, histName,  dySF,  dyErr);
    hDYdd->SetFillColor(kGreen+col_depth);
    hDYdd->GetXaxis()->SetRangeUser(xmin,xmax);
    //create same dir to the data driven qcd file
    std::string histPath = std::string(chDir+baseDir+histDir_);
    TDirectory *d = fDYdd->GetDirectory(histPath.c_str());
    if(!d) fDYdd->mkdir(histPath.c_str());
    fDYdd->cd(histPath.c_str());
    //hDY->Draw();
    hDYdd->Write();
    hStack->Add(hDYdd);
    hMC->Add(hDYdd);
  }
  else hDYdd = MyESHist.stackHisto(fDY, chDir, baseDir, histDir, histName, kGreen +col_depth, 1,   hStack, hMC); 

  if(isData) c1->cd(1);
  else c1->cd();
   gPad->SetTopMargin(0.10);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.05);
  hStack->Draw("HIST");
  hStack->SetMinimum(1.0);
  if(histDir.Contains("ZTag")) hStack->SetMinimum(0.1);
  hStack ->GetXaxis()->SetRangeUser(xmin, xmax);
  //cout<<hStack->GetMaximum()<<endl;
  if(isData){
    hStack->GetYaxis()->SetTitleOffset(0.70);
    hStack->GetYaxis()->SetTitleSize(0.10);   
    hStack->GetYaxis()->SetLabelSize(0.07);   
    hStack->GetYaxis()->SetTickLength(0.04); 
    hStack->GetYaxis()->SetTitle("Events");
    hStack->GetXaxis()->SetTitleOffset(1.20);
  }
  else{
  hStack->GetYaxis()->SetTitle("Events");
  hStack->GetXaxis()->SetTitle(xTitle);
  hStack->GetXaxis()->SetTitleSize(0.07);
  hStack->GetXaxis()->SetLabelSize(0.05);   
  hStack->GetXaxis()->SetTickLength(0.05); 
  //hStack->GetYaxis()->SetNdivisions(5);
  hStack->GetYaxis()->SetTitleSize(0.08);   
  hStack->GetYaxis()->SetLabelSize(0.05);   
  hStack->GetYaxis()->SetTickLength(0.04); 
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.15);
  hStack->GetYaxis()->SetTitleOffset(0.80);
  hStack->GetXaxis()->SetTitleOffset(0.90);
  }

  //-------------------------------------///
  //unc band
  //-------------------------------------///
  TGraphAsymmErrors *UncBand;
  if(unc){
  UncBand = MyESHist.UNCGRAPH(
            MyESHist.addHistoForUnc("base/", 	 baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JESPlus/",      baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JESMinus/",     baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JERPlus/",      baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JERMinus/",     baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("bTagPlus/",     baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("bTagMinus/",    baseDir, histDir, histName),
	    true, false);
  UncBand->SetFillColor(17);
  UncBand->SetFillStyle(3008);
  UncBand->Draw(" E2 same");
  }

  //-------------------------------
  //Data
  //-------------------------------
  TH1F* hData;
  if(isMuChannel) hData = MyESHist.getHisto(fMuData, chDir, baseDir, histDir, histName);
  if(isEleChannel) hData = MyESHist.getHisto(fEleData, chDir, baseDir, histDir, histName);
  ///MyESHist.decorateHisto(hData, "", xTitle, "Events");
  hData->SetFillColor(kBlack);
  hData->SetMarkerStyle(20); hData->SetMarkerSize(1.2);
  if(isData)hData->Draw("Esame"); 

  //-------------------------------
  //Signal 
  //-------------------------------
  TH1F* hSig;
  if(isMuChannel) hSig = MyESHist.getHisto(fSigMuMuZ, chDir, baseDir, histDir, histName);
  if(isEleChannel) hSig = MyESHist.getHisto(fSigEEZ, chDir, baseDir, histDir, histName);
  ///MyESHist.decorateHisto(hSig, "", xTitle, "Events");
  hSig->SetLineColor(kRed); hSig->SetLineStyle(2);
  hSig->SetLineWidth(3); hSig->SetFillColor(0);
  if(isSig)hSig->Draw("HISTSAME"); 

  //-------------------------------
  //Legends
  //-------------------------------
  TLegend* leg = new TLegend(0.7218792,0.3061504,0.9212081,0.8798861,NULL,"brNDC");
  //TLegend* leg = new TLegend(0.7618792,0.3061504,0.9712081,0.8798861,NULL,"brNDC");
  if(hist_name.find("pt") != string::npos || hist_name.find("mt") != string::npos || hist_name.find("Fit") != string::npos ||hist_name.find("RelIso") != string::npos){
    leg = new TLegend(0.6018792,0.6061504,0.9512081,0.8898861,NULL,"brNDC");
    leg->SetNColumns(2);
  }
  if(hist_name.find("mlZ") != string::npos){
    leg = new TLegend(0.6518792,0.3061504,0.8512081,0.8798861,NULL,"brNDC");
  }
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(kBlack);
  leg->SetTextFont(42);
  leg->SetTextAngle(0);
  leg->SetTextSize(0.04);
  leg->SetTextAlign(12);
  //leg->AddEntry(hSig, "#splitline{Signal}{M_{H^{+}} = 120 GeV}","L");
  if(isData)leg->AddEntry(hData,"Data","PE");
  leg->AddEntry(hTT,"t#bar{t} + jets","F");
  leg->AddEntry(hDYdd,"DY","F");
  leg->AddEntry(hVV,"VV","F");
  leg->AddEntry(hWJ,"W + jets","F");
  //if(isSig)leg->AddEntry(hSig, "Signal","L");
  if(unc)leg->AddEntry(UncBand, "Uncertainty","F");
  if(isSig)leg->AddEntry(hSig, "Signal","L");
  if(isSig)leg->AddEntry((TObject*)0, "(M_{l^{*}} = 250 GeV)","");
  leg->Draw();

  double yMax = 0;
  if(hData->GetMaximum() > hSig->GetMaximum()) yMax = hData->GetMaximum();
  else yMax = hSig->GetMaximum();
  if(yMax < hMC->GetMaximum()) yMax = hMC->GetMaximum();

  if(isData) hStack->SetMaximum(4.0*hStack->GetMaximum());
  else hStack->SetMaximum(1.1*yMax);
  TPaveText *cct = MyESHist.paveText(0.30,0.8454,0.40,0.8462, 0, 19, 1, 0, 132);
  cct->SetTextSize(0.055);
  cct->AddText("CMS, Preliminary");

  //-------------------------------------///
  //  Draw Pave Text 
  //-------------------------------------///
  //hist name
  TPaveText *hLable = MyESHist.paveText(0.6513423,0.7754898,0.6010067,0.8962187, 0, 19, 1, 0, 132);
  hLable->SetTextSize(0.07);
  hLable->AddText(xTitle);
  
  //channel
  TPaveText *ch = MyESHist.paveText(0.823,0.9154898,0.9210067,0.9762187, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.10);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  //CMS prili
  TPaveText *pt = MyESHist.paveText(0.01,0.9554,0.82,0.9562, 0, 19, 1, 0, 132);
  if(isData) pt->SetTextSize(0.080);
  else pt->SetTextSize(0.05);
  pt->AddText(histDir_+": 35.9 fb^{-1} (13 TeV)");
  //TText *text = pt->AddText(histDir+": CMS Preliminary, #sqrt{s} = 13 TeV, 35.9 fb^{-1}");
  pt->Draw();
  if(isSig) cct->Draw();
  ch->Draw();
  //hLable->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  if(isData){
    c1->cd(2);
    gPad->SetTopMargin(0); gPad->SetBottomMargin(0.5); //gPad->SetGridy();
    if(histDir=="") gPad->SetBottomMargin(0.55);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)hData->Clone("hRatio");
    hRatio->Reset();
    hRatio->Add(hData);
    hRatio->Divide(hMC); 
    MyESHist.decorateHisto(hRatio, "", xTitle, "#frac{Data}{Bkg}");
    hRatio->SetFillColor(kBlack);
    if(histDir.Contains("ZTag")) hRatio->GetYaxis()->SetRangeUser(0, 2);
    else hRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
    //else hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetYaxis()->SetTitleOffset(0.40);
    hRatio->GetXaxis()->SetTitleOffset(0.90);
    hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(1.2);
    hRatio->GetYaxis()->SetTitleSize(0.15); 
    hRatio->GetXaxis()->SetTitleSize(0.15);
    hRatio->GetXaxis()->SetLabelSize(0.10); 
    hRatio->GetYaxis()->SetLabelSize(0.10); 
    if(hist_name.find("mjj") != string::npos){
      hRatio->GetXaxis()->SetTitleSize(0.05); 
      hRatio->GetXaxis()->SetTitleOffset(1.40);
    }
    //lable x-axis, for cutflow
    if(histName=="cutflow"){
      vector<string >cut_label;
      if(isEleChannel){
        cut_label.push_back("Ele trigger");
      }
      if(isMuChannel){
        cut_label.push_back("Mu trigger");
      }
      cut_label.push_back("Control Sel");
      cut_label.push_back("b-jet veto");
      cut_label.push_back("PreSel");
      cut_label.push_back("ZTag Sel");
      for(int istep=0; istep<cut_label.size(); istep++ ){
       hRatio->GetXaxis()->SetBinLabel(istep+1, cut_label[istep].c_str());
      }
      hRatio->GetXaxis()->LabelsOption("v");
      hRatio->GetXaxis()->SetTickLength(0.08); 
      hRatio->GetXaxis()->SetLabelOffset(0.03);
      hRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    }
    //unc band
    hRatio->Draw("E"); // use "P" or "AP"
    if(unc){
    TGraphAsymmErrors *UncBand_Ratio;
    UncBand_Ratio = MyESHist.UNCGRAPH(
	    MyESHist.addHistoForUnc("base/", 	baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JESPlus/",     baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JESMinus/",    baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JERPlus/",     baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("JERMinus/",    baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("bTagPlus/",    baseDir, histDir, histName),
      	    MyESHist.addHistoForUnc("bTagMinus/",   baseDir, histDir, histName),
	    false, true);
    UncBand_Ratio->SetFillColor(17);
    //UncBand_Ratio->SetFillStyle(19);
    UncBand_Ratio->Draw("E2 same");
    }
    hRatio->Draw("E same"); // use "P" or "AP"
    //base line at 1
    TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
    baseLine->SetLineColor(kBlack);
    baseLine->Draw("SAME");
    c1->Update();
  }
  if(isSaveHisto){
    mkdir(chDir, S_IRWXU);
    mkdir(chDir+histDir, S_IRWXU);
    //mkdir(histDir_, S_IRWXU);
    TString outFile("$PWD/");
    outFile += chDir+histDir+histName;
    if(isMuChannel) outFile += "_mu"+histDir_+".pdf";
    if(isEleChannel) outFile += "_ele"+histDir_+".pdf";
    c1->SaveAs(outFile);
    //c1->Close();
  }
}
void stackAllHisto(TString chDir, TString histDir, bool isDYdd){
  TString baseDir = "base/";
  bool isDataMjj= false;
  bool isData = true;
  //flags
  bool isSig = false;
  bool isUnc = false;
  bool isDDtmp = false;
  //---------------------------------------
  bool isMuChannel = false;
  bool isEleChannel = false;
  if(chDir.Contains("Muon")) isMuChannel = true;
  if(chDir.Contains("Electron")) isEleChannel = true;
  //---------------------------------------
  //plotStackedHisto(chDir, baseDir, "", "cutflow","cutflow", isData,  isSig,  0.5, 7.5, false);
  if(isMuChannel){
    plotStackedHisto(chDir, baseDir, histDir, "mll", "M^{#mu_{1}#mu_{2}} [GeV]", isData, isSig,  50, 1000, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "pt_Z", "Pt^{#mu_{1}#mu_{2}}[GeV]", isData, isSig,  0, 1000, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "eta_Z", "#eta^{#mu_{1}#mu_{2}}", isData, isSig,  -3.0, 4.5, isUnc);
  }
  if(isEleChannel){
    plotStackedHisto(chDir, baseDir, histDir, "mll", "M^{e_{1}e_{2}} [GeV]", isData, isSig,  50, 1000, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "pt_Z", "Pt^{e_{1}e_{2}}[GeV]", isData, isSig,  0, 1000, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "eta_Z", "#eta^{e_{1}e_{2}}", isData, isSig,  -3.0, 4.5, isUnc);
  }
}

void MyEtherStackHisto(){
  bool isDYdd = false;
  stackAllHisto("Muon/", "ControlP2/",       isDYdd);
  //stackAllHisto("Electron/", "ControlP/",   isDYdd);
}

