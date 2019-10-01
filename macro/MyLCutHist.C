
//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;

//INPUT FILES
TFile* fData = TFile::Open("all_Data.root");
TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY.root");
//TFile* fDY	= TFile::Open("all_DY_M50.root");
//TFile* fDY	= TFile::Open("all_DY_Pt.root");
TFile* fTT	= TFile::Open("all_TT.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile *fAllBkg     = TFile::Open("all_Bkg.root");
TFile *fSig250     = TFile::Open("all_ExLep_M250.root");

string doubleToStr(double val){
     ostringstream convert;
     string result("");
     convert <<std::setprecision(4)<<val;
     result = convert.str();
  return result;
}

void MyLCutHist(TString histName = "/base/Iso/ZTag/mlZmin_mlZmax", TString xaxis_title="mlZmin_mlZmax"){
  gStyle->SetFrameLineWidth(2);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  
  TCanvas *c = new TCanvas();
  gPad->SetBottomMargin(2.);
  gPad->SetLeftMargin(0.11);
  
  TH2F *hAllBkg = (TH2F*)fAllBkg->Get(histName);
  hAllBkg->SetLineWidth(2);
  hAllBkg->SetMarkerColor(kBlack);
  hAllBkg->SetLineColor(kBlack);
  if(isMuChannel)
    hAllBkg->SetTitle("#mu#mu*->#mu#muZ->#mu#muqq");
  else
    hAllBkg->SetTitle("ee*->eeZ->eeqq");
  hAllBkg->GetXaxis()->CenterTitle();
  hAllBkg->GetYaxis()->CenterTitle();
  hAllBkg->GetXaxis()->SetTitle("M_{lZ}^{min} [GeV]");
  hAllBkg->GetYaxis()->SetTitle("M_{lZ}^{max} [GeV]");
  hAllBkg->GetXaxis()->SetRangeUser(0, 5500);
  hAllBkg->GetXaxis()->SetTitleOffset(1.2);
  hAllBkg->GetYaxis()->SetTitleOffset(1.5);
  hAllBkg->GetYaxis()->SetRangeUser(0, 8000);
  hAllBkg->Draw();
  
  TLegend *leg = new TLegend(0.15,0.75,0.45,0.89,"");
  leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->SetFillColor(0);
  //mass vector
  vector<double> sigMass;
  vector<double> lCutMin;
  vector<double> lCutMax;
  sigMass.push_back(250);   lCutMax.push_back(540);  lCutMin.push_back(300);
  sigMass.push_back(750);   lCutMax.push_back(800);  lCutMin.push_back(900);
  sigMass.push_back(1250);  lCutMax.push_back(1300); lCutMin.push_back(1370);
  sigMass.push_back(1750);  lCutMax.push_back(1400); lCutMin.push_back(1950);
  sigMass.push_back(2500);  lCutMax.push_back(1400); lCutMin.push_back(2700);
  sigMass.push_back(3500);  lCutMax.push_back(1400); lCutMin.push_back(3700);
  sigMass.push_back(4000);  lCutMax.push_back(1400); lCutMin.push_back(4200);
  sigMass.push_back(5000);  lCutMax.push_back(1400); lCutMin.push_back(5200);
  //arrow
  TArrow *arX = new TArrow(0, 0, 0 , 0 ,0.02,"|");
  arX->SetLineWidth(2);
  //for loop
  for(int i=0; i<sigMass.size(); i++){
    TString mass = TString(doubleToStr(sigMass[i]));
    TFile *fSig    = TFile::Open("all_ExLep_M"+mass+".root");
    TH2F *hSig = (TH2F*)fSig->Get(histName);
    hSig->SetMarkerColor(i+2);
    hSig->SetLineColor(i+2);
    leg->AddEntry(hSig,"M_{l^{*}}= "+mass+" GeV","lp");
    hSig->Draw("same");
    //draw arrow
    double min = lCutMax[i];
    double max = lCutMin[i];
    arX->SetLineColor(i+2);
    arX->DrawArrow(max, min, max, 8000, 0.02, "|");
    arX->Draw("same");
    arX->DrawArrow(0, min, max, min, 0.02, "|");
    arX->Draw("same");
  }
  leg->AddEntry(hAllBkg,"AllBkg","lp");
  
  TH2F *hData = (TH2F*)fData->Get(histName);
  hData->SetLineWidth(2);
  hData->SetMarkerColor(kYellow +2);
  hData->SetLineColor(kYellow +2);
  //hData->Draw("same");
  //leg->AddEntry(hData,"Data","lp");
  leg->Draw();
  if(isMuChannel) c->SaveAs("lCutHist_mu.png");
  else c->SaveAs("lCutHist_ele.png");
}
