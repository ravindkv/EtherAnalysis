
//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;

//INPUT FILES
TFile* fData = TFile::Open("all_Data.root");
TFile* fVV	= TFile::Open("all_VV.root");
TFile* fDY	= TFile::Open("all_DY_M50.root");
//TFile* fDY	= TFile::Open("all_DY_Pt.root");
TFile* fTT	= TFile::Open("all_TT.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
//TFile *fSig     = TFile::Open("all_ExLepEEZ_M250.root");
TFile *fSig250     = TFile::Open("all_ExLep_M250.root");
TFile *fSig500     = TFile::Open("all_ExLep_M500.root");
TFile *fSig750     = TFile::Open("all_ExLep_M750.root");
TFile *fSig1000     = TFile::Open("all_ExLep_M1000.root");
TFile *fSig1250     = TFile::Open("all_ExLep_M1250.root");
TFile *fSig1500     = TFile::Open("all_ExLep_M1500.root");
TFile *fSig1750     = TFile::Open("all_ExLep_M1750.root");
TFile *fSig2000     = TFile::Open("all_ExLep_M2000.root");
TFile *fSig2500     = TFile::Open("all_ExLep_M2500.root");
TFile *fSig3000     = TFile::Open("all_ExLep_M3000.root");
TFile *fSig3500     = TFile::Open("all_ExLep_M3500.root");
TFile *fSig4000     = TFile::Open("all_ExLep_M4000.root");
TFile *fSig4500     = TFile::Open("all_ExLep_M4500.root");
TFile *fSig5000     = TFile::Open("all_ExLep_M5000.root");

TH1F* getHisto(TFile *histFile, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  TString fullPath = histPath+"/"+dir+"/"+histName;
  if(!(histFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(histFile->Get(fullPath))->Clone(histName);
  return hist;
}

TH1F* getHistoSig(TString fileName, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  TFile * histFile = TFile::Open(fileName);
  TString fullPath = histPath+"/"+dir+"/"+histName;
  if(!(histFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(histFile->Get(fullPath))->Clone(histName);
  return hist;
}

double getBinContent(TH1F *hist, double tau21){
  double passedEvt = 0.0;
  TAxis *xaxis = hist->GetXaxis();
  int binX = 1;
  binX = xaxis->FindBin(tau21);
  passedEvt = hist->GetBinContent(binX);
  return passedEvt; 
}
double getSigEff(TString fileName, double tau21){
  TH1F * histTotal = (TH1F*) getHistoSig(fileName, "base/Iso", "ZTag", "totalEvt");
  TH1F * histPassed = (TH1F*) getHistoSig(fileName, "base/Iso", "ZTag", "passedEvt");
  double totalEvt = histTotal->Integral();
  double passedEvt = getBinContent(histPassed, tau21);
  double effSig = 0.0;
  effSig = 100*passedEvt/totalEvt;
  return effSig;
}


double getBkgEff(double tau21){
  TH1F * hTotalVV = (TH1F*) getHisto(fVV, "base/Iso", "ZTag", "totalEvt");
  TH1F * hTotalTT = (TH1F*) getHisto(fTT, "base/Iso", "ZTag", "totalEvt");
  TH1F * hTotalWJ = (TH1F*) getHisto(fWJ, "base/Iso", "ZTag", "totalEvt");
  TH1F * hTotalDY = (TH1F*) getHisto(fDY, "base/Iso", "ZTag", "totalEvt");
  TH1F * histTotal = (TH1F*) hTotalVV->Clone();
  histTotal->Add(hTotalTT);
  histTotal->Add(hTotalWJ);
  histTotal->Add(hTotalDY);

  TH1F * hPassedVV = (TH1F*) getHisto(fVV, "base/Iso", "ZTag", "passedEvt");
  TH1F * hPassedTT = (TH1F*) getHisto(fTT, "base/Iso", "ZTag", "passedEvt");
  TH1F * hPassedWJ = (TH1F*) getHisto(fWJ, "base/Iso", "ZTag", "passedEvt");
  TH1F * hPassedDY = (TH1F*) getHisto(fDY, "base/Iso", "ZTag", "passedEvt");
  TH1F * histPassed = (TH1F*) hPassedVV->Clone();
  histPassed->Add(hPassedTT);
  histPassed->Add(hPassedWJ);
  histPassed->Add(hPassedDY);

  double totalEvt = histTotal->Integral();
  double passedEvt = getBinContent(histPassed, tau21);
  double effBkg = 0.0;
  effBkg = 100*passedEvt/totalEvt;
  return effBkg;
  //return passedEvt;
}

double getRatioEff(TString fileName, double tau21){
  double eSig = getSigEff(fileName, tau21);
  double eBkg = getBkgEff(tau21);
  //double ratio = eSig/sqrt(eBkg);
  double ratio = eSig/eBkg;
  //cout<<setprecision(4)<<tau21<<"\t"<<eSig<<"\t"<<eBkg<<"\t"<<ratio<<endl;
  cout<<setprecision(4)<<tau21<<"\t"<<eSig<<"\t"<<sqrt(eBkg)<<"\t"<<ratio<<endl;
  //return eBkg;
  return ratio;
}

//--------------------------------------------
//function to make graph from two array
//--------------------------------------------
TGraph* makeGraph(vector<double> Xarray, vector<double> Yarray){
  const int n = Xarray.size();
  double x[n], y[n];
   for (Int_t i=0;i<n;i++) {
     x[i] = Xarray[i];
     y[i] = Yarray[i] ;
   }
   TGraph* gr = new TGraph(n,x,y);
   return gr;
}

//--------------------------------------------
//function to decorate Graph  
//--------------------------------------------
TGraph* decorateGraph(TGraph *graph, TString xTitle, TString yTitle, TString myTitle, double yMin, double yMax, int color){
  graph->SetTitle(myTitle);
  graph->GetYaxis()->SetTitleOffset(1.30);
  graph->GetXaxis()->SetTitle(xTitle);
  graph->GetXaxis()->SetTitleSize(0.04);
  graph->GetXaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetTitle(yTitle);
  graph->GetYaxis()->SetTitleSize(0.04);
  graph->GetYaxis()->SetRangeUser(yMin, yMax);
  graph->GetXaxis()->SetRangeUser(0.44, 0.72);
  graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetXaxis()->SetNdivisions(14);
  graph->GetYaxis()->SetNdivisions(10);
  graph->GetYaxis()->CenterTitle();
  graph->SetLineColor(color);
  graph->SetLineWidth(3);
  graph->SetMarkerSize(20);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(color);
  return graph;
}

//-----------------------------------------------
//overlay all the graphs in one canvas
//-----------------------------------------------
void overlayGraphs(std::vector<std::string> diode_pos, string colX, string colY, TString xTit, TString yTit, TString myTit, int yMin, int yMax, TCanvas *c1){
}

void MyOptimizeTau21(){
  vector<string>massPoints;
  massPoints.push_back("all_ExLep_M250.root");
  massPoints.push_back("all_ExLep_M500.root");
  massPoints.push_back("all_ExLep_M750.root");
  massPoints.push_back("all_ExLep_M1000.root");
  massPoints.push_back("all_ExLep_M1250.root");
  massPoints.push_back("all_ExLep_M1500.root");
  massPoints.push_back("all_ExLep_M1750.root");
  massPoints.push_back("all_ExLep_M2000.root");
  massPoints.push_back("all_ExLep_M2500.root");
  massPoints.push_back("all_ExLep_M3000.root");
  massPoints.push_back("all_ExLep_M3500.root");
  massPoints.push_back("all_ExLep_M4000.root");
  massPoints.push_back("all_ExLep_M4500.root");
  massPoints.push_back("all_ExLep_M5000.root");
  

  vector<string>massLables;
  massLables.push_back("M_{l*} = 250 GeV");
  massLables.push_back("M_{l*} = 500 GeV");
  massLables.push_back("M_{l*} = 750 GeV");
  massLables.push_back("M_{l*} = 1000 GeV");
  massLables.push_back("M_{l*} = 1250 GeV");
  massLables.push_back("M_{l*} = 1500 GeV");
  massLables.push_back("M_{l*} = 1750 GeV");
  massLables.push_back("M_{l*} = 2000 GeV");
  massLables.push_back("M_{l*} = 2500 GeV");
  massLables.push_back("M_{l*} = 3000 GeV");
  massLables.push_back("M_{l*} = 3500 GeV");
  massLables.push_back("M_{l*} = 4000 GeV");
  massLables.push_back("M_{l*} = 4500 GeV");
  massLables.push_back("M_{l*} = 5000 GeV");
  TCanvas * can = new TCanvas();
  TLegend* leg = new TLegend(0.50,0.15,0.85,0.45,NULL,"brNDC");
  gPad->SetGridy();
  gPad->SetGridx();
  leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  for(int i = 0; i<massPoints.size(); i++){
    cout<<endl;
    cout<<"==============="<<massPoints[i]<<"============="<<endl; 
    vector<double> tauVec;
    vector<double> ratioEffVec;
    //cout<<"Tau21"<<"\t"<<"effSig"<<"\t"<<"effBkg"<<"\t"<<"effSig/effBkg"<<endl;
    cout<<"Tau21"<<"\t"<<"effSig"<<"\t"<<"sqrt(nBkg)"<<"\t"<<"effSig/sqrt(nBkg)"<<endl;
    cout<<""<<"\t"<<"(%)"<<"\t"<<""<<"\t"<<"(%)"<<endl;
    //cout<<""<<"\t"<<"(%)"<<"\t"<<"(%)"<<"\t"<<""<<endl;
    for(int t = 0; t<26; t++){
      double tau = 0.45+t/100.0;
      double ratioEff = getRatioEff(TString(massPoints[i]), tau);
      tauVec.push_back(tau);
      ratioEffVec.push_back(ratioEff);
    }
   can->cd();
   TGraph * myGraph = makeGraph(tauVec, ratioEffVec);
   int color = i+1;
   if(i==9) color = 433;
   //decorateGraph(myGraph, "Tau21", "N_{Bkg}", "optimization of tau21: #mu -channel", 200, 450, color);
   //decorateGraph(myGraph, "Tau21", "100* effSig", "optimization of tau21: #mu -channel", 45, 100, color);
   //decorateGraph(myGraph, "Tau21", "100* effSig/sqrt(nBkg)", "optimization of tau21: #mu -channel", 3.0, 5.1, color);
   if(isMuChannel)decorateGraph(myGraph, "#tau_{21}", "#epsilon_{Signal}/#epsilon_{Bkg}", "optimization of #tau_{21} :  #mu-channel", 0.75, 1.15, color);
   if(isEleChannel)decorateGraph(myGraph, "#tau_{21}", "#epsilon_{Signal}/#epsilon_{Bkg}", "optimization of #tau_{21} :  e-channel", 0.75, 1.15, color);
   if(i==0) myGraph->Draw("ALP");
   else myGraph->Draw("LPsame");
   leg->AddEntry(myGraph, TString(massLables[i]),"PL");
  }
  leg->Draw();
  if(isMuChannel)can->SaveAs("tau21_mu.pdf");
  if(isEleChannel)can->SaveAs("tau21_ele.pdf");
}


