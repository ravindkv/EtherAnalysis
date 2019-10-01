
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
TFile * fBkg = TFile::Open("all_Bkg.root");
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
  graph->GetXaxis()->SetLabelSize(0.03);
  graph->GetYaxis()->SetTitle(yTitle);
  graph->GetYaxis()->SetTitleSize(0.04);
  graph->GetYaxis()->SetRangeUser(yMin, yMax);
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
void printOutput(TFile *fBkg, TString sigFile,  double mass){
  TH1F * hOneTrigBkg = (TH1F*)fBkg->Get("base/oneTrig");
  TH1F * hTowTrigBkg = (TH1F*)fBkg->Get("base/twoTrig");
  double oneTrigBkg = hOneTrigBkg->Integral();
  double twoTrigBkg = hTowTrigBkg->Integral();
  
  TFile *fSig = TFile::Open(sigFile);
  TH1F * hOneTrigSig = (TH1F*)fSig->Get("base/oneTrig");
  TH1F * hTowTrigSig = (TH1F*)fSig->Get("base/twoTrig");
  double oneTrigSig = hOneTrigSig->Integral();
  double twoTrigSig = hTowTrigSig->Integral();
  cout<<"--------------------------------------------------------------"<<endl;
  cout<<std::setprecision(4)<<endl;
  cout<<setw(20)<<"Trigger Selection"<<setw(15)<<"Siganl,M"<<mass<<setw(15)<<"AllBkg"<<setw(15)<<"Sig/AllBkg"<<endl;
  cout<<setw(20)<<"One Trigger Sel"<<setw(15)<<oneTrigSig<<setw(15)<<oneTrigBkg<<setw(15)<<oneTrigSig/oneTrigBkg<<endl;
  cout<<setw(20)<<"Tow Trigger  Sel"<<setw(15)<<twoTrigSig<<setw(15)<<twoTrigBkg<<setw(15)<<twoTrigSig/twoTrigBkg<<endl;
}

double getRatio(TFile *fBkg, TString sigFile, TString histName){
  TFile *fSig = TFile::Open(sigFile);
  TH1F * hBkg = (TH1F*)fBkg->Get("base/"+histName);
  TH1F * hSig = (TH1F*)fSig->Get("base/"+histName);
  double bkg = hBkg->Integral();
  double sig = hSig->Integral();
  double ratio = 1.0;
  ratio = sig/bkg;
  return ratio;
}
void MyOptimizeMuTrig(){
  vector<string>massFiles;
  massFiles.push_back("all_ExLep_M250.root");
  massFiles.push_back("all_ExLep_M500.root");
  massFiles.push_back("all_ExLep_M750.root");
  massFiles.push_back("all_ExLep_M1000.root");
  massFiles.push_back("all_ExLep_M1250.root");
  massFiles.push_back("all_ExLep_M1500.root");
  massFiles.push_back("all_ExLep_M1750.root");
  massFiles.push_back("all_ExLep_M2000.root");
  massFiles.push_back("all_ExLep_M2500.root");
  massFiles.push_back("all_ExLep_M3000.root");
  massFiles.push_back("all_ExLep_M3500.root");
  massFiles.push_back("all_ExLep_M4000.root");
  massFiles.push_back("all_ExLep_M4500.root");
  massFiles.push_back("all_ExLep_M5000.root");

  vector<double>massPoints;
  massPoints.push_back(250);
  massPoints.push_back(500);
  massPoints.push_back(750);
  massPoints.push_back(1000);
  massPoints.push_back(1250);
  massPoints.push_back(1500);
  massPoints.push_back(1750);
  massPoints.push_back(2000);
  massPoints.push_back(2500);
  massPoints.push_back(3000);
  massPoints.push_back(3500);
  massPoints.push_back(4000);
  massPoints.push_back(4500);
  massPoints.push_back(5000);
  
  //get info regarding trig selection
  vector<double> ratioOneTrig;
  vector<double> ratioTwoTrig;
  for(int i = 0; i<massPoints.size(); i++){
    double ratioOne = 1000*getRatio(fBkg, TString(massFiles[i]), "oneTrig");
    double ratioTwo = 1000*getRatio(fBkg, TString(massFiles[i]), "twoTrig");
    cout<<ratioOne<<"\t"<<ratioTwo<<endl;
    ratioOneTrig.push_back(ratioOne);
    ratioTwoTrig.push_back(ratioTwo);
    printOutput(fBkg, TString(massFiles[i]), massPoints[i]);
  }
  //create graph for trig selection
  TGraph * graphOneTrig = makeGraph(massPoints, ratioOneTrig);
  TGraph * graphTwoTrig = makeGraph(massPoints, ratioTwoTrig);
  //decorate graph
  TString chName = "muon";
  double yMin = 0.0;
  double yMax = 0.0;
  if(isMuChannel){
    chName = "#mu -channel";
    yMin = 3.5;
    yMax = 5.0;
  }
  else{
    chName = "e -channel";
    yMin = 0.1;
    yMax = 2.0;
  }
  decorateGraph(graphOneTrig, "M_{l^{*}}", "1000*Sig/Bkg", chName, yMin, yMax, kGreen);
  decorateGraph(graphTwoTrig, "M_{l^{*}}", "1000*Sig/Bkg", chName, yMin, yMax, kRed);
  //draw graphs in a canvas
  TCanvas * canCharge = new TCanvas();
  canCharge->cd();
  gPad->SetGridy();
  gPad->SetGridx();
  graphOneTrig->Draw();
  graphTwoTrig->Draw("same");
  //draw legend
  TLegend* leg = new TLegend(0.45,0.15,0.70,0.30,NULL,"brNDC");
  leg->AddEntry(graphOneTrig, "HLT_IsoMu24","PL");
  leg->AddEntry(graphTwoTrig, "HLT_IsoMu24 || HLT_IsoTkMu","PL");
  //leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->Draw();
  if(isMuChannel)canCharge->SaveAs("trigSel_mu.png");
  if(isEleChannel)canCharge->SaveAs("trigSel_ele.png");
}

