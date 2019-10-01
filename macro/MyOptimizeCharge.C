#include <cmath>
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
  //graph->GetYaxis()->SetRangeUser(yMin, yMax);
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
double getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);  
  //double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;
  double statUnc = fabs(sError);
  return statUnc;
}

void printOutput(TFile *fBkg, TString sigFile,  string mass){
  TH1F * hNoChBkg = (TH1F*)fBkg->Get("base/noCharge");
  TH1F * hOppChBkg = (TH1F*)fBkg->Get("base/oppCharge");
  TH1F * hSameChBkg = (TH1F*)fBkg->Get("base/sameCharge");
  double noChBkg = hNoChBkg->Integral();
  double oppChBkg = hOppChBkg->Integral();
  double sameChBkg = hSameChBkg->Integral();
  double noChBkgErr   = 0.0;
  double oppChBkgErr  = 0.0;
  double sameChBkgErr = 0.0;
  noChBkgErr   = getStatUnc(hNoChBkg, 0.0);
  oppChBkgErr  = getStatUnc(hOppChBkg, 0.0);
  sameChBkgErr = getStatUnc(hSameChBkg, 0.0);
    
  TFile *fSig = TFile::Open(sigFile);
  TH1F * hNoChSig = (TH1F*)fSig->Get("base/noCharge");
  TH1F * hOppChSig = (TH1F*)fSig->Get("base/oppCharge");
  TH1F * hSameChSig = (TH1F*)fSig->Get("base/sameCharge");
  double noChSig = hNoChSig->Integral();
  double oppChSig = hOppChSig->Integral();
  double sameChSig = hSameChSig->Integral();
  double noChSigErr   = 0.0;
  double oppChSigErr  = 0.0;
  double sameChSigErr = 0.0;
  noChSigErr   = getStatUnc(hNoChSig, 0.0);
  oppChSigErr  = getStatUnc(hOppChSig, 0.0);
  sameChSigErr = getStatUnc(hSameChSig, 0.0);

  TH1F* hRatioNoCh = (TH1F*)hNoChSig->Clone("hRatioNoCh");  
  hRatioNoCh->Divide(hNoChBkg);
  TH1F* hRatioOppCh = (TH1F*)hOppChSig->Clone("hRatioOppCh");  
  hRatioOppCh->Divide(hOppChBkg);
  TH1F* hRatioSameCh = (TH1F*)hSameChSig->Clone("hRatioSameCh");  
  hRatioSameCh->Divide(hSameChBkg);
  double rNoChSig   = hRatioNoCh->Integral();
  double rOppChSig  = hRatioOppCh->Integral();
  double rSameChSig = hRatioSameCh->Integral();
  double rNoChSigErr   = 0.0;
  double rOppChSigErr  = 0.0;
  double rSameChSigErr = 0.0;
  rNoChSigErr   = getStatUnc(hRatioNoCh, 0.0);
  rOppChSigErr  = getStatUnc(hRatioOppCh, 0.0);
  rSameChSigErr = getStatUnc(hRatioSameCh, 0.0);

  cout<<"--------------------------------------------------------------"<<endl;
  cout<<std::setprecision(3)<<endl;
  //cout<<setw(20)<<"Charge Selection"<<setw(15)<<"Siganl,M"<<mass<<setw(15)<<"AllBkg"<<setw(15)<<"Sig/AllBkg"<<endl;
  //cout<<setw(20)<<"NoChargeSel"<<setw(15)<<noChSig<<setw(15)<<noChBkg<<setw(15)<<noChSig/noChBkg<<endl;
  //cout<<setw(20)<<"OppChargeSel"<<setw(15)<<oppChSig<<setw(15)<<oppChBkg<<setw(15)<<oppChSig/oppChBkg<<endl;
  //cout<<setw(20)<<"SameChargeSel"<<setw(15)<<sameChSig<<setw(15)<<sameChBkg<<setw(15)<<sameChSig/sameChBkg<<endl;
cout<<setw(15)<<"Selection"<<setw(15)<<"Nsig(M"+mass+")"<<setw(10)<<"Esig"<<setw(10)<<"Nbkg"<<setw(10)<<"Ebkg"<<setw(10)<<"Nsig/Nbkg"<<setw(10)<<"E(S/B)"<<endl;
cout<<setw(15)<<"NoChargeSel"<<setw(15)<<noChSig<<setw(10)<<noChSigErr<<setw(10)<<noChBkg<<setw(10)<<noChBkgErr<<setw(10)<<rNoChSig<<setw(10)<<rNoChSigErr<<endl;
cout<<setw(15)<<"OppChargeSel"<<setw(15)<<oppChSig<<setw(10)<<oppChSigErr<<setw(10)<<oppChBkg<<setw(10)<<oppChBkgErr<<setw(10)<<rOppChSig<<setw(10)<<rOppChSigErr<<endl;
cout<<setw(15)<<"SameChargeSel"<<setw(15)<<sameChSig<<setw(10)<<sameChSigErr<<setw(10)<<sameChBkg<<setw(10)<<sameChBkgErr<<setw(10)<<rSameChSig<<setw(10)<<rSameChSigErr<<endl;
}
double getRatio(TFile *fBkg, TString sigFile, TString histName){
  TFile *fSig = TFile::Open(sigFile);
  TH1F * hBkg = (TH1F*)fBkg->Get("base/"+histName);
  TH1F * hSig = (TH1F*)fSig->Get("base/"+histName);
  double bkg = hBkg->Integral();
  double sig = hSig->Integral();
  double ratio = 1.0;
  //cout<<sig<<", "<<bkg<<endl;
  ratio = sig/bkg;
  return ratio;
}

void MyOptimizeCharge(){
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

  vector<string>massLabel;
  massLabel.push_back("250");
  massLabel.push_back("500");
  massLabel.push_back("750");
  massLabel.push_back("1000");
  massLabel.push_back("1250");
  massLabel.push_back("1500");
  massLabel.push_back("1750");
  massLabel.push_back("2000");
  massLabel.push_back("2500");
  massLabel.push_back("3000");
  massLabel.push_back("3500");
  massLabel.push_back("4000");
  massLabel.push_back("4500");
  massLabel.push_back("5000");
  
  //get info regarding charge selection
  vector<double> ratioOppCharge;
  vector<double> ratioNoCharge;
  for(int i = 0; i<massPoints.size(); i++){
    double ratioOpp = getRatio(fBkg, TString(massFiles[i]), "oppCharge");
    double ratioNo = getRatio(fBkg, TString(massFiles[i]), "noCharge");
    //cout<<ratioOpp<<"\t"<<ratioNo<<endl;
    ratioOppCharge.push_back(ratioOpp);
    ratioNoCharge.push_back(ratioNo);
    printOutput(fBkg, TString(massFiles[i]), massLabel[i]);
  }
  //create graph for charge selection
  TGraph * graphOppCharge = makeGraph(massPoints, ratioOppCharge);
  TGraph * graphNoCharge = makeGraph(massPoints, ratioNoCharge);
  //decorate graph
  TString chName = "muon";
  double yMin = 0.0;
  double yMax = 0.0;
  if(isMuChannel){
    chName = "#mu -channel";
    yMin = 1.7;
    yMax = 2.2;
  }
  else{
    chName = "e -channel";
    yMin = 0.1;
    yMax = 2.0;
  }
  decorateGraph(graphOppCharge, "M_{l^{*}}", "Sig/Bkg", chName, yMin, yMax, kGreen);
  decorateGraph(graphNoCharge, "M_{l^{*}}", "Sig/Bkg", chName, yMin, yMax, kBlue);
  //draw graphs in a canvas
  TCanvas * canCharge = new TCanvas();
  canCharge->cd();
  gPad->SetGridy();
  gPad->SetGridx();
  graphOppCharge->Draw();
  graphNoCharge->Draw("same");
  //draw legend
  TLegend* leg = new TLegend(0.50,0.15,0.70,0.30,NULL,"brNDC");
  leg->AddEntry(graphNoCharge, "no charge selection","PL");
  leg->AddEntry(graphOppCharge, "opposite charge selection","PL");
  //leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->Draw();
  if(isMuChannel)canCharge->SaveAs("chargeSel_mu.png");
  if(isEleChannel)canCharge->SaveAs("chargeSel_ele.png");
}

