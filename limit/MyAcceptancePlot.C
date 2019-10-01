
//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;
TFile* fTT	= TFile::Open("all_TT.root");

TH1F* getHistoSig(TString fileName, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  TFile * histFile = TFile::Open(fileName);
  TString fullPath = histPath+"/"+dir+"/"+histName;
  //cout<<fileName<<"/"<<fullPath<<endl;
  if(!(histFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(histFile->Get(fullPath))->Clone(histName);
  return hist;
}

double getSigEff(TString fileName, TString dirIn, TString dirOut, TString histNameIn, TString histNameOut, TString mass){
  TH1F * histIn = (TH1F*) getHistoSig(fileName, "base/Iso", dirIn, histNameIn);
  TH1F * histOut = (TH1F*) getHistoSig(fileName, "base/Iso", dirOut, histNameOut);
  double evtIn = histIn->Integral();
  double evtOut = histOut->Integral();
  double effSig = 0.0;
  effSig = evtOut/evtIn;
  printf("eventIn =%.5f\teventOut = %.5f\teff= %.5f\t ", evtIn, evtOut, effSig);
  return effSig;
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
  graph->GetXaxis()->SetLabelSize(0.04);
  graph->GetYaxis()->SetTitle(yTitle);
  graph->GetYaxis()->SetTitleSize(0.04);
  graph->GetYaxis()->SetRangeUser(yMin, yMax);
  //graph->GetXaxis()->SetRangeUser(0.44, 0.72);
  graph->GetYaxis()->SetLabelSize(0.05);
  graph->GetXaxis()->SetNdivisions(14);
  graph->GetYaxis()->SetNdivisions(5);
  graph->GetYaxis()->CenterTitle();
  graph->SetLineColor(color);
  graph->SetLineWidth(4);
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

void MyAcceptancePlot(){
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
  
  vector<double> sigMass;         vector<string> massLabel; 
  sigMass.push_back(250);         massLabel.push_back("250"); 
  sigMass.push_back(500);         massLabel.push_back("500"); 
  sigMass.push_back(750);         massLabel.push_back("750"); 
  sigMass.push_back(1000);        massLabel.push_back("1000");
  sigMass.push_back(1250);        massLabel.push_back("1250");
  sigMass.push_back(1500);        massLabel.push_back("1500");
  sigMass.push_back(1750);        massLabel.push_back("1750");
  sigMass.push_back(2000);        massLabel.push_back("2000");
  sigMass.push_back(2500);        massLabel.push_back("2500");
  sigMass.push_back(3000);        massLabel.push_back("3000");
  sigMass.push_back(3500);        massLabel.push_back("3500");
  sigMass.push_back(4000);        massLabel.push_back("4000");
  sigMass.push_back(4500);        massLabel.push_back("4500");
  sigMass.push_back(5000);        massLabel.push_back("5000");
  vector<double>sigXss;
  if(isMuChannel){
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
  if(isEleChannel){
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

  double lumi = 35900;
  TCanvas * can = new TCanvas();
  gPad->SetGridy();
  gPad->SetGridx();
  vector<double> massVec;
  vector<double> effVec;
  vector<double> accVec;
  vector<double> accEffVec;
  TString dirIn = "PreSel";
  TString dirOut = "ZTag";
  TString histNameIn = "mll";
  TString histNameOut = "mll";

  for(int i = 0; i<massFiles.size(); i++){
    cout<<endl;
    cout<<"-------------------"<<endl;
    cout<<massFiles[i]<<endl; 
    TH1F * histOut = (TH1F*) getHistoSig(massFiles[i], "base/Iso", dirOut, histNameOut);
    double evtOut = histOut->Integral();
    double acc = evtOut/(lumi*sigXss[i]);
    double eff = getSigEff(TString(massFiles[i]), dirIn, dirOut, histNameIn, histNameOut, TString(massLabel[i]));
    massVec.push_back(sigMass[i]);
    effVec.push_back(eff);
    accVec.push_back(acc);
    accEffVec.push_back(acc*eff);
    cout<<endl;
    cout<<"eff = "<<eff<<endl;
    cout<<"acc = "<<acc<<endl;
    cout<<"acc*eff = "<<acc*eff<<endl;
  }
  TString xAxis = "M_{lep^{*}} [GeV]";
  if(isMuChannel) xAxis = "M_{#mu^{*}} [GeV]";
  if(isEleChannel)xAxis = "M_{e^{*}} [GeV]";
  
  can->cd();
  TGraph * graphEff = makeGraph(massVec, effVec);
  decorateGraph(graphEff, xAxis, "For PreSel and ZTag", "", 0, 1, kRed+1);
  graphEff->Draw("");
  TGraph * graphAcc = makeGraph(massVec, accVec);
  decorateGraph(graphAcc, xAxis, "For PreSel and ZTag", "", 0, 1, kGreen+1);
  graphAcc->Draw("same");
  TGraph * graphAccEff = makeGraph(massVec, accEffVec);
  decorateGraph(graphAccEff, xAxis, "For PreSel and ZTag", "", 0, 1, kBlue+1);
  graphAccEff->Draw("same");

  //pave text channel box
  TPaveText *ch = new TPaveText(0.75,0.9154898,0.80,0.9762187,"brNDC");
  ch->SetFillColor(19);
  ch->SetFillStyle(0);
  ch->SetLineColor(0);
  ch->SetTextSize(0.08);
  ch->SetBorderSize(1);
  if(isMuChannel)ch->AddText("#mu-channel");
  if(isEleChannel) ch->AddText("e-channel");
  ch->Draw("SAME");
  //pave text CMS box
  TPaveText *pt = new TPaveText(0.20,0.9354,0.50,0.9362, "brNDC"); // good_v1
  pt->SetBorderSize(1);
  pt->SetFillColor(19);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.06);
  pt->SetLineColor(0);
  pt->SetTextFont(132);
  TText *text = pt->AddText("13 TeV (35.9 fb^{-1}) ");
  text->SetTextAlign(11);
  pt->Draw("SAME");
  if(isMuChannel)can->SaveAs("accept_mu.png");
  if(isEleChannel)can->SaveAs("accept_ele.pdf");
}


