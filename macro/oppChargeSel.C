
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
TFile *fSig250     = TFile::Open("all_ExLepMuMuZ_M250.root");
TFile *fSig500     = TFile::Open("all_ExLepMuMuZ_M500.root");
TFile *fSig750     = TFile::Open("all_ExLepMuMuZ_M750.root");
TFile *fSig1000     = TFile::Open("all_ExLepMuMuZ_M1000.root");
TFile *fSig1250     = TFile::Open("all_ExLepMuMuZ_M1250.root");
TFile *fSig1500     = TFile::Open("all_ExLepMuMuZ_M1500.root");
TFile *fSig1750     = TFile::Open("all_ExLepMuMuZ_M1750.root");
TFile *fSig2000     = TFile::Open("all_ExLepMuMuZ_M2000.root");
TFile *fSig2500     = TFile::Open("all_ExLepMuMuZ_M2500.root");
TFile *fSig3000     = TFile::Open("all_ExLepMuMuZ_M3000.root");
TFile *fSig3500     = TFile::Open("all_ExLepMuMuZ_M3500.root");
TFile *fSig4000     = TFile::Open("all_ExLepMuMuZ_M4000.root");
TFile *fSig4500     = TFile::Open("all_ExLepMuMuZ_M4500.root");
TFile *fSig5000     = TFile::Open("all_ExLepMuMuZ_M5000.root");

TH1F* getHisto(TFile *histFile, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  TString fullPath = histPath+"/"+dir+"/"+histName;
  if(!(histFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(histFile->Get(fullPath))->Clone(histName);
  return hist;
}

void printOutput(TFile *fBkg, TFile *fSig,  string mass){
  TH1F * hNoChBkg = fBkg->Get("base/noCharge");
  TH1F * hOppChBkg = fBkg->Get("base/oppCharge");
  TH1F * hSameChBkg = fBkg->Get("base/sameCharge");
  double noChBkg = hNoChBkg->Integral();
  double oppChBkg = hOppChBkg->Integral();
  double sameChBkg = hSameChBkg->Integral();

  TH1F * hNoChSig = fSig->Get("base/noCharge");
  TH1F * hOppChSig = fSig->Get("base/oppCharge");
  TH1F * hSameChSig = fSig->Get("base/sameCharge");
  double noChSig = hNoChSig->Integral();
  double oppChSig = hOppChSig->Integral();
  double sameChSig = hSameChSig->Integral();
  cout<<"--------------------------------------------------------------"<<endl;
  cout<<std::setprecision(4)<<endl;
  cout<<setw(20)<<"Charge Selection"<<setw(15)<<"Siganl,M"+mass<<setw(15)<<"AllBkg"<<setw(15)<<"Sig/AllBkg"<<endl;
  cout<<setw(20)<<"NoChargeSel"<<setw(15)<<noChSig<<setw(15)<<noChBkg<<setw(15)<<noChSig/noChBkg<<endl;
  cout<<setw(20)<<"OppChargeSel"<<setw(15)<<oppChSig<<setw(15)<<oppChBkg<<setw(15)<<oppChSig/oppChBkg<<endl;
  cout<<setw(20)<<"SameChargeSel"<<setw(15)<<sameChSig<<setw(15)<<sameChBkg<<setw(15)<<sameChSig/sameChBkg<<endl;
}
void oppChargeSel(){
  printOutput(fBkg, fSig250, "250");
  printOutput(fBkg, fSig500, "500");
  printOutput(fBkg, fSig750, "750");
  printOutput(fBkg, fSig1000, "1000");
  printOutput(fBkg, fSig1250, "1250");
  printOutput(fBkg, fSig1500, "1500");
  printOutput(fBkg, fSig1750, "1750");
  printOutput(fBkg, fSig2000, "2000");
  printOutput(fBkg, fSig2500, "2500");
  printOutput(fBkg, fSig3000, "3000");
  printOutput(fBkg, fSig3500, "3500");
  printOutput(fBkg, fSig4000, "4000");
  printOutput(fBkg, fSig4500, "4500");
  printOutput(fBkg, fSig5000, "5000");
}



