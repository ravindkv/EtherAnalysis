#include <iostream>
#include <fstream>
#include <iomanip>

//INPUT FILES
TFile* fMuData = TFile::Open("inputFiles/all_muData.root");
TFile* fEleData = TFile::Open("inputFiles/all_eleData.root");
TFile* fVV	= TFile::Open("inputFiles/all_VV.root");
TFile* fDY	= TFile::Open("inputFiles/all_DY.root");
TFile* fTT	= TFile::Open("inputFiles/all_TT.root");
TFile* fWJ	= TFile::Open("inputFiles/all_WJets.root");
TFile *fSigMuMuZ250      = TFile::Open("inputFiles/all_ExLepMuMuZ_M250.root");
TFile *fSigMuMuZ500      = TFile::Open("inputFiles/all_ExLepMuMuZ_M500.root");
TFile *fSigMuMuZ750      = TFile::Open("inputFiles/all_ExLepMuMuZ_M750.root");
TFile *fSigMuMuZ1000     = TFile::Open("inputFiles/all_ExLepMuMuZ_M1000.root");
TFile *fSigMuMuZ1250     = TFile::Open("inputFiles/all_ExLepMuMuZ_M1250.root");
TFile *fSigMuMuZ1500     = TFile::Open("inputFiles/all_ExLepMuMuZ_M1500.root");
TFile *fSigMuMuZ1750     = TFile::Open("inputFiles/all_ExLepMuMuZ_M1750.root");
TFile *fSigMuMuZ2000     = TFile::Open("inputFiles/all_ExLepMuMuZ_M2000.root");
TFile *fSigMuMuZ2500     = TFile::Open("inputFiles/all_ExLepMuMuZ_M2500.root");
TFile *fSigMuMuZ3000     = TFile::Open("inputFiles/all_ExLepMuMuZ_M3000.root");
TFile *fSigMuMuZ3500     = TFile::Open("inputFiles/all_ExLepMuMuZ_M3500.root");
TFile *fSigMuMuZ4000     = TFile::Open("inputFiles/all_ExLepMuMuZ_M4000.root");
TFile *fSigMuMuZ4500     = TFile::Open("inputFiles/all_ExLepMuMuZ_M4500.root");
TFile *fSigMuMuZ5000     = TFile::Open("inputFiles/all_ExLepMuMuZ_M5000.root");

TFile *fSigEEZ250      = TFile::Open("inputFiles/all_ExLepEEZ_M250.root");
TFile *fSigEEZ500      = TFile::Open("inputFiles/all_ExLepEEZ_M500.root");
TFile *fSigEEZ750      = TFile::Open("inputFiles/all_ExLepEEZ_M750.root");
TFile *fSigEEZ1000     = TFile::Open("inputFiles/all_ExLepEEZ_M1000.root");
TFile *fSigEEZ1250     = TFile::Open("inputFiles/all_ExLepEEZ_M1250.root");
TFile *fSigEEZ1500     = TFile::Open("inputFiles/all_ExLepEEZ_M1500.root");
TFile *fSigEEZ1750     = TFile::Open("inputFiles/all_ExLepEEZ_M1750.root");
TFile *fSigEEZ2000     = TFile::Open("inputFiles/all_ExLepEEZ_M2000.root");
TFile *fSigEEZ2500     = TFile::Open("inputFiles/all_ExLepEEZ_M2500.root");
TFile *fSigEEZ3000     = TFile::Open("inputFiles/all_ExLepEEZ_M3000.root");
TFile *fSigEEZ3500     = TFile::Open("inputFiles/all_ExLepEEZ_M3500.root");
TFile *fSigEEZ4000     = TFile::Open("inputFiles/all_ExLepEEZ_M4000.root");
TFile *fSigEEZ4500     = TFile::Open("inputFiles/all_ExLepEEZ_M4500.root");
TFile *fSigEEZ5000     = TFile::Open("inputFiles/all_ExLepEEZ_M5000.root");

string doubleToStr(double val){
     ostringstream convert;
     string result("");
     convert <<std::setprecision(4)<<val;
     result = convert.str();
  return result;
}

string makeCutFlowRow(TFile *file, string rowLable, bool isMu, bool isEle){
  TString chName = "";
  if(isMu) chName = "Muon";
  if(isEle) chName = "Electron";
  TH1F * hCutFlow = (TH1F*)file->Get(chName+"/base/cutflow");
  TH1F * hTotalEvents = (TH1F*)file->Get(chName+"/base/totalYield");
  double totalEvents = hTotalEvents->Integral();
  cout<<rowLable<<"\t"<<totalEvents<<endl;
  string row = rowLable;
  for(int i=1; i<hCutFlow->GetSize()-1; i++){
    if(i<=5){
    //if(hCutFlow->GetBinContent(i)!=0){
      double preYield = 1.0; //previous yield
      if (i==1) preYield = totalEvents; 
      else preYield = hCutFlow->GetBinContent(i-1);
      double currYield = 1.0; //current yueld
      currYield = hCutFlow->GetBinContent(i);
      double eff = 100*currYield/preYield;
      row += " & "+ doubleToStr(currYield)+" ("+doubleToStr(eff)+")";
    }
  }
  row += "\\\\";
  //cout<<row<<endl;
  return row;
}
string makeCutFlowRowRatio(TFile *fData, string rowLableBkg, string rowLableRatio, bool isRatio, bool isMu, bool isEle){
  TString chName = "";
  if(isMu) chName = "Muon";
  if(isEle) chName = "Electron";
  TH1F * hCutFlowData = (TH1F*)fData->Get(chName+"/base/cutflow");
  //add all bkgs
  TH1F * hEventsTT = (TH1F*)fTT->Get(chName+"/base/totalYield");
  TH1F * hEventsDY = (TH1F*)fDY->Get(chName+"/base/totalYield");
  TH1F * hEventsWJ = (TH1F*)fWJ->Get(chName+"/base/totalYield");
  TH1F * hEventsVV = (TH1F*)fVV->Get(chName+"/base/totalYield");
  TH1F * hEventsBkg = (TH1F*)hEventsTT->Clone();
  hEventsBkg->Add(hEventsDY);
  hEventsBkg->Add(hEventsWJ);
  hEventsBkg->Add(hEventsVV);
  
  TH1F * hCutFlowTT = (TH1F*)fTT->Get(chName+"/base/cutflow");
  TH1F * hCutFlowDY = (TH1F*)fDY->Get(chName+"/base/cutflow");
  TH1F * hCutFlowWJ = (TH1F*)fWJ->Get(chName+"/base/cutflow");
  TH1F * hCutFlowVV = (TH1F*)fVV->Get(chName+"/base/cutflow");
  TH1F * hCutFlowBkg = (TH1F*)hCutFlowTT->Clone();
  hCutFlowBkg->Add(hCutFlowDY);
  hCutFlowBkg->Add(hCutFlowWJ);
  hCutFlowBkg->Add(hCutFlowVV);

  //print row for all bkgs
  double totalEvents = hEventsBkg->Integral();
  string rowBkg = rowLableBkg;
  for(int i=1; i<hCutFlowBkg->GetSize()-1; i++){
    if(i<=5){
    //if(hCutFlowBkg->GetBinContent(i)!=0){
      double preYield = 1.0; //previous yield
      if (i==1) preYield = totalEvents; 
      else preYield = hCutFlowBkg->GetBinContent(i-1);
      double currYield = 1.0; //current yueld
      currYield = hCutFlowBkg->GetBinContent(i);
      double eff = 100*currYield/preYield;
      rowBkg += " & "+ doubleToStr(currYield)+" ("+doubleToStr(eff)+")";
    }
  }
  rowBkg += "\\\\";

  //print row for data/MC ratio
  string rowRatio = rowLableRatio;
  for(int i=1; i<hCutFlowBkg->GetSize()-1; i++){
    if(i<=5){
      double yieldBkg = 1.0; //current yueld
      double yieldData = 1.0; //current yueld
      yieldBkg = hCutFlowBkg->GetBinContent(i);
      yieldData = hCutFlowData->GetBinContent(i);
      rowRatio += " & "+ doubleToStr(yieldData/yieldBkg);
    }
  }
  rowRatio += "\\\\";

  //return 
  if(isRatio){
    cout<<rowRatio<<endl;
    return rowRatio;
  }
  else{
    //cout<<rowBkg<<endl;
    return rowBkg;
  }
}

void makeTable(ofstream & outFile, bool isMu, bool isEle){
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outFile<<"\\begin{table}"<<endl;
  outFile<<"\\begin{center}"<<endl;  
  //outFile<<"\\begin{LARGE}"<<endl;  
  outFile<<"\\begin{tabular}{ |c|c|c|c|c|c|}"<<endl;  
  outFile<<"\\hline "<<endl; 
  if(isMu) outFile<<" & "<<"Muon trigger & "<<"ControlSel & "<<"b-jet Veto & "<<"PreSel & "<<"ZTag \\\\"<<endl;
  if(isEle) outFile<<" & "<<"Electron trigger & "<<"ControlSel & "<<"b-jet Veto & "<<"PreSel & "<<"ZTag \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"Process & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ & "<<"$N_{evt}(100\\times \\epsilon)$ \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  if(isMu){
    outFile<<makeCutFlowRow(fSigMuMuZ250 , "$M_{l^*} $ = 250  GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ500 , "$M_{l^*} $ = 500  GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ750 , "$M_{l^*} $ = 750  GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ1000, "$M_{l^*} $ = 1000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ1250, "$M_{l^*} $ = 1250 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ1500, "$M_{l^*} $ = 1500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ1750, "$M_{l^*} $ = 1750 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ2000, "$M_{l^*} $ = 2000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ2500, "$M_{l^*} $ = 2500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ3000, "$M_{l^*} $ = 3000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ3500, "$M_{l^*} $ = 3500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ4000, "$M_{l^*} $ = 4000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ4500, "$M_{l^*} $ = 4500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigMuMuZ5000, "$M_{l^*} $ = 5000 GeV", isMu, isEle)<<endl;
  }
  if(isEle){
    outFile<<makeCutFlowRow(fSigEEZ250 , "$M_{l^*} $ = 250  GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ500 , "$M_{l^*} $ = 500  GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ750 , "$M_{l^*} $ = 750  GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ1000, "$M_{l^*} $ = 1000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ1250, "$M_{l^*} $ = 1250 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ1500, "$M_{l^*} $ = 1500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ1750, "$M_{l^*} $ = 1750 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ2000, "$M_{l^*} $ = 2000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ2500, "$M_{l^*} $ = 2500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ3000, "$M_{l^*} $ = 3000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ3500, "$M_{l^*} $ = 3500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ4000, "$M_{l^*} $ = 4000 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ4500, "$M_{l^*} $ = 4500 GeV", isMu, isEle)<<endl;
    outFile<<makeCutFlowRow(fSigEEZ5000, "$M_{l^*} $ = 5000 GeV", isMu, isEle)<<endl;
  }
  outFile<<"\\hline "<<endl;
  outFile<<makeCutFlowRow(fTT, "$t\\bar{t}$ + jets", isMu, isEle)<<endl;
  outFile<<makeCutFlowRow(fDY, "$Z/\\gamma$+ jets",  isMu, isEle)<<endl;
  outFile<<makeCutFlowRow(fVV, "$VV$", isMu, isEle)<<endl;
  outFile<<makeCutFlowRow(fWJ, "$W$ + jets", isMu, isEle)<<endl;
  outFile<<"\\hline "<<endl;
  if(isMu)outFile<<makeCutFlowRowRatio(fMuData, "All Bkg", "Data/MC", false, isMu,isEle)<<endl;
  if(isEle)outFile<<makeCutFlowRowRatio(fEleData, "All Bkg", "Data/MC", false,isMu,isEle)<<endl;
  outFile<<"\\hline "<<endl;
  if(isMu)outFile<<makeCutFlowRow(fMuData, "Data", isMu, isEle)<<endl;
  if(isEle)outFile<<makeCutFlowRow(fEleData, "Data", isMu, isEle)<<endl;
  outFile<<"\\hline "<<endl;
  if(isMu)outFile<<makeCutFlowRowRatio(fMuData,"All Bkg","Data/MC",true,isMu, isEle)<<endl;
  if(isEle)outFile<<makeCutFlowRowRatio(fEleData,"All Bkg","Data/MC",true,isMu, isEle)<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  //outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;
  string chName = "muon";
  if(isEle) chName = "electron";
  outFile<<"\\caption{Cutflow table for "+chName+" channel. The $\\epsilon$ is the efficiency $(N_{after}/N_{before})$ and $N_{evt} (=N_{after})$ is the number of survived events after different cuts.}"<<endl;
  outFile<<"\\end{table}"<<endl;

}

void MyCutFlowTable(){  
  ofstream outFile; 
  outFile.open("outputFiles/cutflowTable.tex"); 
  outFile<<"\\documentclass[]{article}"<<endl;  
  outFile<<"\\pagestyle{empty}"<<endl;  
  outFile<<"\\usepackage{epsfig}"<<endl;  
  outFile<<"\\usepackage{amsmath}"<<endl;  
  outFile<<"\\usepackage{array}"<<endl;  
  outFile<<"\\usepackage{multirow}"<<endl;  
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;  
  outFile<<"\\textheight = 8.in"<<endl;  
  outFile<<"\\textwidth 7.0in"<<endl;  
  outFile<<"\\begin{document}"<<endl;  
  outFile<<""<<endl;
  makeTable(outFile, true, false);
  makeTable(outFile, false, true);
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
} 
