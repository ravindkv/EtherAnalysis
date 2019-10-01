#include <iostream>
#include <fstream>
#include <iomanip>

//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;
//TString zTagDir = "ZTag1";

//INPUT FILES
TFile* fData = TFile::Open("all_Data.root");
TFile* fVV	= TFile::Open("all_VV.root");
//TFile* fDY	= TFile::Open("all_DY_dd.root");
TFile* fDY	= TFile::Open("all_DY.root");
TFile* fTT	= TFile::Open("all_TT.root");
TFile* fWJ	= TFile::Open("all_WJets.root");
TFile* fBkg	= TFile::Open("all_NonDYBkg.root");
//TFile* fBkg	= TFile::Open("all_Bkg.root");


TH1F* getHisto(TFile *histFile, TString histPath, TString dir, TString histName){
  TH1F* hist; 
  TString fullPath = histPath+"/"+dir+"/"+histName;
  if(!(histFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get("base/Iso/ZTag/mll")->Clone(histName));
    hist->Reset();
  }else hist = (TH1F*)(histFile->Get(fullPath))->Clone(histName);
  return hist;
}

//get statistical uncertainity
double getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  //double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;
  double statUnc = fabs(sError);
  return statUnc;
}

string doubleToStr(double val){
     ostringstream convert;
     string result("");
     //convert <<std::setprecision(2)<<std::scientific<<val;
     convert <<std::setprecision(3)<<val;
     result = convert.str();
  return result;
}

string makeOneRow(string procName, TFile *inFile, TString zTagDir, TString mass){
  TH1F * hA = (TH1F*)getHisto(inFile, "base/Iso", zTagDir, "A_mlZ_max_sig"+mass);
  TH1F * hB = (TH1F*)getHisto(inFile, "base/Iso", zTagDir, "B_mlZ_max_sig"+mass);
  TH1F * hC = (TH1F*)getHisto(inFile, "base/Iso", zTagDir, "C_mlZ_max_sig"+mass);
  TH1F * hD = (TH1F*)getHisto(inFile, "base/Iso", zTagDir, "D_mlZ_max_sig"+mass);
  double yieldA = hA->Integral();
  double yieldB = hB->Integral();
  double yieldC = hC->Integral();
  double yieldD = hD->Integral();
  double errorA = getStatUnc(hA, 0.0);
  double errorB = getStatUnc(hB, 0.0);
  double errorC = getStatUnc(hC, 0.0);
  double errorD = getStatUnc(hD, 0.0);
  double yieldBkg = yieldB*yieldC/yieldD;
  double errorBkg = 0.0;
  return procName+" & "+doubleToStr(yieldA)+" $\\pm$ "+doubleToStr(errorA)+" & "+
	  doubleToStr(yieldB)+" $\\pm $ "+doubleToStr(errorB)+" & "+
	  doubleToStr(yieldC)+" $\\pm $ "+doubleToStr(errorC)+" & "+
	  doubleToStr(yieldD)+" $\\pm $ "+doubleToStr(errorD);
	  //+" & "+ doubleToStr(yieldBkg)+" $\\pm $ "+doubleToStr(errorBkg);
}

void makeOneTable(ofstream & outFile, TString zTagDir, TString mass){  
  outFile<<"\\begin{table}"<<endl;
  outFile<<"\\begin{center}"<<endl;  
  //outFile<<"\\begin{LARGE}"<<endl;  
  outFile<<"\\begin{tabular}{cccccc}"<<endl;  
  outFile<<"\\hline "<<endl; 
  outFile<<"Process & $N_A \\pm stat$ &  $N_B \\pm stat$ & $N_C \\pm stat$ & $N_D \\pm stat$ \\\\"<<endl; // & $\\frac{N_B\\times N_C}{N_D}$\\\\"<<endl;
  outFile<<" & ($>200, <0.60$) & ($<200, <0.60$) &($>200, >0.60$) & ($<200, >0.60$)  \\\\"<<endl; 
  outFile<<"\\hline  \\noalign{\\vskip 0.1cm}"<<endl; 
  outFile<<makeOneRow("MC $Z/\\gamma$ + jets", fDY, zTagDir, mass)<<" \\\\ \\noalign{\\vskip 0.1cm}"<<endl;
  outFile<<"\\hline  \\noalign{\\vskip 0.1cm}"<<endl; 
  outFile<<makeOneRow("$ t\\bar{t}$ + jets ", fTT, zTagDir, mass)<<" \\\\ \\noalign{\\vskip 0.1cm}"<<endl;
  outFile<<makeOneRow("VV + jets", fVV, zTagDir, mass)<<" \\\\ \\noalign{\\vskip 0.1cm}"<<endl;
  outFile<<makeOneRow("Non DY Bkg", fBkg, zTagDir, mass)<<" \\\\ \\noalign{\\vskip 0.1cm}"<<endl;
  outFile<<"\\hline  \\noalign{\\vskip 0.1cm}"<<endl; 
  outFile<<makeOneRow(" Data", fData, zTagDir, mass)<<" \\\\ \\noalign{\\vskip 0.1cm}"<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  //outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;
  string chName = "muon";
  if(isEleChannel) chName = "electron";
  outFile<<"\\caption{Event yields in the ABCD region formed by ($M_{l_{1}l_{2}} (\\rm{GeV)}, \\tau_{21})$ for the excited lepton mass  " +mass+ " GeV, after L-cut for "+chName+" channel.}"<<endl;
  outFile<<"\\end{table}"<<endl;
} 

void MyBkgEstTable(){  
  ofstream outFile; 
  outFile.open("bkgEstTable.tex"); 
  outFile<<"\\documentclass[]{article}"<<endl;  
  outFile<<"\\pagestyle{empty}"<<endl;  
  outFile<<"\\usepackage{epsfig}"<<endl;  
  outFile<<"\\usepackage{amsmath}"<<endl;  
  outFile<<"\\begin{document}"<<endl;  
  outFile<<""<<endl;
  makeOneTable(outFile, "ZTag3", "250");
  makeOneTable(outFile, "ZTag3", "500");
  makeOneTable(outFile, "ZTag3", "750");
  makeOneTable(outFile, "ZTag3", "1000");
  makeOneTable(outFile, "ZTag3", "1250");
  makeOneTable(outFile, "ZTag3", "1500");
  makeOneTable(outFile, "ZTag3", "1750");
  makeOneTable(outFile, "ZTag3", "2000");
  makeOneTable(outFile, "ZTag3", "2500");
  makeOneTable(outFile, "ZTag3", "3000");
  makeOneTable(outFile, "ZTag3", "3500");
  makeOneTable(outFile, "ZTag3", "4000");
  makeOneTable(outFile, "ZTag3", "4500");
  makeOneTable(outFile, "ZTag3", "5000");
  outFile<<"\\end{document}"<<endl; 
  outFile.close(); 
} 


