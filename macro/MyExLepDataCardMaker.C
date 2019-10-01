#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> 
#include "MyExLepDataCardMaker.h"

//----------------------------------------//
//make data card for each mass
//----------------------------------------//
MyExLepDataCardMaker DC;
void MyExLepDataCardMaker(TString inFileDir="stack_tmp_Mu_Sys", 
        TString histSubDir_="ZTag1", 
        TString histName="mlZ_max_sig250", 
        TString channelName="mu", 
        int mass = 250, 
        TString label="ExLep250", 
        TString hPlusFileName="all_ExLep_M250.root")
  {
  TString histSubDir = "Iso/"+histSubDir_+"/";
  bool isMuChannel = false; 
  if(channelName=="mu") isMuChannel = true;
  
  TFile *fLstar  = TFile::Open(inFileDir+"/"+hPlusFileName);
  
  cout<<" ======> mass point: "<<mass<<endl;
  //ttbar
  double sf_ttbar = 1; 
  TH1F* hTTbar = DC.readHisto(fTT, "base/", histSubDir, histName, sf_ttbar, "TTbar"); 
  TH1F* hTTbar_JESUp = DC.readHisto(fTT, "JESPlus/", histSubDir, histName, sf_ttbar, "TTbar");
  TH1F* hTTbar_JESDown = DC.readHisto(fTT, "JESMinus/", histSubDir, histName, sf_ttbar, "TTbar"); 
  TH1F* hTTbar_JERUp = DC.readHisto(fTT, "JERPlus/", histSubDir, histName, sf_ttbar, "TTbar"); 
  TH1F* hTTbar_JERDown = DC.readHisto(fTT, "JERMinus/", histSubDir, histName, sf_ttbar, "TTbar"); 
  TH1F* hTTbar_bTagUp = DC.readHisto(fTT, "bTagPlus/", histSubDir, histName, sf_ttbar, "TTbar");
  TH1F* hTTbar_bTagDown = DC.readHisto(fTT, "bTagMinus/", histSubDir, histName, sf_ttbar, "TTbar"); 
  
  //Z+Jets
  double sf_dyjet = 1;
  TH1F* hDYJet = DC.readHisto(fDY, "base/", histSubDir, histName, sf_dyjet, "DYJets"); 
  TH1F* hDYJet_JESUp = DC.readHisto(fDY, "JESPlus/", histSubDir, histName, sf_dyjet, "DYJets"); 
  TH1F* hDYJet_JESDown = DC.readHisto(fDY, "JESMinus/", histSubDir, histName, sf_dyjet, "DYJets"); 
  TH1F* hDYJet_JERUp = DC.readHisto(fDY, "JERPlus/", histSubDir, histName, sf_dyjet, "DYJets"); 
  TH1F* hDYJet_JERDown = DC.readHisto(fDY, "JERMinus/", histSubDir, histName, sf_dyjet, "DYJets"); 
  TH1F* hDYJet_bTagUp = DC.readHisto(fDY, "bTagPlus/", histSubDir, histName, sf_dyjet, "DYJets"); 
  TH1F* hDYJet_bTagDown = DC.readHisto(fDY, "bTagMinus/", histSubDir, histName, sf_dyjet, "DYJets"); 
  
  //Dibosons
  double sf_vv = 1;
  TH1F* hVV = DC.readHisto(fVV, "base/", histSubDir, histName, sf_vv, "VV"); 
  TH1F* hVV_JESUp = DC.readHisto(fVV, "JESPlus/", histSubDir, histName, sf_vv, "VV"); 
  TH1F* hVV_JESDown = DC.readHisto(fVV, "JESMinus/", histSubDir, histName, sf_vv, "VV"); 
  TH1F* hVV_JERUp = DC.readHisto(fVV, "JERPlus/", histSubDir, histName, sf_vv, "VV"); 
  TH1F* hVV_JERDown = DC.readHisto(fVV, "JERMinus/", histSubDir, histName, sf_vv, "VV"); 
  TH1F* hVV_bTagUp = DC.readHisto(fVV, "bTagPlus/", histSubDir, histName, sf_vv, "VV"); 
  TH1F* hVV_bTagDown = DC.readHisto(fVV, "bTagMinus/", histSubDir, histName, sf_vv, "VV"); 
  //Data
  double sf_data = 1; //should be 1, always
  TH1F* hData = DC.readHisto(fData, "base/", histSubDir, histName, sf_data, "Data"); 

  //wh
  //lstar
  double sf_lstar = 1; 
  TH1F* hLstar = DC.readHisto(fLstar, "base/", histSubDir, histName, sf_lstar, "Signal"); 
  TH1F* hLstar_JESUp = DC.readHisto(fLstar, "JESPlus/", histSubDir, histName, sf_lstar, "Signal"); 
  TH1F* hLstar_JESDown = DC.readHisto(fLstar, "JESMinus/", histSubDir, histName, sf_lstar, "Signal"); 
  TH1F* hLstar_JERUp = DC.readHisto(fLstar, "JERPlus/", histSubDir, histName, sf_lstar, "Signal"); 
  TH1F* hLstar_JERDown = DC.readHisto(fLstar, "JERMinus/", histSubDir, histName, sf_lstar, "Signal"); 
  TH1F* hLstar_bTagUp = DC.readHisto(fLstar, "bTagPlus/", histSubDir, histName, sf_lstar, "Signal"); 
  TH1F* hLstar_bTagDown = DC.readHisto(fLstar, "bTagMinus/", histSubDir, histName, sf_lstar, "Signal"); 

  //open input template data card
  ifstream in;
  char* c = new char[1000];
  in.open("MyTemplateDataCard.txt");
  //create output data card for 13 TeV
  string outDataCard = "datacard_llstar_llZ_llq_13TeV_ExLep.txt";
  string histName_str(histSubDir_+TString("_")+histName);
  //if(isMuChannel) outDataCard = "datacard_llstar_llZ_llq_13TeV_mu_"+histName_str+"_ExLep%d.txt"; 
  //else outDataCard = "datacard_llstar_llZ_llq_13TeV_ele_"+histName_str+"_ExLep%d.txt";
  if(isMuChannel) outDataCard = "datacard_llstar_llZ_llq_13TeV_mu_"+histName_str+".txt"; 
  else outDataCard = "datacard_llstar_llZ_llq_13TeV_ele_"+histName_str+".txt";
  ofstream out(Form(outDataCard.c_str(), mass));
  out.precision(8);

  time_t secs=time(0);
  tm *t=localtime(&secs);
  while (in.good()){
    in.getline(c,1000,'\n');
    if (in.good()){
      string line(c);
      if(line.find("Date")!=string::npos){
        string day = string(Form("%d",t->tm_mday));
        string month = string(Form("%d",t->tm_mon+1));
        string year = string(Form("%d",t->tm_year+1900));
        line.replace( line.find("XXX") , 3 , day+"/"+month+"/"+year);
        out << line << endl;
      }
      else if(line.find("Description")!=string::npos){
        line.replace( line.find("YYY") , 3 , string(Form("%d", mass)) );
        line.replace( line.find("ZZZ") , 3 , string(Form("%f", totLumi)) ); 
        out << line << endl;
      }
      else if(line.find("Observation")!=string::npos){
        line.replace( line.find("XXX") , 3 , string(Form("%.0f", hData->Integral())));
        out << line << endl;
      }
      else if(line.find("process")!=string::npos && line.find("lstar")!=string::npos){
        line.replace( line.find("YYY") , 3 , string(Form("%d", mass)) );
        out << line << endl;
      }
      else if(line.find("rate")!=string::npos){
        string rate = "rate               ";  
        string space = "     ";
        out << rate ;
	out << space << hLstar->Integral()
            << space << hTTbar->Integral()
            << space << hDYJet->Integral()
            << space << hVV->Integral()
            << endl;
      }
      else if(line.find("CMS_stat_lstar")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", DC.getStatUnc(hLstar,  0))));  
        out << line << endl;
      } 
      else if(line.find("CMS_stat_tt")!=string::npos){  
	line.replace( line.find("XXXX") , 4 , string(Form("%.2f", DC.getStatUnc(hTTbar,  0))));   
        out << line << endl;
      }  
      else if(line.find("CMS_stat_zjet")!=string::npos){ 
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", DC.getStatUnc(hDYJet,  0))));  
        out << line << endl; 
      }
      else if(line.find("CMS_stat_vv")!=string::npos){  
        line.replace( line.find("XXXX") , 4 , string(Form("%.2f", DC.getStatUnc(hVV,  0))));   
        out << line << endl;  
      }
     else if(line.find("CMS_scale_j")!=string::npos){
        float JESUnc_lstar = (hLstar->Integral() > 0) ? DC.getSysUnc(hLstar, hLstar_JESUp, hLstar_JESDown) : 1.00;
        line.replace( line.find("LLLL") , 4 , string(Form("%.3f", JESUnc_lstar)) );

        float JESUnc_ttbar = (hTTbar->Integral() > 0) ? DC.getSysUnc(hTTbar, hTTbar_JESUp, hTTbar_JESDown) : 1.00;
        line.replace( line.find("TTTT") , 4 , string(Form("%.3f", JESUnc_ttbar)) );

        float JESUnc_zjet = (hDYJet->Integral() > 0) ? DC.getSysUnc(hDYJet, hDYJet_JESUp, hDYJet_JESDown) : 1.00;
        line.replace( line.find("DDDD") , 4 , string(Form("%.3f", JESUnc_zjet)) );

        float JESUnc_vv = (hVV->Integral() > 0) ? DC.getSysUnc(hVV, hVV_JESUp, hVV_JESDown) : 1.00;
        line.replace( line.find("VVVV") , 4 , string(Form("%.3f", JESUnc_vv)) );
        out << line << endl;
      }
     else if(line.find("CMS_res_j")!=string::npos){
        float JERUnc_hLstar = (hLstar->Integral() > 0) ? DC.getSysUnc(hLstar, hLstar_JERUp, hLstar_JERDown) : 1.00;
        line.replace( line.find("LLLL") , 4 , string(Form("%.3f", JERUnc_hLstar)) );

        float JERUnc_ttbar = (hTTbar->Integral() > 0) ? DC.getSysUnc(hTTbar, hTTbar_JERUp, hTTbar_JERDown) : 1.00;
        line.replace( line.find("TTTT") , 4 , string(Form("%.3f", JERUnc_ttbar)) );

        float JERUnc_zjet = (hDYJet->Integral() > 0) ? DC.getSysUnc(hDYJet, hDYJet_JERUp, hDYJet_JERDown) : 1.00;
        line.replace( line.find("DDDD") , 4 , string(Form("%.3f", JERUnc_zjet)) );

        float JERUnc_vv = (hVV->Integral() > 0) ? DC.getSysUnc(hVV, hVV_JERUp, hVV_JERDown) : 1.00;
        line.replace( line.find("VVVV") , 4 , string(Form("%.3f", JERUnc_vv)) );
        out << line << endl;
      }
     else if(line.find("CMS_eff_b")!=string::npos){
        float bTagUnc_hLstar = (hLstar->Integral() > 0) ? DC.getSysUnc(hLstar, hLstar_bTagUp, hLstar_bTagDown) : 1.00;
        line.replace( line.find("LLLL") , 4 , string(Form("%.3f", bTagUnc_hLstar)) );

        float bTagUnc_ttbar = (hTTbar->Integral() > 0) ? DC.getSysUnc(hTTbar, hTTbar_bTagUp, hTTbar_bTagDown) : 1.00;
        line.replace( line.find("TTTT") , 4 , string(Form("%.3f", bTagUnc_ttbar)) );

        float bTagUnc_zjet = (hDYJet->Integral() > 0) ? DC.getSysUnc(hDYJet, hDYJet_bTagUp, hDYJet_bTagDown) : 1.00;
        line.replace( line.find("DDDD") , 4 , string(Form("%.3f", bTagUnc_zjet)) );

        float bTagUnc_vv = (hVV->Integral() > 0) ? DC.getSysUnc(hVV, hVV_bTagUp, hVV_bTagDown) : 1.00;
        line.replace( line.find("VVVV") , 4 , string(Form("%.3f", bTagUnc_vv)) );
        out << line << endl;
      }
      else{ //default without changes
        out << line << endl;
      }
    }
  } 
  out.close();
  in.close();
}
