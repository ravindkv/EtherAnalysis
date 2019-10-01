//#include "TROOT.h"    
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <cmath>
#include <math.h>
#include<string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <fstream>
#include <map>
#include <iomanip>
#include <iostream>

void MyCheckProEvent(){

  bool is2016Data = true;
  bool is2017Data = false;

  if(is2016Data){	
    bool isMC = true;
    bool isMuMC = true;
    bool isMuData = true;
    bool isEleMC = false;
    bool isEleData = false;
  
    std::map<std::string, double> mcEvtDBS;
    mcEvtDBS["DYJetsToLL_M50_Merged.root"]         = 118692000; 
    mcEvtDBS["DYJetsToLL_Pt100To250_Merged.root"]              = 2040596 ; 
    mcEvtDBS["DYJetsToLL_Pt250To400_Merged.root"]              = 423976  ;
    mcEvtDBS["DYJetsToLL_Pt400To650_Merged.root"]              = 432056  ;
    mcEvtDBS["DYJetsToLL_Pt50To100_Merged.root"]               = 21890432;
    mcEvtDBS["DYJetsToLL_Pt650ToInf_Merged.root"]              = 430691  ;
    mcEvtDBS["TT_Merged.root"] 		           = 77081156;
    mcEvtDBS["WJetsToLNu_Merged.root"]	           = 29017600;
    mcEvtDBS["WW_Merged.root"]		           = 994012;
    mcEvtDBS["WZ_Merged.root"]		           = 1000000;
    mcEvtDBS["ZZ_Merged.root"] 		           = 990064;
  						
    if(isMuMC){
       mcEvtDBS["ExLepMuMuZ_M1000_Merged.root"]        = 191300;
       mcEvtDBS["ExLepMuMuZ_M1250_Merged.root"]        = 181400;
       mcEvtDBS["ExLepMuMuZ_M1500_Merged.root"]        = 190900;
       mcEvtDBS["ExLepMuMuZ_M1750_Merged.root"]        = 192200;
       mcEvtDBS["ExLepMuMuZ_M2000_Merged.root"]        = 200000;
       mcEvtDBS["ExLepMuMuZ_M250_Merged.root"]         = 198200;
       mcEvtDBS["ExLepMuMuZ_M2500_Merged.root"]        = 200000;
       mcEvtDBS["ExLepMuMuZ_M3000_Merged.root"]        = 200000;
       mcEvtDBS["ExLepMuMuZ_M3500_Merged.root"]        = 200000;
       mcEvtDBS["ExLepMuMuZ_M4000_Merged.root"]        = 192300;
       mcEvtDBS["ExLepMuMuZ_M4500_Merged.root"]        = 200000;
       mcEvtDBS["ExLepMuMuZ_M500_Merged.root"]         = 200000;
       mcEvtDBS["ExLepMuMuZ_M5000_Merged.root"]        = 192500;
       mcEvtDBS["ExLepMuMuZ_M750_Merged.root"]         = 191400;
    }
    if(isEleMC){
      mcEvtDBS["ExLepEEZ_M1000_Merged.root"]            = 197600;
      mcEvtDBS["ExLepEEZ_M1250_Merged.root"]            = 200000;
      mcEvtDBS["ExLepEEZ_M1500_Merged.root"]            = 41900 ;
      mcEvtDBS["ExLepEEZ_M1750_Merged.root"]            = 192400;
      mcEvtDBS["ExLepEEZ_M2000_Merged.root"]            = 189000;
      mcEvtDBS["ExLepEEZ_M250_Merged.root"]            = 182100;
      mcEvtDBS["ExLepEEZ_M2500_Merged.root"]            = 187400;
      mcEvtDBS["ExLepEEZ_M3000_Merged.root"]            = 192800;
      mcEvtDBS["ExLepEEZ_M3500_Merged.root"]            = 170000;
      mcEvtDBS["ExLepEEZ_M4000_Merged.root"]            = 199800;
      mcEvtDBS["ExLepEEZ_M4500_Merged.root"]            = 187400;
      mcEvtDBS["ExLepEEZ_M500_Merged.root"]             = 200000;
      mcEvtDBS["ExLepEEZ_M5000_Merged.root"]            = 167000;
      mcEvtDBS["ExLepEEZ_M750_Merged.root"]             = 192800;
    }
   
    std::map<std::string, double>::iterator itr_mc;
    if(isMC){ 
      cout<<"=============================="<<endl;
      cout<<"        ALL MC SAMPLES        "<<endl;
      cout<<"=============================="<<endl;
      for(itr_mc = mcEvtDBS.begin(); itr_mc != mcEvtDBS.end(); ++itr_mc){
	//cout<<itr_mc->first<<endl;
        TString inFile(itr_mc->first);
        TFile* ttbar= new TFile(inFile);
        TString path= "base/totalEvents";
        TH1F* hist= (TH1F*)(ttbar->Get(path));
        int entries= hist->GetBinContent(1);
        double mean= hist->GetMean();
        double event_cond = entries*mean;//events from condor submission
        double event_dbs = itr_mc->second;//events at data base server
        double ratio = event_dbs/event_cond;
        cout<<setw(30)<<inFile<<setw(15)<<event_dbs<<setw(15)<<event_cond<<setw(15)<<ratio<<endl;
      }
    }
  
    std::map<std::string, double> muDataEvtDBS;
    muDataEvtDBS["MuRunB2v2_Merged.root"] =77914283 ;  
    muDataEvtDBS["MuRunCv1_Merged.root"]  =26691984 ;
    muDataEvtDBS["MuRunDv1_Merged.root"]  =33279413 ;
    muDataEvtDBS["MuRunEv1_Merged.root"]  =26454838 ;
    muDataEvtDBS["MuRunFv1_Merged.root"]  =20178544 ;
    muDataEvtDBS["MuRunGv1_Merged.root"]  =44581284 ;
    muDataEvtDBS["MuRunH2v1_Merged.root"] =46463367 ; 
    muDataEvtDBS["MuRunH3v1_Merged.root"] =1218674 ;
    
    std::map<std::string, double>::iterator itr_mudata;
    if(isMuData){ 
      cout<<"=============================="<<endl;
      cout<<"     ALL Muon Data SAMPLES    "<<endl;
      cout<<"=============================="<<endl;
      for(itr_mudata = muDataEvtDBS.begin(); itr_mudata != muDataEvtDBS.end(); ++itr_mudata){
        TString inFile(itr_mudata->first);
        TFile* ttbar= new TFile(inFile);
        TString path= "base/totalEvents";
        TH1F* hist= (TH1F*)(ttbar->Get(path));
        int entries= hist->GetBinContent(1);
        double mean= hist->GetMean();
        double event_cond = entries*mean;//events from condor submission
        double event_dbs = itr_mudata->second;//events at data base server
        double ratio = event_dbs/event_cond;
        cout<<setw(30)<<inFile<<setw(15)<<event_dbs<<setw(15)<<event_cond<<setw(15)<<ratio<<endl;
      }
    }
  
    std::map<std::string, double> eleDataEvtDBS;
    eleDataEvtDBS["EleRunBver2v2_Merged.root"]  =135276391 ;  
    eleDataEvtDBS["EleRunCv1_Merged.root"    ]  =45523434  ;
    eleDataEvtDBS["EleRunDv1_Merged.root"    ]  =52423951  ;
    eleDataEvtDBS["EleRunEv1_Merged.root"    ]  =47328841  ;
    eleDataEvtDBS["EleRunFv1_Merged.root"    ]  =33972078  ;
    eleDataEvtDBS["EleRunGv1_Merged.root"    ]  =76538905  ;
    eleDataEvtDBS["EleRunHver2v1_Merged.root"]  =79557798  ; 

    std::map<std::string, double>::iterator itr_eledata;
    if(isEleData){ 
      cout<<"=============================="<<endl;
      cout<<" ALL Electron Data SAMPLES    "<<endl;
      cout<<"=============================="<<endl;
      for(itr_eledata = eleDataEvtDBS.begin(); itr_eledata != eleDataEvtDBS.end(); ++itr_eledata){
        TString inFile(itr_eledata->first);
        TFile* ttbar= new TFile(inFile);
        TString path= "base/totalEvents";
        TH1F* hist= (TH1F*)(ttbar->Get(path));
        int entries= hist->GetBinContent(1);
        double mean= hist->GetMean();
        double event_cond = entries*mean;//events from condor submission
        double event_dbs = itr_eledata->second;//events at data base server
        double ratio = event_dbs/event_cond;
        cout<<setw(30)<<inFile<<setw(15)<<event_dbs<<setw(15)<<event_cond<<setw(15)<<ratio<<endl;
      }
    }
  }//for2016Data

  if(is2017Data){
    bool isMC = true;
    bool isMuMC = true;
    bool isMuData = true;
    bool isEleMC = false;
    bool isEleData = false;

    std::map<std::string, double> mcEvtDBS;
    //mcEvtDBS["DY1JetsToLL_Merged.root"]       =  62079400;
    //mcEvtDBS["DY2JetsToLL_Merged.root"]       =  19970551;
    //mcEvtDBS["DY3JetsToLL_Merged.root"]       =  5856110;
    //mcEvtDBS["DY4JetsToLL_Merged.root"]       =  4197868;
    //mcEvtDBS["DYJetsToLLamcatnlo_Merged.root"]=  115626000;
    mcEvtDBS["DYJetsToLL_M50_Merged.root"]      =  182359906;
    mcEvtDBS["TT_Merged.root"]                  =  153531390;   
    //mcEvtDBS["W1JetsToLNu_Merged.root"]       =  44813600;
    //mcEvtDBS["W2JetsToLNu_Merged.root"]       =  29878415;
    //mcEvtDBS["W3JetsToLNu_Merged.root"]       =  19798117;
    //mcEvtDBS["W4JetsToLNu_Merged.root"]       =  9170576;
    mcEvtDBS["WJetsToLNu_Merged.root"]        =  44767978;
    mcEvtDBS["WW_Merged.root"]                =  7791498;
    mcEvtDBS["WZ_Merged.root"]                =  3928630;
    mcEvtDBS["ZZ_Merged.root"]                =  1949768; 

    if(isMuMC){
      mcEvtDBS["ExLepMuMuZ_M250_Merged.root"]    =  188000;
      mcEvtDBS["ExLepMuMuZ_M500_Merged.root"]    =  200000;
      mcEvtDBS["ExLepMuMuZ_M750_Merged.root"]    =  196000;
      mcEvtDBS["ExLepMuMuZ_M1000_Merged.root"]   =  200000;
      mcEvtDBS["ExLepMuMuZ_M1250_Merged.root"]   =  200000;
      mcEvtDBS["ExLepMuMuZ_M1500_Merged.root"]   =  200000;
      mcEvtDBS["ExLepMuMuZ_M1750_Merged.root"]   =  188000;
      mcEvtDBS["ExLepMuMuZ_M2000_Merged.root"]   =  200000;
      mcEvtDBS["ExLepMuMuZ_M2500_Merged.root"]   =  200000;
      mcEvtDBS["ExLepMuMuZ_M3000_Merged.root"]   =  190000;
      mcEvtDBS["ExLepMuMuZ_M3500_Merged.root"]   =  197000;
      mcEvtDBS["ExLepMuMuZ_M4000_Merged.root"]   =  194000;
      mcEvtDBS["ExLepMuMuZ_M4500_Merged.root"]   =  190000;
      mcEvtDBS["ExLepMuMuZ_M5000_Merged.root"]   =  194000;
    }
    if(isEleMC){
      mcEvtDBS["ExLepEEZ_M1000_Merged.root"]     = 192000;
      mcEvtDBS["ExLepEEZ_M1250_Merged.root"]     = 194000;
      mcEvtDBS["ExLepEEZ_M1500_Merged.root"]     = 198000;
      mcEvtDBS["ExLepEEZ_M1750_Merged.root"]     = 194000;
      mcEvtDBS["ExLepEEZ_M2000_Merged.root"]     = 186000;
      mcEvtDBS["ExLepEEZ_M250_Merged.root"]      = 198000;
      mcEvtDBS["ExLepEEZ_M2500_Merged.root"]     = 187000;
      mcEvtDBS["ExLepEEZ_M3000_Merged.root"]     = 200000;
      mcEvtDBS["ExLepEEZ_M3500_Merged.root"]     = 188000;
      mcEvtDBS["ExLepEEZ_M4000_Merged.root"]     = 200000;
      mcEvtDBS["ExLepEEZ_M4500_Merged.root"]     = 200000;
      mcEvtDBS["ExLepEEZ_M500_Merged.root"]      = 200000;
      mcEvtDBS["ExLepEEZ_M5000_Merged.root"]     = 200000;
      mcEvtDBS["ExLepEEZ_M750_Merged.root"]      = 189000;
    }
    
    std::map<std::string, double>::iterator itr_mc;
    if(isMC){ 
      cout<<"=============================="<<endl;
      cout<<"        ALL MC SAMPLES        "<<endl;
      cout<<"=============================="<<endl;
      for(itr_mc = mcEvtDBS.begin(); itr_mc != mcEvtDBS.end(); ++itr_mc){
        TString inFile(itr_mc->first);
        TFile* ttbar= new TFile(inFile);
        TString path= "base/totalEvents";
        TH1F* hist= (TH1F*)(ttbar->Get(path));
        int entries= hist->GetBinContent(1);
        double mean= hist->GetMean();
        double event_cond = entries*mean;//events from condor submission
        double event_dbs = itr_mc->second;//events at data base server
        double ratio = event_dbs/event_cond;
        cout<<setw(30)<<inFile<<setw(15)<<event_dbs<<setw(15)<<event_cond<<setw(15)<<ratio<<endl;
      }
    }

    std::map<std::string, double> muDataEvtDBS;
    muDataEvtDBS["MuRunB2v2_Merged.root"] = 13572040;  
    muDataEvtDBS["MuRunCv1_Merged.root"] = 46273101;
    muDataEvtDBS["MuRunDv1_Merged.root"] = 22646081;
    muDataEvtDBS["MuRunEv1_Merged.root"] = 50596536;
    muDataEvtDBS["MuRunFv1_Merged.root"]  =75301738;
    muDataEvtDBS["MuRunGv1_Merged.root"]  =44581284 ;
    muDataEvtDBS["MuRunH2v1_Merged.root"] =46463367 ; 
    muDataEvtDBS["MuRunH3v1_Merged.root"] =1218674 ;
    
    std::map<std::string, double>::iterator itr_mudata;
    if(isMuData){ 
      cout<<"=============================="<<endl;
      cout<<"     ALL Muon Data SAMPLES    "<<endl;
      cout<<"=============================="<<endl;
      for(itr_mudata = muDataEvtDBS.begin(); itr_mudata != muDataEvtDBS.end(); ++itr_mudata){
        TString inFile(itr_mudata->first);
        TFile* ttbar= new TFile(inFile);
        TString path= "base/totalEvents";
        TH1F* hist= (TH1F*)(ttbar->Get(path));
        int entries= hist->GetBinContent(1);
        double mean= hist->GetMean();
        double event_cond = entries*mean;//events from condor submission
        double event_dbs = itr_mudata->second;//events at data base server
        double ratio = event_dbs/event_cond;
        cout<<setw(30)<<inFile<<setw(15)<<event_dbs<<setw(15)<<event_cond<<setw(15)<<ratio<<endl;
      }
    }

    std::map<std::string, double> eleDataEvtDBS;
    eleDataEvtDBS["EleRunBver2v2_Merged.root"]  =134354864 ;  
    eleDataEvtDBS["EleRunCv1_Merged.root"    ]  =45523434  ;
    eleDataEvtDBS["EleRunDv1_Merged.root"    ]  =52423951  ;
    eleDataEvtDBS["EleRunEv1_Merged.root"    ]  =47328841  ;
    eleDataEvtDBS["EleRunFv1_Merged.root"    ]  =33972078  ;
    eleDataEvtDBS["EleRunGv1_Merged.root"    ]  =76538905  ;
    eleDataEvtDBS["EleRunHver2v1_Merged.root"]  =79557798  ; 
    
    std::map<std::string, double>::iterator itr_eledata;
    if(isEleData){ 
      cout<<"=============================="<<endl;
      cout<<" ALL Electron Data SAMPLES    "<<endl;
      cout<<"=============================="<<endl;
      for(itr_eledata = eleDataEvtDBS.begin(); itr_eledata != eleDataEvtDBS.end(); ++itr_eledata){
        TString inFile(itr_eledata->first);
        TFile* ttbar= new TFile(inFile);
        TString path= "base/totalEvents";
        TH1F* hist= (TH1F*)(ttbar->Get(path));
        int entries= hist->GetBinContent(1);
        double mean= hist->GetMean();
        double event_cond = entries*mean;//events from condor submission
        double event_dbs = itr_eledata->second;//events at data base server
        double ratio = event_dbs/event_cond;
        cout<<setw(30)<<inFile<<setw(15)<<event_dbs<<setw(15)<<event_cond<<setw(15)<<ratio<<endl;
      }
    }
  }//for2017Data
}

