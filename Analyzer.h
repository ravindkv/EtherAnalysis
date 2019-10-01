
#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TTimeStamp.h"
#include "Math/VectorUtil.h"

#include "interface/Reader.h"
#include "interface/ObjectSelector.hh"
#include "interface/MomentumVec.h"
#include "interface/LumiReweighting.h"
#include "interface/UncertaintyComputer.hh"
#include "interface/HistogramPlotter.hh"
#include "interface/BTagCalibrationStandalone.h"


//---------------------------------------------------//
// Main Class
//---------------------------------------------------//
class Analyzer : public ObjectSelector, HistogramPlotter
{
public :
  Analyzer() : ObjectSelector(), HistogramPlotter()
  {
    //---------------------------------------------------//
    //Pileup reweigting 
    //---------------------------------------------------//
    //PU info for Data: 
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
    //PU info for MC:
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pileup_MC_Information
    LumiWeights_ = reweight::LumiReWeighting("data/lumiRewgt/trueInTimePU_mcDY.root","data/lumiRewgt/trueMinBiasPU_dataMu.root", "pileup", "pileup");
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    //---------------------------------------------------//
    //MC cross sections at 13 TeV 
    //---------------------------------------------------//
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV.py
    //https://github.com/BristolTopGroup/AnalysisSoftware/blob/master/python/DataSetInfo_13TeV_25ns.py
    //https://indico.cern.ch/event/617002/contributions/2490586/attachments/1419016/2173704/update_27022017.pdf
    //evtDBS= event at Data Base Server i.e in DAS (https://cmsweb.cern.ch/das/).
    //Backgrounds
    //negative_weight * nEvents
    double evt_M50   = 0.6720* 122055388; 
    double evt_Pt50  = 0.3710* 21890432;
    double evt_Pt100 = 0.3611* 2040596;
    double evt_Pt250 = 0.3717* 423976;
    double evt_Pt400 = 0.3878* 432056;
    double evt_Pt650 = 0.4168* 430691;
    
    double xss_M50   = 5758.4 ; 
    double xss_Pt50  = 354.8  ;
    double xss_Pt100 = 81.22  ;
    double xss_Pt250 = 2.991  ;
    double xss_Pt400 = 0.3882 ;
    double xss_Pt650 = 0.03737;

    double lumiDY_M50  	= evt_M50/xss_M50;
    double lumiDY_Pt50  = evt_Pt50/xss_Pt50 ;
    double lumiDY_Pt100 = evt_Pt100/xss_Pt100;
    double lumiDY_Pt250 = evt_Pt250/xss_Pt250;
    double lumiDY_Pt400 = evt_Pt400/xss_Pt400;
    double lumiDY_Pt650 = evt_Pt650/xss_Pt650;

    double normDY_Pt50  = 1/(lumiDY_M50 + lumiDY_Pt50 );
    double normDY_Pt100 = 1/(lumiDY_M50 + lumiDY_Pt100);
    double normDY_Pt250 = 1/(lumiDY_M50 + lumiDY_Pt250);
    double normDY_Pt400 = 1/(lumiDY_M50 + lumiDY_Pt400);
    double normDY_Pt650 = 1/(lumiDY_M50 + lumiDY_Pt650);
    xss["DYJetsToLL_M50"]            = xss_M50;      		
    xss["DYJetsToLL_Pt50To100"]      = normDY_Pt50 * xss_Pt50 ;    
    xss["DYJetsToLL_Pt100To250"]     = normDY_Pt100* xss_Pt100;    
    xss["DYJetsToLL_Pt250To400"]     = normDY_Pt250* xss_Pt250;    
    xss["DYJetsToLL_Pt400To650"]     = normDY_Pt400* xss_Pt400;    
    xss["DYJetsToLL_Pt650ToInf"]     = normDY_Pt650* xss_Pt650;    
    evtDBS["DYJetsToLL_Pt650ToInf"]  = evt_Pt650;
    evtDBS["DYJetsToLL_M50"]         = evt_M50  ; 
    evtDBS["DYJetsToLL_Pt50To100"]   = evt_Pt50 ;
    evtDBS["DYJetsToLL_Pt100To250"]  = evt_Pt100; 
    evtDBS["DYJetsToLL_Pt250To400"]  = evt_Pt250;
    evtDBS["DYJetsToLL_Pt400To650"]  = evt_Pt400;

    xss["TT"] 			= 831.76;         
    xss["WJetsToLNu"]	= 50690;          
    xss["WW"]			= 118.7;          
    xss["WZ"]			= 46.74;	    
    xss["ZZ"] 			= 17.72;          
    evtDBS["TT"] 		= 77081156;
    evtDBS["WJetsToLNu"]= 29705748;
    evtDBS["WW"]		= 994012;
    evtDBS["WZ"]		= 1000000;
    evtDBS["ZZ"] 		= 990064;

    //Muon signal
    xss["ExLepMuMuZ_M250"] = 0.00427    ;   evtDBS["ExLepMuMuZ_M250"] = 198200;
    xss["ExLepMuMuZ_M500"] = 0.00291    ;   evtDBS["ExLepMuMuZ_M500"] = 200000;
    xss["ExLepMuMuZ_M750"] = 0.001761   ;   evtDBS["ExLepMuMuZ_M750"] = 191400;
    xss["ExLepMuMuZ_M1000"]= 0.001177   ;   evtDBS["ExLepMuMuZ_M1000"]= 191300;
    xss["ExLepMuMuZ_M1250"]= 0.0007263  ;   evtDBS["ExLepMuMuZ_M1250"]= 181400;
    xss["ExLepMuMuZ_M1500"]= 0.0004267  ;   evtDBS["ExLepMuMuZ_M1500"]= 190900;
    xss["ExLepMuMuZ_M1750"]= 0.0002654  ;   evtDBS["ExLepMuMuZ_M1750"]= 192200;
    xss["ExLepMuMuZ_M2000"]= 0.0002021  ;   evtDBS["ExLepMuMuZ_M2000"]= 200000;
    xss["ExLepMuMuZ_M2500"]= 0.00006755 ;   evtDBS["ExLepMuMuZ_M2500"]= 200000;
    xss["ExLepMuMuZ_M3000"]= 0.00002223 ;   evtDBS["ExLepMuMuZ_M3000"]= 200000;
    xss["ExLepMuMuZ_M3500"]= 0.000008166;   evtDBS["ExLepMuMuZ_M3500"]= 200000;
    xss["ExLepMuMuZ_M4000"]= 0.000003209;   evtDBS["ExLepMuMuZ_M4000"]= 192300;
    xss["ExLepMuMuZ_M4500"]= 0.000001133;   evtDBS["ExLepMuMuZ_M4500"]= 200000;
    xss["ExLepMuMuZ_M5000"]= 0.000000452;   evtDBS["ExLepMuMuZ_M5000"]= 192500;

    //Electron signal
    xss["ExLepEEZ_M250"]   = 0.004407   ;    evtDBS["ExLepEEZ_M250"]  = 182100;
    xss["ExLepEEZ_M500"]   = 0.002916   ;    evtDBS["ExLepEEZ_M500"]  = 200000;
    xss["ExLepEEZ_M750"]   = 0.001607   ;    evtDBS["ExLepEEZ_M750"]  = 192800;
    xss["ExLepEEZ_M1000"]  = 0.001079   ;    evtDBS["ExLepEEZ_M1000"] = 197600;
    xss["ExLepEEZ_M1250"]  = 0.0007171  ;    evtDBS["ExLepEEZ_M1250"] = 200000;
    xss["ExLepEEZ_M1500"]  = 0.0004856  ;    evtDBS["ExLepEEZ_M1500"] = 41900 ;
    xss["ExLepEEZ_M1750"]  = 0.0002907  ;    evtDBS["ExLepEEZ_M1750"] = 192400;
    xss["ExLepEEZ_M2000"]  = 0.0001908  ;    evtDBS["ExLepEEZ_M2000"] = 189000;
    xss["ExLepEEZ_M2500"]  = 0.00005731 ;    evtDBS["ExLepEEZ_M2500"] = 187400;
    xss["ExLepEEZ_M3000"]  = 0.00002497 ;    evtDBS["ExLepEEZ_M3000"] = 192800;
    xss["ExLepEEZ_M3500"]  = 0.000009023;    evtDBS["ExLepEEZ_M3500"] = 170000;
    xss["ExLepEEZ_M4000"]  = 0.000003897;    evtDBS["ExLepEEZ_M4000"] = 199800;
    xss["ExLepEEZ_M4500"]  = 0.000001243;    evtDBS["ExLepEEZ_M4500"] = 187400;
    xss["ExLepEEZ_M5000"]  = 0.000000461;    evtDBS["ExLepEEZ_M5000"] = 167000;

    //Dummy sample
    xss["sampCode_"]       =  1;  	         evtDBS["sampCode_"]      =  1; 

  };
  ~Analyzer() {
    delete evR;
  };
  
  void CutFlowAnalysis(TString url, bool isMuChannel, bool isEleChannel, TFile *outFile_);
  void CutFlowProcessor(TString url,  TString cutflowType="base", TFile *outFile_=0);
  //void CreateAnalHistos(TString flowType, TFile* outFile_);
  void processEvents();
private :
  double DRMIN_JET, DRMIN_ELE, METCUT_;
  Reader *evR;
  
  reweight::LumiReWeighting LumiWeights_;
  reweight::PoissonMeanShifter PShiftUp_;   //pileup syst up
  reweight::PoissonMeanShifter PShiftDown_; //pileup syst down 
  std::map<string, double> xss;
  std::map<string, double> evtDBS;
  std::map<string, double> muSF;
  std::map<string, double> eleSF;
  
  BTagCalibrationReader readCSV(const std::string &filename,const std::string &tagger, BTagEntry::OperatingPoint op, const std::string & measurementType, const std::string & sysType, const std::vector<std::string> & otherSysTypes, BTagEntry::JetFlavor jf);
  ofstream outfile_;
  Double_t getMuonSF(TH2D *h2, double eta, double pt);
  Double_t getMuonTrigSF(TH2D *h2, double eta, double pt);
  Double_t getMuonTrackSF(TGraphAsymmErrors *tg, double eta);
  Double_t getEleSF(TH2D *h2, double etaSC, double pt);
  Double_t getEleTrigSF(TH2D *h2, double etaSC, double pt);
  Double_t getEleHeep1SF(TGraphAsymmErrors *tg, double pt);
  Double_t getEleHeep2SF(TGraphAsymmErrors *tg, double eta);
  double deltaPhi12(double phi1, double phi2);
  double phi0to2pi(double phi);
};

BTagCalibrationReader Analyzer::readCSV(const std::string &filename, const std::string &tagger,  
		BTagEntry::OperatingPoint op, 
		const std::string & measurementType,
		const std::string & sysType, 
		const std::vector<std::string> & otherSysTypes, 
		BTagEntry::JetFlavor jf
		)
{ 
  BTagCalibration calib(tagger, filename);
  BTagCalibrationReader reader(op, sysType, otherSysTypes);      
  reader.load(calib, jf, measurementType); 
  return reader;
}
//https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
Double_t Analyzer::getMuonSF(TH2D *h2, double eta, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <120
  //for Pt >120, we use SF of Pt = 120
  if(pt<=120){
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(120);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}
Double_t Analyzer::getMuonTrigSF(TH2D *h2, double eta, double pt){
  
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <120
  //for Pt >120, we use SF of Pt = 120
  if(pt<=500){
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(eta));
    Int_t binY = yaxis->FindBin(500);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}

Double_t Analyzer::getMuonTrackSF(TGraphAsymmErrors *tg, double eta){
 
  Double_t *eta_array = tg->GetX();
  Double_t *sf_array = tg->GetY();
  Int_t n_points = tg->GetN();

  double SF = 1.0;
  // eta < eta_array[0]
  if(abs(eta)<eta_array[0]) SF = sf_array[0];
  
  // eta_array[0]<eta<eta_array[n_points -1]
  for(Int_t i = 0; i < n_points-1; i++){
    if(abs(eta) >= eta_array[i] && abs(eta) < eta_array[i+1]) SF = sf_array[i+1];
  }
  // eta > eta_array[n_points -]
  if(abs(eta)>eta_array[n_points-1]) SF = sf_array[n_points -1];
  return SF;
}

Double_t Analyzer::getEleSF(TH2D *h2, double etaSC, double pt){
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <500
  //for Pt >500, we use SF of Pt = 500
  if(pt<=500){
    Int_t binX = xaxis->FindBin(abs(etaSC));
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(abs(etaSC));
    Int_t binY = yaxis->FindBin(500);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}
Double_t Analyzer::getEleHeep1SF(TGraphAsymmErrors *tg, double pt){
  Double_t *pt_array = tg->GetX();
  Double_t *sf_array = tg->GetY();
  Int_t n_points = tg->GetN();
  double SF = 1.0;
if(pt < pt_array[0]) SF = sf_array[0];
for(Int_t i = 0; i < n_points-1; i++){
  if(pt >= pt_array[i] && pt < pt_array[i+1]) SF = sf_array[i+1];
  }
  if(pt > pt_array[n_points-1]) SF = sf_array[n_points-1];
  return SF;
}

Double_t Analyzer::getEleHeep2SF(TGraphAsymmErrors *tg, double etaSC){

  Double_t *eta_array = tg->GetX();
  Double_t *sf_array = tg->GetY();
  Int_t n_points = tg->GetN();

double SF = 1.0;
if(abs(etaSC)<eta_array[0]) SF = sf_array[0];
for(Int_t i = 0; i < n_points-1; i++){
  if(abs(etaSC) >= eta_array[i] && abs(etaSC) < eta_array[i+1]) SF = sf_array[i+1];
  }
  if(abs(etaSC)>eta_array[n_points-1]) SF = sf_array[n_points -1];
return SF;
}

Double_t Analyzer::getEleTrigSF(TH2D *h2, double etaSC, double pt){
  TAxis *xaxis = h2->GetXaxis();
  TAxis *yaxis = h2->GetYaxis();
  //since the Pt range of 2D histo is <500
  //for Pt >500, we use SF of Pt = 500
  if(pt<=500){
    Int_t binX = xaxis->FindBin(etaSC);
    Int_t binY = yaxis->FindBin(pt);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }
  else{
    Int_t binX = xaxis->FindBin(etaSC);
    Int_t binY = yaxis->FindBin(500);
    double sf = h2->GetBinContent(binX, binY);
    double err = h2->GetBinError(binX, binY);
    if(sf!=0) return sf;
    else return 1.0;
  }	  
}

double phi0to2pi(double phi){
    double pi = 3.1415926535;
    while (phi >= 2.*pi) phi -= 2.*pi;
    while (phi < 0.) phi += 2.*pi;
    return phi;
}

double deltaPhi12(double phi1_, double phi2_){
    // build the delta Phi angle between the two vectors
    double pi = 3.1415926535;
    double phi1 = phi0to2pi(phi1_);
    double phi2 = phi0to2pi(phi2_);
    double dPhi = phi0to2pi(phi1 - phi2);
    dPhi = (dPhi > (2*pi - dPhi)) ? 2*pi - dPhi : dPhi;
    return dPhi;
}

//---------------------------------------------------//
//muon scale factors from 2D histograms 
//---------------------------------------------------//      
//https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
//https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2 
//Trigger SF
TFile *f_trigSF_BCDEF 	= new TFile("data/muonSF/triggreSF_BCDEF.root");
TFile *f_trigSF_GH 		= new TFile("data/muonSF/triggreSF_GH.root");
TH2D *h2_trigSF_BCDEF 	= (TH2D*)f_trigSF_BCDEF->Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio");
TH2D *h2_trigSF_GH 		= (TH2D*)f_trigSF_GH->Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio");
//Identification SF
TFile *f_idSF_BCDEF 		= new TFile("data/muonSF/idSF_BCDEF.root");
TFile *f_idSF_GH 		= new TFile("data/muonSF/idSF_GH.root");
TH2D *h2_idSF_BCDEF 		= (TH2D*)f_idSF_BCDEF->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/abseta_pair_ne_ratio");
TH2D *h2_idSF_GH 		= (TH2D*)f_idSF_GH->Get("MC_NUM_HighPtID_DEN_genTracks_PAR_newpt_eta/abseta_pair_ne_ratio");
//Isolation SF
TFile *f_isoSF_BCDEF 		= new TFile("data/muonSF/isoSF_BCDEF.root");
TFile *f_isoSF_GH 		= new TFile("data/muonSF/isoSF_GH.root");
TH2D *h2_isoSF_BCDEF 		= (TH2D*)f_isoSF_BCDEF->Get("tkLooseISO_highptID_newpt_eta/abseta_pair_ne_ratio");
TH2D *h2_isoSF_GH 		= (TH2D*)f_isoSF_GH->Get("tkLooseISO_highptID_newpt_eta/abseta_pair_ne_ratio");

//---------------------------------------------------//
//Electron scale factors from 2D histograms 
//---------------------------------------------------//      
//https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
//Reconstruction SF
TFile *f_ele_recoSF 	  	= new TFile("data/eleSF/ele_recoSF.root");
TH2D *h2_ele_recoSF 		= (TH2D*)f_ele_recoSF->Get("EGamma_SF2D");
//Identification SF
TFile *f_ele_veto_idSF 	= new TFile("data/eleSF/ele_veto_idSF.root");
TH2D *h2_ele_veto_idSF 	= (TH2D*)f_ele_veto_idSF->Get("EGamma_SF2D");
TFile *f_ele_heep_SF_EE     = new TFile("data/eleSF/HEEP_SF.root");
TFile *f_ele_heep_SF_EB     = new TFile("data/eleSF/HEEP_SF.root");
TFile *f_ele_heep_SF            = new TFile("data/eleSF/HEEP_SF.root");
TGraphAsymmErrors *tg_heep_SF_EE      = (TGraphAsymmErrors*)f_ele_heep_SF_EE->Get("SF_Et_Endcap");
TGraphAsymmErrors *tg_heep_SF_EB      = (TGraphAsymmErrors*)f_ele_heep_SF_EB->Get("SF_Et_Barrel");
TGraphAsymmErrors *tg_heep_SF             = (TGraphAsymmErrors*)f_ele_heep_SF->Get("SF_Eta");

//Trigger scale factors
//https://indico.cern.ch/event/604912/
//TFile *f_ele_trigSF 		= new TFile("data/eleSF/ele_trigSF_Run2016All_v1.root");
//TH2D *h2_ele_trigSF 		= (TH2D*)f_ele_trigSF->Get("Ele27_WPTight_Gsf");
TFile *f_ele_trigSF 		= new TFile("data/eleSF/MW_2ndleg_EGM2D.root");
TH2D *h2_ele_trigSF 		= (TH2D*)f_ele_trigSF->Get("EGamma_SF2D");
