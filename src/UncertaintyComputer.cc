#include "interface/UncertaintyComputer.hh"
#include <iostream>
#include <iomanip>
#include "TRandom3.h"

ClassImp(UncertaintyComputer)

using namespace std;
double UncertaintyComputer::getJERSF(double eta, int jer){
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
  std::vector<double> JEREtaMap;
  std::vector<double> JERSF;
  std::vector<double> JERSFUp;
  std::vector<double> JERSFDown;
  //eta values
  JEREtaMap.push_back(0.000);
  JEREtaMap.push_back(0.522);
  JEREtaMap.push_back(0.783);
  JEREtaMap.push_back(1.131);
  JEREtaMap.push_back(1.305);
  JEREtaMap.push_back(1.740);
  JEREtaMap.push_back(1.930);
  JEREtaMap.push_back(2.043);
  JEREtaMap.push_back(2.322);
  JEREtaMap.push_back(2.500);
  JEREtaMap.push_back(2.853);
  JEREtaMap.push_back(2.964);
  JEREtaMap.push_back(3.139);
  JEREtaMap.push_back(5.191);
  //SF values
  JERSF.push_back(1.1595); JERSFDown.push_back(1.095);  JERSFUp.push_back(1.224);
  JERSF.push_back(1.1948); JERSFDown.push_back(1.1296); JERSFUp.push_back(1.26);
  JERSF.push_back(1.1464); JERSFDown.push_back(1.0832); JERSFUp.push_back(1.2096);
  JERSF.push_back(1.1609); JERSFDown.push_back(1.0584); JERSFUp.push_back(1.2634);
  JERSF.push_back(1.1278); JERSFDown.push_back(1.0292); JERSFUp.push_back(1.2264);
  JERSF.push_back(1.1000); JERSFDown.push_back(0.9921); JERSFUp.push_back(1.2079);
  JERSF.push_back(1.1426); JERSFDown.push_back(1.0212); JERSFUp.push_back(1.264);
  JERSF.push_back(1.1512); JERSFDown.push_back(1.0372); JERSFUp.push_back(1.2652);
  JERSF.push_back(1.2963); JERSFDown.push_back(1.0592); JERSFUp.push_back(1.5334);
  JERSF.push_back(1.3418); JERSFDown.push_back(1.1327); JERSFUp.push_back(1.5509);
  JERSF.push_back(1.7788); JERSFDown.push_back(1.578);  JERSFUp.push_back(1.9796);
  JERSF.push_back(1.1869); JERSFDown.push_back(1.0626); JERSFUp.push_back(1.3112);
  JERSF.push_back(1.1922); JERSFDown.push_back(1.0434); JERSFUp.push_back(1.341);
  double SF = 1.0;
  for(size_t i = 0; i < 13; i++){
    if(TMath::Abs(eta) >= JEREtaMap[i] && TMath::Abs(eta) < JEREtaMap[i+1]){
      if(jer == 0)SF = JERSF[i];
      else if (jer == 1) SF = JERSFUp[i];
      else if(jer == -1) SF = JERSFDown[i];  
    }
  }
  return SF;
}

double UncertaintyComputer::jetPtWithJESJER(MyJet jet, int jes, int jer){
  double gen_pt = jet.Genp4.pt();  
  double jet_pt = jet.p4.pt();
  double sigmaJER = jet.resolution ;
  //apply JES uncert scaling 
  jet_pt *= (1+(jet.JECUncertainty*double(jes)));
  //apply JER uncert, scaling
  double delR = DeltaR(jet.Genp4, jet.p4);
  double rCone = 0.4;
  if(gen_pt> 0 && delR<rCone/2 && abs(jet_pt -gen_pt)<3*sigmaJER*jet_pt ){
  //if(gen_pt > 0){
    double SF = getJERSF(jet.p4.eta(), jer);
    double ptscale = max(0.0, 1.0 + (SF - 1)*(jet_pt - gen_pt)/ jet_pt);
    jet_pt *= ptscale;
  }
  return jet_pt;
}

//bottom mistagging, by event re-weighting 
double UncertaintyComputer::getBTagPmcSys(TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet){
  double csv =jet.bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
  double pMC = 1.0; 
  pMC = btsf->getBTagPmc(h2_qTagEff_Num, h2_qTagEff_Denom, jet.p4.eta(), jet.p4.pt(), csv);
  return pMC;
}
double UncertaintyComputer::getBTagPdataSys(BTagCalibrationReader &reader, TH2D *h2_qTagEff_Num, TH2D *h2_qTagEff_Denom, MyJet jet, int scale){
  double pData = 1.0;
  double csv =jet.bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
  double eta = jet.p4.eta();
  double pt = jet.p4.pt();
  int flavor = abs(jet.partonFlavour);
  if(scale == 0) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,0);
  else if(scale == 1) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,1);
  else if(scale == -1) pData = btsf->getBTagPdata(reader, h2_qTagEff_Num, h2_qTagEff_Denom, eta, pt, csv, flavor ,-1);
  return pData;
}

double UncertaintyComputer::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  double delR = sqrt(deta*deta + dphi*dphi);
  return delR;
}
