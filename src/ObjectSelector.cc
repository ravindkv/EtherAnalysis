#include "interface/ObjectSelector.hh"
#include <iostream>
#include <iomanip>
ClassImp(ObjectSelector)

using namespace std;
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
//Electron ID: veto
bool ObjectSelector::cutBasedElectronID_Summer16_80X_V1_veto(const MyElectron *e)
{
  bool passID = false;
  //barrel
  if(abs(e->eleSCEta) <= 1.479 
     && e->sigmaIetaIeta    < 0.0115     
     && abs(e->dEtaInSeed)  < 0.00749   
     && abs(e->dPhiIn)      < 0.228 
     && e->hadOverEm        < 0.356 
     && e->relCombPFIsoEA   < 0.175 
     && abs(e->iEminusiP)   < 0.299 
     && e->nInnerHits           <= 2
     && e->passConversionVeto  
    )passID = true;
  //endcap 
  if(abs(e->eleSCEta) > 1.479 
     && e->sigmaIetaIeta    < 0.037 
     && abs(e->dEtaInSeed)  < 0.00895    
     && abs(e->dPhiIn)      < 0.213 
     && e->hadOverEm        < 0.211 
     && e->relCombPFIsoEA   < 0.159 
     && abs(e->iEminusiP)   < 0.15  
     && e->nInnerHits           <= 3
     && e->passConversionVeto  
     )passID = true;
  return passID;
}
//Electron HEEP ID
  bool ObjectSelector::heepElectronID_HEEPV70(const MyElectron *e, MyVertex & vertex)
{
  bool passID = false;
  float energy2x5Overenergy5x5 = e->energy2x5/e->energy5x5;
  //for barrel
  if(abs(e->eleSCEta)                <=1.4442
     && e->isEcalDriven              == 1 
     && abs(e->dEtaInSeed)           < 0.004 
     && abs(e->dPhiIn)               < 0.06 
     && e->hadOverEm                 < 1/e->p4.E() + 0.05 
     && energy2x5Overenergy5x5       > 0.94 
     && e->GsfEleEmHadD1IsoRhoCut    < 2+0.03*e->p4.pt()+0.28*e->eleRho 
     && e->eleTrkPt                  < 5  
     && e->nInnerHits                <=1   //Inner Lost Hits
     && abs(e->D0)                   < 0.02
     )passID = true;
  //endcap
  double HadDepth = 0.0;
  if(e->p4.E() < 50) HadDepth = 2.5 +0.28*e->eleRho;
  else HadDepth = 2.5+0.03*(e->p4.E()-50) +0.28*e->eleRho; 
  if(abs(e->eleSCEta) > 1.566 && abs(e->eleSCEta) < 2.5
     && e->isEcalDriven                == 1
     && abs(e->dEtaInSeed)             < 0.006
     && abs(e->dPhiIn)                 < 0.06 
     && e->hadOverEm                   < 5/e->p4.E() + 0.05 
     && e->sigmaIetaIeta               <0.03 
     && e->GsfEleEmHadD1IsoRhoCut      < HadDepth
     && e->eleTrkPt                    < 5
     && e->nInnerHits                  <=1
     && abs(e->D0)                     < 0.05
     )passID = true;
  return passID; 
}

void ObjectSelector::preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex){
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
  for(unsigned int i=0;i<vE.size();i++){
    const MyElectron * e   = &vE[i];
    double ePt     	   = TMath::Abs(e->p4.pt());
    double eEta     	   = TMath::Abs(e->p4.eta());
    bool passID = heepElectronID_HEEPV70(e, vertex); 
    if(passID && ePt >35 && eEta <2.5 ){e_i->push_back(i);}
  }
}

//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
//Jonas talk:
//https://indico.cern.ch/event/848326/contributions/3564004/attachments/1910537/3156831/talk.pdf
bool ObjectSelector::isHighPtMuon(const MyMuon * m){
  bool isHighPt(false);
  bool globORTuneHits = false;
  if(m->nMuonHits > 0 || m->nMuonHitsTuneP >0) globORTuneHits = true;
  bool for10_4_X = false;
  if(m->nMatchedStations==1){
    if(m->isTrackerMuon && (m->expectedMatchedStations <2 || 
             (m->nStationMask!=1 && m->nStationMask!=16) || m->nRPCLayers >2)){
       for10_4_X = true; 
    }
  }
  bool isMatchedStations = false;
  if(m->nMatchedStations >1 || for10_4_X) isMatchedStations = true;
  isHighPt = m->isGlobalMuon &&
	      globORTuneHits && isMatchedStations &&
          m->bestMuPtErr/m->bestMuPtTrack < 0.3 &&
          m->nPixelHits > 0 &&
          m->nTrackerLayers > 5;
return isHighPt;
}
void ObjectSelector::preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isData){
  for( int i=0;i< (int) vM.size();i++){
    const MyMuon * m = &vM[i];
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt = TMath::Abs(m->p4.pt());
    bool passID = isHighPtMuon(m);
    double iso = m->trkRelIso;
    if(passID && mPt > 35 && mEta < 2.4 && 
		    abs(m->D0) < 0.2 && abs(m->Dz) < 0.5 && iso < 0.10){ 
      m_i->push_back(i);
    }
  }
}

void ObjectSelector::preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, int jes, int jer){
  for(unsigned int i=0;i<vJ.size();i++){
    const MyJet *jet = &vJ[i];
    double jetEta     = TMath::Abs(jet->p4.eta());
    double jetPt      = jetPtWithJESJER(vJ[i], jes, jer); 
    double neutralHadEnFrac = jet->neutralHadronEnergyFraction;
    double neutralEmEnFrac = jet->neutralEmEnergyFraction;
    double chargedHadEnFrac = jet->chargedHadronEnergyFraction;
    double chargedEmFrac = jet->chargedEmEnergyFraction;
    if(jetEta < 2.4){
      if(jetPt > 170 
        && neutralHadEnFrac < 0.90
        && neutralEmEnFrac  < 0.90
	&& jet->NumConst > 1
        && chargedHadEnFrac > 0
	&& chargedEmFrac < 0.99
	&& jet->chargedMultiplicity > 0){
        j_i->push_back(i);
      }
    }
    if(jetEta >= 2.4 && jetEta < 2.5){
      if(jetPt > 170 
        && neutralHadEnFrac < 0.90
        && neutralEmEnFrac  < 0.90
	&& jet->NumConst > 1){
        j_i->push_back(i);
      }
    }
  }
}


bool ObjectSelector::looseMuonVeto( int selMu1, int selMu2, const vector<MyMuon> & vM){
  bool looseVeto(false);
  for(int i=0;i< (int)vM.size();i++){
    const MyMuon * m = &vM[i];
    bool isGlobalMuon = m->isGlobalMuon;
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mRelIso  = m->pfRelIso;
    if(! isGlobalMuon) continue;
    if( i!=selMu1 && i!=selMu2 && mEta<2.4  && mPt> 15.0 && mRelIso < 0.25){ looseVeto = true; }
  }
  return looseVeto;
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
bool ObjectSelector::looseElectronVeto(unsigned long selEle1,unsigned long selEle2, const vector<MyElectron> & vE){
  bool looseVeto(false);
  for(unsigned long i=0;i<vE.size();i++){
    const MyElectron * e   = &vE[i];
    double ePt         = TMath::Abs(e->p4.pt());
    bool passID = cutBasedElectronID_Summer16_80X_V1_veto(e);
    double dxy = abs(e->D0);
    double dz  = abs(e->Dz);
    if(i!=selEle1 && i!=selEle2 && passID && ePt >15.0 
    && dxy<0.05 && dz < 0.1){looseVeto = true;}
  }
  return looseVeto;
}

double ObjectSelector::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  double delR = sqrt(deta*deta + dphi*dphi);
  return delR;
}

