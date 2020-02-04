
///////////////////////
// Muon Channel
///////////////////////

#include "Analyzer.h"
#include <map>

using namespace std;
void Analyzer::CutFlowAnalysis(TString url, bool isMuChannel, bool isEleChannel, TFile *outFile_){
  //check if the input file is MC or Data  
  Reader *evR_;  
  evR_ = new Reader();
  TFile *f_ = TFile::Open(url);
  int nEntries = evR_->AssignEventTreeFrom(f_);
  MyEvent *ev_;
  ev_ = evR_->GetNewEvent(1);
  TString chName = "Lepton";
  if(isMuChannel) chName = "Muon";
  if(isEleChannel) chName = "Electron";

  CutFlowProcessor(url,  chName+"/base", outFile_);
  /*
  //---------------------------------------------------//
  //for systematics (all sys in one go)
  //---------------------------------------------------//  
  if(!ev_->isData){ 
    CutFlowProcessor(url,  chName+"/JESPlus",       outFile_);
    CutFlowProcessor(url,  chName+"/JESMinus",      outFile_);
    CutFlowProcessor(url,  chName+"/JERPlus",       outFile_);
    CutFlowProcessor(url,  chName+"/JERMinus",      outFile_);
    CutFlowProcessor(url,  chName+"/TopPtPlus", 	outFile_);
    CutFlowProcessor(url,  chName+"/TopPtMinus", 	outFile_);
  }
  */
  f_->Close();
  delete f_;
}

//---------------------------------------------------//
//Process the cuts, event by event
//---------------------------------------------------//  
void Analyzer::CutFlowProcessor(TString url, TString cutflowType, TFile *outFile_){
  int input_count_PreSel = 0;
  int input_count_ZTag = 0;
  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  bool isMuChannel = false;
  bool isEleChannel = false;
  if(cutflowType.Contains("Muon")) isMuChannel = true;
  if(cutflowType.Contains("Electron")) isEleChannel = true;

  //Uncertainty variations, JES, JER, MET unclustered, bTag
  int jes = 0, jer = 0, metuc = 0, bScale = 0;
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }
  
  //---------------------------------------------------//
  //get initial number of events, from ntuples
  //store initial informations, in a txt file
  //---------------------------------------------------//
  double lumiTotal = 35860;
  int nEntries = evR->AssignEventTreeFrom(f);
  if(nEntries == 0) {return; }
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"));
  double initialEvents = inputcf->GetBinContent(1);
  cout<<"\033[01;32m input file: \033[00m"<<url<<"\n"<<endl;
  fillHisto(outFile_, cutflowType, "", "totalEvents", 10, 0, 10000000000, initialEvents, 1 );
  MyEvent *ev;
 
  //---------------------------------------------------//
  //loop over each event, of the ntuple
  //---------------------------------------------------//
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
    //if(i > 10000) break; 
    //---------------------------------------------------//
    //apply lumi, k factor and pileup weight
    //---------------------------------------------------//
    double evtWeight = 1.0;
    double genWeight = 0.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      TString sampleName_(sampleName);
      if(sampleName_.Contains("DYJetsToLL")){
        double sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
        evtWeight *= sampleWeight;
        evtWeight *= ev->sampleInfo.gen_weight;
	    genWeight = ev->sampleInfo.gen_weight;
        fillHisto(outFile_, cutflowType, "", "lumiWeight", 10, 0, 1000, sampleWeight, 1 );
      }	  
      //lumi weight
      else {
        double sampleWeight(1.0);
        sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
        evtWeight *= sampleWeight; 
        fillHisto(outFile_, cutflowType, "", "lumiWeight", 10, 0, 1000, sampleWeight, 1 );
      } 
      //pileup weight
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
        float npu = pu[0];
        double weightPU = LumiWeights_.weight(npu);
        evtWeight *= weightPU;  
        fillHisto(outFile_, cutflowType, "", "puWeight", 1000, 0, 100, weightPU, 1 );
      }
      if(i==0){
        double sampleWeight(1.0);
        sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
	fillHisto(outFile_, cutflowType, "", "totalYield", 10, 0, 2, 1, initialEvents*sampleWeight);
      }
    } 
    else{ 
      if(i==0)fillHisto(outFile_, cutflowType, "", "totalYield", 10, 0, 2, 1, initialEvents);
    }
    double topPtWt = 1.0;
    if(!ev->isData){
      string sampleName = ev->sampleInfo.sampleName;
      if(sampleName.find("TT") != string::npos){
        vector<double>topptweights = ev->sampleInfo.topPtWeights;
        if(topptweights.size() > 0){
          topPtWt = topptweights[0]; 
          if(cutflowType.Contains("TopPtPlus")){
            topPtWt = topptweights[0];
            topPtWt = topPtWt*topPtWt;
          }
          else if(cutflowType.Contains("TopPtMinus"))
            topPtWt = 1.0;
        }
      }
    }
    fillHisto(outFile_, cutflowType, "", "SF_topPtWeights", 1000, 0, 3, topPtWt, 1 );
    evtWeight *= topPtWt; //Multiply to the total weights
        
    //---------------------------------------------------//
    //apply muon triggers
    //---------------------------------------------------//
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5
    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6
    bool passTrig = false;
    bool passOneTrig = false;
    vector<string> trig = ev->hlt;
    //Muon channel
    if(isMuChannel){
      for(size_t it = 0; it < trig.size(); it++){
        if(trig[it].find("HLT_Mu50") != string::npos || 
      	      trig[it].find("HLT_TkMu50") != string::npos) {
          passTrig = true;
        }
      } 
      for(size_t it = 0; it < trig.size(); it++){
        if(trig[it].find("HLT_Mu50") != string::npos) passOneTrig = true;
      } 
    }
    //Electron channel
    if(isEleChannel){
      for(size_t it = 0; it < trig.size(); it++){
        if(trig[it].find("HLT_DoubleEle33_CaloIdL_MW") != string::npos) {
          passTrig = true;
        }
      }
    }
    if(!passTrig) continue;
    double nCutPass = 1.0;
    fillHisto(outFile_, cutflowType, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );
   
    //---------------------------------------------------//
    //get all objets e.g. leptons, jets, vertices etc.
    //---------------------------------------------------//
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    if(Vertices.size() <= 0){
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);

    //preselect objects 
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], ev->isData);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0]);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    
    //---------------------------------------------------//
    //apply selection cuts on leptons
    //---------------------------------------------------//
    int nLepton = 0;
    if(isMuChannel) nLepton = m_init.size();
    if(isEleChannel) nLepton = e_init.size();
    if(nLepton < 2)continue;
    double pri_vtxs = Vertices[0].totVtx;
    int charge1 = 0;
    int charge2 = 0;
    int lepton1 = 0;
    int lepton2 = 0;
    if(isMuChannel){
      lepton1 = m_init[0];
      lepton2 = m_init[1];
      charge1 = pfMuons[lepton1].charge;
      charge2 = pfMuons[lepton2].charge;
    }
    if(isEleChannel){
      lepton1 = e_init[0];
      lepton2 = e_init[1];
      charge1 = pfElectrons[lepton1].charge;
      charge2 = pfElectrons[lepton2].charge;
    }
    //both lepton should have opposite charge
    if(charge1 == charge2) continue;
    //veto third loose lepton
    bool isVeto = false;
    if(isMuChannel){
        isVeto = looseElectronVeto(-1, -1, pfElectrons) || 
            looseMuonVeto(lepton1, lepton2, pfMuons) ;
    }
    if(isEleChannel){
        isVeto = looseElectronVeto(lepton1, lepton2, pfElectrons) || 
            looseMuonVeto(-1, -1, pfMuons) ;
    }
    if(isVeto) continue;
    if(isMuChannel && pfMuons[lepton1].p4.pt() < 53) continue;

    //---------------------------------------------------//
    //apply lepton SF to eventWeights 
    //---------------------------------------------------//
    double leptonSF = 1.0;
    if(isMuChannel && !ev->isData){
      double lumi_BCDEF = 19711; double lumi_GH = 16138;	
      double lumi = lumi_BCDEF + lumi_GH;
      //get muon scale factor for fist muon
      //trigger 	
      double muSFtrig_BCDEF1 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFtrig_GH1 	= getMuonTrigSF(h2_trigSF_GH, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFtrig1 		= (muSFtrig_BCDEF1*lumi_BCDEF + muSFtrig_GH1*lumi_GH)/lumi; 

      //identification
      double muSFid_BCDEF1 	= getMuonSF(h2_idSF_BCDEF, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFid_GH1 		= getMuonSF(h2_idSF_GH, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFid1 		= (muSFid_BCDEF1*lumi_BCDEF + muSFid_GH1*lumi_GH)/lumi; 
      //isolation 
      double muSFiso_BCDEF1 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFiso_GH1 		= getMuonSF(h2_isoSF_GH, pfMuons[lepton1].p4.eta(), pfMuons[lepton1].p4.pt());
      double muSFiso1 		= (muSFiso_BCDEF1*lumi_BCDEF + muSFiso_GH1*lumi_GH)/lumi; 
      double muSF1 = muSFtrig1*muSFid1*muSFiso1;	
      
      //get muon scale factor for 2nd muon
      //trigger 	
      double muSFtrig_BCDEF2 	= getMuonTrigSF(h2_trigSF_BCDEF, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFtrig_GH2 	= getMuonTrigSF(h2_trigSF_GH, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFtrig2 		= (muSFtrig_BCDEF2*lumi_BCDEF + muSFtrig_GH2*lumi_GH)/lumi; 
      //identification
      double muSFid_BCDEF2 	= getMuonSF(h2_idSF_BCDEF, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFid_GH2 		= getMuonSF(h2_idSF_GH, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFid2 		= (muSFid_BCDEF2*lumi_BCDEF + muSFid_GH2*lumi_GH)/lumi; 
      //isolation 
      double muSFiso_BCDEF2 	= getMuonSF(h2_isoSF_BCDEF, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFiso_GH2 		= getMuonSF(h2_isoSF_GH, pfMuons[lepton2].p4.eta(), pfMuons[lepton2].p4.pt());
      double muSFiso2 		= (muSFiso_BCDEF2*lumi_BCDEF + muSFiso_GH2*lumi_GH)/lumi; 
      //tracking 
      double muSF2 = muSFtrig2*muSFid2*muSFiso2;	
      leptonSF = muSF1*muSF2;
    } 
    if(isEleChannel && !ev->isData){
      double eleSF1 =0;
      double ele_recoSF1       = getEleSF(h2_ele_recoSF, pfElectrons[lepton1].eleSCEta, pfElectrons[lepton1].p4.pt());
      //This is cut-based ID, we are using Heep ID
      //double ele_medium_idSF1  = getEleSF(h2_ele_medium_idSF, pfElectrons[lepton1].eleSCEta, pfElectrons[lepton1].p4.pt());
      double ele_trigSF1       = getEleTrigSF(h2_ele_trigSF, pfElectrons[lepton1].eleSCEta, pfElectrons[lepton1].p4.pt());
      double ele_heep_SF1      = getEleHeep2SF(tg_heep_SF, pfElectrons[lepton1].eleSCEta);
      eleSF1 = ele_recoSF1*ele_trigSF1*ele_heep_SF1;  

      double eleSF2 =0;
      double ele_recoSF2       = getEleSF(h2_ele_recoSF, pfElectrons[lepton2].eleSCEta, pfElectrons[lepton2].p4.pt());
      //double ele_medium_idSF2  = getEleSF(h2_ele_medium_idSF, pfElectrons[lepton2].eleSCEta, pfElectrons[lepton2].p4.pt());
      double ele_trigSF2       = getEleTrigSF(h2_ele_trigSF, pfElectrons[lepton2].eleSCEta,  pfElectrons[lepton2].p4.pt());
      double ele_heep_SF2      = getEleHeep2SF(tg_heep_SF, pfElectrons[lepton2].eleSCEta);
      eleSF2 = ele_recoSF2*ele_trigSF2*ele_heep_SF2;
      double eleSF = 1.0;
      leptonSF = eleSF1*eleSF2;
    }
    evtWeight *= leptonSF;
    fillHisto(outFile_, cutflowType, "", "leptonSF", 1000, 0, 100, leptonSF, 1 );
    string cutflowType_(cutflowType);
    cutflowType_ = cutflowType;
    
    //---------------------------------------------------//
    //get 4 vector for Z boson
    //---------------------------------------------------//
    MyLorentzVector vZ; 
    if(isMuChannel) vZ = pfMuons[lepton1].p4 + pfMuons[lepton2].p4;
    if(isEleChannel) vZ = pfElectrons[lepton1].p4 + pfElectrons[lepton2].p4;
    //if(j_init.size()==0) continue;
    //---------------------------------------------------//
    //Fill histos with for Control Plots
    //---------------------------------------------------//
    //fill histos for muon
    double leptonPt1 = 0.0;
    double leptonPt2 = 0.0; 
    double etaLep1 = 0.0;
    double etaLep2 =0.0;
    if(isMuChannel){
      leptonPt1 = pfMuons[lepton1].p4.pt();
      leptonPt2 = pfMuons[lepton2].p4.pt(); 
      etaLep1   =  pfMuons[lepton1].p4.eta();
      etaLep2   =  pfMuons[lepton2].p4.eta();
    }
    if(isEleChannel){
      leptonPt1 = pfElectrons[lepton1].p4.pt();
      leptonPt2 = pfElectrons[lepton2].p4.pt(); 
      etaLep1   =  pfElectrons[lepton1].p4.eta();
      etaLep2   =  pfElectrons[lepton2].p4.eta();
    }
    //fill histos for jets
    double count_CSVT_SF = 0.0;
    for(size_t ijet = 0; ijet < j_init.size(); ijet++){
      int ind_jet = j_init[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer);
      fillHisto(outFile_, cutflowType_, "ControlP","pt_jet", 500, 0, 10000, jetPt, evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","eta_jet", 50, -5, 5, pfJets[ind_jet].p4.eta(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","phi_jet", 50, -5, 5, pfJets[ind_jet].p4.phi(), evtWeight );
      fillHisto(outFile_, cutflowType_, "ControlP","mass_jet", 500, 0, 5000, pfJets[ind_jet].p4.M(), evtWeight );
    }
    fillHisto(outFile_, cutflowType_, "ControlP","final_multi_jet", 15, 0, 15, j_init.size(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","multiLep",  15, 0.5, 15.5, nLepton, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt1stLep", 500, 0, 10000, leptonPt1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pt2ndLep", 500, 0, 10000, leptonPt2, evtWeight );
    fillHisto2D(outFile_, cutflowType_,"ControlP", "ptLep1_ptLep2", 500, 0, 10000, leptonPt1,500, 0, 10000, leptonPt2, 1);
    fillHisto(outFile_, cutflowType_, "ControlP","eta1stLep", 50, -5, 5, etaLep1, evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta2ndLep", 50, -5, 5, etaLep2, evtWeight );
    //fill histos for Z boson
    fillHisto(outFile_, cutflowType_, "ControlP","pt_Z",  200, 0, 1000, vZ.Pt(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","pz_Z",  200, 0, 1000, vZ.Pz(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","eta_Z", 50, -5, 5, vZ.Rapidity(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","phi_Z", 50, -5, 5, vZ.Phi(), evtWeight );
    fillHisto(outFile_, cutflowType_, "ControlP","mll",   1000, 0, 200, vZ.M(), evtWeight );
    //fill histos for nvtx
    fillHisto(outFile_, cutflowType_, "ControlP","nvtx", 100, 0, 100, pri_vtxs, evtWeight );
    for(std::size_t n=0; n<Vertices.size(); n++){
      fillHisto(outFile_, cutflowType_, "ControlP","rhoAll", 100, 0, 100, Vertices[n].rhoAll, evtWeight );
    }
    nCutPass++;
    fillHisto(outFile_, cutflowType_, "", "cutflow", 20, 0.5, 20.5, nCutPass, evtWeight );

  }//event loop
  cout<<"Total events  = "<<nEntries<<endl;
  f->Close(); 
  delete f;
}

void Analyzer::processEvents(){ 
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+"outputFile_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel(9);
  
  //Local 
  TString pathLocal = "13TeV/inputDir/DYJetsToLL_Ntuple.root";
  CutFlowAnalysis(pathLocal, true, false, outFile_); 
  //CutFlowAnalysis(pathLocal, false, true, outFile_); 
  
  //T2
  //TString pathT2 = "/cms/store/user/rverma/ntuple_for2016MC_20190922/MC_20190922/DYJetsToLL_M50_MC_20190922/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M50_MC_20190922/190922_144942/0000/DYJetsToLL_M50_MC_20190922_Ntuple_1.root";
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/"+pathT2, true, false, outFile_);
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/"+pathT2, false, true, outFile_);
  //================
  //condor submission
  //================
  /*
  CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", true, false, outFile_);
  CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", false, true, outFile_);
  */
  outFile_->Write(); 
  outFile_->Close();
} 
