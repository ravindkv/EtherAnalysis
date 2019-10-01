#include <iostream>
#include <fstream>
#include <iomanip>

void makeLimitTable(TString CHANNEL, TString zTagDir, TString label,
         bool obs, bool isOut, ofstream & outFile)
  {

  float X[] = {250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
  float obsY[]      = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY[]      = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
		                                              
  float expX1sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expX1sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sL_[]  = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY1sH_[]  = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
                                                              
		                                              
  float expX2sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expX2sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY2sL[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};
  float expY2sH[]   = {0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0};

  int nMassPoints = 14;
  TString ch_hist = CHANNEL;
  TString massFiles [14] = {
  "Mass250/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH250.root",
  "Mass500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH500.root",
  "Mass750/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH750.root",
  "Mass1000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1000.root",
  "Mass1250/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1250.root",
  "Mass1500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1500.root",
  "Mass1750/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH1750.root",
  "Mass2000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH2000.root",
  "Mass2500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH2500.root",
  "Mass3000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH3000.root",
  "Mass3500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH3500.root",
  "Mass4000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH4000.root",
  "Mass4500/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH4500.root",
  "Mass5000/higgsCombine_llstar_llZ_llq_13TeV_"+ch_hist+".AsymptoticLimits.mH5000.root"};
  
  double maxY = 1.0;
  for(unsigned int i = 0 ; i < nMassPoints; i++){
    TFile f("output/"+CHANNEL+"/"+zTagDir+"/"+massFiles[i],"READ"); 
    if(f.IsZombie()){
      cout << "Cannot open file for " << string(CHANNEL.Data()) << " and mass " << X[i] << endl;
      continue;
    }
    vector<double>sigXss;
    if(CHANNEL=="mu"){
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
    if(CHANNEL=="ele"){
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
    if(CHANNEL=="mu_ele"){
      sigXss.push_back(0.00427     + 0.004407)  ;
      sigXss.push_back(0.00291     + 0.002916)  ;
      sigXss.push_back(0.001761    + 0.001607)  ;
      sigXss.push_back(0.001177    + 0.001079)  ;
      sigXss.push_back(0.0007263   + 0.0007171) ;
      sigXss.push_back(0.0004267   + 0.0004856) ;
      sigXss.push_back(0.0002654   + 0.0002907) ;
      sigXss.push_back(0.0002021   + 0.0001908) ;
      sigXss.push_back(0.00006755  + 0.00005731);
      sigXss.push_back(0.00002223  + 0.00002497);
      sigXss.push_back(0.000008166 + 0.000009023);
      sigXss.push_back(0.000003209 + 0.000003897);
      sigXss.push_back(0.000001133 + 0.000001243);
      sigXss.push_back(0.000000452 + 0.000000461);
    }

    Double_t r;
    TTree* limit = (TTree*)f.Get("limit");
    limit->SetBranchAddress("limit",&r);
    
    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);
      //multiply by xss
      r = 1000*r*sigXss[i];
      if(k==0) expY2sL[i] = r;
      if(k==1) expY1sL[i] = r;
      if(k==1) expY1sL_[i] = r;
      if(k==2) expY[i]    = r;
      if(k==3) expY1sH[i] = r;
      if(k==3) expY1sH_[i] = r;
      if(k==4) expY2sH[i] = r;
      if(k==5) obsY[i]    = r;
    }
    if(massFiles[i].Contains("5000")) maxY = expY2sH[i];
  }
  cout<<std::setprecision(4)<<endl;
  cout<<"===========: " <<zTagDir<<endl;
  cout<<"Mass:"<<setw(15)<<"base value"<<setw(15)<<"-2 #sigma"<<setw(15)<<"-1 #sigma"<<setw(15)<<"+1 #sigma"<<setw(15)<<"+2 #sigma"<<endl; 

  //make table
  outFile<<"\\begin{table}"<<endl; 
  outFile<<"\\begin{center}"<<endl; 
  outFile<<"\\begin{tabular}{cccccc}"<<endl; 
  //outFile<<"\\begin{tabular}{ |l|l|l|l|l|l|}"<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"$M_{l^*}$"<<" & "<< "-2$\\sigma$ & -1$\\sigma$ & median & +1$\\sigma$ & +1$\\sigma$\\\\ "<<endl;
  outFile<<"(GeV) & "<< "(fb) & "<< "(fb) & "<< "(fb) & "<<"(fb) & "<<"(fb) "<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  for(int i1 = 0 ; i1 < nMassPoints ; i1++){
  cout<<X[i1]<<setw(15)<<expY[i1]<<setw(15)<<expY2sL[i1]<<setw(15)<< expY1sL[i1]<<setw(15)<<expY1sH[i1]<<setw(15)<<expY2sH[i1]<<endl; 
  outFile<<X[i1]<<std::setprecision(3)<<" & "<<expY2sL[i1]<<" & "<< expY1sL[i1]<<" & "<<expY[i1]<<" & "<<expY1sH[i1]<<" & "<<expY2sH[i1]<<"\\\\"<<endl; 
  }
  outFile<<"\\hline "<<endl;
  outFile<<"\\end{tabular}"<<endl; 
  TString ch_name = "mu";
  if(CHANNEL=="mu") ch_name = "muon";
  if(CHANNEL=="ele") ch_name = "electron";
  if(CHANNEL=="mu_ele") ch_name = "lepton";
  ///outFile<<"\\caption{95\\% CL exclusion limit for "+ch_name+" channel.}"<<endl; 
  outFile<<"\\caption{95\\% CL exclusion limit after different L-cut $(M_{lZ}^{min} < X, M_{lZ}^{max} > Y + "+label+")$ for "+CHANNEL+" channel.}"<<endl;
  outFile<<"\\label{tab:limit_"+ch_name+"}"<<endl;
  outFile<<"\\end{center}"<<endl; 
  outFile<<"\\end{table}"<<endl; 
  outFile<<endl;
  
}

void MyMakeLimitTable(){
  ofstream outFile; 
  outFile.open("output/limitTable.tex"); 
  //outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;  
  outFile<<"\\documentclass[]{article}"<<endl;  
  outFile<<"\\pagestyle{empty}"<<endl;  
  outFile<<"\\usepackage{epsfig}"<<endl;  
  outFile<<"\\usepackage{amsmath}"<<endl;  
  outFile<<"\\usepackage{array}"<<endl;  
  outFile<<"\\usepackage{multirow}"<<endl;  
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;  
  outFile<<"\\begin{document}"<<endl;  
  outFile<<""<<endl;
  //muon
  //makeLimitTable("mu", "ZTag1", "0",   false, true, outFile);
  //makeLimitTable("mu", "ZTag2", "100", false, true, outFile);
  //makeLimitTable("mu", "ZTag3", "200", false, true, outFile);
  //makeLimitTable("mu", "ZTag4", "350", false, true, outFile);
  //makeLimitTable("mu", "ZTag5", "400", false, true, outFile);
  //makeLimitTable("mu", "ZTag6", "500", false, true, outFile);
  //makeLimitTable("mu", "ZTag7", "600", false, true, outFile);
  //makeLimitTable("mu", "ZTag8", "700", false, true, outFile);
  //makeLimitTable("mu", "ZTag9", "800", false, true, outFile);
  
  //electron 
  makeLimitTable("ele", "ZTag1", "0",   false, true, outFile);
  //makeLimitTable("ele", "ZTag2", "100", false, true, outFile);
  //makeLimitTable("ele", "ZTag3", "200", false, true, outFile);
  //makeLimitTable("ele", "ZTag4", "350", false, true, outFile);
  //makeLimitTable("ele", "ZTag5", "400", false, true, outFile);
  //makeLimitTable("ele", "ZTag6", "500", false, true, outFile);
  //makeLimitTable("ele", "ZTag7", "600", false, true, outFile);
  //makeLimitTable("ele", "ZTag8", "700", false, true, outFile);
  //makeLimitTable("ele", "ZTag9", "800", false, true, outFile);
  //makeLimitTable("ele", "ZTag9",  false, true, outFile);
  
  //lepton 
  //makeLimitTable("mu_ele", false, true, outFile);
  outFile<<"\\end{document}"<<endl;  
  outFile.close(); 
}

