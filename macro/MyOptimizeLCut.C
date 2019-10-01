
//CHANNEL
bool isMuChannel = true;
bool isEleChannel = false;

//--------------------------------------------
//function to make graph from two array
//--------------------------------------------
TGraph* makeGraph(vector<double> Xarray, vector<double> Yarray){
  const int n = Xarray.size();
  double x[n], y[n];
   for (Int_t i=0;i<n;i++) {
     x[i] = Xarray[i];
     y[i] = Yarray[i] ;
   }
   TGraph* gr = new TGraph(n,x,y);
   return gr;
}

//--------------------------------------------
//function to decorate Graph  
//--------------------------------------------
TGraph* decorateGraph(TGraph *graph, TString xTitle, TString yTitle, TString myTitle, double yMin, double yMax, int color){
  graph->SetTitle(myTitle);
  graph->GetYaxis()->SetTitleOffset(1.30);
  graph->GetXaxis()->SetTitle(xTitle);
  graph->GetXaxis()->SetTitleSize(0.04);
  graph->GetXaxis()->SetLabelSize(0.035);
  graph->GetYaxis()->SetTitle(yTitle);
  graph->GetYaxis()->SetTitleSize(0.04);
  graph->GetYaxis()->SetRangeUser(yMin, yMax);
  //graph->GetXaxis()->SetRangeUser(0.44, 0.72);
  graph->GetYaxis()->SetMoreLogLabels(true);
  graph->GetXaxis()->SetMoreLogLabels(true);
  graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetXaxis()->SetNdivisions(14);
  graph->GetYaxis()->SetNdivisions(10);
  graph->GetYaxis()->CenterTitle();
  graph->SetLineColor(color);
  graph->SetLineWidth(3);
  graph->SetMarkerSize(20);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(color);
  return graph;
}

void MyOptimizeLCut(){
  vector<double> massVec;
  vector<double>limitVec1;
  vector<double>limitVec2;
  vector<double>limitVec3;
  vector<double>limitVec4;
  vector<double>limitVec5;
  vector<double>limitVec6;
  massVec.push_back(250 );       
  massVec.push_back(500 );       
  massVec.push_back(750 );       
  massVec.push_back(1000);       
  massVec.push_back(1250);       
  massVec.push_back(1500);       
  massVec.push_back(1750);       
  massVec.push_back(2000);       
  massVec.push_back(2500);       
  massVec.push_back(3000);       
  massVec.push_back(3500);       
  massVec.push_back(4000);       
  massVec.push_back(4500);       
  massVec.push_back(5000);       
  
  if(isMuChannel){
    limitVec1.push_back(7.639 );       limitVec2.push_back(20.35 );       
    limitVec1.push_back(4.842 );       limitVec2.push_back(8.093 );      
    limitVec1.push_back(2.215 );       limitVec2.push_back(2.146 );      
    limitVec1.push_back( 1.209);       limitVec2.push_back( 1.563);      
    limitVec1.push_back(0.7008);       limitVec2.push_back(0.8029);      
    limitVec1.push_back(0.6501);       limitVec2.push_back(0.6084);      
    limitVec1.push_back(0.6428);       limitVec2.push_back(0.5992);      
    limitVec1.push_back(  0.66);       limitVec2.push_back(0.5874);      
    limitVec1.push_back(0.6628);       limitVec2.push_back(0.6185);      
    limitVec1.push_back(0.6655);       limitVec2.push_back(0.6211);      
    limitVec1.push_back(0.6717);       limitVec2.push_back(0.6247);      
    limitVec1.push_back(0.6819);       limitVec2.push_back(0.6314);      
    limitVec1.push_back(0.6832);       limitVec2.push_back(0.6311);      
    limitVec1.push_back( 0.687);       limitVec2.push_back(0.6342);      
    
    limitVec3.push_back(9.307 );      limitVec4.push_back(10.94 ); 
    limitVec3.push_back(5.388 );      limitVec4.push_back( 5.57 ); 
    limitVec3.push_back(2.428 );      limitVec4.push_back(2.717 ); 
    limitVec3.push_back( 1.237);      limitVec4.push_back( 1.462); 
    limitVec3.push_back(0.7433);      limitVec4.push_back(0.7887); 
    limitVec3.push_back(0.6801);      limitVec4.push_back(0.7067); 
    limitVec3.push_back(0.6614);      limitVec4.push_back(0.6884); 
    limitVec3.push_back(0.6852);      limitVec4.push_back(0.7105); 
    limitVec3.push_back(0.6861);      limitVec4.push_back(0.7114); 
    limitVec3.push_back(0.6891);      limitVec4.push_back(0.7141); 
    limitVec3.push_back(0.6962);      limitVec4.push_back(0.7206); 
    limitVec3.push_back(0.7076);      limitVec4.push_back(0.7317); 
    limitVec3.push_back(0.7093);      limitVec4.push_back(0.7342); 
    limitVec3.push_back(0.7146);      limitVec4.push_back(0.7399); 
    
    limitVec5.push_back(22.35 );      limitVec6.push_back(7.039 );      
    limitVec5.push_back(9.094 );      limitVec6.push_back(4.308 );      
    limitVec5.push_back(2.215 );      limitVec6.push_back(2.146 );      
    limitVec5.push_back( 1.522);      limitVec6.push_back( 1.214);      
    limitVec5.push_back( 0.749);      limitVec6.push_back(0.7462);      
    limitVec5.push_back(0.6534);      limitVec6.push_back( 0.605);      
    limitVec5.push_back(0.6469);      limitVec6.push_back(0.5951);      
    limitVec5.push_back( 0.641);      limitVec6.push_back(0.6047);      
    limitVec5.push_back(0.6713);      limitVec6.push_back(0.6101);      
    limitVec5.push_back(0.6752);      limitVec6.push_back(0.6113);      
    limitVec5.push_back(0.6798);      limitVec6.push_back(0.6186);      
    limitVec5.push_back(0.6883);      limitVec6.push_back(0.6274);      
    limitVec5.push_back(0.6866);      limitVec6.push_back(0.6277);      
    limitVec5.push_back(0.6898);      limitVec6.push_back(0.6319);      
  }
  if(isEleChannel){
    limitVec1.push_back(13.84 );       limitVec2.push_back(37.46 );       
    limitVec1.push_back(9.113 );       limitVec2.push_back(15.31 );      
    limitVec1.push_back(4.281 );       limitVec2.push_back(4.143 );      
    limitVec1.push_back( 2.276);       limitVec2.push_back( 2.875);      
    limitVec1.push_back( 1.317);       limitVec2.push_back( 1.451);      
    limitVec1.push_back( 1.222);       limitVec2.push_back( 1.134);      
    limitVec1.push_back( 1.222);       limitVec2.push_back( 1.131);      
    limitVec1.push_back( 1.258);       limitVec2.push_back( 1.115);      
    limitVec1.push_back( 1.257);       limitVec2.push_back( 1.161);      
    limitVec1.push_back( 1.248);       limitVec2.push_back( 1.155);      
    limitVec1.push_back( 1.254);       limitVec2.push_back( 1.159);      
    limitVec1.push_back(  1.28);       limitVec2.push_back( 1.183);      
    limitVec1.push_back( 1.285);       limitVec2.push_back( 1.186);      
    limitVec1.push_back( 1.291);       limitVec2.push_back( 1.191);      
    
    limitVec3.push_back(16.94 );      limitVec4.push_back( 19.9 ); 
    limitVec3.push_back(10.16 );      limitVec4.push_back(10.52 ); 
    limitVec3.push_back(4.708 );      limitVec4.push_back(5.298 ); 
    limitVec3.push_back( 2.369);      limitVec4.push_back( 2.816); 
    limitVec3.push_back( 1.417);      limitVec4.push_back( 1.513); 
    limitVec3.push_back( 1.282);      limitVec4.push_back( 1.332); 
    limitVec3.push_back( 1.258);      limitVec4.push_back( 1.308); 
    limitVec3.push_back( 1.318);      limitVec4.push_back( 1.365); 
    limitVec3.push_back( 1.311);      limitVec4.push_back( 1.358); 
    limitVec3.push_back( 1.305);      limitVec4.push_back( 1.355); 
    limitVec3.push_back( 1.306);      limitVec4.push_back( 1.356); 
    limitVec3.push_back( 1.335);      limitVec4.push_back( 1.383); 
    limitVec3.push_back( 1.339);      limitVec4.push_back( 1.387); 
    limitVec3.push_back( 1.344);      limitVec4.push_back( 1.393); 
    
    limitVec5.push_back(41.04 );      limitVec6.push_back(12.81  );      
    limitVec5.push_back(17.27 );      limitVec6.push_back( 8.11  );      
    limitVec5.push_back(4.281 );      limitVec6.push_back(4.156  );      
    limitVec5.push_back( 2.807);      limitVec6.push_back( 2.284 );      
    limitVec5.push_back( 1.361);      limitVec6.push_back( 1.401 );      
    limitVec5.push_back( 1.225);      limitVec6.push_back( 1.134 );      
    limitVec5.push_back( 1.222);      limitVec6.push_back( 1.126 );      
    limitVec5.push_back( 1.216);      limitVec6.push_back(  1.16 );      
    limitVec5.push_back( 1.261);      limitVec6.push_back( 1.153 );      
    limitVec5.push_back( 1.255);      limitVec6.push_back( 1.149 );      
    limitVec5.push_back( 1.261);      limitVec6.push_back( 1.153 );      
    limitVec5.push_back( 1.288);      limitVec6.push_back( 1.177 );      
    limitVec5.push_back(  1.29);      limitVec6.push_back( 1.181 );      
    limitVec5.push_back( 1.295);      limitVec6.push_back( 1.187 );      
  }
  
  vector<vector<double>> limitVecAll;
  limitVecAll.push_back(limitVec1);
  limitVecAll.push_back(limitVec2);
  limitVecAll.push_back(limitVec3);
  limitVecAll.push_back(limitVec4);
  limitVecAll.push_back(limitVec5);
  limitVecAll.push_back(limitVec6);

  vector<string>lCutVec;
  lCutVec.push_back("L-cut: set1");
  lCutVec.push_back("L-cut: set2");
  lCutVec.push_back("L-cut: set3");
  lCutVec.push_back("L-cut: set4");
  lCutVec.push_back("L-cut: set5");
  lCutVec.push_back("L-cut: set6");
  TCanvas * can = new TCanvas();
  TLegend* leg = new TLegend(0.50,0.55,0.85,0.85,NULL,"brNDC");
  gPad->SetGridy();
  gPad->SetGridx();
  gPad->SetLogy();
  gPad->SetLogx();
  leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  for(int i = 0; i<limitVecAll.size(); i++){
   can->cd();
   TGraph * myGraph = makeGraph(massVec, limitVecAll.at(i));
   int color = i+1;
   if(i==9) color = 433;
   if(isMuChannel){
     TString fullProcess = "#mu#mu*#rightarrow #mu#muZ #rightarrow 2#mu2q";   
     decorateGraph(myGraph, "M_{l^{*}} [GeV]", "#sigma (fb) #times BR("+fullProcess+")", "Expected limits for different L-cut :  #mu-channel", 0.1, 50, color);
   }
   if(isEleChannel){
     TString fullProcess = "ee*#rightarrow eeZ #rightarrow 2e2q";   
     decorateGraph(myGraph, "M_{l^{*}} [GeV]", "#sigma (fb) #times BR("+fullProcess+")", "Expected limits for different L-cut :  e-channel", 0.1, 50, color);
   }
   if(i==0) myGraph->Draw("ALP");
   else myGraph->Draw("LPsame");
   leg->AddEntry(myGraph, TString(lCutVec[i]),"PL");
  }
  leg->Draw();
  if(isMuChannel)can->SaveAs("limit_lCut_mu.pdf");
  if(isEleChannel)can->SaveAs("limit_lCut_ele.pdf");
}

