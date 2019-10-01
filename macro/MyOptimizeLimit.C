
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

void MyOptimizeLimit(){
  vector<double> sigMass;
  vector<double> limitTau45;
  vector<double> limitTau55;
  vector<double> limitTau60;
  vector<double> limitTau65;
  
  sigMass.push_back(250);    
  sigMass.push_back(500);    
  sigMass.push_back(750);    
  sigMass.push_back(1000);   
  sigMass.push_back(1250);   
  sigMass.push_back(1500);   
  sigMass.push_back(1750);   
  sigMass.push_back(2000);   
  sigMass.push_back(2500);   
  sigMass.push_back(3000);   
  sigMass.push_back(3500);   
  sigMass.push_back(4000);   
  sigMass.push_back(4500);   
  sigMass.push_back(5000);   
  
  limitTau45.push_back(22.95);
  limitTau45.push_back(10.82);
  limitTau45.push_back(3.591);
  limitTau45.push_back(1.223);
  limitTau45.push_back(0.754);
  limitTau45.push_back(0.543);
  limitTau45.push_back(0.551);
  limitTau45.push_back(0.570);
  limitTau45.push_back(0.601);
  limitTau45.push_back(0.630);
  limitTau45.push_back(0.665);
  limitTau45.push_back(0.704);
  limitTau45.push_back(0.728);
  limitTau45.push_back(0.759);
  
  limitTau55.push_back(21.42);
  limitTau55.push_back(10.50);
  limitTau55.push_back(3.563);
  limitTau55.push_back(1.439);
  limitTau55.push_back(0.859);
  limitTau55.push_back(0.763);
  limitTau55.push_back(0.765);
  limitTau55.push_back(0.757);
  limitTau55.push_back(0.787);
  limitTau55.push_back(0.808);
  limitTau55.push_back(0.837);
  limitTau55.push_back(0.871);
  limitTau55.push_back(0.891);
  limitTau55.push_back(0.915);
  
  limitTau60.push_back(21.28);     
  limitTau60.push_back(10.50);     
  limitTau60.push_back(3.605);     
  limitTau60.push_back(1.554);     
  limitTau60.push_back(0.854);     
  limitTau60.push_back(0.770);     
  limitTau60.push_back(0.765);     
  limitTau60.push_back(0.792);     
  limitTau60.push_back(0.814);     
  limitTau60.push_back(0.836);     
  limitTau60.push_back(0.857);     
  limitTau60.push_back(0.885);     
  limitTau60.push_back(0.904);     
  limitTau60.push_back(0.923);     
  
  limitTau65.push_back(21.02);
  limitTau65.push_back(10.73);
  limitTau65.push_back(3.742);
  limitTau65.push_back(1.531);
  limitTau65.push_back(0.896);
  limitTau65.push_back(0.813);
  limitTau65.push_back(0.810);
  limitTau65.push_back(0.824);
  limitTau65.push_back(0.844);
  limitTau65.push_back(0.861);
  limitTau65.push_back(0.881);
  limitTau65.push_back(0.906);
  limitTau65.push_back(0.922);
  limitTau65.push_back(0.938);

  TCanvas * can = new TCanvas();
  TLegend* leg = new TLegend(0.50,0.65,0.85,0.85,NULL,"brNDC");
  gStyle->SetFrameLineWidth(2);
  gPad->SetGridy();
  gPad->SetGridx();
  gPad->SetLogy();
  gPad->SetLogx();
  //leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  
  TGraph * graphTau45 = makeGraph(sigMass, limitTau45);
  TGraph * graphTau55 = makeGraph(sigMass, limitTau55);
  TGraph * graphTau60 = makeGraph(sigMass, limitTau60);
  TGraph * graphTau65 = makeGraph(sigMass, limitTau65);
   
  TString yTitle = " 95% CL limit on #sigma (fb) #times BR(#mu#mu*#rightarrow #mu#muZ #rightarrow 2#mu2q)";
  decorateGraph(graphTau45, "M_{l^{*}}", yTitle, "#mu-channel", 0.5, 25, kBlue);
  decorateGraph(graphTau55, "M_{l^{*}}", yTitle, "#mu-channel", 0.5, 25, kBlack);
  decorateGraph(graphTau60, "M_{l^{*}}", yTitle, "#mu-channel", 0.5, 25, kGreen);
  decorateGraph(graphTau65, "M_{l^{*}}", yTitle, "#mu-channel", 0.5, 25, kRed);

  graphTau45->Draw();
  graphTau55->Draw("same");
  graphTau60->Draw("same");
  graphTau65->Draw("same");

  leg->AddEntry(graphTau45, "Tau21 < 0.45","PL");
  leg->AddEntry(graphTau55, "Tau21 < 0.55","PL");
  leg->AddEntry(graphTau60, "Tau21 < 0.60","PL");
  leg->AddEntry(graphTau65, "Tau21 < 0.65","PL");
  leg->Draw();
  can->SaveAs("tau21_limit.png");
}


