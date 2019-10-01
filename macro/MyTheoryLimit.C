
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
  graph->GetYaxis()->SetMoreLogLabels(true);
  graph->GetXaxis()->SetMoreLogLabels(true);
  //graph->GetYaxis()->SetRangeUser(yMin, yMax);
  //graph->GetXaxis()->SetRangeUser(0.44, 0.72);
  graph->GetYaxis()->SetLabelSize(0.035);
  graph->GetXaxis()->SetNdivisions(14);
  graph->GetYaxis()->SetNdivisions(10);
  graph->GetYaxis()->CenterTitle();
  graph->SetLineColor(color);
  graph->SetLineWidth(3);
  //graph->SetMarkerStyle(20);
  //graph->SetMarkerSize(1);
  graph->SetMarkerColor(color);
  return graph;
}

void MyTheoryLimit(){
  vector<double>lambVec250;             vector<double>lambVec500;
  lambVec250.push_back(861.4	);      lambVec500.push_back(167.4    );   
  lambVec250.push_back(98.63	);      lambVec500.push_back(35.9     );
  lambVec250.push_back(11.87	);      lambVec500.push_back(8.921    );
  lambVec250.push_back(2.455	);      lambVec500.push_back(2.247    );
  lambVec250.push_back(0.784	);      lambVec500.push_back(0.6843   );
  lambVec250.push_back(0.2793	);      lambVec500.push_back(0.2766   );
  lambVec250.push_back(0.1133	);      lambVec500.push_back(0.1373   );
  lambVec250.push_back(0.06922	);      lambVec500.push_back(0.08491  );
  lambVec250.push_back(0.01524	);      lambVec500.push_back(0.02411  );
  lambVec250.push_back(0.006653	);      lambVec500.push_back(0.009437 );
  lambVec250.push_back(0.002478	);      lambVec500.push_back(0.004609 );
  lambVec250.push_back(0.001263	);      lambVec500.push_back(0.002353 );
  lambVec250.push_back(0.0006454	);      lambVec500.push_back(0.001109 );
  lambVec250.push_back(0.0003316	);      lambVec500.push_back(0.000748 );

  vector<double>lambVec750;             vector<double>lambVec1000;
  lambVec750.push_back(59.4     );       lambVec1000.push_back(23.51     );   
  lambVec750.push_back(16.27    );       lambVec1000.push_back(7.798     );
  lambVec750.push_back(4.485    );       lambVec1000.push_back(2.778     );
  lambVec750.push_back(1.646    );       lambVec1000.push_back(1.084     );
  lambVec750.push_back(0.5934   );       lambVec1000.push_back(0.4112    );
  lambVec750.push_back(0.2133   );       lambVec1000.push_back(0.1598    );
  lambVec750.push_back(0.09433  );       lambVec1000.push_back(0.06778   );
  lambVec750.push_back(0.04673  );       lambVec1000.push_back(0.03527   );
  lambVec750.push_back(0.02302  );       lambVec1000.push_back(0.01221   );
  lambVec750.push_back(0.007783 );       lambVec1000.push_back(0.005419  );
  lambVec750.push_back(0.004525 );       lambVec1000.push_back(0.002615  );
  lambVec750.push_back(0.002142 );       lambVec1000.push_back(0.001787  );
  lambVec750.push_back(0.001279 );       lambVec1000.push_back(0.0009694 );
  lambVec750.push_back(0.0008643);       lambVec1000.push_back(0.0005933 );

  vector<double>lambVec1250;            vector<double>lambVec1500;        
  lambVec1250.push_back(11.83    );      lambVec1500.push_back(6.69     );
  lambVec1250.push_back(4.825    );      lambVec1500.push_back(2.705    );
  lambVec1250.push_back(1.579    );      lambVec1500.push_back(1.11     );
  lambVec1250.push_back(0.6337   );      lambVec1500.push_back(0.3795   );
  lambVec1250.push_back(0.2394   );      lambVec1500.push_back(0.1995   );
  lambVec1250.push_back(0.1085   );      lambVec1500.push_back(0.0855   );
  lambVec1250.push_back(0.05432  );      lambVec1500.push_back(0.0394   );
  lambVec1250.push_back(0.02816  );      lambVec1500.push_back(0.02071  );
  lambVec1250.push_back(0.008043 );      lambVec1500.push_back(0.005673 );
  lambVec1250.push_back(0.003708 );      lambVec1500.push_back(0.002339 );
  lambVec1250.push_back(0.001767 );      lambVec1500.push_back(0.001348 );
  lambVec1250.push_back(0.001057 );      lambVec1500.push_back(0.000795 );
  lambVec1250.push_back(0.0007054);      lambVec1500.push_back(0.0004033);
  lambVec1250.push_back(0.0004234);      lambVec1500.push_back(0.0002709);

  vector<double>lambVec1750;            vector<double>lambVec2000;       
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(2.276     );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(1.107     );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.4566    );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.1832    );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.09486   );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.04627   );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.0223    );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.01141   );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.003424  );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.001149  );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.0005401 );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.0002838 );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.0001541 );
  lambVec1750.push_back(0.0000 );        lambVec2000.push_back(0.00009769);

  vector<double>lambVec2500;             vector<double>lambVec3000;       
  lambVec2500.push_back(1.086     );     lambVec3000.push_back(0.5149    );
  lambVec2500.push_back(0.5275    );     lambVec3000.push_back(0.2822    );
  lambVec2500.push_back(0.2277    );     lambVec3000.push_back(0.1424    );
  lambVec2500.push_back(0.1155    );     lambVec3000.push_back(0.07183   );
  lambVec2500.push_back(0.05286   );     lambVec3000.push_back(0.03289   );
  lambVec2500.push_back(0.03009   );     lambVec3000.push_back(0.01738   );
  lambVec2500.push_back(0.01253   );     lambVec3000.push_back(0.008909  );
  lambVec2500.push_back(0.008244  );     lambVec3000.push_back(0.004948  );
  lambVec2500.push_back(0.00212   );     lambVec3000.push_back(0.001507  );
  lambVec2500.push_back(0.0006827 );     lambVec3000.push_back(0.0004379 );
  lambVec2500.push_back(0.000284  );     lambVec3000.push_back(0.0001703 );
  lambVec2500.push_back(0.0001393 );     lambVec3000.push_back(0.00008872);
  lambVec2500.push_back(0.00007738);     lambVec3000.push_back(0.00003888);
  lambVec2500.push_back(0.00005921);     lambVec3000.push_back(0.00002767);
  
  vector<string>lambVecLabel;
  lambVecLabel.push_back("#lambda = 250  GeV");
  lambVecLabel.push_back("#lambda = 500  GeV");
  lambVecLabel.push_back("#lambda = 750  GeV");
  lambVecLabel.push_back("#lambda = 1000 GeV");
  lambVecLabel.push_back("#lambda = 1250 GeV");
  lambVecLabel.push_back("#lambda = 1500 GeV");
  lambVecLabel.push_back("#lambda = 1750 GeV");
  lambVecLabel.push_back("#lambda = 2000 GeV");
  lambVecLabel.push_back("#lambda = 2500 GeV");
  lambVecLabel.push_back("#lambda = 3000 GeV");
  lambVecLabel.push_back("#lambda = 3500 GeV");
  lambVecLabel.push_back("#lambda = 4000 GeV");
  lambVecLabel.push_back("#lambda = 4500 GeV");
  lambVecLabel.push_back("#lambda = 5000 GeV");

  vector<vector<double>> lambdaVecAll;
  lambdaVecAll.push_back(lambVec250);
  lambdaVecAll.push_back(lambVec500);
  lambdaVecAll.push_back(lambVec750);
  lambdaVecAll.push_back(lambVec1000);
  lambdaVecAll.push_back(lambVec1250);
  lambdaVecAll.push_back(lambVec1500);
  lambdaVecAll.push_back(lambVec1750);
  lambdaVecAll.push_back(lambVec2000);
  lambdaVecAll.push_back(lambVec2500);
  lambdaVecAll.push_back(lambVec3000);

  vector<double> massVec;
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
  
  TFile *outFile = new TFile("theoryXss.root", "RECREATE");
  TCanvas * can = new TCanvas();
  can->cd();
  gPad->SetGridy();
  gPad->SetGridx();
  gPad->SetLogy();
  gPad->SetLogx();
  TLegend* leg = new TLegend(0.15,0.78,0.90,0.88,NULL,"brNDC");
  leg->SetNColumns(3);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  for(int i = 0; i<lambdaVecAll.size(); i++){
    TGraph * myGraph = makeGraph(massVec, lambdaVecAll.at(i));
    myGraph->Write();
    if(isMuChannel)decorateGraph(myGraph, "M_{l^{*}} [GeV]", "#sigma (fb) #times BR", "Theoretical cross section", 0.75, 1.15, i+1);
    if(i==0) myGraph->Draw("ALP");
    else myGraph->Draw("LPsame");
    leg->AddEntry(myGraph, TString(lambVecLabel.at(i)),"PL");
  }
  outFile->Close();
  leg->Draw();
  can->SaveAs("theoryXss.pdf");
}


