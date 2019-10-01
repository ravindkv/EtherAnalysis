
TH2F * decorateHist2D(TH2F* hist2D, TString title, int color){
hist2D->SetLineWidth(2);
hist2D->SetMarkerColor(kRed);
hist2D->SetTitle(title);
hist2D->GetXaxis()->CenterTitle();
hist2D->GetYaxis()->CenterTitle();
hist2D->GetXaxis()->SetTitle("#tau_{21}");
hist2D->GetYaxis()->SetTitle("M_{l1l2} [GeV]");
hist2D->GetXaxis()->SetRangeUser(0, 1);
hist2D->GetXaxis()->SetTitleOffset(1.2);
hist2D->GetYaxis()->SetTitleOffset(1.5);
hist2D->GetYaxis()->SetRangeUser(50, 500);
return hist2D;
}

void MyBkgEstHist2D(TString histName = "/base/Iso/ZTag/ml1l2_tau21", TString xaxis_title="ml1l2_tau21"){
gStyle->SetFrameLineWidth(2);
gStyle->SetOptStat(0);
gROOT->ForceStyle();
TFile *fDY = TFile::Open("all_DY.root");
TFile *fTT = TFile::Open("all_TT.root");
TFile *fVV = TFile::Open("all_VV.root");
TFile *fData = TFile::Open("all_Data.root");

TCanvas *c = new TCanvas();
c->Divide(2,2);
gPad->SetBottomMargin(2.);
gPad->SetLeftMargin(0.11);

c->cd(1);
TH2F *hDY = (TH2F*)fDY->Get(histName);
decorateHist2D(hDY, "DY + jets", 1);
hDY->Draw("colz");

c->cd(2);
TH2F *hTT = (TH2F*)fTT->Get(histName);
decorateHist2D(hTT, "ttbar + jets" , 1);
hTT->Draw("colz");

c->cd(3);
TH2F *hVV = (TH2F*)fVV->Get(histName);
decorateHist2D(hVV, "VV" ,1);
hVV->Draw("colz");

c->cd(4);
TH2F *hAllBkg = (TH2F*)hDY->Clone();
hAllBkg->Add(hTT);
hAllBkg->Add(hVV);
decorateHist2D(hAllBkg, "All Bkg" ,1);
hAllBkg->Draw("colz");
c->SaveAs("bkgEstHist2D.pdf");
}
