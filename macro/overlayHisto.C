
void overlayHisto(TString histName = "/base/Iso/ZTag/mlZmin_mlZmax", TString xaxis_title="mlZmin_mlZmax"){
gStyle->SetFrameLineWidth(3);
gStyle->SetOptStat(0);
gROOT->ForceStyle();
TFile *fSig250 = TFile::Open("all_ExLepMuMuZ_M250.root", "READ");
TFile *fSig2000 = TFile::Open("all_ExLepMuMuZ_M2000.root");
TFile *fSig4000 = TFile::Open("all_ExLepMuMuZ_M4000.root");
TFile *fDY = TFile::Open("all_DY_all.root");
TFile *fTT = TFile::Open("all_TT.root");
TFile *fVV = TFile::Open("all_VV.root");
TFile *fData = TFile::Open("all_muData.root");

TCanvas *c = new TCanvas();
gPad->SetBottomMargin(2.);
gPad->SetLeftMargin(0.11);

TH2F *hSig250 = (TH2F*)fSig250->Get(histName);
hSig250->SetLineWidth(2);
hSig250->SetMarkerColor(kRed);
hSig250->SetTitle("#mu#mu*->#mu#muZ->#mu#muqq");
hSig250->GetXaxis()->CenterTitle();
hSig250->GetYaxis()->CenterTitle();
hSig250->GetXaxis()->SetTitle("M_{lZ}^{min} [GeV]");
hSig250->GetYaxis()->SetTitle("M_{lZ}^{max} [GeV]");
hSig250->GetXaxis()->SetRangeUser(0, 5000);
hSig250->GetXaxis()->SetTitleOffset(1.2);
hSig250->GetYaxis()->SetTitleOffset(1.5);
hSig250->GetYaxis()->SetRangeUser(0, 7000);
hSig250->Draw();

TH2F *hSig2000 = (TH2F*)fSig2000->Get(histName);
hSig2000->SetLineWidth(2);
hSig2000->SetMarkerColor(kGreen);
hSig2000->Draw("same");

TH2F *hSig4000 = (TH2F*)fSig4000->Get(histName);
hSig4000->SetLineWidth(2);
hSig4000->SetMarkerColor(kBlue);
hSig4000->Draw("same");

TH2F *hDY = (TH2F*)fDY->Get(histName);
TH2F *hTT = (TH2F*)fTT->Get(histName);
TH2F *hVV = (TH2F*)fVV->Get(histName);
TH2F *hAllBkg = (TH2F*)hDY->Clone();
hAllBkg->Add(hTT);
hAllBkg->Add(hVV);
hAllBkg->SetLineWidth(2);
hAllBkg->SetMarkerColor(kCyan);
hAllBkg->Draw("same");

TH2F *hData = (TH2F*)fData->Get(histName);
hData->SetLineWidth(2);
hData->SetMarkerColor(kYellow +2);
//hData->Draw("same");

TLegend *leg = new TLegend( 0.60,0.20,0.80,0.40,"");
hSig250->SetLineColor(kRed);
hSig2000->SetLineColor(kGreen);
hSig4000->SetLineColor(kBlue);
hAllBkg->SetLineColor(kCyan);
hData->SetLineColor(kYellow +2);
leg->SetTextSize(0.03);
leg->AddEntry(hSig250,"Signal 250 GeV","lp");
leg->AddEntry(hSig2000,"Signal 2000 GeV","lp");
leg->AddEntry(hSig4000,"Signal 4000 GeV","lp");
leg->AddEntry(hAllBkg,"AllBkg","lp");
//leg->AddEntry(hData,"Data","lp");
leg->Draw();
}
