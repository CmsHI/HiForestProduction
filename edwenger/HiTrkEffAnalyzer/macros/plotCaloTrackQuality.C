//-----------------------------------------------------------------------------
// Andre S. Yoon, Nov 22 2010
//
// Simple macro that plots the track qualities of track-calo compatibility 
// rejected and accepted tracks in comparison for |eta| < eta' 
// 
//-----------------------------------------------------------------------------


void plotCaloTrackQuality() {

  

  TFile *f = new TFile("../root_files/data/trkhistDATA_nov19_trkhist_cAll_v2.root");

  TString hDir = "pfCandidateAnalyzer";

  TString hN_nHits_acc; TString hN_ptErr_acc; TString hN_nAlgo_acc; TString hN_D0_acc; TString hN_DZ_acc; TString hN_D0Err_acc; 
  TString hN_DZErr_acc; TString hN_D0PerErr_acc; TString hN_DZPerErr_acc;

  TString hN_nHits_rej;TString hN_ptErr_rej; TString hN_nAlgo_rej; TString hN_D0_rej; TString hN_DZ_rej; TString hN_D0Err_rej;
  TString hN_DZErr_rej;TString hN_D0PerErr_rej; TString hN_DZPerErr_rej;

  
  hN_nHits_acc = "hHitsEtaAccept"; hN_ptErr_acc = "hPtErrEtaAccept"; hN_nAlgo_acc = "hAlgoEtaAccept"; hN_D0_acc = "hD0EtaAccept"; 
  hN_DZ_acc = "hDZEtaAccept"; hN_D0Err_acc = "hD0ErrEtaAccept"; hN_DZErr_acc = "hDZErrEtaAccept"; hN_D0PerErr_acc = "hD0PerErrEtaAccept"; 
  hN_DZPerErr_acc = "hDZPerErrEtaAccept";

  hN_nHits_rej = "hHitsEtaReject"; hN_ptErr_rej= "hPtErrEtaReject"; hN_nAlgo_rej = "hAlgoEtaReject"; hN_D0_rej = "hD0EtaReject";
  hN_DZ_rej ="hDZEtaReject"; hN_D0Err_rej = "hD0ErrEtaReject"; hN_DZErr_rej = "hDZErrEtaReject"; hN_D0PerErr_rej = "hD0PerErrEtaReject";
  hN_DZPerErr_rej = "hDZPerErrEtaReject";

  // accepted tracks
  TH2F *hNhitsEta_acc = (TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_nHits_acc.Data()));
  TH2F *hPtErrEta_acc = (TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_ptErr_acc.Data()));
  TH2F *hD0PerErrEta_acc = (TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_D0PerErr_acc.Data()));
  TH2F *hDZPerErrEta_acc = (TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_DZPerErr_acc.Data()));
  
  // rejected tracks
  TH2F *hNhitsEta_rej = (TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_nHits_rej.Data()));
  TH2F *hPtErrEta_rej =(TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_ptErr_rej.Data()));
  TH2F *hD0PerErrEta_rej = (TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_D0PerErr_rej.Data()));
  TH2F *hDZPerErrEta_rej = (TH2F*) f->Get(Form("%s/%s",hDir.Data(),hN_DZPerErr_rej.Data()));
  

   // bin ranges (eta) 
   double small = 0.01;

   int eta_min_full = hNhitsEta_acc->GetXaxis()->FindBin(-2.4+small);
   int eta_max_full = hNhitsEta_acc->GetXaxis()->FindBin(2.4-small);

   int eta_min_narr = hNhitsEta_acc->GetXaxis()->FindBin(-1.0+small);
   int eta_max_narr = hNhitsEta_acc->GetXaxis()->FindBin(1.0-small);

   // projection eta < 2.4  and normalization
   TH1D *hNhits_acc_full = (TH1D*) hNhitsEta_acc->ProjectionY("hNhits_acc_full",eta_min_full,eta_max_full,"e"); 
   TH1D *hPtErr_acc_full = (TH1D*) hPtErrEta_acc->ProjectionY("hPtErr_acc_full",eta_min_full,eta_max_full,"e");
   TH1D *hD0PerErr_acc_full  = (TH1D*) hD0PerErrEta_acc->ProjectionY("hD0PerErr_acc_full",eta_min_full,eta_max_full,"e");
   TH1D *hDZPerErr_acc_full  = (TH1D*) hDZPerErrEta_acc->ProjectionY("hDZPerErr_acc_full",eta_min_full,eta_max_full,"e");

   TH1D *hNhits_rej_full = (TH1D*) hNhitsEta_rej->ProjectionY("hNhits_rej_full",eta_min_full,eta_max_full,"e");
   TH1D *hPtErr_rej_full = (TH1D*) hPtErrEta_rej->ProjectionY("hPtErr_rej_full",eta_min_full,eta_max_full,"e");
   TH1D *hD0PerErr_rej_full  = (TH1D*) hD0PerErrEta_rej->ProjectionY("hD0PerErr_rej_full",eta_min_full,eta_max_full,"e");
   TH1D *hDZPerErr_rej_full  = (TH1D*) hDZPerErrEta_rej->ProjectionY("hDZPerErr_rej_full",eta_min_full,eta_max_full,"e");

   hNhits_acc_full->Scale(1./hNhits_acc_full->GetEntries());
   hPtErr_acc_full->Scale(1./hPtErr_acc_full->GetEntries());
   hD0PerErr_acc_full->Scale(1./hD0PerErr_acc_full->GetEntries());
   hDZPerErr_acc_full->Scale(1./hDZPerErr_acc_full->GetEntries());

   hNhits_rej_full->Scale(1./hNhits_rej_full->GetEntries());
   hPtErr_rej_full->Scale(1./hPtErr_rej_full->GetEntries());
   hD0PerErr_rej_full->Scale(1./hD0PerErr_rej_full->GetEntries());
   hDZPerErr_rej_full->Scale(1./hDZPerErr_rej_full->GetEntries());


   // projection eta < 1.0
   TH1D *hNhits_acc_narr = (TH1D*) hNhitsEta_acc->ProjectionY("hNhits_acc_narr",eta_min_narr,eta_max_narr,"e");
   TH1D *hPtErr_acc_narr = (TH1D*) hPtErrEta_acc->ProjectionY("hPtErr_acc_narr",eta_min_narr,eta_max_narr,"e");
   TH1D *hD0PerErr_acc_narr  = (TH1D*) hD0PerErrEta_acc->ProjectionY("hD0PerErr_acc_narr",eta_min_narr,eta_max_narr,"e");
   TH1D *hDZPerErr_acc_narr  = (TH1D*) hDZPerErrEta_acc->ProjectionY("hDZPerErr_acc_narr",eta_min_narr,eta_max_narr,"e");

   TH1D *hNhits_rej_narr = (TH1D*) hNhitsEta_rej->ProjectionY("hNhits_rej_narr",eta_min_narr,eta_max_narr,"e");
   TH1D *hPtErr_rej_narr = (TH1D*) hPtErrEta_rej->ProjectionY("hPtErr_rej_narr",eta_min_narr,eta_max_narr,"e");
   TH1D *hD0PerErr_rej_narr  = (TH1D*) hD0PerErrEta_rej->ProjectionY("hD0PerErr_rej_narr",eta_min_narr,eta_max_narr,"e");
   TH1D *hDZPerErr_rej_narr  = (TH1D*) hDZPerErrEta_rej->ProjectionY("hDZPerErr_rej_narr",eta_min_narr,eta_max_narr,"e");



   // dum histogram 
   float xnhitsmin = 8; float xnhitsmax = 28;
   TH1D *hDumNhits = new TH1D("hDumNhits","hDumNhits",100,xnhitsmin,xnhitsmax); hDumNhits->SetStats(0);
   
   float xpterrmin = 0; float xpterrmax = 0.05;
   TH1D *hDumPtErr = new TH1D("hDumPtErr","hDumPtErr",100,xpterrmin,xpterrmax); hDumPtErr->SetStats(0);

   float xnd0min = 0; float xnd0max = 3;
   TH1D *hDumND0 = new TH1D("hDumND0","hDumND0",100,xnd0min,xnd0max); hDumND0->SetStats(0);
   
   float xndzmin = 0; float xndzmax = 3;
   TH1D *hDumNDZ = new TH1D("hDumNDZ","hDumNDZ",100,xndzmin,xndzmax); hDumNDZ->SetStats(0);

   // canvas 1
   TCanvas *c1 = new TCanvas("c1","Track Qualities Distribution 1(|#eta|<2.4)",1000,500);
   c1->Divide(2,1);
   
   TLegend *legNhits = new TLegend(0.56,0.69,0.82,0.90);
   legNhits->SetFillColor(0); legNhits->SetBorderSize(0);
   legNhits->SetHeader("  |#eta|<2.4");
   legNhits->AddEntry(hNhits_acc_full,"Accepted","lp");
   legNhits->AddEntry(hNhits_rej_full,"Rejected","lp");
   
   c1->cd(1); gPad->SetRightMargin(0.12);
   hDumNhits->SetMaximum(0.4);
   hDumNhits->GetXaxis()->SetTitle("Number of valid hits"); hDumNhits->GetYaxis()->SetTitle("A.U. (Normalized to N_{track})");
   hDumNhits->GetXaxis()->CenterTitle(); hDumNhits->GetYaxis()->CenterTitle(); 

   hDumNhits->Draw("");
   hNhits_acc_full->SetMarkerStyle(20); hNhits_acc_full->SetMarkerColor(1); hNhits_acc_full->SetLineColor(1); hNhits_acc_full->Draw("pzsame"); 
   hNhits_rej_full->SetMarkerStyle(24); hNhits_rej_full->SetMarkerColor(2); hNhits_rej_full->SetLineColor(2); hNhits_rej_full->Draw("pzsame");
   hNhits_acc_full->SetMarkerSize(1.2); hNhits_rej_full->SetMarkerSize(1.2);
   legNhits->Draw("");

   c1->cd(2); gPad->SetRightMargin(0.12);
   hDumPtErr->SetMaximum(0.5);
   hDumPtErr->GetXaxis()->SetTitle("#Delta p_{T}/p_{T}"); hDumPtErr->GetYaxis()->SetTitle("A.U. (Normalized to N_{track})");
   hDumPtErr->GetXaxis()->CenterTitle(); hDumPtErr->GetYaxis()->CenterTitle();

   hDumPtErr->Draw("");
   hPtErr_acc_full->SetMarkerStyle(20); hPtErr_acc_full->SetMarkerColor(1); hPtErr_acc_full->SetLineColor(1); hPtErr_acc_full->Draw("pzsame");
   hPtErr_rej_full->SetMarkerStyle(24); hPtErr_rej_full->SetMarkerColor(2); hPtErr_rej_full->SetLineColor(2); hPtErr_rej_full->Draw("pzsame");
   hPtErr_acc_full->SetMarkerSize(1.2); hPtErr_rej_full->SetMarkerSize(1.2);
   legNhits->Draw("");


   // canvas 2 
   TCanvas *c2 = new TCanvas("c2","Track Qualities Distribution 2 (|#eta|<2.4)",1000,500);
   c2->Divide(2,1);

   TLegend *legD0PerErr = new TLegend(0.56,0.69,0.82,0.90);
   legD0PerErr->SetFillColor(0); legD0PerErr->SetBorderSize(0);
   legD0PerErr->SetHeader("  |#eta|<2.4");
   legD0PerErr->AddEntry(hD0PerErr_acc_full,"Accepted","lp");
   legD0PerErr->AddEntry(hD0PerErr_rej_full,"Rejected","lp");

   c2->cd(1); gPad->SetRightMargin(0.12);
   hDumND0->SetMaximum(0.25);
   hDumND0->GetXaxis()->SetTitle("#||{d_{0}/#Delta d_{0}}"); hDumND0->GetYaxis()->SetTitle("A.U. (Normalized to N_{track})");
   hDumND0->GetXaxis()->CenterTitle(); hDumND0->GetYaxis()->CenterTitle();

   hDumND0->Draw("");
   hD0PerErr_acc_full->SetMarkerStyle(20); hD0PerErr_acc_full->SetMarkerColor(1); hD0PerErr_acc_full->SetLineColor(1); hD0PerErr_acc_full->Draw("pzsame");
   hD0PerErr_rej_full->SetMarkerStyle(24); hD0PerErr_rej_full->SetMarkerColor(2); hD0PerErr_rej_full->SetLineColor(2); hD0PerErr_rej_full->Draw("pzsame");
   hD0PerErr_acc_full->SetMarkerSize(1.2); hD0PerErr_rej_full->SetMarkerSize(1.2);
   legD0PerErr->Draw("");

   c2->cd(2); gPad->SetRightMargin(0.12);
   hDumNDZ->SetMaximum(0.16);
   hDumNDZ->GetXaxis()->SetTitle("#||{d_{Z}/#Delta d_{Z}}"); hDumNDZ->GetYaxis()->SetTitle("A.U. (Normalized to N_{track})");
   hDumNDZ->GetXaxis()->CenterTitle(); hDumNDZ->GetYaxis()->CenterTitle();

   hDumNDZ->Draw("");
   hDZPerErr_acc_full->SetMarkerStyle(20); hDZPerErr_acc_full->SetMarkerColor(1); hDZPerErr_acc_full->SetLineColor(1); hDZPerErr_acc_full->Draw("pzsame");
   hDZPerErr_rej_full->SetMarkerStyle(24); hDZPerErr_rej_full->SetMarkerColor(2); hDZPerErr_rej_full->SetLineColor(2); hDZPerErr_rej_full->Draw("pzsame");
   hDZPerErr_acc_full->SetMarkerSize(1.2); hDZPerErr_rej_full->SetMarkerSize(1.2);
   legD0PerErr->Draw("");




   


}



void set_plot_style() {

  // nicer colz plots                                                                                                                                                    

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}
