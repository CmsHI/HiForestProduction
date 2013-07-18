void plotEvtSelEff() {

  gStyle->SetOptStat(0);

  TFile *f = new TFile("../test/ROOTupleMC_HighPurity.root");
  TH1F *hNSD = (TH1F*) f->Get("preTrgAna/hGenMultNSD");
  TH1F *hNSDtrg = (TH1F*) f->Get("postTrgAna/hGenMultNSD");
  TH1F *hNSDevt = (TH1F*) f->Get("postEvtSelAna/hGenMultNSD");
  TH1F *hNSDvtx = (TH1F*) f->Get("postVtxAna/hGenMultNSD");
  TH1F *hNSDtrkvtx = (TH1F*) f->Get("postTrkVtxAna/hGenMultNSD");

  TGraphAsymmErrors *gNSDtrg = new TGraphAsymmErrors();
  TGraphAsymmErrors *gNSDevt = new TGraphAsymmErrors();
  TGraphAsymmErrors *gNSDvtx = new TGraphAsymmErrors();
  TGraphAsymmErrors *gNSDtrkvtx = new TGraphAsymmErrors();

  gNSDtrg->BayesDivide(hNSDtrg,hNSD);
  gNSDevt->BayesDivide(hNSDevt,hNSD);
  gNSDvtx->BayesDivide(hNSDvtx,hNSD);
  gNSDtrkvtx->BayesDivide(hNSDtrkvtx,hNSD);

  TH1F *dum1 = new TH1F("dum1",";Charged-particle multiplicity;Fraction of events",100,0,100);
  dum1->SetMaximum(0.05);
  dum1->GetYaxis()->SetTitleOffset(1.8);
  TH1F *dum2 = new TH1F("dum2",";Charged-particle multiplicity;Selection efficiency",100,0,60);

  TCanvas *c1 = new TCanvas("c1","Event Selection",900,500);
  c1->Divide(2,1);

  c1->cd(1);
  dum1->Draw();
  hNSDtrkvtx->Sumw2();
  hNSDtrkvtx->Scale(1./hNSDtrkvtx->GetEntries());
  hNSDtrkvtx->SetMarkerStyle(24);
  hNSDtrkvtx->Draw("pzsame");
  hNSDvtx->Sumw2();
  hNSDvtx->Scale(1./hNSDvtx->GetEntries());
  hNSDvtx->SetMarkerStyle(20);
  hNSDvtx->Draw("pzsame");

  c1->cd(2);
  dum2->Draw();
  gNSDtrkvtx->SetMarkerStyle(24);
  gNSDtrkvtx->Draw("pzsame");
  gNSDvtx->SetMarkerStyle(20);
  gNSDvtx->Draw("pzsame");
}
