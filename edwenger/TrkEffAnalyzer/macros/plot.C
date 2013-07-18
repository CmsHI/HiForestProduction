void plot() {

  gStyle->SetOptStat(0);

  set_plot_style();

  //TFile *f = new TFile("../test/hists.root");
  TChain *tsim = new TChain("trkEffAnalyzer/simTrackTree");
  TChain *trec = new TChain("trkEffAnalyzer/recTrackTree");
  
  int nfiles = tsim->Add("dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/edwenger/MinBias/MB-S10-357-TRKANASKIM-v1/4addbb78ddff117c63c16d9a48f17d13/ROOTupleMC_HighPurity_1*_*.root");
  
  cout << "added " << nfiles  << " file(s) to SIM chain" << endl;
  
  int nfiles2 = trec->Add("dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/edwenger/MinBias/MB-S10-357-TRKANASKIM-v1/4addbb78ddff117c63c16d9a48f17d13/ROOTupleMC_HighPurity_1*_*.root");
  
  cout << "added " << nfiles2  << " file(s) to RECO chain" << endl;

  TString trksel = "1";
  trksel += " && abs(d0/d0err)<3.0 && abs(dz/dzerr)<3.0";
  trksel += " && (pterr/ptr)<0.05";
  trksel += " && algo<=7";
  trksel += " && hitr>=5";
  trksel += " && abs(dz)<0.2 && abs(d0)<0.2";

  // sim-to-reco hists
  TH2F *hSim = new TH2F("hSim","Primary Charged Sim Tracks; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  tsim->Draw("pts:etas>>hSim","","goff");
  TH2F *hAcc = new TH2F("hAcc","Tracking Acceptance; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  tsim->Draw("pts:etas>>hAcc","acc","goff");
  TH2F *hEff = new TH2F("hEff","Tracking Efficiency; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  tsim->Draw("pts:etas>>hEff",Form("acc && nrec>0 && %s",trksel.Data()),"goff");
  TH2F *hMul = new TH2F("hMul","Multiple Reconstruction; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  tsim->Draw("pts:etas>>hMul",Form("acc && nrec>1 && %s",trksel.Data()),"goff");

  // reco-to-sim hists
  TH2F *hRec = new TH2F("hRec","Reconstructed Tracks; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  trec->Draw("ptr:etar>>hRec",Form("%s",trksel.Data()),"goff");
  TH2F *hFak = new TH2F("hFak","Fake Tracks; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  trec->Draw("ptr:etar>>hFak",Form("nsim==0 && %s",trksel.Data()),"goff");
  TH2F *hSec = new TH2F("hSec","Secondary Tracks; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  trec->Draw("ptr:etar>>hSec",Form("nsim==1 && status<1 && %s",trksel.Data()),"goff");
  //TH2F *hDec = new TH2F("hDec","Weak-decay Tracks; #eta;p_{T} [GeV/c]",60,-3.,3.,80,0.,10.);
  //trec->Draw("ptr:etar>>hDec",Form("parids!=0 && %s",trksel.Data()),"goff");

  // ratio histograms
  TH2F *rAcc = (TH2F*) hAcc->Clone("rAcc");
  TH2F *rEff = (TH2F*) hEff->Clone("rEff");
  TH2F *rMul = (TH2F*) hMul->Clone("rMul");
  TH2F *rFak = (TH2F*) hFak->Clone("rFak");
  TH2F *rSec = (TH2F*) hSec->Clone("rSec");
  //TH2F *rDec = (TH2F*) hDec->Clone("rDec");

  //---------------------------------------------

  // acceptance fraction
  TCanvas *c1 = new TCanvas("c1","Acceptance Fraction",600,500);
  gPad->SetRightMargin(0.15);
  rAcc->Divide(hAcc,hSim,1,1,"B");
  rAcc->SetStats(0);
  rAcc->SetMaximum(1.0); rAcc->SetMinimum(0.0);
  rAcc->Draw("colz");

  // reco efficiency fraction
  TCanvas *c2 = new TCanvas("c2","Reco Efficiency Fraction",600,500);
  gPad->SetRightMargin(0.15);
  rEff->Divide(hEff,hAcc,1,1,"B");
  rEff->SetStats(0);
  rEff->SetMaximum(1.0); rEff->SetMinimum(0.0);
  rEff->Draw("colz");

  // multiple reco fraction
  TCanvas *c3 = new TCanvas("c3","Multiple Reco Fraction",600,500);
  gPad->SetRightMargin(0.15);
  rMul->Divide(hMul,hAcc,1,1,"B");
  rMul->SetStats(0);
  rMul->SetMaximum(1.0); rMul->SetMinimum(0.0);
  rMul->Draw("colz");

  // fake reco fraction
  TCanvas *c4 = new TCanvas("c4","Fake Reco Fraction",600,500);
  gPad->SetRightMargin(0.15);
  rFak->Divide(hFak,hRec,1,1,"B");
  rFak->SetStats(0);
  rFak->SetMaximum(1.0); rFak->SetMinimum(0.0);
  rFak->Draw("colz");

  // secondary reco fraction
  TCanvas *c5 = new TCanvas("c5","Secondary Fraction",600,500);  
  gPad->SetRightMargin(0.15);
  rSec->Divide(hSec,hRec,1,1,"B");
  rSec->SetStats(0);
  rSec->SetMaximum(1.0); rSec->SetMinimum(0.0);
  rSec->Draw("colz");

  //---------------------------------------------

  // projected hists: pt > 0.5 GeV/c
  TH1D* hSimEta = (TH1D*) hSim->ProjectionX("hSimEta",5,80,"e");
  TH1D* hAccEta = (TH1D*) hAcc->ProjectionX("hAccEta",5,80,"e");
  TH1D* hEffEta = (TH1D*) hEff->ProjectionX("hEffEta",5,80,"e");
  TH1D* hMulEta = (TH1D*) hMul->ProjectionX("hMulEta",5,80,"e");
  TH1D* hRecEta = (TH1D*) hRec->ProjectionX("hRecEta",5,80,"e");
  TH1D* hFakEta = (TH1D*) hFak->ProjectionX("hFakEta",5,80,"e");
  TH1D* hSecEta = (TH1D*) hSec->ProjectionX("hSecEta",5,80,"e");
  TH1D* hDumEta = new TH1D("hDumEta","p_{T} > 500 MeV/c;#eta",60,-3.0,3.0); hDumEta->SetMaximum(1.0);
  TH1D* hDumEta2 = (TH1D*) hDumEta->Clone("hDumEta2"); hDumEta2->SetMaximum(0.05); 

  // projected hists: abs(eta) < 2.4
  TH1D* hSimPt  = (TH1D*) hSim->ProjectionY("hSimPt",7,54,"e");
  TH1D* hAccPt  = (TH1D*) hAcc->ProjectionY("hAccPt",7,54,"e");
  TH1D* hEffPt  = (TH1D*) hEff->ProjectionY("hEffPt",7,54,"e");
  TH1D* hMulPt  = (TH1D*) hMul->ProjectionY("hMulPt",7,54,"e");
  TH1D* hRecPt  = (TH1D*) hRec->ProjectionY("hRecPt",7,54,"e");
  TH1D* hFakPt  = (TH1D*) hFak->ProjectionY("hFakPt",7,54,"e");
  TH1D* hSecPt  = (TH1D*) hSec->ProjectionY("hSecPt",7,54,"e");
  TH1D* hDumPt = new TH1D("hDumPt","|#eta| < 2.4;p_{T}",80,0.0,10.0); hDumPt->SetMaximum(1.0);
  TH1D* hDumPt2 = (TH1D*) hDumPt->Clone("hDumPt2"); hDumPt2->SetMaximum(0.05); 

  // Acceptance
  TGraphAsymmErrors *gAccEta = new TGraphAsymmErrors();
  gAccEta->BayesDivide(hAccEta,hSimEta);
  gAccEta->SetMarkerStyle(20);

  TGraphAsymmErrors *gAccPt = new TGraphAsymmErrors();
  gAccPt->BayesDivide(hAccPt,hSimPt);
  gAccPt->SetMarkerStyle(20);

  TCanvas *c6 = new TCanvas("c6","Acceptance Fraction",900,500);
  c6->Divide(2,1);
  c6->cd(1); hDumEta->Draw(); gAccEta->Draw("pz");
  c6->cd(2); hDumPt->Draw(); gAccPt->Draw("pz");

  // Efficiency
  TGraphAsymmErrors *gEffEta = new TGraphAsymmErrors();
  gEffEta->BayesDivide(hEffEta,hAccEta);
  gEffEta->SetMarkerStyle(20);

  TGraphAsymmErrors *gEffPt = new TGraphAsymmErrors();
  gEffPt->BayesDivide(hEffPt,hAccPt);
  gEffPt->SetMarkerStyle(20);

  TCanvas *c7 = new TCanvas("c7","Efficiency Fraction",900,500);
  c7->Divide(2,1);
  c7->cd(1); hDumEta->Draw(); gEffEta->Draw("pz");
  c7->cd(2); hDumPt->Draw(); gEffPt->Draw("pz");

  // Multiple Reco
  TGraphAsymmErrors *gMulEta = new TGraphAsymmErrors();
  gMulEta->BayesDivide(hMulEta,hAccEta);
  gMulEta->SetMarkerStyle(20);

  TGraphAsymmErrors *gMulPt = new TGraphAsymmErrors();
  gMulPt->BayesDivide(hMulPt,hAccPt);
  gMulPt->SetMarkerStyle(20);

  TCanvas *c8 = new TCanvas("c8","Multiple Fraction",900,500);
  c8->Divide(2,1);
  c8->cd(1); hDumEta2->Draw(); gMulEta->Draw("pz");
  c8->cd(2); hDumPt2->Draw(); gMulPt->Draw("pz");

  // Fakes
  TGraphAsymmErrors *gFakEta = new TGraphAsymmErrors();
  gFakEta->BayesDivide(hFakEta,hRecEta);
  gFakEta->SetMarkerStyle(20);

  TGraphAsymmErrors *gFakPt = new TGraphAsymmErrors();
  gFakPt->BayesDivide(hFakPt,hRecPt);
  gFakPt->SetMarkerStyle(20);

  TCanvas *c9 = new TCanvas("c9","Fake Fraction",900,500);
  c9->Divide(2,1);
  c9->cd(1); hDumEta2->Draw(); gFakEta->Draw("pz");
  c9->cd(2); hDumPt2->Draw(); gFakPt->Draw("pz");

  // Secondaries
  TGraphAsymmErrors *gSecEta = new TGraphAsymmErrors();
  gSecEta->BayesDivide(hSecEta,hRecEta);
  gSecEta->SetMarkerStyle(20);

  TGraphAsymmErrors *gSecPt = new TGraphAsymmErrors();
  gSecPt->BayesDivide(hSecPt,hRecPt);
  gSecPt->SetMarkerStyle(20);

  TCanvas *c10 = new TCanvas("c10","Secondary Fraction",900,500);
  c10->Divide(2,1);
  c10->cd(1); hDumEta2->Draw(); gSecEta->Draw("pz");
  c10->cd(2); hDumPt2->Draw(); gSecPt->Draw("pz");

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
