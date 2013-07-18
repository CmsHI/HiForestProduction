//-----------------------------------------------------------------------------
//  Andre S. Yoon, Nov 16 2010 
//
//  Simple macro that plots various distribution/correlation of/between
//  track and calo activity.
//
//  
//-----------------------------------------------------------------------------


void plotCaloTrack() {

  set_plot_style();

  
  TFile *f = new TFile("../root_files/data/trkhistDATA_nov19_trkhist_cAll_v2.root");

  bool isData = true;

  //TString hDir = "pfCandidateAnalyzer_test";
  TString hDir = "pfCandidateAnalyzer";

  TString hName_ecal_2d; TString hName_hcal_2d; TString hName_calo_2d;
  TString hName_ecal_2df; TString hName_hcal_2df; TString hName_calo_2df;

  if(isData){ // for data, no fake
    hName_ecal_2d = "hTrkPtEtaEcalEtSum"; hName_hcal_2d = "hTrkPtEtaHcalEtSum"; hName_calo_2d = "hTrkPtEtaCaloEtSum";
    hName_ecal_2df = "hTrkPtEtaEcalEtSum"; hName_hcal_2df = "hTrkPtEtaHcalEtSum"; hName_calo_2df = "hTrkPtEtaCaloEtSum";
  }else{
    hName_ecal_2d = "hTrkPtEtaEcalEtSum_real"; hName_hcal_2d = "hTrkPtEtaHcalEtSum_real"; hName_calo_2d = "hTrkPtEtaCaloEtSum_real";
    hName_ecal_2df = "hTrkPtEtaEcalEtSum_fake"; hName_hcal_2df = "hTrkPtEtaHcalEtSum_fake"; hName_calo_2df = "hTrkPtEtaCaloEtSum_fake";
  }

  // eta vs pt vs calo
  TH3F *hEtaPtEcal = (TH3F*) f->Get(Form("%s/%s",hDir.Data(),hName_ecal_2d.Data()));
  TH3F *hEtaPtHcal = (TH3F*) f->Get(Form("%s/%s",hDir.Data(),hName_hcal_2d.Data()));
  TH3F *hEtaPtCalo = (TH3F*) f->Get(Form("%s/%s",hDir.Data(),hName_calo_2d.Data()));

  TH3F *hEtaPtEcal_fake = (TH3F*) f->Get(Form("%s/%s",hDir.Data(),hName_ecal_2df.Data()));
  TH3F *hEtaPtHcal_fake = (TH3F*) f->Get(Form("%s/%s",hDir.Data(),hName_hcal_2df.Data()));
  TH3F *hEtaPtCalo_fake = (TH3F*) f->Get(Form("%s/%s",hDir.Data(),hName_calo_2df.Data()));


  // pt vs calo vs nhits

  // variable range
  double small = 0.01;
  
  int eta_min_full = hEtaPtEcal->GetXaxis()->FindBin(-2.4+small);
  int eta_max_full = hEtaPtEcal->GetXaxis()->FindBin(2.4-small);

  int eta_min_narr = hEtaPtEcal->GetXaxis()->FindBin(-1.0+small);
  int eta_max_narr = hEtaPtEcal->GetXaxis()->FindBin(1.0-small);

  //int nhit_min_low = hPtEcalNhit->GetZaxis()->GetFirst();
  //int nhit_max_low = hPtEcalNhit->GetZaxis()->FindBin(15-small);

  //  int nhit_min_high = nhit_min_low + 1;
  //int nhit_max_high = hPtEcalNhit->GetZaxis()->GetLast();


  // projections

  // projection eta < 2.4
  hEtaPtEcal->GetXaxis()->SetRange(eta_min_full,eta_max_full);
  hEtaPtHcal->GetXaxis()->SetRange(eta_min_full,eta_max_full);
  hEtaPtCalo->GetXaxis()->SetRange(eta_min_full,eta_max_full);
  
  TH2D *hEtaPtEcal_full = (TH2D*) hEtaPtEcal->Project3D("zy"); hEtaPtEcal_full->SetName("hEtaPtEcal_full"); 
  TH2D *hEtaPtHcal_full = (TH2D*) hEtaPtHcal->Project3D("zy"); hEtaPtHcal_full->SetName("hEtaPtHcal_full");
  TH2D *hEtaPtCalo_full = (TH2D*) hEtaPtCalo->Project3D("zy"); hEtaPtCalo_full->SetName("hEtaPtCalo_full");

  hEtaPtEcal_fake->GetXaxis()->SetRange(eta_min_full,eta_max_full);
  hEtaPtHcal_fake->GetXaxis()->SetRange(eta_min_full,eta_max_full);
  hEtaPtCalo_fake->GetXaxis()->SetRange(eta_min_full,eta_max_full);

  TH2D *hEtaPtEcal_fake_full = (TH2D*) hEtaPtEcal_fake->Project3D("zy"); hEtaPtEcal_fake_full->SetName("hEtaPtEcal_fake_full");
  TH2D *hEtaPtHcal_fake_full = (TH2D*) hEtaPtHcal_fake->Project3D("zy"); hEtaPtHcal_fake_full->SetName("hEtaPtHcal_fake_full");
  TH2D *hEtaPtCalo_fake_full = (TH2D*) hEtaPtCalo_fake->Project3D("zy"); hEtaPtCalo_fake_full->SetName("hEtaPtCalo_fake_full");

  // projection eta < 1.0
  hEtaPtEcal->GetXaxis()->SetRange(eta_min_narr,eta_max_narr);
  hEtaPtHcal->GetXaxis()->SetRange(eta_min_narr,eta_max_narr);
  hEtaPtCalo->GetXaxis()->SetRange(eta_min_narr,eta_max_narr);

  TH2D *hEtaPtEcal_narr = (TH2D*) hEtaPtEcal->Project3D("zy"); hEtaPtEcal_narr->SetName("hEtaPtEcal_narr");
  TH2D *hEtaPtHcal_narr = (TH2D*) hEtaPtHcal->Project3D("zy"); hEtaPtHcal_narr->SetName("hEtaPtHcal_narr");
  TH2D *hEtaPtCalo_narr = (TH2D*) hEtaPtCalo->Project3D("zy"); hEtaPtCalo_narr->SetName("hEtaPtCalo_narr");

  hEtaPtEcal_fake->GetXaxis()->SetRange(eta_min_narr,eta_max_narr);
  hEtaPtHcal_fake->GetXaxis()->SetRange(eta_min_narr,eta_max_narr);
  hEtaPtCalo_fake->GetXaxis()->SetRange(eta_min_narr,eta_max_narr);

  TH2D *hEtaPtEcal_fake_narr = (TH2D*) hEtaPtEcal_fake->Project3D("zy"); hEtaPtEcal_fake_narr->SetName("hEtaPtEcal_fake_narr");
  TH2D *hEtaPtHcal_fake_narr = (TH2D*) hEtaPtHcal_fake->Project3D("zy"); hEtaPtHcal_fake_narr->SetName("hEtaPtHcal_fake_narr");
  TH2D *hEtaPtCalo_fake_narr = (TH2D*) hEtaPtCalo_fake->Project3D("zy"); hEtaPtCalo_fake_narr->SetName("hEtaPtCalo_fake_narr");


  // projection nhits < 15

  // dum histogram 1
  float xmin=12; float xmax=60;
  float ymin=0; float ymax=80;

  TH1D *hDum1D = new TH1D("hDum1D","hDum1D",100,xmin,xmax); hDum1D->SetStats(0);
  TH2D *hDumEta2D = new TH2D("hDumEta2D","hDumEta2D",100,xmin,xmax,100,ymin,ymax); hDumEta2D->SetStats(0);

  // functions 
  //TF1 *fFakeReject_t1 = new TF1("","4.0+(x>20)*(x-20.)",0,xmax);
  //TF1 *fFakeReject_t1 = new TF1("","4.0 + 5./(2.0+exp(-1.0*(x-20.))) + (x>25)*x",0,xmax);
  //TF1 *fFakeReject_t1 = new TF1("","4.0 + 5./(2.0+exp(-1.0*(x-20.)))",0,xmax); 
  //TF1 *fFakeReject_t1 = new TF1("","3.0 + (x>20)*5./(2.0+exp(x)) + 0.015*(x>30)*pow(x,1.9)",0,xmax);
  //TF1 *fFakeReject_t1 = new TF1("","3.5 + (x<30)*(x>20)*0.45*(x-15) +0.015*(x>30)*pow(x,1.9)",0,xmax); 
  //TF1 *fFakeReject_t1 = new TF1("","3.0 + 0.015*pow(x,1.9)*(1/(1+0.02*exp(0.5*x)))",0,xmax);  
  //TF1 *fFakeReject_t1 = new TF1("","2.5 + (x<20)*(x>15)*1.5*(x-15) +0.02*(x>20)*pow(x,2)",0,xmax);  
  //TF1 *fFakeReject_t1 = new TF1("","2.5 + (x<20)*(x>15)*1.5*(x-15) + (0.02*(x>20)*pow(x,2)) * 1./(exp(0.009*(x-30)))",0,xmax); 
  //TF1 *fFakeReject_t1 = new TF1("","2.5 + (x<20)*(x>15)*1.5*(x-15) + (1.5*(x>20)*pow(x,8/9)) * 1./(exp(0.009*(x-30)))",0,xmax);
  //TF1 *fFakeReject_t1 = new TF1("","(x>14)*(1.2*pow((x-10),8.7/9)) * (2+1./(exp(-1*(x-13))-1))",0,xmax); 
  //TF1 *fFakeReject_t1 = new TF1("","(x>14)*(1.2*pow(abs(x-10),8.7/9)) * (2+1./(exp(-1*(x-13))-1))",0,xmax);
  //TF1 *fFakeReject_t1 = new TF1("","(x>14)*(1.2*(TMath::abs(x-10))^(8.7/9)) * (2+1./(exp(-1*(x-13))-1))",0,xmax); 
  //TF1 *fFakeReject_t1 = new TF1("","(x>14) * (1.2*pow(TMath::Abs(x-10),8.7/9)) * (2+1./(exp(-1*(x-13))-1))",0,xmax);     

  //0~10 %
  TF1 *fFakeReject_t1 = new TF1("","(x>14) * (1.2*pow(TMath::Abs(x-10),8.5/9)) * (2+1./(exp(-1*(x-13))-1))",0,xmax);

  //30~50%
  //TF1 *fFakeReject_t1 = new TF1("","(x>14) * (1.1*pow(TMath::Abs(x-10),8.3/9)) * (2+1./(exp(-1*(x-13))-1))",0,xmax);

  //50~100%
  //TF1 *fFakeReject_t1 = new TF1("","(x>15) * (1.1*pow(TMath::Abs(x-10),8.3/9)) * (2+1./(exp(-1*(x-13))-1))",0,xmax); 
  

  fFakeReject_t1->SetLineColor(2); fFakeReject_t1->SetLineWidth(2.0); fFakeReject_t1->SetLineStyle(9);  

  // canvase setting up
  
  // legend
  TLatex *labelEta24 = new TLatex(0.7*xmax,0.9*ymax,"|#eta| < 2.4");
  TLatex *labelEta10 = new TLatex(0.7*xmax,0.9*ymax,"|#eta| < 1.0");

  // pt vs calo full eta range

  TCanvas *c1 = new TCanvas("c1","p_{T} vs E_{T} (|#eta|<2.4)",1300,500);
  c1->Divide(3,1); 
  hDum1D->SetMaximum(100);
  hDum1D->GetXaxis()->SetTitle("p_{T} (GeV)"); hDum1D->GetYaxis()->SetTitle("ECAL sum E_{T} (GeV)");
  hDum1D->GetXaxis()->CenterTitle(); hDum1D->GetYaxis()->CenterTitle();

  TH1D *hDum1De = (TH1D*) hDum1D->Clone("hDum1De");
  hDum1De->GetYaxis()->SetTitle("ECAL sum E_{T} (GeV)");

  TH1D *hDum1Dh = (TH1D*) hDum1D->Clone("hDum1Dh");
  hDum1Dh->GetYaxis()->SetTitle("HCAL sum E_{T} (GeV)");
  
  TH1D *hDum1Dc = (TH1D*) hDum1D->Clone("hDum1Dc");
  hDum1Dc->GetYaxis()->SetTitle("ECAL+HCAL sum E_{T} (GeV)");

  c1->cd(1); gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1De->Draw(""); labelEta24->Draw(); hEtaPtEcal_full->Draw("COLZsame"); 
  if(!isData) hEtaPtEcal_fake_full->Draw("BOXsame");


  c1->cd(2); gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1Dh->Draw(""); labelEta24->Draw(); hEtaPtHcal_full->Draw("COLZsame"); 
  if(!isData) hEtaPtHcal_fake_full->Draw("BOXsame");

  c1->cd(3); gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1Dc->Draw(""); labelEta24->Draw(); hEtaPtCalo_full->Draw("COLZsame"); 
  if(!isData) hEtaPtCalo_fake_full->Draw("BOXsame");
  
  TCanvas *c1_1 = new TCanvas("c1_1","p_{T} vs E_{T} (|#eta|<2.4)",450,500);
  gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1Dc->Draw(""); labelEta24->Draw(); hEtaPtCalo_full->Draw("COLZsame");
  if(!isData) hEtaPtCalo_fake_full->Draw("BOXsame");
  fFakeReject_t1->Draw("same");


  // pt vs calo |eta|<1.0 

  TCanvas *c2 = new TCanvas("c2","p_{T} vs E_{T} (|#eta|<1.0)",1300,500);
  c2->Divide(3,1);
  hDum1D->SetMaximum(100);
  hDum1D->GetXaxis()->SetTitle("p_{T} (GeV)"); hDum1D->GetYaxis()->SetTitle("E_{T} (GeV)");
  hDum1D->GetXaxis()->CenterTitle(); hDum1D->GetYaxis()->CenterTitle();

  c2->cd(1); gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1De->Draw(""); labelEta10->Draw(); hEtaPtEcal_narr->Draw("COLZsame"); 
  if(!isData) hEtaPtEcal_fake_narr->Draw("BOXsame");

  c2->cd(2); gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1Dh->Draw(""); labelEta10->Draw(); hEtaPtHcal_narr->Draw("COLZsame"); 
  if(!isData) hEtaPtHcal_fake_narr->Draw("BOXsame");

  c2->cd(3); gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1Dc->Draw(""); labelEta10->Draw(); hEtaPtCalo_narr->Draw("COLZsame"); 
  if(!isData)hEtaPtCalo_fake_narr->Draw("BOXsame");

  TCanvas *c2_1 = new TCanvas("c2_1","p_{T} vs E_{T} (|#eta|<1.0)",450,500);
  gPad->SetRightMargin(0.12); gPad->SetLogz();
  hDum1Dc->Draw(""); labelEta10->Draw(); hEtaPtCalo_narr->Draw("COLZsame");
  if(!isData)hEtaPtCalo_fake_narr->Draw("BOXsame");
  fFakeReject_t1->Draw("same");

  // pt vs calo for nhits < 15


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


