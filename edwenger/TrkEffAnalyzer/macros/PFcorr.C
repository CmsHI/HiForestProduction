void PFcorr(bool isMC=false) {
  
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.8,"Y");

  TString filename = "../test/pftuple";
  if(isMC) filename += "MC";
  filename += ".root";

  TFile *f = new TFile(filename.Data());
  TNtuple *nt = (TNtuple*) f->Get("pfCandidateAnalyzer/nt");
  
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);

  TString hadron = "tkptmax>10";
  //hadron += " && type==1";
  nt->SetAlias("trksel","(algo<=7 && nhits>=5 && relpterr<=0.05 && abs(nd0)<3 && abs(ndz)<3)");
  TString seltrk = " && trksel";
  TString fakehadron = hadron + " && fake";

  Float_t axisrange=500.;
  if(!isMC) axisrange=350.;

  c1->cd(1);
  nt->Draw(Form("tkptsum:eetsum+hetsum>>h1(50,0,%f,50,0,%f)",axisrange,axisrange),hadron.Data(),"goff");
  h1->SetTitle("PF candidate with p_{T}>10 GeV/c track element; PF calo cluster E_{T} sum (ECAL+HCAL); Track p_{T} sum");
  h1->Draw("colz");
  if(isMC) {
    nt->Draw(Form("tkptsum:eetsum+hetsum>>h1f(50,0,%f,50,0,%f)",axisrange,axisrange),fakehadron.Data(),"goff");
    TH2F* h1 = h1;
    h1f->Divide(h1); h1f->Scale(100.);
    h1f->Draw("boxsame");
  }
  gPad->SetLogz();

  c1->cd(2);
  nt->Draw(Form("tkptsum:eetsum+hetsum>>h2(50,0,%f,50,0,%f)",axisrange,axisrange),(hadron+seltrk).Data(),"goff");
  h2->SetTitle("PF candidate with p_{T}>10 GeV/c track element; PF calo cluster E_{T} sum (ECAL+HCAL); SELECTED Track p_{T} sum");
  h2->Draw("colz");
  if(isMC) {
    nt->Draw(Form("tkptsum:eetsum+hetsum>>h2f(50,0,%f,50,0,%f)",axisrange,axisrange),(fakehadron+seltrk).Data(),"goff");
    h2f->Draw("boxsame");
  }
  gPad->SetLogz();

  c1->cd(3);
  nt->Draw(Form("tkptmax:eetsum+hetsum>>h3(50,0,%f,50,0,%f)",axisrange,axisrange),hadron.Data(),"goff");
  h3->SetTitle("PF candidate with p_{T}>10 GeV/c track element; PF calo cluster E_{T} sum (ECAL+HCAL); Track with MAX p_{T}");
  h3->Draw("colz");
  if(isMC) {
    nt->Draw(Form("tkptmax:eetsum+hetsum>>h3f(50,0,%f,50,0,%f)",axisrange,axisrange),fakehadron.Data(),"goff");
    h3f->Draw("boxsame");
  }
  gPad->SetLogz();

  c1->cd(4);
  nt->Draw(Form("tkptmax:eetsum+hetsum>>h4(50,0,%f,50,0,%f)",axisrange,axisrange),(hadron+seltrk).Data(),"goff");
  h4->SetTitle("PF candidate with p_{T}>10 GeV/c track element; PF calo cluster E_{T} sum (ECAL+HCAL); SELECTED Track with MAX p_{T}");
  h4->Draw("colz");
  if(isMC) {
    nt->Draw(Form("tkptmax:eetsum+hetsum>>h4f(50,0,%f,50,0,%f)",axisrange,axisrange),(fakehadron+seltrk).Data(),"goff");
    h4f->Draw("boxsame");
  }
  gPad->SetLogz();

}
