#include "weightMix.C"
#include "TProfile.h"
#include "Slices.h"

static int centMode = 1;
static string centName[] = {"N_{trk}^{offline}","HF+","N_{side}","N_{trk}^{low p_{T}}"};

TH1D* getMult(int isample = 2, string var = "xreco"){

  TH1::SetDefaultSumw2();

  TChain* ntdijet = new TChain("ntdijet");
  TChain* ntmatch = new TChain("ntmatch");
  TChain* ntevt = new TChain("ntevt");
  TChain* ntw;

  const char* pPbFile = "/d101/yetkin/analysis/d0204/ntuple_data_pPb_akPu3PF_forest71_20130204_01.root";
  const char* PbPbFile = "/d101/yetkin/analysis/d0204/ntuple_data_PbPb_akPu3PF_forest71_20130204_01.root";
  const char* HydjetFile = "/d101/yetkin/analysis/d0204/hydjet.root";
  const char* HijingFile = "/d101/yetkin/analysis/d0204/hijing.root";


  pPbFile = "ntuple_data_pPb_akPu3PF_forest71_20130204_01.root";
  //  PbPbFile = "ntuple_data_PbPb2011_akPu3PF_forest71_20130204_01.root";
  PbPbFile = "ntuple_data_PbPb_akPu3PF_forest71_20130205_01.root";

  HijingFile = "hijing.root";
  HydjetFile = "hydjet.root";


  if(isample == 0){
    ntdijet->Add(pPbFile);
    ntevt->Add(pPbFile);
  }

  if(isample == 1){
    ntdijet->Add(PbPbFile);
    ntevt->Add(PbPbFile);
  }

  if(isample == 2){
    ntw = new TChain("ntw");
    ntmatch->Add(HijingFile);
    ntdijet->Add(HijingFile);
    ntevt->Add(HijingFile);
    ntw->Add(HijingFile);
  }

  if(isample == 3){
    ntw = new TChain("ntw");
    ntmatch->Add(HydjetFile);
    ntdijet->Add(HydjetFile);
    ntevt->Add(HydjetFile);
    ntw->Add(HydjetFile);
  }

  string varName = "NNNNN";
  if(var == "ngen30") varName = "N_{GenJet}(p_{T}>30)";
  if(var == "njet30") varName = "N_{RecoJet}(p_{T}>30)";
  if(var == "xreco") varName = "#LTp_{T,2}/p_{T,1}#GT";
  if(var == "dphi") varName = "#LT#Delta#phi_{1,2}#GT";

  if(var == "xgen"){
    varName = "True (GenJet) #LTp_{T,2}/p_{T,1}#GT";    
  }

  ntdijet->AddFriend(ntevt);
  if(isample > 1){
    ntdijet->AddFriend(ntw);
    ntdijet->AddFriend(ntmatch);

  }

  TH2D* h[20];
  TH1D* p[20];

  int nmax = 300;
  if(centMode == 1) nmax = 50;
  if(centMode == 2) nmax = 50;

  int Nbins = 50;
  double ymax = 10;

  if(var == "dphi"){
    ymax = pi;
  }
  if(var == "xreco" || var == "xgen"){
    ymax = 1;
  }

  //  if(isample == 2) nmax = 500;

  for(int i = 0; i < 5; ++i){
    h[i] = new TH2D(Form("h%d",i),"",Nbins,0,nmax,100,0,ymax);
  }
  
  TH2D* hn = new TH2D(Form("hn%d",isample),"",Nbins,0,nmax,100,0,2.);

  TCut dijet("jtpt1 > 120 && jtpt2 > 30");

  if(var != "dphi") dijet = dijet&&"acos(cos(jtphi1-jtphi2))> 2.0944";
  TCut weight("weight*(pthat > 80)");
  if(isample < 2) weight = "jtpt1 > -999";

  if(centMode == 0){
    ntdijet->SetAlias("cent","ntrk");
    ntdijet->SetAlias("gens","nps");
    ntdijet->SetAlias("genb","npb");
  }

  if(centMode == 1){
    ntdijet->SetAlias("cent","hfp");
    ntdijet->SetAlias("gens","hfs");
    ntdijet->SetAlias("genb","hfb");
  }

  if(centMode == 2){
    ntdijet->SetAlias("cent","nside");
    ntdijet->SetAlias("gens","nls");
    ntdijet->SetAlias("genb","nlb");
  }

  if(centMode == 3){
    ntdijet->SetAlias("cent","nls+nlb");
    ntdijet->SetAlias("gens","nls");
    ntdijet->SetAlias("genb","nlb");
  }

  ntdijet->SetAlias("gen","gens+genb");


  ntdijet->SetAlias("xgen","pt2/pt1");
  ntdijet->SetAlias("xreco","jtpt2/jtpt1");
  ntdijet->SetAlias("dphi","acos(cos(jtphi1-jtphi2))");
  ntdijet->SetAlias("etaDijet","(jteta1+jteta2)/2");


  ntdijet->SetAlias("var",var.data());


  if(var == "jtpt2/jtpt1") var = "xreco";
  if(var == "pt2/pt1") var = "xgen";
  if(var == "acos(cos(jtphi1-jtphi2))") var = "dphi";


  TCanvas* c2 = new TCanvas("c2","",600,600);
  ntdijet->Draw(Form("var:cent>>h%d",0),weight*dijet,"colz");
  c2->Print(Form("figure_c2_%s_%d_centMode%d.pdf",var.data(),isample,centMode));


  if(isample > 1){

    TCanvas* c1 = new TCanvas("c1","",600,600);
    ntdijet->Draw(Form("gens/genb:cent>>hn%d",isample),weight*dijet,"colz");
    hn->SetTitle(Form(";%s;N_{signal particle}/N_{bkg particle}",centName[centMode].data()));
    hn->GetXaxis()->CenterTitle();
    hn->GetYaxis()->CenterTitle();

    hn->Draw("colz");


    hn->ProfileX(Form("phn_%d",isample))->Draw("same");
    c1->Print(Form("figure_c1_%s_%d_centMode%d.pdf",var.data(),isample,centMode));
    
    TCanvas* c3 = new TCanvas("c3","",600,600);
    ntdijet->Draw(Form("var:genb+gens>>h%d",1),weight*dijet,"colz");
    c3->Print(Form("figure_c3_%s_%d_centMode%d.pdf",var.data(),isample,centMode));

    
    TCanvas* c4 = new TCanvas("c4","",600,600);
    ntdijet->Draw(Form("var:gens>>h%d",2),weight*dijet,"colz");
    c4->Print(Form("figure_c4_%s_%d_centMode%d.pdf",var.data(),isample,centMode));

    TCanvas* c5 = new TCanvas("c5","",600,600);
    ntdijet->Draw(Form("var:genb>>h%d",3),weight*dijet,"colz");
    c5->Print(Form("figure_c5_%s_%d_centMode%d.pdf",var.data(),isample,centMode));
  
  }

  TCanvas* cc1 = new TCanvas("cc1","",600,600);

  for(int i = 0; i < 4; ++i){
    //    h[i]->RebinX(2);

    fitSlices(h[i]);
    p[i] = (TH1D*)gDirectory->Get(Form("%s_1",h[i]->GetName()));
    //    p[i] = h[i]->ProfileX();

    p[i]->SetName(Form("p%d_%d",i,isample));
  }

  p[2]->SetMarkerColor(2);
  p[2]->SetLineColor(2);

  p[3]->SetMarkerColor(4);
  p[3]->SetLineColor(4);

  p[1]->SetMaximum(5);
  p[1]->SetMinimum(0);

  if(var == "xreco" || var == "xgen"){
    p[1]->SetMaximum(1);
  }
  if(var == "dphi"){
    p[1]->SetMaximum(pi);
    p[1]->SetMinimum(pi-pi/3);

  }
  if(var == "nps" || var == "npb"){
    p[1]->SetMaximum(50);
  }

  p[1]->GetXaxis()->SetTitle("N_{?}");
  p[0]->GetXaxis()->SetTitle(centName[centMode].data());

  for(int i = 0; i < 4; ++i){
    p[i]->GetYaxis()->SetTitle(varName.data());
    p[i]->GetXaxis()->CenterTitle();
    p[i]->GetYaxis()->CenterTitle();
  }


  p[1]->Draw("hist");
  p[0]->Draw("same");

  p[2]->Draw("same");
  p[3]->Draw("same");

  TLegend *t3=new TLegend(0.26,0.2,0.53,0.45);
  t3->SetFillColor(0);
  t3->SetBorderSize(0);
  t3->SetFillStyle(0);
  t3->SetTextFont(63);
  t3->SetTextSize(20);
  if(centMode == 0){
    t3->AddEntry(p[1],Form("N = %s","N_{gen}"),"l");
    t3->AddEntry(p[0],Form("N = %s","N_{trk}"),"p");
    t3->AddEntry(p[2],Form("N = %s","N_{signal}"),"p");
    t3->AddEntry(p[3],Form("N = %s","N_{bkg}"),"p");
  }
  if(centMode == 1){
    t3->AddEntry(p[1],Form("N = %s","E_{T}^{HF}_{gen}"),"l");
    t3->AddEntry(p[0],Form("N = %s","E_{T}^{HF}_{reco}"),"p");
    t3->AddEntry(p[2],Form("N = %s","E_{T}^{HF}_{signal}"),"p");
    t3->AddEntry(p[3],Form("N = %s","E_{T}^{HF}_{bkg}"),"p");
  }

  t3->Draw();

  cc1->Print(Form("figure_cc1_%s_%d_centMode%d.pdf",var.data(),isample,centMode));


  return p[0];

}


void plotVariable(string var = "jtpt2/jtpt1"){

   //   var = "nps";
   //   var = "npb";

  TCanvas* cc2 = new TCanvas("cc2","",600,600);
  TFile* outf = new TFile("nothing.root","recreate");
  TH1D* p[4];

  for(int i = 0; i < 4; ++i){
    p[i] = getMult(i,var);
    p[i]->SetName(Form("pp%d",i));
    p[i]->SetDirectory(0);
  }

  p[1]->GetXaxis()->SetTitle(centName[centMode].data());

  p[1]->SetMarkerStyle(25);
  p[1]->SetMarkerColor(4);
  p[1]->SetLineColor(4);

  p[3]->SetLineColor(4);

  p[1]->SetMaximum(5);
  p[1]->SetMinimum(0);
  if(var == "xreco" || var == "xgen"){
    p[1]->SetMaximum(1);
  }
  if(var == "dphi"){
    p[1]->SetMaximum(pi);
    p[1]->SetMinimum(pi-pi/3);

  }
  if(var == "nps" || var == "npb"){
    p[1]->SetMaximum(250);
  }

  cc2->cd();
  p[1]->Draw("");
  p[0]->Draw("same");
  p[2]->Draw("hist same");
  p[3]->Draw("hist same");

  TLegend *t3=new TLegend(0.26,0.2,0.53,0.45);
  t3->SetFillColor(0);
  t3->SetBorderSize(0);
  t3->SetFillStyle(0);
  t3->SetTextFont(63);
  t3->SetTextSize(20);
  t3->AddEntry(p[0],"pPb data","p");
  t3->AddEntry(p[1],"PbPb data","p");
  t3->AddEntry(p[2],"pPb Hijing + Pythia","l");
  t3->AddEntry(p[3],"PbPb Hydjet + Pythia","l");

  t3->Draw();

  if(var == "jtpt2/jtpt1") var = "xreco";
  if(var == "pt2/pt1") var = "xgen";
  if(var == "acos(cos(jtphi1-jtphi2))") var = "dphi";

  cc2->Print(Form("figure_cc0_%s_centMode%d.pdf",var.data(),centMode));


  TCanvas* cc3 = new TCanvas("cc3","",600,600);

  p[0]->GetXaxis()->SetTitle(centName[centMode].data());

  p[0]->Draw("");
  p[2]->Draw("hist same");

  t3=new TLegend(0.26,0.2,0.53,0.45);
  t3->SetFillColor(0);
  t3->SetBorderSize(0);
  t3->SetFillStyle(0);
  t3->SetTextFont(63);
  t3->SetTextSize(20);
  t3->AddEntry(p[0],"pPb data","p");
  t3->AddEntry(p[2],"pPb Hijing + Pythia","l");

  t3->Draw();

  cc3->Print(Form("figure_cc3_pPb_%s_centMode%d.pdf",var.data(),centMode));

  TCanvas* cc4 = new TCanvas("cc4","",600,600);
  if(centMode == 0){
    p[1]->GetXaxis()->SetTitle("N_{trk}^{heavy ion tracking}");
  }
  if(centMode == 1){
    p[1]->GetXaxis()->SetTitle("E_{T}^{HF+}_{heavy ion reco}");
  }

  p[1]->Draw("");
  p[3]->Draw("hist same");

  t3=new TLegend(0.26,0.2,0.53,0.45);
  t3->SetFillColor(0);
  t3->SetBorderSize(0);
  t3->SetFillStyle(0);
  t3->SetTextFont(63);
  t3->SetTextSize(20);
  t3->AddEntry(p[1],"PbPb data","p");
  t3->AddEntry(p[3],"PbPb Hydjet + Pythia","l");

  t3->Draw();

  cc4->Print(Form("figure_cc4_PbPb_%s_centMode%d.pdf",var.data(),centMode));

  outf->Write();
  outf->Close();


}



void plotSignalMultiplicity(){

  string var = "";

  var = "etaDijet";
  plotVariable(var);

  return;

  var = "xreco";
  plotVariable(var);

  var = "xgen";
  plotVariable(var);

  var = "dphi";
  plotVariable(var);

  var = "ngen30";
  plotVariable(var);

  var = "njet30";
  plotVariable(var);

  var = "nps";
  plotVariable(var);
  var = "npb";
  plotVariable(var);


}



