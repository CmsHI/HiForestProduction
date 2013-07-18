#include "CommonParameters.h"
#include "TLine.h"
#include <string>
#include "Slices.h"


using namespace std;


void analyzePuRejection(string var = "etaDijet"){


  TH1::SetDefaultSumw2();


  TString infname = "/d101/yetkin/analysis/d0204/ntuple_data_pPb_akPu3PF_forest71_20130204_01.root";
  TString refname = "/d101/yetkin/analysis/d0204/ntuple_data_PbPb_akPu3PF_forest71_20130204_01.root";
  TString mixname = "/d101/yetkin/analysis/d0204/hijing.root";

  TFile *inf = new TFile(infname.Data());
  TTree *nt =(TTree*)inf->FindObjectAny("ntdijet");
  TTree *ntevt =(TTree*)inf->FindObjectAny("ntevt");
  nt->AddFriend(ntevt);

  TFile *infMix = new TFile(mixname.Data());
  TTree *ntMix =(TTree*)infMix->FindObjectAny("ntdijet");
  TTree *ntevtMix =(TTree*)infMix->FindObjectAny("ntevt");
  TTree *ntw =(TTree*)infMix->FindObjectAny("ntw");
  ntMix->AddFriend(ntw);
  ntMix->AddFriend(ntevtMix);

  nt->SetAlias("pt1","jtpt1");
  nt->SetAlias("pt2","jtpt2");
  nt->SetAlias("eta1","jteta1");
  nt->SetAlias("eta2","jteta2");
  nt->SetAlias("phi1","jtphi1");
  nt->SetAlias("phi2","jtphi2");

  ntMix->SetAlias("pt1","jtpt1");
  ntMix->SetAlias("pt2","jtpt2");
  ntMix->SetAlias("eta1","jteta1");
  ntMix->SetAlias("eta2","jteta2");
  ntMix->SetAlias("phi1","jtphi1");
  ntMix->SetAlias("phi2","jtphi2");

  nt->SetAlias("cent","hfp");
  ntMix->SetAlias("cent","hfp");

  nt->SetAlias("etaDijet","(eta1+eta2)/2");
  ntMix->SetAlias("etaDijet","(eta1+eta2)/2");

  nt->SetAlias("x","pt2/pt1");
  ntMix->SetAlias("x","pt2/pt1");

  nt->SetAlias("dphi","acos(cos(phi1-phi2))");
  ntMix->SetAlias("dphi","acos(cos(phi1-phi2))");


  nt->SetAlias("var",var.data());
  ntMix->SetAlias("var",var.data());

  nt->SetAlias("w","pt1> -999");
  ntMix->SetAlias("w","weight*(pthat > 50)");

  TCut noises[2] = {"noise == -1","noise == 1"};

  TTree* t[2] = {nt,ntMix};

  TH2D* h[2];
  TH1D* hMean[3];
  TH1D* hSigma[3];

  double yMax = 250;
  double yMin = 0;

  string varName = "N_{trk}^{offline}";

  if(var == "etaDijet"){
    yMax = 2;
    yMin = -2;
    varName = "#eta_{dijet}=(#eta_{1}+#eta_{2})/2";
  }

  if(var == "x"){
    yMax = 1;
    yMin = 0;
    varName = "p_{T,2}/p_{T,1}";
  }

  if(var == "dphi"){
    yMax = pi;
    yMin = pi-pi/3;
    varName = "#Delta#phi_{1,2}";
  }


  TCanvas* c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);

  for(int i = 0; i < 2; ++i){
    c1->cd(i+1)->SetLogz();
    h[i] = new TH2D(Form("h%d",i),Form(";E_{T}^{HF[#eta > 4]};%s",varName.data()),14,0,70,250,yMin,yMax);
    h[i]->GetXaxis()->CenterTitle();
    h[i]->GetYaxis()->CenterTitle();

    t[0]->Draw(Form("var:cent>>%s",h[i]->GetName()),dijet&&deltaPhi&&jetID&&vtx&&noises[i],"colz");
 
    fitSlices(h[i]);

    hMean[i] = (TH1D*)gDirectory->Get(Form("%s_1",h[i]->GetName()));
    hSigma[i] = (TH1D*)gDirectory->Get(Form("%s_2",h[i]->GetName()));
    hMean[i]->SetTitle("");
    hSigma[i]->SetTitle("");

    hMean[i]->SetLineColor(i+1);
    hSigma[i]->SetLineColor(i+1);
    hMean[i]->SetMarkerColor(i+1);
    hSigma[i]->SetMarkerColor(i+1);
    hMean[i]->GetXaxis()->CenterTitle();
    hMean[i]->GetYaxis()->CenterTitle();
    hSigma[i]->GetXaxis()->CenterTitle();
    hSigma[i]->GetYaxis()->CenterTitle();

    hMean[i]->GetYaxis()->SetTitle(Form("#LT%s#GT",varName.data()));
    hSigma[i]->GetYaxis()->SetTitle(Form("#sigma(%s)",varName.data()));


 }

  c1->Print(Form("figure_%s_HF_correlations_dijets_pu.gif",var.data()));
  c1->Print(Form("figure_%s_HF_correlations_dijets_pu.pdf",var.data()));


  TCanvas* c2 = new TCanvas("c2","",600,600);

  if(var == "etaDijet"){
    hMean[0]->SetMinimum(-0.3);
    hMean[0]->SetMaximum(0.3);
  }

  if(var == "x"){
    hMean[0]->SetMinimum(0.6);
    hMean[0]->SetMaximum(0.9);
  }


  hMean[0]->Draw("hist");
  hMean[1]->Draw("same");

  hMean[2] = (TH1D*)hMean[1]->Clone("hr0");
  hMean[2]->Divide(hMean[0]);

  TLegend *t3=new TLegend(0.22,0.72,0.49,0.88);
  t3->SetFillColor(0);
  t3->SetBorderSize(0);
  t3->SetFillStyle(0);
  t3->SetTextFont(63);
  t3->SetTextSize(20);
  t3->AddEntry(hMean[0],Form("pVertexFilterCutGplus"),"l");
  t3->AddEntry(hMean[1],Form("!pVertexFilterCutGplus"),"p");
  t3->Draw();
  c2->Print(Form("figure_%s_HF_mean_dijets_pu.gif",var.data()));
  c2->Print(Form("figure_%s_HF_mean_dijets_pu.pdf",var.data()));


  TCanvas* cr2 = new TCanvas("cr2","",600,600);
  TLine* line = new TLine(0,1,70,1);
  line->SetLineStyle(2);
  hMean[2]->SetMinimum(0);
  hMean[2]->SetMaximum(2);

  hMean[2]->Draw("");
  line->Draw("same");
  t3->Draw();

  cr2->Print(Form("figure_%s_HF_ratio_dijets_pu.gif",var.data()));
  cr2->Print(Form("figure_%s_HF_ratio_dijets_pu.pdf",var.data()));


  TCanvas* c3 = new TCanvas("c3","",600,600);

  hSigma[1]->SetMinimum(0);
  hSigma[1]->SetMaximum(120);

  hSigma[0]->Draw("hist");
  hSigma[1]->Draw("same");
  t3->Draw();


  c3->Print(Form("figure_%s_HF_sigma_dijets_pu.gif",var.data()));
  c3->Print(Form("figure_%s_HF_sigma_dijets_pu.pdf",var.data()));





  



}

