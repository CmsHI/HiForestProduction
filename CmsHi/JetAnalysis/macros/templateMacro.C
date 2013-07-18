#include "CommonParameters.h"
#include "TLine.h"
#include <string>
#include "Slices.h"


using namespace std;




void templateMacro(){

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


  nt->SetAlias("w","pt1> -999");
  ntMix->SetAlias("w","weight*(pthat > 50)");

  TCut weight("w");

  TTree* t[2] = {nt,ntMix};

  TH2D* h[2];
  TH1D* hMean[3];
  TH1D* hSigma[3];


  TCanvas* c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);

  for(int i = 0; i < 2; ++i){
    c1->cd(i+1)->SetLogz();
    h[i] = new TH2D(Form("h%d",i),";E_{T}^{HF[#eta > 4]};N_{trk}^{offline}",35,0,70,250,0,250);
    h[i]->GetXaxis()->CenterTitle();
    h[i]->GetYaxis()->CenterTitle();

    t[i]->Draw(Form("ntrk:cent>>%s",h[i]->GetName()),weight*dijet&&deltaPhi&&jetID&&vtx&&noise,"colz");

    fitSlices(h[i]);
    hMean[i] = (TH1D*)gDirectory->Get(Form("%s_1",h[i]->GetName()));
    hSigma[i] = (TH1D*)gDirectory->Get(Form("%s_2",h[i]->GetName()));

  }


}

