#include "CommonParameters.h"
#include "TLine.h"
#include <string>
#include "Slices.h"


using namespace std;


void plotSummary(int iPlot = 0){

  bool showPbPb = 0;
  if(iPlot == 0 || iPlot == 3) showPbPb = 1;


  TH1::SetDefaultSumw2();

  TString infname = "/d101/yetkin/analysis/d0204/ntuple_data_pPb_akPu3PF_forest71_20130204_01.root";
  TString refname = "/d101/yetkin/analysis/d0205/ntuple_data_PbPb_akPu3PF_forest71_20130205_01.root";
  TString mixname = "/d101/yetkin/analysis/d0204/hijing.root";


  infname = "ntuple_data_pPb_akPu3PF_forest71_20130204_01.root";
  refname = "ntuple_data_PbPb_akPu3PF_forest71_20130205_01.root";
  mixname = "hijing.root";


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

  TFile *infReference = new TFile(refname.Data());
  TTree *ntReference;
  TTree *ntevtReference;

  ntReference = (TTree*) infReference->FindObjectAny("ntdijet");
  ntevtReference = (TTree*) infReference->FindObjectAny("ntevt");
  ntReference->AddFriend(ntevtReference);


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


  ntReference->SetAlias("pt1","jtpt1");
  ntReference->SetAlias("pt2","jtpt2");
  ntReference->SetAlias("eta1","jteta1");
  ntReference->SetAlias("eta2","jteta2");
  ntReference->SetAlias("phi1","jtphi1");
  ntReference->SetAlias("phi2","jtphi2");



  nt->SetAlias("cent","hfp");
  ntMix->SetAlias("cent","hfp");


  nt->SetAlias("w","pt1> -999");
  ntReference->SetAlias("w","pt1> -999");

  ntMix->SetAlias("w","weight*(pthat > 50)");

  TCut weight("w");

  TTree* t[3] = {nt,ntMix,ntReference};

  TH2D* h[3];
  TH1D* hMean[3];
  TH1D* hSigma[3];

  string varNames[] = {
    "p_{T,2}/p_{T,1}",
    "#eta_{dijet} = (#eta_{1}+#eta_{2})/2",
    "",
    "#Delta #phi_{1,2}",
    "",

    "","","","","",
    "","","","","",
    "","","","","",

    "",
    "#eta_{1}",
    "#eta_{1}"
  };

  string figures[] = {"imbalance","etaDijet","",    "dphi","",
		      "","",        "", "pudifference", "ntrk",
		      "","pu1","pu2","trkMax1","trkMax2",
		      "","","","","",
		      "","eta1","eta2"
  };


  nt->SetAlias("dphi","acos(cos(phi1-phi2))");
  ntMix->SetAlias("dphi","acos(cos(phi1-phi2))");
  ntReference->SetAlias("dphi","acos(cos(phi1-phi2))");

  if(iPlot == 0){
    nt->SetAlias("var","pt2/pt1");
    ntReference->SetAlias("var","pt2/pt1");
    ntMix->SetAlias("var","pt2/pt1");
  }

  if(iPlot == 1){
    nt->SetAlias("var","(eta1+eta2)/2");
    ntReference->SetAlias("var","(eta1+eta2)/2");
    ntMix->SetAlias("var","(eta1+eta2)/2");
  }

  if(iPlot == 3){
    nt->SetAlias("var","acos(cos(phi1-phi2))");
    ntReference->SetAlias("var","acos(cos(phi1-phi2))");
    ntMix->SetAlias("var","acos(cos(phi1-phi2))");
  }


  int nBins = 12;
  double xBins[] = {0,5,10,15,20,25,30,35,40,45,50,60,70};

  for(int i = 0; i < 3; ++i){
    h[i] = new TH2D(Form("h%d",i),Form(";E_{T}^{HF[#eta > 4]};%s"),nBins,xBins,100,0,1);
    h[i]->GetXaxis()->CenterTitle();
    h[i]->GetYaxis()->CenterTitle();

    if(i == 1) t[i]->Draw(Form("var:hfp>>%s",h[i]->GetName()),weight*dijet&&deltaPhi&&jetID&&vtx,"colz");
    else t[i]->Draw(Form("var:hfp>>%s",h[i]->GetName()),weight*dijet&&deltaPhi&&jetID&&vtx&&noise,"colz");

    fitSlices(h[i]);
    hMean[i] = (TH1D*)gDirectory->Get(Form("%s_1",h[i]->GetName()));
    hSigma[i] = (TH1D*)gDirectory->Get(Form("%s_2",h[i]->GetName()));
  }


  TCanvas* c1 = new TCanvas("c1","",600,600);

  hMean[1]->SetMarkerColor(2);
  hMean[1]->SetLineColor(2);
  hMean[2]->SetMarkerColor(4);
  hMean[2]->SetLineColor(4);
  hMean[2]->SetMarkerStyle(25);

  hMean[1]->SetMinimum(0.5);
  hMean[1]->SetMaximum(0.9);

  hMean[1]->Draw("hist");
  hMean[0]->Draw("same");
  if(showPbPb) hMean[2]->Draw("same");


  const char* date = "20130205";

  string formats[] = {"gif","pdf","eps","C"};

  for(int it = 0; it < 4; ++it){
    c1->Print(Form("./fig/figure_mean_%s_%s%s.%s",figures[iPlot].data(),date,subtract?"_subt":"",formats[it].data()));
  }



}

