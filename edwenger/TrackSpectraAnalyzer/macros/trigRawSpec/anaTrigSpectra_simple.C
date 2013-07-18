#include <TROOT.h>                   // access to gROOT, entry point to ROOT system
#include <TStyle.h>                  // class to handle plot styles in ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "TCut.h"
#include "TSystem.h"

#include "edwenger/TrackSpectraAnalyzer/macros/trigRawSpec/CPlot.h"           // helper class for plots
#include "edwenger/TrackSpectraAnalyzer/macros/trigRawSpec/HisTGroup.h"
#include "edwenger/TrackSpectraAnalyzer/macros/trigRawSpec/tgraphTools.C"

#include <map>
#include <cassert>
using namespace std;

Double_t GetEffJetPt(TGraph *eff, Double_t minEff)
{
  assert(eff);
  cout << "Begin: " << eff->GetName() << endl;
  Double_t minpt = 0;
  for(Int_t ig = 0; ig<eff->GetN(); ++ig){
    Double_t x,y;
    eff->GetPoint(ig,x,y);
    //cout << "Point " << ig << " x,y: " << x << "," << y << endl;
    if(y>minEff){
      minpt = x;
      cout<<"jet pt with eff > "<<minEff<<" = "<<minpt<< " GeV" << endl;
      break;
    }
  }
  return  minpt;
}

void anaTrigSpectra_simple(TString sampleName="Data",
    const char * inFileName = "plots/MB-C10-PR9-MBskim-v0_p0628_a2/anahlt.root",
    const char * ntFiles = "/net/hibat0003/d00/scratch/frankma/data/MinimumBias/MB-C10-PR9-MBskim-v0_proc0628_trkAnaNoFilter/trkhists_trkAnaSkimAOD_9*.root",
    TString outdir="plots/MB-C10-PR9-MBskim-v0_p0628_a2",
    TString histDir = "trackAna/",
    TString trkNtName = "nt_jettrack",
    Double_t minHltJetEff = 0.99)
{
  CPlot::sOutDir = outdir;
  Float_t histJetEtMax = 300;
  Int_t numPtBins=75;
  TFile * inFile = new TFile(inFileName);
  TH1::SetDefaultSumw2();

  TGraphAsymmErrors *gHltEff_HltJet6U, *gHltEff_HltJet15U, *gHltEff_HltJet30U, *gHltEff_HltJet50U;
  inFile->GetObject("gHltEff_HltJet15U",gHltEff_HltJet15U);
  inFile->GetObject("gHltEff_HltJet30U",gHltEff_HltJet30U);
  inFile->GetObject("gHltEff_HltJet50U",gHltEff_HltJet50U);

  map<TString, TGraphAsymmErrors* > gAEs;
  gAEs["gHltEff_HltJet15U"] = gHltEff_HltJet15U;
  gAEs["gHltEff_HltJet30U"] = gHltEff_HltJet30U;
  gAEs["gHltEff_HltJet50U"] = gHltEff_HltJet50U;
  map<TString, TGraphAsymmErrors*>::iterator ig;
  for (ig=gAEs.begin(); ig != gAEs.end(); ++ig) {
    assert(ig->second);
    cout << ig->first << " " << ig->second->GetName() << endl;
    clearXErrorBar(ig->second);
  }
  // check
  TCanvas * cHltEff = new TCanvas("cHltEff","cHltEff",500,500);
  CPlot cpHltEff("HltEff","Hlt Eff","p_{T}^{corr jet} [GeV/c]","#epsilon_{HLT}");
  cpHltEff.SetXRange(0,histJetEtMax);
  cpHltEff.SetYRange(0,1.1);
  cpHltEff.AddGraph(gAEs["gHltEff_HltJet15U"],"HLT_Jet15U","pz",kGreen-3);
  cpHltEff.AddGraph(gAEs["gHltEff_HltJet30U"],"HLT_Jet30U","pz",kOrange-5);
  cpHltEff.AddGraph(gAEs["gHltEff_HltJet50U"],"HLT_Jet50U","pz",kRed-2);
  cpHltEff.SetLegendHeader(sampleName);
  cpHltEff.SetLegend(0.58,0.21,0.91,0.50);
  cpHltEff.Draw(cHltEff,false);

  // Find Efficiency > 99% jet et
  map<TString, Double_t> effJetPt;
  for (ig=gAEs.begin(); ig != gAEs.end(); ++ig) {
    effJetPt[ig->first] = GetEffJetPt(ig->second,minHltJetEff);
  }
  // check
  if (//effJetPt["gHltEff_HltJet6U"] >= effJetPt["gHltEff_HltJet15U"] ||
      effJetPt["gHltEff_HltJet15U"] >= effJetPt["gHltEff_HltJet30U"] ||
      effJetPt["gHltEff_HltJet30U"] >= effJetPt["gHltEff_HltJet50U"]) {
    cout << "Problem: turn on point of a lower jet trigger is higher, check HLT efficiency" << endl;
    exit(1);
  }


  // ================== Triggered Spectra Ana ==================
  TChain * nt_jettrk = new TChain(histDir+trkNtName,"ntuple: jets and tracks");
  nt_jettrk->Add(ntFiles);
  cout << "Ana: " << histDir << trkNtName << ": " << nt_jettrk->GetEntries() << " tracks" << endl;;
  HisTGroup<TH1D> hgTrigSpec("TrigSpec",100,0,100);
  hgTrigSpec.Add1D("MB");
  hgTrigSpec.Add1D("MBJet");
  hgTrigSpec.Add1D("15U");
  hgTrigSpec.Add1D("30U");
  hgTrigSpec.Add1D("50U");

  TCut selSpecMBJet(Form("mb && jet<%1.1f",effJetPt["gHltEff_HltJet15U"]));
  TCut selSpec15U(Form("jet15 && jet>=%1.1f && jet<%1.1f",effJetPt["gHltEff_HltJet15U"],effJetPt["gHltEff_HltJet30U"]));
  TCut selSpec30U(Form("jet30 && jet>=%1.1f && jet<%1.1f",effJetPt["gHltEff_HltJet30U"],effJetPt["gHltEff_HltJet50U"]));
  TCut selSpec50U(Form("jet50 && jet>=%1.1f",effJetPt["gHltEff_HltJet50U"]));

  cout << "hSpecMBJet sel: " << TString(selSpecMBJet) << endl;
  cout << "hSpec15U sel: " << TString(selSpec15U) << endl;
  cout << "hSpec30U sel: " << TString(selSpec30U) << endl;
  cout << "hSpec50U sel: " << TString(selSpec50U) << endl;

  nt_jettrk->Draw(Form("pt>>%s",hgTrigSpec.GetH("MB")->GetName()),"mb","goff");
  nt_jettrk->Draw(Form("pt>>%s",hgTrigSpec.GetH("MBJet")->GetName()),selSpecMBJet,"goff");
  nt_jettrk->Draw(Form("pt>>%s",hgTrigSpec.GetH("15U")->GetName()),selSpec15U,"goff");
  nt_jettrk->Draw(Form("pt>>%s",hgTrigSpec.GetH("30U")->GetName()),selSpec30U,"goff");
  nt_jettrk->Draw(Form("pt>>%s",hgTrigSpec.GetH("50U")->GetName()),selSpec50U,"goff");

  // Draw Final Plots
  TCanvas * cTrigSpec = new TCanvas("cTrigSpec","cTrigSpec",500,500);
  CPlot cpTrigSpec("TrigSpec","Jet triggered spectra","p_{T}^{trk} [GeV/c]","# evt");
  cpTrigSpec.SetLogy(1);
  cpTrigSpec.SetXRange(0,hgTrigSpec.xmax_);
  cpTrigSpec.AddHist1D(hgTrigSpec.H("MB"),"MB","hist",kBlack);
  cpTrigSpec.AddHist1D(hgTrigSpec.H("MBJet"),"MB+Jet","",kViolet+2);
  cpTrigSpec.AddHist1D(hgTrigSpec.H("15U"),"HLT_Jet15U","",kGreen-3);
  cpTrigSpec.AddHist1D(hgTrigSpec.H("30U"),"HLT_Jet30U","",kOrange-5);
  cpTrigSpec.AddHist1D(hgTrigSpec.H("50U"),"HLT_Jet50U","",kRed-2);
  cpTrigSpec.SetLegendHeader(sampleName);
  cpTrigSpec.Draw(cTrigSpec,false);
  
  // All done, save hists
  gSystem->mkdir(outdir,true);
  TFile * outf = new TFile(Form("%s/anaspec.root",outdir.Data()),"RECREATE");
  hgTrigSpec.Save();
}
