#include <TROOT.h>                   // access to gROOT, entry point to ROOT system
#include <TStyle.h>                  // class to handle plot styles in ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "TCut.h"
#include "TSystem.h"

#include "edwenger/TrackSpectraAnalyzer/macros/hlt_eff_simple/CPlot.h"           // helper class for plots
#include "edwenger/TrackSpectraAnalyzer/macros/hlt_eff_simple/HisGroup.h"
#include "edwenger/TrackSpectraAnalyzer/macros/hlt_eff_simple/tgraphTools.C"

#include <map>
#include <cassert>
using namespace std;

void compareTrigSpectra(TString recoSampleName = "PythiaD6T",
    const char * recoFileName = "plots/mcproc0607_v3/all_mcrecotrk/trigSpec/anaspec.root",
    TString genSampleName = "PythiaD6T",
    const char * genFileName = "plots/mcproc0607_v3/all_mcgentrk/trigSpec/anaspec.root",
    TString outdir="plots/mcproc0607_v3/compare")
{
  CPlot::sOutDir = outdir;
  Float_t histJetEtMax = 300;
  Int_t numPtBins=75;
  TFile * recoFile = new TFile(recoFileName);
  cout << recoFileName << endl;
  recoFile->ls();
  TFile * genFile = new TFile(genFileName);
  cout << genFileName << endl;
  genFile->ls();
  TH1::SetDefaultSumw2();

  // === reco ===
  HisGroup psTrigSpec("scTrigSpec");
  psTrigSpec.Add(recoFile,"hSpecMBJet_scaled");
  psTrigSpec.Add(recoFile,"hSpec6U_scaled");
  psTrigSpec.Add(recoFile,"hSpec15U_scaled");
  psTrigSpec.Add(recoFile,"hSpec30U_scaled");
  psTrigSpec.Add(recoFile,"hSpec50U_scaled");
  psTrigSpec.Sum();

  // Find spectra ratio of combined triggered spectra
  // and all MB spectra
  TH1D * hSpecRatio = (TH1D*)psTrigSpec.hSum_->Clone("hSpecRatio");
  TH1D * hSpecMBHF_scaled = (TH1D*)recoFile->FindObjectAny("hSpecMBHF_scaled");
  hSpecRatio->Divide(hSpecMBHF_scaled);

  // Final plots
  // Spectra comparison
  TCanvas * cPSTrigSpec = new TCanvas("cPSTrigSpec","cPSTrigSpec",510,640);
  CPlot cpPSTrigSpec("finalPSTrigSpec","Jet triggered spectra","p_{T}^{trk} [GeV/c]","# Events/GeV");
  cpPSTrigSpec.SetLogy(1);
  cpPSTrigSpec.SetXRange(0,70);
  cpPSTrigSpec.AddHist1D(recoFile,"hSpecMBHF_scaled","MinBias All","E",kBlack,kFullCircle);
  cpPSTrigSpec.AddHist1D(psTrigSpec.hm_["hSpecMBJet_scaled"],"MinBias All and E_{T}^{jet1} < 30GeV","E",kViolet+2);
  cpPSTrigSpec.AddHist1D(psTrigSpec.hm_["hSpec6U_scaled"],"HLT: Jet6^{Raw} and E_{T}^{jet1} = 30-42GeV","E",kAzure+6);
  cpPSTrigSpec.AddHist1D(psTrigSpec.hm_["hSpec15U_scaled"],"HLT: Jet15^{Raw} and E_{T}^{jet1} = 42-78GeV","E",kGreen-3);
  cpPSTrigSpec.AddHist1D(psTrigSpec.hm_["hSpec30U_scaled"],"HLT: Jet30^{Raw} and E_{T}^{jet1} = 78-94GeV","E",kOrange-5);
  cpPSTrigSpec.AddHist1D(psTrigSpec.hm_["hSpec50U_scaled"],"HLT: Jet50^{Raw} and E_{T}^{jet1} > 94GeV","E",kRed-2);
  cpPSTrigSpec.AddHist1D(psTrigSpec.hSum_,"Triggered Classes Combined","hist",kRed,0);
  cpPSTrigSpec.SetLegendHeader(recoSampleName);
  cpPSTrigSpec.SetLegend(0.35,0.66,0.86,0.93);
  cpPSTrigSpec.Draw(cPSTrigSpec,true,"all");

  // Ratio Plot
  TCanvas * cMBTrigSpecRatio = new TCanvas("cMBTrigSpecRatio","cMBTrigSpecRatio",500,500);
  CPlot cpMBTrigSpecRatio("finalMBTrigSpecRatio","Jet triggered spectra","p_{T}^{trk} [GeV/c]","# evt");
  cpMBTrigSpecRatio.SetXRange(0,70);
  cpMBTrigSpecRatio.SetYRange(0,1.1);
  cpMBTrigSpecRatio.AddHist1D(hSpecRatio,"#frac{Triggered Classes Combined}{MB_HF}","E",kViolet+2);
  cpMBTrigSpecRatio.SetLegend(0.18,0.18,0.72,0.36);
  cpMBTrigSpecRatio.Draw(cMBTrigSpecRatio,true,"gif");
  
  // === gen ===
  HisGroup psTrigGenSpec("scTrigGenSpec");
  psTrigGenSpec.Add(genFile,"hSpecMBJet_scaled");
  psTrigGenSpec.Add(genFile,"hSpec6U_scaled");
  psTrigGenSpec.Add(genFile,"hSpec15U_scaled");
  psTrigGenSpec.Add(genFile,"hSpec30U_scaled");
  psTrigGenSpec.Add(genFile,"hSpec50U_scaled");
  psTrigGenSpec.Sum();

  // Find spectra ratio of combined triggered spectra
  // and all MB spectra
  TH1D * hGenSpecRatio = (TH1D*)psTrigGenSpec.hSum_->Clone("hGenSpecRatio");
  TH1D * hGenSpecMBHF_scaled = (TH1D*)genFile->FindObjectAny("hSpecMBHF_scaled");
  hGenSpecRatio->Divide(hGenSpecMBHF_scaled);

  // Final plots
  // Spectra comparison
  TCanvas * cPSTrigGenSpec = new TCanvas("cPSTrigGenSpec","cPSTrigGenSpec",510,640);
  CPlot cpPSTrigGenSpec("finalPSTrigGenSpec","Jet triggered spectra","p_{T}^{gentrk} [GeV/c]","# Events/GeV");
  cpPSTrigGenSpec.SetLogy(1);
  cpPSTrigGenSpec.SetXRange(0,70);
  cpPSTrigGenSpec.AddHist1D(genFile,"hSpecMBHF_scaled","MB_HF","E",kBlack,kFullCircle);
  cpPSTrigGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpecMBJet_scaled"],"MB_HF, E_{T}^{jet1} < 30GeV","E",kViolet+2);
  cpPSTrigGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpec6U_scaled"],"HLT_L1Jet6U, E_{T}^{jet1} = 30-42GeV","E",kAzure+6);
  cpPSTrigGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpec15U_scaled"],"HLT_Jet15U, E_{T}^{jet1} = 42-78GeV","E",kGreen-3);
  cpPSTrigGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpec30U_scaled"],"HLT_Jet30U, E_{T}^{jet1} = 78-94GeV","E",kOrange-5);
  cpPSTrigGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpec50U_scaled"],"HLT_Jet50U, E_{T}^{jet1} > 94GeV","E",kRed-2);
  cpPSTrigGenSpec.AddHist1D(psTrigGenSpec.hSum_,"Triggered Classes Combined","hist",kRed,0);
  cpPSTrigGenSpec.SetLegendHeader(genSampleName);
  cpPSTrigGenSpec.SetLegend(0.35,0.66,0.86,0.93);
  cpPSTrigGenSpec.Draw(cPSTrigGenSpec,true,"all");

  // Ratio Plot
  TCanvas * cMBTrigGenSpecRatio = new TCanvas("cMBTrigGenSpecRatio","cMBTrigGenSpecRatio",500,500);
  CPlot cpMBTrigGenSpecRatio("finalMBTrigGenSpecRatio","Jet triggered spectra","p_{T}^{gentrk} [GeV/c]","# evt");
  cpMBTrigGenSpecRatio.SetXRange(0,70);
  cpMBTrigGenSpecRatio.SetYRange(0,1.1);
  cpMBTrigGenSpecRatio.AddHist1D(hGenSpecRatio,"#frac{Triggered Classes Combined}{MB_HF}","E",kViolet+2);
  cpMBTrigGenSpecRatio.SetLegend(0.18,0.18,0.72,0.36);
  cpMBTrigGenSpecRatio.Draw(cMBTrigGenSpecRatio,true,"gif");

  // === Reco vs Gen Trk comparison ===
  // Spectra comparison
  TCanvas * cPSTrigRecoGenSpec = new TCanvas("cPSTrigRecoGenSpec","cPSTrigRecoGenSpec",510,640);
  CPlot cpPSTrigRecoGenSpec("finalPSTrigRecoGenSpec","Jet triggered spectra","p_{T}^{trk} [GeV/c]","# Events/GeV");
  cpPSTrigRecoGenSpec.SetLogy(1);
  cpPSTrigRecoGenSpec.SetXRange(0,50);
  cpPSTrigRecoGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpecMBJet_scaled"],"MinBias All and E_{T}^{jet1} < 38GeV, GenTrk","E",kViolet+4,kOpenSquare);
  cpPSTrigRecoGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpec15U_scaled"],"HLT: Jet15^{Raw} and E_{T}^{jet1} = 46-70GeV, GenTrk","E",kGreen-1,kOpenSquare);
  //cpPSTrigRecoGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpec30U_scaled"],"HLT: Jet30^{Raw} and E_{T}^{jet1} = 70-94GeV, GenTrk","E",kOrange-3,kOpenSquare);
  cpPSTrigRecoGenSpec.AddHist1D(psTrigGenSpec.hm_["hSpec50U_scaled"],"HLT: Jet50^{Raw} and E_{T}^{jet1} > 94GeV, GenTrk","E",kRed+3,kOpenSquare);
  cpPSTrigRecoGenSpec.AddHist1D(psTrigSpec.hm_["hSpecMBJet_scaled"],"MinBias All and E_{T}^{jet1} < 38GeV, RecoTrk","E",kViolet+2);
  cpPSTrigRecoGenSpec.AddHist1D(psTrigSpec.hm_["hSpec15U_scaled"],"HLT: Jet15^{Raw} and E_{T}^{jet1} = 46-70GeV, RecoTrk","E",kGreen-3);
  //cpPSTrigRecoGenSpec.AddHist1D(psTrigSpec.hm_["hSpec30U_scaled"],"HLT: Jet30^{Raw} E_{T}^{jet1} = 70-94GeV, RecoTrk","E",kOrange-5);
  cpPSTrigRecoGenSpec.AddHist1D(psTrigSpec.hm_["hSpec50U_scaled"],"HLT: Jet50^{Raw} and E_{T}^{jet1} > 94GeV, RecoTrk","E",kRed-2);
  cpPSTrigRecoGenSpec.SetLegendHeader(recoSampleName);
  cpPSTrigRecoGenSpec.SetLegend(0.26,0.70,0.70,0.93);
  cpPSTrigRecoGenSpec.Draw(cPSTrigRecoGenSpec,true,"all");

  // All done, save hists
  gSystem->mkdir(outdir,true);
  TFile * outf = new TFile(Form("%s/comparespec.root",outdir.Data()),"RECREATE");
  psTrigSpec.Save();
  psTrigGenSpec.Save();
}
