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

void finalTrigSpectra(TString sampleName = "Data",
    const char * inFileName = "plots/MB-C10-PR9-MBskim-v0_p0628_a2/anaspec.root",
    TString outdir="plots/MB-C10-PR9-MBskim-v0_p0628_a2")
{
  CPlot::sOutDir = outdir;
  CPlot::sPlotStyle = 50;
  Float_t histJetEtMax = 300;
  Int_t numPtBins=75;
  TFile * inFile = new TFile(inFileName);
  TH1::SetDefaultSumw2();

  HisTGroup<TH1D> hgTrigSpec("TrigSpec");
  hgTrigSpec.Add(inFile,"hTrigSpecMBJet","MBJet");
  hgTrigSpec.Add(inFile,"hTrigSpec15U","15U");
  hgTrigSpec.Add(inFile,"hTrigSpec30U","30U");
  hgTrigSpec.Add(inFile,"hTrigSpec50U","50U");
  hgTrigSpec.Sum();

  // Find spectra ratio of combined triggered spectra
  // and all MB spectra
  HisTGroup<TH1D> hgCompSpec("CompSpec");
  hgCompSpec.Add(inFile,"hTrigSpecMB","MB");
  hgCompSpec.Add(hgTrigSpec.R("Sum"),"Comb");
  hgCompSpec.Divide("Comb","MB");

  // Final plots
  // Spectra comparison
  TCanvas * cPSTrigSpec = new TCanvas("cPSTrigSpec","cPSTrigSpec",510,640);
  CPlot cpPSTrigSpec("finalPSTrigSpec","Jet triggered spectra","p_{T}^{trk} [GeV/c]","# Events/GeV");
  cpPSTrigSpec.SetLogy(1);
  cpPSTrigSpec.SetXRange(0,70);
  cpPSTrigSpec.AddHist1D(inFile,"hTrigSpecMB","MinBias All","E",kBlack,kFullCircle);
  cpPSTrigSpec.AddHist1D(hgTrigSpec.H("MBJet"),"MinBias All and E_{T}^{jet1} < 38GeV","E",kViolet+2);
  cpPSTrigSpec.AddHist1D(hgTrigSpec.H("15U"),"HLT: Jet15^{Raw} and E_{T}^{jet1} = 38-62GeV","E",kGreen-3);
  cpPSTrigSpec.AddHist1D(hgTrigSpec.H("30U"),"HLT: Jet30^{Raw} and E_{T}^{jet1} = 62-90GeV","E",kOrange-5);
  cpPSTrigSpec.AddHist1D(hgTrigSpec.H("50U"),"HLT: Jet50^{Raw} and E_{T}^{jet1} > 90GeV","E",kRed-2);
  cpPSTrigSpec.AddHist1D(hgTrigSpec.R("Sum"),"Triggered Classes Combined","hist",kRed,0);
  cpPSTrigSpec.SetLegendHeader(sampleName);
  cpPSTrigSpec.SetLegend(0.35,0.60,0.86,0.87);
  cpPSTrigSpec.Draw(cPSTrigSpec,true,"all");

  // Ratio Plot
  TCanvas * cMBTrigSpecRatio = new TCanvas("cMBTrigSpecRatio","cMBTrigSpecRatio",500,500);
  CPlot cpMBTrigSpecRatio("finalMBTrigSpecRatio","Jet triggered spectra","p_{T}^{trk} [GeV/c]","# evt");
  cpMBTrigSpecRatio.SetXRange(0,70);
  cpMBTrigSpecRatio.SetYRange(0,1.1);
  cpMBTrigSpecRatio.AddHist1D(hgCompSpec.R("CombDivMB"),"#frac{Triggered Classes Combined}{MB_HF}","E",kViolet+2);
  cpMBTrigSpecRatio.SetLegend(0.18,0.18,0.72,0.36);
  cpMBTrigSpecRatio.Draw(cMBTrigSpecRatio,true,"gif");
}
