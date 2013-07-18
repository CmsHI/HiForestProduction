#include <TROOT.h>                   // access to gROOT, entry point to ROOT system
#include <TStyle.h>                  // class to handle plot styles in ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "TCut.h"

#include "CPlot.h"           // helper class for plots
#include "HisTGroup.h"
#include "tgraphTools.C"

#include <map>
#include <cassert>
using namespace std;

Double_t countEvt(const char * inFileName, TString want="All", TString histDir="trackAna/", TString trigHistName="hJet0Pt_HLT_Jet15U")
{
  HisTGroup<TH1F> hgCt("Ct");
  TFile * inFile = new TFile(inFileName);
  hgCt.Add(inFile,histDir+trigHistName,"Trig0");
  hgCt.Add(inFile,histDir+"hNumJets","NJ");
  hgCt.Add(inFile,histDir+"hJet0Pt","JEt");

  // checks
  Double_t numEvt = hgCt.GetH("NJ")->GetEntries();
  Double_t numJets = hgCt.GetH("NJ")->Integral(2,1000);
  Double_t numGt60GeV = hgCt.GetH("JEt")->Integral(121,1000);
  cout << "Total Number of Events: " << (Int_t)numEvt << endl;
  cout << "Total Number of Selected Jets: " << (Int_t)numJets << endl;
  cout << "Total Number of Selected Jets Above 60GeV: " << (Int_t)numGt60GeV << endl;

  // count
  Double_t numGt60GeVTrig = hgCt.GetH("Trig0")->Integral(121,1000);
  inFile->Close();

  Double_t ans=0;
  if (want.Contains("All")) ans = numEvt;
  if (want.Contains("Gt60Trig")) ans = numGt60GeVTrig;
  return ans;
}

void anaHltJetEff_simple(TString sampleName="#sqrt{s} = 7 TeV",
    const char * inFileName = "~frankma/scratch01/data/MinimumBias/MB-C10-PR9-MBskim-v0_proc0710_trkAnaNoFilter/trkhists_trkAnaSkimAOD_*.root",
    TString outdir="plots/MB-C10-PR9-MnJ-v0_p0710_a3",
    const char * mergedFileName = "~frankma/scratch01/data/MinimumBias/MB-C10-PR9-MBskim-v0_proc0710_trkAnaNoFilter/mergeAll/trkhists_sub0.root",
    const char * trigFileName = "~frankma/scratch01/data/MinimumBias/MB-C10-PR9-JMTskim-v0_proc0710_trkAnaNoFilter/trkhists_trkAnaSkimAOD_*.root",
    TString j50FileName = "~frankma/scratch01/data/JetMETTau/JM-R10A-PR4-Jet50Uskim-v0_proc0111_JetVtx/trkhists_*",
    TString histDir = "trackAna/",
    Bool_t useHist = false)
{
  CPlot::sOutDir = outdir;
  Float_t histJetEtMax = 300;
  Int_t numPtBins=75;
  Float_t histJEtBinWidth = histJetEtMax/numPtBins;
  TH1::SetDefaultSumw2();

  HisTGroup<TH1D> hgJet0Et("Jet0Et",numPtBins,0,histJetEtMax);
  HisTGroup<TH1D> hgTrigJet0Et("TrigJet0Et",numPtBins,0,histJetEtMax);
  TH1D * hTrigJet50UJet0Et=0;
  Double_t numSelMBEvt=0, numJet15UGt60MB=0, fracJet15UGt60MB=0, numJet15UGt60Trig=0;
  Double_t numJet15UGt120TrigJet15U=0, numJet50UGt120TrigJet50U=0;
  // === Inputs ===
  // Plot Jet Pt distributions from various triggers
  if (useHist) {
    TFile * inFile = new TFile(inFileName);
    return;
  } else {
    // === HLT Min Bias ===
    TChain * nt_jet = new TChain(histDir+"nt_jet","ntuple: jets");
    nt_jet->Add(inFileName);
    //nt_jet->Print();
    hgJet0Et.Add1D("MB");
    hgJet0Et.Add1D("15U");
    hgJet0Et.Add1D("30U");
    hgJet0Et.Add1D("50U");
    cout << "===== HLT Min Bias =====" << endl;
    cout << inFileName << endl;
    cout << "Tree Analysis on " << nt_jet->GetEntries() << " events" << endl;
    cout << "Tree Analysis - Jets > 60GeV: " << nt_jet->GetEntries("jet>60") << " events" << endl;
    numJet15UGt60MB = nt_jet->GetEntries("jet>60 && jet15");
    cout << "Tree Analysis 15U - Jets > 60GeV: " << numJet15UGt60MB << " events" << endl;
    // Get Histograms
    TCut baseJetSel="mb";
    nt_jet->Draw(Form("jet>>%s",hgJet0Et.GetH("MB")->GetName()),baseJetSel,"goff");
    nt_jet->Draw(Form("jet>>%s",hgJet0Et.GetH("15U")->GetName()),baseJetSel&&"jet15 && jet>10","goff");
    baseJetSel="jet15";
    nt_jet->Draw(Form("jet>>%s",hgJet0Et.GetH("30U")->GetName()),baseJetSel&&"jet30 && jet>20","goff");
    nt_jet->Draw(Form("jet>>%s",hgJet0Et.GetH("50U")->GetName()),baseJetSel&&"jet50 && jet>30","goff");
    // check number
    cout << "Hist 15U check: # of Jets above 60GeV: " << hgJet0Et.GetH("15U")->Integral(60./histJEtBinWidth+1,1000) << endl;
    // Get MB Event Normalization
    cout << " -- Normalization --" << endl;
    numSelMBEvt = countEvt(mergedFileName,"All",histDir);
    fracJet15UGt60MB=numJet15UGt60MB/numSelMBEvt;
    cout << "Tree Analysis 15U - frac Jets > 60GeV: " << fracJet15UGt60MB << " events" << endl;

    // === HLT Jet15U ===
    TChain * nt_trigjet = new TChain(histDir+"nt_jet","ntuple: jets in triggered events");
    nt_trigjet->Add(trigFileName);
    //nt_trigjet->Print();
    hgTrigJet0Et.Add1D("15U");
    hgTrigJet0Et.Add1D("30U");
    hgTrigJet0Et.Add1D("50U");
    cout << endl << "===== HLT Jet15U =====" << endl;
    cout << trigFileName << endl;
    cout << "Tree Analysis on " << nt_trigjet->GetEntries() << " events" << endl;
    numJet15UGt60Trig = nt_trigjet->GetEntries("jet>60 && jet15");
    numJet15UGt120TrigJet15U = nt_trigjet->GetEntries("jet>120 && jet15");
    cout << "Tree Analysis 15U - Jets > 60GeV: " << numJet15UGt60Trig << " events" << endl;
    cout << "Tree Analysis 15U - Jets > 120GeV: " << numJet15UGt120TrigJet15U << " events" << endl;
    // Get Histograms
    TCut baseTrigJetSel="jet15";
    nt_trigjet->Draw(Form("jet>>%s",hgTrigJet0Et.GetH("15U")->GetName()),baseTrigJetSel&&"jet15 && jet>10","goff");
    nt_trigjet->Draw(Form("jet>>%s",hgTrigJet0Et.GetH("30U")->GetName()),baseTrigJetSel&&"jet30 && jet>20","goff");
    nt_trigjet->Draw(Form("jet>>%s",hgTrigJet0Et.GetH("50U")->GetName()),baseTrigJetSel&&"jet50 && jet>30","goff");
    // check number
    cout << "Hist 15U check: # of Jets above 60GeV: " << hgTrigJet0Et.GetH("15U")->Integral(60./histJEtBinWidth+1,1000) << endl;

    // === Jet50U ===
    TChain * djj50 = new TChain("djcalo/djTree");
    djj50->Add(j50FileName);
    hTrigJet50UJet0Et = new TH1D("hTrigJet50UJet0Et","",numPtBins,0,histJetEtMax);
    cout << endl << "===== HLT Jet50U =====" << endl;
    cout << j50FileName << endl;
    cout << "Tree Analysis on " << djj50->GetEntries() << " events" << endl;
    djj50->Project("hTrigJet50UJet0Et","nljet","nljet>50 && hlt[2]");
    numJet50UGt120TrigJet50U = djj50->GetEntries("nljet>120 && hlt[2]");
    cout << "Tree Analysis 50U - Jets > 120GeV: " << numJet50UGt120TrigJet50U << " events" << endl;
  }

  // === Define Output ===
  gSystem->mkdir(outdir.Data(),true);
  TFile * outf = new TFile(Form("%s/anahlt.root",outdir.Data()),"RECREATE");

  // ======= Begin Jet Turnon Ana =======
  cout << endl << "===== Begin Jet Turn-on Ana =====" << endl;
  // Scale Histograms
  HisTGroup<TH1D> hgScJet0Et("ScJet0Et",numPtBins,0,histJetEtMax);
  cout << "MB Trig Data: Normalize # of events with: " << numSelMBEvt << endl;
  hgScJet0Et.Add(hgJet0Et.H("MB"),"MB",1./(numSelMBEvt*hgJet0Et.H("MB")->GetBinWidth(1)));
  hgScJet0Et.Add(hgJet0Et.H("15U"),"15U",1./(numSelMBEvt*hgJet0Et.H("15U")->GetBinWidth(1)));
  hgScJet0Et.Add(hgJet0Et.H("30U"),"30U",1./(numSelMBEvt*hgJet0Et.H("30U")->GetBinWidth(1)));
  hgScJet0Et.Add(hgJet0Et.H("50U"),"50U",1./(numSelMBEvt*hgJet0Et.H("50U")->GetBinWidth(1)));
  // check number
  cout << "15U: Frac of Jets above 60GeV: " << hgScJet0Et.GetH("15U")->Integral(60./histJEtBinWidth+1,1000)*hgScJet0Et.GetH("15U")->GetBinWidth(1) << endl;

  Double_t numSelTrigEvt = numJet15UGt60Trig/fracJet15UGt60MB;
  HisTGroup<TH1D> hgScTrigJet0Et("ScTrigJet0Et",numPtBins,0,histJetEtMax);
  cout << "Jet15U Trig Data: Normalize # of events with: " << numSelTrigEvt << endl;
  hgScTrigJet0Et.Add(hgTrigJet0Et.H("15U"),"15U",1./(numSelTrigEvt*hgTrigJet0Et.H("15U")->GetBinWidth(1)));
  hgScTrigJet0Et.Add(hgTrigJet0Et.H("30U"),"30U",1./(numSelTrigEvt*hgTrigJet0Et.H("30U")->GetBinWidth(1)));
  hgScTrigJet0Et.Add(hgTrigJet0Et.H("50U"),"50U",1./(numSelTrigEvt*hgTrigJet0Et.H("50U")->GetBinWidth(1)));
  // check number
  cout << "15U: Frac of Jets above 60GeV: " << hgScTrigJet0Et.GetH("15U")->Integral(60./histJEtBinWidth+1,1000)*hgScTrigJet0Et.GetH("15U")->GetBinWidth(1) << endl;

  // Jet50U
  hTrigJet50UJet0Et->Scale((numJet15UGt120TrigJet15U/numJet50UGt120TrigJet50U)/numSelTrigEvt/hTrigJet50UJet0Et->GetBinWidth(1));

  // === Check Histograms ===
  TCanvas * cJet0Et = new TCanvas("cJet0Et","cJet0Et",510,510);
  CPlot cpJet0Et("Jet0Et","Jet Et","E_{T}^{corr jet} [GeV/c]","# Events/4 GeV");
  cpJet0Et.SetLogy(1);
  cpJet0Et.SetXRange(0,histJetEtMax);
  cpJet0Et.AddHist1D(hgJet0Et.H("MB"),"MinBias All","E",kViolet);
  cpJet0Et.AddHist1D(hgJet0Et.H("15U"),"HLT: Jet15^{Raw}","E",kGreen-3);
  cpJet0Et.AddHist1D(hgJet0Et.H("30U"),"HLT: Jet30^{Raw}","E",kOrange-5);
  cpJet0Et.AddHist1D(hgJet0Et.H("50U"),"HLT: Jet50^{Raw}","E",kRed-2);
  cpJet0Et.SetLegendHeader(sampleName);
  cpJet0Et.SetLegend(0.58,0.54,0.98,0.82);
  cpJet0Et.Draw(cJet0Et,false);

  TCanvas * cTrigJet0Et = new TCanvas("cTrigJet0Et","cTrigJet0Et",510,510);
  CPlot cpTrigJet0Et("TrigJet0Et","Jet Et","E_{T}^{corr jet} [GeV/c]","# Events/4 GeV");
  cpTrigJet0Et.SetLogy(1);
  cpTrigJet0Et.SetXRange(0,histJetEtMax);
  cpTrigJet0Et.AddHist1D(hgTrigJet0Et.H("15U"),"HLT: Jet15^{Raw}","E",kGreen-3);
  cpTrigJet0Et.AddHist1D(hgTrigJet0Et.H("30U"),"HLT: Jet30^{Raw}","E",kOrange-5);
  cpTrigJet0Et.AddHist1D(hgTrigJet0Et.H("50U"),"HLT: Jet50^{Raw}","E",kRed-2);
  cpTrigJet0Et.SetLegendHeader(sampleName);
  cpTrigJet0Et.SetLegend(0.58,0.54,0.98,0.82);
  cpTrigJet0Et.Draw(cTrigJet0Et,false);

  // ===== Final Plot =====
  // ===Plotting Pad ===
  TCanvas * cJetTurnOn = new TCanvas("cJetTurnOn","cJetTurnOn",510,640);
  cJetTurnOn->cd();
  TPad * pUpper = new TPad("pUpper","pUpper",0,0.34,1,1,0,0,0);
  pUpper->SetBottomMargin(0);
  //pUpper->SetLeftMargin(0.2);
  pUpper->Draw();
  pUpper->cd();
  pUpper->SetNumber(1);
  cJetTurnOn->cd();
  TPad * pLower = new TPad("pLower","pLower",0,0,1,0.34,0,0,0);
  pLower->SetTopMargin(0);
  pLower->SetBottomMargin(0.1*(1./0.34));
  //pLower->SetLeftMargin(0.2);
  pLower->Draw();
  pLower->cd();
  pLower->SetNumber(2);

  // === Begin Ana ===
  CPlot::sPlotStyle = 50;
  pUpper->cd();
  CPlot cpScJet0Et("ScJet0Et","Lead Jet Et","E_{T}^{corr jet} [GeV/c]","#frac{1}{N_{MB}^{Evt}} #frac{dN^{Leading Jet}}{dE_{T}}");
  cpScJet0Et.SetLogy(1);
  cpScJet0Et.SetXRange(0,histJetEtMax);
  cpScJet0Et.SetYRange(1.1e-10,1e-1);
  TH1D * hTemp = (TH1D*)hgScJet0Et.GetH("MB")->Clone("hTemp");
  hTemp->Scale(0);
  hTemp->SetTitle("");
  cpScJet0Et.AddHist1D(hTemp,"MinBias Triggered","E",0,0);
  cpScJet0Et.AddHist1D(hgScJet0Et.H("MB"),"  HLT MinBias","E",kViolet,kOpenSquare,1.1);
  cpScJet0Et.AddHist1D(hgScJet0Et.H("15U"),"  HLT Jet15U","E",kGreen-1,kOpenSquare,1.1);
  //cpScJet0Et.AddHist1D(hgScJet0Et.H("30U"),"MB: HLT Jet30U","E",kOrange-5);
  //cpScJet0Et.AddHist1D(hgScJet0Et.H("50U"),"MB: HLT Jet50U","E",kRed-2);
  cpScJet0Et.AddHist1D(hTemp,"","E",0,0);
  cpScJet0Et.AddHist1D(hTemp,"Jet15U Triggered","E",0,0);
  cpScJet0Et.AddHist1D(hgScTrigJet0Et.H("15U"),"  HLT Jet15U","E",kGreen-3,kFullCircle);
  //cpScJet0Et.AddHist1D(hgScTrigJet0Et.H("30U"),"Triggered: HLT Jet30U","E",kOrange-5,kFullCircle);
  cpScJet0Et.AddHist1D(hgScTrigJet0Et.H("50U"),"  HLT Jet50U","E",kRed-2,kFullCircle);
  cpScJet0Et.AddHist1D(hTemp,"","E",0,0);
  cpScJet0Et.AddHist1D(hTemp,"Jet50U Triggered","E",0,0);
  cpScJet0Et.AddHist1D(hTrigJet50UJet0Et,"  HLT Jet50U","E",kRed,kOpenStar);
  cpScJet0Et.SetLegend(0.518,0.413,0.919,0.847);
  cpScJet0Et.SetLegendStyle(0.045);
  cpScJet0Et.SetAxisLabeling(15,63,18,63,4,2.5);
  cpScJet0Et.Draw(pUpper,false);
  cpScJet0Et.AddText(sampleName,0.33,0.884);
  cpScJet0Et.AddText("#intLdt = 10.2 nb^{-1}",0.285,0.824);

  // === HLT Eff Ana ===
  map<TString, TGraphAsymmErrors* > gAEs;
  gAEs["gHltEff_HltJet15U"] = new TGraphAsymmErrors();
  gAEs["gHltEff_HltJet30U"] = new TGraphAsymmErrors();
  gAEs["gHltEff_HltJet50U"] = new TGraphAsymmErrors();
  gAEs["gTrigHltEff_HltJet30U"] = new TGraphAsymmErrors();
  gAEs["gTrigHltEff_HltJet50U"] = new TGraphAsymmErrors();
  // Calc Hlt Eff
  gAEs["gHltEff_HltJet15U"]->BayesDivide(hgJet0Et.H("15U"),hgJet0Et.H("MB"));
  gAEs["gHltEff_HltJet30U"]->BayesDivide(hgJet0Et.H("30U"),hgJet0Et.H("15U"));
  gAEs["gHltEff_HltJet50U"]->BayesDivide(hgJet0Et.H("50U"),hgJet0Et.H("15U"));
  gAEs["gTrigHltEff_HltJet30U"]->BayesDivide(hgTrigJet0Et.H("30U"),hgTrigJet0Et.H("15U"));
  gAEs["gTrigHltEff_HltJet50U"]->BayesDivide(hgTrigJet0Et.H("50U"),hgTrigJet0Et.H("15U"));
  map<TString, TGraphAsymmErrors*>::iterator ig;
  for (ig=gAEs.begin(); ig != gAEs.end(); ++ig) {
    // check
    assert(ig->second);
    ig->second->SetName(ig->first);
    cout << ig->first << " " << ig->second->GetName() << endl;
    clearXErrorBar(ig->second);
  }

  pLower->cd();
  CPlot::sPlotStyle = 0;
  CPlot cpHltEff("HltEff","Hlt Eff","E_{T}^{corr jet} [GeV/c]","#frac{Trigger}{Lower Trigger}");
  cpHltEff.SetXRange(0,histJetEtMax);
  cpHltEff.SetYRange(0,1.1);
  cpHltEff.AddGraph(gAEs["gHltEff_HltJet15U"],"HLT: Jet15U","pz",kGreen-1,kOpenSquare,1.1);
  //cpHltEff.AddGraph(gAEs["gHltEff_HltJet30U"],"HLT: Jet30U","pz",kOrange-5);
  //cpHltEff.AddGraph(gAEs["gHltEff_HltJet50U"],"HLT: Jet50U","pz",kRed-2);
  //cpHltEff.AddGraph(gAEs["gTrigHltEff_HltJet30U"],"HLT: Jet30U","pz",kOrange-5,kFullCircle);
  cpHltEff.AddGraph(gAEs["gTrigHltEff_HltJet50U"],"HLT: Jet50U","pz",kRed-2,kFullCircle);
  cpHltEff.ShowLegend(0);
  cpHltEff.SetAxisLabeling(15,63,18,63,4,2.5);
  cpHltEff.Draw(pLower,false);

  // === Xcheck comparison ===
  HisTGroup<TH1D> hgCompJetSc("CompJetSc");
  hgCompJetSc.Add(hgScJet0Et.H("15U"),"TrigMB");
  hgCompJetSc.Add(hgScTrigJet0Et.H("15U"),"TrigJet");
  hgCompJetSc.Divide("TrigMB","TrigJet");
  TCanvas * cCompJetSc = new TCanvas("cCompJetSc","cCompJetSc",510,510);
  hgCompJetSc.R("TrigMBDivTrigJet")->Draw();

  // All done, save hists
  cJetTurnOn->Print(Form("%s/cJetTurnOn.gif",outdir.Data()));
  cJetTurnOn->Print(Form("%s/cJetTurnOn.pdf",outdir.Data()));
  hgJet0Et.Save();
  for (map<TString, TGraphAsymmErrors*>::iterator
      ig=gAEs.begin(); ig != gAEs.end(); ++ig) {
    ig->second->Write();
  }
}
