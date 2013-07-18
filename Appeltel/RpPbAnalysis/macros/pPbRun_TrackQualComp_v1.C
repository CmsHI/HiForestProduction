// Rylan Conway, Jan 2013
// This macro plots the various track quality variables from the 2013 pPb run and compares them with the 
// data from MC and the pilot run in September 2012.


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

using namespace std;

void drawTrackPlots(TString runFile = "trackAnalysis_testA_v2.root")
{

  TString etaString = "etaFull";

  TFile * f_pPb = TFile::Open(runFile);
  TFile * f_MC = TFile::Open("trackAnalysis_epos_v2.root");
  TFile * f_pilot = TFile::Open("trackAnalysis_pilot202792_v2.root");

  double mSize = 0.85;

  ////////////////////////////////////////////////////////
  // Read in the data from the corresponding data files //
  ////////////////////////////////////////////////////////

  //pPb run 2013
  // int nEvents_pPb = f_pPb->Get("trkAna/events");
  // int nTracks_pPb = f_pPb->Get("trkAna/tracks");
  // int trackseta_pPb = f_pPb->Get("trkAna/trackseta");
  // int nVertices_pPb = f_pPb->Get("trkAna/vertices");

  // TH1D * evtNtrk_pPb = (TH1D*)f_pPb->Get("trkAna/evtNtrk");
  // TH1D * evtNvtx_pPb = (TH1D*)f_pPb->Get("trkAna/evtNvtx");
  // TH1D * vtxX_pPb = (TH1D*)f_pPb->Get("trkAna/vtxX");
  // TH1D * vtxY_pPb = (TH1D*)f_pPb->Get("trkAna/vtxY");
  // TH1D * vtxZ_pPb = (TH1D*)f_pPb->Get("trkAna/vtxZ");
  // TH2D * vtxNtrk2D_pPb = (TH2D*)f_pPb->Get("trkAna/vtxNtrk2D");
  // TH1D * trkNhit_pPb = (TH1D*)f_pPb->Get("trkAna/trkNhit");
  // TH1D * trkPt_pPb = (TH1D*)f_pPb->Get("trkAna/trkPt");
  // //  TH1D * trkPtErr_pPb = (TH1D*)f_MC->Get(Form("trkAna_%s/trkPterr",MCetaString.Data()));
  // TH1D * trkPtErr_pPb = new TH1D();
  // trkPtErr_pPb->SetName("trkPtErr_pPb");

  // TH1D * trkEta_pPb = (TH1D*)f_pPb->Get("trkAna/trkEta");
  // TH1D * trkPhi_pPb = (TH1D*)f_pPb->Get("trkAna/trkPhi");
  // //  TH1D * trkChi2_pPb = (TH1D*)f_MC->Get(Form("trkAna_%s/trkChi2",MCetaString.Data()));
  // TH1D * trkChi2_pPb = new TH1D();
  // trkChi2_pPb->SetName("trkChi2_pPb");

  // TH1D * trkDxyErr_pPb = (TH1D*)f_pPb->Get("trkAna/trkDxyErr");
  // TH1D * trkDzErr_pPb = (TH1D*)f_pPb->Get("trkAna/trkDzErr");
  // TH2D * trkEtaPhi_pPb = (TH2D*)f_pPb->Get("trkAna/trkNhit");

  int nEvents_pPb = f_pPb->Get(Form("trkAna_%s/events",etaString.Data()));
  int nTracks_pPb = f_pPb->Get(Form("trkAna_%s/tracks",etaString.Data()));
  int trackseta_pPb = f_pPb->Get(Form("trkAna_%s/trackseta",etaString.Data()));
  int nVertices_pPb = f_pPb->Get(Form("trkAna_%s/vertices",etaString.Data()));

  TH1D * evtNtrk_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/evtNtrk",etaString.Data()));
  TH1D * evtNvtx_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/evtNvtx",etaString.Data()));
  TH1D * vtxX_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/vtxX",etaString.Data()));
  TH1D * vtxY_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/vtxY",etaString.Data()));
  TH1D * vtxZ_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/vtxZ",etaString.Data()));
  TH2D * vtxNtrk2D_pPb = (TH2D*)f_pPb->Get(Form("trkAna_%s/vtxNtrk2D",etaString.Data()));
  TH1D * trkNhit_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkNhit",etaString.Data()));
  TH1D * trkPt_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkPt",etaString.Data()));
  TH1D * trkPtErr_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkPterr",etaString.Data()));
  TH1D * trkEta_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkEta",etaString.Data()));
  TH1D * trkPhi_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkPhi",etaString.Data()));
  TH1D * trkChi2_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkChi2",etaString.Data()));
  TH1D * trkDxyErr_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkDxyErr",etaString.Data()));
  TH1D * trkDzErr_pPb = (TH1D*)f_pPb->Get(Form("trkAna_%s/trkDzErr",etaString.Data()));
  TH2D * trkEtaPhi_pPb = (TH2D*)f_pPb->Get(Form("trkAna_%s/trkEtaPhi",etaString.Data()));



  //pPb MC run
  int nEvents_MC = f_MC->Get(Form("trkAna_%s/events",etaString.Data()));
  int nTracks_MC = f_MC->Get(Form("trkAna_%s/tracks",etaString.Data()));
  int trackseta_MC = f_MC->Get(Form("trkAna_%s/trackseta",etaString.Data()));
  int nVertices_MC = f_MC->Get(Form("trkAna_%s/vertices",etaString.Data()));

  TH1D * evtNtrk_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/evtNtrk",etaString.Data()));
  TH1D * evtNvtx_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/evtNvtx",etaString.Data()));
  TH1D * vtxX_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/vtxX",etaString.Data()));
  TH1D * vtxY_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/vtxY",etaString.Data()));
  TH1D * vtxZ_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/vtxZ",etaString.Data()));
  TH2D * vtxNtrk2D_MC = (TH2D*)f_MC->Get(Form("trkAna_%s/vtxNtrk2D",etaString.Data()));
  TH1D * trkNhit_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkNhit",etaString.Data()));
  TH1D * trkPt_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkPt",etaString.Data()));
  TH1D * trkPtErr_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkPterr",etaString.Data()));
  TH1D * trkEta_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkEta",etaString.Data()));
  TH1D * trkPhi_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkPhi",etaString.Data()));
  TH1D * trkChi2_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkChi2",etaString.Data()));
  TH1D * trkDxyErr_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkDxyErr",etaString.Data()));
  TH1D * trkDzErr_MC = (TH1D*)f_MC->Get(Form("trkAna_%s/trkDzErr",etaString.Data()));
  TH2D * trkEtaPhi_MC = (TH2D*)f_MC->Get(Form("trkAna_%s/trkEtaPhi",etaString.Data()));


  //pPb pilot run
  int nEvents_pilot = f_pilot->Get(Form("trkAna_%s/events",etaString.Data()));
  int nTracks_pilot = f_pilot->Get(Form("trkAna_%s/tracks",etaString.Data()));
  int trackseta_pilot = f_pilot->Get(Form("trkAna_%s/trackseta",etaString));
  int nVertices_pilot = f_pilot->Get(Form("trkAna_%s/vertices",etaString.Data()));

  TH1D * evtNtrk_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/evtNtrk",etaString.Data()));
  TH1D * evtNvtx_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/evtNvtx",etaString.Data()));
  TH1D * vtxX_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/vtxX",etaString.Data()));
  TH1D * vtxY_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/vtxY",etaString.Data()));
  TH1D * vtxZ_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/vtxZ",etaString.Data()));
  TH2D * vtxNtrk2D_pilot = (TH2D*)f_pilot->Get(Form("trkAna_%s/vtxNtrk2D",etaString.Data()));
  TH1D * trkNhit_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkNhit",etaString.Data()));
  TH1D * trkPt_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkPt",etaString.Data()));
  TH1D * trkPtErr_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkPterr",etaString.Data()));
  TH1D * trkEta_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkEta",etaString.Data()));
  TH1D * trkPhi_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkPhi",etaString.Data()));
  TH1D * trkChi2_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkChi2",etaString.Data()));
  TH1D * trkDxyErr_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkDxyErr",etaString.Data()));
  TH1D * trkDzErr_pilot = (TH1D*)f_pilot->Get(Form("trkAna_%s/trkDzErr",etaString.Data()));
  TH2D * trkEtaPhi_pilot = (TH2D*)f_pilot->Get(Form("trkAna_%s/trkEtaPhi",etaString.Data()));





  cout << endl << "Total # of events is: " << nEvents_pPb << endl;
  cout << "Total # of tracks is: " << nTracks_pPb << endl;
  cout << "Total # of vertices is: " << nVertices_pPb << endl;
  cout << "trackseta is: " << trackseta_pPb << endl << endl;

  cout << endl << "Total # of MC events is: " << nEvents_MC << endl;
  cout << "Total # of tracks is: " << nTracks_MC << endl;
  cout << "Total # of vertices is: " << nVertices_MC << endl;
  cout << "trackseta is: " << trackseta_MC << endl << endl;

  cout << endl << "Total # of pilot events is: " << nEvents_pilot << endl;
  cout << "Total # of tracks is: " << nTracks_pilot << endl;
  cout << "Total # of vertices is: " << nVertices_pilot << endl;
  cout << "trackseta is: " << trackseta_pilot << endl << endl;

  /////////////////////////////////
  // Normalize the distributions //
  /////////////////////////////////

  evtNtrk_pPb->Sumw2();
  evtNtrk_pPb->Scale(1./evtNtrk_pPb->Integral());
  evtNvtx_pPb->Sumw2();
  evtNvtx_pPb->Scale(1./evtNvtx_pPb->Integral());

  vtxX_pPb->Sumw2();
  vtxX_pPb->Scale(1./vtxX_pPb->Integral());
  vtxY_pPb->Sumw2();
  vtxY_pPb->Scale(1./vtxY_pPb->Integral());
  vtxZ_pPb->Sumw2();
  vtxZ_pPb->Scale(1./vtxZ_pPb->Integral());

  trkNhit_pPb->Sumw2();
  trkNhit_pPb->Scale(1./trkNhit_pPb->Integral());
  trkPt_pPb->Sumw2();
  trkPt_pPb->Scale(1./trkPt_pPb->Integral());
  trkPtErr_pPb->Sumw2();
  trkPtErr_pPb->Scale(1./trkPtErr_pPb->Integral());
  trkEta_pPb->Sumw2();
  trkEta_pPb->Scale(1./trkEta_pPb->Integral());
  trkPhi_pPb->Sumw2();
  trkPhi_pPb->Scale(1./trkPhi_pPb->Integral());
  trkChi2_pPb->Sumw2();
  trkChi2_pPb->Scale(1./trkChi2_pPb->Integral());
  trkDxyErr_pPb->Sumw2();
  trkDxyErr_pPb->Scale(1./trkDxyErr_pPb->Integral());
  trkDzErr_pPb->Sumw2();
  trkDzErr_pPb->Scale(1./trkDzErr_pPb->Integral());



  evtNtrk_pilot->Sumw2();
  evtNtrk_pilot->Scale(1./evtNtrk_pilot->Integral());
  evtNvtx_pilot->Sumw2();
  evtNvtx_pilot->Scale(1./evtNvtx_pilot->Integral());

  vtxX_pilot->Sumw2();
  vtxX_pilot->Scale(1./vtxX_pilot->Integral());
  vtxY_pilot->Sumw2();
  vtxY_pilot->Scale(1./vtxY_pilot->Integral());
  vtxZ_pilot->Sumw2();
  vtxZ_pilot->Scale(1./vtxZ_pilot->Integral());

  trkNhit_pilot->Sumw2();
  trkNhit_pilot->Scale(1./trkNhit_pilot->Integral());
  trkPt_pilot->Sumw2();
  trkPt_pilot->Scale(1./trkPt_pilot->Integral());
  trkPtErr_pilot->Sumw2();
  trkPtErr_pilot->Scale(1./trkPtErr_pilot->Integral());
  trkEta_pilot->Sumw2();
  trkEta_pilot->Scale(1./trkEta_pilot->Integral());
  trkPhi_pilot->Sumw2();
  trkPhi_pilot->Scale(1./trkPhi_pilot->Integral());
  trkChi2_pilot->Sumw2();
  trkChi2_pilot->Scale(1./trkChi2_pilot->Integral());
  trkDxyErr_pilot->Sumw2();
  trkDxyErr_pilot->Scale(1./trkDxyErr_pilot->Integral());
  trkDzErr_pilot->Sumw2();
  trkDzErr_pilot->Scale(1./trkDzErr_pilot->Integral());


  
  evtNtrk_MC->Sumw2();
  evtNtrk_MC->Scale(1./evtNtrk_MC->Integral());
  evtNvtx_MC->Sumw2();
  evtNvtx_MC->Scale(1./evtNvtx_MC->Integral());

  vtxX_MC->Sumw2();
  vtxX_MC->Scale(1./vtxX_MC->Integral());
  vtxY_MC->Sumw2();
  vtxY_MC->Scale(1./vtxY_MC->Integral());
  vtxZ_MC->Sumw2();
  vtxZ_MC->Scale(1./vtxZ_MC->Integral());

  trkNhit_MC->Sumw2();
  trkNhit_MC->Scale(1./trkNhit_MC->Integral());
  trkPt_MC->Sumw2();
  trkPt_MC->Scale(1./trkPt_MC->Integral());
  trkPtErr_MC->Sumw2();
  trkPtErr_MC->Scale(1./trkPtErr_MC->Integral());
  trkEta_MC->Sumw2();
  trkEta_MC->Scale(1./trkEta_MC->Integral());
  trkPhi_MC->Sumw2();
  trkPhi_MC->Scale(1./trkPhi_MC->Integral());
  trkChi2_MC->Sumw2();
  trkChi2_MC->Scale(1./trkChi2_MC->Integral());
  trkDxyErr_MC->Sumw2();
  trkDxyErr_MC->Scale(1./trkDxyErr_MC->Integral());
  trkDzErr_MC->Sumw2();
  trkDzErr_MC->Scale(1./trkDzErr_MC->Integral());



  ////////////////////////
  //  Make ratio hists  //
  ////////////////////////

  TH1D* trkPtRatio = (TH1D*)trkPt_pPb->Clone("trkPtRatio");
  trkPtRatio->Divide(trkPt_pilot);
  TH1D* trkPhiRatio = (TH1D*)trkPhi_pPb->Clone("trkPhiRatio");
  trkPhiRatio->Divide(trkPhi_pilot);
  TH1D* trkEtaRatio = (TH1D*)trkEta_pPb->Clone("trkEtaRatio");
  trkEtaRatio->Divide(trkEta_pilot);

  TH1D* trkPtErrRatio = (TH1D*)trkPtErr_pPb->Clone("trkPtErrRatio");
  trkPtErrRatio->Divide(trkPtErr_pilot);
  TH1D* trkDxyErrRatio = (TH1D*)trkDxyErr_pPb->Clone("trkDxyErrRatio");
  trkDxyErrRatio->Divide(trkDxyErr_pilot);
  TH1D* trkDzErrRatio = (TH1D*)trkDzErr_pPb->Clone("trkDzErrRatio");
  trkDzErrRatio->Divide(trkDzErr_pilot);

  TH1D* trkNhitRatio = (TH1D*)trkNhit_pPb->Clone("trkNhitRatio");
  trkNhitRatio->Divide(trkNhit_pilot);
  TH1D* trkChi2Ratio = (TH1D*)trkChi2_pPb->Clone("trkChi2Ratio");
  trkChi2Ratio->Divide(trkChi2_pilot);


  ////////////////////////
  //  Make pretty plots //
  ////////////////////////


  //vertex plot
  TCanvas* c0 = new TCanvas("c0","Vertex Informatin", 850,750);
  c0->Divide(2,2);

  c0->cd(1);
  evtNvtx_pilot->SetXTitle("Number of vertices per event");
  evtNvtx_pilot->SetMarkerSize(mSize);
  evtNvtx_pilot->SetLineWidth(mSize);
  evtNvtx_pPb->SetMarkerStyle(33);
  evtNvtx_pPb->SetMarkerColor(kRed);
  evtNvtx_pPb->SetMarkerSize(mSize);
  evtNvtx_pPb->SetLineWidth(mSize);
  evtNvtx_MC->SetMarkerStyle(21);
  evtNvtx_MC->SetMarkerColor(kBlue);
  evtNvtx_MC->SetMarkerSize(mSize);
  evtNvtx_MC->SetLineWidth(mSize);
  evtNvtx_pilot->Draw("PE");
  evtNvtx_MC->Draw("PSAME");
  evtNvtx_pilot->Draw("PSAME");
  evtNvtx_pPb->Draw("PSAME");

  evtNvtx_pPb->SetXTitle("Number of vertices per event");
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.85);
  leg->AddEntry(evtNvtx_pPb, "Run #XXXX", "p");
  leg->AddEntry(evtNvtx_pilot, "Pilot Run", "p");
  leg->AddEntry(evtNvtx_MC, "EPOS", "p");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw("SAME");

  c0->cd(2);
  //  vtxX_pPb->Draw("");
  vtxX_pilot->SetXTitle("V_{x}");
  vtxX_pilot->SetMarkerSize(mSize);
  vtxX_pilot->SetLineWidth(mSize);
  vtxX_pPb->SetMarkerStyle(33);
  vtxX_pPb->SetMarkerColor(kRed);
  vtxX_pPb->SetMarkerSize(mSize);
  vtxX_pPb->SetLineWidth(mSize);
  vtxX_pPb->SetLineColor(kRed);
  vtxX_MC->SetMarkerStyle(25);
  vtxX_MC->SetMarkerColor(kBlue);
  vtxX_MC->SetMarkerSize(mSize);
  vtxX_MC->SetLineWidth(mSize);
  vtxX_MC->SetLineColor(kBlue);
  vtxX_MC->Draw("HIST");
  vtxX_pilot->Draw("HISTSAME");
  vtxX_MC->Draw("HISTSAME");
  vtxX_pPb->Draw("HISTSAME");
  gPad->Update();
  vtxX_pPb->SetStats(1);

  c0->cd(3);
  vtxY_pilot->SetXTitle("V_{y}");
  vtxY_pilot->SetMarkerSize(mSize);
  vtxY_pilot->SetLineWidth(mSize);
  vtxY_pPb->SetMarkerStyle(33);
  vtxY_pPb->SetMarkerColor(kRed);
  vtxY_pPb->SetMarkerSize(mSize);
  vtxY_pPb->SetLineWidth(mSize);
  vtxY_pPb->SetLineColor(kRed);
  vtxY_MC->SetMarkerStyle(25);
  vtxY_MC->SetMarkerColor(kBlue);
  vtxY_MC->SetMarkerSize(mSize);
  vtxY_MC->SetLineWidth(mSize);
  vtxY_MC->SetLineColor(kBlue);
  vtxY_MC->Draw("HIST");
  vtxY_pilot->Draw("HISTSAME");
  vtxY_pPb->Draw("HISTSAME");

  c0->cd(4);
  vtxZ_pilot->SetXTitle("V_{z}");
  vtxZ_pilot->SetMarkerSize(mSize);
  vtxZ_pilot->SetLineWidth(mSize);
  vtxZ_pPb->SetMarkerStyle(33);
  vtxZ_pPb->SetMarkerColor(kRed);
  vtxZ_pPb->SetMarkerSize(mSize);
  vtxZ_pPb->SetLineWidth(mSize);
  vtxZ_MC->SetMarkerStyle(21);
  vtxZ_MC->SetMarkerColor(kBlue);
  vtxZ_MC->SetMarkerSize(mSize);
  vtxZ_MC->SetLineWidth(mSize);
  vtxZ_pilot->Draw("PE");
  vtxZ_MC->Draw("PESAME");
  vtxZ_pilot->Draw("PESAME");
  vtxZ_pPb->Draw("PESAME");







  //pt, eta, phi
  TCanvas* c1 = new TCanvas("c1","Track Variables", 1100,750);
  c1->Divide(3,2);

  c1->cd(1);
  trkPt_pPb->SetMarkerStyle(33);
  trkPt_pPb->SetMarkerColor(kRed);
  trkPt_pPb->SetMarkerSize(mSize);
  trkPt_pPb->SetLineWidth(mSize);
  trkPt_MC->SetMarkerStyle(21);
  trkPt_MC->SetMarkerColor(kBlue);
  trkPt_MC->SetMarkerSize(mSize);
  trkPt_MC->SetLineWidth(mSize);
  trkPt_pilot->SetMarkerSize(mSize);
  trkPt_pilot->SetLineWidth(mSize);
  trkPt_MC->SetAxisRange(0,0.16,"Y");
  trkPt_MC->Draw("PE");
  trkPt_pilot->Draw("PESAME");
  trkPt_pPb->Draw("PESAME");

  TLegend* leg1 = new TLegend(0.6, 0.7, 0.9, 0.85);
  leg1->AddEntry(trkPt_pPb, "Run #XXXX", "p");
  leg1->AddEntry(trkPt_pilot, "Pilot Run", "p");
  leg1->AddEntry(trkPt_MC, "EPOS", "p");
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->Draw("SAME");


  c1->cd(2);
  trkPhi_pPb->SetMarkerStyle(33);
  trkPhi_pPb->SetMarkerColor(kRed);
  trkPhi_pPb->SetMarkerSize(mSize);
  trkPhi_pPb->SetLineWidth(mSize);
  trkPhi_MC->SetMarkerStyle(21);
  trkPhi_MC->SetMarkerColor(kBlue);
  trkPhi_MC->SetMarkerSize(mSize);
  trkPhi_MC->SetLineWidth(mSize);
  trkPhi_pilot->SetMarkerSize(mSize);
  trkPhi_pilot->SetLineWidth(mSize);
  trkPhi_pPb->SetAxisRange(0.0075,0.0115,"Y");
  trkPhi_pPb->Draw("PE");
  trkPhi_MC->Draw("PESAME");
  trkPhi_pilot->Draw("PESAME");
  trkPhi_pPb->Draw("PESAME");


  c1->cd(3);
  trkEta_pPb->SetMarkerStyle(33);
  trkEta_pPb->SetMarkerColor(kRed);
  trkEta_pPb->SetMarkerSize(mSize);
  trkEta_pPb->SetLineWidth(mSize);
  trkEta_MC->SetMarkerStyle(21);
  trkEta_MC->SetMarkerColor(kBlue);
  trkEta_MC->SetMarkerSize(mSize);
  trkEta_MC->SetLineWidth(mSize);
  trkEta_MC->SetAxisRange(0.0075,trkEta_MC->GetMaximum() + 0.002,"Y");
  trkEta_pilot->SetMarkerSize(mSize);
  trkEta_pilot->SetLineWidth(mSize);
  trkEta_MC->Draw("PE");
  trkEta_pilot->Draw("PESAME");
  trkEta_pPb->Draw("PESAME");


  c1->cd(4);
  trkPtRatio->SetAxisRange(0.8,1.2,"Y");
  trkPtRatio->SetMarkerSize(mSize);
  trkPtRatio->SetLineWidth(mSize);
  trkPtRatio->Draw("HIST P");
  TLegend* legR = new TLegend(0.34, 0.7, 0.75, 0.9);
  legR->AddEntry(trkPtRatio, "Ratio: 2013/pilot", "p");
  legR->SetBorderSize(0);
  legR->SetFillColor(0);
  legR->Draw("PESAME");


  c1->cd(5);
  trkPhiRatio->SetAxisRange(0.8,1.2,"Y");
  trkPhiRatio->SetMarkerSize(mSize);
  trkPhiRatio->SetLineWidth(mSize);
  trkPhiRatio->Draw("PE");


  c1->cd(6);
  trkEtaRatio->SetAxisRange(0.8,1.2,"Y");
  trkEtaRatio->SetMarkerSize(mSize);
  trkEtaRatio->SetLineWidth(mSize);
  trkEtaRatio->Draw("PE");










 //ptErr, DxyErr, DzErr
  TCanvas* c2 = new TCanvas("c2","Track Variables 2", 1100,750);
  c2->Divide(3,2);

  c2->cd(1);
  trkPtErr_pPb->SetMarkerStyle(33);
  trkPtErr_pPb->SetMarkerColor(kRed);
  trkPtErr_pPb->SetMarkerSize(mSize);
  trkPtErr_pPb->SetLineWidth(mSize);
  trkPtErr_MC->SetMarkerStyle(21);
  trkPtErr_MC->SetMarkerColor(kBlue);
  trkPtErr_MC->SetMarkerSize(mSize);
  trkPtErr_MC->SetLineWidth(mSize);
  trkPtErr_pilot->SetMarkerSize(mSize);
  trkPtErr_pilot->SetLineWidth(mSize);
  trkPtErr_MC->SetAxisRange(0,0.16,"Y");
  trkPtErr_pPb->Draw("PE");
  trkPtErr_MC->Draw("PESAME");
  trkPtErr_pilot->Draw("PESAME");
  trkPtErr_pPb->Draw("PESAME");

  TLegend* leg2 = new TLegend(0.6, 0.7, 0.9, 0.85);
  leg2->AddEntry(trkPtErr_pPb, "Run #XXXX", "p");
  leg2->AddEntry(trkPtErr_pilot, "Pilot Run", "p");
  leg2->AddEntry(trkPtErr_MC, "EPOS", "p");
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->Draw("SAME");


  c2->cd(2);
  trkDxyErr_pPb->SetMarkerStyle(33);
  trkDxyErr_pPb->SetMarkerColor(kRed);
  trkDxyErr_pPb->SetMarkerSize(mSize);
  trkDxyErr_pPb->SetLineWidth(mSize);
  trkDxyErr_MC->SetMarkerStyle(21);
  trkDxyErr_MC->SetMarkerColor(kBlue);
  trkDxyErr_MC->SetMarkerSize(mSize);
  trkDxyErr_MC->SetLineWidth(mSize);
  trkDxyErr_pilot->SetMarkerSize(mSize);
  trkDxyErr_pilot->SetLineWidth(mSize);
  //  trkDxyErr_pPb->SetAxisRange(0.0075,0.0115,"Y");
  trkDxyErr_pPb->Draw("PE");
  trkDxyErr_MC->Draw("PESAME");
  trkDxyErr_pilot->Draw("PESAME");
  trkDxyErr_pPb->Draw("PESAME");


  c2->cd(3);
  trkDzErr_pPb->SetMarkerStyle(33);
  trkDzErr_pPb->SetMarkerColor(kRed);
  trkDzErr_pPb->SetMarkerSize(mSize);
  trkDzErr_pPb->SetLineWidth(mSize);
  trkDzErr_MC->SetMarkerStyle(21);
  trkDzErr_MC->SetMarkerColor(kBlue);
  trkDzErr_MC->SetMarkerSize(mSize);
  trkDzErr_MC->SetLineWidth(mSize);
  //  trkDzErr_MC->SetAxisRange(0.0075,trkDzErr_MC->GetMaximum() + 0.002,"Y");
  trkDzErr_pilot->SetMarkerSize(mSize);
  trkDzErr_pilot->SetLineWidth(mSize);
  trkDzErr_MC->Draw("PE");
  trkDzErr_pilot->Draw("PESAME");
  trkDzErr_pPb->Draw("PESAME");


  c2->cd(4);
  trkPtErrRatio->SetAxisRange(0.8,1.2,"Y");
  trkPtErrRatio->SetMarkerSize(mSize);
  trkPtErrRatio->SetLineWidth(mSize);
  trkPtErrRatio->Draw("PE");
  TLegend* legR2 = new TLegend(0.34, 0.7, 0.75, 0.9);
  legR2->AddEntry(trkPtErrRatio, "Ratio: 2013/pilot", "p");
  legR2->SetBorderSize(0);
  legR2->SetFillColor(0);
  legR2->Draw("SAME");


  c2->cd(5);
  trkDxyErrRatio->SetAxisRange(0.8,1.2,"Y");
  trkDxyErrRatio->SetMarkerSize(mSize);
  trkDxyErrRatio->SetLineWidth(mSize);
  trkDxyErrRatio->Draw("PE");


  c2->cd(6);
  trkDzErrRatio->SetAxisRange(0.8,1.2,"Y");
  trkDzErrRatio->SetMarkerSize(mSize);
  trkDzErrRatio->SetLineWidth(mSize);
  trkDzErrRatio->Draw("PE");
  








  //Nhit, Chi2, etaPhi
  TCanvas* c3 = new TCanvas("c3","Track Variables 3", 1100,750);
  c3->Divide(3,2);

  c3->cd(1);
  trkNhit_pPb->SetMarkerSize(mSize);
  trkNhit_pPb->SetMarkerStyle(33);
  trkNhit_pPb->SetMarkerColor(kRed);
  trkNhit_pPb->SetLineWidth(1);
  trkNhit_MC->SetMarkerSize(mSize);
  trkNhit_MC->SetMarkerStyle(21);
  trkNhit_MC->SetMarkerColor(kBlue);
  trkNhit_MC->SetLineWidth(0);
  trkNhit_pilot->SetMarkerSize(mSize);
  trkNhit_pilot->SetLineWidth(0);
  trkNhit_MC->SetAxisRange(0,0.16,"Y");
  trkNhit_pPb->Draw("PE");
  trkNhit_MC->Draw("PESAME");
  trkNhit_pilot->Draw("PESAME");
  trkNhit_pPb->Draw("PESAME");

  TLegend* leg3 = new TLegend(0.6, 0.7, 0.9, 0.85);
  leg3->AddEntry(trkNhit_pPb, "Run #XXXX", "p");
  leg3->AddEntry(trkNhit_pilot, "Pilot Run", "p");
  leg3->AddEntry(trkNhit_MC, "EPOS", "p");
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->Draw("SAME");


  c3->cd(2);
  trkChi2_pPb->SetMarkerStyle(33);
  trkChi2_pPb->SetMarkerColor(kRed);
  trkChi2_pPb->SetMarkerSize(mSize);
  trkChi2_pPb->SetLineWidth(mSize);
  trkChi2_MC->SetMarkerStyle(21);
  trkChi2_MC->SetMarkerColor(kBlue);
  trkChi2_MC->SetMarkerSize(mSize);
  trkChi2_MC->SetLineWidth(mSize);
  trkChi2_pilot->SetLineWidth(mSize);
  trkChi2_pilot->SetMarkerSize(mSize);
  //  trkChi2_pPb->SetAxisRange(0.0075,0.0115,"Y");
  trkChi2_pPb->Draw("PE");
  trkChi2_MC->Draw("PESAME");
  trkChi2_pilot->Draw("PESAME");
  trkChi2_pPb->Draw("PESAME");


  c3->cd(3);
  trkEtaPhi_pPb->Draw("colz");
  trkEtaPhi_pPb->SetXTitle("Run #XXXX    #eta");
  // TLegend* legEtaPhi= new TLegend(0.34, 0.7, 0.65, 0.85);
  // legEtaPhi->AddEntry((TObject*)0,"Run #XXXX","");
  // legEtaPhi->Draw("same");

  c3->cd(4);
  trkNhitRatio->SetAxisRange(0.8,1.2,"Y");
  trkNhitRatio->SetMarkerSize(mSize);
  trkNhitRatio->SetLineWidth(mSize);
  trkNhitRatio->Draw("PE");
  TLegend* legR3 = new TLegend(0.34, 0.7, 0.75, 0.9);
  legR3->AddEntry(trkNhitRatio, "Ratio: 2013/pilot", "p");
  legR3->SetBorderSize(0);
  legR3->SetFillColor(0);
  legR3->Draw("SAME");


  c3->cd(5);
  trkChi2Ratio->SetAxisRange(0.85,1.15,"Y");
  trkChi2Ratio->SetMarkerSize(mSize);
  trkChi2Ratio->SetLineWidth(mSize);
  trkChi2Ratio->Draw("PE");


  c3->cd(6);
  //  trkEtaPhiRatio->SetAxisRange(0.85,1.15,"Y");
  trkEtaPhi_pilot->Draw("colz");
  trkEtaPhi_pilot->SetXTitle("Pilot Run     #eta");
  // TLegend* legEtaPhi2= new TLegend(0.34, 0.7, 0.65, 0.85);
  // legEtaPhi2->AddEntry((TObject*)0,"Pilot Run","");
  // legEtaPhi2->Draw("same");






  return;
}
