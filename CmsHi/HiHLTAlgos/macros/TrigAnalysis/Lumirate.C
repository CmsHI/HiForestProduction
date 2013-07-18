//June 30, 2011. Plots the rate of trigger firing as a function of Lumi Section.


#include <iostream>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"

using namespace std;



void Lumirate(TString infile="/net/hidsk0001/d00/scratch/frankma/data/HAZS/hazsv2-run152791to152957-unfiltered-jetoh-v1/merge/all.root") 
{
  TChain * tree = new TChain("icPu5JetAnalyzer/t"); 
  tree->Add(infile);
        tree->AddFriend("hltanalysis/HltTree", infile);
  cout << "Total: " << tree->GetEntries() << endl;

//TCut evtSel("Run==152791&&hiBin<40"); // Another Run Possibility. 
  TCut evtSel("Run==152957&&hiBin<40");  
  
  //Creating Space of available histograms.
  TH1D * hLumi1 = new TH1D("hLumi1",";Lumi Section; dN^{Evt}/dt (Hz)", 50, -100, 1050);
  TH1D * hLumi2 = new TH1D("hLumi2",";Lumi Section; dN^{Evt}/dt (Hz)", 50, -100, 1050);
  TH1D * hLumi3 = new TH1D("hLumi3",";Lumi Section; dN^{Evt}/dt (Hz)", 50, -100, 1050);
  TH1D * hLumi4 = new TH1D("hLumi4",";Lumi Section; dN^{Evt}/dt (Hz)", 50, -100, 1050);
  TH1D * hLumi5 = new TH1D("hLumi5",";Lumi Section; dN^{Evt}/dt (Hz)", 50, -100, 1050);

  //Assinging event selection to histograms.
  tree->Draw("LumiBlock>>hLumi1", evtSel&&"HLT_HIMinBiasHfOrBSC","goff");
  tree->Draw("LumiBlock>>hLumi2", evtSel&&"HLT_HIMinBiasHfOrBSC","goff");
  tree->Draw("LumiBlock>>hLumi3", evtSel&&"HLT_HIJet35U","goff");
  tree->Draw("LumiBlock>>hLumi4", evtSel&&"HLT_HIL2DoubleMu0","goff");
  tree->Draw("LumiBlock>>hLumi5", evtSel&&"HLT_HIPhoton15","goff");

  //Making Graph properties attractive.
  hLumi2->SetLineColor(kRed);
  hLumi3->SetLineColor(kOrange+1);
  hLumi4->SetLineColor(kGreen+2);
  hLumi5->SetLineColor(kBlue);
  hLumi2->SetMarkerColor(kRed);
  hLumi3->SetMarkerColor(kOrange+1);
  hLumi4->SetMarkerColor(kGreen+2);
  hLumi5->SetMarkerColor(kBlue);
  hLumi2->SetMarkerStyle(kOpenCircle);
  hLumi3->SetMarkerStyle(kOpenCircle);
  hLumi4->SetMarkerStyle(kOpenCircle);
  hLumi5->SetMarkerStyle(kOpenCircle);

  //Scaling histogram entries to convert lumi section to rate. 
  hLumi1->Scale(1./hLumi1->GetBinWidth(1));
  hLumi1->Scale(1./23); // Based on the fact that there are 23 seconds in a single Lumisection. 
  hLumi1->Scale(1*3); //Prescale is 3
  hLumi2->Scale(1./hLumi2->GetBinWidth(1));
  hLumi2->Scale(1./23);
  hLumi3->Scale(1./hLumi3->GetBinWidth(1));
  hLumi3->Scale(1./23);
  hLumi4->Scale(1./hLumi4->GetBinWidth(1));
  hLumi4->Scale(1./23);
  hLumi5->Scale(1./hLumi5->GetBinWidth(1));
  hLumi5->Scale(1./23);

  hLumi1->SetAxisRange(1,600,"Y");

  //Draw histograms
  TCanvas * c4 = new TCanvas("c4", "Event Rate for Lumi Section", 800, 750);
  c4->SetLogy();
  hLumi1->SetAxisRange(1,700,"Y");
  hLumi1->Draw("");
  hLumi2->Draw("same");
  hLumi3->Draw("same");
  hLumi4->Draw("same");
  hLumi5->Draw("same");

  //Create Lengend
  TLegend *t = new TLegend(0.595477,0.7915477,0.918342,0.921008);
   t->SetHeader("AllPhysics, All #eta and Centrality");
   t->SetBorderSize(0);
   t->SetFillStyle(0);
   t->AddEntry(hLumi1,"HLT_HIMinBiasHfOrBSC x Prescale","pl");
   t->AddEntry(hLumi2,"HLT_HIMinBiasHfOrBSC","pl");
   t->AddEntry(hLumi3,"HLT_HIJet35U","pl");
   t->AddEntry(hLumi4,"HLT_HIL2DoubleMu0","pl");
   t->AddEntry(hLumi5,"HLT_HIPhoton15","pl");
   t->Draw();

  //hLumi1->Draw("");
  
  }

