//July 11, 2011: Gives Accept Fraction as a function of uncorrected and corrected Jet pt. Macro is 
//modified to give interaction Rate as a function of pt.

#include "TChain.h"
#include "TH1.h"
#include "TCut.h"
#include "TCanvas.h"
#include <iostream>


void RateEstimateVpt(TString infile="alldatacleaned.root") //Contains both cleaned and uncleaned data
{
  //TChain * tree = new TChain("icPu5JetAnalyzer_clean/t");
  TChain * tree = new TChain("icPu5JetAnalyzer/t");
  tree->Add(infile);
  //tree->AddFriend("hltanalysis_clean/HltTree",infile);
  tree->AddFriend("hltanalysis/HltTree",infile);
  //  cout << "Total: " << tree->GetEntries() << endl; 

  TCut evtSel("Run==152957&&HLT_HIMinBiasHfOrBSC&&hiBin<40&&abs(jteta[0])<2");


  //Canvas definition and choosing uncorrected and corrected jet pt variables 
  TCanvas *e1 = new TCanvas("e1","Rate Estimate vs. Threshold",200,10,700,500);
  e1->SetLogy();
  TH1D * RJetPt = new TH1D("RJetPt",";Leading Jet p_{T} (GeV/c);count;",300,0,300);
  TH1D * JetPt = new TH1D("JetPt",";Leading Jet p_{T} (GeV/c);count;",300,0,300);
  tree->Draw("rawpt[0]>>RJetPt",evtSel,"goff");
  tree->Draw("jtpt[0]>>JetPt",evtSel,"goff");
  cout << "Integral: " << JetPt->Integral(0,101) << endl;

  //Defines total # of entries for normalization and max rate
  TH1F * h4 = new TH1F("h4","",50,0,300);
  tree->Draw("LumiBlock>>h4",evtSel,"goff");
  Float_t R = 1 ; //Defines Maximum Rate
  Float_t TA = h4->GetEntries();

  //Computes Accept Fraction as a function of pt 
  const Int_t n = 61;
  Float_t x[n], y[n], z[n], w[n], t[n];
  for (Int_t i=0;i<n;i++) {
    x[i] = 300-i*5; 
    y[i] = RJetPt->Integral(x[i],301);
    z[i] = R*y[i]/TA;
    w[i] = JetPt->Integral(x[i],301);
    t[i] = R*w[i]/TA;
    printf(" i %i %f %f %f \n",i,x[i],y[i],z[i]);
  }

  //Defines graphs and axis. 
  gr = new TGraph(n,x,z);
  gs = new TGraph(n,x,t);
  gr->GetYaxis()->SetRangeUser(0.001,5000);

  // Attractive features for graphs
  gr->SetLineColor(kRed);
  gs->SetLineColor(kBlue);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(kOpenCircle);
  gs->SetMarkerColor(kBlue);
  gs->SetMarkerStyle(kOpenCircle);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Jet Pt Threshold (GeV/c)");
  gr->GetYaxis()->SetTitle("Rate Estimate (Hz)");

  //Draw graphs
  gr->Draw("LAP");
  gs->Draw("LP");

  //Defines Legened
  TLegend *L = new TLegend(0.52292,0.675103,0.932471,0.919831);
  L->SetHeader("AllPhysics, |#eta| < 2, All Centrality" );
  L->SetBorderSize(0);
  L->SetFillStyle(0);
  L->AddEntry(gr,"Uncorrected p_{T}","pl");
  L->AddEntry(gs,"Corrected p_{T}","pl");
  L->Draw();

}
