// 18 July, 2011. This macro provides a comparison between the pt spectrum of a cleaned and not cleaned HLT_HIJet35U trigger data sample. 

#include "TChain.h"
#include "TH1.h"
#include "TCut.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

void ptspectrum(TString infile="alldatacleaned.root") //contains both cleaned and uncleaned data. 
{
   TChain * tree = new TChain("icPu5JetAnalyzer_clean/t");
  tree->Add(infile);
	tree->AddFriend("hltanalysis_clean/HltTree",infile);
	//  cout << "Total: " << tree->GetEntries() << endl; 
   TChain * tree2 = new TChain("icPu5JetAnalyzer/t");
  tree2->Add(infile);
	tree2->AddFriend("hltanalysis/HltTree",infile);
	
	//Trigger Selection
  TCut evtSel("HLT_HIJet35U&&hiBin<40");

  //Information assignment
  TH1D * hJetPt = new TH1D("hJetPt",";Leading Jet p_{T} (GeV/c); Events;",85,0,500);
     TCanvas * c4 = new TCanvas("c4","Leading Jet p_{T}" ,600,600);
     c4->SetLogy();
     tree->Draw("jtpt[0]>>hJetPt",evtSel&&"abs(jteta[0])<2");
  TH1D * hJetPt2 = new TH1D("hJetPt2",";Leading Jet p_{T} (GeV/c);Events;",85,0,500);
     TCanvas * c5 = new TCanvas("c5","Leading Jet p_{T}" ,600,600);
     c5->SetLogy();
     tree2->Draw("jtpt[0]>>hJetPt2",evtSel&&"abs(jteta[0])<2");

    //Set Colors
    hJetPt->SetLineColor(kOrange+10);
    hJetPt->SetMarkerColor(kOrange+10);
    hJetPt2->SetLineColor(kGreen+2);
    hJetPt2->SetMarkerColor(kGreen+2);
    
    //Make Canvas
   TCanvas * g4 = new TCanvas("g4", "pt spectrum", 800, 500);
     g4->SetLogy();
  // hJetPt->SetAxisRange(1,6000,"Y");
  hJetPt->Draw("");
  hJetPt2->Draw("same");

  //Make Legened
TLegend *t = new TLegend(0.431208,0.609756,0.919463,0.921603);
   t->SetHeader("AllPhysics, |#eta|<2, All centrality ");
   t->SetBorderSize(0);
   t->SetFillStyle(0);
   t->AddEntry(hJetPt,"Cleaned","pl");
   t->AddEntry(hJetPt2,"Not Cleaned","pl");
   t->Draw();
       


}
