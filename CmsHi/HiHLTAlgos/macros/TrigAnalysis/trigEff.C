// 18 July 2011. This macro calculates the efficiency for various triggers for data and monte carlo. 

#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include "TChain.h"
#include <TH1F.h>
#include <TCut.h>
#include <TLegend.h>
#include <TLine.h>
#include <TCanvas.h>
#include <iostream> 

#include "commonUtility.h"
#include "tgraphTools.C"

using namespace std;


// Chooses bin number and bin content range. 
const int nBin=25; 
Float_t bin[nBin+1]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,140,160,200,240,300}; 

//Defines Basic efficiency Function
TGraphAsymmErrors *calcEff1(TH1* h1, TH1* h2,char *hName="hEff") // We use BayesDivide to ensure that the error bars don't do something nonsensical such as giving this result a value above 1.
{  
  TGraphAsymmErrors *gEfficiency = new TGraphAsymmErrors()
   gEfficiency->BayesDivide(h2,h1);
   return gEfficiency;
}



// Defines more complex efficiency function with more output information.
TGraphAsymmErrors *eff(TTree * t,TString var="nljet",TCut sel="",TCut cut="",TString tag="") 
{
  cout << "var: " << var << endl;
  cout << "Evt Sel: " << TString(sel) << ": " << t->GetEntries(sel) << endl;
  cout << "Trigger: " << TString(sel&&cut) << ": " << t->GetEntries(sel&&cut) << endl;


  TH1F *h = new TH1F("h"+tag,"",nBin,bin);
  TH1F *hCut = new TH1F("hCut"+tag,"",nBin,bin);
  t->Draw(var+">>h"+tag,sel,"");
  t->Draw(var+">>hCut"+tag,sel&&cut,"");
  TGraphAsymmErrors *g = calcEff1(h,hCut);
   clearXErrorBar(g);
  return g;
}



// Inputs tree information and graphs efficiency as a function of pt. 
TChain * trigEff(TString infile="alldata.root")
// TChain * trigEff(TString infile="allgen.root")
{
  TChain * tree = new TChain("hltanalysis/HltTree");
  tree->Add(infile);
	tree->AddFriend("icPu5JetAnalyzer/t",infile);
  cout << "Total: " << tree->GetEntries() << endl;

  TCut evtSel("HLT_HIMinBiasHfOrBSC&&hiBin<40&&abs(jteta[0])<2");
  
  TGraphAsymmErrors *g0=0,*g1=0,*g2=0,*g3=0,*g4=0;
     g0 = eff(tree,"jtpt[0]",evtSel,"HLT_HIJet35U","Jet35U");
  // g0 = eff(tree,"genpt[0]",evtSel,"HLT_HIJet35U","Jet35U");
  // g0 = eff(tree,"jtpt[0]",evtSel,"L1_SingleJet30_BptxAND","Single Jet 30");
  // g0 = eff(tree,"genpt[0]",evtSel,"L1_SingleJet30","Single Jet 30");
  
     //Draws Canvas
  TCanvas *cTrigEff2 = new TCanvas("cTrigEff2","cTrigEff2",600,550);
  TH1F *hTmp = new TH1F("hTmp","",nBin,bin);
    hTmp->SetAxisRange(0,1.3,"Y");
   hTmp->SetXTitle("Leading Jet p_{T} (GeV/c)"); 
    //hTmp->SetXTitle("Gen Jet p_{T} (GeV/c)"); // Use for Monte Carlo
  hTmp->SetYTitle("HLT  Eff. (Trigger/MB)");
   // hTmp->SetYTitle(" HLT Eff. (Trigger/Generated Jets)"); // Use for Monte Carlo
    handsomeTH1(hTmp);
  hTmp->Draw();
  g0->Draw("p");

  //Draws Legend
  TLine *l = new TLine(0,1,bin[nBin],1);
  l->SetLineStyle(2);
  l->Draw();
   TLegend *t = new TLegend(0.251678,0.742366,.927852,0.9561);
  //TLegend *t = new TLegend(0.236577,0.742366,1,0.956107); //Use for Monte Carlo
  t->SetHeader("All Physics |#eta|<2, All Centrality");
  //("PYTHIA+HYDJET |#eta|<2, All Centrality"); //Use for Monte Carlo
  t->SetBorderSize(0);
  t->SetFillStyle(0);
  t->Draw();

  cTrigEff2->Print("trigeff(1).pdf"); // Must change name after each change in trigger selection. 

  return tree;
}
