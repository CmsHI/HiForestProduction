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
#include "HisMath.C"

using namespace std;


// Chooses bin number and bin content range. 
// const int nBin=25; 
// Float_t bin[nBin+1]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,140,160,200,240,300}; 
const int nBin=33;
Float_t bin[nBin+1]={20,22,24,26,28,30,32,34,36,38,40,42,46,48,50,52,54,56,58,60,62,64,68,72,76,80,85,90,95,100,110,120,140,200};

//Defines Basic efficiency Function
TGraphAsymmErrors *calcEff1(TH1* h1, TH1* h2,char *hName="hEff") // We use BayesDivide to ensure that the error bars don't do something nonsensical such as giving this result a value above 1.
{  
  TGraphAsymmErrors *gEfficiency = new TGraphAsymmErrors();
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
void AcceptFracvPt()
{
   TString infile="/net/hisrv0001/home/yenjie/scratch/hiForest/prod/production/CMSSW_4_4_2_patch6/test/pPb/test_pPb.root";
   TChain * tree = new TChain("hltanalysis/HltTree");
   tree->Add(infile);
   tree->AddFriend("icPu5JetAnalyzer/t",infile);
   tree->AddFriend("hiEvtAnalyzer/HiTree",infile);
   cout << "Total: " << tree->GetEntries() << endl;
   
//    TCut evtSel("HLT_HIMinBiasHfOrBSC_v1&&hiBin<40&&abs(jteta[0])<2");
   TCut evtSel("hiBin<40");

   TH1F *hJetPt = new TH1F("hJetPt",";leading jet p_{T} (GeV/c)",50,0,200);
   float nSel = tree->Draw("jtpt[0]>>hJetPt",evtSel);
   
   hJetPt->Scale(1./nSel);
   
   TCanvas * c2 = new TCanvas("c2","",500,500);
   hJetPt->SetYTitle("a.n.");
   hJetPt->Draw();
   c2->Print("fig/pB/jet_pt.gif");

   TCanvas * c3 = new TCanvas("c3","",500,500);
   c3->SetLogy();
   TH1D * hIntPt = integrateHist(hJetPt,false,-1);
   hIntPt->SetTitle(";Jet p_{T} Threshold (GeV/c);Accept Fraction (from MB)");
   hIntPt->Draw();
   c3->Print("fig/pB/accept_fraction.gif");
}
