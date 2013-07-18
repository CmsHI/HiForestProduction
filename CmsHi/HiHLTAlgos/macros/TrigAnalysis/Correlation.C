// 19 July, 2011. Caluclates the overlaps and efficiencies for various triggers. 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TPaletteAxis.h"
#include "TCanvas.h"
#include "TSystem.h"

using namespace std;


vector <string>* triggers = new vector<string>;
vector <string>* triggerCuts = new vector<string>;
vector< vector <double*> >* results = new vector< vector <double*> >;
TString goutdir;
TString anaMode;


//Defines error function
double calcEffErr(double nEvt,double nAcc)
{
  double eff =  nAcc / nEvt;
  return sqrt(eff*(1-eff)/nEvt);
}


//Defines efficiency function
double* calcEff(TTree* HltTree, const char *title, double nEvt,const char *cut,int flag=1)
{
  double nAcc = HltTree->GetEntries(cut);
  double effErr = calcEffErr(nEvt,nAcc);
  double eff = nAcc/nEvt;
  double *result = new double[2];
  result[0]=eff;
  result[1]=effErr;


  // print out or not
  if (flag) {
    triggers->push_back(string(title));
    triggerCuts->push_back(string(cut));
    cout <<" | "<<setw(20)<<title;
    cout <<" | "<<setw(10)<<fixed<<setprecision(3)<<eff*100.<<"% +-";
    cout        <<setw(10)<<fixed<<setprecision(3)<<effErr*100.<<"%";
    cout <<" | "<<endl;
  }

  return result;
}

void printEff(TTree* HltTree,const char *cut,const char *title, char *projectTitle)
{
  cout <<"   * "<<title<<":"<<endl;
  cout <<"      * Efficiencies:"<<endl;
  double nEvt = HltTree->GetEntries(cut);

  cout <<" | "<<setw(20)<<"HLT Path";
  cout <<" | "<<setw(25)<<"Efficiency";
  cout <<" | "<<endl;

  triggers->clear();
  triggerCuts->clear();

  vector <double*> effs;

  // Chooses the triggers for the correlation matrix   
  effs.push_back(calcEff(HltTree,"HLT_HIJet35U",nEvt,Form("(%s)&&HLT_HIJet35U",cut)));
  cout << cut << endl;
  if (anaMode=="HLTMB") {
    effs.push_back(calcEff(HltTree,"HLT_HIJet50U",nEvt,Form("(%s)&&HLT_HIJet50U",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIJet75U",nEvt,Form("(%s)&&HLT_HIJet75U",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIPhoton15",nEvt,Form("(%s)&&HLT_HIPhoton15",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIPhoton20",nEvt,Form("(%s)&&HLT_HIPhoton20",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIPhoton30" ,nEvt,Form("(%s)&&HLT_HIPhoton30",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIL2DoubleMu0",nEvt,Form("(%s)&&HLT_HIL2DoubleMu0",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIL2Mu3",nEvt,Form("(%s)&&HLT_HIL2Mu3",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIL2DoubleMu3",nEvt,Form("(%s)&&HLT_HIL2DoubleMu3",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIL2Mu20",nEvt,Form("(%s)&&HLT_HIL2Mu20",cut)));
    effs.push_back(calcEff(HltTree,"HLT_HIL2Mu5Tight",nEvt,Form("(%s)&&HLT_HIL2Mu5Tight",cut)));
  }
  results->push_back(effs);

  cout <<"      * Correlation Matrix:"<<endl;
  int tsize = (int)triggers->size();
  TH2D *h = new TH2D(Form("h%s",title),"",tsize,0,tsize,tsize,0,tsize);
  TH2D *hct = new TH2D(Form("h%s_ct",title),"",tsize,0,tsize,tsize,0,tsize);

  for (int i=tsize-1;i>=0;i--){
    int nEvtAfterCut = HltTree->GetEntries((*triggerCuts)[i].c_str());
    h->GetXaxis()->SetBinLabel(i+1,(*triggers)[i].c_str());
    h->GetYaxis()->SetBinLabel(i+1,(*triggers)[i].c_str());
    hct->GetXaxis()->SetBinLabel(i+1,(*triggers)[i].c_str());
    hct->GetYaxis()->SetBinLabel(i+1,(*triggers)[i].c_str());
    for (int j=0;j<tsize;j++){
      string cut ="("+(*triggerCuts)[i]+")&&("+(*triggerCuts)[j]+")";
      double* eff = calcEff(HltTree,"",nEvtAfterCut,Form("%s",cut.c_str()),0);
      if (nEvtAfterCut==0) eff[0]=0;
      if (i!=j) {
	h->SetBinContent(i+1,j+1,int(eff[0]*1000)/10.);
	hct->SetBinContent(i+1,j+1,HltTree->GetEntries(cut.c_str()));
      }
    }
  }
  h->GetXaxis()->LabelsOption("v");
  h->SetLabelSize(0.05,"XY");
  hct->GetXaxis()->LabelsOption("v");

  TCanvas *c1 = new TCanvas(Form("c%s",title), Form("c_%s",title),1000,800);
  c1->Range(-3.609756,-1.910995,12.7561,10.60209);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetTickx();
  c1->SetTicky();
  c1->SetLeftMargin(0.37);
  c1->SetRightMargin(0.1684054);
  c1->SetTopMargin(0.02);
  c1->SetBottomMargin(0.4);
  c1->SetFrameLineColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameLineColor(0);
  c1->SetFrameBorderMode(0);

  h->Draw("col text");

  TPaletteAxis *palette = new TPaletteAxis(tsize*1.02,0,tsize*1.1,tsize,h);
  palette->SetLabelColor(1);
  palette->SetLabelFont(42);
  palette->SetLabelOffset(0.005);
  palette->SetLabelSize(0.045);
  palette->SetTitleOffset(1);
  palette->SetTitleSize(0.04);
  palette->SetFillColor(100);
  palette->SetFillStyle(1001);
  h->GetListOfFunctions()->Add(palette,"br");
  h->Draw("col text z");

  string fname(Form("%s/trigCorr_%s_%s.gif",goutdir.Data(),projectTitle,title));
  c1->SaveAs(fname.c_str());
  c1->SaveAs(Form("%s/trigCorr_%s_%s.C",goutdir.Data(),projectTitle,title));

  cout <<"<img src=\"%ATTACHURLPATH%/"<<fname<<"\" alt=\""<<fname<<"\" width='671'   height='478' />"<<endl;   
}

//Correlation function definition
void Correlation(
    TString mode="HLTMB", // L1MB, L1Algo, HLTMB, HLTAlgo
    Int_t runNum = 152957,
    Int_t goodLumiStart = 0,
    TString inFile0Name="/net/hidsk0001/d00/scratch/frankma/data/HAZS/hazsv2-run152791to152957-jetoh-v2/merge/all.root", 
    TString outdir="out/HR10AllPR2",
    char *projectTitle = "HIAllPhy2010",
    string source="data")
{
  // Load input
  TChain * HltTree = new TChain("hltanalysis/HltTree","HI OpenHLT Tree");
  HltTree->Add(inFile0Name);
  cout << inFile0Name << endl;
  cout << " # entries: " << HltTree->GetEntries() << endl;
  anaMode=mode;

  // Define Output
  TString evtSel("HLT_HIMinBiasHfOrBSC");//("L1Tech_BPTX_plus_AND_minus.v0_5bx"); 
  outdir+=Form("/run%d_%s_%s",runNum,anaMode.Data(),evtSel.Data());
  gSystem->mkdir(outdir.Data(),kTRUE);
  goutdir=outdir;
  TFile *outf = new TFile(Form("%s/%s_hist.root",outdir.Data(),projectTitle),"RECREATE");

  // define event types
    vector<string> evtType;
  vector<string> evtTypeCut;
  evtType.push_back("BptxAND"); evtTypeCut.push_back(evtSel.Data());
  if (source=="mc") {
  }
  else if (source=="data") {
    for (UInt_t i=0; i<evtTypeCut.size(); ++i) evtTypeCut[i]+=Form("&&Run==%d&&LumiBlock>%d",runNum,goodLumiStart);
  }

  // Print out event type fractions
  cout <<"<pre>"<<endl;
  for (unsigned int i=0; i<evtType.size(); ++i) {
    int nSel = HltTree->GetEntries(evtTypeCut[i].c_str());
    cout <<std::right<<std::setw(20)<<evtType[i] << ": " <<nSel<<endl;
  }
  cout <<"</pre>"<<endl;

  // Calculate efficiencies
  for (unsigned int i=0; i<evtType.size(); ++i) {
    printEff(HltTree,evtTypeCut[i].c_str(),evtType[i].c_str(),projectTitle);
  }

  // Print efficiency results
  cout <<" | "<<setw(20)<<" ";
  for (unsigned int i=0;i<evtType.size();i++) {
    cout <<" | "<<setw(8)<<evtType[i];
  }
  cout <<" | " <<endl;

  for (int i=0;i<(int)triggers->size();i++){
    cout <<" | " <<setw(20)<<(*triggers)[i];
    for (int j=0;j<(int)results->size();j++) {
      double eff=(((*results)[j])[i])[0];
      double effErr=(((*results)[j])[i])[1];
      cout <<" | "<<setw(10)<<fixed<<setprecision(3)<<eff*100.<<"% +-";
      cout        <<setw(10)<<fixed<<setprecision(3)<<effErr*100.<<"% ";
    }
    cout <<" | "<<endl;
  }   

  // save
  outf->Write();
}
