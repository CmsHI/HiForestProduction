#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TInterpreter.h>
#include <TGraphAsymmErrors.h>
#endif

void Example2(){
    gROOT->Reset();
    gROOT->Clear();
    gROOT->SetStyle("Plain");

    gStyle->SetTextSize(0.01908148);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.55);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    // select the one of the versus : Hf energy, Centrality, # Charged , # Tracks
    int Vselect = 1;

    // number of Trigger you will use , 
    const int Ntr = 2; 
    //if you want to use more than two trigger please specify the Ntr2 for dividivg the canvas
    const int Ntr2 = 1;
    TFile* inf = new TFile("openhlt2.root");

    string triggers[6] = {"L1Tech_BSC_minBias_threshold1.v0","L1Tech_BSC_minBias_threshold2.v0","L1_SingleJet30","L1_TripleJet30","L1_QuadJet15","L1_DoubleJet70"};
    string vers[4] = {"hiHF","hiBin","Ncharged","hiNtracks"};

    double bins [4] = {170,40,240,140};
    double limits[4] = {170000,40,24000,1400};

    TCanvas* c1 = new TCanvas("c1","c1",800,400);
    c1->Divide(Ntr,Ntr2);

    TNtuple* nt = (TNtuple*)inf->Get("HltTree");

    TProfile* p1[Ntr];

    for(unsigned int i =0; i<Ntr ; i++){
        c1->cd(i+1);

        p1[i] = new TProfile("p1",Form(";%s;%s",vers[Vselect].data(),triggers[i].data()),bins[1],0,limits[1]);
        nt->SetAlias("trigger",triggers[i].data());

        nt->SetAlias("versus",vers[Vselect].data());
        nt->Draw("trigger:versus>>p1","","prof");

    }

    c1->Print(Form("%s_vsCent.gif",triggers[Vselect].data()));

}
