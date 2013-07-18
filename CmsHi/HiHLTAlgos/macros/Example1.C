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

void Example1(){
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
    // please select the trigger from triggers arrays
    //selected trigger will be plotted versus Hf Energy, Centrality, Charged Particle multiplicity and TrackMultiplicity
    int TrSelect = 3; 
    //specify the roottuple path
    TFile* inf = new TFile("openhlt.root");  

    string triggers[4] = {"L1_SingleJet30","L1_TripleJet30","L1_QuadJet15","L1_DoubleJet70"};
    string vers[4] = {"hiHF","hiBin","Ncharged","hiNtracks"};

    double bins [4] = {170,40,240,140};
    double limits[4] = {170000,40,24000,1400};

    TCanvas* c1 = new TCanvas("c1","c1",800,800);

    c1->Divide(2,2);

    TNtuple* nt = (TNtuple*)inf->Get("HltTree");

    TProfile* p1[4];


    for(unsigned int i =0; i<4 ; i++){
        c1->cd(i+1);

        p1[i] = new TProfile("p1",Form(";%s;%s",vers[i].data(),triggers[TrSelect].data()),bins[i],0,limits[i]);
        nt->SetAlias("trigger",triggers[TrSelect].data());

        nt->SetAlias("versus",vers[i].data());
        nt->Draw("trigger:versus>>p1","","prof");

    }

    // save the plot current directory

    c1->Print(Form("%s_vs4.gif",triggers[TrSelect].data()));

}

