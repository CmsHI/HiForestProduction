#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
using namespace std;

void analyze_tree(){
  TFile* inf = new TFile("tree.root");
  TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));

  int zbins = 1;
  int hitbins = 20;
  int nbins = zbins*hitbins;
  
  vector<TH1D*> h;
  h.reserve(nbins);
  for(int i = 0; i< nbins; ++i){
    h[i] = new TH1D(Form("dNdEtaHits_%02d",i),"dNdEta Hits Layer 1",100,-3,3);
  }
  
  TH1D* h1 = new TH1D("h1","dNdEta Hits Layer 1",100,-3,3);
  TH1D* h2 = new TH1D("h2","dNdEta Hits Layer 2",100,-3,3);

  TH1D* hm1 = new TH1D("hm1","Number of Hits Layer 1",50,0,50);

  float eta1[1000],phi1[1000],r1[1000],eta2[1000],phi2[1000],r2[1000],vz[100];

  int nhits1,nhits2,mult,nv;

  t->SetBranchAddress("eta1",eta1);
  t->SetBranchAddress("phi1",phi1);
  t->SetBranchAddress("r1",r1);
  t->SetBranchAddress("eta2",eta2);
  t->SetBranchAddress("phi2",phi2);
  t->SetBranchAddress("r2",r2);
  t->SetBranchAddress("nhits1",&nhits1);
  t->SetBranchAddress("nhits2",&nhits2);
  t->SetBranchAddress("mult",&mult);
  t->SetBranchAddress("vz",vz);
  t->SetBranchAddress("nv",&nv);

  for(int i = 0; i< t->GetEntries(); ++i){
    t->GetEntry(i);
    std::cout<<"Number of Hits : "<<nhits1<<endl;
    std::cout<<"Hit Multiplicity : "<<mult<<endl;

    hm1->Fill(mult);
    for(int ihit = 0; ihit < nhits1; ++ihit){
      std::cout<<"Hit Eta : "<<eta1[ihit]<<endl;
      h1->Fill(eta1[ihit]);
    }

    for(int ihit2 = 0; ihit2 < nhits2; ++ihit2){
      std::cout<<"Hit Eta : "<<eta1[ihit2]<<endl;
      h2->Fill(eta2[ihit2]);
    }

  }
  h1->Draw();
  hm1->Draw();

  t->Draw("nhits1");

}
