#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>

void formatHist(TH1* h, int col = 1, double norm = 1);
void saveCanvas(TCanvas* c, int date = 20080916);

void analyze_tracklet_deltaEta(const char* infile = "tracklet0915.root", double etaMax = 2.){
  const char* outfile = "output.root";

  /// Parameters
  double normRange = 0.6;
  double deltaCut = 0.2;
  int etaBins = 8;
  int nBins = 100;

  gROOT->Reset();
  gROOT->ProcessLine(".x rootlogon.C");
  gStyle->SetErrorX(0.);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  gStyle->SetPadLeftMargin(0.132);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.02);

  TH1::SetDefaultSumw2();

  TFile *f = new TFile(infile);
  f->cd("ana;1");
  TNtuple * ntparticle= dynamic_cast<TNtuple *>(f->Get("ana/ntparticle"));
  TNtuple * ntmatched= dynamic_cast<TNtuple *>(f->Get("ana/ntmatched"));
  TNtuple * ntInvMatched= dynamic_cast<TNtuple *>(f->Get("ana/ntInvMatched"));

  TNtuple * ntgen = dynamic_cast<TNtuple *>(f->Get("ana/ntgen"));
  TNtuple * ntevent = dynamic_cast<TNtuple *>(f->Get("ana/ntevent"));

  TFile *of = new TFile(outfile,"recreate");

  TH1F * h1 = new TH1F("h1",Form("Everything;#Delta#eta;#_{pixel pairs}/event/%.2f",1/(double)nBins),nBins,0,1); 
  TH1F * h2 = new TH1F("h2","Signal",nBins,0,1);
  TH1F * h3 = new TH1F("h3","Background",nBins,0,1);
  TH1F * h4 = new TH1F("h4","Normalized Reproduced Background",nBins,0,1);
  TH1F * h5 = new TH1F("h5","",nBins,0,1);
  TH1F * h6 = new TH1F("h6","Reproduced Background Subtracted",nBins,0,1);
  TH1F * h7 = new TH1F("h7","",nBins,0,1);
  TH1F * h8 = new TH1F("h8","",nBins,0,1);
  TH1F * h9 = new TH1F("h9","",nBins,0,1);
  TH1F * h10 = new TH1F("h10","",nBins,0,1);

  TProfile * dNdEtaHadron = new TProfile("dNdEtaHadron","",32,-2.1,2.1);
  TProfile * dNdEtaLepton = new TProfile("dNdEtaLepton","",32,-2.1,2.1);
  TProfile * dNdEtaTracklet = new TProfile("dNdEtaTracklet","",32,-2.1,2.1);

  // TH2D * dNdEtaHadron = new TH2D("dNdEtaHadron","",32,-2.1,2.1,50,0,5);
  // TH2D * dNdEtaLepton = new TH2D("dNdEtaLepton","",32,-2.1,2.1,50,0,5);

  TH2D * corr = new TH2D("correlation","; #_{hadrons}; #_{tracklets}",100,0,50,100,0,50);

  Float_t matchedeta1;
  Float_t matchedeta2;
  Float_t matchedinveta2;
  Float_t signalcheck;
  Float_t layer1hits;
  Float_t charge;
  Float_t ntpartlayer1hits;
 
  Float_t eta1;
  Float_t eta2;

  ntmatched->SetBranchAddress("eta1",&matchedeta1);
  ntmatched->SetBranchAddress("matchedeta",&matchedeta2);
  ntInvMatched->SetBranchAddress("matchedeta",&matchedinveta2);
  ntmatched->SetBranchAddress("signalCheck",&signalcheck);
  //  ntmatched->SetBranchAddress("layer1hits",&layer1hits);
  ntparticle->SetBranchAddress("eta1",&eta1);
  ntparticle->SetBranchAddress("eta2",&eta2);
  ntparticle->SetBranchAddress("charge",&charge);
  ntparticle->SetBranchAddress("layer1hits",&ntpartlayer1hits);


  int nevents = ntgen->GetEntries();
  int partentries = ntparticle->GetEntries();
  int matchedentries = ntmatched->GetEntries();

  for(int i = 0; i<matchedentries;i++){
    ntmatched->GetEntry(i);
    ntInvMatched->GetEntry(i);

    if(fabs(matchedeta1)>etaMax) continue;
    if(fabs(matchedeta2)>etaMax) continue;
    if(fabs(matchedinveta2)>etaMax) continue;
    //    if(layer1hits>20) continue;
    //    if(layer1hits<0) continue;    
    h1->Fill(fabs(matchedeta1-matchedeta2));
    h6->Fill(fabs(matchedeta1-matchedeta2));
    if(signalcheck==1){
      h2->Fill(fabs(matchedeta1-matchedeta2));
    }
    if(signalcheck==0){
      h3->Fill(fabs(matchedeta1-matchedeta2));
    } 
    h4->Fill(fabs(matchedeta1-matchedinveta2));
  }
  
  for(int i = 0; i<partentries;i++){
    ntparticle->GetEntry(i);
    if(fabs(eta1)>etaMax || fabs(eta2)>etaMax) continue;
    if(charge==0) continue;
    h5->Fill(fabs(eta1-eta2));
  }

  formatHist(h1,1,nevents);
  formatHist(h2,2,nevents);
  formatHist(h3,3,nevents);
  formatHist(h4,4,nevents);
  formatHist(h5,5,nevents);
  formatHist(h6,6,nevents);

  //// Normalization of background
  Float_t sc = ((h1->Integral((int)(normRange*nBins),nBins+1,"width"))/(h4->Integral((int)(normRange*nBins),nBins+1,"width")));
  cout<<"background normalization: "<<sc<<endl;
  h4->Scale(sc);
  h6->Add(h4,-1);

  //// Determination of correction factor beta
  double beta = 1-((h2->Integral(0,(int)(deltaCut*nBins),"width"))/(h6->Integral(0,(int)(deltaCut*nBins),"width")));
  cout<<"beta: "<<beta<<endl;

  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->SetLogy();

  TLegend * leg1 = new TLegend(0.25,0.66,0.56,0.84);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(h1,"Everything","l");
  leg1->AddEntry(h2,"Signal","l");
  leg1->AddEntry(h3,"Background","l");
  leg1->AddEntry(h4,"Normalized Reproduced Background","l");
  leg1->AddEntry(h6,"Reproduced Background Subtracted","l");

  h1->Draw("hist");
  h2->Draw("hist same");
  h3->Draw("hist same");
  h4->Draw("hist same");
  h6->Draw("hist same");
  leg1->Draw();

  // Might be necessary to increase the sizes here
  float had[8] = {0,0,0,0,0,0,0,0};
  float lep[8] = {0,0,0,0,0,0,0,0};
  float trt[8] = {0,0,0,0,0,0,0,0};

  //Important information:
  //
  //The ntuple is filled in event-by-event basis
  //therefore the total number of particles in a certain
  //eta bin is recorded by the hardcoded variable names
  //had1,had2,had3 etc...
  //Here this is corrected to be from 0 to 7
  //The way to interpret the etabins is this:
  //had[0] = had1 : -2, -1.5
  //had[1] = had2 : -1.5, 1
  //...
  //had[6] = had7 : 1, 1.5
  //had[7] = had8 : 1.5, 2
  // had: number of hadrons
  // lep: number of leptons
  // trt: number of tracklets

  for(int ig = 0; ig < etaBins; ++ig){
    ntgen->SetBranchAddress(Form("had%d",ig+1),&(had[ig]));
    ntgen->SetBranchAddress(Form("lep%d",ig+1),&(lep[ig]));
    ntevent->SetBranchAddress(Form("trt%d",ig+1),&(trt[ig]));
    //ntevent->SetBranchAddress(Form("had%d",ig+1),&(trt[ig])); // in the old run
  }

  for(int iev = 0; iev < nevents; ++iev){
    ntgen->GetEntry(iev);
    ntevent->GetEntry(iev);
    double nhad = 0;
    double ntrt = 0;
    for(int ieta = 0; ieta < etaBins; ++ieta){
      dNdEtaHadron->Fill((0.5*(double)ieta-1.75),had[ieta]);
      dNdEtaLepton->Fill((0.5*(double)ieta-1.75),lep[ieta]);
      dNdEtaTracklet->Fill((0.5*(double)ieta-1.75),trt[ieta]);
      nhad += had[ieta];
      ntrt += trt[ieta];
    }
    corr->Fill(nhad,ntrt);
  }

  TCanvas* c2 = new TCanvas("Correlation","",700,700);
  corr->Draw("col");

  TCanvas* c3 = new TCanvas("dNdEta","",700,400);

  TLegend * leg2 = new TLegend(0.42,0.26,0.74,0.56);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(dNdEtaTracklet,"Tracklets including Background","lp");
  leg2->AddEntry(dNdEtaTracklet,"#Delta#eta < 0.02","");
  leg2->AddEntry(dNdEtaHadron,"Hadrons","lp");
  leg2->AddEntry(dNdEtaLepton,"Leptons","lp");

  formatHist(dNdEtaTracklet,1);
  formatHist(dNdEtaHadron,4);
  formatHist(dNdEtaLepton,3);

  dNdEtaTracklet->SetMaximum(1.8);
  dNdEtaTracklet->Draw("");
  dNdEtaHadron->Draw("same");
  dNdEtaLepton->Draw("same");
  leg2->Draw();

  TCanvas* c4 = new TCanvas("Alpha","",700,400);
  
  TH1D* histTrt = dNdEtaTracklet->ProjectionX("dNdEtaRECO");
  TH1D* histAlpha = dNdEtaHadron->ProjectionX("alpha");

  histAlpha->SetTitle(";#eta;#alpha = #_{primaries}/#_{tracklets}");

  formatHist(histTrt);
  formatHist(histAlpha);
  
  histTrt->Scale(1.-beta);

  histAlpha->SetMaximum(1.4);
  histAlpha->Divide(histTrt);

  histAlpha->Draw();

  saveCanvas(c1);
  saveCanvas(c2);
  saveCanvas(c3);
  saveCanvas(c4);

}

void formatHist(TH1* h, int col, double norm){

  h->Scale(1/norm);
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

}

void saveCanvas(TCanvas* c, int date){
  c->Write();
  c->Draw();
  c->Print(Form("./figures/%s_d%d.gif",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.eps",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.C",c->GetName(),date));

}
