#include <iostream>
#include <math.h>
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
#include <TLine.h>
#include <TMath.h>
#include "../../PixelTracklet/interface/TrackletCorrections.h"

struct TrackletData
{
  TrackletData(
               TNtuple *nt1=0,
               TNtuple *nt2=0,
               TNtuple *nt3=0,
               TNtuple *nt4=0,
               TNtuple *nt5=0,
               TNtuple *nt6=0,
	       TrackletCorrections* tc = 0
               ) :
    ntparticle(nt1),
       ntmatched(nt2),
       ntInvMatched(nt3),
       ntgen(nt4),
       ntevent(nt5),
       ntvertex(nt6),
       corr(tc)
  {;}
  
  TrackletData(TFile* f, TrackletCorrections* tc = 0){

    corr = tc;  // It wouldn't copy anything, right?
    f->cd("ana;1");
    ntparticle= dynamic_cast<TNtuple *>(f->Get("anasim/ntparticle"));
    ntmatched= dynamic_cast<TNtuple *>(f->Get("ana/ntmatched"));
    ntInvMatched= dynamic_cast<TNtuple *>(f->Get("ana/ntInvMatched"));
    ntgen = dynamic_cast<TNtuple *>(f->Get("ana/ntgen"));
    ntevent = dynamic_cast<TNtuple *>(f->Get("ana/ntevent"));
    ntvertex = dynamic_cast<TNtuple *>(f->Get("ana/ntvertex"));
    ntcorr = new TNtuple("ntcorr","","nhad:ntrt:nstrt:vz:nvtx:vntrk");

    counter = 0;

    ntvertex->SetBranchAddress("z",&vz);
    ntvertex->SetBranchAddress("n",&vntrk);
    ntvertex->SetBranchAddress("nvtx",&nvtx);

    ntmatched->SetBranchAddress("eta1",&matchedeta1);
    ntmatched->SetBranchAddress("phi1",&matchedphi1);
    ntmatched->SetBranchAddress("matchedeta",&matchedeta2);
    ntmatched->SetBranchAddress("matchedphi",&matchedphi2);
    ntmatched->SetBranchAddress("signalCheck",&signalcheck);

    //    ntmatched->SetBranchAddress("sid",&signalcheck);

    ntmatched->SetBranchAddress("deta",&deta);
    ntmatched->SetBranchAddress("dphi",&dphi);
    ntmatched->SetBranchAddress("nhit1",&layer1hits);

    ntInvMatched->SetBranchAddress("eta1",&inveta1);
    ntInvMatched->SetBranchAddress("phi1",&invphi1);
    ntInvMatched->SetBranchAddress("matchedeta",&inveta2);
    ntInvMatched->SetBranchAddress("matchedphi",&invphi2);
    ntInvMatched->SetBranchAddress("deta",&invdeta);
    ntInvMatched->SetBranchAddress("dphi",&invdphi);

    ntparticle->SetBranchAddress("eta1",&eta1);
    ntparticle->SetBranchAddress("eta2",&eta2);
    ntparticle->SetBranchAddress("charge",&charge);
    ntparticle->SetBranchAddress("layer1hits",&ntpartlayer1hits);



  }

  ~TrackletData(){;}
  TNtuple * ntparticle;
  TNtuple * ntmatched;
  TNtuple * ntInvMatched;
  TNtuple * ntgen;
  TNtuple * ntevent;
  TNtuple * ntvertex;
  TNtuple * ntcorr;

  TrackletCorrections* corr;

  Float_t matchedeta1;
  Float_t matchedeta2;
  Float_t matchedinveta2;
  Float_t matchedinvphi2;
  Float_t signalcheck;
  Float_t layer1hits;
  Float_t charge;
  Float_t ntpartlayer1hits;
  Float_t matchedphi1;
  Float_t matchedphi2;

  Float_t eta1;
  Float_t eta2;
  Float_t deta;
  Float_t dphi;
  Float_t invdeta;
  Float_t invdphi;
  Float_t inveta1;
  Float_t inveta2;
  Float_t invphi1;
  Float_t invphi2;
  Float_t vz;
  Float_t vntrk;
  Float_t nvtx;

  double getBeta(int bin, bool saveplots = false);
  double getAlpha(int bin, bool saveplots = false);
  double deltaR(double eta1, double phi1, double eta2, double phi2);
  double deltaR2(double eta1, double phi1, double eta2, double phi2);

  int counter;
};


double TrackletData::deltaR(double eta1, double phi1, double eta2, double phi2){
  return sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

double TrackletData::deltaR2(double eta1, double phi1, double eta2, double phi2){
  double c = this->corr->getCPhi();
  double pi = TMath::Pi();
  double eta = eta1-eta2;
  double phi = phi1-phi2;
  if(phi > 2*pi) phi -= 2*pi;
  if(phi > pi) phi = 2*pi-phi;
  return eta*eta+c*c*phi*phi;
}

void formatHist(TH1* h, int col, double norm){
  h->Scale(1/norm);
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->Write();
}

void saveCanvas(TCanvas* c, int date){
  c->Write();
  c->Draw();
  c->Print(Form("./figures/%s_d%d.gif",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.eps",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.C",c->GetName(),date));
}

