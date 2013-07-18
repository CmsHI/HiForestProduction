#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream.h>
#include <TLine.h>
#include <math.h>

void formatHist(TH1* h, int col = 1, double norm = 1);
void saveCanvas(TCanvas* c, int date = 20080829);

double* getBeta(const char* infile, double etaMax, int MinHit, int MaxHit);

void analyze_beta(const char* infile = "p0829.root", double etaMax = 1.){

   int di=4;
   TH1F *beta = new TH1F("beta","",48/di,0,48);

   for (int i=0;i<48;i+=di) {
      double* myBeta = getBeta(infile,etaMax,i,i+di);
      cout <<myBeta[0]<<" "<<myBeta[1]<<endl;
      beta->SetBinContent(i/di+1,myBeta[0]);
      beta->SetBinError(i/di+1,myBeta[1]);
   }

   TCanvas *c = new TCanvas ("c","",400,400);
   beta->Fit("pol6","m");
   beta->SetXTitle("N_{Hits}");
   beta->SetYTitle("#beta");
   beta->SetAxisRange(0,0.3,"Y");
   beta->Draw("p");
}

double* getBeta(const char* infile, double etaMax, int MinHit, int MaxHit)
{
  /// Parameters
  double normRange = 0.5;
  int nBins = 40;
  double deltaCut = 0.5 / (double)nBins * 2;
  int etaBins = 8;

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
  TNtuple * ntmatched= dynamic_cast<TNtuple *>(f->Get("ana/ntmatched"));
  TNtuple * ntInvMatched= dynamic_cast<TNtuple *>(f->Get("ana/ntInvMatched"));
  TNtuple * ntgen = dynamic_cast<TNtuple *>(f->Get("ana/ntgen"));
  TNtuple * ntevent = dynamic_cast<TNtuple *>(f->Get("ana/ntevent"));
  TNtuple * ntvertex = dynamic_cast<TNtuple *>(f->Get("ana/ntvertex"));
 
  TNtuple * ntcorr = new TNtuple("ntcorr","","nhad:ntrt:nstrt:vz:nvtx:vntrk");
  TH1F * h1 = new TH1F("h1",Form("Everything;D;#_{pixel pairs}/event/%.2f",1/(double)nBins),nBins,0,2); 
  TH1F * h2 = new TH1F("h2","Signal",nBins,0,2);
  TH1F * h3 = new TH1F("h3","Background",nBins,0,2);
  TH1F * h4 = new TH1F("h4","Normalized Reproduced Background",nBins,0,2);
  TH1F * h5 = new TH1F("h5","",nBins,0,2);
  TH1F * h6 = new TH1F("h6","Reproduced Background Subtracted",nBins,0,2);

  TProfile * dNdEtaHadron = new TProfile("dNdEtaHadron","",32,-2.1,2.1);
  TProfile * dNdEtaLepton = new TProfile("dNdEtaLepton","",32,-2.1,2.1);
  TProfile * dNdEtaTracklet = new TProfile("dNdEtaTracklet","",32,-2.1,2.1);

  // TH2D * dNdEtaHadron = new TH2D("dNdEtaHadron","",32,-2.1,2.1,50,0,5);
  // TH2D * dNdEtaLepton = new TH2D("dNdEtaLepton","",32,-2.1,2.1,50,0,5);

  TH2D * corr = new TH2D("correlation","; #_{hadrons}; #_{tracklets}",100,0,50,100,0,50);

  Float_t matchedeta1;
  Float_t matchedeta2;
  Float_t matchedinveta2;
  Float_t matchedinvphi2;
  Float_t signalcheck;
  Float_t layer1hits;
  Float_t invlayer1hits;
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


  ntvertex->SetBranchAddress("z",&vz);
  ntvertex->SetBranchAddress("n",&vntrk);
  ntvertex->SetBranchAddress("nvtx",&nvtx);

  ntmatched->SetBranchAddress("eta1",&matchedeta1);
  ntmatched->SetBranchAddress("phi1",&matchedphi1);
  ntmatched->SetBranchAddress("matchedeta",&matchedeta2);
  ntmatched->SetBranchAddress("matchedphi",&matchedphi2);
  ntmatched->SetBranchAddress("signalCheck",&signalcheck);
  ntmatched->SetBranchAddress("deta",&deta);
  ntmatched->SetBranchAddress("dphi",&dphi);
  ntmatched->SetBranchAddress("nhit1",&layer1hits);

  ntInvMatched->SetBranchAddress("eta1",&inveta1);
  ntInvMatched->SetBranchAddress("phi1",&invphi1);
  ntInvMatched->SetBranchAddress("matchedeta",&inveta2);
  ntInvMatched->SetBranchAddress("matchedphi",&invphi2);
  ntInvMatched->SetBranchAddress("deta",&invdeta);
  ntInvMatched->SetBranchAddress("dphi",&invdphi);
  ntInvMatched->SetBranchAddress("nhit1",&invlayer1hits);
  
  int nevents = ntgen->GetEntries();
  int matchedentries = ntmatched->GetEntries();


  for(int i = 0; i<matchedentries;i++){
    ntmatched->GetEntry(i);
    if(fabs(matchedeta1)>=etaMax) continue;
    if(fabs(matchedeta2)>etaMax) continue;
    if(fabs(matchedeta1)<0.1) continue;
    if(layer1hits>MaxHit) continue;
    if(layer1hits<MinHit) continue;    
    float dR= sqrt(deta*deta+dphi*dphi/43./43.);


    h1->Fill(fabs(dR));
    h6->Fill(fabs(dR));
    if(signalcheck==1){
      h2->Fill(fabs(dR));
    }
    if(signalcheck==0){
      h3->Fill(fabs(dR));
    } 
  }

  int invmatchedentries = ntInvMatched->GetEntries();

  for(int i = 0; i<invmatchedentries;i++){
    ntInvMatched->GetEntry(i);
    if(fabs(inveta1)>etaMax) continue;
    if(fabs(inveta1)<0.1) continue;
    if(invlayer1hits>MaxHit) continue;
    if(invlayer1hits<MinHit) continue;    
    

    float dR= sqrt(invdeta*invdeta+invdphi/43.*invdphi/43.);
    h4->Fill(dR);
  }


  formatHist(h1,1,1);
  formatHist(h2,2,1);
  formatHist(h3,3,1);
  formatHist(h4,4,1);
  formatHist(h6,6,1);

  //// Normalization of background

  Float_t nSideBand = h1->Integral((int)(normRange*nBins),nBins,"width");
  Float_t nSBErr = sqrt(nSideBand);
  Float_t nReproducedSideBand = (h4->Integral((int)(normRange*nBins),nBins,"width"));
  Float_t nRSBErr = sqrt(nReproducedSideBand);
  Float_t sc = (nSideBand/nReproducedSideBand);
  Float_t scerr = sc * sqrt((nRSBErr/nReproducedSideBand)*(nRSBErr/nReproducedSideBand)
                           +(nSBErr/nSideBand)*(nSBErr /nSideBand));

  cout<<"nSideBand: "<<nSideBand<<" +- "<<nSBErr<<endl;
  cout<<"nReproducedSideBand: "<<nReproducedSideBand<<" +- "<<nRSBErr<<endl;

  cout<<"background normalization: "<<sc<<endl;
  h4->Scale(sc);
  h6->Add(h4,-1);
  
  //// Determination of correction factor beta
  double beta,betaErr=0 ;
  
  
  if ((h6->Integral(0,(int)(deltaCut*nBins),"width"))==0) {
     beta=0;
     betaErr=0;
  } else {
     Float_t nBck = (h4->Integral(0,(int)(deltaCut*nBins),"width"))/sc;
     Float_t nTotal = (h1->Integral(0,(int)(deltaCut*nBins),"width"));
     Float_t nBckErr = sqrt(1/nBck+(scerr/sc)*(scerr/sc))*nBck*sc;
     Float_t nTErr = sqrt(nTotal);
     beta=(nBck*sc)/(nTotal);
     betaErr=sqrt((nBckErr/nBck/sc)*(nBckErr/nBck/sc)+(nTErr/nTotal)*(nTErr/nTotal))*beta;
  }
  cout<<"sc: "<<sc<<" +- "<<scerr<<endl;

  cout<<"beta: "<<beta<<" +- "<<betaErr<<endl;

  double* val = new double[2];
  val[0]=beta;
  val[1]=betaErr;
  return val;
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
}
