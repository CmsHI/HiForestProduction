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
#include <TF1.h>

void formatHist(TH1* h, int col = 1, double norm = 1);
void saveCanvas(TCanvas* c, int date = 20080829);

double* getBeta(TFile *e, double etaMin, double etaMax, int MinHit, int MaxHit);
void prepareHist(TFile* f,TFile* g, double etaMin, double etaMax);

/// Parameters
double normRange = 0.2;
int nBins = 25;
double deltaCut = 0.1 / (double)nBins * 25;
int etaBins = 8;
int maxHit = 100;
int di=1;
double vzCut=10;



TH1F *h1[200];
TH1F *h2[200];
TH1F *h3[200];
TH1F *h4[200];
TH1F *h5[200];
TH1F *h6[200];;


void analyze_betasFromTree(const char* infile , double etaMin = 0, double etaMax = 1.)
{

   TFile *f = new TFile(infile);
   f->cd("ana;1");
   TFile *g = new TFile("pseudoExp.root");

   prepareHist(f,g,etaMin,etaMax);

   TH1F *beta = new TH1F(Form("beta%d",(int)(etaMin*10+30)),"",maxHit/di,0,maxHit);
   formatHist(beta,2,1);

   for (int i=0;i<maxHit;i+=di) {
      double* myBeta = getBeta(f,etaMin, etaMax,i,i+di);
      //cout <<myBeta[0]<<" "<<myBeta[1]<<endl;
      if (myBeta[1]!=100) {
         beta->SetBinContent(i/di+1,myBeta[0]);
         beta->SetBinError(i/di+1,myBeta[1]);
      }
   }

   TF1 *fun = new TF1("fun","[2]*x+[1]/x+[0]",0,100);


   TCanvas *c = new TCanvas (Form("c%.0f",etaMin*100),"",400,400);
   beta->Fit("pol6","m");
   beta->SetXTitle("N_{Hits}");
   beta->SetYTitle(Form("#beta %.1f < #eta < %.1f",etaMin,etaMax));
   beta->SetAxisRange(0,1,"Y");
   beta->Draw("p");
   TF1 *f1 = beta->GetFunction("pol6");
   cout <<"#   if((eta)>="<<etaMin<<"&&(eta)<="<<etaMax<<") {"<<endl;
   cout <<"#      beta     = ("<<f1->GetParameter(0)<<")+"<<endl;
   cout <<"#                 ("<<f1->GetParameter(1)<<") * cnhit +"<<endl;
   cout <<"#                 ("<<f1->GetParameter(2)<<") * cnhit * cnhit +"<<endl;
   cout <<"#                 ("<<f1->GetParameter(3)<<") * cnhit * cnhit * cnhit +"<<endl;
   cout <<"#                 ("<<f1->GetParameter(4)<<") * cnhit * cnhit * cnhit * cnhit +"<<endl;
   cout <<"#                 ("<<f1->GetParameter(5)<<") * cnhit * cnhit * cnhit * cnhit * cnhit +"<<endl;
   cout <<"#                 ("<<f1->GetParameter(6)<<") * cnhit * cnhit * cnhit * cnhit * cnhit * cnhit;"<<endl;
//   cout <<"#                 ("<<f1->GetParameter(7)<<") * cnhit * cnhit * cnhit * cnhit * cnhit * cnhit * cnhit;"<<endl;
//   cout <<"#                 ("<<f1->GetParameter(8)<<") * cnhit * cnhit * cnhit * cnhit * cnhit * cnhit * cnhit * cnhit;"<<endl;
   cout <<"#   }    "<<endl;

   beta->Fit("fun","m");
   beta->SetXTitle("N_{Hits}");
   beta->SetYTitle(Form("#beta %.1f < #eta < %.1f",etaMin,etaMax));
   beta->SetAxisRange(0,1,"Y");
   beta->Draw("p");
   f1 = beta->GetFunction("fun");
   cout <<"@   if((eta)>="<<etaMin<<"&&(eta)<="<<etaMax<<") {"<<endl;
   cout <<"@      beta     = ("<<f1->GetParameter(0)<<")+"<<endl;
   cout <<"@                 ("<<f1->GetParameter(1)<<") / cnhit +"<<endl;
   cout <<"@                 ("<<f1->GetParameter(2)<<") * cnhit;"<<endl;
   cout <<"@   }    "<<endl;


}

void prepareHist(double etaMin, double etaMax)
{
   TFile *f = new TFile("output.root");
   TFile *g = new TFile("pseudoExp.root");
   prepareHist(f,g,etaMin,etaMax);
}
void prepareHist(TFile* f, TFile *g, double etaMin, double etaMax)
{
  gROOT->Reset();
//  gROOT->ProcessLine(".x rootlogon.C");
  gStyle->SetErrorX(0.);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  gStyle->SetPadLeftMargin(0.132);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.02);

  TH1::SetDefaultSumw2();


  TNtuple * ntmatched= dynamic_cast<TNtuple *>(f->Get("ntmatched"));
  TNtuple * ntInvMatched= dynamic_cast<TNtuple *>(g->Get("ntmatched"));
 
  TNtuple * ntcorr = new TNtuple("ntcorr","","nhad:ntrt:nstrt:vz:nvtx:vntrk");


  for (int i=0;i<200;i++) {
     h1[i] = new TH1F(Form("h1hit%d",i),Form("Everything;D;#_{pixel pairs}/%.2f",1/(double)nBins),nBins,0,2); 
     h2[i] = new TH1F(Form("h2hit%d",i),"Signal",nBins,0,2);
     h3[i] = new TH1F(Form("h3hit%d",i),"Background",nBins,0,2);
     h4[i] = new TH1F(Form("h4hit%d",i),"Normalized Reproduced Background",nBins,0,2);
     h5[i] = new TH1F(Form("h5hit%d",i),"",nBins,0,2);
     h6[i] = new TH1F(Form("h6hit%d",i),"Reproduced Background Subtracted",nBins,0,2); 
  }
  
  Float_t matchedeta1;
  Float_t matchedeta2;
  Float_t matchedinveta2;
  Float_t matchedinvphi2;
  Float_t signalcheck;
  Float_t matchedLayer1hits;
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

  // Tracklet Ntuple
  Float_t evtId;

  ntmatched->SetBranchAddress("eta1",&matchedeta1);
  ntmatched->SetBranchAddress("phi1",&matchedphi1);
  ntmatched->SetBranchAddress("matchedeta",&matchedeta2);
  ntmatched->SetBranchAddress("matchedphi",&matchedphi2);
  ntmatched->SetBranchAddress("signalCheck",&signalcheck);
  ntmatched->SetBranchAddress("deta",&deta);
  ntmatched->SetBranchAddress("dphi",&dphi);
  ntmatched->SetBranchAddress("evtid",&evtId);
  ntmatched->SetBranchAddress("nhit1",&matchedLayer1hits);

  Float_t evtId2;

  ntInvMatched->SetBranchAddress("eta1",&inveta1);
  ntInvMatched->SetBranchAddress("phi1",&invphi1);
  ntInvMatched->SetBranchAddress("matchedeta",&inveta2);
  ntInvMatched->SetBranchAddress("matchedphi",&invphi2);
  ntInvMatched->SetBranchAddress("deta",&invdeta);
  ntInvMatched->SetBranchAddress("dphi",&invdphi);
  ntInvMatched->SetBranchAddress("nhit1",&invlayer1hits);
  ntInvMatched->SetBranchAddress("evtid",&evtId2);
  
  int evts=0;
  int npass=0;
 
  int matchedentries = ntmatched->GetEntries();

  // Get the first eventId  
  for(int i = 0; i<matchedentries;i++){
    ntmatched->GetEntry(i);
    
//    float dR= sqrt(deta*deta+dphi*dphi);
    float dR= sqrt(deta*deta);
//    if(fabs(dphi)>0.05) continue;
    if((matchedeta1)>=etaMax) continue;
    if((matchedeta1)<=etaMin) continue;
    
    int bin = (int)matchedLayer1hits - (int) matchedLayer1hits % di;
    h1[bin]->Fill(fabs(dR));
    h6[bin]->Fill(fabs(dR));

    if(signalcheck==1){
      h2[bin]->Fill(fabs(dR));
    }
    if(signalcheck==0){
      h3[bin]->Fill(fabs(dR));
    } 
  }

  // Get the first eventId
  npass=0;
  evts=0;
  int invmatchedentries = ntInvMatched->GetEntries();
  for(int i = 0; i<invmatchedentries;i++){
    ntInvMatched->GetEntry(i);
//    if(fabs(invdphi)>0.05) continue;
    if((inveta1)>=etaMax) continue;
    if((inveta1)<=etaMin) continue;

    float dR= sqrt(invdeta*invdeta);
//    float dR= sqrt(invdeta*invdeta+invdphi*invdphi);
    int bin = (int)invlayer1hits - (int) invlayer1hits % di;

    h4[bin]->Fill(dR);
  }

  for (int i=0; i<maxHit;i++) {
     formatHist(h1[i],1,1);
     formatHist(h2[i],2,1);
     formatHist(h3[i],3,1);
     formatHist(h4[i],4,1);
     formatHist(h6[i],6,1);
  }
  
  cout <<"Histograms prepared!"<<endl;
}

double* getBeta(TFile* f, double etaMin, double etaMax, int MinHit, int MaxHit)
{

  //// Normalization of background

  Float_t nSideBand = h1[MinHit]->Integral((int)(normRange*nBins),nBins);
  Float_t nSBErr = sqrt(nSideBand);
  Float_t nReproducedSideBand = (h4[MinHit]->Integral((int)(normRange*nBins),nBins));
  Float_t nRSBErr = sqrt(nReproducedSideBand);
  Float_t sc = (nSideBand/nReproducedSideBand);
  Float_t scerr = sc * sqrt((nRSBErr/nReproducedSideBand)*(nRSBErr/nReproducedSideBand)
                           +(nSBErr/nSideBand)*(nSBErr /nSideBand));

  cout<<"nSideBand: "<<nSideBand<<" +- "<<nSBErr<<endl;
  cout<<"nReproducedSideBand: "<<nReproducedSideBand<<" +- "<<nRSBErr<<endl;

  cout<<"background normalization: "<<sc<<endl;
  h4[MinHit]->Scale(sc);
  h6[MinHit]->Add(h4[MinHit],-1);
  
  //// Determination of correction factor beta
  double beta=0,betaErr=100 ;
  
  
  if (nSideBand==0||nReproducedSideBand==0) {
     beta=0;
     betaErr=100;
  } else {
     Float_t nBck = (h4[MinHit]->Integral(0,(int)(deltaCut*nBins)))/sc;
     Float_t nTotal = (h1[MinHit]->Integral(0,(int)(deltaCut*nBins)));
     Float_t nBckErr = sqrt(1/nBck+(scerr/sc)*(scerr/sc))*nBck*sc;
     Float_t nTErr = sqrt(nTotal);
     cout<<"nBck: "<<nBck<<" +- "<<nBckErr<<endl;
     cout<<"nTotal: "<<nTotal<<" +- "<<nTErr<<endl;
     if (nTotal!=0&&nBck!=0) {
        beta=(nBck*sc)/(nTotal);
        betaErr=sqrt((nBckErr/nBck/sc)*(nBckErr/nBck/sc)+(nTErr/nTotal)*(nTErr/nTotal))*beta;
     }
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
  h->SetMarkerSize(0.7);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

}

void saveCanvas(TCanvas* c, int date){
  c->Write();
  c->Draw();
}

void analyze_betas(const char* infile)
{
   for (double i = -3; i<3; i+=0.5) {
      analyze_beta(infile,i,i+0.5);
   }   
}

void analyze(TH1F* h, TH1F* g,int bin)
{
    double n1=h->Integral(bin,10000);
    double n2=g->Integral(bin,10000);
    TH1F *x = (TH1F* )g->Clone();
    TH1F *y = (TH1F* )h->Clone();
    x->Scale(1./n2*n1);
    h->Draw();
    x->Draw("same");
    y->Add(x,-1);
    y->SetMarkerColor(2);
    y->SetLineColor(2);
    y->Draw("same");
    TLegend * leg1 = new TLegend(0.25,0.66,0.56,0.84);
    leg1->SetFillStyle(0);  
    leg1->SetFillColor(0); 
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h,"Everything","l");
    leg1->AddEntry(g,"Reproduced Combinatorial Background","l");
    leg1->AddEntry(y,"Reproduced Combinatorial Background Subtracted","l");

    leg1->Draw();

}
