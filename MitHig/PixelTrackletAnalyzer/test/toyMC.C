#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <vector>
#include <TF1.h>
#include <TRandom.h>

#include "tracklet.h"

#define PI 3.14159265358979

TFile *f;

void test(int nhit1);
void genLayer1(int nhit1,vector<double> &hits1);
void genLayer2(int nhit2,vector<double> &hits2);

bool verbose_=false;
bool useDeltaPhi_=false;
bool checkSecondLayer_=true;
bool compareDeltaR(Tracklet a,Tracklet b) { return fabs(a.dR2())<fabs(b.dR2());}
bool compareDeltaEta(Tracklet a,Tracklet b) {return fabs(a.deta())<fabs(b.deta());}
vector<Tracklet> cleanTracklets(vector<Tracklet> input, int matchNumber);
void toyMC(int run,int nevt);

void toyMCs(int start,int end)
{
   for (int i=start;i<=end;i++){
      toyMC(i,10000);   
   }
}
void toyMC(int run,int nevt)
{
    f = new TFile("output.root");

    // Random Seed
    gRandom->SetSeed( run );

    TH3F *nhits = (TH3F*) f->FindObjectAny("nhits");

    TFile *outf = new TFile (Form("exp-%05d.root",run),"recreate");
    TNtuple *ntmatched = new TNtuple("ntmatched","","eta1:matchedeta:phi1:matchedphi:deta:dphi:signalCheck:tid:r1id:r2id:evtid:nhit1:sid:ptype");
    TNtuple *ntHit1 = new TNtuple("ntHit1","","eta1:phi1:nhit1");
    TNtuple *ntHit2 = new TNtuple("ntHit2","","eta2:phi2:nhit1");

    for (int i=0;i<nevt;i++) {
       vector<double> hits1;
       vector<double> hits2;
       vector<Tracklet> protoTracklets;
       vector<Tracklet> recoTracklets;
       double mult,nhit1,nhit2;
       nhits->GetRandom3(mult,nhit1,nhit2);
       hits1.clear();
       hits2.clear();
       if (i% 1000 == 000 ) cout <<"Run: "<<run<<" Event "<<i<<" "<<(int)mult<<" "<<(int)nhit1<<" "<<hits1.size()<<" "<<(int)nhit2<<" "<<hits2.size()<<endl;
       genLayer1((int)nhit1,hits1);
       genLayer2((int)nhit2,hits2);
       
       for (int j=0;j<(int) hits1.size();j+=2)
       {
          ntHit1->Fill(hits1[j],hits1[j+1],mult);
          for (int k=0;k<(int)hits2.size();k+=2)
	  {
             ntHit2->Fill(hits2[k],hits2[k+1],mult);
	     Tracklet mytracklet(hits1[j],hits2[k],hits1[j+1],hits2[k+1]);
	     mytracklet.setIt1(j/2);
	     mytracklet.setIt2(k/2);
	     protoTracklets.push_back(mytracklet); 
	  }
       }
       recoTracklets = cleanTracklets(protoTracklets,0);
       for (int j=0;j<(int)recoTracklets.size();j++)
       {
           float var[100];
	   var[0] = recoTracklets[j].eta1();
	   var[1] = recoTracklets[j].eta2();
	   var[2] = recoTracklets[j].phi1();
	   var[3] = recoTracklets[j].phi2();
	   var[4] = recoTracklets[j].deta();
	   var[5] = recoTracklets[j].dphi();
	   var[6] = 0;
	   var[7] = recoTracklets[j].getId();
	   var[8] = recoTracklets[j].getId1();
	   var[9] = recoTracklets[j].getId2();
	   var[10] = i;
	   var[11] = (int)mult;
	   var[12] = recoTracklets[j].getSId();
	   var[13] = recoTracklets[j].getType();
	   ntmatched->Fill(var);
       }
    }

    ntmatched->Write();
    ntHit1->Write();
    ntHit2->Write();
    outf->Close();
}

void genLayer1(int nhit1,vector<double> &hits1)
{
    TH1F *h = (TH1F*) f->FindObjectAny(Form("dNdEtaHits1_%02d",nhit1));
    TH1F *h2 = (TH1F*) f->FindObjectAny(Form("dNdPhiHits1_%02d",nhit1));
    hits1.clear();
    double eta,phi;
    int nGen=0;
    for (int i=0;nGen<nhit1;i++){
       eta = h->GetRandom();
       phi = h2->GetRandom();
       if (fabs(eta)<1000) nGen++;
       hits1.push_back(eta); 
       hits1.push_back(phi); 
    }
}

void genLayer2(int nhit2,vector<double> &hits2)
{
    TH1F *h = (TH1F*) f->FindObjectAny(Form("dNdEtaHits2_%02d",nhit2));
    TH1F *h2 = (TH1F*) f->FindObjectAny(Form("dNdPhiHits2_%02d",nhit2));
    hits2.clear();
    double eta,phi;
    int nGen=0;
    for (int i=0;nGen<nhit2;i++){
       eta = h->GetRandom();
       phi = h2->GetRandom();
       if (fabs(eta)<1000) nGen++;
       hits2.push_back(eta); 
       hits2.push_back(phi); 
    }
}

void test(int nhit1)
{
    TH1F *h = (TH1F*) f->FindObjectAny(Form("dNdEtaHits1_%02d",nhit1));
    TH1F *g = new TH1F("toymc","",100,-3,3);

    for (int i=0;i<h->GetEntries();i++){
       g->Fill(h->GetRandom());
    }
    
    g->SetLineColor(2);
    h->Draw();
    g->Draw("same");
}

vector<Tracklet> cleanTracklets(vector<Tracklet> input, int matchNumber)
{
   vector<Tracklet> output;

   if(useDeltaPhi_)
      sort( input.begin() , input.end() , compareDeltaR);
   else
      sort( input.begin() , input.end() , compareDeltaEta);

   if (verbose_) {
      for (unsigned int i = 0; i < input.size(); i++)
      {
         cout <<input[i].deta()<<" "<<input[i].getIt1()<<" "<<input[i].getIt2()<<endl;
      }
   }

   int used1[1000];
   int used2[1000];

   for (int i=0;i<1000;i++) { 
      used1[i]=0;
      used2[i]=0;
   } 

   if (verbose_) cout<<"Printing Hits"<<endl;
   
   for (unsigned int i = 0; i < input.size(); i++){
      

      if(useDeltaPhi_)
	 if (verbose_) cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta R : "<<input[i].dR()<<endl;
      else
	 if (verbose_) cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta Eta : "<<input[i].deta()<<endl; 
      
      int i1=input[i].getIt1();
      int i2=input[i].getIt2();

      if (used1[i1]==0&&used2[i2]==matchNumber) {
	 Tracklet tmp = input[i];
	 output.push_back(tmp);
	 used1[i1]++;
	 if (checkSecondLayer_) used2[i2]++;
      }
      if (used1[i1]==0&&used2[i2]<matchNumber) {
	 if (checkSecondLayer_) used2[i2]++;
      }
   }
   if (verbose_) {
      cout <<"Output:"<<endl;
      for (unsigned int i = 0; i < output.size(); i++)
      {
         cout <<output[i].deta()<<" "<<output[i].getIt1()<<" "<<output[i].getIt2()<<endl;
      }
   }
   
   return output;
}
