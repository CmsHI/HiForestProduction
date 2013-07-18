#include <iostream>
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
#include "TrackletData.h"
#include "TrackletCorrections.h"

using namespace std;

void formatHist(TH1* h, int col = 1, double norm = 1);
void saveCanvas(TCanvas* c, int date = 20080829);

void create_beta(const char* infile = //"/server/03a/yetkin/data/tracklet-vtxFlat-20081006.root",
		 "/server/03a/yetkin/data/pythia_mb_900GeV_vtxFlat_d20081001.hist.root", 
		 const char* outfile = "histograms.root"){

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

  int hitbins = 1;
  int etabins = 1;
  int zbins = 10;

  TrackletCorrections* corr = new TrackletCorrections(hitbins, etabins,zbins);

  corr->setMidEtaCut(0.1); // 0.1
  corr->setDeltaRCut(0.25);
  corr->setCPhi(1/43);
  corr->setNormDRMin(0.75);
  corr->setNormDRMax(2.75);

  corr->setHistBins(40);
  corr->setHistMax(5);

  corr->setHitMax(15);
  corr->setEtaMax(3);
  corr->setZMax(20);

  corr->start();

  int size  =  corr->size();
  cout<<"Corrections size : "<<size<<endl;

  TFile *f = new TFile(infile);
  TrackletData data(f,corr);

  TFile *of = new TFile(outfile,"recreate");
  for(int ibin = 0; ibin < size; ++ibin){
    corr->setBeta(ibin,data.getBeta(ibin,true));
  }

  of->Write();
  of->Close();

  corr->save("corrections20081007_z.root");

}

double TrackletData::getBeta(int bin,bool saveplots)
{

  /// Parameters

  corr->getBin(bin);

  //  double etaMin = corr->binEtaMin();
  //  double etaMax = corr->binEtaMax();

  double etaMax = corr->getEtaMax();
  double etaMin = -etaMax;

  double MinHit = corr->binHitMin();
  double MaxHit = corr->binHitMax();
  double MinZ = corr->binZMin();
  double MaxZ = corr->binZMax();

  double midEtaCut = corr->getMidEtaCut();
  
  cout<<"Bin : "<<bin<<endl;
  cout<<"      "<<endl;

  cout<<"Hit Min = "<<MinHit<<endl;
  cout<<"Hit Max = "<<MaxHit<<endl;
  cout<<"      "<<endl;

  cout<<"Eta Min = "<<etaMin<<endl;
  cout<<"Eta Max = "<<etaMax<<endl;
  cout<<"      "<<endl;

  cout<<"Z Min   = "<<MinZ<<endl;
  cout<<"Z Max   = "<<MaxZ<<endl;
  cout<<"      "<<endl;

  double deltaCut = corr->getDeltaRCut();
  int etaBins = corr->getEtaBins();
  int nBins = corr->getHistBins(); //20;
  double max = corr->getHistMax(); //20;

  double normMin = corr->getNormDRMin();
  double normMax = corr->getNormDRMax();

  TH1F * h1 = new TH1F(Form("h1_%d",counter),Form("Everything;D;#_{pixel pairs}/event/%.2f",1/(double)nBins),nBins,0,max);
  TH1F * h2 = new TH1F(Form("h2_%d",counter),"Signal",nBins,0,max);
  TH1F * h3 = new TH1F(Form("h3_%d",counter),"Background",nBins,0,max);
  TH1F * h4 = new TH1F(Form("h4_%d",counter),"Normalized Reproduced Background",nBins,0,max);
  TH1F * h5 = new TH1F(Form("h5_%d",counter),"Un-Normalized Reproduced Background",nBins,0,max);
  TH1F * h6 = new TH1F(Form("h6_%d",counter),"Reproduced Background Subtracted",nBins,0,max);
  TH1F * h7 = new TH1F(Form("h7_%d",counter),"",nBins,0,max);
  TH1F * h8 = new TH1F(Form("h8_%d",counter),"",nBins,0,max);
  TH1F * h9 = new TH1F(Form("h9_%d",counter),"",nBins,0,max);
  TH1F * h10 = new TH1F(Form("h10_%d",counter),"",nBins,0,max);

  TProfile * dNdEtaHadron = new TProfile(Form("dNdEtaHadron_%d",counter),"",32,-2.1,2.1);
  TProfile * dNdEtaLepton = new TProfile(Form("dNdEtaLepton_%d",counter),"",32,-2.1,2.1);
  TProfile * dNdEtaTracklet = new TProfile(Form("dNdEtaTracklet_%d",counter),"",32,-2.1,2.1);

  // TH2D * dNdEtaHadron = new TH2D(Form("dNdEtaHadron_%d",counter),"",32,-2.1,2.1,50,0,5);
  // TH2D * dNdEtaLepton = new TH2D(Form("dNdEtaLepton_%d",counter),"",32,-2.1,2.1,50,0,5);

  // TH2D * correlation = new TH2D(Form("correlation_%d",counter),"; #_{hadrons}; #_{tracklets}",100,0,50,100,0,50);

  counter++;

  int nevents = ntgen->GetEntries();
  //  int partentries = ntparticle->GetEntries();
  int matchedentries = ntmatched->GetEntries();

  ntevent->GetEntry(0);
  ntvertex->GetEntry(0);

  int nPassCounter=0;
  int oldEvtId=0;
  int evts=0;

  double layer1hits = 0;
  for(int ig =0; ig < 12; ++ig) layer1hits += hits[ig];

  for(int i = 0; i<matchedentries;i++){
    ntmatched->GetEntry(i);

    if (evtId != (int)eventEvtId) {
      while (evtId!=(int)eventEvtId && evts < ntevent->GetEntries()){
        // Check if the event is empty

	//	cout <<evts<<" "<<eventEvtId<<" "<<evtId<<" "<<nPassCounter<<endl;
        nPassCounter = 0;
        oldEvtId++;
        evts++;
        ntevent->GetEntry(evts);
        ntvertex->GetEntry(evts);

	layer1hits = 0;
	for(int ig =0; ig < 8; ++ig) layer1hits += hits[ig];
      }
    }

    if(matchedeta1>=etaMax) continue;
    if(matchedeta2>=etaMax) continue;
    if(matchedeta1<etaMin) continue;
    if(matchedeta2<etaMin) continue;

    if(z>=MaxZ) continue;
    if(z<MinZ) continue;

    if(layer1hits>=MaxHit && layer1hits<corr->getHitMax()) continue;
    if(layer1hits<MinHit) continue;    

    if(fabs(matchedeta1)<midEtaCut) continue;

    //    float dR= sqrt(deta*deta+dphi*dphi/43./43.);
    // This is not the standard delta R - It is with the phi normalized with corr->getCPhi();
    double dR = deltaR(matchedeta1,matchedphi1,matchedeta2,matchedphi2);

    h1->Fill(dR);
    h6->Fill(dR);
    if(signalcheck==1){
      h2->Fill(dR);
    }
    if(signalcheck==0){
      h3->Fill(dR);
    } 
  }


  ntevent->GetEntry(0);
  ntvertex->GetEntry(0);
  int nInvPassCounter=0;
  int oldInvEvtId=0;
  evts=0;
  double layer1InvHits = 0;
  for(int ig =0; ig < 8; ++ig) layer1InvHits += hits[ig];

  int invmatchedentries = ntInvMatched->GetEntries();
  for(int i = 0; i<invmatchedentries;i++){
    ntInvMatched->GetEntry(i);

    if (evtInvId != (int)eventEvtId) {
      while (evtInvId!=(int)eventEvtId && evts < ntevent->GetEntries()){
        // Check if the event is empty

	//cout <<evts<<" "<<eventEvtId<<" "<<evtInvId<<" "<<nInvPassCounter<<endl;
        nInvPassCounter = 0;
        oldInvEvtId++;
        evts++;
        ntevent->GetEntry(evts);
        ntvertex->GetEntry(evts);

        layer1InvHits = 0;
        for(int ig =0; ig < 8; ++ig) layer1InvHits += hits[ig];
      }
    }


    if(z>=MaxZ) continue;
    if(z<MinZ) continue;

    if(inveta2>=etaMax) continue;
    if(inveta2<etaMin) continue;
    if(inveta1>=etaMax) continue;
    if(inveta1<etaMin) continue;

    //    if(layer1InvHits>=MaxHit && layer1InvHits<corr->getHitMax()) continue;
    //    if(layer1InvHits<MinHit) continue;

    if(fabs(inveta1)<midEtaCut) continue;
    
    //    float dR= sqrt(invdeta*invdeta+invdphi*invdphi);
    // This is not the standard delta R - It is with the phi normalized with corr->getCPhi();
    double dR = deltaR(inveta1,invphi1,inveta2,invphi2);
    h4->Fill(dR);


    //    h4->Fill(fabs(invdeta));

  }

  /*  
  for(int i = 0; i<partentries;i++){
    ntparticle->GetEntry(i);
    if(eta1>etaMax || eta2>etaMax) continue;
    if(eta1<=etaMin || eta2<=etaMin) continue;

    if(charge==0) continue;
    h7->Fill(fabs(eta1-eta2));
  }
  */


  if(saveplots){
  formatHist(h1,1,nevents);
  formatHist(h2,2,nevents);
  formatHist(h3,3,nevents);
  formatHist(h4,4,nevents);
  formatHist(h5,5,nevents);
  formatHist(h6,6,nevents);
  }

  //// Normalization of background

  Float_t sc = ((h1->Integral((int)(normMin*nBins/max),(int)(normMax*nBins/max),"width"))/(h4->Integral((int)(normMin*nBins/max),(int)(normMax*nBins/max),"width")));

  cout<<"background normalization: "<<sc<<endl;
  h5->Add(h4);
  h4->Scale(sc);
  h6->Add(h4,-1);

  //// Determination of correction factor beta
  double beta = (h4->Integral(0,(int)(deltaCut*nBins/max),"width")/(h1->Integral(0,(int)(deltaCut*nBins/max),"width")));
  cout<<"beta: "<<beta<<endl;

  cout<<"      "<<endl;
  cout<<"      "<<endl;

  cout<<"-------------"<<endl;
  cout<<"      "<<endl;
  cout<<"      "<<endl;

  cout<<"      "<<endl;


  return beta;

}

double TrackletData::getAlpha(int bin, bool saveplots){
  return 0;
}






