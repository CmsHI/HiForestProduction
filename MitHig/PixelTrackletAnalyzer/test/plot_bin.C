#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "../../PixelTracklet/interface/TrackletCorrections.h"
using namespace std;

void plot_bin(int bin = 0, const char* histfile = "histograms.root", const char* corrfile = "corrections.root"){

  TFile* corrf = new TFile(corrfile,"read");
  TrackletCorrections* corr = new TrackletCorrections(corrf);

  corr->getBin(bin);

  double etaMax = corr->binEtaMax();
  double etaMin = corr->binEtaMin();
  double hitMin = corr->binHitMin();
  double hitMax = corr->binHitMax();
  double zMin = corr->binZMin();
  double zMax = corr->binZMax();

  double beta = corr->beta(bin);


  TFile* hists = new TFile(histfile,"read");
  TH1F* h1 = dynamic_cast<TH1F*>(hists->Get(Form("h1_%d",bin)));
  TH1F* h2 = dynamic_cast<TH1F*>(hists->Get(Form("h2_%d",bin)));
  TH1F* h3 = dynamic_cast<TH1F*>(hists->Get(Form("h3_%d",bin)));
  TH1F* h4 = dynamic_cast<TH1F*>(hists->Get(Form("h4_%d",bin)));
  TH1F* h5 = dynamic_cast<TH1F*>(hists->Get(Form("h5_%d",bin)));
  TH1F* h6 = dynamic_cast<TH1F*>(hists->Get(Form("h6_%d",bin)));


  TCanvas* c1 = new TCanvas("c1","Backgrounds",500,400);
  c1->SetLogy();

  h1->Draw();
  h3->Draw("same");
  h4->Draw("same");


  TLegend* leg1 = new TLegend(0.50,0.60,0.80,0.87);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(h1,"Everything","lp");
  //  leg1->AddEntry(h2,"Signal","lp");
  leg1->AddEntry(h3,"Background","lp");
  leg1->AddEntry(h4,"Normalized Reproduced Background","lp");
  //  leg1->AddEntry(h6,"Reproduced Background Subtracted","lp");
  leg1->AddEntry(h1,Form("%d < #_{hits} < %d",(int)hitMin,(int)hitMax),"");
  leg1->AddEntry(h1,Form("%0.1f < #eta < %0.1f",etaMin,etaMax),"");
  leg1->AddEntry(h1,Form("%d < z_{vtx} < %d",(int)zMin,(int)zMax),"");
  leg1->AddEntry(h1,Form("#beta = %0.4f",beta),"");

  leg1->Draw();


  cout<<"Bin : "<<bin<<endl;
  cout<<"      "<<endl;


  cout<<"Hit Min = "<<hitMin<<endl;
  cout<<"Hit Max = "<<hitMax<<endl;
  cout<<"      "<<endl;

  cout<<"Eta Min = "<<etaMin<<endl;
  cout<<"Eta Max = "<<etaMax<<endl;
  cout<<"      "<<endl;

  cout<<"Z Min   = "<<zMin<<endl;
  cout<<"Z Max   = "<<zMax<<endl;
  cout<<"      "<<endl;

  cout<<"beta    = "<<beta<<endl;

  c1->Print(Form("BackgroundRecoForBin%d.gif",bin));

}
