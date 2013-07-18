#include "TFile.h"
#include "TString.h"
#include "TH3.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;
TH1D* momResHist(TFile*,int);

void momres() {

  int color[12]={kViolet+2,kBlue,kAzure+6,
		 kGreen-3,kOrange-5,kOrange-3,
		 kOrange+4,kRed-3,kRed,kRed+2,1,1};

  TString dir = "/home/sungho/sctch101/mc/7TeV/crab/trkhistsMC_june09";

  vector<TFile*> f;
  f.push_back(new TFile(Form("%s/TrkHistMC_june09_qcdMB.root",dir.Data())));
  f.push_back(new TFile(Form("%s/TrkHistMC_june09_qcdPt15.root",dir.Data())));
  f.push_back(new TFile(Form("%s/TrkHistMC_june09_qcdPt30v3.root",dir.Data())));
  f.push_back(new TFile(Form("%s/TrkHistMC_june09_qcdPt80v3.root",dir.Data())));
  f.push_back(new TFile(Form("%s/TrkHistMC_june09_qcdPt170v3.root",dir.Data())));
  f.push_back(new TFile(Form("%s/TrkHistMC_june09_qcdPt470.root",dir.Data())));

  vector<TH1D*> h;
  for(unsigned i=0; i<f.size(); i++) {
    h.push_back(momResHist(f[i],i));
    h[i]->SetMarkerColor(color[i]);
    for(int j=1; j<=h[i]->GetNbinsX(); j++) {
      if(h[i]->GetBinError(j)>0.002) {
	h[i]->SetBinContent(j,0.0);
	h[i]->SetBinError(j,0.0);
      }
    }
  }

  TCanvas *c3 = new TCanvas("c3","c3",500,500);
  c3->cd();
  for(unsigned i=0; i<h.size(); i++) {
    if(i==0) {h[i]->GetXaxis()->SetRangeUser(0.5,160); h[i]->Draw("p"); }
    else h[i]->Draw("psame");
  }


  TLegend *leg = new TLegend(0.6,0.6,0.85,0.85,"MC samples");
  leg->SetFillColor(0);
  leg->AddEntry(h[0],"MB","lp");
  leg->AddEntry(h[1],"Pt15","lp");
  leg->AddEntry(h[2],"Pt30","lp");
  leg->AddEntry(h[3],"Pt80","lp");
  leg->AddEntry(h[4],"Pt170","lp");
  leg->AddEntry(h[5],"Pt470","lp");
  leg->Draw();


  //gPad->SetLogx();

}

TH1D* momResHist(TFile* f, int i) {

  Double_t xbins[33] = {0,0.2,0.4,0.6,0.8,
			1.0,1.6,2.2,2.8,3.4,
			4.4,5.4,6.0,8.0,10.0,
			12.0,14.0,18.0,22.0,26.0,
			30.0,36.0,42.0,50.0,55.0,
			60.0,70.0,80.0,100.0,120.0,
			140.0,160.0,200.0};

  TH3F *hRes3D = (TH3F*) f->Get("trkEffAnalyzer/hresStoR3D");
  
  double feta = 2.4;
  int binMaxEta = hRes3D->GetXaxis()->FindBin(feta);
  int binMinEta = hRes3D->GetXaxis()->FindBin(-1.0*feta);

  hRes3D->GetXaxis()->SetRange(binMinEta,binMaxEta);
  TH2D *hRes2D = (TH2D*) hRes3D->Project3D("zy");  // zy: rec-> y axis, sim -> x axis
  hRes2D->SetName("hRes2D");

  hRes2D->FitSlicesY(0,0,-1,10);

  TH1D *h1 = (TH1D*)gDirectory->Get("hRes2D_1");
  TH1D *h2 = (TH1D*)gDirectory->Get("hRes2D_2");

  h1->Rebin(32,Form("h1new%d",i),xbins);
  h2->Rebin(32,Form("h2new%d",i),xbins);
  TH1D* h1new = (TH1D*) gDirectory->Get(Form("h1new%d",i));
  TH1D* h2new = (TH1D*) gDirectory->Get(Form("h2new%d",i));


  TH1D *momres = (TH1D*) h2new->Clone("momres");
  momres->SetName(Form("momres_%d",i));
  momres->SetMaximum(0.15); momres->SetMinimum(0.0);
  momres->SetTitle("Momentum resolution;p_{T} [GeV/c]");
  momres->SetStats(0);
  momres->SetMarkerStyle(20);
  momres->Divide(h1new);

  return momres;

}
