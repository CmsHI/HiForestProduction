#include "TF2.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

const double small = 1e-3;
vector<double> ptBins;
void define_variable_pt_bins();

void testMomRes(bool useRealSmearing=true, float momres=0.03, bool doVariableBins=false) {

  // parameterization of dN/dpt
  TCanvas * cDnDpt = new TCanvas("cDnDpt","cDnDpt",500,500);
  cDnDpt->SetLogy();
  TF1 *fDnDpt = new TF1("dndpt","[0]*([5]+(x/[1])+(pow(x,[2])/[3]))^[4]",0,250);
  fDnDpt->SetParameters(2.26959e+08,3.50661,1.56865,3.12002,-3.85260,6.86711e-01);
  fDnDpt->SetMinimum(1e-5);
  fDnDpt->SetMaximum(1e8);
  fDnDpt->Draw();
  // 1st iteration
  TF1 * fResCorr = new TF1("fResCorr","[0]+[1]*x+[2]*x^2",0,250);
  fResCorr->SetParameters(1.00379,2.62093e-05,7.21870e-07);
  fResCorr->Draw("same");
  fDnDpt = new TF1("dndptCorr0","dndpt/fResCorr",0,250);
  fDnDpt->SetLineColor(kRed);
  fDnDpt->Draw("same");

  // Either fill true vs. reco momentum by simple Gaussian with width from function argument
  TF2 *fMomRes = new TF2("fMomRes","[0]*exp(-0.5*(1-y/x)*(1-y/x)/[1]/[1])",0,250,0,250);
  fMomRes->SetParameters(1,momres);
  TH2D *hMomRes=0;
  if(!useRealSmearing) {
    if(!doVariableBins) {
      hMomRes = new TH2D("hMomRes","smearing factor",5000,0,250,5000,0,250);
    } else {
      define_variable_pt_bins();
      hMomRes = new TH2D("hMomRes","smearing factor", ptBins.size()-1, &ptBins[0], ptBins.size()-1, &ptBins[0]);
    }
    hMomRes->FillRandom("fMomRes",1e+06);
    hMomRes->Sumw2();
  }
  
  // ... or get true vs. reco momentum from QCD_Pt470 MC
  TString dir = "../root_files/trkhistsMC_june09";
  TFile *f = new TFile(Form("%s/TrkHistMC_june09_qcdPt470.root",dir.Data()));
  TH3F *hRes3D = (TH3F*) f->Get("trkEffAnalyzer/hresStoR3D");
  double feta = 2.4;
  int binMaxEta = hRes3D->GetXaxis()->FindBin(feta);
  int binMinEta = hRes3D->GetXaxis()->FindBin(-1.0*feta);
  hRes3D->GetXaxis()->SetRange(binMinEta,binMaxEta);
  if(useRealSmearing) {
    hMomRes = (TH2D*) hRes3D->Project3D("zy");  // zy: rec-> y axis, sim -> x axis
    hMomRes->SetName("hMomRes");
  }
  
  TH2D *hNormMomRes = (TH2D*) hMomRes->Clone("hNormMomRes");
  hNormMomRes->SetTitle("normalized smearing factor");
  TH2D *hWtMomRes = (TH2D*) hMomRes->Clone("hWtMomRes");
  hWtMomRes->SetTitle("weighted smearing factor");
  
  for(unsigned ptbin=hMomRes->GetNbinsX(); ptbin>0; ptbin--) {
    
    // get pt at center of bin
    float pt = hMomRes->GetXaxis()->GetBinCenter(ptbin);
    
    // get integral of column in true pt
    float integ = hMomRes->Integral(ptbin,ptbin, 1, hMomRes->GetNbinsY());
    
    // get the value of the dN/dpt fit at the true pt value of the bin center
    float wt = fDnDpt->Eval(pt);
    
    if (integ>0.0) {
      for (unsigned resbin=hMomRes->GetNbinsY(); resbin>0; resbin--) {
	double cont  = hMomRes->GetCellContent(ptbin,resbin);
	double err  = hMomRes->GetCellError(ptbin,resbin);
	cont /= integ;
	err /= integ;
	
	// weight each cell so that every column of true pt has unit integral
	hNormMomRes->SetCellContent(ptbin, resbin, cont);
	hNormMomRes->SetCellError(ptbin, resbin, err);
	
	// weight the normalized histogram by the dN/dpt fit value at the center of the bin
	hWtMomRes->SetCellContent(ptbin, resbin, wt * cont);
	hWtMomRes->SetCellError(ptbin, resbin, wt * err);
	
      }
    }
    
  }

  // project weighted true vs. reco 2-d hist to get smeared dN/dpt
  TH1D *hDnDptSmeared = (TH1D*) hWtMomRes->ProjectionY("hDnDptSmeared",1,hMomRes->GetNbinsY(),"e");
  
  // unsmeared histogram filled randomly from fit
  TH1D *hDnDptUnsmeared = (TH1D*) hDnDptSmeared->Clone("hDnDptUnsmeared");
  hDnDptUnsmeared->Reset();
  hDnDptUnsmeared->SetTitle("unsmeared spectra");
  for(int i=1; i<=hDnDptUnsmeared->GetNbinsX(); i++) {
    hDnDptUnsmeared->SetBinContent(i,fDnDpt->Eval(hDnDptUnsmeared->GetBinCenter(i)));
  }
  
  // divide smeared vs. unsmeared spectra
  TH1D *hSmearCorr = (TH1D*) hDnDptSmeared->Clone("hSmearCorr");
  //hSmearCorr->SetTitle(Form("Correction for %2.0f%% momentum resolution smearing;p_{T} GeV/c",100*momres));
  hSmearCorr->SetTitle("Correction for momentum resolution smearing from QCD Pt470 sample ;p_{T} GeV/c");
  hSmearCorr->Divide(hDnDptUnsmeared);

  // rebin if you like, but be careful because the very lowest bins are a little crazy
  int nrebin=5;
  hSmearCorr->Rebin(nrebin); hSmearCorr->Scale(1.0/nrebin);

  // cosmetics
  TCanvas * cMomResCorr = new TCanvas("cMomResCorr","cMomResCorr",500,500);
  hSmearCorr->SetMaximum(1.1);
  hSmearCorr->SetMinimum(1.0);
  hSmearCorr->SetMarkerStyle(20);
  hSmearCorr->GetXaxis()->SetRangeUser(0,150);
  hSmearCorr->Draw("p");

  hSmearCorr->Fit("pol2","","",2,170);
  hSmearCorr->GetFunction("pol2")->SetLineColor(kRed);

  float ptMin=0,ptMax=170;
  TF1 * fResCorrUpper = new TF1("fResCorrUpper","[0]+[1]*x+[2]*x^2",ptMin,ptMax);
  fResCorrUpper->SetParameters(1.00379,6.58589e-05,1.3465e-06);
  fResCorrUpper->SetLineColor(kBlue);
  fResCorrUpper->Draw("same");

  TF1 * fResCorrLower = new TF1("fResCorrLower","[0]+[1]*x+[2]*x^2",ptMin,ptMax);
  fResCorrLower->SetParameters(1.00379,-2.58589e-05,5.3465e-07);
  fResCorrLower->SetLineColor(kGreen+2);
  fResCorrLower->Draw("same");
}


void define_variable_pt_bins() {

  // define variable pt bins
  double ptb;
  
  for(ptb =   0  ; ptb <   1.2-small; ptb +=  0.05) ptBins.push_back(ptb); // 24 bins
  for(ptb =   1.2; ptb <   2.4-small; ptb +=  0.1 ) ptBins.push_back(ptb); // 12 bins
  for(ptb =   2.4; ptb <   7.2-small; ptb +=  0.2 ) ptBins.push_back(ptb); // 24 bins
  for(ptb =   7.2; ptb <  13.2-small; ptb +=  0.5 ) ptBins.push_back(ptb); // 12 bins
  for(ptb =  13.2; ptb <  25.2-small; ptb +=  1.0 ) ptBins.push_back(ptb); // 12 bins
  for(ptb =  25.2; ptb <  61.2-small; ptb +=  3.0 ) ptBins.push_back(ptb); // 12 bins
  for(ptb =  61.2; ptb < 121.2-small; ptb +=  5.0 ) ptBins.push_back(ptb); // 12 bins
  for(ptb = 121.2; ptb < 241.2-small; ptb += 10.0 ) ptBins.push_back(ptb); // 12 bins

}
