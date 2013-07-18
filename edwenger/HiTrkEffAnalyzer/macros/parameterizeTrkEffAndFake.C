//--------------------------------------------------------------------
//  Andre S. Yoon, December 03 2010
//
//  A root macro that does the parameterization of tracking eff/fake
//  for a given eta, pt, z bins (where z can be centralities, and etc)
//  At the end, summary of all parameters for each eta, pt, z ranges
//  is given. 
//
//  * definition of efficiency = algo eff. x acceptance
//--------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TStyle.h>
#include "TError.h"

#include "TFile.h"
#include "TCanvas.h"

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TF1.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TObjArray.h"

#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TChain.h"
#include "TGraph.h"

#include "TKey.h"
#include "TMath.h"
#include "TRandom.h"

#include "TLatex.h"

using namespace std;

#endif

TFile *fMC=0;

TH3F *hSim3D=0;
TH3F *hAcc3D=0;
TH3F *hEff3D=0;

TH3F *hRec3D=0;
TH3F *hFak3D=0;

TH3F *rAcc3D=0;
TH3F *rEff3D=0;
TH3F *rFak3D=0;

vector<TH1D*> hSim;
vector<TH1D*> hRec;

vector<TH1D*> hEff;
vector<TH1D*> hFak;

vector<vector<TH1D*> > hSimVV;
vector<vector<TH1D*> > hRecVV;

vector<vector<TH1D*> > hEffVV;
vector<vector<TH1D*> > hFakVV;

vector<TGraphAsymmErrors*> tgEff;
vector<TGraphAsymmErrors*> tgFak;

vector<vector<TGraphAsymmErrors*> > tgEffVV;
vector<vector<TGraphAsymmErrors*> > tgFakVV;

vector<TCanvas*> tcan;
vector<TCanvas*> tcan2;

TF1* fit1=0;
TF1* fit2=0;
TF1* fit3=0;
TF1* fit_tot=0;

TF1* fitfake1=0;
TF1* fitfake2=0;
TF1* fitfake3=0;

vector<TF1*> fitV1;
vector<TF1*> fitV2;
vector<TF1*> fitV3;

vector<double> etaBins;
vector<double> ptBins;
vector<double> zvarBins;


bool drawing = true;
bool fulltrk = true;
float ptmax = 0;
float fit1_min = 0, fit1_max = 0;
float fit2_min = 0, fit2_max = 0;
float fit3_min = 0, fit3_max = 0;

TString infile;
ofstream fdata;

//--------------------------------------------------
void prepareThings();
void drawText(const char *text, float xp, float yp);
TH1D *hackedProjectionY(TH3F *me, const char *name, int ixmin, int ixmax, int izmin, int izmax);
void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns,
                          const Int_t rows, const Float_t leftOffset=0.,
                          const Float_t bottomOffset=0.,
                          const Float_t leftMargin=0.2,
                          const Float_t bottomMargin=0.2,
                          const Float_t edge=0.05);
//--------------------------------------------------

void prepareThings(){

   // input file
   if(fulltrk) infile = "validation3D_bass_hiGoodTracks";
   else infile = "validation3D_bass_hiConformalPixelTracks";

   // output data file
   fdata.open(Form("%s_fulltrk%d.dat",infile.Data(),fulltrk));

   // eta range 1 (-2.4 ~ -1.0)
   static float etaMin1   = -2.4;
   static float etaMax1   = -1.0;
   static float etaWidth1 =  0.2;
   
   for(double eta = etaMin1; eta < etaMax1 + etaWidth1/2; eta += etaWidth1)
      etaBins.push_back(eta);
   
   // eta range 2 (-1.0 ~ 1.0)
   etaBins.push_back(-0.6);
   etaBins.push_back(-0.2);
   etaBins.push_back(0.2);
   etaBins.push_back(0.6);
   
  // eta range 1 (1.0 ~ 2.4) 
  static float etaMin2   = 1.0;
  static float etaMax2   = 2.4;
  static float etaWidth2 =  0.2;

  for(double eta = etaMin2; eta < etaMax2 + etaWidth2/2; eta += etaWidth2)
    etaBins.push_back(eta);
  

  static float zvarMin   =  0.0;
  static float zvarMax   =  40.0;
  static float zvarWidth =  20;

  //for(double zvar = zvarMin; zvar < zvarMax + zvarWidth/2; zvar += zvarWidth)
  //zvarBins.push_back(zvar);

  zvarBins.push_back(0);
  zvarBins.push_back(4);
  zvarBins.push_back(12);
  zvarBins.push_back(40);

  if(fulltrk) ptmax = 10.0;
  else ptmax = 3.0;
  drawing = true;
  
  fit1_min = 0.2, fit1_max = 1.5;
  fit2_min = 1.5, fit2_max = 1.80;
  fit3_min = 1.80, fit3_max = 100;

  
  if(fulltrk){
     // fit function hiGoodTracks
     fit1 = new TF1("fit1","[0]*pow(x,[3])/(1+exp([1]*(x+[2]))) + [4]*pow(x,[5])",0.7,200);
     fit1->SetParLimits(3,0.1,1.0);
     fit1->SetParLimits(4,-6.0E-2,-1.0E-3);
     fit1->SetParLimits(5,0.2,1.6);
     
     fitfake1 = new TF1("fitfake1","[0]*pow(x,[3])/(1+exp([1]*(x+[2]))) + [4]*pow(x,[5])",0.7,200);  
     fitfake1->SetParLimits(4,8.5E-1,1.0);
     fitfake1->SetParLimits(5,-1.0,1.0); // -1.0 can be lowered 
  }else{
     // fit function hiConformalPixelTracks
     fit1 = new TF1("fit1","[0]*pow(x,[3])/(1+exp([1]*(x+[2]))) + [4]*pow(x,[5])",0.2,200);
     fit1->SetParLimits(3,0.1,1.0);
     fit1->SetParLimits(4,-6.0E-2,-1.0E-3);
     fit1->SetParLimits(5,0.2,1.6);

     fitfake1 = new TF1("fitfake1","[0]*pow(x,[3])/(1+exp([1]*(x+[2]))) + [4]*pow(x,[5])",0.2,200);  
     fitfake1->SetParLimits(4,8.5E-1,1.0);
     fitfake1->SetParLimits(5,-1E-2,1.0);
  }


}

void parameterizeTrkEffAndFake(){

   // prepare binning, etc 
   prepareThings();

   fMC = new TFile(Form("./root_files/%s.root",infile.Data()));
   cout<<"Input file loaded = "<<infile.Data()<<endl;

  // sim-to-reco hists
  hSim3D = (TH3F*) fMC->Get("hitrkPixelEffAnalyzer/hsim3D"); 
  hAcc3D = (TH3F*) fMC->Get("hitrkPixelEffAnalyzer/hacc3D"); 
  hEff3D = (TH3F*) fMC->Get("hitrkPixelEffAnalyzer/heff3D"); 

  // reco-to-sim hists
  hRec3D = (TH3F*) fMC->Get("hitrkPixelEffAnalyzer/hrec3D"); 
  hFak3D = (TH3F*) fMC->Get("hitrkPixelEffAnalyzer/hfak3D"); 

  //loop over and slice the 3D hists in bins of eta, cent
  int counter = 0;

  cout<<"[Loop begins]==============================================" <<endl;
  for(int i=0;i<zvarBins.size()-1;i++){  // -1 because bin ranges has one more number 
    
    cout<<"[zvar slicing].................."<<endl;
    cout<<" number of bins = "<<zvarBins.size()-1<<endl;

    float small = 0.01;
    int zvar_min = 0.0;
    int zvar_max = 0.0;

    zvar_min  = hSim3D->GetZaxis()->FindBin(zvarBins[i]+small);
    zvar_max  = hSim3D->GetZaxis()->FindBin(zvarBins[i+1]-small);

    cout<<"  input zvar ranges: "<<zvarBins[i]<<" to "<<zvarBins[i+1]<<endl;
    cout<<"    found zvar bin ranges: "<<zvar_min<<" to "<<zvar_max<<endl; 
    cout<<"    found zvar low edge: "<<hSim3D->GetZaxis()->GetBinLowEdge(zvar_min)
	<<" up edge: "<<hSim3D->GetZaxis()->GetBinUpEdge(zvar_max)<<endl;
  
    cout<<"  [eta slicing].................."<<endl;
    cout<<"   number of bins = "<<etaBins.size()-1<<endl;

    // first push_back ith row with empty vector
    hSimVV.push_back(hSim), hEffVV.push_back(hEff);

    hRecVV.push_back(hRec), hFakVV.push_back(hFak);

    tgEffVV.push_back(tgEff), tgFakVV.push_back(tgFak); 

    for(int j=0;j<etaBins.size()-1;j++){

      int eta_min = 0.0;
      int eta_max = 0.0;

      eta_min  = hSim3D->GetXaxis()->FindBin(etaBins[j]+small);
      eta_max  = hSim3D->GetXaxis()->FindBin(etaBins[j+1]-small);

      cout<<"   input eta ranges: "<<etaBins[j]<<" to "<<etaBins[j+1]<<endl;
      cout<<"     found eta bin ranges: "<<eta_min<<" to "<<eta_max<<endl;
      cout<<"     found eta low edge: "<<hSim3D->GetXaxis()->GetBinLowEdge(eta_min)
	  <<" up edge: "<<hSim3D->GetXaxis()->GetBinUpEdge(eta_max)<<endl;

      hSimVV[i].push_back(hackedProjectionY(hSim3D,Form("hsim_proj_%d",100*i+j),eta_min,eta_max,zvar_min,zvar_max));
      hEffVV[i].push_back(hackedProjectionY(hEff3D,Form("heff_proj_%d",100*i+j),eta_min,eta_max,zvar_min,zvar_max));
      tgEffVV[i].push_back(new TGraphAsymmErrors()), tgEffVV[i][j]->SetName(Form("tgEff_%d",100*i+j));

      hRecVV[i].push_back(hackedProjectionY(hRec3D,Form("hrec_proj_%d",100*i+j),eta_min,eta_max,zvar_min,zvar_max));
      hFakVV[i].push_back(hackedProjectionY(hFak3D,Form("hfak_proj_%d",100*i+j),eta_min,eta_max,zvar_min,zvar_max));
      tgFakVV[i].push_back(new TGraphAsymmErrors()), tgFakVV[i][j]->SetName(Form("tgFak_%d",100*i+j));

      counter++;
    }

    cout<<" counter = "<<counter<<endl;
  }
  cout<<"[Loop ended]==============================================" <<endl;
  

   // BayesDivide
   cout<<"[BayesDivide] size of tgEffW (number of bins in zvar) = "<<tgEffVV.size()<<endl;

   for(int i=0; i<tgEffVV.size();i++){
     cout<<"[BayesDivide] size of tgEffW[i] (num of bins in eta) = "<<tgEffVV[i].size()<<endl;
     for(int j=0; j<tgEffVV[i].size();j++){

	//efficiency
	tgEffVV[i][j]->BayesDivide(hEffVV[i][j],hSimVV[i][j]);
	tgEffVV[i][j]->SetMarkerStyle(25), tgEffVV[i][j]->SetLineStyle(2);
	tgEffVV[i][j]->SetLineColor(2), tgEffVV[i][j]->SetMarkerColor(2);
	tgEffVV[i][j]->SetMarkerSize(1.0);
	
	//fake rate
	tgFakVV[i][j]->BayesDivide(hFakVV[i][j],hRecVV[i][j]);
        tgFakVV[i][j]->SetMarkerStyle(25), tgFakVV[i][j]->SetLineStyle(2);
        tgFakVV[i][j]->SetLineColor(4), tgFakVV[i][j]->SetMarkerColor(4);
	tgFakVV[i][j]->SetMarkerSize(1.0);
     }

   }

  // Dum hist
  TH1F *hDum = new TH1F("hDum","",100,0,ptmax);
  hDum->SetStats(0), hDum->SetMinimum(0.0), hDum->SetMaximum(1.0);
  hDum->GetXaxis()->SetTitle("p_{T}"), hDum->GetYaxis()->SetTitle("Eff (or Fake)");
  hDum->GetXaxis()->CenterTitle(), hDum->GetYaxis()->CenterTitle();

  hDum->GetXaxis()->SetLabelSize(17), hDum->GetXaxis()->SetLabelFont(43);
  hDum->GetXaxis()->SetTitleSize(18), hDum->GetXaxis()->SetTitleFont(43);
  hDum->GetXaxis()->SetTitleOffset(3.0), hDum->GetXaxis()->CenterTitle();
  
  hDum->GetYaxis()->SetLabelSize(17), hDum->GetYaxis()->SetLabelFont(43);
  hDum->GetYaxis()->SetTitleSize(18), hDum->GetYaxis()->SetTitleFont(43);
  hDum->GetYaxis()->SetTitleOffset(3.2), hDum->GetYaxis()->CenterTitle();
  
  hDum->GetYaxis()->SetNdivisions(905,true);

  TH1F *hDum2 = (TH1F*) hDum->Clone("hDum2");
  hDum2->SetMinimum(-0.1), hDum2->SetMaximum(0.5);

  // Make canvas for the number of zvar
  for(int i=0; i<zvarBins.size()-1; i++){
    tcan.push_back(new TCanvas(Form("canvas_%d",i),"",1200,700)), makeMultiPanelCanvas(tcan[i],5,4,0.0,0.0,0.2,0.2,0.05);
    tcan2.push_back(new TCanvas(Form("canvas2_%d",i),"",1200,700)), makeMultiPanelCanvas(tcan2[i],5,4,0.0,0.0,0.2,0.2,0.05);
  }

  gStyle->SetOptFit(0000);

  // canvas 1 - efficiency
  cout<<"[Efficiency table] ===============================================\n"<<endl;
  for(int i=0; i<tgEffVV.size(); i++){
    tcan[i]->cd();
    for(int j=0; j<tgEffVV[i].size(); j++){
       double par[7];
       tcan[i]->cd(j+1);
       hDum->Draw("");
       tgEffVV[i][j]->Draw("pc");
       tgEffVV[i][j]->Fit(fit1,"R");
       drawText(Form("%1.1f ~ %1.1f",etaBins[j],etaBins[j+1]),0.24,0.24);
       
       fdata<<" 0 "<<" "<<fulltrk<<" "<<j<<" "<<i<<" "<<fit1->GetParameter(0)<<" "<<fit1->GetParameter(1)
	   <<" "<<fit1->GetParameter(2)<<" "<<fit1->GetParameter(3)<<" "<<fit1->GetParameter(4)
	   <<" "<<fit1->GetParameter(5)<<" "<<endl;

    }
  }
  cout<<"[Efficiency table] ===============================================\n"<<endl;

  // canvas 2 - fake rate
  cout<<"[Fake rate table] ===============================================\n"<<endl;
  for(int i=0; i<tgFakVV.size(); i++){
     tcan2[i]->cd();
     for(int j=0; j<tgFakVV[i].size(); j++){
	tcan2[i]->cd(j+1);
	hDum2->Draw("");
	tgFakVV[i][j]->Draw("pc"); 
	tgFakVV[i][j]->Fit("fitfake1","R");
	drawText(Form("%1.1f ~ %1.1f",etaBins[j],etaBins[j+1]),0.24,0.5);
	
	fdata<<" 1 "<<" "<<fulltrk<<" "<<j<<" "<<i<<" "<<fitfake1->GetParameter(0)<<" "<<fitfake1->GetParameter(1)
	    <<" "<<fitfake1->GetParameter(2)<<" "<<fitfake1->GetParameter(3)<<" "<<fitfake1->GetParameter(4)
	    <<" "<<fitfake1->GetParameter(5)<<" "<<endl;
     }
  }
  cout<<"[Fake rate table] ===============================================\n"<<endl;


  fdata.close();
  cout<<Form("%s_fulltrk%d.dat is created",infile.Data(),fulltrk)<<endl;
}


TH1D *hackedProjectionY(TH3F *me, const char *name, int ixmin, int ixmax, int izmin, int izmax){

  int pixmin = me->GetXaxis()->GetFirst();
  int pixmax = me->GetXaxis()->GetLast();
  int pizmin = me->GetZaxis()->GetFirst();
  int pizmax = me->GetZaxis()->GetLast();

  me->GetXaxis()->SetRange(ixmin,ixmax);
  me->GetZaxis()->SetRange(izmin,izmax);

  TH1D * h1 =  (TH1D*) me->Project3D("ye");
  h1->SetName( name );

  // restore axis range 
  me->GetXaxis()->SetRange(pixmin,pixmax);
  me->GetZaxis()->SetRange(pizmin,pizmax);

  return h1;

}


void drawText(const char *text, float xp, float yp){
   TLatex *tex = new TLatex(xp,yp,text);
   tex->SetTextFont(63);
   tex->SetTextSize(15);
   tex->SetTextColor(kBlack);
   tex->SetLineWidth(1);
   tex->SetNDC();
   tex->Draw();
}


void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
		      (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			(1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
  
  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
  
  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }

  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
			   Xlow[i],Ylow[j],Xup[i],Yup[j]);
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
      else pad[i][j]->SetLeftMargin(0);

      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);

      if(j==0) pad[i][j]->SetTopMargin(edge);
      else pad[i][j]->SetTopMargin(0);

      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}
