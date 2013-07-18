#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <utility>
#include <TROOT.h>
#include <TStyle.h>

#include "TFile.h"
#include "TCanvas.h"

#include "TH1F.h"
#include "TH1D.h"

#include "TH2F.h"
#include "TH2D.h"

#include "TF1.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TChain.h"
#include "TGraph.h"
#include "TSpline.h"

//#include "cms7TeV_lowpt.h"
//#include "cms900GeV_lowpt.h"

using namespace std;


//const double small = 1e-3;
//double ptb;
//std::vector<double> ptBins;


// In case any of these are called in this file =====================================
TCanvas *GetDummyCanvas();
TH1D* ratio_hist_to_func(TH1D* num, TF1* f3);
TH1D* ratio_hist_to_func(TH1D* num, TF1* f3, double minx, double maxx);
TH1D* ratio_hist_to_tspline(TH1D* num, TSpline3* f3, double minx, double maxx);
TH1D* ratio_hist_to_hist(TH1D* num, TH1D* den);
TGraphErrors* ratio_func_to_func(TGraphErrors* num, TF1* f3);
TGraphErrors* ratio_func_to_func(TGraphErrors* num, TF1* f3, double minx, double maxx);
TGraphErrors* ratio_func_to_tspline(TGraphErrors* num, TSpline3* f3, double minx, double maxx);
TH1D* HistIt(TGraphErrors* tg, TH1D* histC);
TGraphErrors* TgraphIt(TH1D* hist);
TGraphErrors* TgraphIt(TH1D* hist, bool xerr);
TH1D* RebinIt(TH1D* hist_pre, bool REBIN, char *name);
TH1D* RebinIt(TH1D* hist_pre, bool REBIN);
TH1D* RebinSimple(TH1D* hRS, char *name);
TH1D* RebinSimple(TH1D* hRS, char *name, bool skipCorr);
TH1D* RebinSimple(TH1D* hRS, char *name, bool varB, int rebFact);
TH1D* GetDummyHist(Float_t xmax, Float_t min, Float_t max,char *xttl,char *yttl);
TH1D* GetDummyHist(Float_t xmin, Float_t xmax, Float_t min, Float_t max,char *xttl,char *yttl, Bool_t stat);
void PrintXsection(TH1D* hist, float scale);
void PrintNameAndTitle(TH1D* hPNAT);
TH1D* FillFromFunc(TH1D* hist, TF1* func);
void CheckNumberOfBinEntries(TH1D* hCNOBE);
void PrepareGlobVar();
//double* rebinned();
// In case any of these are called in this file =====================================    
/*
const double small = 1e-3;
double ptb, ptbV;
double fact = 2.0;
double fact2 = 3.0;
std::vector<double> ptBins;
std::vector<double> ptBinsV;


void PrepareGlobVar(){

   double minW = 0.2; //minimum qunata 
   
   double nB = 0;
   double mF = 0;
   double iPt = 0;

   for(ptb =   0; ptb <    1-small; ptb +=  minW) ptBins.push_back(ptb);
   for(ptb =   1; ptb <    2-small; ptb +=  minW*2) ptBins.push_back(ptb);
   for(ptb =   2; ptb <    8-small; ptb +=  minW*6) ptBins.push_back(ptb);
   for(ptb =   8; ptb <   16-small; ptb +=  minW*10) ptBins.push_back(ptb);
   for(ptb =   16; ptb <  24-small; ptb +=  minW*20) ptBins.push_back(ptb);
   for(ptb =   24; ptb <  44-small; ptb +=  minW*50) ptBins.push_back(ptb);
   for(ptb =   44; ptb <  76-small; ptb +=  minW*80) ptBins.push_back(ptb);
   for(ptb =   76; ptb < 100-small; ptb +=  minW*120) ptBins.push_back(ptb);
   for(ptb =  100; ptb < 180-small; ptb +=  minW*200) ptBins.push_back(ptb);
   for(ptb =  180; ptb < 480-small; ptb +=  minW*500) ptBins.push_back(ptb);


   for(ptbV =   0; ptbV <   1-small; ptbV +=  0.05*fact ) ptBinsV.push_back(ptbV);
   for(ptbV =   1; ptbV <   2-small; ptbV +=  0.1*fact  ) ptBinsV.push_back(ptbV);
   for(ptbV =   2; ptbV <   5-small; ptbV +=  0.2*fact2 ) ptBinsV.push_back(ptbV);
   for(ptbV =   5; ptbV <  10-small; ptbV +=  0.5*fact  ) ptBinsV.push_back(ptbV);
   for(ptbV =  10; ptbV <  20-small; ptbV +=  1.0*fact  ) ptBinsV.push_back(ptbV);
   for(ptbV =  20; ptbV <  50-small; ptbV +=  2.0*fact2 ) ptBinsV.push_back(ptbV);
   for(ptbV =  50; ptbV < 100-small; ptbV +=  5.0*fact  ) ptBinsV.push_back(ptbV);
   for(ptbV = 100; ptbV < 200-small; ptbV += 10.0*fact  ) ptBinsV.push_back(ptbV);
   for(ptbV = 200; ptbV < 500-small; ptbV += 20.0*fact2 ) ptBinsV.push_back(ptbV);

}
*/

/*
TH1D* RebinSimple(TH1D* hRS, char *name, bool skipCorr){

   double Odbin = hRS->GetBinWidth(1); // assume const bin width
   cout<<"[RebinSimple]----------------------"<<endl;
   cout<<"original bin width = "<<Odbin<<endl;
   cout<<"number of rebins = "<<ptBins.size()-1<<endl;
   
   TH1D* tempH = (TH1D*) hRS->Rebin(ptBinsV.size()-1,name,&ptBinsV[0]);
   
   // correction due to rebinning  
   if(skipCorr){
      for(int i=0;i<ptBins.size()-1;i++){
	 double dbin = tempH->GetBinWidth(i+1);
	 double ratio = dbin/Odbin;
	 tempH->SetBinContent(i+1,tempH->GetBinContent(i+1)/ratio);
	 tempH->SetBinError(i+1,tempH->GetBinError(i+1)/ratio);
      }
   }
   cout<<"pointer returned = "<<tempH<<endl;
   cout<<"[RebinSimple]----------------------"<<endl;
   return tempH;
}

TH1D* RebinSimple(TH1D* hRS, char *name){

   //TH1D *hRS = (TH1D*) hRS_pre->Clone("hRS"); 
   double Odbin = hRS->GetBinWidth(1); // assume const bin width

   cout<<"[RebinSimple]----------------------"<<endl;
   cout<<"original bin width = "<<Odbin<<endl;
   cout<<"number of rebins = "<<ptBins.size()-1<<endl;

   TH1D* tempH = (TH1D*) hRS->Rebin(ptBins.size()-1,name,&ptBins[0]);

   // correction due to rebinning
   for(int i=0;i<ptBins.size()-1;i++){
      double dbin = tempH->GetBinWidth(i+1);
      double ratio = dbin/Odbin;

      tempH->SetBinContent(i+1,tempH->GetBinContent(i+1)/ratio);
      tempH->SetBinError(i+1,tempH->GetBinError(i+1)/ratio);
   }
   cout<<"pointer returned = "<<tempH<<endl;
   cout<<"[RebinSimple]----------------------"<<endl;
   return tempH;
}

TH1D* RebinSimple(TH1D* hRS, char *name, bool varB, int rebFact){

   double Odbin = hRS->GetBinWidth(1); // assume const bin width

   cout<<"[RebinSimple]----------------------"<<endl;
   cout<<"original bin width = "<<Odbin<<endl;
   cout<<"number of rebins = "<<ptBins.size()-1<<endl;

   if(!varB) TH1D* tempH = (TH1D*) hRS->Rebin(ptBins.size()-1,name,&ptBins[0]);
   TH1D* tempH = (TH1D*) hRS->Rebin(rebFact,name);

   // correction due to rebinning                                                                                                                           
   for(int i=0;i<ptBins.size()-1;i++){
      double dbin = tempH->GetBinWidth(i+1);
      double ratio = 1;
      if(!varB) ratio = dbin/Odbin;
      else ratio = rebFact;
      tempH->SetBinContent(i+1,tempH->GetBinContent(i+1)/ratio);
      tempH->SetBinError(i+1,tempH->GetBinError(i+1)/ratio);
   }
   cout<<"pointer returned = "<<tempH<<endl;
   cout<<"[RebinSimple]----------------------"<<endl;
   return tempH;
}
*/

void CheckNumberOfBinEntries(TH1D* hCNOBE){
   int nbin = hCNOBE->GetNbinsX();
   for(int i=0;i<nbin;i++){
      cout<<"x value : "<<hCNOBE->GetBinCenter(i+1)<<" entries (from (err*err)) : "<<hCNOBE->GetBinError(i+1)*hCNOBE->GetBinError(i+1)<<endl;
   }
}

void PrintNameAndTitle(TH1D* hPNAT){
   cout<<"name = "<<hPNAT->GetName()<<" && title = "<<hPNAT->GetTitle()<<endl;
}

TH1D* FillFromFunc(TH1D* hFFF_pre, TF1* fFFF){

   TH1D* hFFF = (TH1D*) hFFF_pre->Clone("hFFF");
   hFFF->Reset();

   int nbin = hFFF->GetNbinsX();
   //hist->Sumw2(); 
   for(int i = 0; i<nbin ;i++){
      double low = hFFF->GetXaxis()->GetBinLowEdge(i+1);
      double high = hFFF->GetXaxis()->GetBinUpEdge(i+1);
      double deltaN = fFFF->Integral(low,high);
      double pt = hFFF->GetBinCenter(i+1);
      double dbin = hFFF->GetBinWidth(i+1);
      hFFF->Fill(pt,deltaN/dbin);
      hFFF->SetBinError(i+1,0);
   }
   return hFFF;
}

TH1D* ratio_hist_to_func(TH1D* num, TF1* f3){

   cout<<"[Ratio to fit used]"<<endl;

   TH1D *hRatio = (TH1D*) num->Clone("hRatio");
   //hRatio->Reset();

   int nbin = hRatio->GetNbinsX();
   //cout<<"[test1]"<<endl;
   for(int i=0;i<nbin;i++){

      double cms_value = (double) f3->Eval(hRatio->GetBinCenter(i+1));
      //double ratio = num->GetBinContent(i+1)/cms_value; 
      //double ratio_err =  (num->GetBinError(i+1)/num->GetBinContent(i+1))*ratio;
      double ratio = hRatio->GetBinContent(i+1)/cms_value;                                                                                                   
      double ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;

      hRatio->SetBinContent(i+1,ratio);
      hRatio->SetBinError(i+1,ratio_err);
   }
   //cout<<"[test2]"<<endl;
   return hRatio;
}



TH1D* ratio_hist_to_func(TH1D* num, TF1* f3, double minx, double maxx){

   cout<<"[Ratio to fit with range used]"<<endl;

   TH1D *hRatio = (TH1D*) num->Clone("hRatio");

   int nbin = num->GetNbinsX();

   for(int i=0;i<nbin;i++){

      double cms_value = (double) f3->Eval(hRatio->GetBinCenter(i+1));

      double ratio = -999;
      double ratio_err = 0;

      if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
         ratio = hRatio->GetBinContent(i+1)/cms_value;
	 ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
      }else{
         ratio = -999;
         ratio_err = 0.0;
      }

      hRatio->SetBinContent(i+1,ratio);
      hRatio->SetBinError(i+1,ratio_err);
   }

   return hRatio;

}


TH1D* ratio_hist_to_tspline(TH1D* num, TSpline3* ts, double minx, double maxx){

  cout<<"[Ratio to fit with range used]"<<endl;

  TH1D *hRatio = (TH1D*) num->Clone("hRatio");

  int nbin = num->GetNbinsX();

  for(int i=0;i<nbin;i++){

    double cms_value = (double) ts->Eval(hRatio->GetBinCenter(i+1));

    double ratio = -999;
    double ratio_err = 0;

    if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
      ratio = hRatio->GetBinContent(i+1)/cms_value;
      ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
    }else{
      ratio = -999;
      ratio_err = 0.0;
    }

    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,ratio_err);
  }

  return hRatio;

}

TH1D* ratio_hist_to_hist(TH1D* num, TH1D* den){

  //-----------------------------------------------------------
  // this is used when binning is identical but a bin 
  // range is different compare the ranges and use narrower one 
  //-----------------------------------------------------------
  
  int nNum = num->GetXaxis()->GetNbins();
  int nDen = den->GetXaxis()->GetNbins();
  int nbin = (nNum>nDen) ? nDen : nNum;

  TH1D* ratio = (nNum>nDen) ? (TH1D*) den->Clone("ratio") : (TH1D*) num->Clone("ratio");   
  ratio->Reset();     
  
  for(int i=0;i<nbin;i++){
    if(den->GetBinContent(i+1)==0) return 0;
    double rat = num->GetBinContent(i+1)/den->GetBinContent(i+1);
    ratio->SetBinContent(i+1,rat);
    ratio->SetBinError(i+1,0.0001*rat); // need to be calculated correctly
  }

  return ratio;

}


TGraphErrors* ratio_func_to_func(TGraphErrors* num, TF1* f3){

   cout<<"[Ratio to fit used]"<<endl;

   TGraphErrors *tg;

   int nbin = num->GetN();

   const int nlines = nbin;
   double pt[nlines], pterr[nlines];
   double xsec[nlines], xsecerr[nlines];
   double ratio[nlines], ratioerr[nlines];

   for(int i=0;i<nbin;i++){

      num->GetPoint(i,pt[i],xsec[i]);
      xsecerr[i] = num->GetErrorY(i);
      double cms_value = (double) f3->Eval(pt[i]);

      ratio[i] = xsec[i]/cms_value;
      //ratioerr[i] = xsecerr[i]/cms_value;
      ratio[i] = (xsecerr[i]/xsec[i])*ratio[i];
      pterr[i] = 0;
   }

   tg = new TGraphErrors(nlines,pt,ratio,pterr,ratioerr);
   return tg;
}

TGraphErrors* ratio_func_to_func(TGraphErrors* num, TF1* f3, double minx, double maxx){

   cout<<"[Ratio to fit used]"<<endl;

   TGraphErrors *tg;

   int nbin = num->GetN();

   const int nlines = nbin;
   double pt[nlines], pterr[nlines];
   double xsec[nlines], xsecerr[nlines];
   double ratio[nlines], ratioerr[nlines];

   for(int i=0;i<nbin;i++){

      num->GetPoint(i,pt[i],xsec[i]);
      xsecerr[i] = num->GetErrorY(i);
      double cms_value = (double) f3->Eval(pt[i]);

      if(pt[i]>minx && pt[i]<maxx){
         ratio[i] = xsec[i]/cms_value;
	 ratioerr[i] = (xsecerr[i]/xsec[i])*ratio[i];
         //ratioerr[i] = xsecerr[i]/cms_value;
      }else{
         ratio[i] = -999;
         ratioerr[i] = 0.0;
      }

      pterr[i] = 0;
   }

   tg = new TGraphErrors(nlines,pt,ratio,pterr,ratioerr);
   return tg;
}


TGraphErrors* ratio_func_to_tspline(TGraphErrors* num, TSpline3* ts, double minx, double maxx){

  cout<<"[Ratio to fit used]"<<endl;

  TGraphErrors *tg;

  int nbin = num->GetN();

  const int nlines = nbin;
  double pt[nlines], pterr[nlines];
  double xsec[nlines], xsecerr[nlines];
  double ratio[nlines], ratioerr[nlines];

  for(int i=0;i<nbin;i++){

    num->GetPoint(i,pt[i],xsec[i]);
    xsecerr[i] = num->GetErrorY(i);
    double cms_value = (double) ts->Eval(pt[i]);

    if(pt[i]>minx && pt[i]<maxx){
      ratio[i] = xsec[i]/cms_value;
      ratioerr[i] = (xsecerr[i]/xsec[i])*ratio[i];
      //ratioerr[i] = xsecerr[i]/cms_value;                                                                                                              
    }else{
      ratio[i] = -999;
      ratioerr[i] = 0.0;
    }

    pterr[i] = 0;
  }

  tg = new TGraphErrors(nlines,pt,ratio,pterr,ratioerr);
  return tg;
}



// convert tgraph  to hist
// assume predefined hist
TH1D* HistIt(TGraphErrors* tg, TH1D* histC){

   //TH1D* histC = (TH1D*) hist->Clone("histC");
   int nbins = tg->GetN();

   const int nlines = nbins;
   
   /*
   double *x[nlines], *y[nlines];
   
   x = (double*) tg->GetX();
   y = (double*) tg->GetY();
   
   double dbin = (x[1]-x[0]);
   double nbinx = nbins;
   
   hist = new TH1D("hist","",nbins,x[0]-(dbin/2),x[nbinx]+(dbin/2));
   hist->Sumw2();
   */

   double pt[nlines], xsec[nlines];
   double pterr[nlines], xsecerr[nlines];

   for(int i = 0; i<nbins; i++ ){
      tg->GetPoint(i,pt[i],xsec[i]);
      pterr[i] = tg->GetErrorX(i);
      xsecerr[i] = tg->GetErrorY(i);
      histC->SetBinContent(i+1,xsec[i]);
      histC->SetBinError(i+1,xsecerr[i]);
   }
   return histC;
}

/*
TH1D* HistIt(TGraphErrors* tg){

  int nbins = tg->GetN();
  const int nlines = nbins;

  double pt[nlines], xsec[nlines];
  double pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    tg->GetPoint(i,pt[i],xsec[i]);
    pterr[i] = tg->GetErrorX(i);
    xsecerr[i] = tg->GetErrorY(i);
  }

  double dbin = pt[1]-pt[0];
  if(dbin!=(pt[2]-pt[1])||dbin==0){
    Error("[HistIt]","variable bins or bin width - 0");
    return;
  }
  double xi = pt[0]-dbin;
  double xf = pt[nbins-1]+dbin;
  
  TH1D *hist = new TH1D("","",nbins,xi,xf);
 
}
*/

// convert hist to tgraph
TGraphErrors* TgraphIt(TH1D* hist){

   TGraphErrors *tg;
   int nbins = hist->GetNbinsX();

   const int nlines = nbins;

   float pt[nlines], xsec[nlines];
   float pterr[nlines], xsecerr[nlines];

   for(int i = 0; i<nbins; i++ ){
      pt[i] = hist->GetBinCenter(i+1);
      xsec[i] = hist->GetBinContent(i+1);
      xsecerr[i] = hist->GetBinError(i+1);
      pterr[i] = 0;
   }

   tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
   return tg;
}

TGraphErrors* TgraphIt(TH1D* hist, bool xerr){

  TGraphErrors *tg;
  int nbins = hist->GetNbinsX();

  const int nlines = nbins;

  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    if(xerr) pterr[i] = 0.5*(hist->GetBinWidth(i+1));
    else pterr[i] = 0;
    //cout<<" pt = "<<pt[i]<<" and bin width = "<<hist->GetBinWidth(i+1)<<endl;
  }

  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}





TH1D* RebinIt(TH1D* hist_pre, bool REBIN, char *name){

   // ########################################### prepare variable bin ######

   TH1D* hRebinIt = (TH1D*) hist_pre->Clone("hRebinIt");
   //hist->Reset();

   Double_t dBin = hRebinIt->GetBinWidth(1);//assume hist has constan bin width 
   
   Int_t NumBins = hRebinIt->GetNbinsX();
   const Int_t Nlines = NumBins;
   
   Int_t binTemp = 0;
   Int_t nbins = 0;
   Double_t binsTemp[Nlines+1];
   Int_t totBins = Nlines;
   
   for(Int_t i = 0; i < totBins ; i++){
      binsTemp[i] = binTemp;
      binTemp++;
   }
   binsTemp[totBins] = totBins;
   
   int bin = 0;
   Double_t bins[Nlines+1];
   Double_t binWidth = hRebinIt->GetBinWidth(1);
   cout<<"bin width of original histogram "<<binWidth<<endl;
   
   if(REBIN){
      while (bin < totBins) {
	 bins[nbins] = binWidth*binsTemp[bin];
	 //cout<<"bins[nbins] = "<<bins[nbins]<<endl;
	 nbins++;
	 // DATA

	 /*
	 if (bin < 3)          bin += 1;
	 else if(bin < 8)      bin += 2;
	 else if (bin < 12)    bin += 4;
	 else if (bin < 30)    bin += 6;
	 else if (bin < 40)    bin += 8;
	 else if (bin < 70)    bin += 25;
	 else if (bin < 100)   bin += 40;
	 else if (bin < 200)   bin += 55;
	 else                  bin += 100;
	 */

	 if (bin < 3)          bin += 1;
         else if(bin < 8)      bin += 2;
         else if (bin < 12)    bin += 4;
         else if (bin < 30)    bin += 6;
         else if (bin < 40)    bin += 8;
         else if (bin < 70)    bin += 25;
         else if (bin < 100)   bin += 45;
         else if (bin < 200)   bin += 65;
         else                  bin += 150;


	 //bin += 2;
      }
      bins[nbins] = binWidth*binsTemp[totBins];
      //cout<<"bins[nbins] = "<<bins[nbins]<<endl;
   }
   
   cout<<"number of bins after rebinned = "<<nbins<<endl;
   
   //################################################# Rebin ###################            
   Char_t hname[500];
   //sprintf(hname,"pre_hRInvX");
   TH1D* pre_hRInvX = (TH1D*) hRebinIt->Rebin(nbins,name,bins);

   const Int_t RNlines = nbins;
   
   Float_t ptR[RNlines], xsecR[RNlines];
   Float_t err_ptR[RNlines], err_xsecR[RNlines];
   
   
   for(Int_t j = 0; j<nbins; j++ ){
      
      // In order to scale the content properly (due to rebinning)
      
      Float_t dbin = pre_hRInvX->GetBinWidth(j+1);
      Float_t ratio = dbin/dBin;
      Float_t content = pre_hRInvX->GetBinContent(j+1);
      Float_t err = pre_hRInvX->GetBinError(j+1);

      pre_hRInvX->SetBinContent(j+1,(content/ratio));
      pre_hRInvX->SetBinError(j+1,(err/ratio));
      
      ptR[j] = pre_hRInvX->GetBinCenter(j+1);
      xsecR[j] = (content/ratio);
      err_ptR[j] = 0;
      err_xsecR[j] = (err/ratio);
      
   }

   delete hist_pre;
   delete hRebinIt;

   return pre_hRInvX;
   
}


TH1D* RebinIt(TH1D* hist_pre, bool REBIN){

   // ########################################### prepare variable bin ######

   TH1D* hRebinIt2 = (TH1D*) hist_pre->Clone("hist");
   //hist->Reset();

   Double_t dBin = hRebinIt2->GetBinWidth(1);//assume hist has constan bin width 
   
   Int_t NumBins = hRebinIt2->GetNbinsX();
   const Int_t Nlines = NumBins;
   
   Int_t binTemp = 0;
   Int_t nbins = 0;
   Double_t binsTemp[Nlines+1];
   Int_t totBins = Nlines;
   
   for(Int_t i = 0; i < totBins ; i++){
      binsTemp[i] = binTemp;
      binTemp++;
   }
   binsTemp[totBins] = totBins;
   
   Int_t bin = 0;
   Double_t bins[Nlines+1];
   Double_t binWidth = hRebinIt2->GetBinWidth(1);
   cout<<"bin width of original histogram "<<binWidth<<endl;
   
   if(REBIN){
      while (bin < totBins) {
	 bins[nbins] = binWidth*binsTemp[bin];
	 //cout<<"bins[nbins] = "<<bins[nbins]<<endl;	
	 nbins++;
	 // DATA

	 /*
	 if (bin < 3)          bin += 1;
	 else if(bin < 8)      bin += 2;
	 else if (bin < 12)    bin += 4;
	 else if (bin < 30)    bin += 6;
	 else if (bin < 40)    bin += 8;
	 else if (bin < 70)    bin += 25;
	 else if (bin < 100)   bin += 40;
	 else if (bin < 200)   bin += 55;
	 else                  bin += 100;
	 */

         if (bin < 3)          bin += 1;
         else if(bin < 8)      bin += 2;
         else if (bin < 12)    bin += 4;
         else if (bin < 30)    bin += 6;
         else if (bin < 40)    bin += 8;
         else if (bin < 70)    bin += 25;
         else if (bin < 100)   bin += 45;
         else if (bin < 200)   bin += 65;
         else                  bin += 150;

	 //bin += 2;
      }
      bins[nbins] = binWidth*binsTemp[totBins];
      //cout<<"bins[nbins] = "<<bins[nbins]<<endl;	
   }
   
   cout<<"number of bins after rebinned = "<<nbins<<endl;
   
   //################################################# Rebin ###################            
   Char_t hname[500];
   sprintf(hname,"pre_hRInvX");
   TH1D* pre_hRInvX = (TH1D*) hRebinIt2->Rebin(nbins,hname,bins);

   const Int_t RNlines = nbins;
   
   Float_t ptR[RNlines], xsecR[RNlines];
   Float_t err_ptR[RNlines], err_xsecR[RNlines];
   
   
   for(Int_t j = 0; j<nbins; j++ ){
      
      // In order to scale the content properly (due to rebinning)
      
      Float_t dbin = pre_hRInvX->GetBinWidth(j+1);
      Float_t ratio = dbin/dBin;
      Float_t content = pre_hRInvX->GetBinContent(j+1);
      Float_t err = pre_hRInvX->GetBinError(j+1);

      pre_hRInvX->SetBinContent(j+1,(content/ratio));
      pre_hRInvX->SetBinError(j+1,(err/ratio));
      
      ptR[j] = pre_hRInvX->GetBinCenter(j+1);
      xsecR[j] = (content/ratio);
      err_ptR[j] = 0;
      err_xsecR[j] = (err/ratio);
      
   }

   delete hist_pre;
   delete hRebinIt2;

   return pre_hRInvX;
   
}

TCanvas *GetDummyCanvas(){
   TCanvas *dum = new TCanvas("call","call",550,410);
   return dum;
}

TH1D* GetDummyHist(Float_t xmin, Float_t xmax, Float_t min, Float_t max,char *xttl,char *yttl, Bool_t stat) {

   //gStyle->SetOptFit(1111);
  //gStyle->SetOptStat(1111);

   TH1D *dum;
   dum = new TH1D("dum","",100,xmin,xmax);

   dum->SetMinimum(min);
   dum->SetMaximum(max);
   dum->SetStats(stat);

   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();

   return dum;

}



TH1D* GetDummyHist(Float_t xmax, Float_t min, Float_t max,Char_t *xttl,Char_t *yttl) {

   TH1D *dum;
   //dum = new TH1D("dum","",100,0.0,xmax);
   //dum = new TH1D("dum","",100,0.5,xmax);  
   dum = new TH1D("dum","",100,0.4,xmax);
   //dum = new TH1D("dum","",100,0.2,xmax);                                                                                                                                                                   

   dum->SetMinimum(min);
   dum->SetMaximum(max);
   dum->SetStats(0);

   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();

   return dum;

}


void PrintXsection(TH1D* hist, float scale){

   cout<<"\n"<<endl;
   cout<<"[print the cross section of "<<hist->GetName()<<"]================="<<endl;
   cout<<"scale factor used : "<<scale<<endl;
   for(int i = 0; i<hist->GetNbinsX(); i++){
      cout<<"      "<<hist->GetBinCenter(i+1)<<"      "
          <<scale*hist->GetBinContent(i+1)<<"       "
          <<scale*hist->GetBinError(i+1)<<endl;
   }
   cout<<"[print the cross section of "<<hist->GetName()<<"]================="<<endl;
   cout<<"\n"<<endl;
}

#endif
