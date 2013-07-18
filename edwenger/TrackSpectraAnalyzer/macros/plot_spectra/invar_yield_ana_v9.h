//------------------------------------------------------------------- 
// June 03 2010, Andre S. Yoon 
// - header file for invar_yield_ana_v9.C
// - efficiency correction is done for tracking
// - efficiency correction is done for event selection 
// - return various histograms/tgraph in struct
//--------------------------------------------------------------------



#include <utility>
#include <iostream>

#include "/home/sungho/plots/common/utilities.h"

#include <TROOT.h>
#include <TStyle.h>

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

#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TChain.h"
#include "TGraph.h"

#include "TKey.h";

using namespace std;

const int NHISTR = 5;

TH1D* FillFromFunc(TH1D* hist, TF1* func);
int GetFirstNonEmptyBin(TH1D* hist);
int GetLastNonEmptyBin(TH1D* hist);
double GetEffFactor(TH3F* num, TH3F* den, double pt, double eta, double jet);
double GetFakFactor(TH3F* num, TH3F* den, double pt, double eta, double jet);
double CorrectNevt(TH1F* histEVT);
double CorrectNevtWithMult(TH1F* histEVT);
double GetZeroBinFraction(TH1F* hist);
double GetOneBinFraction(TH1F* hist);
double GetTwoBinFraction(TH1F* hist);
double GetThreeBinFraction(TH1F* hist);
double GetZOTTBinFraction(TH1F* hist);
double GetOTTBinFraction(TH1F* hist);
double CorrectNevtV2(TH1F* histDen, TH1F* histNum);
std::pair<TH3F*,TH3F*> getEffHists(char *file, char *dirC, char *histN, char *histD);
void PrintDetailsOfHist(TH3F* hSpectra);
std::pair<TH1D*,TH1D*> CorrectForMomRes(char** file, TH1D* hist_pre, bool rebOnly, double ieta, double feta, double xmin, double xmax);


std::pair<TH1D*,TH1D*> CorrectForMomRes(char* file, TH1D* hist_pre, bool rebOnly, double ieta, double feta, double xmin, double xmax){ 
//std::pair<TH1D*,TH1D*> CorrectForMomRes(char** afile, TH1D* hist_pre, bool rebOnly, double ieta, double feta, double xmin, double xmax){ 


   cout<<"\n"<<endl;
   cout<<"[Mom resolution and binning correction ]========================================"<<endl;
   
   double pts_1st, pts_2nd, pts_3rd, pts_4th;
   pts_1st = 10, pts_2nd = 20, pts_3rd = 30, pts_4th = 70;

   
   TH1D *hCFMR = (TH1D*) hist_pre->Clone("hCFMR");

   //==============  obtain fit 
   cout<<"fitting from the spectra"<<endl;
   TF1* fitF = new TF1("fitpythia","[0]*(1+(x/[1])+(pow(x,2)/[2])+(pow(x,3)/[3])+(pow(x,4)/[4])+(pow(x,5)/[5]))^[6]",xmin,xmax);
   fitF->SetParameters(2.96361e+08,2.66339e-01,1.34086e-01,2.96428e-01,1.74015e+00,5.16931e+00,-1.70643e+00);
   fitF->SetLineColor(2);
   fitF->SetLineWidth(2);

   hCFMR->Fit(fitF,"NRO");
   TH1D* hFit = (TH1D*) FillFromFunc(hCFMR,fitF);


   //==============  obtain mom res histograms
   TH3F *hRes3D = 0;
   TH2D *hRes2D = 0;

   if(!rebOnly){
      cout<<"file to be processed = "<<file<<endl;
      TFile *fEVT = new TFile(file);
      hRes3D = (TH3F*) fEVT->Get("trkEffAnalyzer/hresStoR3D");
      
      int binmaxeta = hRes3D->GetXaxis()->FindBin(feta);
      int binmineta = hRes3D->GetXaxis()->FindBin(-1.0*feta);
      
      hRes3D->GetXaxis()->SetRange(binmineta,binmaxeta);
      hRes2D = (TH2D*) hRes3D->Project3D("zy"); 
   }

   //==============  convolute with res!
   TH1D *hRes1D_temp = 0;

   if(!rebOnly){

      int nbins = hCFMR->GetNbinsX();
      cout<<"nbins = "<<nbins<<endl;
      
      int index = 0;

      for(int i=0;i<nbins;i++){
	 
	 bool nullArea = false; // this is needed for emtpy mom res bin

	 double pts = hCFMR->GetBinCenter(i);

	 hRes1D_temp = (TH1D*) hRes2D->ProjectionY("",i+1,i+1,"e");

	 // check if works..
	 int firstB = GetFirstNonEmptyBin(hRes1D_temp);
	 int lastB = GetLastNonEmptyBin(hRes1D_temp);
	 
	 double area = hRes1D_temp->Integral();
	 
	 if(area==0) nullArea =true;
	 else hRes1D_temp->Scale(1./area);
	 
	 double corr;

	 for(int j=firstB;j<lastB+1;j++){
	    
	    if(!nullArea) corr = hRes1D_temp->GetBinContent(j);
	    else corr = 1.0; // no information on res-> no correction
	    
	    double ptr = hRes1D_temp->GetBinCenter(j);
	    
	    int rpt_bin = hRes2D->GetYaxis()->FindBin(ptr); // should be same as j!
	    int spt_bin = i+1;
	    //double smeared = corr*(fitF->Eval(ptr)); // weighted by the fit 
	    double smeared = corr*(hFit->GetBinContent(rpt_bin));// weighted by filled histogram (from fit)  
	    

	    hRes2D->SetBinContent(spt_bin,rpt_bin,smeared);
	    hRes2D->SetBinError(spt_bin,rpt_bin,0);
	 }
	 delete hRes1D_temp;
      }
   }

   TH1D* histCorr = 0;
   TH1D* histC = 0;
   TH1D* histCR = 0;

   // correction sould be like corr1 (res) * corr2 (rebin) 
   // instead of separate!
   if(!rebOnly){
      histCorr = (TH1D*) hRes2D->ProjectionY("",1,hRes2D->GetXaxis()->GetLast(),"e"); // convoluted and projected!
      histC = (TH1D*) histCorr->Clone("histC");
      histCR = (TH1D*) RebinIt(histCorr,true,"histCR");
   }else{
      histCorr = (TH1D*) FillFromFunc(hCFMR,fitF);
      histC = (TH1D*) histCorr->Clone("histC"); 
      histCR = (TH1D*) RebinIt(histCorr,true,"histCR");
   }

   TH1D* hSPEC = (TH1D*) hCFMR->Clone("hSPEC");
   TH1D* hSPECR = (TH1D*) RebinIt(hCFMR,true,"hSPECR");

   //cout<<"name = "<<hist->GetName()<<" && title = "<<hist->GetTitle()<<endl;  
   //PrintNameAndTitle(hist);
   TH1D* histC_ratio = (TH1D*) ratio_hist_to_func(histC,fitF);
   TH1D* histCR_ratio = (TH1D*) ratio_hist_to_func(histCR,fitF);

   for(int i=0;i<histC_ratio->GetNbinsX();i++){
      double yield    = hSPEC->GetBinContent(i+1);
      double yielderr = hSPEC->GetBinError(i+1);
      double corrR = histC_ratio->GetBinContent(i+1);
      if(corrR==0) corrR = 1; // corrR=0 when yield =0
      hSPEC->SetBinContent(i+1,yield/corrR);
      hSPEC->SetBinError(i+1,yielderr/corrR);
   }

   for(int i=0;i<histCR_ratio->GetNbinsX();i++){
      double yield    = hSPECR->GetBinContent(i+1);
      double yielderr = hSPECR->GetBinError(i+1);
      double corrR = histCR_ratio->GetBinContent(i+1);
      if(corrR==0) corrR = 1;
      hSPECR->SetBinContent(i+1,yield/corrR);
      hSPECR->SetBinError(i+1,yielderr/corrR);
   }

   
   //To see how well spectra is fit and the correction factor
   //TCanvas *cIn = new TCanvas("cIn","cIn",500,450);
   //cIn->cd();
   //cIn->SetLogy();
   
   //TH1D *dummyH = new TH1D("dummyH","",100,0.0,200);
   //dummyH->Draw("");
   //hCFMR->SetMarkerSize(0.8);

   //cout<<"name = "<<hist->GetName()<<" && title = "<<hist->GetTitle()<<endl;
   //hCFMR->Draw("histsame");
   //histC_ratio->Draw("histsame");
   //fitF->Draw("pzsame");

   /*
   TCanvas *cIn2 = new TCanvas("cIn2","cIn2",500,450);
   cIn2->cd();
   histC_ratio->Draw("hist");
   histCR_ratio->SetLineColor(kRed);
   histCR_ratio->Draw("histsame");
   */   

   cout<<"[Mom resolution and binning correction ]========================================"<<endl;
   cout<<"\n"<<endl;

   return  std::pair<TH1D*,TH1D*>(hSPEC,hSPECR);

}

/*
TH1D* FillFromFunc(TH1D* hist, TF1* func){
   int nbin = hist->GetNbinsX();
   //hist->Sumw2();
   for(int i = 0; i<nbin ;i++){
      double low = hist->GetXaxis()->GetBinLowEdge(i+1);
      double high = hist->GetXaxis()->GetBinUpEdge(i+1);
      double deltaN = func->Integral(low,high);
      double pt = hist->GetBinCenter(i+1);
      double dbin = hist->GetBinWidth(i+1);
      hist->Fill(pt,deltaN/dbin);
      hist->SetBinError(i+1,0);
   }

   return hist;
}
*/

int GetFirstNonEmptyBin(TH1D* hist){
   int bin = 0;
   for(int i=0;i<hist->GetNbinsX();i++){
      bin = i+1;
      //cout<<hist->GetBinContent(bin)<<endl;                                                                                                                           
      if(hist->GetBinContent(bin)>0) break;
   }
   //cout<<"bin = "<<bin<<endl;                                                                                                                                         
   return bin;
}

int GetLastNonEmptyBin(TH1D* hist){
   int bin = 0;
   for(int i=0;i<hist->GetNbinsX();i++){
      bin = (hist->GetNbinsX())-i;
      if(hist->GetBinContent(bin)>0) break;
   }
   return bin;
}


double GetEffFactor(TH3F* num, TH3F* den, double pt, double eta, double jet){

   int num_bin = num->FindBin(eta,pt,jet);
   int den_bin = den->FindBin(eta,pt,jet);

   //double num_bin_width = num->GetYaxis()->GetBinWidth(num_bin);
   //double den_bin_width = den->GetYaxis()->GetBinWidth(den_bin);

   double n_num = num->GetBinContent(num_bin);
   double n_den = den->GetBinContent(den_bin);
   if(n_den == 0  || n_num == 0 ) return 1;
   else return n_num/n_den;
}


double GetFakFactor(TH3F* num, TH3F* den, double pt, double eta, double jet){

   int num_bin = num->FindBin(eta,pt,jet);
   int den_bin = den->FindBin(eta,pt,jet);

   //double num_bin_width = num->GetYaxis()->GetBinWidth(num_bin);
   //double den_bin_width = den->GetYaxis()->GetBinWidth(den_bin);

   double n_num = num->GetBinContent(num_bin);
   double n_den = den->GetBinContent(den_bin);
   if(n_den == 0) return 1; // meaing nothing is reconstructed.                                                                                                           
   if(n_num == 0) return 0;
   else return n_num/n_den; // be careful with def, with MTV, its 1-(n_num/n_den)                                                                                         
}


double CorrectNevt(TH1F* histEVT){
   double nevt = 0;
   for(int i=0;i<histEVT->GetNbinsX();i++){
      nevt = nevt + (histEVT->GetBinContent(i+1)/histEVT->GetBinCenter(i+1));
   }
   return nevt;
}

double CorrectNevtWithMult(TH1F* histEVT){
   return (double) histEVT->Integral();
}



// fraction is defined as a fraction of event with mult =x out of N_{tot}
// not out of N_{i} (where i = 4,5,6.. (summed implicitly))

double GetZeroBinFraction(TH1F* hist){

   double total = hist->Integral();
   //double non_zerobin = hist->Integral(2,hist->GetNbinsX());
   //double fraction = (total-non_zerobin)/total;
   double  zerobin = hist->Integral(1,1);
   double fraction = (zerobin)/total;
   return fraction;
}

double GetOneBinFraction(TH1F* hist){

   double total = hist->Integral();
   double zerobin = hist->Integral(2,2);
   double fraction = (zerobin)/total;
   return fraction;
}

double GetTwoBinFraction(TH1F* hist){

   double total = hist->Integral();
   double zerobin = hist->Integral(3,3);
   double fraction = (zerobin)/total;
   return fraction;
}

double GetThreeBinFraction(TH1F* hist){

   double total = hist->Integral();
   double zerobin = hist->Integral(4,4);
   double fraction = (zerobin)/total;
   return fraction;
}


double GetZOTTBinFraction(TH1F* hist){

   double total = hist->Integral();
   double zerobin = hist->Integral(1,4);
   double fraction = (zerobin)/total;
   return fraction;
}

double GetOTTBinFraction(TH1F* hist){

   double total = hist->Integral(2,hist->GetNbinsX());//except zero bin 
   //double non_zerobin = hist->Integral(5,hist->GetNbinsX());
   //double fraction = (total-non_zerobin)/total;

   double zerobin = hist->Integral(1,4);
   double fraction = (zerobin)/total;

   return fraction;
}

double CorrectNevtV2(TH1F* histDen, TH1F* histNum){

   double nevt = 0;
   TH1F* histNum_dum = (TH1F*) histNum->Clone("histNum_dum");
   TGraphAsymmErrors *gNSDvtx = new TGraphAsymmErrors();

   gNSDvtx->BayesDivide(histNum,histDen);

   int n_bins = gNSDvtx->GetN();
   const int NBINS = n_bins;
   double x[NBINS], y[NBINS];

   for(int i=0;i<n_bins;i++){
      gNSDvtx->GetPoint(i,x[i],y[i]);
   }

   for(int i=0;i<histNum->GetNbinsX();i++){
      histNum_dum->SetBinContent(i+1,(1/y[i])*histNum->GetBinContent(i+1));
   }

   return histNum_dum->Integral();
}


std::pair<TH3F*,TH3F*> getEffHists(char *file, char *dirC, char *histN, char *histD){

   TFile *efile = new TFile(file,"read");
   TDirectoryFile *efileDIR = (TDirectoryFile*) efile->GetDirectory(dirC);

   TIter next(efileDIR->GetListOfKeys());
   TKey *key;

   Char_t name[100];

   TH3F *hpneff3d=0;
   TH3F *hpdeff3d=0;

   while ((key=(TKey*)next())) {
      sprintf(name,"%s",key->GetName());
      if(strcmp((key->GetName()),(histN))==0){
         //cout<<"[getEffHists] Your numerator for Eff "<<name<<endl;
         hpneff3d = (TH3F*) efileDIR->Get(name);
      }
      if(strcmp((key->GetName()),(histD))==0){
         //cout<<"[getEffHists] Your denominator for Eff "<<name<<endl; 
         hpdeff3d = (TH3F*) efileDIR->Get(name);
      }
   }

   //efileDIR->Close();
   //efile->Close();

   return std::pair<TH3F*,TH3F*>(hpneff3d,hpdeff3d);
}


void PrintDetailsOfHist(TH3F* hSpectra){

   int  nbinX = hSpectra->GetNbinsX();
   int  nbinY = hSpectra->GetNbinsY();
   int  nbinZ = hSpectra->GetNbinsZ();

   cout<<"Number of bins in X (eta) : "<<nbinX
       <<" in Y (pt) : "<<nbinY
       <<" in Z (Jet ET) : "<<nbinZ<<endl;

   int binXmin = hSpectra->GetXaxis()->GetBinLowEdge(1);
   int binXmax = hSpectra->GetXaxis()->GetBinUpEdge(nbinX);

   int binYmin = hSpectra->GetYaxis()->GetBinLowEdge(1);
   int binYmax = hSpectra->GetYaxis()->GetBinUpEdge(nbinY);

   int binZmin = hSpectra->GetZaxis()->GetBinLowEdge(1);
   int binZmax = hSpectra->GetZaxis()->GetBinUpEdge(nbinZ);

   cout<<" eta ragnes[bin] : "<<binXmin<<" to " <<binXmax<<endl;
   cout<<" pt ragnes[bin] : "<<binYmin<<" to " <<binYmax<<endl;
   cout<<" jet ragnes[bin] : "<<binZmin<<" to " <<binZmax<<endl;

   /*
   int binMaxEta = hSpectra->GetXaxis()->FindBin(feta);
   int binMinEta = hSpectra->GetXaxis()->FindBin(-1.0*feta);

   cout<<"[DEBUG]-----"<<endl;
   cout<<"for input eta "<<feta<<" found max eta (bin center): "<<hSpectra->GetBinCenter(binMaxEta)<<endl;
   cout<<"for input eta "<<-1.0*feta<<" found min eta (bin center): "<<hSpectra->GetBinCenter(binMinEta)<<endl;
   cout<<"[DEBUG]-----"<<endl;

   float maxEta = hSpectra->GetXaxis()->GetBinCenter(binMaxEta);
   float minEta =hSpectra->GetXaxis()->GetBinCenter(binMinEta);
   */
}



