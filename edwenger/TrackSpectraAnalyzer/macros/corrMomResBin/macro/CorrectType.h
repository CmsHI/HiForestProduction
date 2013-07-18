#ifndef CORRECTTYPE_H
#define CORRECTTYPE_H


#include <utility>
#include <iostream>

#include <TROOT.h>
#include <TStyle.h>
#include "TError.h"

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TF1.h"

#include "TFile.h"
#include "TMath.h"

using namespace std;

//----------------------------------------------------------------------------------
TFile *loadFile(TFile *file, char *cfile);
TFile *loadFile(TFile *file, char *cfile, char *option);
TH1D *makeItInvHist(TH1D *hmIIH, float deta);
//----------------------------------------------------------------------------------


const double small = 1e-3;
double ptb, ptbV;
double fact = 4.0;
double fact2 = 5.0;
std::vector<double> ptBins;
std::vector<double> ptBinsV;

void prepareVarBin(){

   // same bin as in TrackSpectraAnalyzer
   for(ptb =   0; ptb <   1-small; ptb +=  0.05) ptBins.push_back(ptb); // 20 bins -> 24 bins?
   for(ptb =   1; ptb <   2-small; ptb +=  0.1 ) ptBins.push_back(ptb); // 10 bins -> 12 bins
   for(ptb =   2; ptb <   6-small; ptb +=  0.2 ) ptBins.push_back(ptb); // 20 bins -> 24 bins?
   for(ptb =   6; ptb <  10-small; ptb +=  0.5 ) ptBins.push_back(ptb); // 8 bins  -> 12 bins?
   for(ptb =  10; ptb <  20-small; ptb +=  1.0 ) ptBins.push_back(ptb); // 10 bins ->12
   for(ptb =  20; ptb <  50-small; ptb +=  3.0 ) ptBins.push_back(ptb); // 10 bins ->12
   for(ptb =  50; ptb < 100-small; ptb +=  5.0 ) ptBins.push_back(ptb); // 10 bins ->12
   for(ptb = 100; ptb < 200-small; ptb += 10.0 ) ptBins.push_back(ptb); // 10 bins ->12
   for(ptb = 200; ptb < 400-small; ptb += 20.0 ) ptBins.push_back(ptb); // 10 bins ->12

   

   // rebinned!
   for(ptbV =   0; ptbV <   1-small; ptbV +=  (0.05*fact)) ptBinsV.push_back(ptbV); // 20 bins
   for(ptbV =   1; ptbV <   2-small; ptbV +=  (0.1*fact2) ) ptBinsV.push_back(ptbV); // 10 bins ->
   for(ptbV =   2; ptbV <   6-small; ptbV +=  (0.2*fact) ) ptBinsV.push_back(ptbV); // 20 bins
   for(ptbV =   6; ptbV <  10-small; ptbV +=  (0.5*fact) ) ptBinsV.push_back(ptbV); // 8 bins 
   for(ptbV =  10; ptbV <  20-small; ptbV +=  (1.0*fact2) ) ptBinsV.push_back(ptbV); // 10 bins
   for(ptbV =  20; ptbV <  50-small; ptbV +=  (3.0*fact2) ) ptBinsV.push_back(ptbV); // 10 bins
   for(ptbV =  50; ptbV < 100-small; ptbV +=  (5.0*fact2) ) ptBinsV.push_back(ptbV); // 10 bins
   for(ptbV = 100; ptbV < 200-small; ptbV += (10.0*fact2) ) ptBinsV.push_back(ptbV); // 10 bins
   for(ptbV = 200; ptbV < 400-small; ptbV += (20.0*fact2) ) ptBinsV.push_back(ptbV); // 10 bins



   cout<<"[prepareVarBin] variable pt bins are prepared!"<<endl;
   cout<<"\n"<<endl;
}

TH1D *rebinWitCorrection(TH1D *hRWC_pre){
   //------
   /*
   double vbins[ptBinsV.size()];
   for(int i=0;i<ptBinsV.size();i++){
      vbins[i] = ptBinsV[i];
      cout<<"vbins = "<<vbins[i]<<endl;
   }
   vbins[ptBinsV.size()] = 400; // needed to close the end
   cout<<"vbins = "<<vbins[ptBinsV.size()]<<endl;
   */
   //--------

   TH1D *hRWC = (TH1D*) hRWC_pre->Clone("hRWC");
   cout<<"[rebinWitCorrection] number of bins = "<<hRWC->GetNbinsX()
       <<" && after rebin = "<<ptBinsV.size()-1<<endl;
   TH1D *hRWC_rebinned = (TH1D*) hRWC->Rebin(ptBinsV.size()-1,"hRWC_rebinned", &ptBinsV[0]);
   //TH1D *hRWC_rebinned = (TH1D*) hRWC->Rebin(ptBinsV.size(),"hRWC_rebinned", vbins); 
   
   for(int i = 0; i<hRWC_rebinned->GetNbinsX();i++){

      double pt = hRWC_rebinned->GetBinCenter(i+1);
      double content = hRWC_rebinned->GetBinContent(i+1);
      double err = hRWC_rebinned->GetBinError(i+1);

      if(pt>=0 && pt< 1) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact), hRWC_rebinned->SetBinError(i+1,err/fact);
      else if(pt>=1 && pt< 2) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact2), hRWC_rebinned->SetBinError(i+1,err/fact2);
      else if(pt>=2 && pt< 6) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact), hRWC_rebinned->SetBinError(i+1,err/fact);
      else if(pt>=6 && pt< 10) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact), hRWC_rebinned->SetBinError(i+1,err/fact);
      else if(pt>=10 && pt< 20) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact2), hRWC_rebinned->SetBinError(i+1,err/fact2);
      else if(pt>=20 && pt< 50) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact2), hRWC_rebinned->SetBinError(i+1,err/fact2);
      else if(pt>=50 && pt< 100) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact2), hRWC_rebinned->SetBinError(i+1,err/fact2); 
      else if(pt>=100 && pt< 200) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact2), hRWC_rebinned->SetBinError(i+1,err/fact2);
      else if(pt>=200 && pt< 400) 
	 hRWC_rebinned->SetBinContent(i+1,content/fact2), hRWC_rebinned->SetBinError(i+1,err/fact2);

   }
   return hRWC_rebinned;
}

TFile *loadFile(TFile *file, char *cfile, char *option){

   if(file!=0) delete file;
   file = new TFile(cfile,option);
   if(file->IsOpen()==false) Error("[FileLoading]","could not open hist file [%s].",cfile);
   else cout<<"[FileLoading] input hist file:"<<cfile<<" with an option "<<option<<endl;
   return file;

}

TFile *loadFile(TFile *file, char *cfile){

   //receive a pointer and return a pointer! 
   if(file!=0) delete file;
   file = new TFile(cfile);
   if(file->IsOpen()==false) Error("[FileLoading]","could not open hist file [%s].",cfile);
   else cout<<"[FileLoading] input hist file:"<<cfile<<endl;
   return file;

}

TH1D *makeItInvHist(TH1D *hmIIH, float deta){
   //TH1D *hmIIH = (TH1D*) hmIIH_pre->Clone("hmIIH"); // is this necessary?

   for(int i=0;i<hmIIH->GetNbinsX();i++){
      double pt = hmIIH->GetBinCenter(i+1);
      double dndpt = hmIIH->GetBinContent(i+1);
      double err = hmIIH->GetBinError(i+1);
      hmIIH->SetBinContent(i+1,dndpt*(1/(2*TMath::Pi()*pt*deta)));
      hmIIH->SetBinError(i+1,err*(1/(2*TMath::Pi()*pt*deta)));
   }
   return hmIIH;
}



#endif
