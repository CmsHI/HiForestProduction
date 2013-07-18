//---------------------------------------------------------------------------------------------
// July 1 2010 Andre S. Yoon   
//
// This macro reads output file of CorrectTypeOneNSave.C and apply the mom. res.,
// rebinning, and 1/pT correction (*pT in the center of bin is not nessarily 
// same as the center of gravity in the bin with the steeply falling distribution)
// (Rebinning+1/pT) collectively bin width correction!
// And save output into the existing input file!
//
// PROCEDURE
//   1. input dN/dpT -> fit -> analytically obtain (1/2pipTdeta)(dN/dpT)    => "TRUTH"
//   2. input dN/dpT -> fit -> fill histogram from fit (same bin as input dN/dpT)*
//                   -> [smear with response matrix (MC)*]                  => Smeared dN/dpT
//                   -> (1/2pipTdeta)(dN/dpT) to the smeared spectra 
//                   -> [rebinned]                                          => "MEASURED" (quote!) 
//   3. correction factor c = "TRUTH"/"MEASURED"                            => bare correction
//   4. c -> fit with polynomial                                            => fitted correction
//
// * steps in [ ] can be skipped 
// 
// CORRECTION LEVEL CONVENTION
//   1. lev1(R) -> bin width only(Rebinning) 
//   2. lev2(R) -> full (Rebinning)
//
// DETAILS - smearing
//   1. reco spectra in true (sim) pT bin is taken from the 2D hist from TrkEffAnalyzer
//   2. nomalize to unity (or any constant)
//   3. weight by the fit 
//   4. project weighted 2D hist on to measured (reco) pT bin 
// 
// Q - binning of response matrix?? smaller?
//   - different pt samples!
// 
// DEBUGGING PLOT
//   1. fit dN/dpT - input dN/pT, fit function, fit stat
//   2. fit invariant - fit invariant, weighted fit invariant, filled hist
//   3. smeared spectra
//   4. correction factors - a. rebin corr +  fit, b. full corr + fit
//   5. possible rec pt distribution in each sim pt (gaussin~)
//---------------------------------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "CorrectType.h"

#include <utility>
#include <iostream>

#include "/home/sungho/plots/common/utilities.h"
#include "/home/sungho/plots/common/commonStyle.h"

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

#include "TMath.h"

#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TChain.h"
#include "TGraph.h"

#include "TLine.h"

#include "TKey.h"
#include "TString.h"

using namespace std;

TF1 *fitF=0,#endif

char inFile[300];
char inFileMom[300];

TFile *spectraF=0;
TFile *momResF=0;

TH1D *hdNdpT=0;
TH1D *hfitF=0;
TH1D *hfitF_inv=0;
TF1 *fitF=0, *fitF_inv=0;

TH1D *hRes1D=0;
TH2D *hRes2D = 0;
TH3F *hRes3D = 0;

TH1D *hLev1Spec=0, *hLev1SpecR=0;
TH1D *hLev2Spec=0, *hLev2SpecR=0;

TH1D *hCorrLev1=0, *hCorrLev1R=0;
TH1D *hCorrLev2=0, *hCorrLev2R=0;

TH1D *hInvCorrOne=0, *hInvCorrOneR=0; // previous spectra!
TH1D *hInvCorrTwoLev1=0, *hInvCorrTwoLev1R=0;
TH1D *hInvCorrTwoLev2=0, *hInvCorrTwoLev2R=0;

TGraphErrors *tInvCorrOne=0, *tInvCorrOneR=0; // previous spectra!
TGraphErrors *tInvCorrTwoLev1=0, *tInvCorrTwoLev1R=0;
TGraphErrors *tInvCorrTwoLev2=0, *tInvCorrTwoLev2R=0;

TString fitFormula;

TCanvas *can1=0, *can2=0, *can3=0, *can4=0, *can5=0, *can6=0, *can7=0;
TH1D *dum1=0, *dum2=0, *dum3=0, *dum4=0, *dum5=0, *dum6=0, *dum7=0;

TDirectory *dTypeTwo=0;

//------------------------------------------------------------------------------
void doSmearing();
void fitdNdpT(float ix, float fx);
void obtainFitInv(float ix, float fx, float deta);
TH3F *makeItdNdpT(TH3F *h3Dhist);
int GetFirstNonEmptyBin(TH1D* hist);
int GetLastNonEmptyBin(TH1D* hist);
TH1D *correctedSpectra(TH1D *inputHist, TH1D *corrF);
//------------------------------------------------------------------------------ 

void CorrectTypeTwoNSave(const char *cDir="../root_files/",
			 const char *cFile="TrkHistMC_july04_qcdMBD6T0704vbinV2",
			 //const char *cFile="MB-C10-M6RR-MBskim-v0_proc0628_trkAnaNoFilter_v2_proc0628",
			 const char *cFileTRK1="TrkHistMC_july04_qcdPt170d0704vbinV2", 
			 //const char *cFileTRK1="TrkHistMC_july04_qcdMBD6T0704vbinV2",
			 const char *dir_ana="trackAna_STD",
                         const char *dir_corr="trkEffAnalyzer",
			 bool isGEN=false, bool varBin=true,
                         float ijet=0, float fjet=2000,
                         float ieta=0, float feta=2.4,
			 bool rebOnly=false,
                         bool drawFig =true, bool saveFile=true
			 )
{
   cout<<"\n"<<"\n"<<endl;

   prepareVarBin();

   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TH3::SetDefaultSumw2();
   
   
   // files and hist loading ---------------------------------------------
   sprintf(inFile,"./rootOutput/CORR_%s_eta_%1.1fto%1.1f_jet_%1.1fto%1.1f_%s_GEN%i_varBin%i.root",cFile,ieta,feta,ijet,fjet,dir_ana,isGEN,varBin);
   cout<<"Input file for spectra will be "<<inFile<<endl;
   cout<<"\n"<<endl;

   spectraF = (TFile*) loadFile(spectraF,inFile,"UPDATE");
   
   char histStr[100], histStr2[100];
   sprintf(histStr,"corrTypeOne/hdndpt_full");
   sprintf(histStr2,"corrTypeOne/hInv_full");
   
   cout<<"\n"<<endl;
   cout<<"spectra hist to be looked at = "<<histStr<<endl;
   cout<<"inv spectra hist to be looked at = "<<histStr2<<endl;

   hdNdpT = (TH1D*) spectraF->Get(histStr);
   hInvCorrOne = (TH1D*) spectraF->Get(histStr2);

   if(hdNdpT==0||hInvCorrOne==0) {
      Error("[HistLoading]","could not get histogram");
      return;
   }


   // fitting -----------------------------------------------------------
   float xmin = 0.3;
   float xmax = 1.1*((float) hdNdpT->GetXaxis()->GetXmax()); //10% more than maximum

   fitdNdpT(xmin, xmax);

   hfitF = FillFromFunc(hdNdpT,fitF);
   hfitF->Draw("histsame");

   can2 = new TCanvas("can2","can2",550,410);
   can2->cd();
   can2->SetLogy();

   dum2 = GetDummyHist(0.5,150,2E-9,2E+8,"p_{T}","d^{3}N/d^{3}p",true);
   dum2->Draw();
   
   obtainFitInv(xmin,xmax,2.*feta);

   hfitF_inv = FillFromFunc(hdNdpT,fitF_inv); // fill from function!   
   //hfitF_inv->Draw("histsame");

   // obtaining 3D hist for mom smearing! ----------------------------------

   if(!rebOnly){
      sprintf(inFileMom,"%s%s.root",cDir,cFileTRK1);
      cout<<"Input file for spectra will be "<<inFileMom<<endl;
      cout<<"\n"<<endl;
      momResF = (TFile*) loadFile(momResF,inFileMom);

      hRes3D = (TH3F*) momResF->Get("trkEffAnalyzer/hresStoR3D");
      hRes3D = (TH3F*) makeItdNdpT(hRes3D);

      int eta_min_bin = (ieta<=hRes3D->GetXaxis()->GetXmin()) ?
	 hRes3D->GetXaxis()->GetFirst() : hRes3D->GetXaxis()->FindBin(-1.0*feta);
      int eta_max_bin = (feta>=hRes3D->GetXaxis()->GetXmax()) ?
	 hRes3D->GetXaxis()->GetLast() : hRes3D->GetXaxis()->FindBin(feta);
      
      hRes3D->GetXaxis()->SetRange(eta_min_bin,eta_max_bin); 
      hRes2D = (TH2D*) hRes3D->Project3D("zy");

      if(hRes2D==0 || hRes3D==0) {
	 Error("[HistLoading]","could not get histogram");
	 return;
      }
   }else{
      cout<<"No momemtum resolution smearing is applied!"<<endl;
   }

   can6 = new TCanvas("can6","can6",550,410);
   can6->cd();
   hRes2D->Draw("COLZ");

   doSmearing();

   // "Measured" spectra -------------------------------------------------------
   hLev1Spec = (TH1D*) hfitF->Clone("hLev1Spec"); // rather than hLev1Spec = (TH1D*) hfitF_inv->Clone("hLev1Spec");
   hLev1Spec = (TH1D*) makeItInvHist(hLev1Spec,2.*feta);

   hLev2Spec = (TH1D*) hRes2D->ProjectionY("hLev2Spec",1,hRes2D->GetXaxis()->GetLast(),"e"); // projection of convoluted hist        
   hLev2Spec = (TH1D*) makeItInvHist(hLev2Spec,2.*feta);

   hLev1SpecR = (TH1D*) rebinWitCorrection(hLev1Spec);
   hLev2SpecR = (TH1D*) rebinWitCorrection(hLev2Spec);

   can2->cd();

   hLev1Spec->SetLineWidth(2), hLev1Spec->SetLineColor(2), hLev1Spec->SetLineStyle(2);
   hLev1Spec->SetMarkerStyle(24), hLev1Spec->SetMarkerColor(2), hLev1Spec->SetMarkerSize(1);
   hLev1Spec->Draw("pzsame");

   hLev1SpecR->SetLineWidth(2), hLev1SpecR->SetLineColor(2), hLev1SpecR->SetLineStyle(2);
   hLev1SpecR->SetMarkerStyle(22), hLev1SpecR->SetMarkerColor(2), hLev1SpecR->SetMarkerSize(1);
   hLev1SpecR->Draw("pzsame");
   
   hLev2Spec->SetLineWidth(2), hLev2Spec->SetLineColor(4), hLev2Spec->SetLineStyle(2);
   hLev2Spec->SetMarkerStyle(24), hLev2Spec->SetMarkerColor(4), hLev2Spec->SetMarkerSize(1);
   hLev2Spec->Draw("pzsame");

   hLev2SpecR->SetLineWidth(2), hLev2SpecR->SetLineColor(4), hLev2SpecR->SetLineStyle(2);
   hLev2SpecR->SetMarkerStyle(22), hLev2SpecR->SetMarkerColor(4), hLev2SpecR->SetMarkerSize(1);
   hLev2SpecR->Draw("pzsame");
   

   
   // Correction factor = "Measured"/"Truth" -------------------------------------
   hCorrLev1 = (TH1D*) ratio_hist_to_func(hLev1Spec,fitF_inv); // "Measured"/"Truth" 
   hCorrLev2 = (TH1D*) ratio_hist_to_func(hLev2Spec,fitF_inv); // "Measured"/"Truth"
   
   hCorrLev1R = (TH1D*) ratio_hist_to_func(hLev1SpecR,fitF_inv); // "Measured(Rebinned)"/"Truth"  
   hCorrLev2R = (TH1D*) ratio_hist_to_func(hLev2SpecR,fitF_inv); // "Measured(Rebinned)"/"Truth"
   
   can4 = new TCanvas("can4","can4",550,410);
   can4->cd();

   dum4 = GetDummyHist(0.5,150,0.3,2.3,"p_{T}^{rec}","Correction factor",true);
   dum4->Draw("");

   TLine *line = new TLine(0.5,1,150,1);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();

   hCorrLev1->SetLineWidth(2), hCorrLev1->SetLineColor(2), hCorrLev1->Draw("pzsame");
   hCorrLev1R->SetMarkerStyle(22), hCorrLev1R->SetLineWidth(2), hCorrLev1R->SetLineColor(2), hCorrLev1R->Draw("pzsame");

   hCorrLev2->SetLineWidth(2), hCorrLev2->SetLineColor(4), hCorrLev2->Draw("pzsame");
   hCorrLev2R->SetMarkerStyle(22),hCorrLev2R->SetLineWidth(2), hCorrLev2R->SetLineColor(4), hCorrLev2R->Draw("pzsame");

   
   // Corrected inv yield spectra!
   hInvCorrTwoLev1 = (TH1D*) correctedSpectra(hInvCorrOne,hCorrLev1);
   hInvCorrTwoLev2 = (TH1D*) correctedSpectra(hInvCorrOne,hCorrLev2);

   hInvCorrOneR = (TH1D*) rebinWitCorrection(hInvCorrOne);
   
   hInvCorrTwoLev1R = (TH1D*) correctedSpectra(hInvCorrOneR,hCorrLev1R);
   hInvCorrTwoLev2R = (TH1D*) correctedSpectra(hInvCorrOneR,hCorrLev2R);
   

   hInvCorrOne->SetName("hInvCorrOne");
   hInvCorrTwoLev1->SetName("hInvCorrTwoLev1");
   hInvCorrTwoLev2->SetName("hInvCorrTwoLev2");
   hInvCorrTwoLev1R->SetName("hInvCorrTwoLev1R");
   hInvCorrTwoLev2R->SetName("hInvCorrTwoLev2R");


   can7 = new TCanvas("can7","can7",550,410);
   can7->cd();
   can7->SetLogy();

   dum7 = GetDummyHist(0.5,150,2E-9,2E+8,"p_{T}","d^{3}N/d^{3}p",true);
   dum7->Draw();
   
   hInvCorrOne->SetMarkerStyle(3), hInvCorrOne->SetMarkerSize(1), hInvCorrOne->Draw("pzsame");
   hInvCorrTwoLev1->SetMarkerStyle(24), hInvCorrTwoLev1->SetMarkerSize(1), hInvCorrTwoLev1->SetMarkerColor(2), hInvCorrTwoLev1->Draw("pzsame");
   hInvCorrTwoLev2->SetMarkerStyle(24), hInvCorrTwoLev2->SetMarkerSize(1), hInvCorrTwoLev2->SetMarkerColor(4), hInvCorrTwoLev2->Draw("pzsame");
   hInvCorrTwoLev1R->SetMarkerStyle(22), hInvCorrTwoLev1R->SetMarkerSize(1), hInvCorrTwoLev1R->SetMarkerColor(2), hInvCorrTwoLev1R->Draw("pzsame");
   hInvCorrTwoLev2R->SetMarkerStyle(22), hInvCorrTwoLev2R->SetMarkerSize(1), hInvCorrTwoLev2R->SetMarkerColor(4), hInvCorrTwoLev2R->Draw("pzsame");

   tInvCorrOne = (TGraphErrors*) TgraphIt(hInvCorrOne);
   tInvCorrTwoLev1 = (TGraphErrors*) TgraphIt(hInvCorrTwoLev1);
   tInvCorrTwoLev2 = (TGraphErrors*) TgraphIt(hInvCorrTwoLev2);
   tInvCorrTwoLev1R = (TGraphErrors*) TgraphIt(hInvCorrTwoLev1R);
   tInvCorrTwoLev2R = (TGraphErrors*) TgraphIt(hInvCorrTwoLev2R);


   if(saveFile){
      dTypeTwo = spectraF->GetDirectory("corrTypeTwo");
      if(dTypeTwo==0){
	 dTypeTwo = spectraF->mkdir("corrTypeTwo");
	 dTypeTwo->cd();
	 cout<<"corrTypeTwo directory is created"<<endl;
      }else{
	 dTypeTwo->cd();
      }

      hInvCorrOne->Write("hInvCorrOne",TObject::kOverwrite);
      hInvCorrTwoLev1->Write("hInvCorrTwoLev1",TObject::kOverwrite);
      hInvCorrTwoLev2->Write("hInvCorrTwoLev2",TObject::kOverwrite);
      hInvCorrTwoLev1R->Write("hInvCorrTwoLev1R",TObject::kOverwrite);
      hInvCorrTwoLev2R->Write("hInvCorrTwoLev2R",TObject::kOverwrite);
      
      tInvCorrOne->Write("tInvCorrOne",TObject::kOverwrite);
      tInvCorrTwoLev1->Write("tInvCorrTwoLev1",TObject::kOverwrite);
      tInvCorrTwoLev2->Write("tInvCorrTwoLev2",TObject::kOverwrite);
      tInvCorrTwoLev1R->Write("tInvCorrTwoLev1R",TObject::kOverwrite);
      tInvCorrTwoLev2R->Write("tInvCorrTwoLev2R",TObject::kOverwrite);

      cout<<inFile<<" is updated"<<endl;
   }


   

}


TH1D *correctedSpectra(TH1D *inputHist_pre, TH1D *corrF_pre){

   TH1D *inputHist = (TH1D*) inputHist_pre->Clone("inputHist");
   TH1D *corrF = (TH1D*) corrF_pre->Clone("corrF");

   cout<<"[correctedSpectra] number of bins in input hist = "<<inputHist->GetNbinsX()
       <<" and in correction hist = "<<corrF->GetNbinsX()<<endl;

   if(inputHist->GetNbinsX()!=corrF->GetNbinsX()){
      Error("[correctSpectra]","input inv spectra and correction factor has different number of bins!!");
      return 0;
   }
   for(int i=0;i<inputHist->GetNbinsX();i++){
      double yield = inputHist->GetBinContent(i+1);
      double yielderr = inputHist->GetBinError(i+1);
      if(fabs(inputHist->GetBinCenter(i+1)-corrF->GetBinCenter(i+1))>0.1){
	 Error("[correctSpectra]","input inv spectra and correction factor has different pt bins!");
	 return 0;
      }
      double corrR = corrF->GetBinContent(i+1);
      //cout<<"corr factor ="<<corrR<<" for pt = "<<corrF->GetBinCenter(i+1)<<endl;
      if(corrR==0) corrR = 1; // corrR=0 when yield = 0;
      inputHist->SetBinContent(i+1,yield/corrR);
      inputHist->SetBinError(i+1,yielderr/corrR);
   }

   return inputHist;
}

void doSmearing(){
   //-----------------------------------------------------------------------
   // Q. does the binning of 2D response matrix have to be same as dN/dpT?
   // Yes. if the final correction is applied in the form of histogram
   // No.  if the final correction is applied as a fit function
   //-----------------------------------------------------------------------
   can3 = new TCanvas("can3","can3",550,410);
   can3->cd();
   can3->SetLogy();

   dum3 = (TH1D*) dum1->Clone("dum3");
   dum3->Draw();

   hRes1D = (TH1D*) hRes2D->ProjectionY("hRes1D",hRes2D->GetXaxis()->GetFirst(),hRes2D->GetXaxis()->GetLast(),"e");
   hRes1D->SetMarkerSize(1), hRes1D->Draw("pzsame"); //projected dn  

   can5 = new TCanvas("can5","can5",550,410);
   can5->cd();

   dum5 = GetDummyHist(0.5,150,0,2,"p_{T}^{rec}","Prob",true);
   dum5->SetMinimum(0.00001);
   dum5->Draw("");
   
   int nbin = hRes1D->GetNbinsX();
   cout<<"[doSmearing] nbin = "<<nbin<<endl;

   TH1D *hRes1D_temp=0;

   for(int i=0;i<nbin;i++){

      bool nullArea = false; // this is needed for emtpy mom res bin 

      //- project from SIM pt to REC pt
      /*
      if(nbin_dndpt!= hRes2D->GetNbinsY()){
	 Error("[doSmearing]","");
         return;
      }
      */

      hRes1D_temp = (TH1D*) hRes2D->ProjectionY("",i+1,i+1,"e");
      //TH1D *hRes1D_temp  = (TH1D*) hRes2D->ProjectionY("hRes1D_temp",i+1,i+1,"e"); 

      //- normalize the distribution
      double area = hRes1D_temp->Integral();

      if(area==0) nullArea =true;
      else hRes1D_temp->Scale(1./area);
      
      hRes1D_temp->SetLineColor((i+1)%10);
      hRes1D_temp->Draw("histsame");

      int firstB = GetFirstNonEmptyBin(hRes1D_temp);
      int lastB = GetLastNonEmptyBin(hRes1D_temp);

      //cout<<"first bin = "<<firstB<<" last bin = "<<lastB<<endl;

      double corr;
      
      //firstB=1, lastB=2;
      for(int j=firstB;j<lastB+1;j++){

	 if(!nullArea) corr = hRes1D_temp->GetBinContent(j);
	 else corr = 1.0; // no information on res-> no correction                                                                                       

	 //corr = 1.0; 

	 double ptr = hRes1D_temp->GetBinCenter(j);

	 int rpt_bin = hRes2D->GetYaxis()->FindBin(ptr); // should be same as j!                                                                         
	 int spt_bin = i+1;

	 int bin_spec = hfitF->GetXaxis()->FindBin(ptr); // binning between Mom Res Hist and Spectra may differ
	 double smeared = corr*(hfitF->GetBinContent(bin_spec));// weighted by filled histogram (from fit)
	 //double smeared = corr*(fitF->Eval(ptr)); // weighted by the fit                                                                             
	 
	 hRes2D->SetBinContent(spt_bin,rpt_bin,smeared);
	 hRes2D->SetBinError(spt_bin,rpt_bin,0);
      }
      delete hRes1D_temp;
   }

}


TH3F *makeItdNdpT(TH3F *h3Dhist){
   // bit computing intensive 
   // -> better do it for 2D?

   //TH3F *h3Dhist = (TH3F*) h3Dhist_pre->Clone("h3Dhist");
   //h3Dhist->Reset();

   int  nbinX = h3Dhist->GetNbinsX();
   int  nbinY = h3Dhist->GetNbinsY();
   int  nbinZ = h3Dhist->GetNbinsZ();
   cout<<"number of x bin = "<<nbinX<<" y bin = "<<nbinY<<" z bin = "<<nbinZ<<endl;

   for(int j=0;j<(nbinX*nbinY*nbinZ);j++){

      int xbin = (j % nbinX);
      int ybin = ((j-xbin)/nbinX) % nbinY;
      int zbin = (j-(xbin+((ybin)*nbinX)))/(nbinX*nbinY);

      // dbin from Y axis creates a bummy structure at boundary!
      float  dbin = h3Dhist->GetZaxis()->GetBinWidth(zbin+1);
      //float  dbin = h3Dhist->GetYaxis()->GetBinWidth(ybin+1);

      float  dn = h3Dhist->GetBinContent(xbin+1,ybin+1,zbin+1);
      float  edn = h3Dhist->GetBinError(xbin+1,ybin+1,zbin+1);

      // dN/dpT
      dn = dn/(dbin);
      edn = edn/(dbin);

      h3Dhist->SetBinContent(xbin+1,ybin+1,zbin+1,dn);
      h3Dhist->SetBinError(xbin+1,ybin+1,zbin+1,edn);
   }

   return h3Dhist;
}


void obtainFitInv(float ix, float fx, float deta){

   fitFormula += "*(1/(2*TMath::Pi()*x*[6]))";
   cout<<"\n"<<"[obtainFitInv] fitting function (invaraint form) = "<<fitFormula.Data()<<"\n"<<endl;

   fitF_inv = new TF1("fitFunctionInv",fitFormula.Data(),ix,fx);
   fitF_inv->SetParameters(fitF->GetParameter(0),fitF->GetParameter(1),fitF->GetParameter(2),
                           fitF->GetParameter(3),fitF->GetParameter(4),fitF->GetParameter(5),
                           deta);

   fitF_inv->SetLineColor(4), fitF_inv->SetLineWidth(2);
   //fitF_inv->Draw("same");

}

void fitdNdpT(float ix, float fx){
   
   can1 = new TCanvas("can1","can1",550,410);
   can1->cd();
   can1->SetLogy();

   
   dum1 = GetDummyHist(0.5,150,2E-4,1E+9,"p_{T}","dN/dp_{T}",true);
   dum1->Draw();
   
   fitFormula = "[0]*([5]+(x/[1])+(pow(x,[2])/[3]))^[4]";

   cout<<"\n"<<"[fitdNdpT] fitting function = "<<fitFormula.Data()<<"\n"<<endl;
   
   fitF = new TF1("fitFunction",fitFormula.Data(),ix,fx);

   fitF->SetParameters(2.96361e+08,2.66339e-01,1,1,-5,1);

   fitF->SetLineColor(2);
   fitF->SetLineWidth(3);

   hdNdpT->Fit(fitF,"NRO");

   hdNdpT->SetMarkerSize(1), hdNdpT->Draw("pzsame");
   fitF->Draw("same");
}

int GetFirstNonEmptyBin(TH1D* hist){
   int bin = 0;
   for(int i=0;i<hist->GetNbinsX();i++){
      bin = i+1;
      if(hist->GetBinContent(bin)>0) break;
   }
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
