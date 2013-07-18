//-------------------------------------------------------------------------------------
// June 30 2010 Andre S. Yoon
//
// This macro reads the output of TrackSpectraAnalyzer and extract dN/dpT for a given
// eta range, and correct for trk inefficiency based on the outputs of TrkEffAnalyzer.
// Various histograms are saved in a seperate root file, which in turn feeds the 
// next type correction, i.e CorrectTypeTwoNSave.C
//
// MC vs DATA
//    MC   -> a) RECO  
//         -> b) GEN (by isGEN flag)
//    DATA -> a) RECO
//
// CORRECTION LEVEL CONVENTION:
//    1. lev1 -> trk
//    2. lev2 -> trk+sec
//    3. full -> trk+sec+mult
//-------------------------------------------------------------------------------------


#if !defined(__CINT__) || defined(__MAKECINT__)

#include "CorrectType.h"

#include <utility>
#include <iostream>

#include "utilities.h"
#include "commonStyle.h"

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

#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TChain.h"
#include "TGraph.h"

#include "TKey.h"
using namespace std;

#endif

TH3F *hdndptdetadet=0;
TH3F *hdndptdetadet_raw=0;
TH3F *hdndptdetadet_lev1=0;
TH3F *hdndptdetadet_lev2=0;
TH3F *hdndptdetadet_full=0;

TH1D *hdndptR_raw_dum=0;

TH1D *hdndpt_raw=0, *hdndptR_raw=0;
TH1D *hdndpt_full=0, *hdndptR_full=0;
TH1D *hdndpt_lev1=0, *hdndptR_lev1=0;
TH1D *hdndpt_lev2=0, *hdndptR_lev2=0;

TH1D *hInv_raw=0, *hInvR_raw=0;
TH1D *hInv_full=0, *hInvR_full=0;
TH1D *hInv_lev1=0, *hInvR_lev1=0;
TH1D *hInv_lev2=0, *hInvR_lev2=0;

TGraphErrors *tdndpt_raw=0, *tdndptR_raw=0;
TGraphErrors *tdndpt_full=0, *tdndptR_full=0;
TGraphErrors *tdndpt_lev1=0, *tdndptR_lev1=0;
TGraphErrors *tdndpt_lev2=0, *tdndptR_lev2=0;

TGraphErrors *tInv_raw=0, *tInvR_raw=0;
TGraphErrors *tInv_full=0, *tInvR_full=0;
TGraphErrors *tInv_lev1=0, *tInvR_lev1=0;
TGraphErrors *tInv_lev2=0, *tInvR_lev2=0;


TFile *spectraF=0;
TFile *trkcorF1=0, *trkcorF2=0,*trkcorF3=0,*trkcorF4=0,*trkcorF5=0;

TFile *outputFile=0;
TDirectory *dTypeOne=0;

const int NHIST = 5;
const char *cFileTRKArray[NHIST];

char hist1[100], hist2[100], hist3[100], hist4[100];
TH3F *hPNEff3D[NHIST], *hPDEff3D[NHIST], *hPNFake3D[NHIST], *hPDFake3D[NHIST];
std::pair<TH3F*,TH3F*> h3DEff[NHIST], h3DFake[NHIST];

char histMult1[100], histMult2[100];
TH3F *hPNMult3D[NHIST], *hPDMult3D[NHIST];
std::pair<TH3F*,TH3F*> h3DMult[NHIST];

char histSec1[100], histSec2[100];
TH3F *hPNSec3D[NHIST], *hPDSec3D[NHIST];
std::pair<TH3F*,TH3F*> h3DSec[NHIST];

char xTitle[100],yTitle[100];

TCanvas *call, *callR, *callRR;
TH1D *dum, *dumR, *dumRR;

char cList[100];
char cListTRK1[100], cListTRK2[100], cListTRK3[100];
char cListTRK4[100], cListTRK5[100];
char outFile[300];

//------------------------------------------------------------------------------
void checkEtaRange(float iEta, float fEta, int EtaMin, int EtaMax);
void checkEtRange(float iJet, float fJet, int EtMin, int EtMax);
void prepareTrkEffCorrection(const char *dirCorr, const char *cTRK1,
			     const char *cTRK2, const char *cTRK3,
			     const char *cTRK4, const char *cTRK5);
std::pair<TH3F*,TH3F*> getEffHists(const char *file, const char *dirC, 
				   const char *histN, const char *histD);
double GetEffFactor(TH3F* num, TH3F* den, double pt, double eta, double jet);
double GetFakFactor(TH3F* num, TH3F* den, double pt, double eta, double jet);
void GetEffFactorDEBUG(TH3F* num, TH3F* den, double pt, double eta, double jet);
void drawDebugPlots();
//-----------------------------------------------------------------------------


void CorrectTypeOneNSave(const char *cDir="../root_files/",
			 const char *cDirOut="./rootOutput/",
                         const char *cFile="TrkHistMC_june09_qcdMB",
                         const char *cFileTRK1="TrkHistMC_june09_qcdMB",
                         const char *cFileTRK2="TrkHistMC_june09_qcdMB",
                         const char *cFileTRK3="TrkHistMC_june09_qcdMB",
                         const char *cFileTRK4="TrkHistMC_june09_qcdMB",
                         const char *cFileTRK5="TrkHistMC_june09_qcdMB",
                         const char *dir_ana="trackAna_STD",
                         const char *dir_corr="trkEffAnalyzer",
                         bool isGEN=false, bool varBin=false,
                         float ijet=0, float fjet=2000,
                         float ieta=0, float feta=2.4,
			 double scale_dum=1,
			 bool debug = true,
                         bool drawFig =true, bool saveFile=false
                         )
{

   cout<<"\n"<<"\n"<<endl;

   //double scaleF = 1.13643e+07;
   //double scaleF = 3.79542505524861872e+08;

   //double scaleF = 1.1805E+07;
   //double scaleF = 4.0786E+08;
   //double scaleF = 5.5917e+08;
   double scaleF = 1;

   //double scaleF = 2.3845e+07; //total N_nsd (MB in pr6 and m6rr)
   //double scaleF = 5.5917e+08; 

   cout<<">>>> scale factor = "<<scaleF<<endl;
   
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TH3::SetDefaultSumw2();


   sprintf(outFile,"%sCORR_TEST1_%s_eta_%1.1fto%1.1f_jet_%1.1fto%1.1f_%s_GEN%i_varBin%i.root",cDirOut,cFile,ieta,feta,ijet,fjet,dir_ana,isGEN,varBin);
   cout<<"Output file will be "<<outFile<<endl;
   cout<<"\n"<<endl;

   sprintf(cList,"%s%s.root",cDir,cFile);
   sprintf(cListTRK1,"%s%s.root",cDir,cFileTRK1);
   sprintf(cListTRK2,"%s%s.root",cDir,cFileTRK2);
   sprintf(cListTRK3,"%s%s.root",cDir,cFileTRK3);
   sprintf(cListTRK4,"%s%s.root",cDir,cFileTRK4);
   sprintf(cListTRK5,"%s%s.root",cDir,cFileTRK5);

   spectraF = (TFile*) loadFile(spectraF,cList);
   trkcorF1 = (TFile*) loadFile(trkcorF1,cListTRK1);
   trkcorF2 = (TFile*) loadFile(trkcorF2,cListTRK2);
   trkcorF3 = (TFile*) loadFile(trkcorF3,cListTRK3);
   trkcorF4 = (TFile*) loadFile(trkcorF4,cListTRK4);
   trkcorF5 = (TFile*) loadFile(trkcorF5,cListTRK5);

   char histName[50];

   if(isGEN) sprintf(histName,"hGenTrkPtEtaJetEt");
   else sprintf(histName,"hTrkPtEtaJetEt");

   if(varBin) strcat(histName,"_vbin");

   char histStr[100];
   sprintf(histStr,"%s/threeDHist/%s",dir_ana,histName);
   cout<<"\n"<<endl;
   cout<<"spectra hist to be looked at = "<<histStr<<endl;

   hdndptdetadet = (TH3F*) spectraF->Get(histStr);
   if(hdndptdetadet==0) {
      Error("[HistLoading]","could not get histogram"); 
      return;
   }


   hdndptdetadet_raw = (TH3F*) hdndptdetadet->Clone("hdndptdetadet_raw");
   hdndptdetadet_raw->Reset();

   hdndptdetadet_lev1 = (TH3F*) hdndptdetadet->Clone("hdndptdetadet_lev1");
   hdndptdetadet_lev1->Reset();

   hdndptdetadet_lev2 = (TH3F*) hdndptdetadet->Clone("hdndptdetadet_lev2");
   hdndptdetadet_lev2->Reset();

   hdndptdetadet_full = (TH3F*) hdndptdetadet->Clone("hdndptdetadet_full");
   hdndptdetadet_full->Reset();

   // First decide eta and jet et range -----------------------------
   float small = 0.01; // when the eta value is at the edge (1.0 and 2.4 is at the edges!)  

   int eta_min_bin = (ieta<=hdndptdetadet->GetXaxis()->GetXmin()) ?
      hdndptdetadet->GetXaxis()->GetFirst() : hdndptdetadet->GetXaxis()->FindBin(-1.0*feta+small);
   int eta_max_bin = (feta>=hdndptdetadet->GetXaxis()->GetXmax()) ?
      hdndptdetadet->GetXaxis()->GetLast() : hdndptdetadet->GetXaxis()->FindBin(feta-small);
   
   checkEtaRange(ieta,feta,eta_min_bin,eta_max_bin);

   int jet_min_bin = (ijet<=hdndptdetadet->GetZaxis()->GetXmin()) ? 
      hdndptdetadet->GetZaxis()->GetFirst() : hdndptdetadet->GetZaxis()->FindBin(ijet+small); 
   int jet_max_bin = (fjet>=hdndptdetadet->GetZaxis()->GetXmax()) ?
      hdndptdetadet->GetZaxis()->GetLast() : hdndptdetadet->GetZaxis()->FindBin(fjet-small);

   checkEtRange(ijet,fjet,jet_min_bin,jet_max_bin);


   // Loop over the 3D histograms and apply correction --------------
   if(!isGEN) prepareTrkEffCorrection(dir_corr, cListTRK1,
			   cListTRK2,cListTRK3,
			   cListTRK4,cListTRK5); // first prepare relevant histograms for corrections

   int  nbinX = hdndptdetadet->GetNbinsX();
   int  nbinY = hdndptdetadet->GetNbinsY();
   int  nbinZ = hdndptdetadet->GetNbinsZ();
   
   cout<<"number of bin in x axis = "<<nbinX<<" in y axis = "<<nbinY<<" in z axis = "<<nbinZ<<endl;

   // jet et range for different samples for trk correction
   // make sure jet et range is quantized with min et range of 20 GeV!
   // so the range has to be at n*20, where n = 1,2,3...
   float jet_1st, jet_2nd, jet_3rd, jet_4th;
   //jet_1st = 41, jet_2nd = 61, jet_3rd = 81, jet_4th = 181;  
   jet_1st = 41, jet_2nd = 61, jet_3rd = 101, jet_4th = 181;

   float var_1st, var_2nd, var_3rd, var_4th; // var can be pt, Et, etc..

   // Below is for the purpose of debugging, i.e. to find how many tracks are
   // not corrected due to the lack of statistics or different fragmentation in the correction MC samples!
   float recTracktotal = 0.0; // total rec tracks with a certain kinematic cut!
   float recTrackNotCorrectedEff = 0.0; // total rec tracks with a certain kinematic cut but not corrected for eff!
   float recTrackNotCorrectedFak = 0.0;
   float pt_thres = 30.0; //

   int recTracktotalBin = 0;
   int recTrackNotCorrectedEffBin = 0;
   int recTrackNotCorrectedFakBin = 0;

   
   for(int j=0;j<(nbinX*nbinY*nbinZ);j++){
      
      int xbin = (j % nbinX);
      int ybin = ((j-xbin)/nbinX) % nbinY;
      int zbin = (j-(xbin+((ybin)*nbinX)))/(nbinX*nbinY);
      
      float  eta = hdndptdetadet->GetXaxis()->GetBinCenter(xbin+1);
      float  pt = hdndptdetadet->GetYaxis()->GetBinCenter(ybin+1);
      float  jet = hdndptdetadet->GetZaxis()->GetBinCenter(zbin+1);
      
      float  dbin = hdndptdetadet->GetYaxis()->GetBinWidth(ybin+1);
   
      float  dn = hdndptdetadet->GetBinContent(xbin+1,ybin+1,zbin+1);
      float  edn = hdndptdetadet->GetBinError(xbin+1,ybin+1,zbin+1);
      
      double eff=0, fak=0;
      double mlt=0, sec=0;
      
      float var_for_corr = 0;
      
      var_for_corr = jet;
      var_1st = jet_1st, var_2nd = jet_2nd, var_3rd = jet_3rd, var_4th = jet_4th;
      
      if(!isGEN){  // for GEN truth, no correction is needed!
	 if(var_for_corr<var_1st){
	    eff = GetEffFactor(hPNEff3D[0],hPDEff3D[0],pt,eta,jet);
	    fak = GetFakFactor(hPNFake3D[0],hPDFake3D[0],pt,eta,jet);
	    mlt = GetFakFactor(hPNMult3D[0],hPDMult3D[0],pt,eta,jet);
	    sec = GetFakFactor(hPNSec3D[0],hPDSec3D[0],pt,eta,jet);
	    if(dn>0 && pt>pt_thres && fabs(eta)<feta) GetEffFactorDEBUG(hPNEff3D[0],hPDEff3D[0],pt,eta,jet);
	 }else if(var_for_corr>=var_1st && var_for_corr<var_2nd){
	    eff = GetEffFactor(hPNEff3D[1],hPDEff3D[1],pt,eta,jet);
	    fak = GetFakFactor(hPNFake3D[1],hPDFake3D[1],pt,eta,jet);
	    mlt = GetFakFactor(hPNMult3D[1],hPDMult3D[1],pt,eta,jet);
	    sec = GetFakFactor(hPNSec3D[1],hPDSec3D[1],pt,eta,jet);
	    if(dn>0 && pt>pt_thres && fabs(eta)<feta) GetEffFactorDEBUG(hPNEff3D[1],hPDEff3D[1],pt,eta,jet);
	 }else if(var_for_corr>=var_2nd && var_for_corr<var_3rd){
	    eff = GetEffFactor(hPNEff3D[2],hPDEff3D[2],pt,eta,jet);
	    fak = GetFakFactor(hPNFake3D[2],hPDFake3D[2],pt,eta,jet);
	    mlt = GetFakFactor(hPNMult3D[2],hPDMult3D[2],pt,eta,jet);
	    sec = GetFakFactor(hPNSec3D[2],hPDSec3D[2],pt,eta,jet);
	    if(dn>0 && pt>pt_thres && fabs(eta)<feta) GetEffFactorDEBUG(hPNEff3D[2],hPDEff3D[2],pt,eta,jet);
	 }else if(var_for_corr>=var_3rd && var_for_corr<var_4th){
	    eff = GetEffFactor(hPNEff3D[3],hPDEff3D[3],pt,eta,jet);
	    fak = GetFakFactor(hPNFake3D[3],hPDFake3D[3],pt,eta,jet);
	    mlt = GetFakFactor(hPNMult3D[3],hPDMult3D[3],pt,eta,jet);
	    sec = GetFakFactor(hPNSec3D[3],hPDSec3D[3],pt,eta,jet);
	    if(dn>0 && pt>pt_thres && fabs(eta)<feta) GetEffFactorDEBUG(hPNEff3D[3],hPDEff3D[3],pt,eta,jet);
	 }else{
	    eff = GetEffFactor(hPNEff3D[4],hPDEff3D[4],pt,eta,jet);
	    fak = GetFakFactor(hPNFake3D[4],hPDFake3D[4],pt,eta,jet);
	    mlt = GetFakFactor(hPNMult3D[4],hPDMult3D[4],pt,eta,jet);
	    sec = GetFakFactor(hPNSec3D[4],hPDSec3D[4],pt,eta,jet);
	    if(dn>0 && pt>pt_thres && fabs(eta)<feta) GetEffFactorDEBUG(hPNEff3D[4],hPDEff3D[4],pt,eta,jet);
	 }
      }

      if(dn>0 && pt>pt_thres && fabs(eta)<feta){
	recTracktotalBin++, recTracktotal = recTracktotal + dn;
	if(eff==1 || eff==0){
	  recTrackNotCorrectedEffBin++, recTrackNotCorrectedEff = recTrackNotCorrectedEff + dn;
	  //if(pt>(jet+10.)) cout<<"[Not corrected for EFF and Pt>Jet] pt = "<<pt<<" and jet = "<<jet<<endl; 
	  //cout<<"[Not corrected for EFF]"<<" pt = "<<pt<<" and eta = "<<eta<<" and jet = "<<jet<<endl; 
	}
	if(fak==1 || fak==0){
	  recTrackNotCorrectedFakBin++, recTrackNotCorrectedFak = recTrackNotCorrectedFak + dn;
	  //if(pt>(jet+10.)) cout<<"[Not corrected for EFF and Pt>Jet] pt = "<<pt<<" and jet = "<<jet<<endl;
	  //cout<<"[Not corrected for FR]"<<" pt = "<<pt<<" and eta = "<<eta<<" and jet = "<<jet<<endl;
	}
      }

      // dn/dpt
      dn = dn/(dbin);
      edn = edn/(dbin);
      
      hdndptdetadet_raw->SetBinContent(xbin+1,ybin+1,zbin+1,dn);
      hdndptdetadet_raw->SetBinError(xbin+1,ybin+1,zbin+1,edn);
      
      if(!isGEN){  // for GEN truth, no correction is needed!
	 
	 // efficiency and fake rate ! 
	 if(eff==0 || fak ==1){
	    dn = dn;
	    edn = edn;
	 }else{
	    dn = dn*((1-fak)/eff);
	    edn = edn*((1-fak)/eff);
	 }
	 
	 hdndptdetadet_lev1->SetBinContent(xbin+1,ybin+1,zbin+1,dn);
	 hdndptdetadet_lev1->SetBinError(xbin+1,ybin+1,zbin+1,edn);
	 
	 // secondary !
	 if(sec==1 || sec==0){
	    dn = dn, edn = edn; // no correction!                                                                                                        
	 }else{
	    dn = dn*(1-sec), edn = edn*(1-sec);
	 }
	 
	 hdndptdetadet_lev2->SetBinContent(xbin+1,ybin+1,zbin+1,dn);
	 hdndptdetadet_lev2->SetBinError(xbin+1,ybin+1,zbin+1,edn);
	 
	 // multiple reconstruction !
	 if(mlt==1 || mlt==0){
	    dn = dn, edn = edn; // no correction! 
	 }else{
	    dn = dn*(1./(1+mlt)), edn = edn*(1./(1+mlt));
	 }
	 
	 hdndptdetadet_full->SetBinContent(xbin+1,ybin+1,zbin+1,dn);
	 hdndptdetadet_full->SetBinError(xbin+1,ybin+1,zbin+1,edn);
      }
   }

   cout<<"\n"<<endl;
   cout<<"[Debugging summary]==================================================="<<endl;
   cout<<"number of total tracks = "<<recTracktotal<<" above pT = "<<pt_thres<<" (GeV/c)"<<endl;
   cout<<"number of total tracks no corrected for eff = "<<recTrackNotCorrectedEff<<endl;
   cout<<"fraction = "<<(recTrackNotCorrectedEff/recTracktotal)<<endl;
   cout<<"number of total tracks no corrected for fak = "<<recTrackNotCorrectedFak<<endl;
   cout<<"fraction = "<<(recTrackNotCorrectedFak/recTracktotal)<<endl;
   cout<<""<<endl;
   cout<<"number of total bins = "<<recTracktotalBin<<" above pT = "<<pt_thres<<" (GeV/c)"<<endl;
   cout<<"number of total bins no corrected for eff = "<<recTrackNotCorrectedEffBin<<endl;
   cout<<"fraction = "<<((float)recTrackNotCorrectedEffBin)/((float)recTracktotalBin)<<endl;
   cout<<"number of total bins no corrected for fak = "<<recTrackNotCorrectedFakBin<<endl;
   cout<<"fraction = "<<((float)recTrackNotCorrectedFakBin)/((float)recTracktotalBin)<<endl;
   cout<<"[Debugging summary]==================================================="<<endl;
   cout<<"\n"<<endl;

			  
   // to avoid root projeciton bug
   // ProjectionY is buggy when x axis full x-axis range (-2.4 to 2.4)
   // However, the problem does not show up when the 3D histogram is reset, which is the case in this macro
   // , but to make sure in any case, hackedProjectionY is used.
   hdndpt_raw = (TH1D*) hackedProjectionY(hdndptdetadet_raw,"hdndpt_raw",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin);
   hdndpt_lev1 = (TH1D*) hackedProjectionY(hdndptdetadet_lev1,"hdndpt_lev1",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin);
   hdndpt_lev2 = (TH1D*) hackedProjectionY(hdndptdetadet_lev2,"hdndpt_lev2",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin);
   hdndpt_full = (TH1D*) hackedProjectionY(hdndptdetadet_full,"hdndpt_full",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin);

   hdndpt_raw->Scale(1./scaleF);
   hdndpt_lev1->Scale(1./scaleF);
   hdndpt_lev2->Scale(1./scaleF);
   hdndpt_full->Scale(1./scaleF);



   hInv_raw = (TH1D*) hdndpt_raw->Clone("hInv_raw");
   hInv_lev1 = (TH1D*) hdndpt_lev1->Clone("hInv_lev1");
   hInv_lev2 = (TH1D*) hdndpt_lev2->Clone("hInv_lev2");
   hInv_full = (TH1D*) hdndpt_full->Clone("hInv_full");
   
   hInv_raw = (TH1D*) makeItInvHist(hInv_raw, 2.*feta);
   hInv_lev1 = (TH1D*) makeItInvHist(hInv_lev1, 2.*feta);
   hInv_lev2 = (TH1D*) makeItInvHist(hInv_lev2, 2.*feta);
   hInv_full = (TH1D*) makeItInvHist(hInv_full, 2.*feta);
   
   tdndpt_raw =  (TGraphErrors*) TgraphIt(hdndpt_raw);
   tdndpt_lev1 = (TGraphErrors*) TgraphIt(hdndpt_lev1);
   tdndpt_lev2 = (TGraphErrors*) TgraphIt(hdndpt_lev2);
   tdndpt_full = (TGraphErrors*) TgraphIt(hdndpt_full);

   tInv_raw =  (TGraphErrors*) TgraphIt(hInv_raw);
   tInv_lev1 = (TGraphErrors*) TgraphIt(hInv_lev1);
   tInv_lev2 = (TGraphErrors*) TgraphIt(hInv_lev2);
   tInv_full = (TGraphErrors*) TgraphIt(hInv_full);

   
   if(drawFig){
      drawDebugPlots();
   }
   
   if(saveFile){
      outputFile = new TFile(outFile,"recreate");   
      dTypeOne = outputFile->mkdir("corrTypeOne");
      dTypeOne->cd();                                

      hdndpt_raw->Write(), hdndpt_lev1->Write();
      hdndpt_lev2->Write(), hdndpt_full->Write();

      hInv_raw->Write(), hInv_lev1->Write();
      hInv_lev2->Write(), hInv_full->Write();

      tdndpt_raw->Write("tdndpt_raw"),tdndpt_lev1->Write("tdndpt_lev1");
      tdndpt_lev2->Write("tdndpt_lev2"),tdndpt_full->Write("tdndpt_full");

      tInv_raw->Write("tInv_raw"),tInv_lev1->Write("tInv_lev1");
      tInv_lev2->Write("tInv_lev2"),tInv_full->Write("tInv_full");

      cout<<outFile<<" is created"<<endl;
   }
}


void drawDebugPlots(){

   call = new TCanvas("call","call",310,470);
   call->cd(), call->SetLogx(), call->SetLogy();
   sprintf(yTitle,"dN/dp_{T}");
   sprintf(xTitle,"p_{T} [GeV/c]");
   dum = GetDummyHist(200,1E-11,1E+2,xTitle,yTitle);
   dum->Draw();

   th1Style1(hdndpt_raw,14,28,1.0,14,1.0,1,1);
   th1Style1(hdndpt_lev1,15,24,1.0,15,1.0,1,1);
   th1Style1(hdndpt_lev2,16,25,1.0,16,1.0,1,1);
   th1Style1(hdndpt_full,17,26,1.0,17,1.0,1,1);

   th1Style1(tdndpt_raw,14,28,1.0,14,1.0,1,2);
   th1Style1(tdndpt_lev1,15,24,1.0,15,1.0,1,2);
   th1Style1(tdndpt_lev2,16,25,1.0,16,1.0,1,2);
   th1Style1(tdndpt_full,17,26,1.0,17,1.0,1,2);

   callR = new TCanvas("callR","callR",310,470);
   callR->cd(), callR->SetLogy();
   dumR = (TH1D*) dum->Clone("dumR");
   dumR->Draw();

   th1Style1(hdndpt_raw,14,28,1.0,14,1.0,1,1);
   th1Style1(hdndpt_lev1,15,24,1.0,15,1.0,1,1);
   th1Style1(hdndpt_lev2,16,25,1.0,16,1.0,1,1);
   th1Style1(hdndpt_full,17,26,1.0,17,1.0,1,1);

   th1Style1(tdndpt_raw,14,28,1.0,14,1.0,1,2);
   th1Style1(tdndpt_lev1,15,24,1.0,15,1.0,1,2);
   th1Style1(tdndpt_lev2,16,25,1.0,16,1.0,1,2);
   th1Style1(tdndpt_full,17,26,1.0,17,1.0,1,2);

}

void prepareTrkEffCorrection(const char *dirCorr, const char *cTRK1,const char *cTRK2,
			     const char *cTRK3,const char *cTRK4,const char *cTRK5)
{
   sprintf(hist1,"heff3D"), sprintf(hist2,"hsim3D"), sprintf(hist3,"hfak3D"), sprintf(hist4,"hrec3D");
   sprintf(histMult1,"hmul3D"), sprintf(histMult2,"hsim3D");
   sprintf(histSec1,"hsec3D"), sprintf(histSec2,"hrec3D");

   for(int i=0;i<NHIST;i++){

      if(i==0) cFileTRKArray[i] = cTRK1;
      else if(i==1) cFileTRKArray[i] = cTRK2;
      else if(i==2) cFileTRKArray[i] = cTRK3;
      else if(i==3) cFileTRKArray[i] = cTRK4;
      else if(i==4) cFileTRKArray[i] = cTRK5;

      // trk eff/fr correction                                                                                                                           
      h3DEff[i] = getEffHists(cFileTRKArray[i],dirCorr,hist1,hist2);
      h3DFake[i] = getEffHists(cFileTRKArray[i],dirCorr,hist3,hist4);

      hPNEff3D[i] = h3DEff[i].first, hPDEff3D[i] = h3DEff[i].second;
      hPNFake3D[i] = h3DFake[i].first, hPDFake3D[i] = h3DFake[i].second;

      // multiply reconstructed!                                                                                                                         
      h3DMult[i] = getEffHists(cFileTRKArray[i],dirCorr,histMult1,histMult2);
      hPNMult3D[i] = h3DMult[i].first, hPDMult3D[i] = h3DMult[i].second;

      // secondary !                                                                                                                                     
      h3DSec[i] = getEffHists(cFileTRKArray[i],dirCorr,histSec1,histSec2);
      hPNSec3D[i] = h3DSec[i].first, hPDSec3D[i] = h3DSec[i].second;
   }
}

/*
TFile *loadFile(TFile *file, char *cfile){

   //receive a pointer and return a pointer!
   if(file!=0) delete file;
   file = new TFile(cfile); 
   if(file->IsOpen()==false) Error("[FileLoading]","could not open hist file [%s].",cfile); 
   else cout<<"[FileLoading] input hist file:"<<cfile<<endl; 
   return file;
}
*/


void checkEtaRange(float iEta, float fEta, int EtaMin, int EtaMax){

   cout<<"\n"<<endl;
   cout<<"[checkEtaRange]-------------------------------------------------"<<endl;
   cout<<" bin number eta min = "<<EtaMin<<" eta max = "<<EtaMax<<endl;
   cout<<"for input eta "<<fEta<<" found max eta (bin center): "<<hdndptdetadet->GetXaxis()->GetBinCenter(EtaMax)<<endl;
   cout<<"for input eta "<<-1.0*fEta<<" found min eta (bin center): "<<hdndptdetadet->GetXaxis()->GetBinCenter(EtaMin)<<endl;
   cout<<"for input eta "<<fEta<<" found max eta (bin up edge): "<<hdndptdetadet->GetXaxis()->GetBinUpEdge(EtaMax)<<endl;
   cout<<"for input eta "<<-1.0*fEta<<" found min eta (bin low edge): "<<hdndptdetadet->GetXaxis()->GetBinLowEdge(EtaMin)<<endl;
   cout<<"Integration range is from "<<hdndptdetadet->GetXaxis()->GetBinLowEdge(EtaMin)<<" to "<<hdndptdetadet->GetXaxis()->GetBinUpEdge(EtaMax)<<endl;
   cout<<"[checkEtaRange]-------------------------------------------------"<<endl;
   cout<<"\n"<<endl;

}

void checkEtRange(float iJet, float fJet, int EtMin, int EtMax){

   cout<<"[checkEtRange]--------------------"<<endl;
   cout<<" bin number jet min = "<<EtMin<<" jet max = "<<EtMax<<endl;
   cout<<"for input min Et = "<<iJet<<endl;
   cout<<"jet min Et (low edge) = "<<hdndptdetadet->GetZaxis()->GetBinLowEdge(EtMin)<<endl;
   cout<<"jet min Et (up edge) = "<<hdndptdetadet->GetZaxis()->GetBinUpEdge(EtMin)<<endl;
   cout<<"for input max Et = "<<fJet<<endl;
   cout<<"jet max Et (low edge) = "<<hdndptdetadet->GetZaxis()->GetBinLowEdge(EtMax)<<endl;
   cout<<"jet max Et (up edge) = "<<hdndptdetadet->GetZaxis()->GetBinUpEdge(EtMax)<<endl;
   cout<<"Integration range is from "<<hdndptdetadet->GetZaxis()->GetBinLowEdge(EtMin)<<" to "<<hdndptdetadet->GetZaxis()->GetBinUpEdge(EtMax)<<endl;
   cout<<"[checkEtRange]--------------------"<<endl;
   cout<<"\n"<<endl;

}


std::pair<TH3F*,TH3F*> getEffHists(const char *file, const char *dirC, 
				   const char *histN, const char *histD){

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


void GetEffFactorDEBUG(TH3F* num, TH3F* den, double pt, double eta, double jet){

  int num_bin = num->FindBin(eta,pt,jet);
  int den_bin = den->FindBin(eta,pt,jet);

  double n_num = num->GetBinContent(num_bin);
  double n_den = den->GetBinContent(den_bin);

  if(n_den == 0) cout<<"[GetEffFactorDEBUG den=0] hist = "<<num->GetName()<<" pt = "<<pt<<" eta = "<<eta<<" jet = "<<jet<<endl; 
  if(n_num == 0) cout<<"[GetEffFactorDEBUG num=0] hist = "<<den->GetName()<<" pt = "<<pt<<" eta = "<<eta<<" jet = "<<jet<<endl;

  //if(n_den == 0  || n_num == 0 ) return 1;
  //else return n_num/n_den;
}


double GetEffFactor(TH3F* num, TH3F* den, double pt, double eta, double jet){

   int num_bin = num->FindBin(eta,pt,jet);
   int den_bin = den->FindBin(eta,pt,jet);

   //double num_bin_width = num->GetYaxis()->GetBinWidth(num_bin); 
   //double den_bin_width = den->GetYaxis()->GetBinWidth(den_bin); 

   double n_num = num->GetBinContent(num_bin);
   double n_den = den->GetBinContent(den_bin);

   /*
   if(n_den == 0  || n_num == 0 ) {
     if(fabs(eta)<2.4){ 
	 cout<<"[GetEffFactor] correction factor 1 will be returned as n_den = "<<n_num
	     <<" and n_den = "<<n_den<<endl;
	 cout<<" pt = "<<pt<<" jet = "<<jet<<" and eta = "<<eta<<endl;
     }
   }
   */

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

