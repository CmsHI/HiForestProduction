//-------------------------------------------------------------------
// June 07 2010, Andre S. Yoon
// - simple GEN spectra extractor for different pT hat sample
// 
// 
//--------------------------------------------------------------------


#include <utility>
#include <iostream>

//#include "invar_yield_ana_GEN_v9.h"

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

using namespace std;


struct invar_yield_ana_GEN_v9_data {
   TH1D *hInvX;
   TH1D *hRInvX; // rebinned   
   TGraphErrors *InvX;
   TGraphErrors *RInvX; // rebinned   
   double integratedLum;
};

const float gpi = 3.14159265358979323846;
//typedef pair<TH3D*, TH3D*> 3DHistPair;

const int GNHIST = 10;

//-------------------------------------------------------------------------------------------------------------------------------
// 1. 
// 2. 
// 3. 
//------------------------------------------------------------------------------------------------------------------------------- 

invar_yield_ana_GEN_v9_data  invar_yield_ana_GEN_v9_graph(char *cFile1, char *cFile2, char *cFile3, char *cFile4,
							  char *cFile5, char *cFile6, char *cFile7, char *cFile8, char *cFile9, char *cFile10,
							  bool cross, double scale, double ieta, double feta)
   
{

   char* cFileArray[GNHIST];
   cFileArray[0] = cFile1, cFileArray[1] = cFile2, cFileArray[2] = cFile3, cFileArray[3] = cFile4;
   cFileArray[4] = cFile5, cFileArray[5] = cFile6, cFileArray[6] = cFile7, cFileArray[7] = cFile8;
   cFileArray[8] = cFile9, cFileArray[9] = cFile10;

   TH3F *hSpect3D_pre[GNHIST];
   TH1D *hSpect1D_pre[GNHIST];
   TH1F *hEVT[GNHIST];

   char nameHist[200];
   sprintf(nameHist,"hGenTrkPtEtaJetEtW");

   float xsection[10];
   xsection[0] = 4.844e+10, xsection[1] = 5.794e+08, xsection[2] = 2.361e+08, xsection[3] = 5.311e+07;
   xsection[4] = 6.358e+06, xsection[5] = 7.849e+05, xsection[6] = 1.151e+05, xsection[7] = 2.014e+04;
   xsection[8] = 4.094e+03, xsection[9] = 9.346e+02;

   cout<<"\n"<<endl;
   cout<<"[GEN spectra extracting..]======================================="<<endl;
   for(int i=0;i<GNHIST;i++){
      TFile *fEVT = new TFile(cFileArray[i]);
      hSpect3D_pre[i] = (TH3F*) fEVT->Get("preTrackAna/threeDHist/hGenTrkPtEtaJetEtW");
      hEVT[i] = (TH1F*) fEVT->Get("preTrackAna/hGenNevt");

      int binMaxEta = hSpect3D_pre[i]->GetXaxis()->FindBin(feta);
      int binMinEta = hSpect3D_pre[i]->GetXaxis()->FindBin(-1.0*feta);
      
      hSpect1D_pre[i] = (TH1D*) hSpect3D_pre[i]->ProjectionY("",binMinEta,binMaxEta,
							     hSpect3D_pre[i]->GetZaxis()->GetFirst(),
							     hSpect3D_pre[i]->GetZaxis()->GetLast(),"e");    
      cout<<i<<" : cross section = "<<xsection[i]<<" and num evt = "<<hEVT[i]->GetEntries()<<endl;
      hSpect1D_pre[i]->Sumw2();
      if(cross) hSpect1D_pre[i]->Scale((xsection[i]*1E-9)/(float)hEVT[i]->GetEntries());//to convert it to mb
      else hSpect1D_pre[i]->Scale(1.0/(float)hEVT[i]->GetEntries());
      if(i>0) hSpect1D_pre[0]->Add(hSpect1D_pre[i]);
   }

   hSpect1D_pre[0]->Scale(1./scale);

   int  nbinX = hSpect1D_pre[0]->GetNbinsX();
   double deta = (feta) - (ieta);
   cout<<"delta Eta = "<<deta<<endl;

   for(int j=0;j<(nbinX);j++){

      int xbin = j;

      double  pt   = hSpect1D_pre[0]->GetBinCenter(xbin+1); 
      double  dn   = hSpect1D_pre[0]->GetBinContent(xbin+1);
      double  edn  = hSpect1D_pre[0]->GetBinError(xbin+1);
      double  dbin = hSpect1D_pre[0]->GetXaxis()->GetBinWidth(xbin+1);

      dn = dn/(2*gpi*deta*dbin);
      edn=  edn/(2*gpi*deta*dbin);
      
      hSpect1D_pre[0]->SetBinContent(xbin+1,dn);
      hSpect1D_pre[0]->SetBinError(xbin+1,edn);
   }

   
   invar_yield_ana_GEN_v9_data data;

   data.hInvX = (TH1D*) hSpect1D_pre[0]->Clone("hInvX");
   TH1D *hInvX_dum_pre = (TH1D*) data.hInvX->Clone("hInvX_dum_pre");
   data.InvX = (TGraphErrors*) TgraphIt(hInvX_dum_pre);

   TH1D *hInvX_dum = (TH1D*) data.hInvX->Clone("hInvX_dum");
   data.hRInvX = (TH1D*) RebinIt(hInvX_dum,true);       

   TH1D *hRInvX_dum = (TH1D*) data.hRInvX->Clone("hRInvX_dum");
   data.RInvX = (TGraphErrors*) TgraphIt(hRInvX_dum);


   cout<<"[GEN spectra extracting..]======================================="<<endl;
   cout<<"\n"<<endl;

   return (data);

}


