//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------


#if !defined(__CINT__) || defined(__MAKECINT__)

#include "GetEffAndFake.h"

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


TFile *trkcorF1=0;
TFile *outputFile=0;

TDirectory *dTrkEff=0;

char cListTRK1[100];
char outFile[900];

char hist1[100], hist2[100], hist3[100], hist4[100];
TH3F *hPNEff3D=0, *hPDEff3D=0, *hPNFake3D=0, *hPDFake3D;
std::pair<TH3F*,TH3F*> h3DEff, h3DFake;

TH1D *hPNEff1D=0, *hPDEff1D=0, *hPNFake1D=0, *hPDFake1D=0;
TH1D *hPNEff1D_pre=0, *hPDEff1D_pre=0, *hPNFake1D_pre=0, *hPDFake1D_pre=0;

TGraphAsymmErrors *gEfficiency=0, *gFakerate=0;
int cbini=0, cbinf=0;

TString cbin;
TString outfile;

//--------------------------------------------------------------------------
void checkEtaRange(float iEta, float fEta, int EtaMin, int EtaMax);
void checkEtRange(float iJet, float fJet, int EtMin, int EtMax);
TH1D *rebinVarBin(TH1D *hist, TString hisname);
std::pair<TH3F*,TH3F*> getEffHists(const char *file, const char *dirC,
                                   const char *histN, const char *histD);
//--------------------------------------------------------------------------

void GetEffAndFakeNSave(const char *cDir="../root_files/",
			const char *cDirOut="./rootOutput/",
			const char *cFileTRK1="HiTrkHistMC_july08_hydjet_small_v2",
			const char *dir_corr="trkEffAnalyzer",
			float ijet=0, float fjet=2000,
			float ieta=0, float feta=1.0,
			float ipt=0,
			bool drawFig =true, bool saveFile=true, 
			int var=1, int rebF=1,
			int cbini_in=0, int cbinf_in=0,
			bool effonly=false
			)
{
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TH3::SetDefaultSumw2();

   cbini = cbini_in;
   cbinf = cbinf_in;


   // var = 1 => as a function of pT
   // var = 2 => as a function of whatever variable in z axis (jet ET, multiplicity, etc..)
   // var = 3 => as a function of eta

   cout<<"\n"<<"\n"<<endl;
   outfile = Form("%sTRKEffnFake_%s_%s_eta_%1.1fto%1.1f_jet_%1.1fto%1.1f_minpT%1.1f_reb%d_var%d_cbin%ito%i.root",                                                  
		  cDirOut,cFileTRK1,dir_corr,ieta,feta,ijet,fjet,ipt,rebF,var,cbini,cbinf);
   cout<<"Output file will be "<<outfile.Data()<<endl;
   cout<<"\n"<<endl;

   sprintf(cListTRK1,"%s%s.root",cDir,cFileTRK1);

   trkcorF1 = (TFile*) loadFile(trkcorF1,cListTRK1);
   sprintf(hist1,"heff3D"), sprintf(hist2,"hsim3D"), sprintf(hist3,"hfak3D"), sprintf(hist4,"hrec3D");



   if(cbini!=cbinf) {
     cbin.Append(Form("_cbin%dto%d",cbini,cbinf));
     cout<<"for histogram binned in centrality"<<cbin.Data()<<endl;
     
     if(cbini!=cbinf) {
       strcat(hist1,cbin.Data());
       strcat(hist2,cbin.Data());
       strcat(hist3,cbin.Data());
       strcat(hist4,cbin.Data());
     }
   }

   h3DEff = getEffHists(cListTRK1,dir_corr,hist1,hist2);
   h3DFake = getEffHists(cListTRK1,dir_corr,hist3,hist4);

   hPNEff3D = h3DEff.first, hPDEff3D = h3DEff.second;
   hPNFake3D = h3DFake.first, hPDFake3D = h3DFake.second;


   // eta and jet range for projection!
   float small = 0.01; // when the eta value is at the edge (1.0 and 2.4 is at the edges!)     

   int eta_min_bin = (ieta<=hPNEff3D->GetXaxis()->GetXmin()) ?
      hPNEff3D->GetXaxis()->GetFirst() : hPNEff3D->GetXaxis()->FindBin(-1.0*feta+small);
   int eta_max_bin = (feta>=hPNEff3D->GetXaxis()->GetXmax()) ?
      hPNEff3D->GetXaxis()->GetLast() : hPNEff3D->GetXaxis()->FindBin(feta-small);

   //eta_min_bin = eta_min_bin+1; //to make it symmetric!

   checkEtaRange(ieta,feta,eta_min_bin,eta_max_bin);

   int jet_min_bin = (ijet<=hPNEff3D->GetZaxis()->GetXmin()) ?
      hPNEff3D->GetZaxis()->GetFirst() : hPNEff3D->GetZaxis()->FindBin(ijet);
   int jet_max_bin = (fjet>=hPNEff3D->GetZaxis()->GetXmax()) ?
      hPNEff3D->GetZaxis()->GetLast() : hPNEff3D->GetZaxis()->FindBin(fjet-small);

   checkEtRange(ijet,fjet,jet_min_bin,jet_max_bin);

   int pt_min_bin = (ipt<=hPNEff3D->GetYaxis()->GetXmin()) ?
      hPNEff3D->GetYaxis()->GetFirst() : hPNEff3D->GetYaxis()->FindBin(ipt);
   int pt_max_bin = hPNEff3D->GetYaxis()->GetLast();


   if(var==1){
      hPNEff1D_pre = (TH1D*) hPNEff3D->ProjectionY("hPNEff1D_pre",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin,"e");
      hPDEff1D_pre = (TH1D*) hPDEff3D->ProjectionY("hPDEff1D_pre",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin,"e");
      hPNFake1D_pre = (TH1D*) hPNFake3D->ProjectionY("hPNFake1D_pre",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin,"e");
      hPDFake1D_pre = (TH1D*) hPDFake3D->ProjectionY("hPDFake1D_pre",eta_min_bin,eta_max_bin,jet_min_bin,jet_max_bin,"e");
   }else if(var==2){
      hPNEff1D_pre = (TH1D*) hPNEff3D->ProjectionZ("hPNEff1D_pre",eta_min_bin,eta_max_bin,pt_min_bin,pt_max_bin,"e");
      hPDEff1D_pre = (TH1D*) hPDEff3D->ProjectionZ("hPDEff1D_pre",eta_min_bin,eta_max_bin,pt_min_bin,pt_max_bin,"e");
      hPNFake1D_pre = (TH1D*) hPNFake3D->ProjectionZ("hPNFake1D_pre",eta_min_bin,eta_max_bin,pt_min_bin,pt_max_bin,"e");
      hPDFake1D_pre = (TH1D*) hPDFake3D->ProjectionZ("hPDFake1D_pre",eta_min_bin,eta_max_bin,pt_min_bin,pt_max_bin,"e");
   }else if(var==3){
     hPNEff1D_pre = (TH1D*) hPNEff3D->ProjectionX("hPNEff1D_pre",pt_min_bin,pt_max_bin,jet_min_bin,jet_max_bin,"e");
     hPDEff1D_pre = (TH1D*) hPDEff3D->ProjectionX("hPDEff1D_pre",pt_min_bin,pt_max_bin,jet_min_bin,jet_max_bin,"e");
     hPNFake1D_pre = (TH1D*) hPNFake3D->ProjectionX("hPNFake1D_pre",pt_min_bin,pt_max_bin,jet_min_bin,jet_max_bin,"e");
     hPDFake1D_pre = (TH1D*) hPDFake3D->ProjectionX("hPDFake1D_pre",pt_min_bin,pt_max_bin,jet_min_bin,jet_max_bin,"e");
   }

   if(rebF==0){ // arbitray variable binning
     hPNEff1D = (TH1D*) rebinVarBin(hPNEff1D_pre,"hPNEff1D");
     hPDEff1D = (TH1D*) rebinVarBin(hPDEff1D_pre,"hPDEff1D");
     hPNFake1D = (TH1D*) rebinVarBin(hPNFake1D_pre,"hPNFake1D");
     hPDFake1D = (TH1D*) rebinVarBin(hPDFake1D_pre,"hPDFake1D");
   }else if(rebF!=1){
     hPNEff1D_pre->Rebin(rebF,"hPNEff1D");
     hPDEff1D_pre->Rebin(rebF,"hPDEff1D");
     hPNFake1D_pre->Rebin(rebF,"hPNFake1D");
     hPDFake1D_pre->Rebin(rebF,"hPDFake1D");

     hPNEff1D = (TH1D*) gDirectory->Get("hPNEff1D");
     hPDEff1D = (TH1D*) gDirectory->Get("hPDEff1D");
     hPNFake1D = (TH1D*) gDirectory->Get("hPNFake1D");
     hPDFake1D = (TH1D*) gDirectory->Get("hPDFake1D");
   }else{
     hPNEff1D = (TH1D*) hPNEff1D_pre->Clone("hPNEff1D");
     hPDEff1D = (TH1D*) hPDEff1D_pre->Clone("hPDEff1D");
     hPNFake1D = (TH1D*) hPNFake1D_pre->Clone("hPNFake1D");
     hPDFake1D = (TH1D*) hPDFake1D_pre->Clone("hPDFake1D");
   }

   

   gEfficiency = new TGraphAsymmErrors();
   gEfficiency->BayesDivide(hPNEff1D,hPDEff1D);
   //gEfficiency->SetPointError(0,0.0,0.0,dYlow,dYhigh);

   if(!effonly){
     gFakerate =  new TGraphAsymmErrors();
     gFakerate->BayesDivide(hPNFake1D,hPDFake1D);
     //gFakerate->SetPointError(0,0.0,0.0,dYlowFR,dYhighFR);
   }

   if(saveFile){
     outputFile = new TFile(outfile.Data(),"recreate");
     dTrkEff = outputFile->mkdir("trkEffNFake");
     dTrkEff->cd();
     
     hPNEff1D->Write(), hPDEff1D->Write();
     if(!effonly) hPNFake1D->Write(),hPDFake1D->Write();
     
     
     gEfficiency->Write("gEfficiency");
     if(!effonly) gFakerate->Write("gFakerate");
     
     cout<<outfile.Data()<<" is created"<<endl;
   }



}

TH1D *rebinVarBin(TH1D *hist, TString hisname){
  
  TH1D *hist_cl1 = (TH1D*) hist->Clone("hist_cl1");
  TH1D *hist_cl2 = (TH1D*) hist->Clone("hist_cl2");
  TH1D *hist_cl3 = (TH1D*) hist->Clone("hist_cl3");

  hist_cl1->Rebin(2); //
  hist_cl2->Rebin(3); //
  hist_cl3->Rebin(6); //
  
  int bin = 1;
  double bins[36];
  int nbins = 0;

  bins[0] = 0.0;

  // for pt<16.8
  for(int i=0;i<hist_cl1->GetNbinsX();i++){
    double pt_up = hist_cl1->GetXaxis()->GetBinUpEdge(i+1);
    if(pt_up<17) {
      bins[i+1] = pt_up;
      nbins++;
    }
  }

  // for 16.8<pt<55.2
  for(int i=0;i<hist_cl2->GetNbinsX();i++){
    double pt_up = hist_cl2->GetXaxis()->GetBinUpEdge(i+1);
    if(pt_up<56 && pt_up>16.8) {
      bins[nbins+1] = pt_up;
      nbins++;
    }
  }

  // for 55.2<pt<inf
  for(int i=0;i<hist_cl3->GetNbinsX();i++){
    double pt_up = hist_cl3->GetXaxis()->GetBinUpEdge(i+1);
    if(pt_up>55.2) {
      bins[nbins+1] = pt_up; 
      nbins++;
    }
  }

  //for(int i=0;i<nbins+1;i++){
  //cout<<"pt bins = "<<bins[i]<<endl;
  //}
  
  
  TH1D *histnew = new TH1D("histnew","rebinned histo",nbins,bins);   
  histnew->SetName(hisname.Data());
  
  for(int i=0;i<histnew->GetNbinsX();i++){
    //if(i==0) cout<<"pt bins = "<<histnew->GetXaxis()->GetBinLowEdge(i+1)<<endl;
    //else cout<<"pt bins = "<<histnew->GetXaxis()->GetBinUpEdge(i+1)<<endl;
    double pt_up = histnew->GetXaxis()->GetBinUpEdge(i+1);
    double pt = histnew->GetXaxis()->GetBinCenter(i+1);
    if(pt_up<17) {
      histnew->SetBinContent(i+1,hist_cl1->GetBinContent(hist_cl1->FindBin(pt)));
      histnew->SetBinError(i+1,hist_cl1->GetBinError(hist_cl1->FindBin(pt)));
    }else if(pt_up<56 && pt_up>16.8){
      histnew->SetBinContent(i+1,hist_cl2->GetBinContent(hist_cl2->FindBin(pt)));
      histnew->SetBinError(i+1,hist_cl2->GetBinError(hist_cl2->FindBin(pt)));
    }else{
      histnew->SetBinContent(i+1,hist_cl3->GetBinContent(hist_cl3->FindBin(pt)));
      histnew->SetBinError(i+1,hist_cl3->GetBinError(hist_cl3->FindBin(pt)));
    }
  }

  return histnew;
  //return hist;

  /*
  Int_t bin = 1;
  Double_t bins[80];
  Int_t nbins = 0;
  TAxis *xaxis = h1->GetXaxis();
  while (bin <100) {
    bins[nbins] = xaxis->GetBinLowEdge(bin);
    nbins++;
    if (bin < 50)      bin += 2;
    else if (bin < 80) bin += 5;
    else               bin += 10;
  }             
  bins[nbins] = xaxis->GetBinUpEdge(100);
   
  TH1F *h2 = new TH1F("h2","rebinned histo",nbins,bins);
  for (bin=1;bin<100;bin++) {
    h2->Fill(xaxis->GetBinCenter(bin),h1->GetBinContent(bin));
  }
*/
}


void checkEtRange(float iJet, float fJet, int EtMin, int EtMax){

   cout<<"[checkEtRange]--------------------"<<endl;
   cout<<"for input min Et = "<<iJet<<endl;
   cout<<"jet min Et (low edge) = "<<hPNEff3D->GetZaxis()->GetBinLowEdge(EtMin)<<endl;
   cout<<"jet min Et (up edge) = "<<hPNEff3D->GetZaxis()->GetBinUpEdge(EtMin)<<endl;
   cout<<"for input max Et = "<<fJet<<endl;
   cout<<"jet min Et (low edge) = "<<hPNEff3D->GetZaxis()->GetBinLowEdge(EtMax)<<endl;
   cout<<"jet min Et (up edge) = "<<hPNEff3D->GetZaxis()->GetBinUpEdge(EtMax)<<endl;
   cout<<"[checkEtRange]--------------------"<<endl;
   cout<<"\n"<<endl;

}

void checkEtaRange(float iEta, float fEta, int EtaMin, int EtaMax){

   cout<<"\n"<<endl;
   cout<<"[checkEtaRange]-------------------------------------------------"<<endl;
   cout<<"for input eta "<<fEta<<" found max eta (bin center): "<<hPNEff3D->GetXaxis()->GetBinCenter(EtaMax)<<endl;
   cout<<"for input eta "<<-1.0*fEta<<" found min eta (bin center): "<<hPNEff3D->GetXaxis()->GetBinCenter(EtaMin)<<endl;
   cout<<"for input eta "<<fEta<<" found max eta (bin up edge): "<<hPNEff3D->GetXaxis()->GetBinUpEdge(EtaMax)<<endl;
   cout<<"for input eta "<<-1.0*fEta<<" found min eta (bin low edge): "<<hPNEff3D->GetXaxis()->GetBinLowEdge(EtaMin)<<endl;
   cout<<"[checkEtaRange]-------------------------------------------------"<<endl;
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
         hpneff3d = (TH3F*) efileDIR->Get(name);
      }
      if(strcmp((key->GetName()),(histD))==0){
         hpdeff3d = (TH3F*) efileDIR->Get(name);
      }
   }

   //efileDIR->Close();                                                                                                                                                                                       
   //efile->Close();                                                                                                                                                                                          
   return std::pair<TH3F*,TH3F*>(hpneff3d,hpdeff3d);
}
