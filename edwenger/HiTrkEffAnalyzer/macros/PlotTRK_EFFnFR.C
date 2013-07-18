


#if !defined(__CINT__) || defined(__MAKECINT__)

#include "GetEffAndFake.h"

#include "../plotting/utilities.h"
#include "../plotting/commonStyle.h"

#include "TLine.h"
#include "TLegend.h"
#endif

TGraphAsymmErrors *gEff_v1=0, *gFake_v1=0;
TGraphAsymmErrors *gEff_v2=0, *gFake_v2=0;
TGraphAsymmErrors *gEff_v3=0, *gFake_v3=0;
TGraphAsymmErrors *gEff_v4=0, *gFake_v4=0;
TGraphAsymmErrors *gEff_v5=0, *gFake_v5=0;

char fileName1[200], fileName2[200], fileName3[200]; 
char fileName4[200], fileName5[200];


TFile *trkf_v1=0, *trkf_v2=0, *trkf_v3=0, *trkf_v4=0, *trkf_v5=0;


TPad *pp1=0, *pp1_1=0;
TH1D *dum1=0, *dum1_1=0;

TLegend *legd=0;

bool lowpt_eta;  // low pt cut for eff(eta)

//-------------------------------------------------
//void setPad(TCanvas *can);
//void setLowerPad(TH1D *dum_low);
//TH1D *rerangedHist(TH1D *hRange, int minpt_bin);
void drawText(const char *text, float xp, float yp);
void drawLegend(int var);
//-------------------------------------------------


void PlotTRK_EFFnFR(int sample=0, int var=1, bool save=false, TString trkcoll="hitrkEffAnalyzer"){


  // var = 1, 2, 3

   bool lowpt=false;
   int logx=0;

   lowpt_eta = true;

   float pt_min, pt_max;
   pt_min = 0.0;
   pt_max = 100;

   //if(!lowpt) pt_max = 1.5;

   float ymin, ymax;
   ymin = -0.03, ymax = 0.85;
   if(var==3) ymin = -0.03, ymax = 0.92;


   TString varNum;

   varNum.Append(Form("%d",var));

   //TString dir = "./rootOutput_trk/may112011_invstg";  //may092011_test
   TString dir = "./rootOutput_trk/may122011";

   TString outf;

   TString stag;
   TString otag;

   if(sample==0){
     stag = "off393";
     otag = "MinBias";
     //pt_max = 25;
     pt_max = 50;
   }else if(sample==1){
     stag = "dijet30";
     otag = "QCD Pt30";
     pt_max = 50;
   }else if(sample==2){
     stag = "dijet50";
     otag = "QCD Pt50";
     pt_max = 60;
   }else if(sample==3){
     stag = "dijet80";
     otag = "QCD Pt80";
     pt_max = 100;
   }else if(sample==4){
     stag = "dijet110";
     otag = "QCD Pt110";
   }else if(sample==5){
     stag = "dijet170";
     otag = "QCD Pt170";
   }

   outf.Append(Form("trk_eff_fake_sample_%s_%s_var%d",stag.Data(),trkcoll.Data(),var));
   if(var==3) outf.Append(Form("_lowptcut%d",lowpt_eta));

   cout<<"Output file name = "<<outf.Data()<<endl;
     
   TString f_s;

   f_s = "TRKEffnFake_trkhist_may052011_rereco_chi2_v0_hydjetBass_off393_hitrkEffAnalyzer_eta_0.0to1.0_jet_0.0to2000.0_minpT0.0_reb1_var1__CBIN_";

   TString cb1, cb2, cb3, cb4, cb5;
   cb1 = "cbin0to1";
   cb2 = "cbin2to3";
   cb3 = "cbin4to11";
   cb4 = "cbin12to19";
   cb5 = "cbin20to35";
   //cb5 = "cbin28to35";


   // file loading -----------------------------------------------
   TString f_s1;

   f_s1.Append(Form("%s",f_s.Data()));
   f_s1.ReplaceAll("_STAG_",stag.Data());
   f_s1.ReplaceAll("_TRKCOLL_",trkcoll.Data());
   f_s1.ReplaceAll("_CBIN_",cb5.Data());
   f_s1.ReplaceAll("_VAR_",varNum.Data());


   trkf_v1 = (TFile*) loadFile(trkf_v1,Form("%s/%s.root",dir.Data(),f_s1.Data()));

   gEff_v1 = (TGraphAsymmErrors*) trkf_v1->Get("trkEffNFake/gEfficiency");
   gFake_v1 = (TGraphAsymmErrors*) trkf_v1->Get("trkEffNFake/gFakerate");


   TString f_s2;

   f_s2.Append(Form("%s",f_s.Data()));
   f_s2.ReplaceAll("_STAG_",stag.Data());
   f_s2.ReplaceAll("_TRKCOLL_",trkcoll.Data());
   f_s2.ReplaceAll("_CBIN_",cb4.Data());
   f_s2.ReplaceAll("_VAR_",varNum.Data());


   trkf_v2 = (TFile*) loadFile(trkf_v2,Form("%s/%s.root",dir.Data(),f_s2.Data()));
   
   gEff_v2 = (TGraphAsymmErrors*) trkf_v2->Get("trkEffNFake/gEfficiency");
   gFake_v2 = (TGraphAsymmErrors*) trkf_v2->Get("trkEffNFake/gFakerate");


   TString f_s3;

   f_s3.Append(Form("%s",f_s.Data()));
   f_s3.ReplaceAll("_STAG_",stag.Data());
   f_s3.ReplaceAll("_TRKCOLL_",trkcoll.Data());
   f_s3.ReplaceAll("_CBIN_",cb3.Data());
   f_s3.ReplaceAll("_VAR_",varNum.Data());


   trkf_v3 = (TFile*) loadFile(trkf_v3,Form("%s/%s.root",dir.Data(),f_s3.Data()));

   gEff_v3 = (TGraphAsymmErrors*) trkf_v3->Get("trkEffNFake/gEfficiency");
   gFake_v3 = (TGraphAsymmErrors*) trkf_v3->Get("trkEffNFake/gFakerate");


   TString f_s4;

   f_s4.Append(Form("%s",f_s.Data()));
   f_s4.ReplaceAll("_STAG_",stag.Data());
   f_s4.ReplaceAll("_TRKCOLL_",trkcoll.Data());
   f_s4.ReplaceAll("_CBIN_",cb2.Data());
   f_s4.ReplaceAll("_VAR_",varNum.Data());


   trkf_v4 = (TFile*) loadFile(trkf_v4,Form("%s/%s.root",dir.Data(),f_s4.Data()));

   gEff_v4 = (TGraphAsymmErrors*) trkf_v4->Get("trkEffNFake/gEfficiency");
   gFake_v4 = (TGraphAsymmErrors*) trkf_v4->Get("trkEffNFake/gFakerate");

   TString f_s5;

   f_s5.Append(Form("%s",f_s.Data()));
   f_s5.ReplaceAll("_STAG_",stag.Data());
   f_s5.ReplaceAll("_TRKCOLL_",trkcoll.Data());
   f_s5.ReplaceAll("_CBIN_",cb1.Data());
   f_s5.ReplaceAll("_VAR_",varNum.Data());

   trkf_v5 = (TFile*) loadFile(trkf_v5,Form("%s/%s.root",dir.Data(),f_s5.Data()));

   gEff_v5 = (TGraphAsymmErrors*) trkf_v5->Get("trkEffNFake/gEfficiency");
   gFake_v5 = (TGraphAsymmErrors*) trkf_v5->Get("trkEffNFake/gFakerate");


   if(gEff_v1==0||gFake_v1==0||
      gEff_v2==0||gFake_v2==0||
      gEff_v3==0||gFake_v3==0||
      gEff_v4==0||gFake_v4==0||
      gEff_v5==0||gFake_v5==0) {
      Error("[TGraph Loading]","could not get histogram");
      return;
   }


   // Basic canvas and dummy histogram
   TCanvas *call = new TCanvas("call","call",510,510);
   call->cd();
   //call->SetLogx();
   
   if(var==1){
     dum1 = GetDummyHist(pt_min,pt_max,ymin,ymax,"p_{T} [GeV/c]","A #times #epsilon_{TR}",false);
     dum1->Draw("");
      //dum1->GetXaxis()->SetNdivisions(599);
   }else if(var==2){
     dum1 = GetDummyHist(0,1.0,-0.02,0.85,"Z variable","A #times #epsilon_{TR}",false);
     dum1->Draw("");
   }else if(var==3){
     dum1 = GetDummyHist(-2.5,2.5,ymin,ymax,"#eta","A #times #epsilon_{TR}",false);
     dum1->Draw("");
   }
   
   th1Style1(gEff_v1,15,20,1.0,15,1.5,1,1);
   th1Style1(gFake_v1,15,24,1.0,15,1.5,1,1);
   

   th1Style1(gEff_v2,16,20,1.0,16,1.5,1,1);
   th1Style1(gFake_v2,16,24,1.0,16,1.5,1,1);

   th1Style1(gEff_v3,17,20,1.0,17,1.5,1,1);
   th1Style1(gFake_v3,17,24,1.0,17,1.5,1,1);

   th1Style1(gEff_v4,18,20,1.0,18,1.5,1,1);
   th1Style1(gFake_v4,18,24,1.0,18,1.5,1,1);

   th1Style1(gEff_v5,19,20,1.0,19,1.5,1,1);
   th1Style1(gFake_v5,19,24,1.0,19,1.5,1,1);


   drawLegend(var);

   drawText(otag.Data(),0.79,0.88);


   // save 
   if(save){
     cout<<"Output file = "<<outf.Data()<<endl;
     printCanvases(call,outf.Data(),logx,1);
   }


}


void drawLegend(int var){

  if(var==1){
    legd = new TLegend(0.22,0.29,0.52,0.56); 
  }else if(var==3){
    legd = new TLegend(0.20,0.65,0.50,0.92);
  }

  legd->SetMargin(0.3);
  legd->SetBorderSize(0);
  legd->SetFillColor(kWhite);
  legd->SetFillStyle(0);
  legd->SetTextSize(0.03);
  legd->SetMargin(0.3);
  
  if(var==1){
   legd->SetHeader("#sqrt{s_{NN}}=2.76 TeV HYDJET MC |#eta|<1.0");
   legd->AddEntry(gEff_v1,"50~100%","pl");
   legd->AddEntry(gEff_v2,"30~50%","pl");
   legd->AddEntry(gEff_v3,"10~30%","pl");
   legd->AddEntry(gEff_v4,"5~10%","pl");
   legd->AddEntry(gEff_v5,"0~5%","pl");
  }else if(var==3){
    if(lowpt_eta) legd->SetHeader("#sqrt{s_{NN}}=2.76 TeV HYDJET MC p_{T}>0.9 GeV/c");
    else legd->SetHeader("#sqrt{s_{NN}}=2.76 TeV HYDJET MC p_{T}>2.0 GeV/c"); 
    legd->AddEntry(gEff_v1,"50~100%","pl");
    legd->AddEntry(gEff_v2,"30~50%","pl");
    legd->AddEntry(gEff_v3,"10~30%","pl");
    legd->AddEntry(gEff_v4,"5~10%","pl");
    legd->AddEntry(gEff_v5,"0~5%","pl");
  }
  legd->Draw();

}

void drawText(const char *text, float xp, float yp){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(16);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
