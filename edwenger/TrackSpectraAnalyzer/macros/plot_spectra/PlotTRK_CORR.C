//------------------------------------------------------------------------------
// July 05 2010    Andre S. Yoon
// compilable plotting macro for validation plot for type 1 correction (i.e. trk)
//------------------------------------------------------------------------------ 

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "CorrectType.h"

#include "/home/sungho/plots/common/utilities.h"
#include "/home/sungho/plots/common/commonStyle.h"

#include "TLine.h"
#include "TLegend.h"
#endif


TFile *specTruthF=0;
TFile *specMC_reco=0;

TH1D *hdndpt_truth=0; 
TH1D *hdndptMC_recoRaw=0;
TH1D *hdndptMC_recoLev1=0;
TH1D *hdndptMC_recoLev2=0;
TH1D *hdndptMC_recoFull=0;
TH1D *hdndptDATA_reco=0;

TH1D *hdndpt_truth_dum=0;
TH1D *hdndptMC_recoRaw_dum=0;
TH1D *hdndptMC_recoLev1_dum=0;
TH1D *hdndptMC_recoLev2_dum=0;
TH1D *hdndptMC_recoFull_dum=0;
TH1D *hdndptDATA_reco_dum=0;

TPad *pp1=0, *pp1_1=0;
TH1D *dum1=0, *dum1_1=0;

TLegend *legd=0;

//------------------------------------------
void setPad(TCanvas *can);
void setLowerPad(TH1D *dum_low);
TH1D *rerangedHist(TH1D *hRange, int minpt_bin);
void drawLegend();
//------------------------------------------

void PlotTRK_CORR(){

   bool lowpT = true;
   
   float pt_min, pt_max;
   float ymin, ymax;
   float ymin_r, ymax_r;

   if(lowpT){
      pt_min = 0.5, pt_max = 12;
      ymin = 6E-4, ymax = 3E+2;
      ymin_r = 0.97, ymax_r = 1.03;
   }else{
      pt_min = 0.5, pt_max = 60;
      ymin = 6E-8, ymax = 3E+2;
      ymin_r = 0.7, ymax_r = 1.3;
   }
   //------------------ file loading

   // GEN truth
   specTruthF = (TFile*) loadFile(specTruthF,
				  "./rootOutput/CORR_TrkHistMC_july04_qcdMBD6T0704vbinV2_eta_0.0to2.4_jet_0.0to2000.0_trackAna_GEN1_varBin1.root");
   
   hdndpt_truth = (TH1D*) specTruthF->Get("corrTypeOne/hdndpt_raw");
   if(hdndpt_truth==0) {
      Error("[HistLoading]","could not get histogram");
      return;
   }
   
   // MC reco
   specMC_reco = (TFile*) loadFile(specMC_reco,
				   "./rootOutput/CORR_TrkHistMC_july04_qcdMBD6T0704vbinV2_eta_0.0to2.4_jet_0.0to2000.0_trackAna_GEN0_varBin1.root");

   hdndptMC_recoRaw = (TH1D*) specMC_reco->Get("corrTypeOne/hdndpt_raw");
   hdndptMC_recoLev1 = (TH1D*) specMC_reco->Get("corrTypeOne/hdndpt_lev1");
   hdndptMC_recoLev2 = (TH1D*) specMC_reco->Get("corrTypeOne/hdndpt_lev2");
   hdndptMC_recoFull = (TH1D*) specMC_reco->Get("corrTypeOne/hdndpt_full");

   if(hdndptMC_recoFull==0||
      hdndptMC_recoLev1==0||
      hdndptMC_recoLev2==0||
      hdndptMC_recoRaw==0) {
      Error("[HistLoading]","could not get histogram");
      return;
   }


   // Normalization
   float N_evt = 6873768;
   hdndpt_truth->Scale(1./N_evt);
   hdndptMC_recoRaw->Scale(1./N_evt);
   hdndptMC_recoLev1->Scale(1./N_evt);
   hdndptMC_recoLev2->Scale(1./N_evt);
   hdndptMC_recoFull->Scale(1./N_evt);


   // Rebin and etc
   hdndpt_truth->Rebin(4);
   hdndptMC_recoRaw->Rebin(4);
   hdndptMC_recoLev1->Rebin(4);
   hdndptMC_recoLev2->Rebin(4);
   hdndptMC_recoFull->Rebin(4);

   //hdndpt_truth = (TH1D*) rerangedHist(hdndpt_truth,4);
   hdndptMC_recoRaw = (TH1D*) rerangedHist(hdndptMC_recoRaw,3);
   hdndptMC_recoLev1 = (TH1D*) rerangedHist(hdndptMC_recoLev1,3);
   hdndptMC_recoLev2 = (TH1D*) rerangedHist(hdndptMC_recoLev2,3);
   hdndptMC_recoFull = (TH1D*) rerangedHist(hdndptMC_recoFull,3);
   

   // Basic canvas and dummy histogram 
   TCanvas *call = new TCanvas("call","call",510,670);       
   setPad(call);
   
   // ----- pad 1
   pp1->cd();
   pp1->SetLogy();

   dum1 = GetDummyHist(0,pt_max,ymin,ymax,"p_{T}","(1/N_{evt}) dN/dp_{T}",false);
   dum1->Draw("");

   th1Style1(hdndpt_truth,12,24,1.0,12,1.5,1,3);
   th1Style1(hdndptMC_recoRaw,4,25,1.0,4,1.5,1,1);
   th1Style1(hdndptMC_recoLev1,6,26,1.0,6,1.5,1,1);
   th1Style1(hdndptMC_recoLev2,2,24,1.0,2,1.5,1,1);
   th1Style1(hdndptMC_recoFull,1,28,1.0,1,1.5,1,1);

   drawLegend();

   // ----- pad 2
   pp1_1->cd();

   dum1_1 = GetDummyHist(0,pt_max,ymin_r,ymax_r,"p_{T}","RECO/GEN",false);
   setLowerPad(dum1_1);
   dum1_1->Draw("");

   TLine *line = new TLine(0,1,pt_max,1);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();

   
   hdndpt_truth_dum = (TH1D*) hdndpt_truth->Clone("hdndpt_truth_dum");
   hdndptMC_recoRaw_dum = (TH1D*) hdndptMC_recoRaw->Clone("hdndptMC_recoRaw_dum");
   hdndptMC_recoLev1_dum = (TH1D*) hdndptMC_recoLev1->Clone("hdndptMC_recoLev1_dum");
   hdndptMC_recoLev2_dum = (TH1D*) hdndptMC_recoLev2->Clone("hdndptMC_recoLev2_dum");
   hdndptMC_recoFull_dum = (TH1D*) hdndptMC_recoFull->Clone("hdndptMC_recoFull_dum");
   
   hdndpt_truth_dum = (TH1D*) rerangedHist(hdndpt_truth_dum,3); 
   hdndptMC_recoRaw_dum = (TH1D*) rerangedHist(hdndptMC_recoRaw_dum,3);
   hdndptMC_recoLev1_dum = (TH1D*) rerangedHist(hdndptMC_recoLev1_dum,3);
   hdndptMC_recoLev2_dum = (TH1D*) rerangedHist(hdndptMC_recoLev2_dum,3);
   hdndptMC_recoFull_dum = (TH1D*) rerangedHist(hdndptMC_recoFull_dum,3);

   hdndptMC_recoFull_dum->Divide(hdndpt_truth_dum);
   hdndptMC_recoRaw_dum->Divide(hdndpt_truth_dum);
   hdndptMC_recoLev1_dum->Divide(hdndpt_truth_dum);
   hdndptMC_recoLev2_dum->Divide(hdndpt_truth_dum);


   th1Style1(hdndptMC_recoFull_dum,1,28,1.0,1,1.5,1,1);
   th1Style1(hdndptMC_recoRaw_dum,4,25,1.0,4,1.5,1,1);
   th1Style1(hdndptMC_recoLev1_dum,6,26,1.0,6,1.5,1,1);
   th1Style1(hdndptMC_recoLev2_dum,2,24,1.0,2,1.5,1,1);


}

void drawLegend(){

   legd = new TLegend(0.37,0.52,0.76,0.88);    
   legd->SetMargin(0.3);
   legd->SetBorderSize(0);
   legd->SetFillColor(kWhite);
   legd->SetFillStyle(0);
   legd->SetTextSize(0.039);
   legd->SetMargin(0.3);

   legd->SetHeader("pp #rightarrow h+X, 7 TeV MC, |#eta|<2.4");

   legd->AddEntry(hdndpt_truth,"PYTHIA D6T (selected)","l"); 
   legd->AddEntry(hdndptMC_recoRaw,"RECO TRK Eff/FR Raw.","pl");
   legd->AddEntry(hdndptMC_recoLev1,"RECO TRK Eff/FR Corr.","pl");
   legd->AddEntry(hdndptMC_recoLev2,"RECO TRK Eff/FR/Sec Corr.","pl");
   legd->AddEntry(hdndptMC_recoFull,"RECO TRK Eff/FR/Mlt/Sec Corr.","pl");

   legd->Draw();
}

TH1D *rerangedHist(TH1D *hRange, int minpt_bin){
   hRange->GetXaxis()->SetRange(minpt_bin,hRange->GetXaxis()->GetLast());
   return hRange;
}

void setLowerPad(TH1D *dum_low){

   //Just to make lower pannel lable appearing as similar size as upper one
   dum_low->GetYaxis()->SetLabelSize(0.07);
   dum_low->GetYaxis()->SetTitleSize(0.08);
   dum_low->GetYaxis()->SetTitleOffset(1.0);

   dum_low->GetXaxis()->SetLabelSize(0.09);
   dum_low->GetXaxis()->SetLabelOffset(0.05);
   dum_low->GetXaxis()->SetTitleSize(0.11);
   dum_low->GetXaxis()->SetTitleOffset(1.47);

   dum_low->GetYaxis()->SetNdivisions(405);
   dum_low->GetYaxis()-> SetDecimals(true);
   dum_low->Draw();

}

void setPad(TCanvas *can){

   can->cd();

   pp1 = new TPad("p1","p1",0,0.34,1,1,0,0,0);
   pp1->SetBottomMargin(0.0);
   pp1->SetTopMargin(0.05*(1/0.72));

   pp1->Draw();
   pp1->cd();
   pp1->SetNumber(1);

   can->cd();

   pp1_1 = new TPad("p1_1","p1_1",0,0.0,1,0.34,0,0,0);
   pp1_1->SetTopMargin(0);
   pp1_1->SetBottomMargin(0.14*(1/0.34));
   pp1_1->Draw();
   pp1_1->cd();
   pp1_1->SetNumber(2);
}
