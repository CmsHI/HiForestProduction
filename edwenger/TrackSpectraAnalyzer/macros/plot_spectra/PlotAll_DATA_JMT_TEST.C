//------------------------------------------------------------------- 
// June 07 2010, Andre S. Yoon 
// - plot  the returned histograms/tgraph
// - evolved from PlotAll_EVT.C
//--------------------------------------------------------------------

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

#include "TChain.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"

#include "/home/sungho/plots/common/commonStyle.h"
#include "/home/sungho/plots/common/utilities.h"
#include "cms7TeV_lowpt.h"

#include "invar_yield_ana_v9.C"
#include "invar_yield_ana_GEN_v9.C"

using namespace std;

void PlotAll_DATA_JMT_TEST(bool print=false)
{

   gROOT->LoadMacro("invar_yield_ana_v9.C+");
   gROOT->LoadMacro("invar_yield_ana_GEN_v9.C+");


   char dir[100];
   char dir_corr[100];

   bool GEN = false;
   bool CORRECT = false;
   bool REBIN = false;
   bool MC = false;
   bool multcrct =false;
   bool seccrct = false;
   bool cross = false;
   bool loose = false;
   bool jet_based_correction = true;
   bool entire_range = true;

   bool zerobin = false;
   bool onetwothreebin = false;
   bool evteffcorr = false;
   bool full_eta =true;
   bool oneoverpt=true;
   bool minpt=true;
   bool rescrct=false;
   bool rebOnly=false;

   bool log_scale = true;

   int mom_index = 4;
   int mult_bin = 0;


   double jet_min=0;
   double jet_max=2000;
   double eta_max;
   float pt_max;
   float ymin, ymax;
   float ymin_r, ymax_r;

   if(entire_range){
      if(MC) pt_max = 55,  ymin = 3e-13, ymax = 9e1, ymin_r = 0.30, ymax_r = 1.65;
      else pt_max = 170, ymin = 3e-16, ymax = 9e1, ymin_r = 0.35, ymax_r = 1.65;
   }else{
      if(MC) pt_max = 10, ymin = 5e-7, ymax = 2e1, ymin_r = 0.7, ymax_r = 1.4;       
      else pt_max = 6.5, ymin = 5e-6, ymax = 5e1, ymin_r = 0.75, ymax_r = 1.25;  
   }

   if(full_eta) eta_max = 2.4;
   else eta_max = 1.0;


   char file1[200], file2[200], file3[200], file4[200], file5[200], file6[200], file7[200];

   //sprintf(file1,"../root_files/MB-C10-A20RR-TRKANASKIM-MERGED.root");
   //sprintf(file1,"../root_files/MB-C10-M6RR-JMTskim-v0_0615.root");
   //sprintf(file1,"../root_files/MB-C10-PR9-JMTskim-v0_v2_0615.root");
   sprintf(file1,"../root_files/MB-C10-M6RR-MBskim-v0_0615.root");
   sprintf(file2,"../root_files/TrkHistMC_june09_qcdMB.root");
   sprintf(file3,"../root_files/TrkHistMC_june09_qcdMB.root");
   sprintf(file4,"../root_files/TrkHistMC_june09_qcdPt15.root");
   sprintf(file5,"../root_files/TrkHistMC_june09_qcdPt30v3.root");
   sprintf(file6,"../root_files/TrkHistMC_june09_qcdPt80v3.root");
   sprintf(file7,"../root_files/TrkHistMC_june09_qcdPt170v3.root");


   //MC = false, GEN = false, CORRECT = true, multcrct = true, seccrct = true, evteffcorr = true, zerobin = true, onetwothreebin = true;
   MC = false, GEN = false, CORRECT = true, multcrct = true, seccrct = true, evteffcorr = true, zerobin = true, onetwothreebin = true, rescrct= true;

   sprintf(dir,"trackAna_STD");
   sprintf(dir_corr,"trkEffAnalyzer");

   invar_yield_ana_v9_data spec_rec_tight =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
			       dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
			       evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,rescrct,rebOnly,mom_index,mult_bin,jet_min,jet_max,0,eta_max);

   MC = false, GEN = false, CORRECT = true, multcrct = true, seccrct = true, evteffcorr = true, zerobin = true, onetwothreebin = true, rescrct= true;

   sprintf(dir,"trackAna_STD");
   sprintf(dir_corr,"trkEffAnalyzer");

   jet_min = 0;
   jet_max = 45;//-->60

   invar_yield_ana_v9_data spec_rec_tight_res =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,rescrct,rebOnly,mom_index,mult_bin,jet_min,jet_max,0,eta_max);
   

   //sprintf(file1,"../root_files/MB-C10-PR9-JMTskim-v0_v2_0615.root");     
   //sprintf(file1,"../root_files/MB-C10-PR9-M6RR-JMTskim-v0_0615.root");
   //sprintf(file1,"../root_files/JMT-R10A-PR1-JMTskim-v0_0615.root");
   sprintf(file1,"../root_files/MB-C10-PR9-M6RR-JMTskim-JMT-R10A-PR1-JMTskim-v0_0615.root");
   MC = false, GEN = false, CORRECT = true, multcrct = true, seccrct = true, evteffcorr = true, zerobin = true, onetwothreebin = true, rescrct= true;

   sprintf(dir,"trackAna_STD");
   sprintf(dir_corr,"trkEffAnalyzer");

   jet_min = 45;//=-> 60
   jet_max = 2000;

   invar_yield_ana_v9_data spec_rec_tight_resR =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,rescrct,rebOnly,mom_index,mult_bin,jet_min,jet_max,0,eta_max);


   //GEN MC
   sprintf(file1,"../root_files/TrkHistMC_june09_qcdMB.root");
   sprintf(file2,"../root_files/TrkHistMC_june09_qcdMB.root");
   sprintf(file3,"../root_files/TrkHistMC_june09_qcdMB.root");
   sprintf(file4,"../root_files/TrkHistMC_june09_qcdPt15.root");
   sprintf(file5,"../root_files/TrkHistMC_june09_qcdPt30v3.root");
   sprintf(file6,"../root_files/TrkHistMC_june09_qcdPt80v3.root");
   sprintf(file7,"../root_files/TrkHistMC_june09_qcdPt170v3.root");

   MC = true, GEN = true, CORRECT = false, multcrct = false, seccrct = false, evteffcorr = false, zerobin = false, onetwothreebin = false, rescrct= false;

   sprintf(dir,"preTrackAna");
   sprintf(dir_corr,"trkEffAnalyzer");
  
   jet_min = 0;
   jet_max = 2000;

   invar_yield_ana_v9_data spec_rec_tight_MC =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,rescrct,rebOnly,mom_index,mult_bin,jet_min,jet_max,0,eta_max);


   
   // P0
   sprintf(file1,"../root_files/mergePY8.root");
   sprintf(file2,"../root_files/mergePY8.root");
   sprintf(file3,"../root_files/mergePY8.root");
   sprintf(file4,"../root_files/mergePY8.root");
   sprintf(file5,"../root_files/mergePY8.root");
   sprintf(file6,"../root_files/mergePY8.root");
   sprintf(file7,"../root_files/mergePY8.root");

   /*
   sprintf(file1,"../root_files/mergeP0.root");
   sprintf(file2,"../root_files/mergeP0.root");
   sprintf(file3,"../root_files/mergeP0.root");
   sprintf(file4,"../root_files/mergeP0.root");
   sprintf(file5,"../root_files/mergeP0.root");
   sprintf(file6,"../root_files/mergeP0.root");
   sprintf(file7,"../root_files/mergeP0.root");
   */

   MC = true, GEN = true, CORRECT = false, multcrct = false, seccrct = false, evteffcorr = false, zerobin = false, onetwothreebin = false, rescrct= false;

   sprintf(dir,"preTrackAna");
   sprintf(dir_corr,"trkEffAnalyzer");

   jet_min = 0;
   jet_max = 2000;

   invar_yield_ana_v9_data spec_rec_tight_MC_p0 =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,rescrct,rebOnly,mom_index,mult_bin,jet_min,jet_max,0,eta_max);


   // Pythia 8
   sprintf(file1,"../root_files/mergePY8.root");
   sprintf(file2,"../root_files/mergePY8.root");
   sprintf(file3,"../root_files/mergePY8.root");
   sprintf(file4,"../root_files/mergePY8.root");
   sprintf(file5,"../root_files/mergePY8.root");
   sprintf(file6,"../root_files/mergePY8.root");
   sprintf(file7,"../root_files/mergePY8.root");

   MC = true, GEN = true, CORRECT = false, multcrct = false, seccrct = false, evteffcorr = false, zerobin = false, onetwothreebin = false, rescrct= false;

   sprintf(dir,"preTrackAna");
   sprintf(dir_corr,"trkEffAnalyzer");

   jet_min = 0;
   jet_max = 2000;

   invar_yield_ana_v9_data spec_rec_tight_MC_py8 =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,rescrct,rebOnly,mom_index,mult_bin,jet_min,jet_max,0,eta_max);


   //--------------------------- CMS measurement ---------------------------
   TGraphErrors* cms_7000GeV = (TGraphErrors*) CMS_7TEV(1);

   //TF1* f2 = new TF1("fitpythia","[0]*(1+(x/[1])+(pow(x,2)/[2])+(pow(x,3)/[3])+(pow(x,4)/[4])+(pow(x,5)/[5]))^[6]",0.2,6.0);
   //f2->SetParameters(1.81426e+01,1.08477e+00,1,1,1,1,-9.48034e+00);
   //cms_7000GeV->Fit(f2,"RN");

   //TF1 *f2 = new TF1("fitTsallis","[0]*(1+(sqrt(0.1396**2+x**2)-0.1396)/([1]*[2]))**(-[2])",0.5,200);
   //f2->SetParameters(1.32522e+01,1.67841e-01,7.26036e+00); // fit the measurement
   //f2->SetParameters(1.82045e+01,1.53540e-01,7.06942e+00); // fit the fit (used in paper)
   //f2->SetLineColor(1);

   TF1* f2 = new TF1("fitPythia","[0]*(1+(x/[1])+(pow(x,2)/[2])+(pow(x,3)/[3])+(pow(x,4)/[4])+(pow(x,5)/[5]))^[6]",0.5,180);
   f2->SetParameters(2.96361e+01,2.66339e-01,1.34086e-01,2.96428e-01,1.74015e+00,5.16931e+00,-1.70643e+00);


   //--------------------------------------------------- PYTHIA truth GEN------------------------
   //gROOT->LoadMacro("invar_yield_ana_GEN_v9.C+");

   char file1_GEN[100], file2_GEN[100], file3_GEN[100], file4_GEN[100], file5_GEN[100], file6_GEN[100], file7_GEN[100];
   char file8_GEN[100], file9_GEN[100], file10_GEN[100];

   sprintf(file1_GEN,"../root_files/TrkHistMC_june04_qcdPt0To15v3cFix.root");
   sprintf(file2_GEN,"../root_files/TrkHistMC_june04_qcdPt15to20v2.root");
   sprintf(file3_GEN,"../root_files/TrkHistMC_june04_qcdPt20to30v2.root");
   sprintf(file4_GEN,"../root_files/TrkHistMC_june04_qcdPt30to50v2.root");
   sprintf(file5_GEN,"../root_files/TrkHistMC_june04_qcdPt50to80v2.root");
   sprintf(file6_GEN,"../root_files/TrkHistMC_june04_qcdPt80to120v2.root");
   sprintf(file7_GEN,"../root_files/TrkHistMC_june04_qcdPt120to170v2.root");
   sprintf(file8_GEN,"../root_files/TrkHistMC_june04_qcdPt170to230v2.root");
   sprintf(file9_GEN,"../root_files/TrkHistMC_june04_qcdPt230to300v2.root");
   sprintf(file10_GEN,"../root_files/TrkHistMC_june04_qcdPt300to380v2.root");

   invar_yield_ana_GEN_v9_data  gen_nsd_spect  =
      invar_yield_ana_GEN_v9_graph(file1_GEN,file2_GEN,file3_GEN,file4_GEN,
				   file5_GEN,file6_GEN,file7_GEN,file8_GEN,
				   file9_GEN,file10_GEN,true,117.,0,eta_max);

   //--------------------------------------------------- PYTHIA truth GEN------------------------ 

   double intLum;
   TH1D *dndpt_rec_tight=0;
   TH1D *dndpt_rec_tight_res=0;
   TH1D *dndpt_rec_tight_resR=0;

   TH1D *dndpt_rec_tight_added=0;

   TH1D *dndpt_rec_tight_MC=0;
   TH1D *dndpt_rec_tight_MC_p0=0;
   TH1D *dndpt_rec_tight_MC_py8=0;
   TH1D *dndpt_gen_nsd=0;

   TGraphErrors *dndpt_gen_nsd_tg=0;

   if(REBIN){
      dndpt_rec_tight = spec_rec_tight.hRInvX;
      intLum = spec_rec_tight.integratedLum;
      dndpt_rec_tight->Scale(1./intLum);
      dndpt_rec_tight_res = spec_rec_tight_res.hRInvX;
      intLum = spec_rec_tight_res.integratedLum;
      dndpt_rec_tight_res->Scale(1./intLum);
      dndpt_rec_tight_MC = spec_rec_tight_MC.hRInvX;
      intLum = spec_rec_tight_MC.integratedLum;
      dndpt_rec_tight_MC->Scale(1./intLum);
      //dndpt_gen_nsd = gen_nsd_spect.hRInvX;
      //dndpt_gen_nsd_tg = gen_nsd_spect.RInvX;
   }else{
      dndpt_rec_tight = spec_rec_tight.hRInvX;
      intLum = spec_rec_tight.integratedLum;
      dndpt_rec_tight->Scale(1./intLum);

      dndpt_rec_tight_res = spec_rec_tight_res.hRInvX;
      intLum = spec_rec_tight_res.integratedLum;
      dndpt_rec_tight_res->Scale(1./intLum);

      dndpt_rec_tight_resR = spec_rec_tight_resR.hRInvX;
      intLum = spec_rec_tight_resR.integratedLum;
      dndpt_rec_tight_resR->Scale(1./intLum);
      float scale = 1.15223531582625199e-04/5.06700000000000000e+04;

      //dndpt_rec_tight_resR->Scale(scale);
      //dndpt_rec_tight_resR->Scale(1./8.5E1);
      //dndpt_rec_tight_resR->Scale(1./7.5E2);
      dndpt_rec_tight_resR->Scale(1./4E2);     
      
      dndpt_rec_tight_MC = spec_rec_tight_MC.hInvX;
      intLum = spec_rec_tight_MC.integratedLum;
      dndpt_rec_tight_MC->Scale(1./intLum);
      
      dndpt_rec_tight_MC_p0 = spec_rec_tight_MC_p0.hRInvX;
      intLum = spec_rec_tight_MC_p0.integratedLum;
      dndpt_rec_tight_MC_p0->Scale(1./intLum);

      dndpt_rec_tight_MC_py8 = spec_rec_tight_MC_py8.hRInvX;
      intLum = spec_rec_tight_MC_py8.integratedLum;
      dndpt_rec_tight_MC_py8->Scale(1./intLum);
      
      //dndpt_rec_tight_MC_py8->Draw("hist");

      dndpt_gen_nsd = gen_nsd_spect.hRInvX;
      dndpt_gen_nsd_tg = gen_nsd_spect.RInvX;

      dndpt_rec_tight_added = (TH1D*) dndpt_rec_tight_res->Clone("dndpt_rec_tight_added");
      dndpt_rec_tight_added->Add(dndpt_rec_tight_resR);
   }
   
   dndpt_rec_tight_added->Fit(f2,"RN");

   
   // Basic canvas and dummy histogram 
   TCanvas *call = new TCanvas("call","call",510,670);       
   call->cd();
   //if(entire_range) call->SetLogx();

   // ------------------------------- Pad settings 
   TPad *pp1 = new TPad("p1","p1",0,0.34,1,1,0,0,0);
   pp1->SetBottomMargin(0.0);
   pp1->SetTopMargin(0.05*(1/0.72));

   pp1->Draw();
   pp1->cd();
   pp1->SetNumber(1);

   call->cd();

   TPad *pp1_1 = new TPad("p1_1","p1_1",0,0.0,1,0.34,0,0,0);
   pp1_1->SetTopMargin(0);
   pp1_1->SetBottomMargin(0.14*(1/0.34));
   pp1_1->Draw();
   pp1_1->cd();
   pp1_1->SetNumber(2);


   // ----- pad 1 
   pp1->cd();
   pp1->SetLogy();
   if(log_scale && entire_range) pp1->SetLogx();

   Char_t xTitle[100],yTitle[100];


   if(!cross) sprintf(yTitle,"Ed^{3}N/d^{3}p [GeV^{-2}c^{3}]");
   else sprintf(yTitle,"Ed^{3}#sigma/d^{3}p [mb GeV^{-2}c^{3}]");
   sprintf(xTitle,"p_{T} [GeV/c]");
   TH1D *dum = GetDummyHist(pt_max,ymin,ymax,xTitle,yTitle);

   dum->Draw();
   dum->GetXaxis()->SetNdivisions(908);
   dum->GetYaxis()->SetNdivisions(407);

   dndpt_rec_tight_MC_py8->Draw("histsame");

   TPad *call_pd = (TPad*)gPad;
   call_pd->cd();
   gPad->SetLogy();

   if(minpt){
      dndpt_rec_tight->GetXaxis()->SetRange(3,dndpt_rec_tight->GetXaxis()->GetLast());
      dndpt_rec_tight_res->GetXaxis()->SetRange(3,dndpt_rec_tight_res->GetXaxis()->GetLast());
      dndpt_rec_tight_resR->GetXaxis()->SetRange(3,dndpt_rec_tight_resR->GetXaxis()->GetLast());
      dndpt_rec_tight_MC->GetXaxis()->SetRange(3,dndpt_rec_tight_MC->GetXaxis()->GetLast());
      dndpt_rec_tight_MC_p0->GetXaxis()->SetRange(3,dndpt_rec_tight_MC_p0->GetXaxis()->GetLast());
      //dndpt_rec_tight_MC_p0->GetXaxis()->SetRange(3,
						  //dndpt_rec_tight_MC_p0->GetXaxis()->FindBin(19));
      //dndpt_rec_tight_MC_py8->GetXaxis()->SetRange(3,
      //dndpt_rec_tight_MC_py8->GetXaxis()->GetLast());
      dndpt_rec_tight_added->GetXaxis()->SetRange(3,dndpt_rec_tight_added->GetXaxis()->GetLast());
      dndpt_gen_nsd->GetXaxis()->SetRange(3,dndpt_gen_nsd->GetXaxis()->GetLast());
   }


   //plotting 

   f2->SetLineStyle(1);
   f2->SetLineWidth(1.5);
   //f2->Draw("same");
   th1Style1(cms_7000GeV,13,30,1.0,13,1,1,1);
   //th1Style1(dndpt_gen_nsd,96,20,1,96,2,1,2);
   th1Style1(dndpt_gen_nsd_tg,96,20,1.0,96,1.5,1,2); 
   //th1Style1(dndpt_rec_tight_MC_p0,2,20,1.0,2,1.5,2,3);
   //th1Style1(dndpt_rec_tight_MC_py8,4,20,1.0,4,1.5,2,3);
   //th1Style1(dndpt_rec_tight,6,25,1.0,6,1.5,1,1);
   //th1Style1(dndpt_rec_tight_res,2,24,1.0,2,1.5,1,1);
   //th1Style1(dndpt_rec_tight_resR,4,24,1.0,4,1.5,1,1);
   th1Style1(dndpt_rec_tight_added,1,20,1.0,1,1.5,1,1);




   TLegend *leg2=0;
   if(log_scale && entire_range) leg2 = new TLegend(0.21,0.09,0.58,0.47);
   else leg2 = new TLegend(0.47,0.47,0.84,0.85);
   leg2->SetMargin(0.3);
   leg2->SetBorderSize(0);
   leg2->SetFillColor(kWhite);
   leg2->SetFillStyle(0);
   //leg2->SetTextSize(0.039);
   leg2->SetTextSize(0.035);
   leg2->SetMargin(0.3);
   

   if(full_eta) leg2->SetHeader("pp #rightarrow h+X, 7 TeV DATA, |#eta|<2.4");
   else leg2->SetHeader("pp #rightarrow h+X, 7 TeV DATA, |#eta|<1.0");               
   //leg2->AddEntry(dndpt_rec_tight,"MB full","pl");
   //leg2->AddEntry(dndpt_rec_tight_res,"MB (E_{T}<60 GeV/c)","pl");
   //leg2->AddEntry(dndpt_rec_tight_resR,"Jet15U (E_{T}>60 GeV/c)","pl");
   leg2->AddEntry(dndpt_rec_tight_added,"MB + Jet15U ","pl");
   leg2->AddEntry(cms_7000GeV,"CMS 7 TeV |#eta|<2.4","pl");
   leg2->AddEntry(dndpt_gen_nsd_tg,"PYTHIA D6T","l");
   leg2->AddEntry(dndpt_rec_tight_MC_p0,"PYTHIA P0","l");
   leg2->AddEntry(dndpt_rec_tight_MC_py8,"PYTHIA 8","l");
   leg2->AddEntry(f2,"power-law fit","l");
   leg2->Draw();
   


   TLatex * tex;
   tex = new TLatex(0.69,0.88,"CMS Preliminary");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   tex->Draw();



   // ------ pad2
   pp1_1->cd();
   if(log_scale && entire_range) pp1_1->SetLogx();

   sprintf(yTitle,"Ratio");

   //TH1D *dum2 = GetDummyHist(-0.1,2.3,xTitle,yTitle);
   TH1D *dum2 = GetDummyHist(pt_max,ymin_r,ymax_r,xTitle,yTitle);    

   dum2->GetYaxis()->SetLabelSize(0.07);
   dum2->GetYaxis()->SetTitleSize(0.10);
   dum2->GetYaxis()->SetTitleOffset(0.75);

   dum2->GetXaxis()->SetLabelSize(0.09);
   dum2->GetXaxis()->SetLabelOffset(0.05);
   dum2->GetXaxis()->SetTitleSize(0.11);
   dum2->GetXaxis()->SetTitleOffset(1.47);

   dum2->GetYaxis()->SetNdivisions(405);
   dum2->GetYaxis()-> SetDecimals(true);
   dum2->Draw();
   
   TLine *line = new TLine(0,1,pt_max,1);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();

   
   TH1D *dndpt_rec_tight_dum = (TH1D*) dndpt_rec_tight->Clone("dndpt_rec_tight_dum"); 
   TH1D *dndpt_rec_tight_res_dum = (TH1D*) dndpt_rec_tight_res->Clone("dndpt_rec_tight_res_dum");
   TH1D *dndpt_rec_tight_resR_dum = (TH1D*) dndpt_rec_tight_resR->Clone("dndpt_rec_tight_resR_dum");
   TH1D *dndpt_rec_tight_MC_dum = (TH1D*) dndpt_rec_tight_MC->Clone("dndpt_rec_tight_MC_dum");
   TH1D *dndpt_rec_tight_MC_p0_dum = (TH1D*) dndpt_rec_tight_MC_p0->Clone("dndpt_rec_tight_MC_p0_dum");
   
   TH1D *dndpt_rec_tight_added_dum = (TH1D*) dndpt_rec_tight_added->Clone("dndpt_rec_tight_added_dum");
   TH1D *dndpt_rec_tight_added_dum2 = (TH1D*) dndpt_rec_tight_added->Clone("dndpt_rec_tight_added_dum2");



   TH1D *dndpt_gen_nsd_dum = (TH1D*) dndpt_gen_nsd->Clone("dndpt_gen_nsd_dum");

   TH1D *dndpt1_dum_div_cms7000 = (TH1D*) ratio_hist_to_func(dndpt_rec_tight_added_dum,f2,0.4,200.0);
   TH1D *dndpt1_dum_div_cms7000_res = (TH1D*) ratio_hist_to_func(dndpt_rec_tight_res_dum,f2,0.4,6.0);
   TH1D *dndpt1_dum_div_cms7000_resR = (TH1D*) ratio_hist_to_func(dndpt_rec_tight_resR_dum,f2,0.4,6.0);
   TGraphErrors *cms7000_div_fit = (TGraphErrors*) ratio_func_to_func(cms_7000GeV,f2,0.4,6.0);

   /*
   TH1D *dndpt1_dum_div_cms7000 = (TH1D*) ratio_hist_to_func(dndpt_rec_tight_dum,f2,0.4,200.0);
   TH1D *dndpt1_dum_div_cms7000_res = (TH1D*) ratio_hist_to_func(dndpt_rec_tight_res_dum,f2,0.4,200.0);
   TH1D *dndpt1_dum_div_cms7000_resR = (TH1D*) ratio_hist_to_func(dndpt_rec_tight_resR_dum,f2,0.4,200.0);
   TGraphErrors *cms7000_div_fit = (TGraphErrors*) ratio_func_to_func(cms_7000GeV,f2,0.4,200.0);
   */
   

   dndpt_gen_nsd_dum->Sumw2();
   dndpt_rec_tight_dum->Sumw2();
   dndpt_rec_tight_res_dum->Sumw2();
   dndpt_rec_tight_resR_dum->Sumw2();
   dndpt_rec_tight_MC_dum->Sumw2();
   dndpt_rec_tight_MC_p0_dum->Sumw2();
   dndpt_rec_tight_added_dum->Sumw2();
   dndpt_rec_tight_added_dum2->Sumw2();


   if(minpt){
      dndpt_rec_tight_dum->GetXaxis()->SetRange(3,dndpt_rec_tight_dum->GetXaxis()->GetLast());
      dndpt_rec_tight_added_dum2->GetXaxis()->SetRange(3,dndpt_rec_tight_added_dum2->GetXaxis()->GetLast());

      dndpt_rec_tight_added_dum->GetXaxis()->SetRange(3,dndpt_rec_tight_added_dum->GetXaxis()->GetLast());
      dndpt_rec_tight_MC_dum->GetXaxis()->SetRange(3,dndpt_rec_tight_MC_dum->GetXaxis()->GetLast());
      dndpt_rec_tight_MC_p0_dum->GetXaxis()->SetRange(3,dndpt_rec_tight_MC_p0_dum->GetXaxis()->GetLast());

   }

   dndpt_rec_tight_added_dum2->Divide(dndpt_gen_nsd_dum);
   //dndpt_rec_tight_dum->Divide(dndpt_rec_tight_MC_dum);
   //dndpt_rec_tight_added_dum->Divide(dndpt_rec_tight_dum);
   dndpt_rec_tight_added_dum->Divide(dndpt_rec_tight_MC_p0_dum);      
   

   th1Style1(dndpt_rec_tight_added_dum,2,20,1.0,2,1.5,2,3);
   th1Style1(dndpt_rec_tight_added_dum2,96,20,1.0,96,1.5,2,3);

   th1Style1(dndpt1_dum_div_cms7000,1,20,1.0,1,1.5,1,1);
   //th1Style1(dndpt1_dum_div_cms7000_res,2,24,1.0,2,1.5,1,1);
   //th1Style1(dndpt1_dum_div_cms7000_resR,6,25,1.0,6,1.5,1,1);
   //th1Style1(cms7000_div_fit,13,30,1.0,13,1.5,1,1);
   //th1Style1(dndpt_rec_tight_added_dum,1,20,1.0,1,1.5,1,1);


   // Print the values!!
   //PrintXsection(dndpt_rec_tight,117.);
   //PrintXsection(dndpt_rec_tight_added,1.); 

   if(print){
      if(entire_range){  
	 if(log_scale){
	    if(full_eta) printCanvases(call,"spectra_DATA_JMT_TEST_fullEta_logx_v4",2,1);
	    else printCanvases(call,"spectra_DATA_JMT_TEST_logx_v4",2,1);
	 }else{
	    if(full_eta) printCanvases(call,"spectra_DATA_JMT_TEST_fullEta_v4",2,1);
	    else printCanvases(call,"spectra_DATA_JMT_TEST_v4",2,1);
	 }
      }else{
	 if(full_eta) printCanvases(call,"spectra_DATA_JMT_TEST_narrow_fullEta_v4",1,1);
	 else printCanvases(call,"spectra_DATA_JMT_TEST_narrow_v4",1,1);
      }
   }

}



