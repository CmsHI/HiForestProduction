//------------------------------------------------------------------- 
// June 03 2010, Andre S. Yoon 
// -plot  the returned histograms/tgraph
//--------------------------------------------------------------------

#include "/home/sungho/plots/common/commonStyle.h"

void PlotAll_EVT()
{

   gROOT->LoadMacro("invar_yield_ana_v9.C");

   char dir[100];
   char dir_corr[100];

   bool GEN = false;
   bool CORRECT = false;
   bool REBIN = true;
   bool MC = false;
   bool multcrct =false;
   bool seccrct = false;
   bool cross = false;
   bool loose = false;
   bool jet_based_correction = true;
   bool entire_range = false;
   bool zerobin = false;
   bool onetwothreebin = false;
   bool evteffcorr = false;
   bool full_eta =false;
   bool oneoverpt=true;
   bool minpt=true;

   double eta_max;
   float pt_max;
   float ymin, ymax;
   float ymin_r, ymax_r;

   if(entire_range){
      if(MC) pt_max = 55,  ymin = 3e-13, ymax = 9e1, ymin_r = 0.30, ymax_r = 1.65;
      else pt_max = 55, ymin = 3e-12, ymax = 9e1, ymin_r = 0.7, ymax_r = 1.3;
   }else{
      if(MC) pt_max = 10, ymin = 5e-7, ymax = 2e1, ymin_r = 0.7, ymax_r = 1.4;       
      //else pt_max = 10, ymin = 5e-7, ymax = 5e1, ymin_r = 0.97, ymax_r = 1.03;
      else pt_max = 10, ymin = 5e-7, ymax = 5e1, ymin_r = 0.5, ymax_r = 1.5;        
   }

   if(full_eta) eta_max = 2.4;
   else eta_max = 1.0;


   char file1[100], file2[100], file3[100], file4[100], file5[100], file6[100], file7[100];

   sprintf(file1,"../root_files/TrkHistMC_june04_10Mv2.root");
   sprintf(file2,"../root_files/TrkHistMC_june04_10Mv2.root");
   sprintf(file3,"../root_files/TrkHistMC_june04_10Mv2.root");
   sprintf(file4,"../root_files/TrkHistMC_june04_10Mv2.root");
   sprintf(file5,"../root_files/TrkHistMC_june04_10Mv2.root");
   sprintf(file6,"../root_files/TrkHistMC_june04_10Mv2.root");
   sprintf(file7,"../root_files/TrkHistMC_june04_10Mv2.root");

   //sprintf(file1,"../root_files/TrkHistMC_june04_10McMultv2.root");            
   //sprintf(file2,"../root_files/TrkHistMC_june04_10McMultv2.root"); 
   /*
   sprintf(file3,"../root_files/TrkHistMCv10_QCD_Pt15v3.root");
   sprintf(file4,"../root_files/TrkHistMCv11_qcdPt30.root");
   sprintf(file5,"../root_files/TrkHistMCv11_qcdPt80.root");
   sprintf(file6,"../root_files/TrkHistMCv11_qcdPt170v2.root");
   sprintf(file7,"../root_files/TrkHistMCv11_qcdPt300v2.root");
   */

   MC = true, GEN = true, CORRECT = false;
   sprintf(dir,"preTrackAna");
   sprintf(dir_corr,"trkEffAnalyzer");

   invar_yield_ana_v9_data spec_gen_tight =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
			       dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
			       evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,0,eta_max);

   MC = true, GEN = true, CORRECT = false;
   sprintf(dir,"trackAna");
   sprintf(dir_corr,"trkEffAnalyzer");

   invar_yield_ana_v9_data spec_gen_tight_v2 =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,0,eta_max);


   MC = true, GEN = true, CORRECT = false;
   sprintf(dir,"looseTrackAna");
   sprintf(dir_corr,"trkEffAnalyzer");

   invar_yield_ana_v9_data spec_gen_tight_v3 =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,0,eta_max);



   MC = true, GEN = false, CORRECT = true, multcrct = true, seccrct = true, evteffcorr = true, zerobin = true, onetwothreebin = false;
   sprintf(dir,"looseTrackAna_STD");
   sprintf(dir_corr,"loosetrkEffAnalyzer");

   invar_yield_ana_v9_data spec_rec_loose =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,0,eta_max);


   MC = true, GEN = false, CORRECT = true, multcrct = true, seccrct = true, evteffcorr = true, zerobin = true, onetwothreebin = true;
   sprintf(dir,"trackAna_STD");
   sprintf(dir_corr,"trkEffAnalyzer");

   invar_yield_ana_v9_data spec_rec_tight =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
			       dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
			       evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,0,eta_max);


   sprintf(file1,"../root_files/TrkHistMC_june04_10McMultv2.root");
   sprintf(file2,"../root_files/TrkHistMC_june04_10McMultv2.root");
   sprintf(file3,"../root_files/TrkHistMC_june04_10McMultv2.root");
   sprintf(file4,"../root_files/TrkHistMC_june04_10McMultv2.root");
   sprintf(file5,"../root_files/TrkHistMC_june04_10McMultv2.root");
   sprintf(file6,"../root_files/TrkHistMC_june04_10McMultv2.root");
   sprintf(file7,"../root_files/TrkHistMC_june04_10McMultv2.root");

   MC = true, GEN = false, CORRECT = true, multcrct = true, seccrct = true, evteffcorr = true, zerobin = true, onetwothreebin = true;
   sprintf(dir,"looseTrackAna_STD");
   sprintf(dir_corr,"loosetrkEffAnalyzer");

   invar_yield_ana_v9_data spec_rec_loose_zott =
      invar_yield_ana_v9_graph(file1,file2,file3,file4,file5,file6,file7,
                               dir,dir_corr,GEN,CORRECT,multcrct,seccrct,MC,jet_based_correction,
                               evteffcorr,zerobin,onetwothreebin,cross,oneoverpt,0,eta_max);




   //--------------------------- CMS measurement ---------------------------
   /*
   TF1 *f2 = new TF1("f2","[0]*(1+(sqrt(0.1396**2+x**2)-0.1396)/([1]*[2]))**(-[2])",0,6);
   f2->SetParameters(2*TMath::Pi(),0.13,7.7);
   cms_7000GeV->Fit(f2,"RN");
   f2->SetLineColor(kRed);
   */

   //--------------------------------------------------- PYTHIA truth GEN------------------------
   /*
   gROOT->LoadMacro("pythia_xsection_extractor_scale.C");

   Bool_t REBIN = true;
   Bool_t CDFBIN = true;
   Bool_t ERROR = false;
   Double_t SCALE = 0.04;
   Int_t IREBIN = 2;

   pythia_xsection_extractor_data hmb7000 =
      pythia_xsection_extractor_graph("../root_files/SPECTRAv7_7TeV_Input_7TeV_QCD_IC7_pCutOff10.00.root",
                                      "dsigdpt",REBIN,ERROR,SCALE);

   TGraphErrors *pythia_7000GeV_hqcdR = hmb7000.RInvX;
   TH1D *pythia_7000GeV_hhqcdR = hmb7000.hRInvX;
   */
   //--------------------------------------------------- PYTHIA truth GEN------------------------ 


   if(REBIN){
      TH1D *dndpt_gen_tight = spec_gen_tight.hRInvX;
      TH1D *dndpt_gen_tight_v2 = spec_gen_tight_v2.hRInvX;
      TH1D *dndpt_gen_tight_v3 = spec_gen_tight_v3.hRInvX;
      TH1D *dndpt_rec_loose = spec_rec_loose.hRInvX;
      TH1D *dndpt_rec_tight = spec_rec_tight.hRInvX;
      TH1D *dndpt_rec_loose_zott = spec_rec_loose_zott.hRInvX;
   }else{
      TH1D *dndpt_gen_tight = spec_gen_tight.hInvX;
      TH1D *dndpt_rec_loose = spec_rec_loose.hInvX;
      TH1D *dndpt_rec_tight = spec_rec_tight.hInvX;
   }
   
   
   // Basic canvas and dummy histogram 
   TCanvas *call = new TCanvas("call","call",510,670);       
   call->cd();
   
   // ------------------------------- Pad settings 
   pp1 = new TPad("p1","p1",0,0.34,1,1,0,0,0);
   pp1->SetBottomMargin(0.0);
   pp1->SetTopMargin(0.05*(1/0.72));

   pp1->Draw();
   pp1->cd();
   pp1->SetNumber(1);

   call->cd();

   pp1_1 = new TPad("p1_1","p1_1",0,0.0,1,0.34,0,0,0);
   pp1_1->SetTopMargin(0);
   pp1_1->SetBottomMargin(0.14*(1/0.34));
   pp1_1->Draw();
   pp1_1->cd();
   pp1_1->SetNumber(2);


   // ----- pad 1 
   pp1->cd();
   pp1->SetLogy();


   Char_t xTitle[100],yTitle[100];


   if(!cross) sprintf(yTitle,"Ed^{3}N/d^{3}p [GeV^{-2}c^{3}]");
   else sprintf(yTitle,"Ed^{3}#sigma/d^{3}p [mb GeV^{-2}c^{3}]");
   sprintf(xTitle,"p_{T} [GeV/c]");
   TH1D *dum = GetDummyHist(pt_max,ymin,ymax,xTitle,yTitle);

   dum->Draw();
   dum->GetXaxis()->SetNdivisions(908);
   dum->GetYaxis()->SetNdivisions(407);

   TPad *call_pd = (TPad*)gPad;
   call_pd->cd();
   gPad->SetLogy();

   if(minpt){
      dndpt_rec_tight->GetXaxis()->SetRange(3,dndpt_rec_tight->GetXaxis()->GetLast());
      dndpt_rec_loose->GetXaxis()->SetRange(3,dndpt_rec_loose->GetXaxis()->GetLast());
      dndpt_rec_loose_zott->GetXaxis()->SetRange(3,dndpt_rec_loose_zott->GetXaxis()->GetLast());
   }


   //plotting 
   th1Style1(dndpt_gen_tight,12,24,1.0,12,1.5,1,3);
   th1Style1(dndpt_gen_tight_v2,12,24,1.0,12,1.5,9,3);
   th1Style1(dndpt_gen_tight_v3,12,24,1.0,12,1.5,2,3);
   th1Style1(dndpt_rec_loose,2,24,1.0,2,1.5,1,1);
   th1Style1(dndpt_rec_tight,4,25,1.0,4,1.5,1,1);
   th1Style1(dndpt_rec_loose_zott,6,26,1.0,6,1.5,1,1);


   TLegend *leg2 = new TLegend(0.37,0.47,0.74,0.85);    
   leg2->SetMargin(0.3);
   leg2->SetBorderSize(0);
   leg2->SetFillColor(kWhite);
   leg2->SetFillStyle(0);
   leg2->SetTextSize(0.039);
   leg2->SetMargin(0.3);
   

   if(full_eta) leg2->SetHeader("pp #rightarrow h+X, 7 TeV MC, |#eta|<2.4");
   else leg2->SetHeader("pp #rightarrow h+X, 7 TeV MC, |#eta|<1.0");               
   leg2->AddEntry(dndpt_gen_tight,"PYTHIA D6T (NSD)","l"); 
   leg2->AddEntry(dndpt_gen_tight_v2,"PYTHIA D6T (select)","l");
   leg2->AddEntry(dndpt_gen_tight_v3,"PYTHIA D6T (loose)","l");
   leg2->AddEntry(dndpt_rec_loose,"RECO Loose.","pl");
   leg2->AddEntry(dndpt_rec_tight,"RECO Tight.","pl");
   leg2->AddEntry(dndpt_rec_loose_zott,"RECO Loose (0~3 excl).","pl");

   
   leg2->Draw();



   TLatex * tex;
   tex = new TLatex(0.69,0.88,"CMS Preliminary");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   tex->Draw();



   // ------ pad2
   pp1_1->cd();
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


   
   TH1D *dndpt_gen_tight_dum = (TH1D*) dndpt_gen_tight->Clone("dndpt_gen_tight_dum");
   TH1D *dndpt_rec_loose_dum = (TH1D*) dndpt_rec_loose->Clone("dndpt_rec_loose_dum");
   TH1D *dndpt_rec_tight_dum = (TH1D*) dndpt_rec_tight->Clone("dndpt_rec_tight_dum"); 
   TH1D *dndpt_rec_loose_zott_dum = (TH1D*) dndpt_rec_loose_zott->Clone("dndpt_rec_loose_zott_dum");


   dndpt_gen_tight_dum->Sumw2();
   dndpt_rec_loose_dum->Sumw2();
   dndpt_rec_tight_dum->Sumw2();
   dndpt_rec_loose_zott_dum->Sumw2();


   dndpt_rec_loose_dum->Divide(dndpt_gen_tight_dum);
   dndpt_rec_tight_dum->Divide(dndpt_gen_tight_dum);
   dndpt_rec_loose_zott_dum->Divide(dndpt_gen_tight_dum);


   if(minpt){
      dndpt_rec_loose_dum->GetXaxis()->SetRange(3,dndpt_rec_loose_dum->GetXaxis()->GetLast());
      dndpt_rec_tight_dum->GetXaxis()->SetRange(3,dndpt_rec_tight_dum->GetXaxis()->GetLast());
      dndpt_rec_loose_zott_dum->GetXaxis()->SetRange(3,dndpt_rec_loose_zott_dum->GetXaxis()->GetLast());
   }

   th1Style1(dndpt_rec_loose_dum,2,24,1.0,2,1.5,1,1);
   th1Style1(dndpt_rec_tight_dum,4,25,1.0,4,1.5,1,1);
   th1Style1(dndpt_rec_loose_zott_dum,6,26,1.0,6,1.5,1,1);


   //th1Style1(dndpt1_dum_div_cms7000,13,30,1.0,13,1.5,1,1);
   

   if(entire_range){  
      if(full_eta) printCanvases(call,"spectra_EVT_fullEta",0,1);
      else printCanvases(call,"spectra_EVT_PAS",0,1);
   }else{
      if(full_eta) printCanvases(call,"spectra_EVT_narrow_fullEta",0,1);
      else printCanvases(call,"spectra_EVT_PAS_narrow",0,1);
   }

}





TH1D* GetDummyHist(Float_t xmax, Float_t min, Float_t max,Char_t *xttl,Char_t *yttl) {

   TH1D *dum;
   dum = new TH1D("dum","",100,0.0,xmax);     
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

// Better to keep this in a seperate file
TGraphErrors* CMS_7TEV(Double_t scale){

   //7 TeV measurement at CMS                                                                                                                                                                                 
   TGraphErrors *gre = new TGraphErrors(33);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,0.15,9.94559);
   gre->SetPointError(0,0,0.4460492);
   gre->SetPoint(1,0.25,6.4274);
   gre->SetPointError(1,0,0.2867568);
   gre->SetPoint(2,0.35,3.69525);
   gre->SetPointError(2,0,0.1654095);
   gre->SetPoint(3,0.45,2.2442);
   gre->SetPointError(3,0,0.1009501);
   gre->SetPoint(4,0.55,1.40735);
   gre->SetPointError(4,0,0.06369867);
   gre->SetPoint(5,0.65,0.920301);
   gre->SetPointError(5,0,0.04193949);
   gre->SetPoint(6,0.75,0.616821);
   gre->SetPointError(6,0,0.0283265);
   gre->SetPoint(7,0.85,0.426732);
   gre->SetPointError(7,0,0.01975661);
   gre->SetPoint(8,0.95,0.299085);
   gre->SetPointError(8,0,0.01397616);
   gre->SetPoint(9,1.1,0.183225);
   gre->SetPointError(9,0,0.00845129);
   gre->SetPoint(10,1.3,0.0995915);
   gre->SetPointError(10,0,0.004670335);
   gre->SetPoint(11,1.5,0.056565);
   gre->SetPointError(11,0,0.002704398);
   gre->SetPoint(12,1.7,0.0337805);
   gre->SetPointError(12,0,0.001650384);
   gre->SetPoint(13,1.9,0.0208494);
   gre->SetPointError(13,0,0.001043639);
   gre->SetPoint(14,2.1,0.01323271);
   gre->SetPointError(14,0,0.0006804352);
   gre->SetPoint(15,2.3,0.00852133);
   gre->SetPointError(15,0,0.0004520424);
   gre->SetPoint(16,2.5,0.00555772);
   gre->SetPointError(16,0,0.0003051332);
   gre->SetPoint(17,2.7,0.00369416);
   gre->SetPointError(17,0,0.0002105437);
   gre->SetPoint(18,2.9,0.00275361);
   gre->SetPointError(18,0,0.0001617788);
   gre->SetPoint(19,3.1,0.001784498);
   gre->SetPointError(19,0,0.0001103706);
   gre->SetPoint(20,3.3,0.001251193);
   gre->SetPointError(20,0,8.105581e-05);
   gre->SetPoint(21,3.5,0.000932109);
   gre->SetPoint(22,3.7,0.00072809);
   gre->SetPointError(22,0,5.106009e-05);
   gre->SetPoint(23,3.9,0.000483879);
   gre->SetPointError(23,0,3.643427e-05);
   gre->SetPoint(24,4.1,0.000397255);
   gre->SetPointError(24,0,3.101453e-05);
   gre->SetPoint(25,4.3,0.000281682);
   gre->SetPointError(25,0,2.426278e-05);
   gre->SetPoint(26,4.5,0.0001898557);
   gre->SetPointError(26,0,1.777722e-05);
   gre->SetPoint(27,4.7,0.0001404168);
   gre->SetPointError(27,0,1.377145e-05);
   gre->SetPoint(28,4.9,0.0001152827);
   gre->SetPointError(28,0,1.228788e-05);
   gre->SetPoint(29,5.1,8.94332e-05);
   gre->SetPointError(29,0,1.021356e-05);
   gre->SetPoint(30,5.3,6.56232e-05);
   gre->SetPointError(30,0,8.759586e-06);
   gre->SetPoint(31,5.5,5.48406e-05);
   gre->SetPointError(31,0,8.112416e-06);
   gre->SetPoint(32,5.7,4.50414e-05);
   gre->SetPointError(32,0,7.518979e-06);

   /*
   gre->SetPoint(33,5.9,3.60633e-05);
   gre->SetPointError(33,0,5.641449e-06);
   gre->SetPoint(34,6.25,2.018451e-05);
   gre->SetPointError(34,0,3.081411e-06);
   gre->SetPoint(35,6.75,1.966985e-05);
   gre->SetPointError(35,0,2.63818e-06);
   gre->SetPoint(36,7.25,1.572381e-05);
   gre->SetPointError(36,0,2.668063e-06);
   gre->SetPoint(37,7.75,8.84664e-06);
   gre->SetPointError(37,0,2.136684e-06);
   gre->SetPoint(38,8.25,1.264757e-05);
   gre->SetPointError(38,0,3.756436e-06);
   gre->SetPoint(39,8.75,6.31772e-06);
   gre->SetPointError(39,0,2.762383e-06);
   gre->SetPoint(40,9.25,5.99949e-06);
   gre->SetPointError(40,0,3.246975e-06);
   gre->SetPoint(41,9.25,5.99949e-06);
   gre->SetPointError(41,0,3.239439e-06);
   /*
   

   /*
   Double_t x[33], y[33];
   Double_t xerr[33], yerr[33];

   for(Int_t i=0; i<33;i++){
      gre->GetPoint(i,x[i],y[i]);
      gre->SetPoint(i,x[i],scale*y[i]);
      //gre->GetPointError(i,x[i],yerr[i]);
      xerr[i] = 0, yerr[i] = 0;
      gre->SetPointError(i,scale*xerr[i],scale*yerr[i]);
   }
   */
   return gre;
}



TH1D* ratio_hist_to_func(TH1D* num, double par0, double par1, double par2){
   cout<<"[Ratio to fit used]"<<endl;

   TH1D *hRatio = (TH1D*) num->Clone("hRatio");

   TF1 *f3 = new TF1("f3","[0]*(1+(sqrt(0.1396**2+x**2)-0.1396)/([1]*[2]))**(-[2])",0,6.5);
   f3->SetParameters(par0,par1,par2);
   //f3->SetLineColor(2);

   int nbin = num->GetEntries();
   for(int i=0;i<nbin;i++){

      double cms_value = (double) f3->Eval(hRatio->GetBinCenter(i+1));

      if(hRatio->GetBinCenter(i+1)>0.4 && hRatio->GetBinCenter(i+1)<6.0){
	 //double ratio = cms_value/hRatio->GetBinContent(i+1);  
	 //double ratio_err = cms_value/hRatio->GetBinError(i+1);
	 double ratio = hRatio->GetBinContent(i+1)/cms_value;  
	 double ratio_err = hRatio->GetBinError(i+1)/cms_value;         
      }else{
	 double ratio = -999;
	 double ratio_err = 0.0;
      }
      //double ratio = hRatio->GetBinContent(i+1)/cms_value;
      //double ratio_err = hRatio->GetBinError(i+1)/cms_value;
      hRatio->SetBinContent(i+1,ratio);
      hRatio->SetBinError(i+1,ratio_err);
   }

   return hRatio;

}



/*
TH1D* ratio_hist_to_func(TH1D* num){

   TH1D *hRatio = (TH1D*) num->Clone("hRatio");

   //TGraphErrors* cms_7000 = (TGraphErrors*) CMS_7TEV(1.0);                                                                                                                                                  

   TF1 *f2 = new TF1("f2","[0]*(1+(sqrt(0.1396**2+x**2)-0.1396)/([1]*[2]))**(-[2])",0,10);
   //f2->SetParameters(2*TMath::Pi(),0.145,6.6); 
   f2->SetParameters(9.89667e+00,0.145,6.6);
   //cms_7000->Fit(fit1,"RNO"); 
   //f2->SetLineColor(2); 
   //gre->Fit(f2); 
   //f2->Draw("same"); 

   int nbin = num->GetEntries();
   for(int i=0;i<nbin;i++){
      double cms_value = (double) f2->Eval(hRatio->GetBinCenter(i+1));
      double ratio = hRatio->GetBinContent(i+1)/cms_value;
      double ratio_err = hRatio->GetBinError(i+1)/cms_value;
      hRatio->SetBinContent(i+1,ratio);
      hRatio->SetBinError(i+1,ratio_err);
   }

   return hRatio;

}
*/
