#include "commonStyle.h"

void plotEvtSelEff_RECO_NORMAL_V3() {

  gStyle->SetOptStat(0);

  //TFile *f = new TFile("../root_files/trkhistsMC_general_5K.root");
  //TFile *f = new TFile("/home/sungho/cmssw/CMSSW_357_SpectraPPV2/src/edwenger/Skims/test/trkhistsMC.root");
  //TFile *f = new TFile("../root_files/TrkHistMCv7_QCD_1Mv2.root");
  //TFile *f = new TFile("../root_files/trkhistsMC_general_50Kv3.root");
  //TFile *f = new TFile("../root_files/trkhistsMC_general_1Kv2.root"); 
  //TFile *f = new TFile("../root_files/trkhistsMC_100K.root");
  //TFile *f = new TFile("../root_files/trkhistsMC_50K.root"); 
  //TFile *f = new TFile("../root_files/trkhistsMC_50Kv2.root");
  //TFile *f = new TFile("../root_files/TrkHistMC_june04_10Mv2.root");
  //TFile *f = new TFile("../root_files/TrkHistMC_june04_QCD_10Mv6.root");
  TFile *f = new TFile("../root_files/TrkHistMC_june09_qcdMB.root");

  bool includeCorr = false;
  bool full = false;


  char chNSD[100];
  if(full) sprintf(chNSD,"preTrgAna/hGenRecMultNSD");
  else sprintf(chNSD,"preTrgAna/hGenRecMultNSD_STD");
  
  char chNSDvtx[100];
  if(full) sprintf(chNSDvtx,"postVtxAna/hRecMult");
  else sprintf(chNSDvtx,"postVtxAna/hRecMult_STD");
  
  char chNSDtrkvtx[100];
  if(full) sprintf(chNSDtrkvtx,"postTrkVtxAna/hRecMult");
  else sprintf(chNSDtrkvtx,"postTrkVtxAna/hRecMult_STD");
  

  TH1F *hNSD_GEN = (TH1F*) f->Get("preTrgAna/hGenMultNSD");
  TH1F *hNSDvtx_CORR = (TH1F*) f->Get("looseTrackAna_STD/hRecMult_STD_corr");
  TH1F *hNSDtrkvtx_CORR = (TH1F*) f->Get("trackAna_STD/hRecMult_STD_corr");


  //------------
  TH1F *hNSD = (TH1F*) f->Get(chNSD);
  TH1F *hSD = (TH1F*) f->Get("preTrgAna/hGenRecMultSD_STD");

  TH1F *hNSDvtx = (TH1F*) f->Get(chNSDvtx);
  TH1F *hNSDtrkvtx = (TH1F*) f->Get(chNSDtrkvtx);
  
  TH1F *hSDvtx = (TH1F*) f->Get("postVtxAna/hGenRecMultSD_STD");
  TH1F *hSDtrkvtx = (TH1F*) f->Get("postTrkVtxAna/hGenRecMultSD_STD");
  
  //  hNSDvtx_CORR->Add(hSDvtx,-1.0);


  //------------
  TH1F *gNSD = (TH1F*) hNSD->Clone("gNSD");
  TH1F *gSD = (TH1F*) hSD->Clone("gSD");
  
  TH1F *gNSDvtx = (TH1F*) hNSDvtx->Clone("gNSDvtx");
  TH1F *gNSDtrkvtx = (TH1F*) hNSDtrkvtx->Clone("gNSDtrkvtx");

  TH1F *gNSDvtx_CORR = (TH1F*) hNSDvtx_CORR->Clone("gNSDvtx_CORR");
  TH1F *gNSDtrkvtx_CORR = (TH1F*) hNSDtrkvtx_CORR->Clone("gNSDtrkvtx_CORR");
  
  TH1F *gSDvtx = (TH1F*) hSDvtx->Clone("gSDvtx");
  TH1F *gSDtrkvtx = (TH1F*) hSDtrkvtx->Clone("gSDtrkvtx");

  TH1F *gSDvtxF = (TH1F*) hSDvtx->Clone("gSDvtx");
  TH1F *gSDtrkvtxF = (TH1F*) hSDtrkvtx->Clone("gSDtrkvtx");


  gSDvtx->Divide(gNSDvtx);
  gSDtrkvtx->Divide(gNSDtrkvtx);

  gNSDvtx->Add(gSDvtxF,-1);
  gNSDtrkvtx->Add(gSDtrkvtxF,-1);
  
  //gNSDvtx_CORR->Add(gSDvtxF,-1);
  //gNSDtrkvtx_CORR->Add(gSDtrkvtxF,-1);

  gNSDvtx->Divide(gNSD);
  gNSDtrkvtx->Divide(gNSD);
  
  gNSDvtx_CORR->Divide(gNSD);
  gNSDtrkvtx_CORR->Divide(gNSD);
  

  cout<<"Total (Truth(reco))" <<hNSD->Integral()<<endl;
  cout<<"Total (corrected)" <<hNSDvtx_CORR->Integral()<<endl;
  

  if(includeCorr){
     TH1F *dum1 = new TH1F("dum1",";Charged particle multiplicity;Fraction of events",100,-2,35);
     TH1F *dum2 = new TH1F("dum2",";Charged particle multiplicity;Selection efficiency",100,0,25);
     dum1->SetMaximum(0.045);
  }else{
     TH1F *dum1 = new TH1F("dum1",";Charged particle multiplicity;Fraction of events",100,-2,100);
     TH1F *dum2 = new TH1F("dum2",";Charged particle multiplicity;Selection efficiency",100,0,55);
  }

  dum1->SetMaximum(0.045);
  dum1->GetYaxis()->SetTitleOffset(1.8);
  dum1->GetYaxis()->CenterTitle();
  dum1->GetXaxis()->CenterTitle();
  dum1->GetYaxis()->SetDecimals();
  


  dum2->SetMaximum(1.2);
  dum2->GetYaxis()->CenterTitle();
  dum2->GetXaxis()->CenterTitle();
  dum2->GetYaxis()->SetDecimals();
  

  cout<<hNSD_GEN->GetName()<<" - Fraction of zero bin = "<<GetZeroBinFraction(hNSD_GEN)<<endl;

  TCanvas *call = new TCanvas("call","call",480,550);
  TCanvas *call2 = new TCanvas("call2","call2",480,550);


  call->cd();
  dum1->Draw();

  hNSD_GEN->Sumw2();
  hNSD_GEN->SetMarkerStyle(3);
  hNSD_GEN->Scale(1./hNSD_GEN->GetEntries());

  int norm = hNSD->GetEntries();
  
  hNSD->Sumw2();
  hNSD->Scale(1./norm);


  hNSDtrkvtx->Sumw2();
  hNSDtrkvtx->Scale(1./norm); 
  //hNSDtrkvtx->Scale(1./hNSDtrkvtx->GetEntries());

  hNSDvtx->Sumw2();
  hNSDvtx->Scale(1./norm);
  //hNSDvtx->Scale(1./hNSDvtx->GetEntries());

  hNSDvtx_CORR->Sumw2();
  hNSDvtx_CORR->Scale(1./norm);

  hNSDtrkvtx_CORR->Sumw2();
  hNSDtrkvtx_CORR->Scale(1./norm);


  //plotting
  th1Style1(hNSD,1,25,1.5,1,1.5,1,1);
  th1Style1(hNSDtrkvtx,1,24,1.5,1,1.5,1,1);
  th1Style1(hNSDvtx,1,20,1.5,1,1.5,1,1);
  if(includeCorr) {
     th1Style1(hNSDvtx_CORR,2,20,1.5,2,2,1,1);
     th1Style1(hNSDtrkvtx_CORR,4,24,1.5,4,2,1,1);
  }

  if(includeCorr) TLegend *leg1 = new TLegend(0.33,0.19,0.71,0.43);
  else TLegend *leg1 = new TLegend(0.49,0.61,0.87,0.83);

  leg1->SetMargin(0.3);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(kWhite);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.038);
  leg1->SetMargin(0.3);
  leg1->SetHeader("  PYTHIA 7 TeV");
  leg1->AddEntry(hNSD,"Before selection","pl");
  leg1->AddEntry(hNSDvtx,"pixel vertex","pl");
  leg1->AddEntry(hNSDtrkvtx,"track vertex","pl");
  if(includeCorr) {
     leg1->AddEntry(hNSDvtx_CORR,"1 pixel vertex (corr.)","pl");
     leg1->AddEntry(hNSDtrkvtx_CORR,"~ 4 track vertex (corr.)","pl");
  }
  leg1->Draw();

  TLatex * tex;
  tex = new TLatex(0.62,0.88,"CMS Preliminary");
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->SetNDC();
  tex->Draw();


  //c1->cd(2);
  call2->cd();
  dum2->Draw();

  if(includeCorr){
     TLine *line = new TLine(0,1,25,1);
     line->SetLineStyle(2);
     line->SetLineWidth(3);
     line->Draw();
  }

  th1Style1(gSDvtx,2,22,1.5,1,2,9,3);
  th1Style1(gSDtrkvtx,4,26,1.5,1,2,2,3);

  th1Style1(gNSDvtx,1,20,1.5,1,2,1,1);
  th1Style1(gNSDtrkvtx,1,24,1.5,1,2,1,1);

  if(includeCorr){
     th1Style1(gNSDvtx_CORR,2,20,1.5,2,2,1,1);
     th1Style1(gNSDtrkvtx_CORR,4,24,1.5,4,2,1,1);
  }


  //TLegend *leg2 = new TLegend(0.50,0.46,0.88,0.65);
  //TLegend *leg2 = new TLegend(0.44,0.26,0.82,0.54);    
  //TLegend *leg2 = new TLegend(0.47,0.26,0.85,0.54);  
  if(includeCorr) TLegend *leg2 = new TLegend(0.42,0.29,0.80,0.63); 
  else TLegend *leg2 = new TLegend(0.42,0.43,0.80,0.71);
  leg2->SetMargin(0.3);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(kWhite);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.038);
  leg2->SetMargin(0.3);
  leg2->SetHeader("  PYTHIA 7 TeV");
  leg2->AddEntry(gNSDvtx,"1 pixel vertex (NSD)","pl");
  leg2->AddEntry(gNSDtrkvtx,"~ 4 track vertex (NSD)","pl");
  leg2->AddEntry(gSDvtx,"1 pixel vertex (SD)","l");
  leg2->AddEntry(gSDtrkvtx,"~ 4 track vertex (SD)","l");
  if(includeCorr) {
     leg2->AddEntry(gNSDvtx_CORR,"1 pixel vertex (corr.)","pl");
     leg2->AddEntry(gNSDtrkvtx_CORR,"~ 4 track vertex (corr.)","pl");
  }
  leg2->Draw();

  tex->Draw();

  
  char filename1[100];
  if(!includeCorr) sprintf(filename1,"evt_sel_eff_v1");
  else sprintf(filename1,"evt_sel_eff_corr_v1");

     
  char filename2[100];
  if(!includeCorr) sprintf(filename2,"evt_mult_dist_v1");
  else sprintf(filename2,"evt_mult_dist_corr_v1");


  printCanvases(call,filename2,0);
  printCanvases(call2,filename1,0);

  

}

double GetZeroBinFraction(TH1F* hist){

   //cout<<"[GetZeroBinFraction] check: bin center of the first bin : "<<hist->GetBinCenter(1)<<endl;
   
   double total = hist->Integral();
   double non_zerobin = hist->Integral(2,hist->GetNbinsX());
   double fraction = (total-non_zerobin)/total;
   cout<<"total event : "<<total<<endl;

   return fraction;
}
