//-------------------------------------------------
// June 02 2010, Andre S. Yoon
// plot track pT and leading jet Et cooreation  
// MC vs DATA
//-------------------------------------------------

#include "commonStyle.h"

void PlotPtvsJetEt(){


   gStyle->SetPadRightMargin(0.13); // to show rainbow scale


   //------------ Canvas && dummy hist 
   TCanvas *cMC = new TCanvas("cMC","cMC",490,550);
   TCanvas *cDATA = new TCanvas("cDATA","cDATA",490,550);

   Char_t xTitle[100],yTitle[100];
   sprintf(yTitle,"Track p_{T} [GeV/c]");
   sprintf(xTitle,"Corrected leading Jet E_{T} [GeV]");

   TH1D *dum_MC = GetDummyHist(0,200,0.2,300,xTitle,yTitle);
   TH1D *dum_DATA = GetDummyHist(0,200,0.2,300,xTitle,yTitle);

   TLatex * tex;
   tex = new TLatex(0.55,0.90,"CMS Preliminary");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();


   //------------- MC
   TFile *fMC = new TFile("./../root_files/TrkHistMC_QCD_10M.root");
   TH3F *hTrkPtEtaJetEt_MC = (TH3F*) fMC->Get("trackAna_mc/threeDHist/hTrkPtEtaJetEt");

   TH1F *hEvt_MC = (TH1F*) fMC->Get("trackAna_mc/hNevt");
   float nevt_mc = hEvt_MC->GetEntries();

   cMC->cd();
   cMC->SetLogz();
   cMC->SetLogy();
   dum_MC->Draw();
   TH2D *hPtEt_MC = (TH2D*) hTrkPtEtaJetEt_MC->Project3D("yze");//->Draw("COLZ");
   //hPtEt_MC->GetYaxis()->SetRangeUser(0.01,300);
   hPtEt_MC->Draw("COLZsame");

   TH2D *hPtEt_MC_clone = (TH2D*) hPtEt_MC->Clone("hPtEt_MC_clone");
   TH2D *hPtEt_MC_clone_2GEV = (TH2D*) hPtEt_MC->Clone("hPtEt_MC_clone_2GEV");

   TH1F *hMean_MC = new TH1F("hMean_MC","hMean_MC",15,0,300);
   TH1F *hMean_MC_2GEV = new TH1F("hMean_MC_2GEV","hMean_MC_2GEV",15,0,300);

   TH1F *hMult_MC = new TH1F("hMult_MC","hMult_MC",15,0,300);

   for(int i=0;i<15;i++){
      hPtEt_MC_clone->GetXaxis()->SetRange(i+1,i+2);
      float mean = hPtEt_MC_clone->GetMean(2);
      float mult = hPtEt_MC_clone->ProjectionY()->GetEntries();
      mult = mult/nevt_mc;
      mult = mult*1000;
      
      hMean_MC->Fill(hMean_MC->GetBinCenter(i+1),mean);
      hMult_MC->Fill(hMult_MC->GetBinCenter(i+1),mult);

      hPtEt_MC_clone_2GEV->GetXaxis()->SetRange(i+1,i+2);
      hPtEt_MC_clone_2GEV->GetYaxis()->SetRange(10,1000);
      float mean_2GEV = hPtEt_MC_clone_2GEV->GetMean(2);
      hMean_MC_2GEV->Fill(hMean_MC_2GEV->GetBinCenter(i+1),mean_2GEV);
   }
   th1Style1(hMean_MC,2,24,1.5,2,3,9,2);
   th1Style1(hMean_MC_2GEV,2,24,1.5,2,3,2,2);
   //th1Style1(hMult_MC,4,24,1.5,4,3,9,2);


   char label_MC[100];
   sprintf(label_MC,"PYTHIA 7 TeV");
   TLegend *leg_MC = (TLegend*) GetDummyLegend(label_MC,hPtEt_MC,hMean_MC,hMean_MC_2GEV);
   leg_MC->Draw();

   tex->Draw();

   //-------------- DATA       
   TFile *fDATA = new TFile("./../root_files/MB-C10-Merged.root");
   TH3F *hTrkPtEtaJetEt_DATA = (TH3F*) fDATA->Get("trackAna/threeDHist/hTrkPtEtaJetEt");
   
   cDATA->cd();
   cDATA->SetLogy();
   cDATA->SetLogz();
   dum_DATA->Draw();
   TH2D *hPtEt_DATA = (TH2D*) hTrkPtEtaJetEt_DATA->Project3D("yze");
   hPtEt_DATA->Draw("COLZsame");

   TH2D *hPtEt_DATA_clone = (TH2D*) hPtEt_DATA->Clone("hPtEt_DATA_clone");
   TH2D *hPtEt_DATA_clone_2GEV = (TH2D*) hPtEt_DATA->Clone("hPtEt_DATA_clone_2GEV");


   TH1F *hMean_DATA = new TH1F("hMean_DATA","hMean_DATA",15,0,300);
   TH1F *hMean_DATA_2GEV = new TH1F("hMean_DATA_2GEV","hMean_DATA_2GEV",15,0,300);

   for(int i=0;i<15;i++){
      hPtEt_DATA_clone->GetXaxis()->SetRange(i+1,i+2);
      float mean = hPtEt_DATA_clone->GetMean(2);
      hMean_DATA->Fill(hMean_DATA->GetBinCenter(i+1),mean);

      hPtEt_DATA_clone_2GEV->GetXaxis()->SetRange(i+1,i+2);
      hPtEt_DATA_clone_2GEV->GetYaxis()->SetRange(10,1000);
      float mean_2GEV = hPtEt_DATA_clone_2GEV->GetMean(2);
      hMean_DATA_2GEV->Fill(hMean_DATA_2GEV->GetBinCenter(i+1),mean_2GEV);
   }
   th1Style1(hMean_DATA,2,24,1.5,2,3,9,2);
   th1Style1(hMean_DATA_2GEV,2,24,1.5,2,3,2,2);




   char label_DATA[100];
   sprintf(label_DATA,"DATA 7 TeV");
   TLegend *leg_DATA = (TLegend*) GetDummyLegend(label_DATA,hPtEt_DATA,hMean_DATA,hMean_DATA_2GEV);
   leg_DATA->Draw();

   tex->Draw();


   //--------------- Save figures
   char filenameMC[100];
   sprintf(filenameMC,"pt_jet_correlation_MC");
   printCanvases(cMC,filenameMC,0);
   
   char filenameDATA[100];
   sprintf(filenameDATA,"pt_jet_correlation_DATA");
   printCanvases(cDATA,filenameDATA,0);



}


TH1D* GetDummyHist(float xmin, float xmax, float ymin, float ymax, char *xttl, char *yttl) {

   TH1D *dum;
   dum = new TH1D("dum","",100,xmin,xmax);

   dum->SetMinimum(ymin);
   dum->SetMaximum(ymax);
   dum->SetStats(0);

   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();

   return dum;

}

TLegend* GetDummyLegend(char *label, TH1* hist, TH1* hist2, TH1* hist3){

   TLegend *leg1 = new TLegend(0.18,0.74,0.55,0.93);
   leg1->SetMargin(0.3);
   leg1->SetBorderSize(0);
   leg1->SetFillColor(kWhite);
   leg1->SetFillStyle(0);
   leg1->SetTextSize(0.038);
   leg1->SetMargin(0.3);
   leg1->AddEntry(hist,label,"");
   leg1->AddEntry(hist2,"mean p_{T}","l");
   leg1->AddEntry(hist3,"mean p_{T} (p_{T}>2.0 GeV/c)","l");
   return leg1;
}
