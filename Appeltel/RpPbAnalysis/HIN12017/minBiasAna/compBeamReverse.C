#include "TrackCorrector.C"

compBeamReverse()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  TrackCorrector t("../corrections/trackCorrections_HIN12017v1_HijingCombined.root");
  t.load("trkCorr_HIN12017");

  //TFile * fmc = new TFile("trackAnalysis_HIN12017v1_HijingJet170closure.root");
  TFile * fpPb = new TFile("trackAnalysis_HIN12017v1_singleTrack210759.root");
  TFile * fPbp = new TFile("trackAnalysis_HIN12017v1_Pbp_singleTrack211607.root");

  double etaMin = -0.16;
  double etaMax = 0.14;
  double etaMinRev = -0.14 ;
  double etaMaxRev = 0.16 ;

  TH3F *hUncorrTrk = (TH3F*) fpPb->Get("trkAna_HIN12017/trkSpectrum");
  hUncorrTrk->SetName("trkSpectrumUncorr");
  TH3F *hTrk = (TH3F*) fpPb->Get("trkAna_HIN12017/trkSpectrum");
  hTrk->SetName("trkSpectrumpPb");
  TH3F *hUncorrTrkRev = (TH3F*) fPbp->Get("trkAna_HIN12017_Pbp/trkSpectrum");
  hUncorrTrkRev->SetName("trkSpectrumUncorrRev");
  TH3F *hTrkRev = (TH3F*) fPbp->Get("trkAna_HIN12017_Pbp/trkSpectrum");
  hTrkRev->SetName("trkSpectrumPbp");

  TH1F * hnevt = fpPb->Get("trkAna_HIN12017/events");
  double nevt = hnevt->Integral();
  TH1F * hnevtRev = fPbp->Get("trkAna_HIN12017_Pbp/events");
  double nevtRev = hnevtRev->Integral();

  // build and apply corrections
  for( int bx=1;bx<=hTrk->GetNbinsX();bx++)
  {
    for( int by=1;by<=hTrk->GetNbinsY();by++)
    {
      for( int bz=1;bz<=hTrk->GetNbinsZ();bz++)
      {
        int bin = hTrk->GetBin(bx,by,bz);
        double x = hTrk->GetXaxis()->GetBinCenter(bx);
        double y = hTrk->GetYaxis()->GetBinCenter(by);
        double z = hTrk->GetZaxis()->GetBinCenter(bz);
        double weight = t.getWeight(y,x,z);
        hTrk->SetBinContent(bin, hTrk->GetBinContent(bin) * weight );
        hTrk->SetBinError(bin, hTrk->GetBinError(bin) * weight );
      }
    }
  }

  for( int bx=1;bx<=hTrkRev->GetNbinsX();bx++)
  {
    for( int by=1;by<=hTrkRev->GetNbinsY();by++)
    {
      for( int bz=1;bz<=hTrkRev->GetNbinsZ();bz++)
      {
        int bin = hTrkRev->GetBin(bx,by,bz);
        double x = hTrkRev->GetXaxis()->GetBinCenter(bx);
        double y = hTrkRev->GetYaxis()->GetBinCenter(by);
        double z = hTrkRev->GetZaxis()->GetBinCenter(bz);
        double weight = t.getWeight(y,x,z);
        hTrkRev->SetBinContent(bin, hTrkRev->GetBinContent(bin) * weight );
        hTrkRev->SetBinError(bin, hTrkRev->GetBinError(bin) * weight );
      }
    }
  }


  // create 1D histos

  TH2F * chTrk = (TH2F*) hTrk->Project3D("yxe");
  TH2F * chUncorrTrk = (TH2F*) hUncorrTrk->Project3D("yxe");
  TH2F * chTrkRev = (TH2F*) hTrkRev->Project3D("yxe");
  TH2F * chUncorrTrkRev = (TH2F*) hUncorrTrkRev->Project3D("yxe");

  Int_t etabin10m=chTrk->GetXaxis()->FindBin(etaMin+0.01);
  Int_t etabin10p=chTrk->GetXaxis()->FindBin(etaMax-0.01);
  Int_t etabin10mRev=chTrkRev->GetXaxis()->FindBin(etaMinRev+0.01);
  Int_t etabin10pRev=chTrkRev->GetXaxis()->FindBin(etaMaxRev-0.01);


  TH1D * hTrkPt = (TH1D *) chTrk->ProjectionY("hTrkPt",etabin10m,etabin10p,"e");
  TH1D * hUncorrTrkPt = (TH1D *) chUncorrTrk->ProjectionY("hUncorrTrkPt",etabin10m,etabin10p,"e");
  TH1D * hTrkPtRev = (TH1D *) chTrkRev->ProjectionY("hTrkPtRev",etabin10mRev,etabin10pRev,"e");
  TH1D * hUncorrTrkPtRev = (TH1D *) chUncorrTrkRev->ProjectionY("hUncorrTrkPtRev",etabin10mRev,etabin10pRev,"e");

  double etaWid = etaMax - etaMin;
  double etaWidRev = etaMaxRev - etaMinRev;

  for(int i = 1; i <= hTrkPt->GetNbinsX(); i++) 
  {
    double binWid = hTrkPt->GetBinWidth(i);
    hTrkPt->SetBinContent(i, hTrkPt->GetBinContent(i)/binWid/etaWid/nevt);
    hTrkPt->SetBinError(i, hTrkPt->GetBinError(i)/binWid/etaWid/nevt);
    hUncorrTrkPt->SetBinContent(i, hUncorrTrkPt->GetBinContent(i)/binWid/etaWid/nevt);
    hUncorrTrkPt->SetBinError(i, hUncorrTrkPt->GetBinError(i)/binWid/etaWid/nevt);
    hTrkPtRev->SetBinContent(i, hTrkPtRev->GetBinContent(i)/binWid/etaWidRev/nevtRev);
    hTrkPtRev->SetBinError(i, hTrkPtRev->GetBinError(i)/binWid/etaWidRev/nevtRev);
    hUncorrTrkPtRev->SetBinContent(i, hUncorrTrkPtRev->GetBinContent(i)/binWid/etaWidRev/nevtRev);
    hUncorrTrkPtRev->SetBinError(i, hUncorrTrkPtRev->GetBinError(i)/binWid/etaWidRev/nevtRev);
  }


  //Tsallis Fit
  TF1 * fit = new TF1("fit","6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.3,1.1);
  fit->SetParameters(10.8327,2.221289,100);
  hTrkPt->Fit("fit","0R");

 
  hTrkPt->SetMarkerStyle(20);
  hUncorrTrkPt->SetMarkerColor(kBlack);
  hUncorrTrkPt->SetMarkerStyle(25);
  hTrkPtRev->SetMarkerStyle(20);
  hTrkPtRev->SetMarkerColor(kRed);
  hUncorrTrkPtRev->SetMarkerColor(kRed);
  hUncorrTrkPtRev->SetMarkerStyle(25);

   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   //c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.16);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.13);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);

   c1->Divide(1,2,0,0);
   c1->cd(1);
   c1->GetPad(1)->SetLogy();
  // c1->GetPad(1)->SetLogx();
  // c1->GetPad(2)->SetLogx();
   

   TH1F * hDum = new TH1F("dum","",10,0.,10.);
   hDum->SetMinimum(0.001);
   hDum->SetMaximum(30);
   hDum->GetYaxis()->SetTitle("Yield");

  TH1D * rTrkPt = (TH1D*) hTrkPt->Clone("ratio");
  TH1D * rTrkPtRev = (TH1D*) hTrkPtRev->Clone("ratioRev");
  rTrkPt->Divide(hTrkPtRev);
  rTrkPtRev->Divide(hTrkPt);
  //TH1D * rTrkPt = (TH1D*) hUncorrTrkPt->Clone("ratio");
  //TH1D * rTrkPtRev = (TH1D*) hUncorrTrkPtRev->Clone("ratioRev");
  //rTrkPt->Divide(hUncorrTrkPtRev);
  //rTrkPtRev->Divide(hUncorrTrkPt);

  hDum->Draw();
  fit->SetLineColor(kGreen);
 // fit->Draw("SAME");
  hTrkPt->Draw("same");
  hUncorrTrkPt->Draw("same");
  hTrkPtRev->Draw("same");
  hUncorrTrkPtRev->Draw("same");
  

  TLegend * leg = new TLegend(0.5,0.6,0.9,0.9);
  leg->AddEntry(hTrkPt,"Corrected Tracks (pPb)","lp");
  leg->AddEntry(hUncorrTrkPt,"Uncorrected Tracks (pPb)","lp");
  leg->AddEntry(hTrkPtRev,"Corrected Tracks (Pbp)","lp");
  leg->AddEntry(hUncorrTrkPtRev,"Uncorrected Tracks (pPbp)","lp");
  leg->SetFillColor(kWhite);
  leg->Draw();
  

  c1->cd(2);

  TH1F * rDum = (TH1F*) hDum->Clone("rdum");
  rDum->SetMinimum(0.9);
  rDum->SetMaximum(1.1);
  rDum->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  rDum->GetXaxis()->SetMoreLogLabels();
  rDum->GetYaxis()->SetTitle("Ratio Pbp / pPb (corrected)");
  rDum->Draw();
  //rTrkPt->Draw("same");
  rTrkPtRev->Draw("same");

  TLatex * tex = new TLatex(0.4,1.05,Form("pPb: %5.3f < #eta_{lab} < %5.3f",etaMin,etaMax));
  tex->SetTextSize(0.043);
  tex->Draw();
  tex = new TLatex(5,1.05,Form("Pbp: %5.3f < #eta_{lab} < %5.3f",etaMinRev,etaMaxRev));
  tex->SetTextSize(0.043);
  tex->Draw();

}

