#include "TrackCorrector.C"

compHIN12016()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  TrackCorrector t("../corrections/trackCorrections_HIN12017v1_HijingCombined.root");
  t.load("trkCorr_HIN12017");
  //TFile * fmc = new TFile("trackAnalysis_HIN12017v1_HijingJet170closure.root");
  //TFile * fmc = new TFile("trackAnalysis_HIN12017v1_singleTrack210759.root");
  TFile * fmc = new TFile("trackAnalysis_HIN12017v1_singleTrack211607_quick.root");

  double etaMin = -0.16;
  double etaMax = 0.14;
//  double etaMin = 0.74;
//  double etaMax = 1.04;


  TH3F *hUncorrTrk = (TH3F*) fmc->Get("trkAna_HIN12017/trkSpectrum");
  hUncorrTrk->SetName("trkSpectrumUncorr");
  TH3F *hTrk = (TH3F*) fmc->Get("trkAna_HIN12017/trkSpectrum");
 // TH3F * hTrk = hUncorrTrk->Clone("correctedTrk");

  TF3 * unity = new TF3("unity","1.0+0.0*x*y*z",-1000.,1000.,-1000.,1000.,-1000.,1000.);

  TH1F * hnevt = fmc->Get("trkAna_HIN12017/events");
  double nevt = hnevt->Integral();

  // apply corrections as scaling factors

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

  // create 1D histos

  TH2F * chTrk = (TH2F*) hTrk->Project3D("yxe");
  TH2F * chUncorrTrk = (TH2F*) hUncorrTrk->Project3D("yxe");

  Int_t etabin10m=chTrk->GetXaxis()->FindBin(etaMin+0.01);
  Int_t etabin10p=chTrk->GetXaxis()->FindBin(etaMax-0.01);


  TH1D * hTrkPt = (TH1D *) chTrk->ProjectionY("hTrkPt",etabin10m,etabin10p,"e");
  TH1D * hUncorrTrkPt = (TH1D *) chUncorrTrk->ProjectionY("hUncorrTrkPt",etabin10m,etabin10p,"e");

  double etaWid = etaMax - etaMin;

  for(int i = 1; i <= hTrkPt->GetNbinsX(); i++) 
  {
    double binWid = hTrkPt->GetBinWidth(i);
    hTrkPt->SetBinContent(i, hTrkPt->GetBinContent(i)/binWid/etaWid/nevt);
    hTrkPt->SetBinError(i, hTrkPt->GetBinError(i)/binWid/etaWid/nevt);
    hUncorrTrkPt->SetBinContent(i, hUncorrTrkPt->GetBinContent(i)/binWid/etaWid/nevt);
    hUncorrTrkPt->SetBinError(i, hUncorrTrkPt->GetBinError(i)/binWid/etaWid/nevt);
  }

  // load the HIN-12-016 data

  double v1,v2,v3,v4;
  double xa[100],xe[100],ya[100],ye[100];  

  ifstream infile;
  infile.open("dndpt_HIN12016.dat");
     int np=0;
  if (infile.is_open())
  {
    std::string line;
    while(std::getline(infile,line,'\n'))
    {
      int res = sscanf(line.c_str(),"%lg %lg %lg",&v1,&v2,&v3);
      if (res!=3) continue;
      xa[np]=v1;
      xe[np]=0.0;
      ya[np]=v2;
      ye[np]=v3;
      np++;
    }
  }
  TGraphErrors *ge = new TGraphErrors(np,&xa[0],&ya[0],&xe[0],&ye[0]);
  ge->SetMarkerStyle(20);
  ge->SetMarkerColor(kRed);

  //Tsallis Fit
  TF1 * fit = new TF1("fit","6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.3,3.1);
  fit->SetParameters(10.8327,2.221289,100);
  ge->Fit("fit","0R");

 
  hTrkPt->SetMarkerStyle(20);
  hUncorrTrkPt->SetMarkerColor(kBlue);
  hUncorrTrkPt->SetMarkerStyle(25);

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
  // c1->GetPad(1)->SetLogy();
  // c1->GetPad(1)->SetLogx();
  // c1->GetPad(2)->SetLogx();

   TH1F * hDum = new TH1F("dum","",10,0.,3.);
   hDum->SetMinimum(0);
   hDum->SetMaximum(30);
   hDum->GetYaxis()->SetTitle("Yield");

  TH1D * rTrkPt = (TH1D*) hTrkPt->Clone("ratio");
//  rTrkPt->Divide(hGenPt);
  rTrkPt->Divide(fit);

  hDum->Draw();
  fit->SetLineColor(kGreen);
  fit->Draw("SAME");
  hTrkPt->Draw("same");
  hUncorrTrkPt->Draw("same");
  ge->Draw("p");
  

  TLegend * leg = new TLegend(0.2,0.2,0.5,0.4);
  leg->AddEntry(hTrkPt,"Corrected Tracks","lp");
  leg->AddEntry(hUncorrTrkPt,"Uncorrected Tracks","lp");
  leg->AddEntry(ge,"HIN-12-016","lp");
  leg->Draw();

  c1->cd(2);

  TH1F * rDum = (TH1F*) hDum->Clone("rdum");
  rDum->SetMinimum(0.9);
  rDum->SetMaximum(1.1);
  rDum->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  rDum->GetXaxis()->SetMoreLogLabels();
  rDum->GetYaxis()->SetTitle("Ratio Trk / Fit");
  rDum->Draw();
  rTrkPt->Draw("same");

  TGraphErrors * rge = ge->Clone("rge");
  for( int i=0;i<np;i++)
  {
    Double_t x,y;
    rge->GetPoint(i,x,y);
    rge->SetPoint(i,x,y / fit->Eval(x) );
  }
  rge->Draw("p");

  TLatex * tex = new TLatex(0.4,1.05,Form("%5.3f < #eta < %5.3f",etaMin,etaMax));
  tex->SetTextSize(0.043);
  tex->Draw();

}

