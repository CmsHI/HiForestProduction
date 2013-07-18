void plotGenRatios() {

  int colorwheel[4] = { kBlack, kRed, kBlue, kGreen };
   
  gStyle->SetOptStat(0);

  set_plot_style();

//  TFile *f = new TFile("validation3D_pion_dataMixed_hiBadgerTracks.root");
//  TFile *f = new TFile("validation3D_bass_25k_superFake.root");

  TFile *f;

  f = new TFile("trackCorrections_jetComp_HIJING_Jet120_50k.root");
//  f = new TFile("trackCorrections_5_1_Y4z.root");

//  char dir[4][256] = {"trkCorr_jet03","trkCorr_jet08","trkCorr_jet03inv","trkCorr_jet08inv"};
//  char dir[4][256] = {"trkCorr_akPu3Calo_Et80dR03","trkCorr_akPu3Calo_Et80dR08",
//                      "trkCorr_akPu3Calo_Et80dR03inv","trkCorr_akPu3Calo_Et80dR08inv"};
  char dir[4][256] = {"trkCorr_akPu3PF_Et80dR03","trkCorr_akPu3PF_Et80dR08",
                      "trkCorr_akPu3PF_Et80dR03inv","trkCorr_akPu3PF_Et80dR08inv"};
//  char dir[4][256] = {"trkCorrPion_akPu3PF_Et80dR03","trkCorrPion_akPu3PF_Et80dR08",
//                      "trkCorrPion_akPu3PF_Et80dR03inv","trkCorrPion_akPu3PF_Et80dR08inv"};

  char jetString[256] = "akPu3PF";
  Double_t jetEtMin = 80.;

  Double_t minCent = 0.;    
  Double_t maxCent = 100.;

  Double_t etaMin = -0.4;
  Double_t etaMax = 0.4;

  Double_t ptMin = 0.2;

  char speciesName[5][256] = { "charged particles","pions","kaons","protons","antiprotons"};
  char histoName[5][256] = {"genPartH","genPartPi","genPartK","genPartp+","genPartp-"};

  TH3F *hSim[5][4];
  TH2F *chSim[5][4];

  for ( int i = 0; i < 4;i++) 
  {
    for( int j=0; j<5; j++)
    {
      hSim[j][i] = (TH3F*) f->Get(Form("%s/%s",dir[i],histoName[j])); 
      hSim[j][i]->GetYaxis()->SetRangeUser(0,100);
      hSim[j][i]->SetName(Form("%s3D%d",histoName[j],i));
      hSim[j][i]->GetZaxis()->SetRangeUser( minCent, maxCent);
      chSim[j][i] = (TH2F*) hSim[j][i]->Project3D("yxe"); 
    }
  }

  // find bins corresponding to projections for below
  Int_t etabin10m=chSim[0][0]->GetXaxis()->FindBin(etaMin+0.01);
  Int_t etabin10p=chSim[0][0]->GetXaxis()->FindBin(etaMax-0.01);

  TH1D * hSimPt[5][4];
  TGraphAsymmErrors * gRat[5][4];

  for( int i = 0; i<4; i++ )
  {
    for( int j=0; j<5; j++)
    {
      hSimPt[j][i] = (TH1D*) chSim[j][i]->ProjectionY(Form("%sPtP%d",histoName[j],i),etabin10m,etabin10p,"e");
    }

    for( int j=1; j<5; j++)
    {
      gRat[j][i] = new TGraphAsymmErrors(); gRat[j][i]->SetName(Form("%sRat%d",histoName[j],i));
      gRat[j][i]->BayesDivide(hSimPt[j][i],hSimPt[0][i]);
      gRat[j][i]->SetMarkerStyle(20);
      gRat[j][i]->SetLineStyle(2);
      gRat[j][i]->SetLineColor(colorwheel[i]);
      gRat[j][i]->SetMarkerColor(colorwheel[i]);
    }
  }

  TLegend *leg = new TLegend(0.22,0.32,0.85,0.75);
  leg->SetFillColor(0); leg->SetBorderSize(1);
  leg->AddEntry(gRat[1][0],"dR < 0.3","lp");
  leg->AddEntry(gRat[1][1],"dR < 0.8","lp");
  leg->AddEntry(gRat[1][2],"dR > 0.3","lp");
  leg->AddEntry(gRat[1][3],"dR > 0.8","lp");
 
  TH1D* hDumPt = new TH1D("hDumPt",";p_{T} [GeV/c]",40,0.0,100.0); hDumPt->SetMaximum(1.0);
  hDumPt->SetMinimum(1e-2);
  hDumPt->GetXaxis()->SetMoreLogLabels(); 
  hDumPt->GetYaxis()->SetMoreLogLabels(); 
  hDumPt->GetYaxis()->SetTitleOffset(1.8);
  hDumPt->GetXaxis()->SetTitleOffset(1.2);

  TH1D * hDumPtC[5];

  TLatex * tex;

  TCanvas *c8 = new TCanvas("c8","Particle Compostion",1100,500);
  c8->Divide(4,1,0,0);
  for( int j=1; j<5;j++)
  {
    c8->cd(j);
    c8->GetPad(j)->SetLogx();
  //  c8->GetPad(j)->SetLogy();
    hDumPtC[j]=(TH1D*) hDumPt->Clone(Form("hDumPt%d",j));
    hDumPtC[j]->GetYaxis()->SetTitle("Ratio");
    hDumPtC[j]->SetTitle("");
    hDumPtC[j]->Draw(); 
    for( int i=0; i<4; i++ )
      gRat[j][i]->Draw("p");
    
    if(j==1) tex = new TLatex(0.5,0.25,Form("%s / %s",speciesName[j],speciesName[0]));
    if(j>1) tex = new TLatex(0.5,0.9,Form("%s / %s",speciesName[j],speciesName[0]));
    tex->SetTextSize(0.055);
    tex->Draw();
    if(j==1) tex = new TLatex(0.5,0.2,Form("%3.1f < #eta <  %3.1f",etaMin,etaMax));
    if(j==1) tex->SetTextSize(0.055);
    if(j==1) tex->Draw();
    if(j==1) tex = new TLatex(0.5,0.12,Form("%s , E_{T}>%4.0f",jetString,jetEtMin));
    if(j==1) tex->SetTextSize(0.055);
    if(j==1) tex->Draw();
    leg->SetTextSize(0.055);
    if(j==4) leg->Draw();
  } 
 
}


void set_plot_style() {

  // nicer colz plots

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

