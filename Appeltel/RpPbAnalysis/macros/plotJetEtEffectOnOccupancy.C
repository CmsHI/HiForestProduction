void plotJetEtEffectOnOccupancy() {

  int colorwheel[7] = { kBlack, kRed, kBlue, kGreen, kViolet, kGray, kOrange };
   
  gStyle->SetOptStat(0);

  set_plot_style();

//  TFile *f = new TFile("validation3D_pion_dataMixed_hiBadgerTracks.root");
//  TFile *f = new TFile("validation3D_bass_25k_superFake.root");

  TFile *f;

  f = new TFile("trackCorrections_jetOcc_HIJING_Jet280_batch6.root");
//  f = new TFile("trackCorrections_5_1_Y4z.root");

//  char dir[4][256] = {"trkCorr_jet03","trkCorr_jet08","trkCorr_jet03inv","trkCorr_jet08inv"};
//  char dir[4][256] = {"trkCorr_akPu3Calo_Et80dR03","trkCorr_akPu3Calo_Et80dR08",
//                      "trkCorr_akPu3Calo_Et80dR03inv","trkCorr_akPu3Calo_Et80dR08inv"};
//  char dir[4][256] = {"trkCorr_akPu3PF_Et80dR03","trkCorr_akPu3PF_Et80dR08",
//                      "trkCorr_akPu3PF_Et80dR03inv","trkCorr_akPu3PF_Et80dR08inv"};
  char dir[256] = "trkCorr_akPu3PF_dR03";

  double jetEtMin[7] = {0.,40.,80.,120.,160.,200.,250.};
  double jetEtMax[7] = {40.,80.,120.,160.,200.,250.,500.};

  Double_t etaMin = -0.4;
  Double_t etaMax = 0.4;

  Double_t ptMin = 0.2;

  // sim-to-reco hists
  TH3F *hSim[7];
  TH3F *hEff[7];
  TH3F *hRec[7];
  TH3F *hFak[7];

  for ( int i = 0; i < 7;i++)
  {
    hSim[i] = (TH3F*) f->Get(Form("%s/hsim3D",dir)); hSim[i]->GetYaxis()->SetRangeUser(0,100);
    hSim[i]->SetName(Form("hsim3D%d",i));
    hEff[i] = (TH3F*) f->Get(Form("%s/heff3D",dir)); hEff[i]->GetYaxis()->SetRangeUser(0,100);
    hEff[i]->SetName(Form("heff3D%d",i));
    hRec[i] = (TH3F*) f->Get(Form("%s/hrec3D",dir)); hRec[i]->GetYaxis()->SetRangeUser(0,100);
    hRec[i]->SetName(Form("hrec3D%d",i));
    hFak[i] = (TH3F*) f->Get(Form("%s/hfak3D",dir)); hFak[i]->GetYaxis()->SetRangeUser(0,100);
    hFak[i]->SetName(Form("hfak3D%d",i));
  }



  // Projection for centrality

  TH2F *chSim[7];
  TH2F *chRec[7];
  TH2F *chFak[7];
  TH2F *chEff[7];
  double iota = 1e-3;
 
  for( int i = 0; i<7; i++)
  {
    hSim[i]->GetZaxis()->SetRangeUser( jetEtMin[i]+iota, jetEtMax[i]-iota);
    chSim[i] = (TH2F*) hSim[i]->Project3D("yxe"); 
    hEff[i]->GetZaxis()->SetRangeUser( jetEtMin[i]+iota, jetEtMax[i]-iota);
    chEff[i] = (TH2F*) hEff[i]->Project3D("yxe"); 
    hRec[i]->GetZaxis()->SetRangeUser( jetEtMin[i]+iota, jetEtMax[i]-iota);
    chRec[i] = (TH2F*) hRec[i]->Project3D("yxe"); 
    hFak[i]->GetZaxis()->SetRangeUser( jetEtMin[i]+iota, jetEtMax[i]-iota);
    chFak[i] = (TH2F*) hFak[i]->Project3D("yxe"); 
  }

  // ratio histograms
  TH2F *rEff[7];  
  TH2F *rRec[7];  
  TH2F *rFak[7];  
  TH2F *rSim[7];  

  for( int i = 0; i<7; i++ )
  {
    rEff[i] = (TH2F*) chEff[i]->Clone(Form("rEff%d",i));
    rRec[i] = (TH2F*) chRec[i]->Clone(Form("rMul%d",i));
    rFak[i] = (TH2F*) chFak[i]->Clone(Form("rFak%d",i));
    rSim[i] = (TH2F*) chSim[i]->Clone(Form("rSec%d",i));
  }
   

  //---------------------------------------------


  // find bins corresponding to projections for below
  Int_t ptbin04=chSim[0]->GetYaxis()->FindBin(ptMin+0.01);
  Int_t ptbins=chSim[0]->GetYaxis()->GetNbins();
  Int_t etabin10m=chSim[0]->GetXaxis()->FindBin(etaMin+0.01);
  Int_t etabin10p=chSim[0]->GetXaxis()->FindBin(etaMax-0.01);

  // projected hists: pt 
  TH1D * hSimEta[7];
  TH1D * hEffEta[7];
  TH1D * hRecEta[7];
  TH1D * hFakEta[7];

  for ( int i = 0; i<7; i++)
  {
    hSimEta[i] = (TH1D*) chSim[i]->ProjectionX(Form("hSimEta%d",i),ptbin04,ptbins,"e");
    hEffEta[i] = (TH1D*) chEff[i]->ProjectionX(Form("hEffEta%d",i),ptbin04,ptbins,"e");
    hRecEta[i] = (TH1D*) chRec[i]->ProjectionX(Form("hRecEta%d",i),ptbin04,ptbins,"e");
    hFakEta[i] = (TH1D*) chFak[i]->ProjectionX(Form("hFakEta%d",i),ptbin04,ptbins,"e");
  }


  TH1D* hDumEta = new TH1D("hDumEta",";#eta",60,-2.4,2.4); hDumEta->SetMaximum(1.0);
   hDumEta->GetXaxis()->CenterTitle(); 
  
  hDumEta->GetYaxis()->SetTitleOffset(1.8);
  TH1D* hDumEta2 = (TH1D*) hDumEta->Clone("hDumEta2"); hDumEta2->SetMaximum(0.2); 
  TH1D* hDumEta3 = (TH1D*) hDumEta->Clone("hDumEta3"); hDumEta3->SetMaximum(0.00049); 


  // projected hists: eta 
  TH1D* hSimPt[7];
  TH1D* hFakPt[7];
  TH1D* hRecPt[7];
  TH1D* hEffPt[7];

  for( int i = 0; i<7; i++ )
  {
    hSimPt[i]  = (TH1D*) chSim[i]->ProjectionY(Form("hSimPtP%d",i),etabin10m,etabin10p,"e");
    hRecPt[i]  = (TH1D*) chRec[i]->ProjectionY(Form("hRecPtP%d",i),etabin10m,etabin10p,"e");
    hFakPt[i]  = (TH1D*) chFak[i]->ProjectionY(Form("hFakPtP%d",i),etabin10m,etabin10p,"e");
    hEffPt[i]  = (TH1D*) chEff[i]->ProjectionY(Form("hEffPtP%d",i),etabin10m,etabin10p,"e");
  }

  
  TH1D* hDumPt = new TH1D("hDumPt",";p_{T} [GeV/c]",40,0.0,100.0); 
  hDumPt->SetMaximum(1.0);
  hDumPt->GetXaxis()->SetMoreLogLabels(); 
  hDumPt->GetYaxis()->SetTitleOffset(1.8);
  hDumPt->GetXaxis()->SetTitleOffset(1.2);

  // Efficiency
  TGraphAsymmErrors *gEffEta[7];
  TGraphAsymmErrors *gEffPt[7];
  
  for( int i = 0; i<7; i++ )
  {
    gEffEta[i] = new TGraphAsymmErrors(); gEffEta[i]->SetName(Form("gEffEta%d",i));
    gEffEta[i]->BayesDivide(hEffEta[i],hSimEta[i]);
    gEffEta[i]->SetMarkerStyle(20);
    gEffEta[i]->SetLineStyle(2);
    gEffEta[i]->SetLineColor(colorwheel[i]);
    gEffEta[i]->SetMarkerColor(colorwheel[i]);

    gEffPt[i] = new TGraphAsymmErrors(); gEffPt[i]->SetName(Form("gEffPt%d",i));
    gEffPt[i]->BayesDivide(hEffPt[i],hSimPt[i]);
    gEffPt[i]->SetMarkerStyle(20);
    gEffPt[i]->SetLineStyle(2);
    gEffPt[i]->SetLineColor(colorwheel[i]);
    gEffPt[i]->SetMarkerColor(colorwheel[i]);

  }


  TLegend *legEta = new TLegend(0.72,0.72,0.95,0.95);
  legEta->SetFillColor(0); legEta->SetBorderSize(1);
  for( int i=0;i<7;i++)
    legEta->AddEntry(gEffEta[i],Form("%4.1f < Jet-E_{T} < %4.1f",
                                     jetEtMin[i], jetEtMax[i] ) ,"lp");

  TLegend *legPt = new TLegend(0.52,0.62,0.85,0.85);
  legPt->SetFillColor(0); legPt->SetBorderSize(1);
  for( int i=0;i<7;i++)
    legPt->AddEntry(gEffPt[i],Form("%4.1f < Jet-E_{T} < %4.1f",
                                     jetEtMin[i], jetEtMax[i] ) ,"lp");

  TCanvas *c7 = new TCanvas("c7","Efficiency Fraction",900,500);
//  c7->Divide(2,1);
  TCanvas *c8 = new TCanvas("c8","Efficiency Fraction",900,500);
  c8->Divide(2,1);
  hDumEtaEff=(TH1F*) hDumEta->Clone("hDumEtaEff");
  hDumEtaEff->GetYaxis()->SetTitle("Absolute efficiency");
  hDumPtEff=(TH1F*) hDumPt->Clone("hDumPtEff");
  hDumPtEff->GetYaxis()->SetTitle("Absolute efficiency");
  hDumPtEff->SetTitle(Form("Absolute Efficiency %4.1f < |#eta| < %4.1f", 
                           etaMin,etaMax));
  hDumEtaEff->SetTitle(Form("Absolute Efficiency p_{T} > %4.2f GeV/c"
                            ,ptMin));
  c7->cd(); hDumEtaEff->Draw(); 
  for( int i=0; i<7; i++ )
    gEffEta[i]->Draw("p");
  
    legEta->Draw();
 
  c8->GetPad(1)->SetLogx();
  c8->GetPad(2)->SetLogx();
 
  c8->cd(1); hDumPtEff->Draw(); 
  for( int i=0; i<7; i++)
    gEffPt[i]->Draw("p"); 

//   legPt->Draw();


  // Fakes

  TGraphAsymmErrors *gFakEta[7];
  TGraphAsymmErrors *gFakPt[7];
  
  for( int i = 0; i<7; i++ )
  {
    gFakEta[i] = new TGraphAsymmErrors(); gFakEta[i]->SetName(Form("gFakEta%d",i));
    gFakEta[i]->BayesDivide(hFakEta[i],hRecEta[i]);
    gFakEta[i]->SetMarkerStyle(20);
    gFakEta[i]->SetLineStyle(2);
    gFakEta[i]->SetLineColor(colorwheel[i]);
    gFakEta[i]->SetMarkerColor(colorwheel[i]);

    gFakPt[i] = new TGraphAsymmErrors(); gFakPt[i]->SetName(Form("gFakPt%d",i));
    gFakPt[i]->BayesDivide(hFakPt[i],hRecPt[i]);
    gFakPt[i]->SetMarkerStyle(20);
    gFakPt[i]->SetLineStyle(2);
    gFakPt[i]->SetLineColor(colorwheel[i]);
    gFakPt[i]->SetMarkerColor(colorwheel[i]);

  }


  TCanvas *c9 = new TCanvas("c9","Fake Fraction",900,500);
  TCanvas *c10 = new TCanvas("c10","Fake Fraction",900,500);
  hDumEtaFak=(TH1F*) hDumEta->Clone("hDumEtaMul");
  hDumEtaFak->GetYaxis()->SetTitle("Fake Reconstruction Fraction");
  hDumPtFak=(TH1F*) hDumPt->Clone("hDumPtMul");
  hDumPtFak->GetYaxis()->SetTitle("Fake Reconstruction Fraction");
  hDumPtFak->SetTitle(Form("Fake Rate %4.1f < |#eta| < %4.1f",
                           etaMin,etaMax));
  hDumEtaFak->SetTitle(Form("Fake Rate p_{T} > %4.2f GeV/c"
                            ,ptMin));
  c9->cd(); hDumEtaFak->Draw(); 
  for( int i=0; i<7; i++)  
    gFakEta[i]->Draw("p"); 
   legEta->Draw();
  c8->cd(2); hDumPtFak->Draw(); 
  for ( int i=0; i<7; i++)
    gFakPt[i]->Draw("p");  
  legPt->Draw();


  /*
  TFile *f = new TFile("trkEffPY8.root","RECREATE");
  gAccPt->Write(); gAccPt2->Write(); gAccEta->Write(); gAccEta2->Write();
  gEffPt->Write(); gEffPt2->Write(); gEffEta->Write(); gEffEta2->Write();
  gMulPt->Write(); gMulPt2->Write(); gMulEta->Write(); gMulEta2->Write();
  gFakPt->Write(); gFakPt2->Write(); gFakEta->Write(); gFakEta2->Write();
  gSecPt->Write(); gSecPt2->Write(); gSecEta->Write(); gSecEta2->Write();
  f->Close();
  */
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

