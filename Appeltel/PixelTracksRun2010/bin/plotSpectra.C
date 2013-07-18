{

  TCanvas * c = new TCanvas("c1","c1",800.,600.);
  c->SetFillColor(kWhite);
  c->SetFrameFillColor(kWhite);
  c->SetBorderMode(0);
  c->SetFrameBorderMode(0);
  c->cd();

//  c->SetLogy();

  double norm = 400.0;

  TFile * f1 = TFile::Open("pixel_PROMPT_200ev.root");
  TH1F * data = (TH1F *) f1->Get("pixel/pixelpt");

  cout << "Data dN/deta is " << data->Integral()/norm << "\n";
  cout << "Data <pt> is " << data->GetMean() * 1000 << " MeV/c\n";

  Double_t data_x[100];
  Double_t data_y[100];
  Double_t data_yerr[100];
  Double_t data_xerr[100];

  for( int i=0;i<100;i++)
  {
    data_x[i] = ((float)i)*0.05 + 0.025;
    data_y[i] = data->GetBinContent(i)/norm;
    data_yerr[i] = data->GetBinError(i)/norm;
    data_xerr[i] = 0.0;
  } 

  TGraphErrors dG( 100, data_x, data_y, data_xerr, data_yerr );

  dG.SetTitle("200 Event Uncorrected Pixel Track Spectrum  |#eta|<1.0");
  dG.GetXaxis()->SetTitle("p_{T} (GeV/c)");
  dG.GetYaxis()->SetTitle("#frac{d^{2}N}{d#eta dp_{T}} ");

  dG.SetMarkerSize(0.7);
  dG.SetMarkerStyle(20);

  dG.Draw("ap");


  TFile * f2 = TFile::Open("pixel_BASS_200ev.root");
  TH1F * bass = (TH1F *) f2->Get("pixel/pixelpt");

  cout << "Bass dN/deta is " << bass->Integral()/norm << "\n";
  cout << "Bass <pt> is " << bass->GetMean() * 1000 << " MeV/c\n";


  for( int i=0;i<100;i++)
  {
    data_x[i] = ((float)i)*0.05 + 0.025;
    data_y[i] = bass->GetBinContent(i)/norm;
    data_yerr[i] = bass->GetBinError(i)/norm;
    data_xerr[i] = 0.0;
  }

  TGraphErrors dGbass( 100, data_x, data_y, data_xerr, data_yerr );

  dGbass.SetMarkerSize(0.7);
  dGbass.SetMarkerStyle(24);
  dGbass.SetMarkerColor(kBlue);
  dGbass.SetLineColor(kBlue);


  dGbass.Draw("p");

  TFile * f3 = TFile::Open("pixel_ORGAN_200ev.root");
  TH1F * organ = (TH1F *) f3->Get("pixel/pixelpt");

  cout << "Organ dN/deta is " << organ->Integral()/norm << "\n";
  cout << "Organ <pt> is " << organ->GetMean() * 1000 << " MeV/c\n";


  for( int i=0;i<100;i++)
  {
    data_x[i] = ((float)i)*0.05 + 0.025;
    data_y[i] = organ->GetBinContent(i)/norm;
    data_yerr[i] = organ->GetBinError(i)/norm;
    data_xerr[i] = 0.0;
  }

  TGraphErrors dGorgan( 100, data_x, data_y, data_xerr, data_yerr );

  dGorgan.SetMarkerSize(0.7);
  dGorgan.SetMarkerStyle(24);
  dGorgan.SetMarkerColor(kRed);
  dGorgan.SetLineColor(kRed);


  dGorgan.Draw("p");

  dG.GetYaxis()->SetRangeUser(0.0, 30.0);

  TLegend leg(0.40,0.65,0.98,0.85);


  TH1F * l1 = new TH1F("l1","l1",10,0.,10.);
  l1->SetLineColor(kBlack);
  l1->SetMarkerColor(kBlack);
  l1->SetMarkerSize(0.7);
  l1->SetMarkerStyle(20);
   

  TH1F * l2 = new TH1F("l2","l2",10,0.,10.);
  l2->SetLineColor(kBlue);
  l2->SetMarkerColor(kBlue);
  l2->SetMarkerSize(0.7);
  l2->SetMarkerStyle(24);

  TH1F * l3 = new TH1F("l3","l3",10,0.,10.);
  l3->SetLineColor(kRed);
  l3->SetMarkerColor(kRed);
  l3->SetMarkerSize(0.7);
  l3->SetMarkerStyle(20);


  leg.AddEntry(l1,"Data Run 150476 HLT_HIMinBiasHF");
  leg.AddEntry(l2,"Hydjet Bass");
  leg.AddEntry(l3,"AMPT Organ");

  leg.Draw();

}

