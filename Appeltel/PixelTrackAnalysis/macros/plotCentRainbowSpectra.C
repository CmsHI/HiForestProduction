{

  TCanvas * can = new TCanvas("c1","c1",800.,600.);
  can->SetFillColor(kWhite);
  can->SetFrameFillColor(kWhite);
  can->SetBorderMode(0);
  can->SetFrameBorderMode(0);
  can->cd();

  can->SetLogy();

  int colors[10] = { kBlack, kOrange+3, kRed, kOrange+1, kYellow, kGreen, kCyan, kBlue, kViolet, kGray+2 };

  TFile * f1 = TFile::Open("pixels_DATA_1000ev.root");

  TH1I * cbins = (TH1I *) f1->Get("pixelTrackAna/centbins");

  TH1F * data[10]; 
  TGraphErrors * dCent[10];

  for( int c=0; c<10; c++)
  {
  
    data[c] = (TH1F *) f1->Get(Form("pixelTrackAna/pixelptcent%d",c));

    double norm = 2 * cbins->Integral(4*c, 4*c+3);

    cout << "cbin = " << c <<  " dN/deta is " << data[c]->Integral()/norm << "\n";
    cout << "cbin = " << c << " <pt> is " << data[c]->GetMean() * 1000 << " MeV/c\n";
  
    Double_t data_x[100];
    Double_t data_y[100];
    Double_t data_yerr[100];
    Double_t data_xerr[100];


    for( int i=0;i<100;i++)
    {
      data_x[i] = ((float)i)*0.05 + 0.025;
      data_y[i] = data[c]->GetBinContent(i)/norm;
      data_yerr[i] = data[c]->GetBinError(i)/norm;
      data_xerr[i] = 0.0;
    }  

    dCent[c] = new TGraphErrors( 100, data_x, data_y, data_xerr, data_yerr );
  
    dCent[c]->SetMarkerSize(0.4);
    dCent[c]->SetMarkerStyle(20);
    dCent[c]->SetMarkerColor(colors[c]);
    dCent[c]->SetLineColor(colors[c]);

    if ( c == 0 )
    {
      dCent[0]->SetTitle("1000 Event Uncorrected Pixel Track Spectrum by Centrality  |#eta|<1.0");
      dCent[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      dCent[0]->GetYaxis()->SetTitle("#frac{d^{2}N}{d#eta dp_{T}} ");
      dCent[0]->Draw("ap");
    } else {
      dCent[c]->Draw("p");
    }


  }



  TLegend leg(0.60,0.35,0.98,0.90);

  TH1F * l[10];

  for( int c = 0; c<10; c++ )
  {
    l[c] = new TH1F(Form("l%d",c),Form("l%d",c),10,0.,10.);
    l[c]->SetLineColor(colors[c]);
    l[c]->SetMarkerColor(colors[c]);
    l[c]->SetMarkerSize(0.8);
    l[c]->SetMarkerStyle(20);

    leg.AddEntry(l[c],Form("Centrality %d0-%d0",c,c+1));
  }    




  leg.Draw();

}

