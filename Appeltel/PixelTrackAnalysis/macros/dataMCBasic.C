dataMCBasic() {

  TCanvas * can1 = new TCanvas("c1","c1",1600.,1200.);
  can1->SetFillColor(kWhite);
  can1->SetFrameFillColor(kWhite);
  can1->SetBorderMode(0);
  can1->SetFrameBorderMode(0);
  can1->Divide(2,2);

  TFile * f[3];

  TFile * f[0] = TFile::Open("pixel_BASS_10000ev.root");
  TFile * f[1] = TFile::Open("pixel_ORGAN_10000ev.root");
  TFile * f[2] = TFile::Open("pixel_RUN150476_10000ev.root");

  TH1F * pt[3];
  TH1F * eta[3];
  TH1F * phi[3];
  TH1F * d0[3];
  TH1F * d0err[3];
  TH1F * dz[3];
  TH1F * dzerr[3];

  TGraphErrors * gpt[3];
  TGraphErrors * geta[3];
  TGraphErrors * gphi[3];
  TGraphErrors * gd0[3];
  TGraphErrors * gd0err[3];
  TGraphErrors * gdz[3];
  TGraphErrors * gdzerr[3];

  int colorwheel[3] = { kBlue, kRed, kBlack };

  for ( int i = 0; i<3; i++)
  {
      pt[i] = (TH1F *) f[i]->Get("pixelTrackAna/pixelpt");  
      phi[i] = (TH1F *) f[i]->Get("pixelTrackAna/pixelphi");  
      eta[i] = (TH1F *) f[i]->Get("pixelTrackAna/pixeleta");  
      d0[i] = (TH1F *) f[i]->Get("pixelTrackAna/pixeld0");  
      d0err[i] = (TH1F *) f[i]->Get("pixelTrackAna/pixeld0err");  
      dz[i] = (TH1F *) f[i]->Get("pixelTrackAna/pixeldz");  
      dzerr[i] = (TH1F *) f[i]->Get("pixelTrackAna/pixeldzerr");  
  }


  TLegend * leg = new TLegend(0.40,0.65,0.98,0.85);
  leg->SetFillColor(0);

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


  leg->AddEntry(l1,"Data Run 150476 HLT_HIMinBiasHF");
  leg->AddEntry(l2,"Hydjet Bass");
  leg->AddEntry(l3,"AMPT Organ");


  for ( int i = 0; i<3; i++)
  {
      gpt[i] = convertToTGraphErrors( pt[i], true, colorwheel[i]);
      gphi[i] = convertToTGraphErrors( phi[i], true, colorwheel[i]);
      geta[i] = convertToTGraphErrors( eta[i], true, colorwheel[i]);
      gd0[i] = convertToTGraphErrors( d0[i], true, colorwheel[i]);
      gd0err[i] = convertToTGraphErrors( d0err[i], true, colorwheel[i]);
      gdz[i] = convertToTGraphErrors( dz[i], true, colorwheel[i]);
      gdzerr[i] = convertToTGraphErrors( dzerr[i], true, colorwheel[i]);
  }

  gpt[0]->SetTitle("");
  gpt[0]->GetYaxis()->SetTitle("frequency");
  gpt[0]->GetXaxis()->SetTitle("p_{T} ( GeV/C )");

  gphi[0]->SetTitle("");
  gphi[0]->GetYaxis()->SetTitle("frequency");
  gphi[0]->GetXaxis()->SetTitle("#phi");

  geta[0]->SetTitle("");
  geta[0]->GetYaxis()->SetTitle("frequency");
  geta[0]->GetXaxis()->SetTitle("#eta");

  gd0[0]->SetTitle("");
  gd0[0]->GetYaxis()->SetTitle("frequency");
  gd0[0]->GetXaxis()->SetTitle("Transverse DCA (cm)");

  gd0err[0]->SetTitle("");
  gd0err[0]->GetYaxis()->SetTitle("frequency");
  gd0err[0]->GetXaxis()->SetTitle("Transverse DCA Significance");

  gdz[0]->SetTitle("");
  gdz[0]->GetYaxis()->SetTitle("frequency");
  gdz[0]->GetXaxis()->SetTitle("Longitudinal DCA (cm)");

  gdzerr[0]->SetTitle("");
  gdzerr[0]->GetYaxis()->SetTitle("frequency");
  gdzerr[0]->GetXaxis()->SetTitle("Longitudinal DCA Significance");

  can1->cd(1);
  gpt[0]->Draw("ap");
  gpt[1]->Draw("p");
  gpt[2]->Draw("p");
  leg->Draw();

  can1->cd(2);
  gphi[0]->Draw("ap");
  gphi[1]->Draw("p");
  gphi[2]->Draw("p");

  can1->cd(3);
  geta[0]->Draw("ap");
  geta[1]->Draw("p");
  geta[2]->Draw("p");





  TCanvas * can2 = new TCanvas("c2","c2",1600.,1200.);  can2->SetFillColor(kWhite);
  can2->SetFrameFillColor(kWhite);
  can2->SetBorderMode(0);
  can2->SetFrameBorderMode(0);  can2->cd();
  can2->Divide(2,2);

  can2->cd(1);
  gd0[0]->Draw("ap");
  gd0[1]->Draw("p");
  gd0[2]->Draw("p");

  can2->cd(2);
  gd0err[0]->Draw("ap");
  gd0err[1]->Draw("p");
  gd0err[2]->Draw("p");
  leg->Draw();

  can2->cd(3);
  gdz[0]->Draw("ap");
  gdz[1]->Draw("p");
  gdz[2]->Draw("p");

  can2->cd(4);
  gdzerr[0]->Draw("ap");
  gdzerr[1]->Draw("p");
  gdzerr[2]->Draw("p");
}


TGraphErrors *  convertToTGraphErrors( TH1 * h, bool normalize, int color )
{

   if (!h) {
      cout << "Pointer to histogram is null\n";
      return NULL;
   }
   if (h->GetDimension() != 1) {
      cout << Form("Histogram must be 1-D; h %s is %d-D\n",h->GetName(),h->GetDimension());
      return NULL;
   }
      
   int npoints = ((TH1*)h)->GetXaxis()->GetNbins();

   if ( npoints  > 10000 )
   {
      cout << "Too many points, fix it in the source or tell Eric to write better macros.\n";
      return NULL;
   }

   double norm = 1.;

   if( normalize ) 
   {
     norm = h->Integral();
   }

   Double_t x[10000], y[10000], xErr[10000], yErr[10000];

   TAxis *xaxis = ((TH1*)h)->GetXaxis();
   for (Int_t i=0;i<npoints;i++) {
      x[i] = xaxis->GetBinCenter(i+1);
      y[i] = h->GetBinContent(i+1) / norm;
      xErr[i] = 0;
      yErr[i] = h->GetBinError(i+1) / norm;
   }

   TGraphErrors * g = new TGraphErrors(npoints, x, y, xErr, yErr );

   h->TAttLine::Copy(*g);
   h->TAttFill::Copy(*g);
   h->TAttMarker::Copy(*g);

   std::string gname = "Graph_from_" + std::string(h->GetName() );
   g->SetName(gname.c_str());
   g->SetTitle(h->GetTitle());

   g->SetLineColor(color);
   g->SetMarkerColor(color);
   g->SetMarkerSize(0.7);
   g->SetMarkerStyle(20);

   return g;

}

