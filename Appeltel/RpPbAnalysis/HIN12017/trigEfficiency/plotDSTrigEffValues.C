{
  gStyle->SetOptStat(0);

  TFile * f = new TFile("trackAnalysis_HIN12017v1_HijingMBevlsel.root");

  TH1F * sel = (TH1F*) f->Get("trkAna_HIN12017_MC/multDS");
  TH1F * nosel = (TH1F*) f->Get("trkAna_HIN12017_MC_noSel/multDS");
  sel->Sumw2(); nosel->Sumw2();
  sel->Divide(nosel);
  sel->SetMarkerColor(kBlue);
  sel->SetLineColor(kBlue);
  sel->SetMarkerStyle(20);
  sel->SetMarkerSize(0.8);


  TH1F * nonDS = (TH1F*) f->Get("trkAna_HIN12017_MC/multNonDS");
  TH1F * all = (TH1F*) f->Get("trkAna_HIN12017_MC/multiplicity");
  sel->Sumw2(); all->Sumw2();
  nonDS->Divide(all);
  nonDS->SetMarkerColor(kBlue);
  nonDS->SetLineColor(kBlue);
  nonDS->SetMarkerStyle(20);
  nonDS->SetMarkerSize(0.8);

  TH1F * hDumE = new TH1F("hDumE","",10,0,50);
  TH1F * hDumC = hDumE->Clone("hDumC");
  
  TCanvas * c1 = new TCanvas("c1","Trigger Efficiency",800,500);
  c1->Divide(2,1);
  c1->cd(1);

  hDumE->SetMaximum(1.2);
  hDumE->SetMinimum(0.0);
  hDumE->GetYaxis()->SetTitle("Trigger Efficiency (DS)");
  hDumE->GetXaxis()->SetTitle("Event Multiplicity (Selected Tracks)");

  hDumE->Draw();
  sel->Draw("same");

  c1->cd(2);

  hDumC->SetMaximum(0.04);
  hDumC->SetMinimum(0.0);
  hDumC->GetYaxis()->SetTitle("Non-DS Contamination");
  hDumC->GetXaxis()->SetTitle("Event Multiplicity (Selected Tracks)");

  hDumC->Draw();
  nonDS->Draw("same");

  cout << "Zero-Multiplicity Event Fraction = " 
       << nosel->Integral(1,1) / nosel->Integral() << "\n\n";

  cout << "Trigger Efficiency e(M) : \n"; 
  for( int i=1; i<31; i++)
  {
    if( i % 8 == 0) cout << "\n        ";
    cout << sel->GetBinContent(i) << ", ";
  }
  cout << "\n\n";
 
  cout << "Non-DS Contamination f(M) : \n";
  for( int i=1; i<31; i++)
  {
    if( i % 8 == 0) cout << "\n        ";
    cout << nonDS->GetBinContent(i) << ", ";
  }
  cout << "\n\n";

}
