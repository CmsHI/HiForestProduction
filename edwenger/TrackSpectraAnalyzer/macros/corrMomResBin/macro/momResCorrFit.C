void momResCorrFit()
{
  float ptMin=0,ptMax=170;
  TCanvas * c2 = new TCanvas("c2","c2",500,500);
  TF1 * fResCorr = new TF1("fResCorr","[0]+[1]*x+[2]*x^2",ptMin,ptMax);
  fResCorr->SetParameters(1.00379,2.62093e-05,7.21870e-07);
  fResCorr->SetMinimum(0.9);
  fResCorr->SetMaximum(1.1);
  fResCorr->Draw();

  TF1 * fResCorr2 = new TF1("fResCorr2","[0]+[1]*x+[2]*x^2",ptMin,ptMax);
  fResCorr2->SetParameters(1.00379,2.58589e-05,7.3465e-07);
  fResCorr2->SetLineColor(2);
  fResCorr2->Draw("same");

  TF1 * fResCorrUpper = new TF1("fResCorrUpper","[0]+[1]*x+[2]*x^2",ptMin,ptMax);
  fResCorrUpper->SetParameters(1.00379,6.58589e-05,1.3465e-06);
  fResCorrUpper->SetLineColor(kBlue);
  fResCorrUpper->Draw("same");

  TF1 * fResCorrLower = new TF1("fResCorrLower","[0]+[1]*x+[2]*x^2",ptMin,ptMax);
  fResCorrLower->SetParameters(1.00379,-2.58589e-05,5.3465e-07);
  fResCorrLower->SetLineColor(kGreen+2);
  fResCorrLower->Draw("same");

  TCanvas * c3 = new TCanvas("c3","c3",500,500);
  TF1 *fRat = new TF1("fRat","fResCorr2/fResCorr",ptMin,ptMax);
  fRat->Draw();
}
