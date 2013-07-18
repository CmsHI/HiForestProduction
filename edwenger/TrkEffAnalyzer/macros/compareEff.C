void drawAll(TString type="Acc") {

  TCanvas *cAll = new TCanvas("cAll","cAll",800,800);
  cAll->Divide(2,2);
  cAll->cd(1);
  compareEff(Form("%sEta",type.Data()));
  cAll->cd(2);
  compareEff(Form("%sEta2",type.Data()));
  cAll->cd(3);
  compareEff(Form("%sPt",type.Data()));
  cAll->cd(4);
  compareEff(Form("%sPt2",type.Data()));

  //TCanvas *cDum = new TCanvas("cDum","cDum",500,500);

}


void compareEff(TString str="AccEta") {


  TFile *f356 = new TFile("trkEffD6T_356.root");
  TFile *fD6T = new TFile("trkEffD6T.root");
  TFile *fP0 = new TFile("trkEffP0.root");
  TFile *fPY8 = new TFile("trkEffPY8.root");

  TGraphAsymmErrors *g356 = (TGraphAsymmErrors*) f356->Get(Form("g%s",str.Data()));
  TGraphAsymmErrors *gD6T = (TGraphAsymmErrors*) fD6T->Get(Form("g%s",str.Data()));
  TGraphAsymmErrors *gP0 = (TGraphAsymmErrors*) fP0->Get(Form("g%s",str.Data()));
  TGraphAsymmErrors *gPY8 = (TGraphAsymmErrors*) fPY8->Get(Form("g%s",str.Data()));

  //g356->Draw("apz");
  //gD6T->Draw("pz");
  //gP0->Draw("pz");
  //gPY8->Draw("pz");

  TString title="Ratio of ";
  if(str.Contains("Eff")) title+= "Efficiencies (";
  else if(str.Contains("Acc")) title+= "Acceptances (";
  else if(str.Contains("Fak")) title+= "Fake Rates (";
  else if(str.Contains("Sec")) title+= "Secondary Rates (";

  if(str.Contains("Eta2")) title+= "p_{T} > 2 GeV/c);#eta";
  else if(str.Contains("Eta")) title+= "p_{T} > 0.4 GeV/c);#eta";
  else if(str.Contains("Pt2")) title+= "|#eta| < 2.4);p_{T}";
  else if(str.Contains("Pt")) title+= "|#eta| < 1.0);p_{T}";

  TH1F *hDumEta = new TH1F(Form("hDumEta_%s",str.Data()),title.Data(),100,-2.4,2.4);
  if(str.Contains("Eff") || str.Contains("Acc")) {hDumEta->SetMinimum(0.95); hDumEta->SetMaximum(1.05);}
  else if(str.Contains("Fak") || str.Contains("Sec")) {hDumEta->SetMinimum(0.5); hDumEta->SetMaximum(1.5);}
  hDumEta->SetStats(0);

  TH1F *hDumPt = new TH1F(Form("hDumPt_%s",str.Data()),title.Data(),100,0,10);
  if(str.Contains("Eff") || str.Contains("Acc")) {hDumPt->SetMinimum(0.95); hDumPt->SetMaximum(1.05);}
  if(str.Contains("Fak") || str.Contains("Sec")) {hDumPt->SetMinimum(0.5); hDumPt->SetMaximum(1.5);}
  hDumPt->SetStats(0);

  if(str.Contains("Eta")) hDumEta->Draw();
  else if(str.Contains("Pt")) hDumPt->Draw();
   
  TGraphAsymmErrors *rD6T = divideGraphs(gD6T,g356,Form("rD6T_%s",str.Data()));
  TGraphAsymmErrors *rP0 = divideGraphs(gP0,g356,Form("rP0_%s",str.Data()));
  TGraphAsymmErrors *rPY8 = divideGraphs(gPY8,g356,Form("rPY8_%s",str.Data()));
  
  rD6T->SetMarkerColor(1);
  rD6T->Draw("pz");
  rP0->SetMarkerColor(2);
  rP0->Draw("pz");
  rPY8->SetMarkerColor(4);
  rPY8->Draw("pz");
  
  TLegend *leg = new TLegend(0.3,0.2,0.7,0.4);
  leg->SetFillColor(0);
  leg->AddEntry(rD6T,"D6T real vtx / D6T wide vtx","lp");
  leg->AddEntry(rP0,"P0 real vtx / D6T wide vtx","lp");
  leg->AddEntry(rPY8,"PY8 real vtx / D6T wide vtx","lp");
  leg->Draw("same");


  TCanvas *cDum = new TCanvas("cDum","cDum",500,500);

}

TGraphAsymmErrors* divideGraphs(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TString name, Double_t offsetx=0) {

  Int_t N = g1->GetN();
  cout << N << " points" << endl;
  Double_t x3[50];
  Double_t y3[50];
  Double_t e3h[50], e3l[50], e3x[50];

  Double_t x1=0.0, y1=0.0, x2=0.0, y2=0.0;
  Double_t e1h=0.0, e1l=0.0, e2h=0.0, e2l=0.0;

  for(Int_t i=0; i<TMath::Min(N,50); i++) {
    g1->GetPoint(i,x1,y1);
    if(y1>0) { e1h = g1->GetErrorYhigh(i); e1l = g1->GetErrorYlow(i); }
    else { e1h=0.0; e1l=0.0; }
    g2->GetPoint(i,x2,y2);
    if(y2>0) { e2h = g2->GetErrorYhigh(i); e2l = g2->GetErrorYlow(i); }
    else { e2h=0.0; e2l=0.0; }
    x3[i]=x1+offsetx;
    if(y2>0) y3[i]=y1/y2;
    else y3[i]=0.0;
    e3l[i]=y3[i]*TMath::Sqrt(e1l*e1l+e2l*e2l);
    e3h[i]=y3[i]*TMath::Sqrt(e1h*e1h+e2h*e2h);
    e3x[i]=0.0;
    //cout << "Point " << i << ": x=" << x1 << ", y=" << y1 << "/" << y2 << "=" << y3[i] << endl;
    //cout << "Error " << i << ": low=" << e3l[i] << ", high=" << e3h[i] << endl;

  }

  TGraphAsymmErrors *g3 = new TGraphAsymmErrors(TMath::Min(N,50),x3,y3,e3x,e3x,e3l,e3h);
  g3->SetName(name.Data());
  g3->SetMarkerStyle(20);
  
  return g3;

}
