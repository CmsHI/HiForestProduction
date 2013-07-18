void fit_pi0(){
  
  gStyle->SetOptStat(0);

  TFile *f=new TFile("merged_JetAnalysisTTrees_hiGoodTightTracks_v1.root");

  TTree *t=f->Get("pi0Analyzer/t");


  TH1F *h6=new TH1F("h6","h6",50,0.,0.3);
  
  //  t->Draw("invMass>>h6","bin>=20&&bin<24&&abs(pi0_eta)<1.4&&pi0_pt>3.0&&pi0_pt<4.0");
  t->Draw("invMass>>h6","bin>=20&&bin<36&&abs(pi0_eta)<0.6&&pi0_pt>2.0");

  TF1 *f4=new TF1("tfit4","[0]*exp(-(x-[1])*(x-[1])/2.0/[2]/[2])+[3]*x+[4]+x*x*[5]",0,0.3);
  f4->SetParameters(50,0.14,0.02,1,1,1);    
    
  h6->SetTitle("#pi^{0} Mass Peak from PF photons");
  h6->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
  h6->SetYTitle("Counts per 6 MeV/c^{2} bin");
  
  h6->Fit(f4,"","",0.05,0.25);

  TF1 *fgaus=new TF1("fgaus","[0]*exp(-(x-[1])*(x-[1])/2.0/[2]/[2])",0,0.3);
  TF1 *fbkg=new TF1("fbkg","[0]*x+[1]+x*x*[2]",0.,0.3);
  fgaus->SetParameters(f4->GetParameter(0),f4->GetParameter(1),f4->GetParameter(2));
  fbkg->SetParameters(f4->GetParameter(3),f4->GetParameter(4),f4->GetParameter(5));
 

  float sigma = f4->GetParameter(2);
  float mean = f4->GetParameter(1);

  float signal = fgaus->Integral(mean-3.*sigma,mean+3.*sigma);
  float bkg = fbkg->Integral(mean-3.*sigma,mean+3.*sigma);
  cout<<" S/B "<< signal/bkg<<endl;

}
