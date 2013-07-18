#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <vector>

void formatHist(TH1* h, int col = 1, double norm = 1);
void saveCanvas(TCanvas* c, int date = 20080829);
TH1F* myDivide(TH1F* a,TH1F* b);

void analyze_cut(const char* infile = "p0829.root", double etaMax = 1.){
  const char* outfile = "output.root";

  /// Parameters
  double normRange = 0.6;
  double deltaCut = 0.2;
  int etaBins = 8;
  int nBins = 20;

  //bool useDR = true;

  gROOT->Reset();
  gROOT->ProcessLine(".x rootlogon.C");
  gStyle->SetErrorX(0.);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  gStyle->SetPadLeftMargin(0.132);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.02);

  TH1::SetDefaultSumw2();

  // get Ntuple from input file
  TFile *f = new TFile(infile);
  f->cd("ana;1");
  TNtuple * ntparticle= dynamic_cast<TNtuple *>(f->Get("ana/ntparticle"));
  TNtuple * ntmatched= dynamic_cast<TNtuple *>(f->Get("ana/ntmatched"));
  TNtuple * ntInvMatched= dynamic_cast<TNtuple *>(f->Get("ana/ntInvMatched"));
  TNtuple * ntgen = dynamic_cast<TNtuple *>(f->Get("ana/ntgen"));
  TNtuple * ntevent = dynamic_cast<TNtuple *>(f->Get("ana/ntevent"));
  TNtuple * ntvertex = dynamic_cast<TNtuple *>(f->Get("ana/ntvertex"));
 

  // set up output file
  TFile *of = new TFile(outfile,"recreate");

  TNtuple * ntcorr = new TNtuple("ntcorr","","nhad:ntrt:nstrt:npass:nhit:vz:nvtx:vntrk:beta:alpha");
  TH1F * h1 = new TH1F("h1",Form("Everything;D;#_{pixel pairs}/event/%.2f",2/(double)nBins),nBins,0,2); 
  TH1F * h2 = new TH1F("h2","Signal",nBins,0,2);
  TH1F * h3 = new TH1F("h3","Background",nBins,0,2);
  TH1F * h4 = new TH1F("h4","Normalized Reproduced Background",nBins,0,2);
  TH1F * h5 = new TH1F("h5","",nBins,0,2);
  TH1F * h6 = new TH1F("h6","Reproduced Background Subtracted",nBins,0,2);
  TH1F * h7 = new TH1F("h7","",nBins,0,2);
  TH1F * h8 = new TH1F("h8","",nBins,0,2);
  TH1F * h9 = new TH1F("h9","",nBins,0,2);
  TH1F * h10 = new TH1F("h10","",nBins,0,2);

  TProfile * dNdEtaHadron = new TProfile("dNdEtaHadron","",32,-2.1,2.1);
  TProfile * dNdEtaLepton = new TProfile("dNdEtaLepton","",32,-2.1,2.1);
  TProfile * dNdEtaTracklet = new TProfile("dNdEtaTracklet","",32,-2.1,2.1);

  TH2D * corr = new TH2D("correlation","; #_{hadrons}; #_{tracklets}",100,0,50,100,0,50);

  Float_t matchedeta1;
  Float_t matchedeta2;
  Float_t matchedinveta2;
  Float_t matchedinvphi2;
  Float_t signalcheck;

  Float_t charge;
  Float_t ntpartlayer1hits;
  Float_t matchedphi1;
  Float_t matchedphi2;
 
  Float_t eta1;
  Float_t eta2;
  Float_t deta;
  Float_t dphi;
  Float_t invdeta;
  Float_t invdphi;
  Float_t inveta1;
  Float_t inveta2;
  Float_t invphi1;
  Float_t invphi2;
  Float_t vz;
  Float_t vntrk;
  Float_t nvtx;
  Float_t evtId;
  Float_t eventEvtId;


  ntvertex->SetBranchAddress("z",&vz);
  ntvertex->SetBranchAddress("n",&vntrk);
  ntvertex->SetBranchAddress("nvtx",&nvtx);

  ntmatched->SetBranchAddress("eta1",&matchedeta1);
  ntmatched->SetBranchAddress("phi1",&matchedphi1);
  ntmatched->SetBranchAddress("matchedeta",&matchedeta2);
  ntmatched->SetBranchAddress("matchedphi",&matchedphi2);
  ntmatched->SetBranchAddress("signalCheck",&signalcheck);
  ntmatched->SetBranchAddress("deta",&deta);
  ntmatched->SetBranchAddress("dphi",&dphi);
  ntmatched->SetBranchAddress("evtid",&evtId);

  ntInvMatched->SetBranchAddress("eta1",&inveta1);
  ntInvMatched->SetBranchAddress("phi1",&invphi1);
  ntInvMatched->SetBranchAddress("matchedeta",&inveta2);
  ntInvMatched->SetBranchAddress("matchedphi",&invphi2);
  ntInvMatched->SetBranchAddress("deta",&invdeta);
  ntInvMatched->SetBranchAddress("dphi",&invdphi);


  ntevent->SetBranchAddress("evtid",&eventEvtId);


  int nevents = ntgen->GetEntries();
  int matchedentries = ntmatched->GetEntries();

  // number of passed tracklets

  vector<int> npass;

  int nPassCounter=0;
  int oldEvtId=0; 
  int evts=0;
  double dRCut = 0.05;

  // Get the first eventId
  ntevent->GetEntry(0);

  for(int i = 0; i<matchedentries;i++){

    ntmatched->GetEntry(i);

    float dR= sqrt(deta*deta+dphi/43.*dphi/43.);

    // Check if this is a tracklet which belong to a new event
    if (evtId != (int)eventEvtId) {
       while (evtId!=(int)eventEvtId){
          // Check if the event is empty
          npass.push_back(nPassCounter); 

          cout <<evts<<" "<<eventEvtId<<" "<<evtId<<" "<<nPassCounter<<endl;
          nPassCounter = 0;
          oldEvtId++;
          evts++;
          ntevent->GetEntry(evts);
       }
    }

    // Cuts
    if(fabs(matchedeta1)>etaMax) continue;
    if (dR<dRCut) nPassCounter++;

    h1->Fill(fabs(dR));
    h6->Fill(fabs(dR));
    if(signalcheck==1){
      h2->Fill(fabs(dR));
    }
    if(signalcheck==0){
      h3->Fill(fabs(dR));
    } 
  }

  npass.push_back(nPassCounter); 


  int invmatchedentries = ntInvMatched->GetEntries();

  for(int i = 0; i<invmatchedentries;i++){
    ntInvMatched->GetEntry(i);
    if(fabs(inveta1)>etaMax) continue;
    if(fabs(inveta1)<0.1) continue;
    float dR= sqrt(invdeta*invdeta+invdphi*invdphi);
    h4->Fill(dR);
  }
   
  formatHist(h1,1,nevents);
  formatHist(h2,2,nevents);
  formatHist(h3,3,nevents);
  formatHist(h4,4,nevents);
  formatHist(h6,6,nevents);

  //// Normalization of background

  Float_t sc = ((h1->Integral((int)(normRange*nBins),nBins,"width"))/(h4->Integral((int)(normRange*nBins),nBins,"width")));


  cout<<"background normalization: "<<sc<<endl;
  h4->Scale(sc);
  h6->Add(h4,-1);

  //// Determination of correction factor beta
  double beta = 1-((h2->Integral(0,(int)(deltaCut*nBins),"width"))/(h6->Integral(0,(int)(deltaCut*nBins),"width")));
  cout<<"beta: "<<beta<<endl;

  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->SetLogy();

  TLegend * leg1 = new TLegend(0.25,0.66,0.56,0.84);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(h1,"Everything","l");
  leg1->AddEntry(h4,"Normalized Reproduced Combinatorial Background","l");
  leg1->AddEntry(h6,"Reproduced Combinatorial Background Subtracted","l");

  h1->Draw("e");
  h2->Draw("e same");
  h4->Draw("e same");
  h6->Draw("e same");
  leg1->Draw();

  // Might be necessary to increase the sizes here
  float had[8] = {0,0,0,0,0,0,0,0};
  float lep[8] = {0,0,0,0,0,0,0,0};
  float trt[8] = {0,0,0,0,0,0,0,0};
  float strt[8] = {0,0,0,0,0,0,0,0};
  float layer1hits[8] = {0,0,0,0,0,0,0,0};

  //Important information:
  //
  //The ntuple is filled in event-by-event basis
  //therefore the total number of particles in a certain
  //eta bin is recorded by the hardcoded variable names
  //had1,had2,had3 etc...
  //Here this is corrected to be from 0 to 7
  //The way to interpret the etabins is this:
  //had[0] = had1 : -2, -1.5
  //had[1] = had2 : -1.5, 1
  //...
  //had[6] = had7 : 1, 1.5
  //had[7] = had8 : 1.5, 2
  // had: number of hadrons
  // lep: number of leptons
  // trt: number of tracklets
  // strt: number of signal tracklets

  for(int ig = 0; ig < etaBins; ++ig){
    ntgen->SetBranchAddress(Form("had%d",ig+1),&(had[ig]));
    ntgen->SetBranchAddress(Form("lep%d",ig+1),&(lep[ig]));
    ntevent->SetBranchAddress(Form("trt%d",ig+1),&(trt[ig]));
    ntevent->SetBranchAddress(Form("strt%d",ig+1),&(strt[ig]));
    ntevent->SetBranchAddress(Form("hit%d",ig+1),&(layer1hits[ig]));
    //ntevent->SetBranchAddress(Form("had%d",ig+1),&(trt[ig])); // in the old run
  }

  if (npass.size()!= (unsigned int)nevents) {
      cout <<"Different number of npasses and nevents!"<<npass.size()<<" "<<nevents<<endl;
  }

  for(int iev = 0; iev < nevents; ++iev){
    ntgen->GetEntry(iev);
    ntevent->GetEntry(iev);
    ntvertex->GetEntry(iev);
    double nhad = 0;
    double ntrt = 0;
    double nstrt = 0;
    double nhit = 0;
    // take |eta|<1
    for(int ieta = 2; ieta < 6; ++ieta){
      dNdEtaHadron->Fill((0.5*(double)ieta-1.75),had[ieta]);
      dNdEtaLepton->Fill((0.5*(double)ieta-1.75),lep[ieta]);
      dNdEtaTracklet->Fill((0.5*(double)ieta-1.75),trt[ieta]);
      nhad += had[ieta];
      ntrt += trt[ieta];
      nstrt += strt[ieta];
      nhit += layer1hits[ieta];
    }

    //// define alpha and beta by modeling, to be replaced with Yetkin's class
    double beta = 0;//0.07111+0.006029848*(nhit)-4.319779e-4*(nhit)*(nhit)+1.265031e-5*(nhit)*(nhit)*(nhit);

    double xx = nhit + 0.5;

    double alpha = 1.141827      +
                  (-5.656305e-2) * xx +
                  (3.800675e-3 ) * xx * xx +
		  (-1.164418e-4) * xx * xx * xx +
		  (1.221609e-6 ) * xx * xx * xx * xx;
		  
    if (vz!=0) {
       corr->Fill(nhad,npass[iev]);
       ntcorr->Fill(nhad,ntrt,nstrt,npass[iev],nhit,vz,nvtx,vntrk,beta,alpha);
    }
  }

  TCanvas* c2 = new TCanvas("Correlation","",400,400);
  corr->Draw("col");

  TCanvas* c5 = new TCanvas("Profile","",400,400);

  TProfile* p1 = new TProfile("p1","",20,0,40);

  ntcorr->Draw("npass*(1-beta)*alpha:nhad>>p1","","prof");
  
  p1->SetXTitle("Number of Hadrons");
  p1->SetYTitle("#alpha (1-#beta) N_{Tracklet}");
  TLine *line = new TLine(0,0,80,80);
  p1->SetAxisRange(0,40,"Y");
  line->SetLineColor(2);
  line->Draw();  

  TCanvas* c6 = new TCanvas("Tracklet Profile","",400,400);

  TProfile* p2 = new TProfile("p2","",20,0,40);

  ntcorr->Draw("ntrt:nhad>>p2","","prof");
  
  p2->SetXTitle("Number of Hadrons");
  p2->SetYTitle("Number of Tracklet");
  p2->Fit("pol2","","");
  TLine *line2 = new TLine(0,0,80,80);
  line2->SetLineColor(2);
  line2->Draw();  

  TCanvas* c7 = new TCanvas("nhadron over npass","",400,400);

  TProfile* p30 = new TProfile("p30","",60,0,60);

  p30->SetMarkerSize(0.7);

  TProfile* p31 = new TProfile("p31","",60,0,60);

  p31->SetMarkerSize(0.7);

  ntcorr->Draw("nhad:nhit>>p30","","prof");
  ntcorr->Draw("npass*(1-beta):nhit>>p31","","prof");
   
  TH1F *p3 = myDivide((TH1F*)p30,(TH1F*)p31);
    
  p3->SetXTitle("Number of Hits in the first layer (|#eta|<1)");
  p3->SetYTitle("N_{Hadron} / (N_{Tracklet}#times(1-#beta)) ");
  p3->Fit("pol4","","",1,40);

  TCanvas* c8 = new TCanvas("nmeasured over npass","",400,400);

  TProfile* p4 = new TProfile("p4","",40,0,40);

  p4->SetMarkerSize(0.7);
  ntcorr->Draw("npass*(1-beta)*alpha/nhad:nhit>>p4","","prof");
  
  p4->SetXTitle("Number of Hits in the first layer (|#eta|<1)");
  p4->SetYTitle("n^{Obs}_{Hadron} / N_{Hadron} ");
  p4->Fit("pol1","","",1,40);
  

  saveCanvas(c1);
  saveCanvas(c2);
  saveCanvas(c5);
  saveCanvas(c6);
  saveCanvas(c7);
  ntcorr->Write();


}

void formatHist(TH1* h, int col, double norm){

  h->Scale(1/norm);
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

}

void saveCanvas(TCanvas* c, int date){
  c->Write();
  c->Draw();
  /*
  c->Print(Form("./figures/%s_d%d.gif",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.eps",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.C",c->GetName(),date));
  */
}

TH1F* myDivide(TH1F* a,TH1F* b)
{
   TH1F *h = new TH1F("h","",a->GetNbinsX(),a->GetXaxis()->GetXmin(),a->GetXaxis()->GetXmax());

   for(int i=0;i<h->GetNbinsX();i++) 
   {
      Double_t X = a->GetBinContent(i);
      Double_t Xerr = a->GetBinError(i);
      Double_t Y = b->GetBinContent(i);
      Double_t Yerr = b->GetBinError(i);
      if (Xerr==0) Xerr=X;
      if (Yerr==0) Yerr=Y;
      
      Double_t R = X/Y;
      Double_t Rerr = R*sqrt((Xerr/X)*(Xerr/X)+(Yerr/Y)*(Yerr/Y));   
      if (X==0||Y==0) {
         R=0;
         Rerr=0;
      }
      cout <<i<<" "<<R<<" "<<Rerr<<" "<<X<<" "<<Y<<endl;
      h->SetBinContent(i,R);
      h->SetBinError(i,Rerr);
   }
   return h;
}
