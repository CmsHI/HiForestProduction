#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TF1.h"                                                                                                             

#include "Tracklet.h"
#include "Math/Vector3D.h"
#include "MitHigFunctions.h"


using namespace std;

void analyze_vertex(char * infile = "/Users/yetkinyilmaz/data/pythia_mb_900GeV_contaminated_vtxFlat_d20081031.root", char * outfile = "vertex900cont.root"){

  gROOT->Reset();
  gROOT->ProcessLine(".x rootlogon.C");
  gStyle->SetErrorX(0.);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetMarkerSize(0.7);


  TFile* inf = new  TFile(infile);
  TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));

  int maxEvents = 1000000000;

  int zbins = 1;
  int hitbins = 100;
  int nbins = zbins*hitbins;
  
  // Selection on Hits and events
  SelectionCriteria cuts;
  cuts.drCut   = 0.4;      // to remove double hit
  cuts.dPhiCut = 0.04;     // to remove double hit
  cuts.vzCut   = 10;       // vertex cut

  cuts.verbose_ = false;
  cuts.useDeltaPhi_ = false;
  cuts.checkSecondLayer_ = true;
  
  // Output PDF
  TFile* outf = new TFile(outfile,"recreate");
  TNtuple *ntmatched = new TNtuple("ntmatched","","eta1:matchedeta:phi1:matchedphi:deta:dphi:signalCheck:tid:r1id:r2id:evtid:nhit1:sid:ptype:vz:eta:phi:pt:id");
  TNtuple *ntmult = new TNtuple("ntmult","","mult:nhit1:nhit2");
  
  TTree *trackletTree = new TTree("TrackletTree","Tree of Reconstructed Tracklets");

  TrackletData tdata;

  trackletTree->Branch("nTracklets",&tdata.nTracklet,"nTracklets/I");
  trackletTree->Branch("nhit1",&tdata.nhit1,"nhit1/I");
  trackletTree->Branch("nhit2",&tdata.nhit2,"nhit2/I");
  trackletTree->Branch("mult",&tdata.mult,"mult/I");
  trackletTree->Branch("nv",&tdata.nv,"nv/I");
  trackletTree->Branch("vz",tdata.vz,"vz[nv]/F");
  trackletTree->Branch("eta1",tdata.eta1,"eta1[nTracklets]/F");
  trackletTree->Branch("phi1",tdata.phi1,"phi1[nTracklets]/F");
  trackletTree->Branch("eta2",tdata.eta2,"eta2[nTracklets]/F");
  trackletTree->Branch("phi2",tdata.phi2,"phi2[nTracklets]/F");
  trackletTree->Branch("deta",tdata.deta,"deta[nTracklets]/F");
  trackletTree->Branch("dphi",tdata.dphi,"dphi[nTracklets]/F");

  trackletTree->Branch("npart",&tdata.npart,"npart/I");
  trackletTree->Branch("eta",tdata.eta,"eta[npart]/F");
  trackletTree->Branch("phi",tdata.phi,"phi[npart]/F");
  trackletTree->Branch("pdg",tdata.pdg,"pdg[npart]/I");
  trackletTree->Branch("chg",tdata.chg,"chg[npart]/I");
  trackletTree->Branch("nhad",tdata.nhad,"nhad[12]/F");

  // Parameters for the tree:
  Parameters par;  

  t->SetBranchAddress("eta1",par.eta1);
  t->SetBranchAddress("phi1",par.phi1);
  t->SetBranchAddress("r1",par.r1);
  t->SetBranchAddress("eta2",par.eta2);
  t->SetBranchAddress("phi2",par.phi2);
  t->SetBranchAddress("r2",par.r2);
  t->SetBranchAddress("nhits1",&par.nhits1);
  t->SetBranchAddress("nhits2",&par.nhits2);
  t->SetBranchAddress("mult",&par.mult);
  t->SetBranchAddress("vz",par.vz);
  t->SetBranchAddress("nv",&par.nv);
  t->SetBranchAddress("npart",&par.npart);
  t->SetBranchAddress("pt",par.pt);
  t->SetBranchAddress("eta",par.eta);
  t->SetBranchAddress("phi",par.phi);
  t->SetBranchAddress("chg",par.chg);
  t->SetBranchAddress("pdg",par.pdg);
  t->SetBranchAddress("gp1",par.gp1);
  t->SetBranchAddress("gp2",par.gp2);


  TH1D* hitsall = new TH1D("hitsall","Events; hits in 1st layer; events",50,0,50);
  TH1D* hitsv1 = new TH1D("hitsv1","Events with vertex found; hits in 1st layer (|#eta|<1); efficiency",50,0,50);
  TH1D* hitsv2 = new TH1D("hitsv2","Events with vertex from clusters found; hits in 1st layer (|#eta|<1); events",50,0,50);

  TH1D* hiteff1 = new TH1D("hiteff1","Vertex Finding Efficiency - Tracks; hits in 1st layer (|#eta|<1); efficiency",50,0,50);
  TH1D* hiteff2 = new TH1D("hiteff2","Vertex Finding Efficiency - Clusters; hits in 1st layer (|#eta|<1); efficiency",50,0,50);

  TH1D* hv1 = new TH1D("hv1","Vertices that are found; z; events",50,-25,25);
  TH1D* hv1no = new TH1D("hv1no","Vertices that are lost; z; events",50,-25,25);

  TH1D* hv2 = new TH1D("hv2","Cluster Vertices that are found; z; events",50,-25,25);
  TH1D* hv2no = new TH1D("hv2no","Cluster Vertices that are lost; z; events",50,-25,25);

  TH1D* zall = new TH1D("zall","vertex distribution; z[cm]; events",50,-25,25);

  TH1D* zeff1 = new TH1D("zeff1","Vertex Finding Efficiency - Tracks; z[cm]; efficiency",50,-25,25);
  TH1D* zeff2 = new TH1D("zeff2","Vertex Finding Efficiency - Clusters; z[cm]; efficiency",50,-25,25);

  TH1D* zeffv1 = new TH1D("zeffv1","Vertex Finding Efficiency - Tracks; z[cm]; efficiency",50,-25,25);


  TH2D* eff1 = new TH2D("eff1","Vertex Finding Efficiency - Tracks; z[cm]; hits in 1st layer (|#eta|<1); efficiency",50,-25,25,70,0,70);
  TH2D* eff2 = new TH2D("eff2","Vertex Finding Efficiency - Clusters; z[cm]; hits in 1st layer (|#eta|<1); efficiency",50,-25,25,70,0,70);

  TH2D* zhits = new TH2D("zhits","Vertex Finding Efficiency - Clusters; z[cm]; hits in 1st layer (|#eta|<1); efficiency",50,-25,25,70,0,70);

  TH2D* fake1 = new TH2D("fake1","Fake Vertex Rate (down to 0.5 cm) - Tracklets; z[cm]; hits in 1st layer (|#eta|<1); fake ratio",50,-25,25,70,0,70);
  TH2D* fake2 = new TH2D("fake2","Fake Vertex Rate (down to 0.5 cm) - Clusters; z[cm]; hits in 1st layer (|#eta|<1); fake ratio",50,-25,25,70,0,70);

  TH2D* mfake2 = new TH2D("mfake2","Fake Vertex Rate (down to 0.5 cm) - For no Track Vtx; z[cm]; hits in 1st layer (|#eta|<1); fake ratio",50,-25,25,70,0,70);
  TH2D* meff2 = new TH2D("meff2","Vertex Finding Efficiency - For no Track Vtx; z[cm]; hits in 1st layer (|#eta|<1); efficiency",50,-25,25,70,0,70);

  TH1D* zr2 = new TH1D("zr2","Vertex Resolution - Clusters; z_{reco}-z{sim}[cm]; events",50,-10,10);
  TH2D* zt2 = new TH2D("zt2","Vertex Scatter - Clusters; z_{sim}[cm]; z_{reco}[cm]; events",50,-25,25,50,-25,25);
  TH2D* ch2 = new TH2D("ch2","Vertex From Clusters; #chi^{2}; |z_{reco}-z_{sim}|[cm]; events",50,0,5000,50,0,50);
  TH2D* th2 = new TH2D("th2","Vertex From Clusters; Contained hits; |z_{reco}-z_{sim}|[cm]; events",50,0,150,50,0,50);

  TH2D* ct2 = new TH2D("ct2","Vertex From Clusters; Contained hits; #chi^{2}; events",50,0,50,50,0,500);
  TH2D* ct2fake = new TH2D("ct2fake","Vertex From Clusters; Contained hits; #chi^{2}; events",50,0,50,50,0,500);

  TH1D* zr1 = new TH1D("zr1","Vertex Resolution - Tracks; z_{reco}-z_{sim}[cm]; events",200,-1.,1.);

  TH1D* meff1 = new TH1D("meff1","Vertex Finding Efficiency - Track Vtx; z[cm]; hits in 1st layer (|#eta|<1); efficiency",70,0,70);


  zeff1->Sumw2();
  zall->Sumw2();

  hitsv1->Sumw2();
  hitsall->Sumw2();

  cout <<"Number of Events: "<<t->GetEntries()<<endl;

  int vtxevents = 0;
  // Main loop
  for(int i = 0; i< t->GetEntries() && i < maxEvents; ++i){    
    t->GetEntry(i);
    if (i % 1000 == 0) cout <<"Event "<<i<<endl;    
    // Selection on Events

    //    if (fabs(par.vz[0])>cuts.vzCut) continue;
    
    if (par.vz[0] == -99) continue;

    // Process the first layer
    vector<RecoHit> layer1 = removeDoubleHits(par, cuts,1);
    double l1hits = 0;
    for(int ihit = 0; ihit< (int)layer1.size(); ++ihit) {
      if(fabs(layer1[ihit].eta)<1) l1hits++;
    }

    if (fabs(par.vz[0]) < cuts.vzCut){
      hitsall->Fill(l1hits);
    }

    zall->Fill(par.vz[0]);
    zhits->Fill(par.vz[0],l1hits);

    double zres = 100; //0.5;

    if(par.vz[1] != -99){
	if (fabs(par.vz[0])< cuts.vzCut){
	  hitsv1->Fill(l1hits);
	}
        hv1->Fill(par.vz[0]);
        zeff1->Fill(par.vz[0]);
        eff1->Fill(par.vz[0],l1hits);

	zr1->Fill(par.vz[1] - par.vz[0]);
    }
 
    if(par.vz[2] != -99){
      if(fabs(par.vz[2] - par.vz[0]) < zres){
        hitsv2->Fill(l1hits);
        hiteff2->Fill(l1hits);
        hv2->Fill(par.vz[0]);
        zeff2->Fill(par.vz[0]);
        eff2->Fill(par.vz[0],l1hits);
        if(par.vz[1] == -99) meff2->Fill(par.vz[0],l1hits);

      }else{

        fake2->Fill(par.vz[0],l1hits);
        if(par.vz[1] == -99) mfake2->Fill(par.vz[0],l1hits);
      }
    }else{
      hv2no->Fill(par.vz[0]);
    }

    if(par.vz[1] == -99 && par.vz[2] == -99) continue;
    vtxevents++;
    
  }

  zr1->Sumw2();

  fake1->Divide(fake1,zhits,1,1,"B");
  fake2->Divide(fake2,zhits,1,1,"B");

  mfake2->Divide(mfake2,zhits,1,1,"B");
  meff2->Divide(meff2,zhits,1,1,"B");

  hiteff1->Divide(hitsv1,hitsall,1,1,"B");

  hiteff1->SetMaximum(1.2);

  hiteff2->Divide(hiteff2,hitsall,1,1,"B");

  eff1->Divide(eff1,zhits,1,1,"B");
  eff2->Divide(eff2,zhits,1,1,"B");


  zeffv1->Divide(zeff1,zall,1,1,"B");
  zeffv1->SetMaximum(1.2);


  TLegend * leg1 = new TLegend(0.17,0.84,0.48,0.96);
  leg1->SetFillStyle(1);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  //  leg1->AddEntry(zr2,"Hits in 1st Layer < 10","");

  TLegend * leg2 = new TLegend(0.45,0.61,0.77,0.74);
  leg2->SetFillStyle(1);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  leg2->AddEntry(hiteff1,"z_{vertex} < 10 cm","");

  TCanvas* cv1 = new TCanvas("cv1","VertexEfficiencyWithTracks",400,400);
  eff1->Draw("colz");
  cv1->Print("VtxEfficiency2D.gif");

  TCanvas* cv3 = new TCanvas("cv3","TrackingVertexResolution",400,400);
  cv3->SetLogy();
  zr1->Draw("P");

  TF1* doublegaus = new TF1("doublegaus","gaus(0)",-10,10);

  doublegaus->SetParameter(0,50000);
  doublegaus->SetParameter(2,0.2);
  doublegaus->SetParameter(1,0.);

  /*
  doublegaus->SetParameter(3,100);
  doublegaus->SetParameter(5,0.6);
  doublegaus->FixParameter(4,0.);
  */

  zr1->Fit(doublegaus);
  doublegaus->SetRange(-0.08,0.08);
  leg1->Draw();

  cv3->Print("VtxResolution.gif");


  TCanvas* cv5 = new TCanvas("cv5","Vertex Efficiency",400,400);
  hiteff1->Draw();
  leg2->Draw();

  cv5->Print("VtxEfficiencyHits.gif");


  TCanvas* cv6 = new TCanvas("cv6","Vertex Efficiency",400,400);
  zeffv1->Draw();
  //  leg2->Draw();
  cv6->Print("VtxEfficiencyZ.gif");



  outf->Write();

}

