#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "Tracklet.h"
#include "Math/Vector3D.h"

using namespace std;

void analyze_trackletTree(char * infile, char * outfile = "output.root"){
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
  trackletTree->Branch("deta",tdata.deta,"data[nTracklets]/F");
  trackletTree->Branch("dphi",tdata.dphi,"data[nTracklets]/F");
  

  trackletTree->Branch("npart",&tdata.npart,"npart/I");
  trackletTree->Branch("eta",tdata.eta,"eta[npart]/F");
  trackletTree->Branch("phi",tdata.phi,"phi[npart]/F");
  trackletTree->Branch("pdg",tdata.pdg,"pdg[npart]/I");
  trackletTree->Branch("chg",tdata.chg,"chg[npart]/I");
  trackletTree->Branch("nhad",tdata.nhad,"nhad[12]/F");


  vector<TH1D*> layer1HitEta;
  layer1HitEta.reserve(nbins);
  vector<TH1D*> layer1HitPhi;
  layer1HitPhi.reserve(nbins);

  vector<TH1D*> layer2HitEta;
  layer2HitEta.reserve(nbins);
  vector<TH1D*> layer2HitPhi;
  layer2HitPhi.reserve(nbins);

  vector<RecoHit> layer;
    
  for(int i = 0; i< nbins; ++i){
    layer1HitEta[i] = new TH1D(Form("dNdEtaHits1_%02d",i),"dNdEta Hits Layer 1",500,-3,3);
    layer2HitEta[i] = new TH1D(Form("dNdEtaHits2_%02d",i),"dNdEta Hits Layer 2",500,-3,3);
    layer1HitPhi[i] = new TH1D(Form("dNdPhiHits1_%02d",i),"dNdPhi Hits Layer 1",500,-3.2,3.2);
    layer2HitPhi[i] = new TH1D(Form("dNdPhiHits2_%02d",i),"dNdPhi Hits Layer 2",500,-3.2,3.2);
  }

  TH3F* nhits = new TH3F("nhits","",100,0,100,100,0,100,100,0,100);
  
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
  t->SetBranchAddress("eta",&par.eta);
  t->SetBranchAddress("phi",&par.phi);
  t->SetBranchAddress("chg",&par.chg);
  t->SetBranchAddress("pdg",&par.pdg);
  t->SetBranchAddress("gp1",par.gp1);
  t->SetBranchAddress("gp2",par.gp2);

  cout <<"Number of Events: "<<t->GetEntries()<<endl;

  // Main loop
  for(int i = 0; i< t->GetEntries() && i < maxEvents; ++i){    
    t->GetEntry(i);
    if (i % 1000 == 0) cout <<"Event "<<i<<endl;    
    // Selection on Events
    if (fabs(par.vz[1])>cuts.vzCut) continue;
    
    // Fill vertex information
    tdata.nv = par.nv;
    for(int j = 0; j<par.nv;j++) {
       tdata.vz[j] = par.vz[j];
    }
    
    // Process the first layer
    vector<RecoHit> layer1 = removeDoubleHits(par, cuts,1);
    double mult = 0;
    for(int ihit = 0; ihit< (int)layer1.size(); ++ihit) {
      int hitbin1 = (int)layer1.size();
      if (hitbin1 > 99) hitbin1 = 99;
      layer1HitEta[hitbin1]->Fill(layer1[ihit].eta);
      layer1HitPhi[hitbin1]->Fill(layer1[ihit].phi);
      RecoHit myhit(layer1[ihit].eta,layer1[ihit].phi,layer1[ihit].r,1);
      myhit.geneta = layer1[ihit].geneta;
      myhit.genphi = layer1[ihit].genphi;
      myhit.genpt = layer1[ihit].genpt;
      myhit.genid = layer1[ihit].genid;
      layer.push_back(myhit);
      if(fabs(layer1[ihit].eta)<1) mult++;
    }

    // Process the second layer
    vector<RecoHit> layer2 = removeDoubleHits(par, cuts,2);

    for(int ihit = 0; ihit< (int)layer2.size(); ++ihit) {
      int hitbin2 = (int)layer2.size();
      if (hitbin2 > 99) hitbin2 = 99;
      layer2HitEta[hitbin2]->Fill(layer2[ihit].eta);
      layer2HitPhi[hitbin2]->Fill(layer2[ihit].phi);
      RecoHit myhit(layer2[ihit].eta,layer2[ihit].phi,layer2[ihit].r,2);
      myhit.geneta = layer2[ihit].geneta;
      myhit.genphi = layer2[ihit].genphi;
      myhit.genpt = layer2[ihit].genpt;
      myhit.genid = layer2[ihit].genid;
      layer.push_back(myhit);
    }

    // Form Tracklets        
//    vector<Tracklet> protoTracklets = recoProtoTracklets(layer1,layer2);
//    vector<Tracklet> recoTracklets = cleanTracklets(protoTracklets,0,cuts);
    vector<Tracklet> recoTracklets = recoFastTracklets(layer);

    // Fill Ntuple
    tdata.nTracklet = recoTracklets.size();
    tdata.nhit1 = layer1.size();
    tdata.nhit2 = layer2.size();
    
    for (int j=0;j<(int)tdata.nTracklet;j++)
    {
        float var[100];
	var[0] = recoTracklets[j].eta1();
	var[1] = recoTracklets[j].eta2();
	var[2] = recoTracklets[j].phi1();
	var[3] = recoTracklets[j].phi2();
	var[4] = recoTracklets[j].deta();
	var[5] = recoTracklets[j].dphi();
	var[6] = 0;
	var[7] = recoTracklets[j].getId();
	var[8] = recoTracklets[j].getId1();
	var[9] = recoTracklets[j].getId2();
	var[10] = i;
	var[11] = (int)mult;
	var[12] = recoTracklets[j].getSId();
	var[13] = recoTracklets[j].getType();
	var[14] = par.vz[1];
	var[15] = (float)recoTracklets[j].eta;
	var[16] = (float)recoTracklets[j].phi;
	var[17] = (float)recoTracklets[j].pt;
	var[18] = (float)recoTracklets[j].genid;

        ntmatched->Fill(var);

        tdata.eta1[j] = recoTracklets[j].eta1();	
        tdata.eta2[j] = recoTracklets[j].eta2();	
        tdata.phi1[j] = recoTracklets[j].phi1();	
        tdata.phi2[j] = recoTracklets[j].phi2();
	tdata.deta[j] = recoTracklets[j].deta();
   	tdata.dphi[j] = recoTracklets[j].dphi();
	tdata.mult = (int)mult;	
    }

    tdata.npart=0;
    for (int j=0;j<12;j++) tdata.nhad[j]=0;

    for(int j=0;j<par.npart;j++)
    {
        if (fabs(par.eta[j])>3||par.chg[j]==0) continue;
	tdata.eta[tdata.npart]=par.eta[j];
	tdata.phi[tdata.npart]=par.phi[j];
	tdata.chg[tdata.npart]=par.chg[j];
	tdata.pdg[tdata.npart]=par.pdg[j];
        tdata.npart++;
	int bin = (int)((par.eta[j]+3)*2);
	int pdg = (int)abs(par.pdg[j]);
	if (pdg== 211 || pdg == 321 || pdg == 2212 || pdg ==3122) tdata.nhad[bin]++;
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());
    trackletTree->Fill();
    layer1.clear();
    layer2.clear();
    layer.clear();
    recoTracklets.clear();    
  }

  outf->Write();
  outf->Close(); 
}

