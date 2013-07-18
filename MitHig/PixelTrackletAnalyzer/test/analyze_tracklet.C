#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
//#include "RecoHit.h"
#include "Tracklet.h"
#include "Math/Vector3D.h"

using namespace std;

void analyze_tracklet(char * infile){
  TFile* inf = new  TFile(infile);
  TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));

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
  TFile* outf = new TFile("output.root","recreate");
  TNtuple *ntmatched = new TNtuple("ntmatched","","eta1:matchedeta:phi1:matchedphi:deta:dphi:signalCheck:tid:r1id:r2id:evtid:nhit1:sid:ptype:vz");
  TNtuple *ntmult = new TNtuple("ntmult","","mult:nhit1:nhit2");

  vector<TH1D*> layer1HitEta;
  layer1HitEta.reserve(nbins);
  vector<TH1D*> layer1HitPhi;
  layer1HitPhi.reserve(nbins);

  vector<TH1D*> layer2HitEta;
  layer2HitEta.reserve(nbins);
  vector<TH1D*> layer2HitPhi;
  layer2HitPhi.reserve(nbins);
    
  for(int i = 0; i< nbins; ++i){
    layer1HitEta[i] = new TH1D(Form("dNdEtaHits1_%02d",i),"dNdEta Hits Layer 1",500,-3,3);
    layer2HitEta[i] = new TH1D(Form("dNdEtaHits2_%02d",i),"dNdEta Hits Layer 2",500,-3,3);
    layer1HitPhi[i] = new TH1D(Form("dNdPhiHits1_%02d",i),"dNdPhi Hits Layer 1",500,-3.2,3.2);
    layer2HitPhi[i] = new TH1D(Form("dNdPhiHits2_%02d",i),"dNdPhi Hits Layer 2",500,-3.2,3.2);
  }
  
  TH1D* hm1 = new TH1D("hm1","Number of Hits Layer 1",50,0,50);

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

  cout <<"Number of Events: "<<t->GetEntries()<<endl;

  // Main loop
  for(int i = 0; i< t->GetEntries(); ++i){    
    t->GetEntry(i);
    if (i % 1000 == 0) cout <<"Event "<<i<<endl;    
    // Selection on Events
    if (fabs(par.vz[1])>cuts.vzCut) continue;

    hm1->Fill(par.mult);

    // Process the first layer
    vector<RecoHit> layer1 = removeDoubleHits(par, cuts,1);
    double mult = 0;
    for(int ihit = 0; ihit< (int)layer1.size(); ++ihit) {
      int hitbin1 = (int)layer1.size();
      if (hitbin1 > 99) hitbin1 = 99;
      layer1HitEta[hitbin1]->Fill(layer1[ihit].eta);
      layer1HitPhi[hitbin1]->Fill(layer1[ihit].phi);
      if(fabs(layer1[ihit].eta)<1) mult++;
    }

    // Process the second layer
    vector<RecoHit> layer2 = removeDoubleHits(par, cuts,2);

    for(int ihit = 0; ihit< (int)layer2.size(); ++ihit) {
      int hitbin2 = (int)layer2.size();
      if (hitbin2 > 99) hitbin2 = 99;
      layer2HitEta[hitbin2]->Fill(layer2[ihit].eta);
      layer2HitPhi[hitbin2]->Fill(layer2[ihit].phi);
    }

    // Form Tracklets        
    vector<Tracklet> protoTracklets = recoProtoTracklets(layer1,layer2);
    vector<Tracklet> recoTracklets = cleanTracklets(protoTracklets,0,cuts);

    // Fill Ntuple
    for (int j=0;j<(int)recoTracklets.size();j++)
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
	var[11] = (int)par.mult;
	var[12] = recoTracklets[j].getSId();
	var[13] = recoTracklets[j].getType();
	var[14] = par.vz[1];
        ntmatched->Fill(var);
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());
  }

  outf->Write();
  outf->Close(); 
}

