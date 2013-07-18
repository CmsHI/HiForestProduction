#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "Tracklet.h"
#include "Math/Vector3D.h"
#include "MitHigFunctions.h"

using namespace std;

void analyze_efficiency(char * infile = "/Users/yetkinyilmaz/data/hist00.root", char * outfile = "output900cont.root"){
  TFile* inf = new  TFile(infile);
  TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));

  TH1::SetDefaultSumw2();

  double etaMatch = 0.1;
  double phiMatch = 1;

  int maxEvents = 100000000;
  bool matchByEtaPhi = false;

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

  double ptmax = 10;
  int ptbins = 20;
  int etabins = 24;

  TH1D* hpt1 = new TH1D("hpt1","; p_{T} [GeV];efficiency ",ptbins,0,ptmax);
  TH1D* hpt2 = new TH1D("hpt2","Layer 2 Hit Efficiency; p_{T} [GeV];efficiency ",ptbins,0,ptmax);
  TH1D* hpt3 = new TH1D("hpt3","Two Hit Efficiency; p_{T} [GeV];efficiency ",ptbins,0,ptmax);
  TH1D* hpt4 = new TH1D("hpt4","Tracklet Efficiency; p_{T} [GeV];efficiency ",ptbins,0,ptmax);
  TH1D* hpt5 = new TH1D("hpt5","Tracklet Efficiency; p_{T} [GeV];efficiency ",ptbins,0,ptmax);
  TH1D* hpt0 = new TH1D("hpt0","Hit Efficiency; p_{T} [GeV];efficiency ",ptbins,0,ptmax);

  TH1D* heta1 = new TH1D("heta1","; #eta;efficiency ",etabins,-3,3);
  TH1D* heta2 = new TH1D("heta2","Layer 2 Hit Efficiency; #eta;efficiency ",etabins,-3,3);
  TH1D* heta3 = new TH1D("heta3","Double Hit Efficiency; #eta;efficiency ",etabins,-3,3);
  TH1D* heta4 = new TH1D("heta4","Tracklet Efficiency; #eta;efficiency ",etabins,-3,3);
  TH1D* heta5 = new TH1D("heta5","Tracklet Efficiency; #eta;efficiency ",etabins,-3,3);

  TH1D* heta0 = new TH1D("heta0","Hit Efficiency; #eta;efficiency ",etabins,-3,3);




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
      if(fabs(layer1[ihit].eta)<1) mult++;
    }

    // Process the second layer
    vector<RecoHit> layer2 = removeDoubleHits(par, cuts,2);

    for(int ihit = 0; ihit< (int)layer2.size(); ++ihit) {
      int hitbin2 = (int)layer2.size();
      if (hitbin2 > 99) hitbin2 = 99;
    }

    // Form Tracklets        
    vector<Tracklet> protoTracklets = recoProtoTracklets(layer1,layer2);
    vector<Tracklet> recoTracklets = cleanTracklets(protoTracklets,0,cuts);

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

	bool trackletMatch = false;
	bool hit1Match = false;
	bool hit2Match = false;

	for (int j3=0;j3<(int)tdata.nTracklet;j3++)
	  {
	    if(matchByEtaPhi){
	      if(fabs(tdata.eta1[j3]-par.eta[j]) > etaMatch) continue;
	      if(deltaPhi(tdata.phi1[j3],par.phi[j]) > phiMatch) continue;
	      trackletMatch = true;
	      break;
	    }else{	    
	      if(par.phi[j] == tdata.phi[j3]){
		trackletMatch = true;
		break;	
	      }
	    }
	  }
	
	for (int j1=0;j1<layer1.size();j1++)
          {
            if(matchByEtaPhi){
	      if(fabs(layer1[j1].eta-par.eta[j]) > etaMatch) continue;
	      if(deltaPhi(layer1[j1].phi,par.phi[j]) > phiMatch) continue;
	      hit1Match = true;
	    }else{
	      if(par.phi[j] == layer1[j1].genphi){
		hit1Match = true;
		break;
	      }
	    }
          }

	for (int j2=0;j2<(int)layer2.size();j2++)
          {
            if(matchByEtaPhi){
	      
	      if(fabs(layer1[j2].eta-par.eta[j]) > etaMatch) continue;
	      if(deltaPhi(layer1[j2].phi,par.phi[j]) > phiMatch) continue;
	      hit2Match = true; 
	    }else{
	      if(par.phi[j] == layer2[j2].genphi){
		hit2Match = true;
		break;
	      }
	    }
          }

	if(abs(par.eta[j])<1) hpt0->Fill(par.pt[j]);
	heta0->Fill(par.eta[j]);
	if(trackletMatch){
	  if(abs(par.eta[j])<1) hpt4->Fill(par.pt[j]);
	  heta4->Fill(par.eta[j]);
	}
	if(hit1Match){
          if(abs(par.eta[j])<1) hpt1->Fill(par.pt[j]);
          heta1->Fill(par.eta[j]);
        }
	if(hit2Match){
          if(abs(par.eta[j])<1) hpt2->Fill(par.pt[j]);
          heta2->Fill(par.eta[j]);
        }

	if(hit1Match && hit2Match){
          if(abs(par.eta[j])<1) hpt3->Fill(par.pt[j]);
          heta3->Fill(par.eta[j]);
        }

	if(hit1Match && hit2Match && trackletMatch){
          if(abs(par.eta[j])<1) hpt5->Fill(par.pt[j]);
          heta5->Fill(par.eta[j]);
        }


    }

    ntmult->Fill(mult,layer1.size(),layer2.size());
    trackletTree->Fill();
  }

  hpt1->Divide(hpt1,hpt0,1,1,"B");
  hpt2->Divide(hpt2,hpt0,1,1,"B");
  hpt3->Divide(hpt3,hpt0,1,1,"B");
  hpt4->Divide(hpt4,hpt0,1,1,"B");
  hpt5->Divide(hpt5,hpt0,1,1,"B");

  heta1->Divide(heta1,heta0,1,1,"B");
  heta2->Divide(heta2,heta0,1,1,"B");
  heta3->Divide(heta3,heta0,1,1,"B");
  heta4->Divide(heta4,heta0,1,1,"B");
  heta5->Divide(heta5,heta0,1,1,"B");

  hpt1->SetMarkerColor(4);
  hpt2->SetMarkerColor(2);
  hpt3->SetMarkerColor(6);
  hpt4->SetMarkerColor(3);
  hpt5->SetMarkerColor(1);


  heta1->SetMarkerColor(4);
  heta2->SetMarkerColor(2);
  heta3->SetMarkerColor(6);
  heta4->SetMarkerColor(3);
  heta5->SetMarkerColor(1);

  TLegend * leg1 = new TLegend(0.316,0.216,0.627,0.395);
  leg1->SetFillStyle(1);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.025);
  leg1->AddEntry(hpt1,"Hit in 1st layer","p"); 
  leg1->AddEntry(hpt2,"Hit in 2nd layer","p");
  leg1->AddEntry(hpt3,"Hit in both layers","p");
  leg1->AddEntry(hpt4,"Tracklet exists","p");
  leg1->AddEntry(hpt5,"Hit in both layers .. ","p");
  leg1->AddEntry(hpt5,".. and Tracklet exists","");

  //  leg1->AddEntry(hpt1,Form("#Delta#eta < %0.2f, #Delta#phi < %0.2f",etaMatch,phiMatch),"");

  TCanvas* ceff2 = new TCanvas("ceff2","",400,400);
  heta1->Draw();
  heta2->Draw("same");
  heta3->Draw("same");
  heta4->Draw("same");
  heta5->Draw("same");

  leg1->Draw();
  ceff2->Print("eff_eta.gif");

  leg1->AddEntry(hpt1,Form("|#eta| < 1"),"");

  TCanvas* ceff1 = new TCanvas("ceff1","",400,400);
  hpt1->Draw();
  hpt2->Draw("same");
  hpt3->Draw("same");
  hpt4->Draw("same");
  hpt5->Draw("same");

  leg1->Draw();
  ceff1->Print("eff_pt.gif");

  outf->Write();
  //  outf->Close(); 
}

