// Usage:
// .x rootlogon.C
// .x analyzeTracks.C++

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include "math.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#endif

void analyzeTracks(bool debug=false){

  // is gen track needed
  const bool isGEN = false;
  debug = true;
  const bool jet = true;
  
  // event cuts
  const unsigned int maxEvents = -1;
  const float hpFracCut = 0.25;
  const unsigned nTrackCut = 10;
  const double hfEThreshold = 3.0;
  const int nTowerThreshold = 1;
  // vertex cuts 
  const double nVtxTrackCut = 3.0;
  // track cuts
  //const char srcTrack[50] = "generalTracks";
  const char srcTrack[50] = "refitTracks";
  const string qualityString = "highPurity";
  const double normD0Cut = 999.0;
  const double normDZCut = 999.0;
  const double ptErrCut = 0.1;
  const unsigned nHitsCut = 8; // at least this many hits on track
  const double ptDebug = 3.0;  // fill debug ntuple for this selection
  // kinematic cuts
  const double etaCut = 2.5;
  // trigger names
  const int nTrigs = 5;
  char *hltNames[nTrigs] = {"HLT_MinBiasBSC","HLT_L1Jet6U","HLT_Jet15U","HLT_Jet30U","HLT_Jet50U"};

  //----- input files (900 GeV data) -----
  vector<string> fileNames;
  string fileDir = "/d101/edwenger/data7TeV/v3";       // data skim
  //string fileDir = "/d101/y_alive/mc/crab/v3";   // mc skim 53
  cout << "directory: '" << fileDir << "'" << endl;
  for(int ifile=1; ifile<=53; ifile++) {
     TString name = Form("trkAnaSkimAOD_%d_1.root",ifile);
     cout << "  adding file: " << name.Data() << endl;
     fileNames.push_back(fileDir + "/" + name.Data());
  }
  fwlite::ChainEvent event(fileNames);
  
  //----- define output hists/trees in directories of output file -----
  char fName[200];
  sprintf(fName,"output_%s_etaMax%1.1f_D0%1.1f_DZ%1.1f_pTerr%1.1f_nHits%u_nVtxTrk%1.1f_%s.root","trkAnaSkimAOD_7TeV_Jets",
	  etaCut,normD0Cut,normDZCut,ptErrCut,nHitsCut,nVtxTrackCut,srcTrack);
  TFile *outFile = new TFile(fName,"recreate");
  TH1D::SetDefaultSumw2();

  // evt hists
  outFile->cd(); outFile->mkdir("evt"); outFile->cd("evt");
  TH2D *hRunLumi    = new TH2D("hRunLumi","Event information; run number; luminosity block",500,123549.5,124049.5,200,-0.5,199.5);
  TH1D *hL1TechBits = new TH1D("hL1TechBits","L1 technical trigger bits before mask",64,-0.5,63.5);
  TH2D *hHPFracNtrk = new TH2D("hHPFracNtrk","High purity fraction vs. # of tracks; number of tracks; highPurity fraction",50,0,500,50,0,1);
  TH2D *hHfTowers   = new TH2D("hHfTowers","Number of HF tower above threshold; positive side; negative side",80,-0.5,79.5,80,-0.5,79.5);
  TH1D *hHLTPaths   = new TH1D("hHLTPaths","HLT Paths",3,0,3);
  hHLTPaths->SetBit(TH1::kCanRebin);

  // vtx hists
  outFile->cd(); outFile->mkdir("vtx"); outFile->cd("vtx");
  TH1D *hVtxSize    = new TH1D("hVtxSize","number of reconstructed pixel vertices",10,-0.5,9.5);
  TH1D *hVtxTrks    = new TH1D("hVtxTrks","number of tracks used to fit pixel vertex",50,-0.5,49.5);
  TH1D *hVtxZ       = new TH1D("hVtxZ","z position of best reconstructed pixel vertex", 80,-20,20);
  TH2D *hBeamXRun   = new TH2D("hBeamXRun","x position of beamspot; run number",500,123549.5,124049.5,80,-0.2,0.2);
  TH2D *hBeamYRun   = new TH2D("hBeamYRun","y position of beamspot; run number",500,123549.5,124049.5,80,-0.2,0.2);
  TH2D *hBeamZRun   = new TH2D("hBeamZRun","z position of beamspot; run number",500,123549.5,124049.5,80,-2,2);

  // track hists
  outFile->cd(); outFile->mkdir("trk"); outFile->cd("trk");
  TH1D *hTrkQual    = new TH1D("hTrkQual","track quality", 20, -0.5, 19.5);
  TH1D *hTrkDxyBeam = new TH1D("hTrkDxyBeam","track dxy from beamspot; dxy [cm]", 80, -2.0, 2.0);
  TH1D *hTrkDzVtx   = new TH1D("hTrkDzVtx","track dz from vertex; dz [cm]", 80, -2.0, 2.0);
  TH2D *hTrkPtErrNhits = new TH2D("hTrkPtErrNhits","track relative pt-error vs. nhits; number of valid hits on track; #sigma(p_{T})/p_{T}",30,-0.5,29.5,80,0.0,0.5);
  TH2D *hTrkPtErrEta = new TH2D("hTrkPtErrEta","track relative pt-error vs. #eta; #eta; #sigma(p_{T})/p_{T}",60,-3,3,80,0.0,0.5);
  TH1D *hTrkNhits   = new TH1D("hTrkNhits", "number of valid hits on track", 30,-0.5,29.5);
  TH1D *hTrkPt      = new TH1D("hTrkPt","track p_{T}; p_{T} [GeV/c]", 100, 0.0, 50.0);
  TH1D *hTrkEta     = new TH1D("hTrkEta","track #eta; #eta", 60, -3.0, 3.0);
  TH1D *hTrkPhi     = new TH1D("hTrkPhi","track #phi; #phi [radians]", 56, -3.5, 3.5);

  // gen track hists
  TH1D *hGenTrkPt=0;      
  if(isGEN) hGenTrkPt = new TH1D("hGenTrkPt","track p_{T}; p_{T} [GeV/c]", 100, 0.0, 50.0);

  // debug ntuple
  outFile->cd();
  TNtuple *nt=0;
  if(debug) nt = new TNtuple("nt","track debug ntuple","pt:eta:phi:qual:hits:pterr:d0:d0err:dz:dzerr:nchi2:jet6:jet15:jet30:jet50");  

  // ntuple for jet
  TNtuple *ntjet=0;
  TNtuple *ntjettrack=0;
  if(jet) ntjet = new TNtuple("ntjet","jet spectra ntuple","jet:jeta:jphi:mb:jet6:jet15:jet30:jet50"); 
  if(jet) ntjettrack = new TNtuple("ntjettrack","jet tracks correlation ntuple","ntrks:pt:eta:phi:jet:jeta:jphi:mb:jet6:jet15:jet30:jet50");

  //----- loop over events -----
  unsigned int iEvent=0;
  for(event.toBegin(); !event.atEnd(); ++event, ++iEvent){

    if( iEvent == maxEvents ) break;
    if( iEvent % 1000 == 0 ) cout << "Processing " << iEvent<< "th event: "
				  << "run " << event.id().run() 
				  << ", lumi " << event.luminosityBlock() 
				  << ", evt " << event.id().event() << endl;

    // fill event info
    hRunLumi->Fill(event.id().run(),event.luminosityBlock());

    // select on L1 trigger bits
    fwlite::Handle<L1GlobalTriggerReadoutRecord> gt;
    gt.getByLabel(event, "gtDigis");
    const TechnicalTriggerWord&  word = gt->technicalTriggerWord(); //before mask
    for(int bit=0; bit<64; bit++) hL1TechBits->Fill(bit,word.at(bit));
    if(!word.at(0)) continue;  // BPTX coincidence
    //if(!word.at(34)) continue; // BSC single-side
    if(!(word.at(40) || word.at(41))) continue; // BSC coincidence
    if(word.at(36) || word.at(37) || word.at(38) || word.at(39)) continue; // BSC halo
    
    // select on coincidence of HF towers above threshold
    fwlite::Handle<CaloTowerCollection> towers;
    towers.getByLabel(event, "towerMaker");
    int nHfTowersN=0, nHfTowersP=0;
    for(CaloTowerCollection::const_iterator calo = towers->begin(); calo != towers->end(); ++calo) {
      if(calo->energy() < hfEThreshold) continue;
      if(calo->eta()>3) nHfTowersP++;
      if(calo->eta()<-3) nHfTowersN++;
    }
    hHfTowers->Fill(nHfTowersP,nHfTowersN);
    //if(nHfTowersP < nTowerThreshold || nHfTowersN < nTowerThreshold) continue;

    // select on high-purity track fraction
    fwlite::Handle<std::vector<reco::Track> > tracks;
    tracks.getByLabel(event, srcTrack);
    int numhighpurity = 0;
    float fraction = 0;
    for(unsigned it=0; it<tracks->size(); ++it)
      if((*tracks)[it].quality(reco::TrackBase::qualityByName(qualityString))) numhighpurity++;
    if(tracks->size() > 0) fraction = (float)numhighpurity/(float)tracks->size();
    hHPFracNtrk->Fill(tracks->size(),fraction);
    if(fraction<hpFracCut && tracks->size()>nTrackCut) continue;

    // get hlt bits
    fwlite::Handle<edm::TriggerResults> triggerResults;
    triggerResults.getByLabel(event, "TriggerResults");
    const edm::TriggerNames triggerNames = event.triggerNames(*triggerResults);
    bool accept[nTrigs];
    for(int i=0; i<nTrigs; i++) {
      accept[i] = triggerResults->accept(triggerNames.triggerIndex(hltNames[i]));
      if(accept[i]) hHLTPaths->Fill(hltNames[i],1);
    }

    // select on requirement of valid vertex
    fwlite::Handle<std::vector<reco::Vertex> > vertices;
    //vertices.getByLabel(event, "pixel3Vertices");        //agglomerative pixel vertex
    vertices.getByLabel(event, "offlinePrimaryVertices");  //full-track primary vertex
    hVtxSize->Fill(vertices->size());
    if(!vertices->size()) continue;
    if(vertices->size()!=1) continue;
    size_t maxtracks=0; double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9, bestNchi2=999.9;
    int numFake=0, numVtxTrk=0;
    for(unsigned it=0; it<vertices->size(); ++it) {
      const reco::Vertex & vtx = (*vertices)[it];
      if(vtx.isFake()) numFake++;
      numVtxTrk = vtx.tracksSize();
      if(vtx.tracksSize() > maxtracks
	 || (vtx.tracksSize() == maxtracks && vtx.normalizedChi2() < bestNchi2) ) {
	maxtracks = vtx.tracksSize();
	bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
	bestNchi2 = vtx.normalizedChi2();
      }	
    }
    if(numVtxTrk<nVtxTrackCut) continue;
    if(numFake>=1) continue;
    hVtxTrks->Fill(maxtracks);
    hVtxZ->Fill(bestvz);

    // get beamspot
    fwlite::Handle<reco::BeamSpot> beamspot;
    beamspot.getByLabel(event, "offlineBeamSpot");
    hBeamXRun->Fill(event.id().run(),beamspot->x0());
    hBeamYRun->Fill(event.id().run(),beamspot->y0());
    hBeamZRun->Fill(event.id().run(),beamspot->z0());
    
    //----- loop over jets and store in a vector for a later use-----
    fwlite::Handle<std::vector<reco::CaloJet> > jets;
    jets.getByLabel(event, "iterativeCone5CaloJets");        

    std::vector<const reco::CaloJet *> sortedJets;

    for(unsigned it=0; it<jets->size(); ++it){
       const reco::CaloJet* jts = &((*jets)[it]);
       sortedJets.push_back( & *jts);
       sortByEtRef (&sortedJets);
    }

    //----- fill ntjet --------
    if(jet){
       for(unsigned it=0; sortedJets.size(); ++it){
	  ntjet->Fill(sortedJets[it]->et(),sortedJets[it]->eta(),sortedJets[it]->phi(),
		      accept[0],accept[1],accept[2],accept[3],accept[4]);
	  break;
       }
    }

    //----- loop over tracks -----
    unsigned nfinalTracks = 0;
    for(unsigned it=0; it<tracks->size(); ++it){
      
      const reco::Track & trk = (*tracks)[it];
      
      // select high purity tracks
      hTrkQual->Fill(trk.qualityMask());
      if(!trk.quality(reco::TrackBase::qualityByName(qualityString))) continue;

      // select tracks based on transverse proximity to beamspot
      double dxybeam = trk.dxy(beamspot->position());
      hTrkDxyBeam->Fill(dxybeam);
      if(fabs(dxybeam/trk.d0Error()) > normD0Cut) continue;

      // select tracks based on z-proximity to best vertex 
      math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
      double dzvtx = trk.dz(bestvtx);
      hTrkDzVtx->Fill(dzvtx);
      if(fabs(dzvtx/trk.dzError()) > normDZCut) continue;

      // select tracks based on relative pt error
      double pterr = trk.ptError()/trk.pt();
      unsigned nhits = trk.numberOfValidHits();
      hTrkPtErrNhits->Fill(nhits,pterr);
      hTrkPtErrEta->Fill(trk.eta(),pterr);
      if(debug && trk.pt() > ptDebug) // fill debug ntuple for selection of tracks
	nt->Fill(trk.pt(),trk.eta(),trk.phi(),trk.qualityMask(),
		 nhits,trk.ptError(),dxybeam,trk.d0Error(),dzvtx,trk.dzError(),trk.normalizedChi2(),
		 accept[1],accept[2],accept[3],accept[4]);
      if(pterr > ptErrCut) continue;

      // select tracks based on number of valid rechits
      hTrkNhits->Fill(nhits);
      if(nhits < nHitsCut) continue;

      // select tracks based on kinematic cuts
      if(abs(trk.eta())>etaCut) continue;
      nfinalTracks++;

      // fill selected track histograms
      hTrkPt->Fill(trk.pt());
      hTrkEta->Fill(trk.eta());
      hTrkPhi->Fill(trk.phi());

      //---loop over jet in track loop---
      if(jet){
	 for(unsigned it=0; sortedJets.size(); ++it){
	    ntjettrack->Fill(nfinalTracks,trk.pt(),trk.eta(),trk.phi(),
			     sortedJets[it]->et(),sortedJets[it]->eta(),sortedJets[it]->phi(),
			     accept[0],accept[1],accept[2],accept[3],accept[4]);       
	    break;
	 }
      }
    }
    
    //---- loop over MC gen level track ---- 
    if(!isGEN) continue;
    fwlite::Handle<std::vector<reco::GenParticle> > genTracks;
    genTracks.getByLabel(event, "genParticles");
    
    for(unsigned ip=0; ip<genTracks->size(); ++ip){
       const reco::GenParticle & p = (*genTracks)[ip];
       if(p.status() != 1) continue;
       if(p.collisionId() != 0) continue;
       if(p.charge() == 0) continue;
       if(fabs(p.eta())>etaCut) continue;
       // fill selected GEN track histograms
       hGenTrkPt->Fill(p.pt());
    }
    
  }
  
  cout << "Number of events processed : " << iEvent << endl;
  cout << "Number passing all event selection cuts: " << hVtxZ->GetEntries() << endl;

  // write to output file
  hHLTPaths->LabelsDeflate();
  outFile->Write();
  outFile->ls();
  outFile->Close();

}

