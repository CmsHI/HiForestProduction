// Usage:
// .x rootlogon.C
// .x analyzeTracks.C++

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <fstream>

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

void makeVtxNtuple(){
  
  // event cuts
  const unsigned int maxEvents = -1;
 
  //----- input files (7 TeV data) -----
  vector<string> fileNames;
  //string fileDir = "/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/edwenger/trkAnaSkim/MB_C10-Apr1ReReco_v2_pptrkana_skim/v3";       // data skim
  string fileDir = "/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/edwenger/trkAnaSkim/MB_C10-PromptReco_v8_pptrkana_skim/v4";       // data skim
  //string fileDir = "/d101/y_alive/mc/crab/v3";       // mc skim 53
  cout << "directory: '" << fileDir << "'" << endl;


  ifstream inf("filelist.txt");
  string word;
  while (inf >> word)
    {
      fileNames.push_back( "dcache:" + fileDir + "/" + word );
    }

  fwlite::ChainEvent event(fileNames);
  
  //----- define output hists/trees in directories of output file -----
  TFile *outFile = new TFile("vtxNtuple.root","recreate");
  TH1D::SetDefaultSumw2();

  // vtx hists
  outFile->cd(); outFile->mkdir("vtx"); outFile->cd("vtx");
  TH1D *hVtxSize    = new TH1D("hVtxSize","number of reconstructed pixel vertices",10,-0.5,9.5);
  TH1D *hVtxTrks    = new TH1D("hVtxTrks","number of tracks used to fit pixel vertex",50,-0.5,49.5);
  TH1D *hVtxZ       = new TH1D("hVtxZ","z position of best reconstructed pixel vertex", 80,-20,20);

  // debug ntuple
  outFile->cd();
  TNtuple *nt = new TNtuple("nt","vtx debug ntuple","run:lumi:npvtx:nvtx:nptrk:pvx:pvy:pvz:nptrk2:pvz2:ntrk:vx:vy:vz");  

  //----- loop over events -----
  unsigned int iEvent=0;
  for(event.toBegin(); !event.atEnd(); ++event, ++iEvent){

    if( iEvent == maxEvents ) break;
    if( iEvent % 1000 == 0 ) cout << "Processing " << iEvent<< "th event: "
				  << "run " << event.id().run() 
				  << ", lumi " << event.luminosityBlock() 
				  << ", evt " << event.id().event() << endl;

    //------------Vertex loop------------------
    fwlite::Handle<std::vector<reco::Vertex> > pxlvertices;
    pxlvertices.getByLabel(event, "pixel3Vertices");        //agglomerative pixel vertex
    
    hVtxSize->Fill(pxlvertices->size());
    
    if(!pxlvertices->size()) continue;
    
    size_t maxtracks=0; double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9, bestNchi2=999.9;
    size_t nexttracks=0; double nextvz=-999.9, nextNchi2=999.9;
 
    for(unsigned it=0; it<pxlvertices->size(); ++it) {
      const reco::Vertex & vtx = (*pxlvertices)[it];
      // if this is the new best vertex...
      if(vtx.tracksSize() > maxtracks
	 || (vtx.tracksSize() == maxtracks && vtx.normalizedChi2() < bestNchi2) ) {
	// move previous best to next best
	nexttracks = maxtracks; nextvz = bestvz; nextNchi2 = bestNchi2;
	// and update new best vertex
	maxtracks = vtx.tracksSize();
	bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
	bestNchi2 = vtx.normalizedChi2();
      }	
      // or if this is the new next-best vertex
      else if (vtx.tracksSize() > nexttracks
	       || (vtx.tracksSize() == nexttracks && vtx.normalizedChi2() < nextNchi2)) {
	// update new next vertex
	nexttracks = vtx.tracksSize();
	nextvz = vtx.z();
	nextNchi2 = vtx.normalizedChi2();
      }
    }

    hVtxTrks->Fill(maxtracks);
    hVtxZ->Fill(bestvz);

    // offline primary vertices
    fwlite::Handle<std::vector<reco::Vertex> > vertices;
    vertices.getByLabel(event, "offlinePrimaryVertices");  //full-track primary vertex

    size_t maxofflinetracks=0;
    double offlineVx = -999.9, offlineVy = -999.9, offlineVz = -999.9, offlineNchi2=999.9;

    for(unsigned it=0; it<vertices->size(); ++it) {
      offlineVx = -999.9; offlineVy = -999.9; offlineVz = -999.9;
      const reco::Vertex & vtx2 = (*vertices)[it];
      if(!vtx2.isFake() && (vtx2.tracksSize()>maxofflinetracks || (vtx2.tracksSize() == maxofflinetracks && vtx2.normalizedChi2() < offlineNchi2)) )
	offlineVx = vtx2.x(); offlineVy = vtx2.y(); offlineVz = vtx2.z(); maxofflinetracks = vtx2.tracksSize();
    
    }
  
    //-----------Beamspot-------------------
    fwlite::Handle<reco::BeamSpot> beamspot;
    beamspot.getByLabel(event, "offlineBeamSpot");
    
    // Fill run, lumi, vtx, beamspot ntuple
    nt->Fill(event.id().run(),
	     event.luminosityBlock(),
	     pxlvertices->size(),
	     vertices->size(),
	     maxtracks,
	     bestvx,
	     bestvy,
	     bestvz,
	     nexttracks,
	     nextvz,
	     maxofflinetracks,
	     offlineVx,
	     offlineVy,
	     offlineVz//,
	     //beamspot->x0(),
	     //beamspot->y0(),
	     //beamspot->z0(),
	     //beamspot->sigmaZ()
	     );

  }
  
  cout << "Number of events processed : " << iEvent << endl;
  cout << "Number passing all event selection cuts: " << hVtxZ->GetEntries() << endl;

  // write to output file
  outFile->Write();
  outFile->ls();
  outFile->Close();

}

