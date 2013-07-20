// Original Author:  Edward Allen Wenger,32 4-A06,+41227676980,
//         Created:  Fri May  7 13:11:39 CEST 2010
// $Id: HiVertexAnalyzer.cc,v 1.6 2011/03/21 16:22:02 sungho Exp $
//

#include "edwenger/HiVertexAnalyzer/interface/HiVertexAnalyzer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "edwenger/HiVertexAnalyzer/interface/HiVertexComparator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

HiVertexAnalyzer::HiVertexAnalyzer(const edm::ParameterSet& iConfig)
:
  vtxlabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxlabel")),
  nTracksBins_(iConfig.getUntrackedParameter<int>("nTracksBins")),
  neededCentBins_(iConfig.getUntrackedParameter<std::vector<int> >("neededCentBins")),
  centrality_(0)
{

}


// ------------ method called to for each event  ------------
void
HiVertexAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  //------ Centrality -------------------
  centrality_ = new CentralityProvider(iSetup);
  centrality_->newEvent(iEvent,iSetup);
  int cbin = centrality_->getBin();


  //------- Vertices ---------------------
  edm::Handle<reco::VertexCollection> vtxsH;
  iEvent.getByLabel(vtxlabel_,vtxsH);

  reco::VertexCollection vtxs = *vtxsH;
  std::sort(vtxs.begin(),vtxs.end(),MoreTracksThenLowerChi2<reco::Vertex>());

  hVtxSize->Fill(vtxs.size());

  if(!vtxs.size()) return;

  hVtxTracks->Fill(vtxs[0].tracksSize());

  hVtxZ->Fill(vtxs[0].z());
  hVtxX->Fill(vtxs[0].x());
  hVtxY->Fill(vtxs[0].y());

  hVtxXnY->Fill(vtxs[0].x(),vtxs[0].y());
  hVtxXnZ->Fill(vtxs[0].x(),vtxs[0].z());
  hVtxYnZ->Fill(vtxs[0].y(),vtxs[0].z());

  for(unsigned i=0;i<neededCentBins_.size()-1;i++){
     if(i==0){
        if(cbin<=neededCentBins_[i+1]){
           hVtxTracks_Cent[i]->Fill(vtxs[0].tracksSize());
	   hVtxZ_Cent[i]->Fill(vtxs[0].z());
	   hVtxX_Cent[i]->Fill(vtxs[0].x());
	   hVtxY_Cent[i]->Fill(vtxs[0].y());
	}
     }else{
        if(cbin>neededCentBins_[i] && cbin<=neededCentBins_[i+1]){
           hVtxTracks_Cent[i]->Fill(vtxs[0].tracksSize());
	   hVtxZ_Cent[i]->Fill(vtxs[0].z());
           hVtxX_Cent[i]->Fill(vtxs[0].x());
           hVtxY_Cent[i]->Fill(vtxs[0].y());
	}
     }
  }


  for(unsigned i=0; i<vtxs.size(); i++)
     hAllVtxZ->Fill(vtxs[i].z());

  if(vtxs.size()<2) return;

  hMultVtxTracks->Fill(vtxs[0].tracksSize(),vtxs[1].tracksSize());
  hMultVtxZ->Fill(vtxs[0].z(),vtxs[1].z());

  for(unsigned i=0; i<vtxs.size(); i++)
    edm::LogVerbatim("HiVertexAnalyzer") << "vtx #" << i << " has " 
		     << vtxs[i].tracksSize() << " tracks";
}


// ------------ method called once each job just before starting event loop  ------------
void 
HiVertexAnalyzer::beginJob()
{

  hVtxSize = f->make<TH1F>("hVtxSize","number of reconstructed vertices",10,-0.5,9.5);
  hVtxTracks = f->make<TH1F>("hVtxTracks","number of tracks fitted to vertex",nTracksBins_,-0.5,(float)nTracksBins_-0.5);

  hVtxZ = f->make<TH1F>("hVtxZ","z position of best reconstructed vertex",120,-30.0,30.0);
  hVtxX = f->make<TH1F>("hVtxX","x position of best reconstructed vertex",120,-10.0,10.0);
  hVtxY = f->make<TH1F>("hVtxY","y position of best reconstructed vertex",120,-10.0,10.0);

  hVtxXnY = f->make<TH2F>("hVtxXnY","x vs y position", 120,-10.0,10.0, 120,-10.0,10.0);
  hVtxXnZ = f->make<TH2F>("hVtxXnZ","x vs z position", 120,-10.0,10.0, 120,-30.0,30.0);
  hVtxYnZ = f->make<TH2F>("hVtxYnZ","y vs z position", 120,-10.0,10.0, 120,-30.0,30.0);

  hMultVtxTracks = f->make<TH2F>("hMultVtxTracks","number of tracks fitted to vertex; most populated; 2nd most populated",nTracksBins_,-0.5,(float)nTracksBins_-0.5,nTracksBins_,-0.5,(float)nTracksBins_-0.5);
  hMultVtxZ = f->make<TH2F>("hMultVtxZ","z position of reconstructed vertex; most populated; 2nd most populated",120,-30.0,30.0,120,-30.0,30.0);
  hAllVtxZ = f->make<TH1F>("hAllVtxZ","z position of all reconstructed vertices",120,-30.0,30.0);

  // centrality binned histogram
  for(unsigned i=0;i<neededCentBins_.size()-1;i++){

     hVtxZ_Cent.push_back(f->make<TH1F>("","z position of best reconstructed vertex", 120,-30.0,30.0));
     hVtxX_Cent.push_back(f->make<TH1F>("","x position of best reconstructed vertex", 120,-10.0,10.0));
     hVtxY_Cent.push_back(f->make<TH1F>("","y position of best reconstructed vertex", 120,-10.0,10.0));

     hVtxTracks_Cent.push_back(f->make<TH1F>("","number of tracks fitted to vertex",nTracksBins_,-0.5,(float)nTracksBins_-0.5));

     if(i==0){
	hVtxZ_Cent[i]->SetName(Form("hVtxZ_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	hVtxY_Cent[i]->SetName(Form("hVtxY_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	hVtxX_Cent[i]->SetName(Form("hVtxX_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
	hVtxTracks_Cent[i]->SetName(Form("hVtxTracks_cbin%dto%d",neededCentBins_[i],neededCentBins_[i+1]));
     }else{
	hVtxZ_Cent[i]->SetName(Form("hVtxZ_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	hVtxY_Cent[i]->SetName(Form("hVtxY_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
	hVtxX_Cent[i]->SetName(Form("hVtxX_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
        hVtxTracks_Cent[i]->SetName(Form("hVtxTracks_cbin%dto%d",neededCentBins_[i]+1,neededCentBins_[i+1]));
     }

  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiVertexAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiVertexAnalyzer);
