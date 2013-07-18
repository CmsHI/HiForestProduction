//
// Original Author:  Edward Wenger
//         Created:  Thu Apr 29 14:31:47 CEST 2010
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Appeltel/PixelTrackAnalysis/interface/HiPixelTrkEffHistograms.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"


// define track efficiency analyzer class
class HiPixelTrkEffAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HiPixelTrkEffAnalyzer(const edm::ParameterSet&);
      ~HiPixelTrkEffAnalyzer(){};


   private:
      virtual void beginJob();
      virtual void beginRun(const edm::Run&, const edm::EventSetup&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      
      SimTrack_t setSimTrack(TrackingParticle&, const reco::Track&, size_t, float);
      RecTrack_t setRecTrack(reco::Track&, const TrackingParticle&, size_t, float);
      bool testVertex(reco::Track&, double&, double&, double&, double&);
      std::pair<bool,bool> isAccepted(TrackingParticle&);
      int getLayerId(const PSimHit&);

      // ----------member data ---------------------------
      edm::InputTag trackTags_; 
      edm::InputTag jetTags_;
      edm::InputTag label_tp_effic_;
      edm::InputTag label_tp_fake_;
      edm::InputTag associatorMap_;
      edm::InputTag vtxTags_;
      edm::InputTag bsTags_;
  
      const TrackerGeometry * theTracker;
      edm::Handle<reco::VertexCollection> vertexCollectionH;
      edm::Handle<reco::BeamSpot> beamSpotH;
      bool doAssociation_;
      bool hasSimInfo_;

      CentralityProvider * centrality_;
     
      HiPixelTrkEffHistograms *histograms;
      edm::Service<TFileService> f;

      enum { BPix1=0, BPix2=1, BPix3=2,
	     FPix1_neg=3, FPix2_neg=4,
	     FPix1_pos=5, FPix2_pos=6,
	     nLayers=7};

};

