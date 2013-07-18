#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "edwenger/TrackSpectraAnalyzer/interface/OnTrackClusterAnalyzer.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/DetId/interface/DetId.h" 
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

OnTrackClusterAnalyzer::OnTrackClusterAnalyzer(const edm::ParameterSet& iConfig) {
   //now do what ever initialization is needed
   src_ = iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("hiGlobalPrimTracks"));
   setQualityBit_ = iConfig.getUntrackedParameter<bool>("setQualityBit", false);

}

// ------------ method called to for each event  ------------
void
OnTrackClusterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
      
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker = tracker.product();

  const std::string qualityString = "highPurity";

  // get track collection    
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent.getByLabel(src_, trackCollection);
  const reco::TrackCollection *tracks = trackCollection.product();

  reco::TrackCollection::const_iterator tciter;  
  for ( tciter=tracks->begin(); tciter!=tracks->end(); tciter++) {

    if(setQualityBit_ && !tciter->quality(reco::TrackBase::qualityByName(qualityString))) continue;

    hTrackEtaPhi->Fill(tciter->eta(),tciter->phi());

    //std::cout << "\nTrack has " << tciter->recHitsSize() << " rechits" << std::endl;

    // loop over recHits on track
    for(trackingRecHit_iterator recHit = tciter->recHitsBegin();recHit!= tciter->recHitsEnd(); recHit++){
      
      //if(!(*recHit)->isValid())
      //continue;
	
      DetId detId = (*recHit)->geographicalId();
      if(!theTracker->idToDet(detId))
	continue;
      
      int layerNumber=0;
      unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());
      
      if ( subdetId ==  PixelSubdetector::PixelBarrel ) 
	{
	  PXBDetId pxbid(detId.rawId());
	  layerNumber = pxbid.layer();
	  //std::cout << "  PXB layer " << layerNumber << std::endl;
	  if((*recHit)->isValid()) {

	    if(recHit==tciter->recHitsBegin() && layerNumber!=1)
	      ntLostHits->Fill(layerNumber,pxbid.ladder(),pxbid.module());

	    GlobalPoint gpos = theTracker->idToDet((*recHit)->geographicalId())->toGlobal((*recHit)->localPosition());
	    //std::cout << "    z=" << gpos.z() << ", phi=" << gpos.phi() << std::endl;
	    if(layerNumber==1) hOnTrackClustersLayer1ZPhi->Fill(gpos.z(),gpos.phi());
	    if(layerNumber==2) hOnTrackClustersLayer2ZPhi->Fill(gpos.z(),gpos.phi());
	    if(layerNumber==3) hOnTrackClustersLayer3ZPhi->Fill(gpos.z(),gpos.phi());
	  } else {
	    //std::cout << "    INVALID HIT on layer " << layerNumber << " and ladder " << pxbid.ladder() << " and module " << pxbid.module() << std::endl;
	    ntInvalidHits->Fill(layerNumber,pxbid.ladder(),pxbid.module());
	    hInvalidHitLadder->Fill(pxbid.ladder());
	  }
	}
      /*
      else if ( subdetId ==  PixelSubdetector::PixelEndcap )
	{
	  PXFDetId pxfid(detId.rawId());
	  layerNumber = pxfid.disk();
	  std::cout << "PXF layer " << layerNumber << std::endl;
	}
      else if ( subdetId == StripSubdetector::TIB)
	{
	  TIBDetId tibid(detId.rawId());
	  layerNumber = tibid.layer();
	  std::cout << "TIB layer " << layerNumber << std::endl;
	}
      else if ( subdetId ==  StripSubdetector::TOB )
	{
	  TOBDetId tobid(detId.rawId());
	  layerNumber = tobid.layer();
	  std::cout << "TOB layer " << layerNumber << std::endl;;
	}
      else if ( subdetId ==  StripSubdetector::TID)
	{
	  TIDDetId tidid(detId.rawId());
	  layerNumber = tidid.wheel();
	  std::cout << "TID layer " << layerNumber << std::endl;
	}
      else if ( subdetId ==  StripSubdetector::TEC )
	{
	  TECDetId tecid(detId.rawId());
	  layerNumber = tecid.wheel();
	  std::cout << "TEC layer " << layerNumber << std::endl;
	}
      */
      else 
	{
	  break;
	}
      
    } // end rechit loop
  } // end track loop
	    
}
// ------------ method called once each job just before starting event loop  ------------
void 
OnTrackClusterAnalyzer::beginJob()
{

  hOnTrackClustersLayer1ZPhi = fs->make<TH2F>("hOnTrackClustersLayer1ZPhi","Layer 1, On-Track, Barrel;z [cm];#phi",100,-30.,30.,100,-3.1416,3.1416);
  hOnTrackClustersLayer2ZPhi = fs->make<TH2F>("hOnTrackClustersLayer2ZPhi","Layer 2, On-Track, Barrel;z [cm];#phi",100,-30.,30.,100,-3.1416,3.1416);
  hOnTrackClustersLayer3ZPhi = fs->make<TH2F>("hOnTrackClustersLayer3ZPhi","Layer 3, On-Track, Barrel;z [cm];#phi",100,-30.,30.,100,-3.1416,3.1416);
  hTrackEtaPhi = fs->make<TH2F>("hTrackEtaPhi","Track distribution;#eta;#phi",100,-3.,3.,100,-3.1416,3.1416);
  hInvalidHitLadder = fs->make<TH1F>("hInvalidHitLadder","Invalid Hits On Track;Ladder",200,0.,200.);
  ntInvalidHits = fs->make<TNtuple>("ntInvalidHits","invalid hits ntuple","layer:ladder:module");
  ntLostHits = fs->make<TNtuple>("ntLostHits","lost PXB1 hits ntuple","layer:ladder:module");
    
}


//define this as a plug-in
DEFINE_FWK_MODULE(OnTrackClusterAnalyzer);
