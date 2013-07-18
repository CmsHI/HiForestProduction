#ifndef SimpleTrackListMerger_h
#define SimpleTrackListMerger_h

// 
// Adapted from the SimpleTrackListMerger
// class originally written by Steve Wagner
//

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace cms
{
  class HiTrackListMerger : public edm::EDProducer
  {
  public:

    explicit HiTrackListMerger(const edm::ParameterSet& conf);

    virtual ~HiTrackListMerger();

    virtual void produce(edm::Event& e, const edm::EventSetup& c);

  private:
    edm::ParameterSet conf_;

    std::auto_ptr<reco::TrackCollection> outputTrks;
    std::auto_ptr<reco::TrackExtraCollection> outputTrkExtras;
    std::auto_ptr< TrackingRecHitCollection>  outputTrkHits;
    std::auto_ptr< std::vector<Trajectory> > outputTrajs;
    std::auto_ptr< TrajTrackAssociationCollection >  outputTTAss;
    std::auto_ptr< TrajectorySeedCollection > outputSeeds;

    reco::TrackRefProd refTrks;
    reco::TrackExtraRefProd refTrkExtras;
    TrackingRecHitRefProd refTrkHits;
    edm::RefProd< std::vector<Trajectory> > refTrajs;
    std::vector<reco::TrackRef> trackRefs;
    edm::RefProd< TrajectorySeedCollection > refTrajSeeds;
    
    bool copyExtras_;
    bool makeReKeyedSeeds_;
  };
}


#endif

