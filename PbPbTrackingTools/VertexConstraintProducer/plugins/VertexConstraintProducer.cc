// -*- C++ -*-
//
// Package:    VertexConstraintProducer
// Class:      VertexConstraintProducer
// 
/**\class VertexConstraintProducer VertexConstraintProducer.cc RecoTracker/ConstraintProducerTest/src/VertexConstraintProducer.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Giuseppe Cerati
//         Created:  Tue Jul 10 15:05:02 CEST 2007
// $Id: VertexConstraintProducer.cc,v 1.2 2010/07/10 21:41:54 sungho Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/TrackConstraintAssociation.h"

// vertex 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// beamspot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"



//
// class decleration
//

class VertexConstraintProducer: public edm::EDProducer {
public:
  explicit VertexConstraintProducer(const edm::ParameterSet&);
  ~VertexConstraintProducer();

private:
  virtual void beginRun(edm::Run & run, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
 const edm::ParameterSet iConfig_;
 
 edm::InputTag tracksrc_;
 edm::InputTag vtxsrc_;
   
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VertexConstraintProducer::VertexConstraintProducer(const edm::ParameterSet& iConfig) : iConfig_(iConfig)
{
  //register your products
  produces<std::vector<VertexConstraint> >();
  produces<TrackVtxConstraintAssociationCollection>();

  //now do what ever other initialization is needed
  tracksrc_ = iConfig.getUntrackedParameter<edm::InputTag>("tracksrc",edm::InputTag("hiSelectedTracks"));
  vtxsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("vtxsrc",edm::InputTag("hiSelectedVertex"));
}


VertexConstraintProducer::~VertexConstraintProducer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void VertexConstraintProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<reco::TrackCollection> theTCollection;
  iEvent.getByLabel(tracksrc_,theTCollection); 
  
  std::auto_ptr<std::vector<VertexConstraint> > pairs(new std::vector<VertexConstraint>);
  std::auto_ptr<TrackVtxConstraintAssociationCollection> output(new TrackVtxConstraintAssociationCollection);
  
  edm::RefProd<std::vector<VertexConstraint> > rPairs = iEvent.getRefBeforePut<std::vector<VertexConstraint> >();

  edm::Handle<reco::VertexCollection> theVCollection;
  iEvent.getByLabel(vtxsrc_,theVCollection);
  
  const reco::Vertex * validVtx = 0;
  
  //bool  validvertex = false;
  for(reco::VertexCollection::const_iterator vtx=theVCollection->begin(); vtx!=theVCollection->end();vtx++){
     // let's keep it simple right now
     validVtx = &(*vtx);
     break;
     /*
     if(vtx->isValid() && (!vtx->isFake())){
	validVtx = &(*vtx);
	validvertex = true;
	break;
     }
     */
  }
  //if(!validvertex) std::cout<<"No real valid vertex found"<<std::endl;

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByLabel("offlineBeamSpot", theBeamSpot);
  const reco::BeamSpot * theBS;
  theBS = theBeamSpot.product();

  //float beamSize = theBS->BeamWidthX();

  int index = 0;
  for (reco::TrackCollection::const_iterator i=theTCollection->begin(); i!=theTCollection->end();i++) {
     VertexConstraint tmp(GlobalPoint(validVtx->position().x(),
				      validVtx->position().y(),
				      validVtx->position().z()),
			  GlobalError(validVtx->covariance(0,0), 0,
				                                 validVtx->covariance(1,1), 0,
				                                 0,validVtx->covariance(2,2))); 
			  //GlobalError(beamSize*beamSize, 0, 
			  //beamSize*beamSize, 0,
			  //0,validVtx->covariance(2,2)));
    pairs->push_back(tmp);
    output->insert(reco::TrackRef(theTCollection,index), edm::Ref<std::vector<VertexConstraint> >(rPairs,index) );
    index++;
  }
  
  iEvent.put(pairs);
  iEvent.put(output);
}

// ------------ method called once each job just before starting event loop  ------------
void VertexConstraintProducer::beginRun(edm::Run & run, const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void VertexConstraintProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexConstraintProducer);
