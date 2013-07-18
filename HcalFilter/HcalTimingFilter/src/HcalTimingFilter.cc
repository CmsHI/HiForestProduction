// -*- C++ -*-
//
// Package:    HcalTimingFilter
// Class:      HcalTimingFilter
// 
/**\class HcalTimingFilter HcalTimingFilter.cc HcalFilter/HcalTimingFilter/src/HcalTimingFilter.cc

 Description: Implement the HcalTiming Filter on an event level (Flagging implements it on a RecHit level)

 Needs: 
addpkg RecoLocalCalo/HcalRecAlgos
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/src/HBHETimingShapedFlag.cc
cvs co -r V00-07-21 RecoLocalCalo/HcalRecAlgos/interface/HBHETimingShapedFlag.h
cvs co -r V00-00-18 JetMETAnalysis/HcalReflagging

 Implementation (cfg.py):
process.load("JetMETAnalysis.HcalReflagging.hbherechitreflaggerJETMET_cfi")

process.hbheReflagNewTimeEnv = process.hbherechitreflaggerJETMET.clone()
process.hbheReflagNewTimeEnv.timingshapedcutsParameters.hbheTimingFlagBit=cms.untracked.int32(8)
process.hcalTimingFilter = cms.EDFilter("HcalTimingFilter",
   hbheHits = cms.untracked.InputTag("hbheReflagNewTimeEnv")
)
process.HFilt = cms.Path(process.hcalTimingFilter*process.output)
# where process.output is an example of the second thing you are doing, the filter should
# come first and it will remove events that are flagged by the Hcal timing filter.

*/
//
// Original Author:  Yen-Jie Lee
//         Created:  Fri Nov 26 07:02:44 EST 2010
// $Id: HcalTimingFilter.cc,v 1.1 2010/12/06 16:21:43 belt Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

//
// class declaration
//

class HcalTimingFilter : public edm::EDFilter {
   public:
      explicit HcalTimingFilter(const edm::ParameterSet&);
      ~HcalTimingFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

      edm::InputTag hbheHits;
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
HcalTimingFilter::HcalTimingFilter(const edm::ParameterSet& iConfig)
{
   hbheHits = iConfig.getUntrackedParameter<edm::InputTag>("hbheHits");

   //now do what ever initialization is needed

}


HcalTimingFilter::~HcalTimingFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HcalTimingFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // get stuff
   Handle<HBHERecHitCollection> hRecHits;
   iEvent.getByLabel(hbheHits, hRecHits);

   for(int i = 0; i < (int)hRecHits->size(); i++)
   {
      HBHERecHit rechit = (*hRecHits)[i];
      if (rechit.flagField(8)) return false;
   }

   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HcalTimingFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HcalTimingFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalTimingFilter);
