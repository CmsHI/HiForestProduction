// -*- C++ -*-
//
// Package:    NjetFilter
// Class:      NjetFilter
// 
/**\class NjetFilter NjetFilter.cc CmsHi/NjetFilter/src/NjetFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yetkin Yilmaz
//         Created:  Fri Sep 10 10:02:12 EDT 2010
// $Id: NjetFilter.cc,v 1.2 2010/10/22 13:34:09 yilmaz Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

//
// class declaration
//

class NjetFilter : public edm::EDFilter {
   public:
      explicit NjetFilter(const edm::ParameterSet&);
      ~NjetFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

   edm::InputTag jetTag_;
   edm::Handle<reco::JetView> jets;

   unsigned int njet_;
   double ptMin_;
   bool equalMode_;
   bool eqGreaterMode_;
   bool lessMode_;
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
NjetFilter::NjetFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   jetTag_ = iConfig.getParameter<edm::InputTag>("src");
   njet_ = iConfig.getParameter<unsigned int>("nJet");
   ptMin_ = iConfig.getParameter<double>("ptMin");
   equalMode_ = iConfig.getParameter<bool>("equalMode");
   eqGreaterMode_ = iConfig.getParameter<bool>("eqGreaterMode");
   lessMode_ = iConfig.getParameter<bool>("lessMode");

}


NjetFilter::~NjetFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
NjetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   iEvent.getByLabel(jetTag_,jets);
   if(equalMode_ && jets->size() == njet_) return true;
   if(eqGreaterMode_ && jets->size() >= njet_) return true;
   if(lessMode_ && jets->size() < njet_) return true;

   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
NjetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NjetFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(NjetFilter);
