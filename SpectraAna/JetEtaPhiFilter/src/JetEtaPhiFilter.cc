// -*- C++ -*-
//
// Package:    JetEtaPhiFilter
// Class:      JetEtaPhiFilter
// 
/**\class JetEtaPhiFilter JetEtaPhiFilter.cc SpectraAna/JetEtaPhiFilter/src/JetEtaPhiFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andre Yoon,32 4-A06,+41227676980,
//         Created:  Wed Jul  6 11:36:18 CEST 2011
// $Id: JetEtaPhiFilter.cc,v 1.1 2011/07/06 11:35:57 sungho Exp $
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

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"

//
// class declaration
//

class JetEtaPhiFilter : public edm::EDFilter {
   public:
      explicit JetEtaPhiFilter(const edm::ParameterSet&);
      ~JetEtaPhiFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag jsrc_;

      double minEta_, maxEta_;   
      double minPhi_, maxPhi_;

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
JetEtaPhiFilter::JetEtaPhiFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   jsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("jsrc",edm::InputTag("iterativeCone5CaloJets"));
   minEta_ = iConfig.getUntrackedParameter<double>("minEta", -6.5);
   maxEta_ = iConfig.getUntrackedParameter<double>("maxEta",  6.5);
   minPhi_ = iConfig.getUntrackedParameter<double>("minPhi", -3.14);
   maxPhi_ = iConfig.getUntrackedParameter<double>("maxPhi",  3.14);
}


JetEtaPhiFilter::~JetEtaPhiFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
JetEtaPhiFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   edm::Handle<reco::CandidateView> jets;
   iEvent.getByLabel(jsrc_,jets);
   
   vector<const reco::Candidate *> sortedJets;
      
   for(unsigned it=0; it<jets->size(); ++it){
      const reco::Candidate* jet = &((*jets)[it]);
      sortedJets.push_back(jet);
      sortByEtRef (&sortedJets);
   }

   double eta = 0.0;
   double phi = 0.0;

   if(sortedJets.size()>0){

      eta = sortedJets[0]->eta();
      phi = sortedJets[0]->phi();

      return !( eta >= minEta_ &&
		eta <= maxEta_ &&
		phi >= minPhi_ &&
		phi <= maxPhi_  );  // if it falls in {etamin,etamax,phimin,phimax} reject it!
   }else{
      return true;  // if no jet, always ture!
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetEtaPhiFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetEtaPhiFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
JetEtaPhiFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
JetEtaPhiFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
JetEtaPhiFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
JetEtaPhiFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetEtaPhiFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(JetEtaPhiFilter);
