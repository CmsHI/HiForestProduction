// -*- C++ -*-
//
// Package:    RechitCount
// Class:      RechitCount
// 
/**\class RechitCount RechitCount.cc CmsHi/RechitCount/src/RechitCount.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yong Kim,32 4-A08,+41227673039,
//         Created:  Fri Oct 29 12:18:14 CEST 2010
// $Id: RechitCount.cc,v 1.1 2011/01/27 23:44:57 kimy Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollectionV.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include <Math/VectorUtil.h>
#include "DataFormats/Math/interface/deltaPhi.h"
//
// class declaration
//

class RechitCount : public edm::EDAnalyzer {
   public:
      explicit RechitCount(const edm::ParameterSet&);
      ~RechitCount();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

   edm::Service<TFileService> fs;

   edm::InputTag photonSrc_;
   edm::InputTag ebAllRecHitCollection_;
   edm::InputTag eeAllRecHitCollection_;
   double etCut_;
   TH1D*  NoE ;
   TH1D* hbar;
   TH1D* hend;
   TTree* theTree;
   int nPho, nBC, nRH; 
   float et;
   float eta;
   float phi;
   float BCet[1000];
   float BCeta[1000];
   float BCphi[1000];
   float RHet[5000];
   float RHdEta[5000];
   float RHdPhi[5000];



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
RechitCount::RechitCount(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   //now do what ever initialization is needed                                                                                                         
   photonSrc_                       = iConfig.getParameter<edm::InputTag>("photonProducer"); // photons                                 
   ebAllRecHitCollection_       = iConfig.getParameter<edm::InputTag>("ebAllRecHitCollection"); //,"reducedEcalRecHitsEB");            
   eeAllRecHitCollection_       = iConfig.getParameter<edm::InputTag>("eeAllRecHitCollection"); //,"reducedEcalRecHitsEE");        
   etCut_                       = iConfig.getUntrackedParameter<double>("etCut",20);
     
}


RechitCount::~RechitCount()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
RechitCount::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //grab the photon collection
   NoE->Fill(0);


   Handle<reco::PhotonCollection> photonColl;
   iEvent.getByLabel(photonSrc_, photonColl);
   const reco::PhotonCollection *photons = photonColl.product();
   reco::PhotonCollection::const_iterator pho;
   reco::PhotonCollection::const_iterator leadingPho;

   //grab rechits                                                                                                                                      
   edm::Handle<EcalRecHitCollection> EBAllRecHits;
   iEvent.getByLabel(ebAllRecHitCollection_, EBAllRecHits);
   const EcalRecHitCollection* rechitsCollectionBarrel = EBAllRecHits.product();
   
   edm::Handle<EcalRecHitCollection> EEAllRecHits;
   iEvent.getByLabel(eeAllRecHitCollection_, EEAllRecHits);
   const EcalRecHitCollection* rechitsCollectionEndcap = EEAllRecHits.product();
  //lazy tool                                                                                                                                     
   EcalClusterLazyTools lazyTool(iEvent, iSetup, ebAllRecHitCollection_, eeAllRecHitCollection_ );
   // get the channel status from the DB                                                                                                     
   edm::ESHandle<EcalChannelStatus> chStatus;
   iSetup.get<EcalChannelStatusRcd>().get(chStatus);

   // get the RecHit geometry
   edm::ESHandle<CaloGeometry> theCaloGeom;
   iSetup.get<CaloGeometryRecord>().get(theCaloGeom);
   const CaloGeometry* caloGeom = theCaloGeom.product();
   const CaloSubdetectorGeometry* subdet_[2];
   subdet_[0] = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
   subdet_[1] = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);


   
   std::vector< std::pair<DetId, float> >::const_iterator rhIt;

   std::auto_ptr<CaloRecHitMetaCollectionV> RecHitsBarrel(0); 
   RecHitsBarrel = std::auto_ptr<CaloRecHitMetaCollectionV>(new EcalRecHitMetaCollection(*rechitsCollectionBarrel));
   std::auto_ptr<CaloRecHitMetaCollectionV> RecHitsEndcap(0);
   RecHitsEndcap = std::auto_ptr<CaloRecHitMetaCollectionV>(new EcalRecHitMetaCollection(*rechitsCollectionEndcap));

   
   double extRadius_ = 100.;
   
   nPho = 0;
   nRH  = 0;
   nBC  = 0;
    for (pho = (*photons).begin(); pho!= (*photons).end(); pho++){
       int valBar(0), valEnd(0);
       et       = (float)pho->et();
       eta      = (float)pho->superCluster()->eta();
       phi      = (float)pho->superCluster()->phi();
       if ( et < etCut_ )  continue;
       
       
       math::XYZPoint theCaloPosition = pho->superCluster()->position();
       GlobalPoint pclu (theCaloPosition.x () , theCaloPosition.y () , theCaloPosition.z () );
       
       CaloSubdetectorGeometry::DetIdSet chosen = subdet_[0]->getCells(pclu,extRadius_);// select cells around cluster
       CaloRecHitMetaCollectionV::const_iterator j=RecHitsBarrel->end();
       CaloSubdetectorGeometry::DetIdSet chosenE = subdet_[1]->getCells(pclu,extRadius_);// select cells around cluster            
       CaloRecHitMetaCollectionV::const_iterator jE=RecHitsEndcap->end();

       for (CaloSubdetectorGeometry::DetIdSet::const_iterator  i = chosen.begin ();i!= chosen.end ();++i){//loop selected cells
	  j=RecHitsBarrel->find(*i); // find selected cell among rechits
	  if( j!=RecHitsBarrel->end()){ // add rechit only if available 
	     const  GlobalPoint & position = theCaloGeom.product()->getPosition(*i);
	     double etarh = position.eta();
	     double phirh = position.phi();
	     double etaDiff = etarh - eta;
	     if ( abs(etaDiff) < 0.5 ) {
		if ( abs(eta)<1.479) valBar++;
		else valEnd++;
	     }
	     
	  }
       }
       
       for (CaloSubdetectorGeometry::DetIdSet::const_iterator  i = chosenE.begin ();i!= chosenE.end ();++i){//loop selected cells               
          jE=RecHitsEndcap->find(*i); // find selected cell among rechits                                                                  
          if( jE!=RecHitsEndcap->end()){ // add rechit only if available                                                                       
             const  GlobalPoint & position = theCaloGeom.product()->getPosition(*i);
             double etarh = position.eta();
             double phirh = position.phi();
             double etaDiff = etarh - eta;
	     if ( abs(etaDiff) < 0.5 ) {
                if ( abs(eta)<1.479) valBar++;
                else valEnd++;
             }
          }
       }
       
       
       if ( abs(eta)<1.479 ) hbar->Fill(valBar);
       else                  hend->Fill(valEnd);
    }
    
    
    
}


// ------------ method called once each job just before starting event loop  ------------
void 
RechitCount::beginJob() 
{
   NoE      = fs->make<TH1D>( "NoE"  , "", 1,  -100., 100. );
   hbar     = fs->make<TH1D>( "hbar" , "", 1000,-.5,10000.5);
   hend    = fs->make<TH1D>( "hend" , "", 1000,-.5,10000.5);
   
   
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RechitCount::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(RechitCount);
