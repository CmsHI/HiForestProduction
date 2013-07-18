// -*- C++ -*-
//
// Package:    IsoConeInspector
// Class:      IsoConeInspector
// 
/**\class IsoConeInspector IsoConeInspector.cc CmsHi/IsoConeInspector/src/IsoConeInspector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yong Kim,32 4-A08,+41227673039,
//         Created:  Fri Oct 29 12:18:14 CEST 2010
// $Id: IsoConeInspector.cc,v 1.10 2011/07/18 15:49:01 kimy Exp $
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
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"


#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"


//
// class declaration
//

class IsoConeInspector : public edm::EDAnalyzer {
   public:
      explicit IsoConeInspector(const edm::ParameterSet&);
      ~IsoConeInspector();
   
   
private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

   edm::Service<TFileService> fs;

  edm::InputTag basicClusterBarrel_;
  edm::InputTag basicClusterEndcap_;
  edm::InputTag hbhe_;
  
  edm::InputTag photonSrc_;
   edm::InputTag ebReducedRecHitCollection_;
   edm::InputTag eeReducedRecHitCollection_;
  const CaloGeometry                 *geometry_;
  bool doSpikeClean_;
   double etCut_;
  double hoeCutMin_;
  double hoeCutMax_;
  
  double sieieCut_;
  double etaCut_;
   TH1D*  NoE ;
   TTree* theTree;
  int nPho, nBC, nRH, nHC; 
  int cBin;
   float et;
   float eta;
   float phi;
   float BCet[1000];
   float BCeta[1000];
   float BCphi[1000];
  float HCet[5000];
  float HCeta[5000];
  float HCphi[5000];
  
   float RHet[70000];
   float RHe[70000];
   float RHdEta[70000];
   float RHdPhi[70000];
  float RHeta[70000];
  float RHphi[70000];
  
  bool doStoreCentrality_;
  CentralityProvider *centrality_;
  const CentralityBins * cbins_;
  

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
IsoConeInspector::IsoConeInspector(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   //now do what ever initialization is needed                                                                                                         
   photonSrc_                       = iConfig.getParameter<edm::InputTag>("photonProducer"); // photons                                 
   ebReducedRecHitCollection_       = iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"); //,"reducedEcalRecHitsEB");            
   eeReducedRecHitCollection_       = iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"); //,"reducedEcalRecHitsEE");        
   doSpikeClean_                    = iConfig.getUntrackedParameter<bool>("doSpikeClean",false);
   etCut_                           = iConfig.getUntrackedParameter<double>("etCut",40);
   etaCut_                           = iConfig.getUntrackedParameter<double>("etaCut",1.479);
   doStoreCentrality_                  = iConfig.getUntrackedParameter<bool>("doStoreCentrality",true);
   hoeCutMin_                           = iConfig.getUntrackedParameter<double>("hoeCutMin",0.0);
   hoeCutMax_                           = iConfig.getUntrackedParameter<double>("hoeCutMax",0.2);
   sieieCut_                           = iConfig.getUntrackedParameter<double>("sieieCut",0.01);

   basicClusterBarrel_              = iConfig.getParameter<edm::InputTag>("basicClusterBarrel");
   hbhe_                            = iConfig.getParameter<edm::InputTag>("bhbe");

   
}


IsoConeInspector::~IsoConeInspector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
IsoConeInspector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //grab the photon collection
   NoE->Fill(0);


   ESHandle<CaloGeometry> geometryHandle;
   iSetup.get<CaloGeometryRecord>().get(geometryHandle);
   if(geometryHandle.isValid())
     geometry_ = geometryHandle.product();
   else
     geometry_ = NULL;


   Handle<reco::PhotonCollection> photonColl;
   iEvent.getByLabel(photonSrc_, photonColl);
   const reco::PhotonCollection *photons = photonColl.product();
   reco::PhotonCollection::const_iterator pho;
   reco::PhotonCollection::const_iterator leadingPho;

   //grab rechits                                                                                                                                      
   edm::Handle<EcalRecHitCollection> EBReducedRecHits;
   iEvent.getByLabel(ebReducedRecHitCollection_, EBReducedRecHits);
   const EcalRecHitCollection* rechitsCollectionBarrel = EBReducedRecHits.product();
   
   edm::Handle<EcalRecHitCollection> EEReducedRecHits;
   iEvent.getByLabel(eeReducedRecHitCollection_, EEReducedRecHits);
   const EcalRecHitCollection* rechitsCollectionEndcap = EEReducedRecHits.product();
  //lazy tool                                                                                                                                     
   EcalClusterLazyTools lazyTool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_ );
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
   
   cBin = -10;
   if (doStoreCentrality_) {
     if(!centrality_) centrality_ = new CentralityProvider(iSetup);
     centrality_->newEvent(iEvent,iSetup);
     cBin = (int)centrality_->getBin();
   }
   
   
   //grab basicClusters                                                                                                                
   Handle<reco::CaloClusterCollection> basicClusterB;
   iEvent.getByLabel(basicClusterBarrel_, basicClusterB);

   reco::CaloClusterCollection myBCs;
   for (reco::CaloClusterCollection::const_iterator bcItr = basicClusterB->begin(); bcItr != basicClusterB->end(); ++bcItr) {
     myBCs.push_back(*bcItr);
   }
   
   
   Handle<HBHERecHitCollection> hehbhandle;
   iEvent.getByLabel(hbhe_, hehbhandle);
   const HBHERecHitCollection         *fHBHERecHits_;
   if(hehbhandle.isValid())
     fHBHERecHits_ = hehbhandle.product();
   else
     fHBHERecHits_ = NULL;



   
   
   double extRadius_ = 100.;
   double etaWidth_ = 2.0;
   nRH  = 0;
   nBC  = 0;
   nHC = 0;
   nPho = (*photons).size();
      
   for (pho = (*photons).begin(); pho!= (*photons).end(); pho++){
     if ( cBin > 4 )
       continue;
     
     et       = (float)pho->et();
     eta      = (float)pho->superCluster()->eta();
     phi      = (float)pho->superCluster()->phi();
     if ( et < etCut_ )  continue;
     if ( fabs(eta) > etaCut_ ) continue;
     if ( pho->hadronicOverEm() > hoeCutMax_) continue;
     if ( pho->hadronicOverEm() < hoeCutMin_) continue;

     if ( pho->sigmaIetaIeta()  > sieieCut_) continue;

      
      math::XYZPoint theCaloPosition = pho->superCluster()->position();
      GlobalPoint pclu (theCaloPosition.x () , theCaloPosition.y () , theCaloPosition.z () );
      
      int subdetnr(0);
      subdetnr = 0;  // barrel
      CaloSubdetectorGeometry::DetIdSet chosen = subdet_[subdetnr]->getCells(pclu,extRadius_);// select cells around cluster
      CaloRecHitMetaCollectionV::const_iterator j=RecHitsBarrel->end();
      
      
      for (reco::CaloClusterCollection::const_iterator bcItr = myBCs.begin(); bcItr != myBCs.end(); ++bcItr) {
        BCet[nBC] = bcItr->energy()/cosh(bcItr->eta());
	BCeta[nBC] = bcItr->eta();
	BCphi[nBC] = bcItr->phi();
	nBC++;
      }

      
      
      nRH = 0;
      
      
      for (CaloSubdetectorGeometry::DetIdSet::const_iterator  i = chosen.begin ();i!= chosen.end ();++i){//loop selected cells
	j=RecHitsBarrel->find(*i); // find selected cell among rechi ts
	if( j!=RecHitsBarrel->end()){ // add rechit only if available 
	  const  GlobalPoint & position = theCaloGeom.product()->getPosition(*i);
	  float etarh = position.eta();
	  float phirh = position.phi();
	  float etaDiff = etarh - eta;
	  float phiDiff= deltaPhi(phirh,phi);
	  float energyrh = j->energy();
	  float etrh    = energyrh/cosh(etarh);
	  
	  if ( fabs(etaDiff) > etaWidth_)  
	    continue;
	  RHet[nRH]   = etrh;
	  RHe[nRH]   = energyrh;
	  RHdEta[nRH] = etaDiff;
	  RHdPhi[nRH] = phiDiff;
	  RHeta[nRH]  = etarh;
	  RHphi[nRH]  = phirh;

	  nRH++;
	}
      }
      
      // Hcal
      nHC = 0;
      for(size_t index = 0; index < fHBHERecHits_->size(); index++) {
	const HBHERecHit &rechit = (*fHBHERecHits_)[index];
	const DetId &detid = rechit.id();
	const GlobalPoint& hitpoint = geometry_->getPosition(detid);
	HCeta[nHC] = hitpoint.eta();
	HCphi[nHC] = hitpoint.phi();
	HCet[nHC]  = rechit.energy()/cosh( hitpoint.eta());
	
	nHC++;
      }
      
      
      
      
      theTree->Fill();
   }
   
   
   
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
IsoConeInspector::beginJob() 
{
  centrality_ = 0;

   NoE      = fs->make<TH1D>( "NoE"  , "", 1,  -100., 100. );
   theTree  = fs->make<TTree>("photon","Tree of Rechits around photon");
   theTree->Branch("nPho",&nPho,"nPho/I");
   theTree->Branch("nBC",&nBC,"nBC/I");
   theTree->Branch("nHC",&nHC,"nHC/I");
   theTree->Branch("nRH",&nRH,"nRH/I");
   theTree->Branch("cBin",&cBin,"cBin/I");
   
   theTree->Branch("et",&et,"et/F");
   theTree->Branch("eta",&eta,"eta/F");
   theTree->Branch("phi",&phi,"phi/F");

   theTree->Branch("BCet",BCet,"BCet[nBC]/F");
   theTree->Branch("BCeta",BCeta,"BCeta[nBC]/F");
   theTree->Branch("BCphi",BCphi,"BCphi[nBC]/F");

   theTree->Branch("rhet",RHet,"rhet[nRH]/F");
   theTree->Branch("rhe",RHe,"rhe[nRH]/F");
   theTree->Branch("rhdEta",RHdEta,"rhdEta[nRH]/F");
   theTree->Branch("rhdPhi",RHdPhi,"rhdPhi[nRH]/F");
   theTree->Branch("rheta",RHeta,"rheta[nRH]/F");
   theTree->Branch("rhphi",RHphi,"rhphi[nRH]/F");
   
   theTree->Branch("HCet",HCet,"HCet[nHC]/F");
   theTree->Branch("HCeta",HCeta,"HCeta[nHC]/F");
   theTree->Branch("HCphi",HCphi,"HCphi[nHC]/F");
   

    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
IsoConeInspector::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoConeInspector);
