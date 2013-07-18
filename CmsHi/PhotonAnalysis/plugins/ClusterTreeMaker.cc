// -*- C++ -*-
//
// Package:    ClusterTreeMaker
// Class:      ClusterTreeMaker
// 
/**\class ClusterTreeMaker ClusterTreeMaker.cc CmsHi/ClusterTreeMaker/src/ClusterTreeMaker.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Yong Kim,32 4-A08,+41227673039,
//         Created:  Fri Oct 29 12:18:14 CEST 2010
// $Id: ClusterTreeMaker.cc,v 1.13 2011/11/12 14:33:41 kimy Exp $
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

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollectionV.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include <Math/VectorUtil.h>
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "TNtuple.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

//
// class declaration
//
float recHitTime(DetId id, const EcalRecHitCollection *recHits);
DetId getMaximumRecHit(const reco::BasicCluster &cluster, const EcalRecHitCollection *recHits);

class ClusterTreeMaker : public edm::EDAnalyzer {
    public:
        explicit ClusterTreeMaker(const edm::ParameterSet&);
        ~ClusterTreeMaker();


    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  edm::InputTag ebReducedRecHitCollection_;
  edm::InputTag eeReducedRecHitCollection_;

  edm::InputTag basicClusterBarrel_;
  edm::InputTag basicClusterEndcap_;

  TTree* theTree;
  TTree* theBC;
  
  bool doRecHit;
  TH1D* hRHetBarrel;
  TH1D* hRHetEndcap;
  TH1D* hRHetEndcap2;
  TH1D* hRHetBarrelCleaned;

  TH1D* hRHetaEndcap;
  TH1D* hRHphiEndcap;
  TH1D* hRHetaBarrel;
  TH1D* hRHphiBarrel;

  std::string mSrc;
  std::string vertexProducer_;      // vertecies producer                                                                                                                                                       
  std::string scProducerB_;
  std::string scProducerE_;
  float et[3000];
  float eta[3000];
  float phi[3000];
  float energy[3000];
  float rawEnergy[3000];
  
  float severity[3000];
  float time[3000];
  float swissCrx[3000];
  int nPar;
  double etCut;
   
   double etCutBC;

  bool doBasicCluster;
  int nBC ; 
  float bcet[3000];
  float bceta[3000];
  float bcphi[3000];

  int nRHBarrel ;
  int nRHEndcap ;

  
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
ClusterTreeMaker::ClusterTreeMaker(const edm::ParameterSet& iConfig)
   
{
   scProducerB_  = iConfig.getUntrackedParameter<std::string>("scTagB","hltRecoHIEcalWithCleaningCandidate");
   scProducerE_  = iConfig.getUntrackedParameter<std::string>("scTagE","hltRecoHIEcalWithCleaningCandidate");
   ebReducedRecHitCollection_       = iConfig.getUntrackedParameter<edm::InputTag>("ebRecHitCollection");
   eeReducedRecHitCollection_       = iConfig.getUntrackedParameter<edm::InputTag>("eeRecHitCollection");
   etCut                            = iConfig.getUntrackedParameter<double>("etCutSC",8.0);
   etCutBC                          = iConfig.getUntrackedParameter<double>("etCutBC",1.0);
   basicClusterBarrel_              = iConfig.getParameter<edm::InputTag>("basicClusterBarrel");
   basicClusterEndcap_              = iConfig.getParameter<edm::InputTag>("basicClusterEndcap");
   doRecHit                         = iConfig.getUntrackedParameter<bool>("doRecHit",true);
   doBasicCluster                   = iConfig.getUntrackedParameter<bool>("doBasicCluster",true);

   
}


ClusterTreeMaker::~ClusterTreeMaker()

{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
    void
    ClusterTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;


   EcalClusterLazyTools lazyTool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_);
   
   Handle<EcalRecHitCollection> barrelRecHitsHandle;
   iEvent.getByLabel(ebReducedRecHitCollection_,barrelRecHitsHandle);
   const EcalRecHitCollection* ecalRecHitsEB = 0;
   if(!barrelRecHitsHandle.isValid()) {
     LogDebug("ClusterTreeMaker") << "Error! Can't get barrelRecHitsHandle product!" << std::endl;
   } else {
     ecalRecHitsEB = barrelRecHitsHandle.product();
   }

   Handle<EcalRecHitCollection> endcapRecHitsHandle;
   iEvent.getByLabel(eeReducedRecHitCollection_,endcapRecHitsHandle);
   const EcalRecHitCollection* ecalRecHitsEE = 0;
   if(!barrelRecHitsHandle.isValid()) {
     LogDebug("ClusterTreeMaker") << "Error! Can't get endcapRecHitsHandle product!" << std::endl;
   } else {
     ecalRecHitsEE = endcapRecHitsHandle.product();
   }





   Handle<reco::SuperClusterCollection> superClusterHandleB;
   iEvent.getByLabel(scProducerB_, superClusterHandleB);
   const reco::SuperClusterCollection* superClusterCollectionB = 0;
   if(!superClusterHandleB.isValid()) {
     LogDebug("SpikeInspector") << "Error! Can't get superClusterHandle product!" << std::endl;
   } else {
     superClusterCollectionB = superClusterHandleB.product();
   }

   Handle<reco::SuperClusterCollection> superClusterHandleE;
   iEvent.getByLabel(scProducerE_, superClusterHandleE);
   const reco::SuperClusterCollection* superClusterCollectionE = 0;

   if(!superClusterHandleE.isValid()) {
     LogDebug("SpikeInspector") << "Error! Can't get superClusterHandle product!" << std::endl;
   } else {
     superClusterCollectionE = superClusterHandleE.product();
   }

   edm::ESHandle<EcalSeverityLevelAlgo> ecalSevLvlAlgoHndl;
   iSetup.get<EcalSeverityLevelAlgoRcd>().get(ecalSevLvlAlgoHndl);

   int maxindex = (int)superClusterCollectionB->size();
   
   nPar = 0;
   // cout << " nPhoton barrel = " << maxindex << endl;
   for(int i = 0; i < maxindex; i++)  {
     const reco::SuperCluster &c1 = (*superClusterCollectionB)[i];
     
     float theEt = c1.energy()/cosh(c1.eta());
     if ( theEt < etCut )  continue;
     
     float theSeverity = -100;
     float theSwissCrx = -100;
     float theTime = -100;
     
     const reco::CaloClusterPtr seed = c1.seed();
     DetId id = lazyTool.getMaximum(*seed).first;
     const EcalRecHitCollection & rechits = *ecalRecHitsEB;
     EcalRecHitCollection::const_iterator it = rechits.find( id );
     if( it != rechits.end() ) {
       theSeverity = ecalSevLvlAlgoHndl->severityLevel(id, rechits);
       theSwissCrx = EcalTools::swissCross   (id, rechits, 0.,true);
       theTime = it->time();
       //      std::cout << "swissCross = " << swissCrx <<std::endl;                                               
       // std::cout << " timing = " << it->time() << std::endl;                                                    
     }

     severity[nPar] = theSeverity;
     time[nPar]    = theTime;
     swissCrx[nPar] =  theSwissCrx;
     et[nPar] = theEt;
     eta[nPar] = (float)c1.eta();
     phi[nPar] = (float)c1.phi();
     energy[nPar] = (float)c1.energy();
     rawEnergy[nPar] = (float)c1.rawEnergy();
     
     nPar++;
   }
   
   maxindex = (int)superClusterCollectionE->size();
   //  cout << " nPhoton Endcap = " << maxindex << endl;
   for(int i = 0; i < maxindex; i++) {
     const reco::SuperCluster &c1 = (*superClusterCollectionE)[i];

     float theEt = c1.energy()/cosh(c1.eta());
     if ( theEt < etCut )  continue;
     float theSeverity = -100;
     float theSwissCrx = -100;
     float theTime = -100;

     const reco::CaloClusterPtr seed = c1.seed();
     DetId id = lazyTool.getMaximum(*seed).first;
     const EcalRecHitCollection & rechits = *ecalRecHitsEE;
     EcalRecHitCollection::const_iterator it = rechits.find( id );
     if( it != rechits.end() ) {
       theSeverity= ecalSevLvlAlgoHndl->severityLevel(id, rechits);
       theSwissCrx = EcalTools::swissCross   (id, rechits, 0.,true);
       theTime = it->time();
       //      std::cout << "swissCross = " << swissCrx <<std::endl;                                                    
       // std::cout << " timing = " << it->time() << std::endl;                                                         
     }

     severity[nPar] = theSeverity;
     time[nPar]    = theTime;
     swissCrx[nPar] =  theSwissCrx;
     et[nPar] = theEt;

     et[nPar] = theEt;
     eta[nPar] = (float)c1.eta();
     phi[nPar] = (float)c1.phi();
     energy[nPar] = (float)c1.energy();
     rawEnergy[nPar] = (float)c1.rawEnergy();
     

     nPar++;
   }


   
      
   // Rechit loop
   //get the rechit geometry
   edm::ESHandle<CaloGeometry> theCaloGeom;
   const CaloGeometry* caloGeom;
   
   //Barrel 
   nRHBarrel=0;
   nRHEndcap=0;

   EcalRecHitCollection::const_iterator rh;
   if ( doRecHit) {
     
     iSetup.get<CaloGeometryRecord>().get(theCaloGeom);
     caloGeom = theCaloGeom.product();
     
     EcalRecHitCollection::const_iterator rh;
     
     for (rh = (*ecalRecHitsEB).begin(); rh!= (*ecalRecHitsEB).end(); rh++){
       DetId id = rh->id();
       double swissCrx = EcalTools::swissCross   (id, *ecalRecHitsEB, 0, true) ; // EcalSeverityLevelAlgo::swissCross(id, *rechitsCollectionBarrel, 0.,true);
       double time     = rh->time();
       
       const GlobalPoint & position = caloGeom->getPosition(rh->id());
       double tempEt = rh->energy()/cosh(position.eta()) ;
       hRHetBarrel->Fill(tempEt);
       hRHetaBarrel->Fill ( position.eta() ) ;
       hRHphiBarrel->Fill ( position.phi() ) ;
       
       if ( (rh->energy()>3) && ((fabs(time) > 3)||(swissCrx > 0.95)) )
	 hRHetBarrelCleaned->Fill(tempEt);
     
       nRHBarrel++;
     }
     //Endcap
     for (rh = (*ecalRecHitsEE).begin(); rh!= (*ecalRecHitsEE).end(); rh++){
       DetId id = rh->id();
       const GlobalPoint & position = caloGeom->getPosition(rh->id());
       double tempEt = rh->energy()/cosh(position.eta()) ;
       hRHetEndcap->Fill(tempEt);
       if ( fabs(position.eta()) < 2 )   
	 hRHetEndcap2->Fill(tempEt);
       
       
       hRHetaEndcap->Fill ( position.eta() ) ;
       hRHphiEndcap->Fill ( position.phi() ) ;
       
       nRHEndcap++;
     }
   }
   
   Handle<reco::CaloClusterCollection> basicClusterB;
   iEvent.getByLabel(basicClusterBarrel_, basicClusterB);
   Handle<reco::CaloClusterCollection> basicClusterE;
   iEvent.getByLabel(basicClusterEndcap_, basicClusterE);
   reco::CaloClusterCollection myBCs;
   for (reco::CaloClusterCollection::const_iterator bcItr = basicClusterB->begin(); bcItr != basicClusterB->end(); ++bcItr) {
     myBCs.push_back(*bcItr);
   }
   for (reco::CaloClusterCollection::const_iterator bcItr = basicClusterE->begin(); bcItr != basicClusterE->end(); ++bcItr) {
     myBCs.push_back(*bcItr);
   }
   
   nBC = 0;
   if (doBasicCluster) {
      for (reco::CaloClusterCollection::const_iterator bcItr = myBCs.begin(); bcItr != myBCs.end(); ++bcItr) {
	 if ( (  bcItr->energy()/cosh(bcItr->eta()) ) < etCutBC ) 
	   continue;
	 bceta[nBC] = bcItr->eta();
	 bcphi[nBC] = bcItr->phi();
	 bcet[nBC] =  bcItr->energy()/cosh(bcItr->eta());
	 nBC++;
      }
      
   }  
   
   theTree->Fill();
   theBC->Fill();
   
}


   // ------------ method called once each job just before starting event loop  ------------
    void 
ClusterTreeMaker::beginJob() 
{
   
  theTree  = fs->make<TTree>("superCluster","Tree of supercluster");
   
   theTree->Branch("nPar",&nPar,"nPar/I");
   theTree->Branch("et",et,"et[nPar]/F");
   theTree->Branch("energy",energy,"energy[nPar]/F");
   theTree->Branch("rawEnergy",rawEnergy,"rawEnergy[nPar]/F");
   theTree->Branch("eta",eta,"eta[nPar]/F");
   theTree->Branch("phi",phi,"phi[nPar]/F");
   theTree->Branch("time",time,"time[nPar]/F");
   theTree->Branch("swissCrx",swissCrx,"swissCrx[nPar]/F");
   theTree->Branch("severity",severity,"severity[nPar]/F");
   theTree->Branch("nRHBarrel",&nRHBarrel,"nRHBarrel/I");
   theTree->Branch("nRHEndcap",&nRHEndcap,"nRHEndcap/I");
 
   theBC  = fs->make<TTree>("basicCluster","Tree of basiccluster");

   theBC->Branch("number",&nBC,"number/I");
   theBC->Branch("et",bcet,"et[number]/F");
   theBC->Branch("eta",bceta,"eta[number]/F");
   theBC->Branch("phi",bcphi,"phi[number]/F");

   

   hRHetBarrel = fs->make<TH1D>( "hRHetBarrel" , "", 4000,-2.5,197.5);;
   hRHetEndcap = fs->make<TH1D>( "hRHetEndcap" , "", 4000,-2.5,197.5);;
   hRHetEndcap2 = fs->make<TH1D>( "hRHetEndcap2" , "", 4000,-2.5,197.5);;
   hRHetaBarrel = fs->make<TH1D>( "hRHetaBarrel" , "", 400,-3,3);
   hRHetaEndcap = fs->make<TH1D>( "hRHetaEndcap" , "", 400,-3,3);
   hRHphiBarrel = fs->make<TH1D>( "hRHphiBarrel" , "", 400,-4,4);
   hRHphiEndcap = fs->make<TH1D>( "hRHphiEndcap" , "", 400,-4,4);
   
   hRHetBarrelCleaned = fs->make<TH1D>( "hRHetBarrelCleaned" , "", 400,-2.5,197.5);;
   
   
   std::cout<<"done beginjob"<<std::endl;
}




// ------------ method called once each job just after ending the event loop  ------------
void 
ClusterTreeMaker::endJob() {
}

float recHitTime(DetId id, const EcalRecHitCollection *recHits) {
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits->find( id );
    if ( it != recHits->end() ) {
      return (*it).time();
    } else {
      //throw cms::Exception("EcalRecHitNotFound") << "The recHit corresponding to the DetId" << id.rawId() << " not found in the EcalRecHitCollection";                                                    
      // the recHit is not in the collection (hopefully zero suppressed)                                                                                                                                    
      return 0;
    }
  }
  return 0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(ClusterTreeMaker);
