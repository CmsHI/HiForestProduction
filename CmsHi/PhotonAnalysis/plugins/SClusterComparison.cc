// -*- C++ -*-
//
// Package:    SClusterComparison
// Class:      SClusterComparison
// 
/**\class SClusterComparison SClusterComparison.cc CmsHi/SClusterComparison/src/SClusterComparison.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Yong Kim,32 4-A08,+41227673039,
//         Created:  Fri Oct 29 12:18:14 CEST 2010
// $Id: SClusterComparison.cc,v 1.3 2010/11/12 12:33:23 kimy Exp $
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

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"


//
// class declaration
//

class SClusterComparison : public edm::EDAnalyzer {
    public:
        explicit SClusterComparison(const edm::ParameterSet&);
        ~SClusterComparison();


    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------

        edm::Service<TFileService> fs;
        edm::InputTag cleanedHybrid_;
        edm::InputTag correctedHybrid_;
        edm::InputTag hybrid_;
        edm::InputTag uncleanedHybrid_;

        edm::InputTag islandBarrel_;
        edm::InputTag correctedIslandBarrel_;

        TTree* theTree;

        int nClean, nCor, nNor, nUncl, nIsl, nCorIsl;

        float cleanedenergy[1000];
        float cleanedrawE[1000];
        float cleanedrawEt[1000];
        float cleanedet[1000];
        float cleanedeta[1000];
        float cleanedphi[1000];


        float cleanedpreshowerenergy[1000];

        float correctedenergy[1000];
        float correctedet[1000];
        float correctedeta[1000];
        float correctedphi[1000];
        float correctedpreshowerenergy[1000];

        float normalenergy[1000];
        float normalet[1000];
        float normaleta[1000];
        float normalphi[1000];
        float normalpreshowerenergy[1000];

        float uncleanedenergy[1000];
        float uncleanedet[1000];
        float uncleanedeta[1000];
        float uncleanedphi[1000];
        float uncleanedpreshowerenergy[1000];

        float islandenergy[1000];
        float islandet[1000];
        float islandrawE[1000];
       float islandrawEt[1000];
float islandeta[1000];
        float islandphi[1000];
        float islandpreshowerenergy[1000];

        float correctedislandenergy[1000];
        float correctedislandet[1000];
        float correctedislandeta[1000];
        float correctedislandphi[1000];
        float correctedislandpreshowerenergy[1000];
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
SClusterComparison::SClusterComparison(const edm::ParameterSet& iConfig)
{
    cleanedHybrid_                  = iConfig.getParameter<edm::InputTag>("cleanedHybridSuperClusters");
    correctedHybrid_                = iConfig.getParameter<edm::InputTag>("correctedHybridSuperClusters");
    hybrid_                         = iConfig.getParameter<edm::InputTag>("hybridSuperClusters");
    uncleanedHybrid_                = iConfig.getParameter<edm::InputTag>("uncleanedHybridSuperClusters");

    islandBarrel_                   = iConfig.getParameter<edm::InputTag>("cleanedIslandBarrelSuperClusters");   // This can be used as the cleanedSC
    correctedIslandBarrel_          = iConfig.getParameter<edm::InputTag>("correctedIslandBarrelSuperClusters");
}


SClusterComparison::~SClusterComparison()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
    void
SClusterComparison::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace reco;

    edm::Handle<SuperClusterCollection> cleanedHybridHandle;
    iEvent.getByLabel(cleanedHybrid_, cleanedHybridHandle);
    const SuperClusterCollection *cleanedHybridCollection = cleanedHybridHandle.product();

    edm::Handle<SuperClusterCollection> correctedHybridHandle;
    iEvent.getByLabel(correctedHybrid_, correctedHybridHandle);
    const SuperClusterCollection *correctedHybridCollection = correctedHybridHandle.product();

    edm::Handle<SuperClusterCollection> HybridHandle;
    iEvent.getByLabel(hybrid_, HybridHandle);
    const SuperClusterCollection *HybridCollection = HybridHandle.product();

    edm::Handle<SuperClusterCollection> uncleanedHybridHandle;
    iEvent.getByLabel(uncleanedHybrid_, uncleanedHybridHandle);
    const SuperClusterCollection *uncleanedHybridCollection = uncleanedHybridHandle.product();

    edm::Handle<SuperClusterCollection> islandBarrelHandle;
    iEvent.getByLabel(islandBarrel_, islandBarrelHandle);
    const SuperClusterCollection *islandBarrelCollection = islandBarrelHandle.product();

    edm::Handle<SuperClusterCollection> correctedIslandBarrelHandle;
    iEvent.getByLabel(correctedIslandBarrel_, correctedIslandBarrelHandle);
    const SuperClusterCollection *correctedIslandBarrelCollection = correctedIslandBarrelHandle.product();

    nClean = 0;
    for(SuperClusterCollection::const_iterator it=cleanedHybridCollection->begin(); it!=cleanedHybridCollection->end(); it++) {
        cleanedenergy[nClean] = it->energy();
        cleanedrawE[nClean] = it->rawEnergy();
        cleanedeta[nClean] = it->eta();
        cleanedphi[nClean] = it->phi();
        cleanedet[nClean] = it->energy()/cosh(it->eta());
        cleanedrawEt[nClean] = it->rawEnergy()/cosh(it->eta());

        nClean++;
    }

    nCor = 0;
    for(SuperClusterCollection::const_iterator it=correctedHybridCollection->begin(); it!=correctedHybridCollection->end(); it++) {
        correctedenergy[nCor] = it->energy();
        correctedeta[nCor] = it->eta();
        correctedphi[nCor] = it->phi();
        correctedet[nCor] = it->energy()/cosh(it->eta());
        nCor++;
    }

    nNor = 0;
    for(SuperClusterCollection::const_iterator it=HybridCollection->begin(); it!=HybridCollection->end(); it++) {
        normalenergy[nNor] = it->energy();
        normaleta[nNor] = it->eta();
        normalphi[nNor] = it->phi();
        normalet[nNor] = it->energy()/cosh(it->eta());
        nNor++;
    }

    nUncl = 0;
    for(SuperClusterCollection::const_iterator it=uncleanedHybridCollection->begin(); it!=uncleanedHybridCollection->end(); it++) {
        uncleanedenergy[nUncl] = it->energy();
        uncleanedeta[nUncl] = it->eta();
        uncleanedphi[nUncl] = it->phi();
        uncleanedet[nUncl] = it->energy()/cosh(it->eta());
        nUncl++;
    }

    nIsl = 0;
    for(SuperClusterCollection::const_iterator it=islandBarrelCollection->begin(); it!=islandBarrelCollection->end(); it++) {
        islandenergy[nIsl] = it->energy();
        islandrawE[nIsl]   = it->rawEnergy();
        islandeta[nIsl] = it->eta();
        islandphi[nIsl] = it->phi();
        islandet[nIsl] = it->energy()/cosh(it->eta());
        islandrawEt[nIsl]   = it->rawEnergy()/cosh(it->eta());
        nIsl++;
    }

    nCorIsl = 0;
    for(SuperClusterCollection::const_iterator it=correctedIslandBarrelCollection->begin(); it!=correctedIslandBarrelCollection->end(); it++) {
        correctedislandenergy[nCorIsl] = it->energy();
        correctedislandeta[nCorIsl] = it->eta();
        correctedislandphi[nCorIsl] = it->phi();
        correctedislandet[nCorIsl] = it->energy()/cosh(it->eta());
        nCorIsl++;
    }

    theTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
    void 
SClusterComparison::beginJob() 
{
    theTree = fs->make<TTree>("clusters","Tree of SuperClusters");

    theTree->Branch("nClean",&nClean,"nClean/I");
    theTree->Branch("cleanedenergy",cleanedenergy,"cleanedenergy[nClean]/F");
    theTree->Branch("cleanedrawE",cleanedrawE,"cleanedrawE[nClean]/F");
    theTree->Branch("cleanedrawEt",cleanedrawEt,"cleanedrawEt[nClean]/F");

    theTree->Branch("cleanedeta",cleanedeta,"cleanedeta[nClean]/F");
    theTree->Branch("cleanedphi",cleanedphi,"cleanedphi[nClean]/F");
    theTree->Branch("cleanedet",cleanedet,"cleanedet[nClean]/F");

    theTree->Branch("nCor",&nCor,"nCor/I");
    theTree->Branch("correctedenergy",correctedenergy,"correctedenergy[nCor]/F");
    theTree->Branch("correctedeta",correctedeta,"correctedeta[nCor]/F");
    theTree->Branch("correctedphi",correctedphi,"correctedphi[nCor]/F");
    theTree->Branch("correctedet",correctedet,"correctedet[nCor]/F");

    theTree->Branch("nNor",&nNor,"nNor/I");
    theTree->Branch("normalenergy",normalenergy,"normalenergy[nNor]/F");
    theTree->Branch("normaleta",normaleta,"normaleta[nNor]/F");
    theTree->Branch("normalphi",normalphi,"normalphi[nNor]/F");
    theTree->Branch("normalet",normalet,"normalet[nNor]/F");

    theTree->Branch("nUncl",&nUncl,"nUncl/I");
    theTree->Branch("uncleanedenergy",uncleanedenergy,"uncleanedenergy[nUncl]/F");
    theTree->Branch("uncleanedeta",uncleanedeta,"uncleanedeta[nUncl]/F");
    theTree->Branch("uncleanedphi",uncleanedphi,"uncleanedphi[nUncl]/F");
    theTree->Branch("uncleanedet",uncleanedet,"uncleanedet[nUncl]/F");

    theTree->Branch("nCleanIsl",&nIsl,"nCleanIsl/I");
    theTree->Branch("cleanIslandenergy",islandenergy,"cleanIslandenergy[nCleanIsl]/F");
    theTree->Branch("cleanIslandrawE",islandrawE,"cleanIslandrawE[nCleanIsl]/F");
    theTree->Branch("cleanIslandrawEt",islandrawEt,"cleanIslandrawEt[nCleanIsl]/F");
    theTree->Branch("cleanIslandeta",islandeta,"cleanIslandeta[nCleanIsl]/F");
    theTree->Branch("cleanIslandphi",islandphi,"cleanIslandphi[nCleanIsl]/F");
    theTree->Branch("cleanIslandet",islandet,"cleanIslandet[nCleanIsl]/F");

    theTree->Branch("nCorIsl",&nCorIsl,"nCorIsl/I");
    theTree->Branch("correctedislandenergy",correctedislandenergy,"correctedislandenergy[nCorIsl]/F");
    theTree->Branch("correctedislandeta",correctedislandeta,"correctedislandeta[nCorIsl]/F");
    theTree->Branch("correctedislandphi",correctedislandphi,"correctedislandphi[nCorIsl]/F");
    theTree->Branch("correctedislandet",correctedislandet,"correctedislandet[nCorIsl]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SClusterComparison::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SClusterComparison);
