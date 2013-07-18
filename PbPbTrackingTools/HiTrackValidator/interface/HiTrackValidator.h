// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// for jet,track,vertices correlation 
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// sim <--> rec track association
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

// Pixel hit information
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

// Calotower matching
#include "DataFormats/Math/interface/deltaR.h"
#include <Math/DistFunc.h>
#include "TMath.h"

// centrality
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"



// ROOT includes
#include <TH2.h>
#include <TMath.h>
#include <TF1.h>

class HiTrackValidator : public edm::EDAnalyzer {
   public:
      explicit HiTrackValidator(const edm::ParameterSet&);
      ~HiTrackValidator(){};

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      void matchByDrAllowReuse(const reco::Track & trk, const edm::Handle<CaloTowerCollection> & towers, double & bestdr, double & bestpt);
      bool hitDeadPXF(const reco::Track&);
      typedef math::XYZPoint Point;

      // ----------member data ---------------------------
      edm::Service<TFileService> f;

      edm::InputTag bslabel_;
      edm::InputTag vtxlabel_;
      edm::InputTag trklabel_;
      edm::InputTag jetlabel_;
      edm::InputTag simtrklabel_;
      edm::InputTag towerlabel_;
      edm::InputTag associatorMap_;

      double etaMax_;
      double jetEtMin_;

      bool hasSimInfo_;
      bool selectFake_;
      bool hasCaloMat_;
      double towerPtMin_;
      bool useQaulityStr_;
      std::string qualityString_;

      bool fiducialCut_;

      // string of functional form
      std::string funcDeltaRTowerMatch_;
      std::string funcCaloComp_;
      TF1 *fDeltaRTowerMatch, *fCaloComp;

      std::vector<int32_t> neededCentBins_;

      CentralityProvider * centrality_;
      
      const TrackerGeometry * theTracker;
      
      // quality cuts distribution 
      TH1F *hVtxSize;
      TH1F *hNlayers;
      TH1F *hNlayers3D;
      TH1F *hNlayersLost;
      TH1F *hNvalidHits;
      TH1F *hChi2n;
      TH1F *hRelPtErr;
      TH1F *hdzPV, *hd0PV;
      TH1F *hdzErr, *hd0Err;
      TH1F *hdzErrPV, *hd0ErrPV;
      TH1F *hvzErr, *hvtErr;
      TH1F *hdzOverdzErr, *hd0Overd0Err;
      TH1F *hdzOverdzErrPV, *hd0Overd0ErrPV;
      TH1F *hEta, *hPhi, *hPt;
      TH1F *hQualityName;

      std::vector<TH1F*> hNlayers_Cent;
      std::vector<TH1F*> hNlayers3D_Cent;
      std::vector<TH1F*> hNlayersLost_Cent;
      std::vector<TH1F*> hNvalidHits_Cent;
      std::vector<TH1F*> hChi2n_Cent;
      std::vector<TH1F*> hRelPtErr_Cent;
      std::vector<TH1F*> hdzPV_Cent;
      std::vector<TH1F*> hd0PV_Cent;
      std::vector<TH1F*> hdzErrPV_Cent;
      std::vector<TH1F*> hd0ErrPV_Cent;
      std::vector<TH1F*> hdzOverdzErrPV_Cent;
      std::vector<TH1F*> hd0Overd0ErrPV_Cent;
      std::vector<TH1F*> hEta_Cent;
      std::vector<TH1F*> hPhi_Cent;
      std::vector<TH1F*> hPt_Cent;


      // as a function of pt
      TH2F *hNlayersdPt;
      TH2F *hNlayers3DdPt;
      TH2F *hNlayersLostdPt;
      TH2F *hNvalidHitsdPt;
      TH2F *hChi2ndPt;
      TH2F *hRelPtErrdPt;
      TH2F *hdzPVdPt, *hd0PVdPt;
      TH2F *hdzErrdPt, *hd0ErrdPt;
      TH2F *hdzErrPVdPt, *hd0ErrPVdPt;
      TH2F *hvzErrdPt, *hvtErrdPt;
      TH2F *hdzOverdzErrdPt, *hd0Overd0ErrdPt;
      TH2F *hdzOverdzErrPVdPt, *hd0Overd0ErrPVdPt;

      std::vector<TH2F*> hNvalidHitsdPt_Cent;
      std::vector<TH2F*> hChi2ndPt_Cent;
      std::vector<TH2F*> hRelPtErrdPt_Cent;
      std::vector<TH2F*> hdzPVdPt_Cent;
      std::vector<TH2F*> hd0PVdPt_Cent;
      std::vector<TH2F*> hdzErrPVdPt_Cent;
      std::vector<TH2F*> hd0ErrPVdPt_Cent;
      std::vector<TH2F*> hdzOverdzErrPVdPt_Cent;
      std::vector<TH2F*> hd0Overd0ErrPVdPt_Cent;

      // correlations
      //TH2D *hdzErrd0Err;
      //TH2D *hdzErrd0ErrPV;
      //TH2D *hdzOverdzErrd0Err;
      //TH2D *hdzOverdzErrd0ErrPV;

      //std::vector<TH2D*> hdzOverdzErrd0Err_Cent;
      //std::vector<TH2D*> hdzOverdzErrd0ErrPV_Cent;
      
      // kinematic distributions
      TH2F *hEtaPhi;
      //std::vector<TH2D*> hEtaPhi_Cent;

      // track-calorimeter matching
      TH1F *hdR;
      TH1F *hPtCaloMat, *hPtReal, *hPtRealCaloMat;
      TH2F *hdRdPt, * hdRdCalE, *hdPtdCalE;

      std::vector<double> ptBins;


};
