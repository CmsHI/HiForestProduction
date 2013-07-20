// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

// ROOT includes
#include "TH2F.h"
#include "TH1F.h"
#include "TNtuple.h"

// define track and jet spectra analyzer
class OnTrackClusterAnalyzer : public edm::EDAnalyzer {
   public:
      explicit OnTrackClusterAnalyzer(const edm::ParameterSet&);
      ~OnTrackClusterAnalyzer(){};


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob(){} ;

   // ----------member data ---------------------------

      TH1F *hInvalidHitLadder;
      TNtuple *ntInvalidHits;
      TNtuple *ntLostHits;
      TH2F *hOnTrackClustersLayer1ZPhi;
      TH2F *hOnTrackClustersLayer2ZPhi;
      TH2F *hOnTrackClustersLayer3ZPhi;
      TH2F *hTrackEtaPhi;
      edm::Service<TFileService> fs;

      edm::InputTag src_;
      bool setQualityBit_;
      
      const TrackerGeometry * theTracker;
      
};
