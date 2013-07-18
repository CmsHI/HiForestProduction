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

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

// ROOT includes
#include <TH2.h>

class HiVertexAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HiVertexAnalyzer(const edm::ParameterSet&);
      ~HiVertexAnalyzer(){};

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag vtxlabel_;
      int nTracksBins_;

      edm::Service<TFileService> f;
      TH1F *hVtxSize;
      TH1F *hVtxTracks;

      TH1F *hVtxZ;
      TH1F *hVtxY;
      TH1F *hVtxX;

      TH2F *hMultVtxTracks;
      TH2F *hMultVtxZ;
      TH1F *hAllVtxZ;
      
      TH2F *hVtxXnY;
      TH2F *hVtxXnZ;
      TH2F *hVtxYnZ;

      std::vector<TH1F*> hVtxZ_Cent;
      std::vector<TH1F*> hVtxX_Cent;
      std::vector<TH1F*> hVtxY_Cent;

      std::vector<TH1F*> hVtxTracks_Cent;

      std::vector<int32_t> neededCentBins_;

      CentralityProvider * centrality_;

};
