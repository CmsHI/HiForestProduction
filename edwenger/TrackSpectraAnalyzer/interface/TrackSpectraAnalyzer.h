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
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

// ROOT includes
#include "TNtuple.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"

// define track and jet spectra analyzer
class TrackSpectraAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackSpectraAnalyzer(const edm::ParameterSet&);
      ~TrackSpectraAnalyzer(){};


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob(){} ;

   // ----------member data ---------------------------

   TNtuple *nt_dndptdeta; 
   TNtuple *nt_gen_dndptdeta;
   TNtuple *nt_jet;
   TNtuple *nt_jettrack;

   TH1F *hNevt;
   TH1F *hNevt_mult1;
   TH1F *hNevt_mult2;
   TH1F *hNevt_mult3;

   TH1F *hRecMult_STD;
   TH1F *hRecMult_STD_corr;
   TH1F *hRecJetEt_STD_corr;

   TH1F *hGenNevt;


   TH2F *hTrkPtEta;
   TH2F *hGenTrkPtEta;
   
   TH3F *hTrkPtEtaJetEt;
   TH3F *hTrkPtEtaJetEtW;
   
   TH3F *hTrkPtEtaJetEt_vbin;
   TH3F *hTrkPtEtaJetEtW_vbin;
   TH3F *hTrkPtEtaJetEtMult0to3_vbin;

   std::vector<TH3F*> hTrkPtEtaJetEt_Trig;
   
   TH3F *hTrkPtEtaJetEtW_mult1;
   TH3F *hTrkPtEtaJetEtW_mult2;
   TH3F *hTrkPtEtaJetEtW_mult3;
   
   TH3F *hGenTrkPtEtaJetEt;
   TH3F *hGenTrkPtEtaJetEtW;

   TH3F *hGenTrkPtEtaJetEt_vbin;
   TH3F *hGenTrkPtEtaJetEtW_vbin;

   //evt sel eff. correction
   TF1 *evtSelEff;

   // Jets
   TH1F *hNumJets;
   TH1F *hJet0Pt;
   TH1F *hJet0Eta;
   TH2F *hJet0PtNvtx;
   TH3F *hJet0EtaPhiEt;
   std::vector<TH1F*> hJet0Pt_Trig;

   edm::Service<TFileService> fs;
   
   bool isGEN_, doJet_, pureGENmode_, nsdOnly_, pythia6_;
   bool histOnly_, includeExtra_;
   bool applyEvtEffCorr_;
   bool setQualityBit_;
   bool mode900GeV_;
   bool trkAcceptedJet_;
   double pthatCut_;

   int32_t evtEffCorrType_, efit_type_, evtMultCut_;

   double  etaMax_, ptMin_;
   double  leadJetEt_, leadJetEta_, leadGJetEt_, leadGJetEta_;

   std::vector<bool> hltAccept_;
   std::vector<double> evtSelEffv_;

   edm::InputTag src_;
   edm::InputTag vsrc_;
   edm::InputTag jsrc_;
   edm::InputTag gsrc_;
   edm::InputTag gjsrc_;
   edm::InputTag src_evtCorr_;
   std::vector<std::string> hltNames_;
   edm::InputTag triglabel_;
   
   std::vector<int32_t> neededTrigSpectra_;
   std::vector<double> etaBins, ptBins, jetBins;
      
};
