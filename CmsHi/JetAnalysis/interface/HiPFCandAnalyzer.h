#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <vector>

const Int_t MAXPARTICLE = 10000;
//
// DiJet ana Event Data Tree definition
//
class TreePFCandEventData
{
 public:
  // ===== Class Methods =====
  void SetDefaults();
  TreePFCandEventData();
  void SetTree(TTree * t) { tree_=t; }
  void SetBranches();
  void Clear();
  bool doJets;
  bool doMC;

  Float_t                 jdphi_;
  // -- particle info --
  Int_t                   nPFpart_, nGENpart_, njets_;
  Int_t                   pfId_[MAXPARTICLE], genPDGId_[MAXPARTICLE];
  Float_t                 pfPt_[MAXPARTICLE], genPt_[MAXPARTICLE],  jetPt_[MAXPARTICLE];
  Float_t                 pfEta_[MAXPARTICLE], genEta_[MAXPARTICLE],  jetEta_[MAXPARTICLE];
  Float_t                 pfPhi_[MAXPARTICLE], genPhi_[MAXPARTICLE],  jetPhi_[MAXPARTICLE];

 private:
  TTree*                 tree_;
};

class HiPFCandAnalyzer : public edm::EDAnalyzer {
  public:
    explicit HiPFCandAnalyzer(const edm::ParameterSet&);
    ~HiPFCandAnalyzer();

    // class methods


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::Service<TFileService> fs;

    // === Ana setup ===

    // Event Info
    edm::InputTag pfCandidateLabel_;
    edm::InputTag genLabel_;
    edm::InputTag jetLabel_;

    TTree	  *pfTree_;
    TreePFCandEventData pfEvt_;

    // cuts
    Double_t        pfPtMin_;
    Double_t        jetPtMin_;
    Double_t        genPtMin_;

   // debug
    Int_t	  verbosity_;

    bool   doJets_;
    bool   doMC_;
    bool   skipCharged_;
};


