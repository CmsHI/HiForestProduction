#ifndef MNguyen_pi0Analyzer_inclusiveJetAnalyzer_
#define MNguyen_pi0Analyzer_inclusiveJetAnalyzer_

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"



#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//

/**\class pi0Analyzer 

\author Matt Nguyen
\date   November 2010
*/




class pi0Analyzer : public edm::EDAnalyzer {
 public:

  explicit pi0Analyzer(const edm::ParameterSet&);

  ~pi0Analyzer();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

  virtual void beginJob();

  

 private:
  
  bool isMerged( const reco::PFCandidate& pfc ) const;

  edm::InputTag   pfCandidatesTag_;


  /// verbose ?
  bool   verbose_;

  bool useCentrality_;
  bool isMC_;


  TTree *t;
  edm::Service<TFileService> fs1;

  CentralityProvider * centrality_;



  static const int MAXPI0S = 9999999;



  struct Pi0{
    
    int run;
    int evt;
    int lumi;
    int bin;
    float vx, vy, vz;
    float b;
    float hf;

    int npi0s;
    float invMass[MAXPI0S], openAng[MAXPI0S];
    float pi0_pt[MAXPI0S], pi0_px[MAXPI0S], pi0_py[MAXPI0S], pi0_pz[MAXPI0S], pi0_E[MAXPI0S], pi0_eta[MAXPI0S];
    float gam1_pt[MAXPI0S], gam1_px[MAXPI0S], gam1_py[MAXPI0S], gam1_pz[MAXPI0S], gam1_E[MAXPI0S], gam1_eta[MAXPI0S];
    float gam2_pt[MAXPI0S], gam2_px[MAXPI0S], gam2_py[MAXPI0S], gam2_pz[MAXPI0S], gam2_E[MAXPI0S], gam2_eta[MAXPI0S];

    int gam1_merged[MAXPI0S], gam2_merged[MAXPI0S];



  };

  Pi0 pi0s_;

};

#endif
