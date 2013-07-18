#ifndef RecoParticleFlow_PFPatProducer_PFCandidateAnalyzer_
#define RecoParticleFlow_PFPatProducer_PFCandidateAnalyzer_

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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//EAW
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//

/**\class PFCandidateAnalyzer 
\brief produces IsolatedPFCandidates from PFCandidates

\author Colin Bernet
\date   february 2008
*/




class PFCandidateAnalyzer : public edm::EDAnalyzer {
 public:

  explicit PFCandidateAnalyzer(const edm::ParameterSet&);

  ~PFCandidateAnalyzer();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

  virtual void beginJob();

 private:
  

  void printElementsInBlocks(const reco::PFCandidate& cand,
			     std::ostream& out=std::cout) const;


  
  /// PFCandidates in which we'll look for pile up particles 
  edm::InputTag   inputTagPFCandidates_;

  /// Vertices for track selection
  edm::InputTag   inputTagVertices_;

  /// SimTracks for track association
  edm::InputTag   inputTagSimTracks_;

  /// Tracks for track association
  edm::InputTag   inputTagTracks_;
  
  /// verbose ?
  bool   verbose_;

  /// contains tracking particles
  bool   hasSimInfo_;

  /// print the blocks associated to a given candidate ?
  bool   printBlocks_;

  double thePtMin_;

  double minHits_, maxPtErr_, maxD0_, maxDZ_, maxD0Norm_, maxDZNorm_;
  bool pixelSeedOnly_;

  TNtuple *nt;
  edm::Service<TFileService> f;

  float cand_type;
  float cand_pt;
  float max_trk;
  float sum_trk;
  float max_ecal;
  float sum_ecal;
  float max_hcal;
  float sum_hcal;

  float max_nhits,max_relpterr,max_algo,max_nd0,max_ndz,max_fake;

  //vector<pair<float,float>> elements; // type, pt 

};

#endif
