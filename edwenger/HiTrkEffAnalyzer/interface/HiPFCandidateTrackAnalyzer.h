#ifndef HiTrkEffAnalyzer_HiPFCandidateTrackAnalyzer
#define HiTrkEffAnalyzer_HiPFCandidateTrackAnalyzer

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

//EAW, ASY
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// centrality
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

/**\class PFCandidateTrackAnalyzer 
\brief produces IsolatedPFCandidates from PFCandidates

\author Colin Bernet
\date   february 2008
--> Modified by Ed Wenger, Matt Nguyen, and Andre Yoon
*/




class HiPFCandidateTrackAnalyzer : public edm::EDAnalyzer {
 public:

  explicit HiPFCandidateTrackAnalyzer(const edm::ParameterSet&);

  ~HiPFCandidateTrackAnalyzer();
  
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

  std::string funcCaloComp_;

  bool applyTrkQCs_;
  
  double minHits_, maxPtErr_, maxD0_, maxDZ_, maxD0Norm_, maxDZNorm_;
  bool pixelSeedOnly_;

  std::vector<int32_t> neededCentBins_;

  CentralityProvider * centrality_;

  bool useQaulityStr_;
  std::string qualityString_;

  // ntuple
  TNtuple *nt;

  // TH1F

  // TH2F
  TH2F *hTrkPtEcalEtSum, *hTrkPtHcalEtSum, *hTrkPtCaloEtSum;
  TH2F *hTrkPtEcalEtSum_fake, *hTrkPtHcalEtSum_fake, *hTrkPtCaloEtSum_fake;
  TH2F *hTrkPtEcalEtSum_real, *hTrkPtHcalEtSum_real, *hTrkPtCaloEtSum_real;

  TH2F *hHitsEtaAccept, *hHitsEtaReject;
  TH2F *hPtErrEtaAccept, *hPtErrEtaReject;
  TH2F *hAlgoEtaAccept, *hAlgoEtaReject;
  TH2F *hD0EtaAccept, *hD0EtaReject;
  TH2F *hDZEtaAccept, *hDZEtaReject;
  TH2F *hD0ErrEtaAccept, *hD0ErrEtaReject; // quadratic sum of d0 err and xy err w.r.t vtx
  TH2F *hDZErrEtaAccept, *hDZErrEtaReject; // quadratic sum of dz err and z err w.r.t vtx 
  TH2F *hD0PerErrEtaAccept, *hD0PerErrEtaReject;
  TH2F *hDZPerErrEtaAccept, *hDZPerErrEtaReject;

  // TH2F binned in centrality
  std::vector<TH2F*> hTrkPtEcalEtSum_Cent;
  std::vector<TH2F*> hTrkPtHcalEtSum_Cent;
  std::vector<TH2F*> hTrkPtCaloEtSum_Cent;

  // TH3F
  TH3F *hTrkPtEtaEcalEtSum, *hTrkPtEtaHcalEtSum, *hTrkPtEtaCaloEtSum;
  TH3F *hTrkPtEtaEcalEtSum_fake, *hTrkPtEtaHcalEtSum_fake, *hTrkPtEtaCaloEtSum_fake;
  TH3F *hTrkPtEtaEcalEtSum_real, *hTrkPtEtaHcalEtSum_real, *hTrkPtEtaCaloEtSum_real;

  // TF1
  TF1 *fCaloComp;

  edm::Service<TFileService> fs;

  float cand_type;
  float cand_pt;
  float cand_eta;
  float max_trk;
  float sum_trk;
  float max_ecal;
  float sum_ecal;
  float max_hcal;
  float sum_hcal;
  float sum_calo;
  float max_calo;

  float max_nhits,max_relpterr,max_algo,max_nd0,max_ndz,max_fake;

  bool prodNtuple_;
  bool isData_;


  std::vector<double> etaBins, ptBins, cEtSumBins, nhitBins;

  //vector<pair<float,float>> elements; // type, pt 

};

#endif
