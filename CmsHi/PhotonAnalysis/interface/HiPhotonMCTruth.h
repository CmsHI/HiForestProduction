#ifndef HiPhotonMCTruth_h
#define HiPhotonMCTruth_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "CmsHi/PhotonAnalysis/interface/HiMCGammaJetSignalDef.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include <vector>

class HiPhotonMCTruth
{
 public:
  HiPhotonMCTruth(edm::Handle<reco::GenParticleCollection> inputHandle);
  bool IsPrompt(const reco::GenParticle &pp);
  bool IsIsolated(const reco::GenParticle &pp);
  bool IsIsolatedPP(const reco::GenParticle &pp,double cone = 0.4, double etCut = 2.0);
  bool IsIsolatedJP(const reco::GenParticle &pp);
 
 private:
  HiMCGammaJetSignalDef mcisocut;
};

#endif

