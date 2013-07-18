#ifndef HiMCGammaJetSignalDef_h
#define HiMCGammaJetSignalDef_h

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

#include "FWCore/Framework/interface/ESHandle.h"
#include <vector>

class HiMCGammaJetSignalDef
{
  public:

  HiMCGammaJetSignalDef();
  HiMCGammaJetSignalDef(const reco::GenParticleCollection *sigPartic);
  bool IsIsolated(const reco::GenParticle &pp);  
  bool IsIsolatedPP(const reco::GenParticle &pp, double cone = 0.4, double etCut = 2.0);
  bool IsIsolatedJP(const reco::GenParticle &pp);

  //  bool IsSignal(const reco::Candidate &pp, double dPhi, bool isIso);
  //  int getIndex(const reco::Candidate &pp);
  double getDeltaR (const reco::Candidate &track1, const reco::Candidate &track2);
  double getDeltaPhi(const reco::Candidate &track1, const reco::Candidate &track2);
  
  
 private:
   const reco::GenParticleCollection        *fSigParticles;

};

#endif

