#include "CmsHi/PhotonAnalysis/interface/HiMCGammaJetSignalDef.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include <vector>

using namespace edm;
using namespace reco;

#define PI 3.141592653589793238462643383279502884197169399375105820974945

HiMCGammaJetSignalDef::HiMCGammaJetSignalDef () :
  fSigParticles(0)
{}
HiMCGammaJetSignalDef::HiMCGammaJetSignalDef (const reco::GenParticleCollection  *sigParticles) :
  fSigParticles(0)
{
  using namespace std;

   fSigParticles =  sigParticles;
 
} 

bool HiMCGammaJetSignalDef::IsIsolated(const reco::GenParticle &pp)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
   // Check if a given particle is isolated.

  double  EtCone = 0;
  double  EtPhoton = 0;
  double  PtMaxHadron = 0;
  double  cone = 0.5;
   
   const int maxindex = (int)fSigParticles->size();
   for(int i=0; i < maxindex ; ++i) {
     
     const GenParticle &p=(*fSigParticles)[i];
     
     if(p.status()!=1) 
       continue; 
     if (p.collisionId() != pp.collisionId())
       continue;
	     
     int apid= abs(p.pdgId());
     if(apid>11 &&  apid<20) 
       continue; //get rid of muons and neutrinos
     
     if(getDeltaR(p,pp)>cone) 
       continue;
     
     double et=p.et();
     double pt = p.pt();
     EtCone+=et;
      
      bool isHadron = false;
      if ( fabs(pp.pdgId())>100 && fabs(pp.pdgId()) != 310)
	isHadron = true;
      
      if(apid>100 && apid!=310 && pt>PtMaxHadron) 
	{
	  if ( (isHadron == true) && (fabs(pp.pt()-pt)<0.001) && (pp.pdgId()==p.pdgId()) )
	    continue;
	  else
	    PtMaxHadron=pt;
	}

   }

   EtPhoton = pp.et();

   // isolation cuts
   if(EtCone-EtPhoton > 5+EtPhoton/20) 
      return kFALSE;
   if(PtMaxHadron > 4.5+EtPhoton/40) 
      return kFALSE;
     
   return kTRUE; 

}



bool HiMCGammaJetSignalDef::IsIsolatedPP(const reco::GenParticle &pp,double cone, double etCut)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

// Check if a given particle is isolated.                                                                                                                                                          

  double  EtCone = 0;
  double  EtPhoton = 0;
  
  const int maxindex = (int)fSigParticles->size();
  for(int i=0; i < maxindex ; ++i) {
    
    const GenParticle &p=(*fSigParticles)[i];
    
  if(p.status()!=1)
    continue;
  if (p.collisionId() != pp.collisionId())
    continue;

   int apid= abs(p.pdgId());
   //  if(apid>11 &&  apid<20)
   //  continue; //get rid of muons and neutrinos                                                                                                                                                     
  if(getDeltaR(p,pp)>cone)
    continue;

  double et=p.et();
  //  double pt = p.pt();
  EtCone+=et;

  }
  
  EtPhoton = pp.et();
  
  // isolation cuts                                                                                                                                                                                  
  if(EtCone-EtPhoton > etCut)
     return kFALSE;
  
  return kTRUE;
  
}

bool HiMCGammaJetSignalDef::IsIsolatedJP(const reco::GenParticle &pp)
{
   using namespace std;
   using namespace edm;
   using namespace reco;

   // Check if a given particle is isolated.                                                                                                                                                                    

   double  EtCone = 0;
   double  EtPhoton = 0;
   double  cone = 0.4;

   const int maxindex = (int)fSigParticles->size();
   for(int i=0; i < maxindex ; ++i) {

      const GenParticle &p=(*fSigParticles)[i];

      if(p.status()!=1)
	 continue;
      if (p.collisionId() != pp.collisionId())
	 continue;
      
      if(getDeltaR(p,pp)>cone)
	 continue;

      double et=p.et();
      //  double pt = p.pt();                                                                                                                                                                                    
      EtCone+=et;

   }

   EtPhoton = pp.et();

   // isolation cuts                                                                                                                                                                                          
   if(EtCone-EtPhoton > 2)
      return kFALSE;
   return kTRUE;

}








/*
  bool HiMCGammaJetSignalDef::IsSignal(const reco::Candidate &pp,                   double dPhi, bool isIso)
  {
  // Check if a given particle is considered to be signal. If isIso=true it is
  // already assumed that it is known that the given particle is isolated.
  
  if(pp.pdgId()!=22)
      return false;
      if(IsIsolated(pp) == false)
      return false;
      
      if(fiNear==-1||fiAway==-1)
      return false;
      
      const reco::Candidate &AwayParton=(*fSigParticles)[fiAway];
      
      if(fabs(getDeltaPhi(pp,AwayParton))<dPhi)
      return false;
      
      return true;
      }
*/


double HiMCGammaJetSignalDef::getDeltaR(const reco::Candidate &track1, const reco::Candidate &track2)
{
    double dEta = track1.eta() - track2.eta();
    double dPhi = track1.phi() - track2.phi();
    
    while(dPhi >= PI)       dPhi -= (2.0*PI);
    while(dPhi < (-1.0*PI)) dPhi += (2.0*PI);
    
    return sqrt(dEta*dEta+dPhi*dPhi);
}

double HiMCGammaJetSignalDef::getDeltaPhi(const reco::Candidate &track1, const reco::Candidate &track2)
{
    double dPhi = track1.phi() - track2.phi();
    
    while(dPhi >= PI)       dPhi -= (2.0*PI);
    while(dPhi < (-1.0*PI)) dPhi += (2.0*PI);
    
    return dPhi;
}
