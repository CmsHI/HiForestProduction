#include "CmsHi/PhotonAnalysis/interface/HiPhotonMCTruth.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"


#include <vector>

using namespace edm;
using namespace reco;

#define PI 3.141592653589793238462643383279502884197169399375105820974945

HiPhotonMCTruth::HiPhotonMCTruth(edm::Handle<GenParticleCollection> inputHandle)
{
  using namespace std;

  const GenParticleCollection *collection1 = inputHandle.product();
   mcisocut = HiMCGammaJetSignalDef(collection1);
  
} 


bool HiPhotonMCTruth::IsIsolated(const reco::GenParticle &pp)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
   // Check if a given particle is isolated.

  return  mcisocut.IsIsolated(pp);
}

bool HiPhotonMCTruth::IsIsolatedPP(const reco::GenParticle &pp, double cone, double etCut)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  // Check if a given particle is isolated.                                                                                                                                                          

  return  mcisocut.IsIsolatedPP(pp,cone,etCut);
}

bool HiPhotonMCTruth::IsIsolatedJP(const reco::GenParticle &pp)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   // Check if a given particle is isolated.                                                                                                                                                                  

   return  mcisocut.IsIsolatedJP(pp);
}





bool HiPhotonMCTruth::IsPrompt(const reco::GenParticle &pp)
{
  using namespace std;
  if ( pp.pdgId() != 22)
    return false;

  if ( pp.mother() ==0 ) 
    {
      //      cout <<    "No mom for this particle.." << endl;
      return false;
    }
  else 
    {
      if (pp.mother()->pdgId() == 22)
	{
	  cout << " found a prompt photon" << endl;
	  return true;
	}
      else
	return false;
    }
  
  
}
