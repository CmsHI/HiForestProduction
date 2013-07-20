/** \class reco::modules::HackedTrackingParticleSelector
 *
 *  Filter to select TrackingParticles according to pt, rapidity, tip, lip, number of hits, pdgId
 *
 *  \author Giuseppe Cerati, INFN (modified to include status selection option, 2010/09/30)
 *
 *  $Date: 2010/09/30 11:39:05 $
 *  $Revision: 1.1 $
 */

#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
//#include "CommonTools/RecoAlgos/interface/HackedTrackingParticleSelector.h"
#include "edwenger/HiTrkEffAnalyzer/interface/HackedTrackingParticleSelector.h"

namespace reco {
  namespace modules {
    typedef SingleObjectSelector<TrackingParticleCollection,::HackedTrackingParticleSelector> 
    HackedTrackingParticleSelector ;

    DEFINE_FWK_MODULE( HackedTrackingParticleSelector );
  }
}
