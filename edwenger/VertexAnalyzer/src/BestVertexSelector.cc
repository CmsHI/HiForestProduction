#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "edwenger/VertexAnalyzer/interface/VertexComparator.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace reco
{
	namespace modules
	{
		
		// define your producer name
		typedef ObjectSelector<
			SortCollectionSelector<
			reco::VertexCollection,
			MoreTracksThenLowerChi2<reco::Vertex>
			> 
			> BestVertexSelection;
		
		// declare the module as plugin
		DEFINE_FWK_MODULE( BestVertexSelection );
	}
}
