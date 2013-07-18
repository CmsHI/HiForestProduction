#ifndef RecoAlgos_HackedAnalyticalTrackSelector_h
#define RecoAlgos_HackedAnalyticalTrackSelector_h
/** \class HackedAnalyticalTrackSelector
 *
 * ***modified to do significance cuts slightly differently***
 * 
 * selects a subset of a track collection, copying extra information on demand
 * 
 * \author Paolo Azzurri, Giovanni Petrucciani 
 *
 * \version $Revision: 1.3 $
 *
 * $Id: HackedAnalyticalTrackSelector.h,v 1.3 2011/07/19 18:27:05 sungho Exp $
 *
 */

#include <utility>
#include <vector>
#include <memory>
#include <algorithm>
#include <map>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"


namespace reco { namespace modules {

    class HackedAnalyticalTrackSelector : public edm::EDProducer {
        private:
        public:
            /// constructor 
            explicit HackedAnalyticalTrackSelector( const edm::ParameterSet & cfg ) ;
            /// destructor
            virtual ~HackedAnalyticalTrackSelector() ;

        private:
            typedef math::XYZPoint Point;
            /// process one event
            void produce( edm::Event& evt, const edm::EventSetup& es ) ;
            /// return class, or -1 if rejected
            bool select (const reco::BeamSpot &vertexBeamSpot, const reco::Track &tk, const std::vector<Point> &points, std::vector<double> &vterr, std::vector<double> &vzerr);
            void selectVertices ( const reco::VertexCollection &vtxs, std::vector<Point> &points, std::vector<double> &vterr, std::vector<double> &vzerr);
            /// source collection label
            edm::InputTag src_;
            edm::InputTag beamspot_;
            edm::InputTag vertices_;
            /// copy only the tracks, not extras and rechits (for AOD)
            bool copyExtras_;
            /// copy also trajectories and trajectory->track associations
            bool copyTrajectories_;

            /// save all the tracks
            bool keepAllTracks_;
            /// do I have to set a quality bit?
            bool setQualityBit_;
            TrackBase::TrackQuality qualityToSet_;

            /// vertex cuts
            int32_t vtxNumber_;
            size_t  vtxTracks_;
            double  vtxChi2Prob_;

			//  parameters for adapted optimal cuts on chi2 and primary vertex compatibility
			std::vector<double> res_par_;
            double  chi2n_par_;
	    std::vector<double> chi2n_par_pt_;
	    std::vector<double> chi2n_par_pt2_; // harder cut for |eta|>0.8

			std::vector<double> d0_par1_;
			std::vector<double> dz_par1_;
			std::vector<double> d0_par2_;
			std::vector<double> dz_par2_;
			// Boolean indicating if adapted primary vertex compatibility cuts are to be applied.
            bool applyAdaptedPVCuts_;
			
            /// Impact parameter absolute cuts
            double max_d0_;
            double max_z0_;

            /// Cuts on numbers of layers with hits/3D hits/lost hits. 
			uint32_t min_layers_;
			uint32_t min_3Dlayers_;
			uint32_t max_lostLayers_;
			
            /// storage
            std::auto_ptr<reco::TrackCollection> selTracks_;
            std::auto_ptr<reco::TrackExtraCollection> selTrackExtras_;
            std::auto_ptr< TrackingRecHitCollection>  selHits_;
            std::auto_ptr< std::vector<Trajectory> > selTrajs_;
            std::auto_ptr< std::vector<const Trajectory *> > selTrajPtrs_;
            std::auto_ptr< TrajTrackAssociationCollection >  selTTAss_;
            reco::TrackRefProd rTracks_;
            reco::TrackExtraRefProd rTrackExtras_;
            TrackingRecHitRefProd rHits_;
            edm::RefProd< std::vector<Trajectory> > rTrajectories_;
            std::vector<reco::TrackRef> trackRefs_;

    };

} }

#endif
