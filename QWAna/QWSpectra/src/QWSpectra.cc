// -*- C++ -*-
//
// Package:    QWSpectra
// Class:      QWSpectra
// 
/**\class QWSpectra QWSpectra.cc QWAna/QWSpectra/src/QWSpectra.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Quan Wang
//         Created:  Tue Apr 23 10:56:56 EDT 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "QWAna/QWSpectra/interface/QWConst.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
//
// class declaration
//

class QWSpectra : public edm::EDAnalyzer {
	public:
		explicit QWSpectra(const edm::ParameterSet&);
		~QWSpectra();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		// ----------member data ---------------------------
		edm::InputTag tracks_; //used to select what tracks to read from configuration file
		edm::InputTag centrality_; // centralityBin
		edm::InputTag vtxCollection_; // centralityBin
		edm::InputTag eventplane_; // ep source
		double minvz_, maxvz_;
		double dzerr_;
		double chi2_;
		double fulldzerr_;
		double fulld0err_;
		double pterr_;
		int nhits_;

		TH1D* hPt[NCENT_PBPB];
		TH2D* hPhiEtaPos[NCENT_PBPB];
		TH2D* hPhiEtaNeg[NCENT_PBPB];
		TH1D* hCent;
		TH1D* hVtxZ[NCENT_PBPB];
		TH2D* hVtxXY[NCENT_PBPB];
		TH1D* hMult[NCENT_PBPB];
		TH1D* hNtrk[NCENT_PBPB];

		TH1D* hEP[hi::NumEPNames][NCENT_PBPB];

		edm::Service<TFileService> fs;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
QWSpectra::QWSpectra(const edm::ParameterSet& iConfig)
	:
		tracks_(iConfig.getUntrackedParameter<edm::InputTag>("tracks_"))

{
	centrality_ = iConfig.getParameter<edm::InputTag>("centrality_");
	vtxCollection_ = iConfig.getParameter<edm::InputTag>("vtxCollection_");
	eventplane_ = iConfig.getParameter<edm::InputTag>("eventplane_");
	minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
	dzerr_ = iConfig.getUntrackedParameter<double>("dzerr_", 10.);
	chi2_ = iConfig.getUntrackedParameter<double>("chi2_", 40.);
	fulldzerr_ = iConfig.getUntrackedParameter<double>("fulldzerr_", 3.);
	fulld0err_ = iConfig.getUntrackedParameter<double>("fulld0err_", 3.);
	pterr_ = iConfig.getUntrackedParameter<double>("pterr_", 0.05);
	nhits_ = iConfig.getUntrackedParameter<int>("pterr_", 12);

	TH1D::SetDefaultSumw2();
	//now do what ever initialization is needed
	for ( int i = 0; i < NCENT_PBPB; i++ ) {
		hPt[i] = fs->make<TH1D>(Form("hPt_%i",i), Form("hPt_%i",i), 1000, 0., 20.);
		hPhiEtaPos[i] = fs->make<TH2D>(Form("hPhiEtaPos_%i", i), Form("hPhiEtaPos_%i", i), 200, -Pi, Pi, 200, -2.5, 2.5);
		hPhiEtaNeg[i] = fs->make<TH2D>(Form("hPhiEtaNeg_%i", i), Form("hPhiEtaNeg_%i", i), 200, -Pi, Pi, 200, -2.5, 2.5);
		hVtxZ[i] = fs->make<TH1D>(Form("hVtxZ_%i", i), Form("hVtxZ_%i", i), 100, -25., 25.);
		hVtxXY[i] = fs->make<TH2D>(Form("hVtxXY_%i", i), Form("hVtxXY_%i", i), 100, -.1, .1, 100, -.1, .1);
		hMult[i] = fs->make<TH1D>(Form("hMult_%i", i), Form("hMult_%i", i), 1500, 0, 1500);
		hNtrk[i] = fs->make<TH1D>(Form("hNtrk_%i", i), Form("hNtrk_%i", i), 1500, 0, 1500);
	}
	hCent = fs->make<TH1D>("hCent", "hCent", 40, 0, 40);

	for ( int i = 0; i < hi::NumEPNames; i++ ) {
		TFileDirectory fep = fs->mkdir(hi::EPNames[i].c_str());
		for ( int j = 0; j < NCENT_PBPB; j++ ) {
			hEP[i][j] = fep.make<TH1D>(Form("hEP_%i",j), Form("hEP_%i",j), 200, -Pi/hi::EPOrder[i], Pi/hi::EPOrder[i]);
		}
	}
}


QWSpectra::~QWSpectra()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
QWSpectra::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

// get centrality
	edm::Handle<int> ch;
	iEvent.getByLabel(centrality_,ch);
	int bin = *(ch.product());
	hCent->Fill(bin);
	if ( bin >= NCENT_PBPB ) {
		std::cout << "!! bin = " << bin << std::endl;
		return;
	}

//	edm::Handle<TrackCollection> tracks;
//	iEvent.getByLabel(trackTags_,tracks);
//	for(TrackCollection::const_iterator itTrack = tracks->begin();
//			itTrack != tracks->end();                      
//			++itTrack) {
//		int charge = 0;
//		charge = itTrack->charge();  
//		histo->Fill( charge );
//	}

// vertex
	edm::Handle<reco::VertexCollection> vtx;
	iEvent.getByLabel(vtxCollection_, vtx);
	const reco::VertexCollection * recoVertices = vtx.product();
	int primaryvtx = 0;
	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
	double vxError = (*recoVertices)[primaryvtx].xError();
	double vyError = (*recoVertices)[primaryvtx].yError();
	double vzError = (*recoVertices)[primaryvtx].zError();

	double vx = (*recoVertices)[primaryvtx].x();
	double vy = (*recoVertices)[primaryvtx].y();
	double vz = (*recoVertices)[primaryvtx].z();
	hVtxZ[bin]->Fill(vz);
	hVtxXY[bin]->Fill(vx, vy);
	if (vz < minvz_ || vz > maxvz_) {
		return;
	}

// track
	edm::Handle<reco::TrackCollection> tracks;
	iEvent.getByLabel(tracks_,tracks);
	int mult = 0;
	hNtrk[bin]->Fill(tracks->size());
	for ( reco::TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {

		if ( itTrack->charge()==0 ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double d0error=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);

		if ( fabs(itTrack->eta()) > 2.4 ) continue;

		if ( itTrack->numberOfValidHits() < 7 ) {
			// pixel track
			// dz significance cut
			if ( fabs(dz/dzerror) > dzerr_ ) continue;
			// chi2/ndof cut
			if ( itTrack->normalizedChi2() > chi2_ ) continue;
		} else {
			// full track
			// dz and d0 significance cuts
			if ( fabs(dz/dzerror) > fulldzerr_ ) continue;
			if ( fabs(d0/d0error) > fulld0err_ ) continue;

			// pt resolution cut
			if ( itTrack->ptError()/itTrack->pt() > pterr_ ) continue;
			// number of valid hits cut
			if ( itTrack->numberOfValidHits() < nhits_ ) continue;
		}
		hPt[bin]->Fill(itTrack->pt());
		if ( itTrack->charge() > 0 ) {
			hPhiEtaPos[bin]->Fill(itTrack->phi(), itTrack->eta());
		} else {
			hPhiEtaNeg[bin]->Fill(itTrack->phi(), itTrack->eta());
		}
		mult++;
	}
	hMult[bin]->Fill(mult);
// EP
	edm::Handle<reco::EvtPlaneCollection> evtPlanes;
	iEvent.getByLabel(eventplane_, evtPlanes);
	if ( !evtPlanes.isValid() ) {
		return;
	}
	int epidx = 0;
	for ( reco::EvtPlaneCollection::const_iterator rp = evtPlanes->begin();rp !=evtPlanes->end(); rp++) {
		hEP[epidx][bin]->Fill(rp->angle());
		epidx++;
	}
}


// ------------ method called once each job just before starting event loop  ------------
	void 
QWSpectra::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
QWSpectra::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
QWSpectra::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
QWSpectra::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
QWSpectra::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
QWSpectra::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QWSpectra::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(QWSpectra);
