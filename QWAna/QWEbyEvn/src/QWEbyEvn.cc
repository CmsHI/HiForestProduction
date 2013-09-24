// -*- C++ -*-
//
// Package:    QWEbyEvn
// Class:      QWEbyEvn
// 
/**\class QWEbyEvn QWEbyEvn.cc QWAna/QWEbyEvn/src/QWEbyEvn.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Quan Wang
//         Created:  Tue Apr 23 04:24:47 EDT 2013
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "QWAna/QWSpectra/interface/QWConst.h"

#include "QWAna/QWEbyEvn/interface/acc_mean_eff.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TRandom3.h"
//
// class declaration
//

class QWEbyEvn : public edm::EDAnalyzer {
	public:
		explicit QWEbyEvn(const edm::ParameterSet&);
		~QWEbyEvn();

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

		double minvz_, maxvz_;
		double dzerr_;
		double chi2_;
		double fulldzerr_;
		double fulld0err_;
		double pterr_;
		int nhits_;
		bool fillntuple_;
		double ptmin_;
		double ptmax_;
		double eta_;
		bool eff_;
		bool acc_;

		int	UseSim_;
		double	simV2_;
		double	simV3_;
		double	SimulatePhi();

		TRandom3 * gRandom;

		TH2D * hVn2Dfull[7][NCENT_PBPB]; 
		TH2D * hVn2Dsub0[7][NCENT_PBPB]; 
		TH2D * hVn2Dsub1[7][NCENT_PBPB]; 
		TH2D * hVn2D0v1[7][NCENT_PBPB]; 

		TH1D * hVnFull[7][NCENT_PBPB]; 
		TH1D * hVnSub0[7][NCENT_PBPB]; 
		TH1D * hVnSub1[7][NCENT_PBPB]; 

		TNtupleD * hNtVn[7][NCENT_PBPB];

		TFile * fnt;
		TFile * feff;

		TH2D * hTot_eff[5];
		int effBin(int);
		bool 	TrackNoffCut(reco::TrackCollection::const_iterator, const reco::VertexCollection *);
		bool 	TrackQualityCut(reco::TrackCollection::const_iterator, const reco::VertexCollection *);
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
QWEbyEvn::QWEbyEvn(const edm::ParameterSet& iConfig)
	:
		tracks_(iConfig.getUntrackedParameter<edm::InputTag>("tracks_"))

{
	centrality_ = iConfig.getParameter<edm::InputTag>("centrality_");
	vtxCollection_ = iConfig.getParameter<edm::InputTag>("vtxCollection_");
	minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
	dzerr_ = iConfig.getUntrackedParameter<double>("dzerr_", 10.);
	chi2_ = iConfig.getUntrackedParameter<double>("chi2_", 40.);
	fulldzerr_ = iConfig.getUntrackedParameter<double>("fulldzerr_", 3.);
	fulld0err_ = iConfig.getUntrackedParameter<double>("fulld0err_", 3.);
	pterr_ = iConfig.getUntrackedParameter<double>("pterr_", 0.05);
	nhits_ = iConfig.getUntrackedParameter<int>("nhits_", 12);
	fillntuple_ = iConfig.getUntrackedParameter<bool>("fillntuple_", false);
	ptmin_ = iConfig.getUntrackedParameter<double>("ptmin_", 0.5);
	ptmax_ = iConfig.getUntrackedParameter<double>("ptmax_", 100.);
	eta_ = iConfig.getUntrackedParameter<double>("eta_", 2.4);
	eff_ = iConfig.getUntrackedParameter<bool>("eff_", false);
	acc_ = iConfig.getUntrackedParameter<bool>("acc_", false);

	UseSim_ = iConfig.getUntrackedParameter<int>("UseSim_", 0);
	simV2_ = iConfig.getUntrackedParameter<double>("simV2_", 0.);
	simV3_ = iConfig.getUntrackedParameter<double>("simV3_", 0.);

	if ( UseSim_ ) gRandom = new TRandom3();
	TH1D::SetDefaultSumw2();

	if ( fillntuple_ ) fnt = new TFile("nt.root", "recreate");
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;
	for ( int c = 0; c < NCENT_PBPB; c++ ) {
		for ( int n = 1; n < 7; n++ ) {
			hVn2Dfull[n][c] = fs->make<TH2D>(Form("hVn2Dfull_%i_%i", n, c), Form("hVn2Dfull_%i_%i", n, c), 150, -0.6, 0.6, 150, -0.6, 0.6);
			hVn2Dsub0[n][c] = fs->make<TH2D>(Form("hVn2Dsub0_%i_%i", n, c), Form("hVn2Dsub0_%i_%i", n, c), 150, -0.6, 0.6, 150, -0.6, 0.6);
			hVn2Dsub1[n][c] = fs->make<TH2D>(Form("hVn2Dsub1_%i_%i", n, c), Form("hVn2Dsub1_%i_%i", n, c), 150, -0.6, 0.6, 150, -0.6, 0.6);
			hVn2D0v1[n][c] = fs->make<TH2D>(Form("hVn2D0v1_%i_%i", n, c), Form("hVn2D0v1_%i_%i", n, c), 150, -0.6, 0.6, 150, -0.6, 0.6);
			hVnFull[n][c] = fs->make<TH1D>(Form("hVnFull_%i_%i", n, c), Form("hVnFull_%i_%i", n, c), 150, 0., 0.6);
			hVnSub0[n][c] = fs->make<TH1D>(Form("hVnSub0_%i_%i", n, c), Form("hVnSub0_%i_%i", n, c), 150, 0., 0.6);
			hVnSub1[n][c] = fs->make<TH1D>(Form("hVnSub1_%i_%i", n, c), Form("hVnSub1_%i_%i", n, c), 150, 0., 0.6);
			if (fillntuple_) {
				fnt->cd();
				hNtVn[n][c] = new TNtupleD(Form("hNtVn_%i_%i", n, c), Form("hNtVn_%i_%i", n, c), "Noff:cfull:sfull:csub0:ssub0:csub1:ssub1:mult_full:mult_sub0:mult_sub1");
			}
		}
	}

	if ( eff_ ) {
		feff = new TFile("../trkEffNew2012_HI_hiGoodTightMerged_xsec_smoothv5true.root");
		if ( feff->IsOpen() ) {
			hTot_eff[0] = (TH2D*) feff->Get("Tot_0");
			hTot_eff[1] = (TH2D*) feff->Get("Tot_1");
			hTot_eff[2] = (TH2D*) feff->Get("Tot_2");
			hTot_eff[3] = (TH2D*) feff->Get("Tot_3");
			hTot_eff[4] = (TH2D*) feff->Get("Tot_4");
		}
	}
}


QWEbyEvn::~QWEbyEvn()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

	if ( fillntuple_ ) {
		fnt->Write();
		fnt->Close();
	}
}


//
// member functions
//

bool
QWEbyEvn::TrackQualityCut(reco::TrackCollection::const_iterator itTrack, const reco::VertexCollection * recoVertices)
{
	if ( itTrack->charge()==0 ) return false;
	if ( itTrack->pt() > ptmax_ or itTrack->pt() < ptmin_ ) return false;
	if ( fabs(itTrack->eta()) > eta_ ) return false;

	int primaryvtx = 0;
	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
	double vxError = (*recoVertices)[primaryvtx].xError();
	double vyError = (*recoVertices)[primaryvtx].yError();
	double vzError = (*recoVertices)[primaryvtx].zError();

	double d0 = -1.* itTrack->dxy(v1);
	double d0error=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
	double dz=itTrack->dz(v1);
	double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);

	if ( itTrack->numberOfValidHits() < 7 ) {
		// pixel track
		// dz significance cut
		if ( fabs(dz/dzerror) > dzerr_ ) return false;
		// chi2/ndof cut
		if ( itTrack->normalizedChi2() > chi2_ ) return false;
	} else {
		// full track
		// dz and d0 significance cuts
		if ( fabs(dz/dzerror) > fulldzerr_ ) return false;
		if ( fabs(d0/d0error) > fulld0err_ ) return false;

		// pt resolution cut
		if ( itTrack->ptError()/itTrack->pt() > pterr_ ) return false;
		// number of valid hits cut
		if ( itTrack->numberOfValidHits() < nhits_ ) return false;
	}

	return true;
}

bool
QWEbyEvn::TrackNoffCut(reco::TrackCollection::const_iterator itTrack, const reco::VertexCollection * recoVertices)
{
	if ( itTrack->charge()==0 ) return false;
	if ( itTrack->pt() < 0.4 ) return false;
	if ( fabs(itTrack->eta()) > 2.4 ) return false;

	int primaryvtx = 0;
	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
	double vxError = (*recoVertices)[primaryvtx].xError();
	double vyError = (*recoVertices)[primaryvtx].yError();
	double vzError = (*recoVertices)[primaryvtx].zError();

	double d0 = -1.* itTrack->dxy(v1);
	double d0error=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
	double dz=itTrack->dz(v1);
	double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);

	if ( itTrack->numberOfValidHits() < 7 ) {
		// pixel track
		// dz significance cut
		if ( fabs(dz/dzerror) > dzerr_ ) return false;
		// chi2/ndof cut
		if ( itTrack->normalizedChi2() > chi2_ ) return false;
	} else {
		// full track
		// dz and d0 significance cuts
		if ( fabs(dz/dzerror) > fulldzerr_ ) return false;
		if ( fabs(d0/d0error) > fulld0err_ ) return false;

		// pt resolution cut
		if ( itTrack->ptError()/itTrack->pt() > pterr_ ) return false;
		// number of valid hits cut
		if ( itTrack->numberOfValidHits() < nhits_ ) return false;
	}

	return true;
}


// ------------ method called for each event  ------------
	void
QWEbyEvn::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
// get centrality
	edm::Handle<int> ch;
	iEvent.getByLabel(centrality_,ch);
	int bin = *(ch.product());
	if ( bin >= NCENT_PBPB ) {
		std::cout << "!! bin = " << bin << std::endl;
		return;
	}

// vertex
	edm::Handle<reco::VertexCollection> vtx;
	iEvent.getByLabel(vtxCollection_, vtx);
	const reco::VertexCollection * recoVertices = vtx.product();
	int primaryvtx = 0;
//	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
//	double vxError = (*recoVertices)[primaryvtx].xError();
//	double vyError = (*recoVertices)[primaryvtx].yError();
//	double vzError = (*recoVertices)[primaryvtx].zError();

//	double vx = (*recoVertices)[primaryvtx].x();
//	double vy = (*recoVertices)[primaryvtx].y();
	double vz = (*recoVertices)[primaryvtx].z();
	if (vz < minvz_ || vz > maxvz_) {
//		return;
	}

// track
	edm::Handle<reco::TrackCollection> tracks;
	iEvent.getByLabel(tracks_,tracks);

	double cfull[7] = { 0,0,0,0,0,0,0 };
	double sfull[7] = { 0,0,0,0,0,0,0 };
	double csub0[7] = { 0,0,0,0,0,0,0 };
	double ssub0[7] = { 0,0,0,0,0,0,0 };
	double csub1[7] = { 0,0,0,0,0,0,0 };
	double ssub1[7] = { 0,0,0,0,0,0,0 };
	double mult_full=0;
	double mult_sub0=0;
	double mult_sub1=0;
	int Noff = 0;
	for ( reco::TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {

		if ( TrackNoffCut(itTrack, recoVertices) ) Noff++;
		if ( !TrackQualityCut(itTrack, recoVertices) ) continue;
		double phi = 0;
		if ( UseSim_ ) {
			phi = SimulatePhi();
		} else {
			phi = itTrack->phi();
		}

		double eff = 1.;

		if ( eff_ ) {
			int effb = effBin(bin);
			int binx = hTot_eff[effb]->GetXaxis()->FindBin(itTrack->eta());
			int biny = hTot_eff[effb]->GetYaxis()->FindBin(itTrack->pt());
			eff = 1./hTot_eff[effb]->GetBinContent(binx, biny);
		}
		for ( int n = 1; n < 7; n++ ) {
			cfull[n] += TMath::Cos(n*phi)*eff;
			sfull[n] += TMath::Sin(n*phi)*eff;
			if ( itTrack->eta() < 0 ) {
				csub0[n] += TMath::Cos(n*phi)*eff;
				ssub0[n] += TMath::Sin(n*phi)*eff;
			} else {
				csub1[n] += TMath::Cos(n*phi)*eff;
				ssub1[n] += TMath::Sin(n*phi)*eff;
			}
		}
		mult_full += eff;
		if ( itTrack->eta() < 0 ) {
			mult_sub0 += eff;
		} else {
			mult_sub1 += eff;
		}
	}
//	if ( (mult_full==0) or (mult_sub0==0) or (mult_sub1==0) ) return;
	for ( int n = 1; n < 7; n++ ) {
		cfull[n] /= mult_full;
		sfull[n] /= mult_full;
		csub0[n] /= mult_sub0;
		ssub0[n] /= mult_sub0;
		csub1[n] /= mult_sub1;
		ssub1[n] /= mult_sub1;

		if ( acc_ ) {
			cfull[n] -= m_cfull[n][bin];
			sfull[n] -= m_sfull[n][bin];
			csub0[n] -= m_csub0[n][bin];
			ssub0[n] -= m_ssub0[n][bin];
			csub1[n] -= m_csub1[n][bin];
			ssub1[n] -= m_ssub1[n][bin];
		}
		hVn2Dfull[n][bin]->Fill(cfull[n], sfull[n]);
		hVn2Dsub0[n][bin]->Fill(csub0[n], ssub0[n]);
		hVn2Dsub1[n][bin]->Fill(csub1[n], ssub1[n]);
		hVn2D0v1[n][bin]->Fill(csub0[n]-csub1[n], ssub0[n]-ssub1[n]);
		hVnFull[n][bin]->Fill(sqrt(cfull[n]*cfull[n] + sfull[n]*sfull[n]));
		hVnSub0[n][bin]->Fill(sqrt(csub0[n]*csub0[n] + ssub0[n]*ssub0[n]));
		hVnSub1[n][bin]->Fill(sqrt(csub1[n]*csub1[n] + ssub1[n]*ssub1[n]));
		if (fillntuple_) hNtVn[n][bin]->Fill(Noff, cfull[n], sfull[n], csub0[n], ssub0[n], csub1[n], ssub1[n], mult_full, mult_sub0, mult_sub1);
	}

}


int
QWEbyEvn::effBin(int bin)
{
	if ( bin < 2 ) return 0;
	if ( bin < 4 ) return 1;
	if ( bin < 12 ) return 2;
	if ( bin < 20 ) return 3;
	return 4;
}

double
QWEbyEvn::SimulatePhi()
{
	double limit = 1 + 2*simV2_ + 2*simV3_;
	double phi=0;
	do {
		phi = gRandom->Uniform(-Pi, Pi);
	} while ( gRandom->Uniform(limit) > 1 + 2*simV2_*TMath::Cos(2*phi) + simV3_*TMath::Cos(3*phi) );
	return phi;
}

// ------------ method called once each job just before starting event loop  ------------
	void 
QWEbyEvn::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
QWEbyEvn::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
QWEbyEvn::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
QWEbyEvn::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
QWEbyEvn::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
QWEbyEvn::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QWEbyEvn::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(QWEbyEvn);
