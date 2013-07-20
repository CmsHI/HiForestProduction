//
// Original Author:  Edward Wenger
//         Created:  Fri May  7 10:33:49 CEST 2010
// $Id: EvtSelAnalyzer.cc,v 1.16 2010/08/20 13:50:25 sungho Exp $
//

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "edwenger/EvtSelAnalyzer/interface/EvtSelAnalyzer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "edwenger/VertexAnalyzer/interface/VertexComparator.h"


EvtSelAnalyzer::EvtSelAnalyzer(const edm::ParameterSet& iConfig)
:
  gtlabel_(iConfig.getUntrackedParameter<edm::InputTag>("gtlabel")),
  towerslabel_(iConfig.getUntrackedParameter<edm::InputTag>("towerslabel")),
  hfEThreshold_(iConfig.getUntrackedParameter<double>("hfEThreshold")),
  trackslabel_(iConfig.getUntrackedParameter<edm::InputTag>("trackslabel")),
  qualityString_(iConfig.getUntrackedParameter<std::string>("qualityString")),
  triglabel_(iConfig.getUntrackedParameter<edm::InputTag>("triglabel")),
  trignames_(iConfig.getUntrackedParameter<std::vector <std::string> >("trignames")),
  vtxlabel_(iConfig.getUntrackedParameter<edm::InputTag>("vtxlabel")),
  isGEN_(iConfig.getUntrackedParameter<bool>("isGEN")),
  includeSelTrk_(iConfig.getUntrackedParameter<bool>("includeSelTrk")),
  etaMaxSpec_(iConfig.getUntrackedParameter<double>("etaMaxSpec", 1.0)),
  ptMin_(iConfig.getUntrackedParameter<double>("ptMin", 0.5)),
  vtxWeight_(iConfig.getUntrackedParameter<bool>("vtxWeight")),
  pvtxG_DATA_(iConfig.getUntrackedParameter< std::vector<double> >("pvtxG_DATA")),
  pvtxG_MC_(iConfig.getUntrackedParameter< std::vector<double> >("pvtxG_MC"))
{

}

// ------------ method called to for each event  ------------
void
EvtSelAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //------- L1 trigger bits --------------------
  edm::Handle<L1GlobalTriggerReadoutRecord> gtH;
  iEvent.getByLabel(gtlabel_,gtH);
  
  // technical bits before mask
  const TechnicalTriggerWord& tword = gtH->technicalTriggerWord();
  for(int tbit=0; tbit<64; tbit++) hL1TechBits->Fill(tbit,tword.at(tbit));
    
  // algo bits before mask
  const DecisionWord& word = gtH->decisionWord();
  for(int bit=0; bit<128; bit++) hL1AlgoBits->Fill(bit,word.at(bit));

  //------- HLT -------------------------------

    edm::Handle<edm::TriggerResults> trigH;
    iEvent.getByLabel(triglabel_,trigH);
    
    // jet triggers
    const edm::TriggerNames names = iEvent.triggerNames(*trigH);
    for(unsigned i=0; i<trignames_.size(); i++) {
      unsigned index = names.triggerIndex(trignames_[i]);
      if(index < trigH->size()) {
	if(trigH->accept(index)) hHLTPaths->Fill(trignames_[i].c_str(),1);  
      } else {
	edm::LogWarning("EvtSelAnalyzer") 
	  << "Index returned by TriggerNames object for trigger '"
	  << trignames_[i]
	  << "' is out of range (" 
	  << index << " >= " << trigH->size() << ")";
      }
    }

  //------- Calotowers ------------------------
  edm::Handle<CaloTowerCollection> towersH;
  iEvent.getByLabel(towerslabel_,towersH);
  
  // HF towers above threshold
  int nHfTowersN=0, nHfTowersP=0;
  for(CaloTowerCollection::const_iterator calo = towersH->begin(); calo != towersH->end(); ++calo) {
    if(calo->energy() < hfEThreshold_) continue;
    if(calo->eta()>3) nHfTowersP++;
    if(calo->eta()<-3) nHfTowersN++;
  }
  hHfTowers->Fill(nHfTowersP,nHfTowersN);
  
  //------- Vertices ---------------------
  float Zvtx = 0;
  float dndvz_data = 0, dndvz_mc = 0;
  float Wvtx = 0;

  if(isGEN_ && vtxWeight_){
     edm::Handle<reco::VertexCollection> vtxsH;
     iEvent.getByLabel(vtxlabel_,vtxsH);
     
     reco::VertexCollection vtxs = *vtxsH;
     std::sort(vtxs.begin(),vtxs.end(),MoreTracksThenLowerChi2<reco::Vertex>());
     
     Zvtx = vtxs[0].z();
     
     dndvz_data = vtxGaussian_DATA->Eval(Zvtx);
     dndvz_mc = vtxGaussian_MC->Eval(Zvtx);
     
     if(dndvz_mc==0) Wvtx = 0;
     else Wvtx = dndvz_data/dndvz_mc;
  }

  //------- Tracks ----------------------------
  edm::Handle< std::vector <reco::Track> > tracksH;
  iEvent.getByLabel(trackslabel_,tracksH);

  // high-purity fraction and track multiplicity
  int numhighpurity=0;
  float fraction=0;
  int nREC=0, nREC_SPEC=0, nREC_STD=0, nREC_AGR=0;

  for(unsigned it=0; it<tracksH->size(); it++){
     const reco::Track & trk = (*tracksH)[it];

     if((*tracksH)[it].quality(reco::TrackBase::qualityByName(qualityString_))) numhighpurity++;
     else continue;
     if(fabs(trk.eta())>2.5) continue;
     nREC++;
     if(trk.pt()>ptMin_ && fabs(trk.eta())<2.4) nREC_STD++; //standard                                                                                                 
     if(trk.pt()>ptMin_ && fabs(trk.eta())<etaMaxSpec_) nREC_SPEC++; // spectra                                                                                        
     if(trk.pt()>0.5 && fabs(trk.eta())<0.8) nREC_AGR++; // agreement b/w experiments        
  }
  if(tracksH->size()) fraction = (float) numhighpurity / (float) tracksH->size();
  hHPFracNtrk->Fill(tracksH->size(),fraction);

  hRecMult->Fill(nREC);
  hRecMult_STD->Fill(nREC_STD);
  hRecMult_SPEC->Fill(nREC_SPEC);
  hRecMult_AGR->Fill(nREC_AGR);

  if(isGEN_ && vtxWeight_) hRecMult_STD_W->Fill(nREC_STD,Wvtx), hGenRecMultNSDvsZvtx_STD->Fill(nREC_STD,Zvtx);
  
  //----- selectTracks--------------------
  if(includeSelTrk_){
     edm::Handle< std::vector <reco::Track> > tracksH_sel;
     iEvent.getByLabel("selectTracks",tracksH_sel);

     int nREC_sel=0, nREC_SPEC_sel=0, nREC_STD_sel=0, nREC_AGR_sel=0;

     for(unsigned it=0; it<tracksH_sel->size(); it++){
	const reco::Track & strk = (*tracksH_sel)[it];

	if(!(*tracksH_sel)[it].quality(reco::TrackBase::qualityByName(qualityString_))) continue;
	if(fabs(strk.eta())>2.5) continue;
	nREC++;
	if(strk.pt()>ptMin_ && fabs(strk.eta())<2.4) nREC_STD_sel++; 
	if(strk.pt()>ptMin_ && fabs(strk.eta())<etaMaxSpec_) nREC_SPEC_sel++; 
	if(strk.pt()>0.5 && fabs(strk.eta())<0.8) nREC_AGR_sel++; 
     }

     hRecMult_sel->Fill(nREC_sel);
     hRecMult_STD_sel->Fill(nREC_STD_sel);
     hRecMult_SPEC_sel->Fill(nREC_SPEC_sel);
     hRecMult_AGR_sel->Fill(nREC_AGR_sel);
  }


  //------- GEN --------------------------------
  if(isGEN_) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    
    int nGEN=0, nGEN_SPEC=0, nGEN_STD=0, nGEN_AGR=0;
    const reco::GenParticleCollection *genCollect = genParticles.product();
    for(unsigned i=0; i<genCollect->size();i++){
      const reco::GenParticle & gen = (*genCollect)[i];
      if(gen.status() != 1) continue;
      //if(gen.collisionId() != 0) continue; // found a bug in pp MC sample, so don't use it
      if(gen.charge() == 0) continue;
      if(fabs(gen.eta())>2.5) continue;
      nGEN++;
      // extra histograms
      if(gen.pt()>ptMin_ && fabs(gen.eta())<2.4) nGEN_STD++; //standard 
      if(gen.pt()>ptMin_ && fabs(gen.eta())<etaMaxSpec_) nGEN_SPEC++; // spectra
      if(gen.pt()>0.5 && fabs(gen.eta())<0.8) nGEN_AGR++; // agreement b/w experiments
    }
    
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator", genEvtInfo);
    int pid = genEvtInfo->signalProcessID();
    edm::LogVerbatim("EvtSelAnalyzer") 
      << "Process ID = " << pid << ", GEN multiplicity = " << nGEN;
    
    hGenMultInel->Fill(nGEN);
    hGenMultInel_STD->Fill(nGEN_STD);
    hGenMultInel_SPEC->Fill(nGEN_SPEC);
    hGenMultInel_AGR->Fill(nGEN_AGR);

    hGenRecMultInel->Fill(nREC);
    hGenRecMultInel_STD->Fill(nREC_STD);
    hGenRecMultInel_SPEC->Fill(nREC_SPEC);
    hGenRecMultInel_AGR->Fill(nREC_AGR);

    switch(pid) {
    case 94:
      hGenMultDD->Fill(nGEN);
      hGenMultNSD->Fill(nGEN);

      hGenMultDD_STD->Fill(nGEN_STD);
      hGenMultNSD_STD->Fill(nGEN_STD);

      hGenMultDD_SPEC->Fill(nGEN_SPEC);
      hGenMultNSD_SPEC->Fill(nGEN_SPEC);
      
      hGenMultDD_AGR->Fill(nGEN_AGR);
      hGenMultNSD_AGR->Fill(nGEN_AGR);
      
      hGenRecMultNSD->Fill(nREC);
      hGenRecMultNSD_STD->Fill(nREC_STD);
      hGenRecMultNSD_SPEC->Fill(nREC_SPEC);
      hGenRecMultNSD_AGR->Fill(nREC_AGR);
      
      if(isGEN_ && vtxWeight_) hGenRecMultNSD_STD_W->Fill(nREC_STD,Wvtx);

      if(nREC_STD==0) hGenToRecZeroBinNSD_STD->Fill(nGEN_STD);
      if(nREC_STD<4) hGenToRec0123BinNSD_STD->Fill(nGEN_STD);
      hGenVsRecMultNSD_STD->Fill(nGEN_STD,nREC_STD);
      break;
    case 92:
    case 93:
      hGenMultSD->Fill(nGEN);
      hGenMultSD_STD->Fill(nGEN_STD);
      hGenMultSD_SPEC->Fill(nGEN_SPEC);
      hGenMultSD_AGR->Fill(nGEN_AGR);

      hGenRecMultSD->Fill(nREC);
      hGenRecMultSD_STD->Fill(nREC_STD);
      hGenRecMultSD_SPEC->Fill(nREC_SPEC);
      hGenRecMultSD_AGR->Fill(nREC_AGR);

      if(isGEN_ && vtxWeight_) hGenRecMultSD_STD_W->Fill(nREC_STD,Wvtx);
      break;
    case 11:
    case 12:
    case 13:
    case 28:
    case 53:
    case 68:
    case 95:
      hGenMultND->Fill(nGEN);
      hGenMultNSD->Fill(nGEN);

      hGenMultND_STD->Fill(nGEN_STD);
      hGenMultNSD_STD->Fill(nGEN_STD);

      hGenMultND_SPEC->Fill(nGEN_SPEC);
      hGenMultNSD_SPEC->Fill(nGEN_SPEC);
      
      hGenMultND_AGR->Fill(nGEN_AGR);
      hGenMultNSD_AGR->Fill(nGEN_AGR);

      hGenRecMultNSD->Fill(nREC);
      hGenRecMultNSD_STD->Fill(nREC_STD);
      hGenRecMultNSD_SPEC->Fill(nREC_SPEC);
      hGenRecMultNSD_AGR->Fill(nREC_AGR);

      if(nREC_STD==0) hGenToRecZeroBinNSD_STD->Fill(nGEN_STD);
      if(nREC_STD<4) hGenToRec0123BinNSD_STD->Fill(nGEN_STD);
      hGenVsRecMultNSD_STD->Fill(nGEN_STD,nREC_STD);
      break;
    default:
      edm::LogWarning("EvtSelAnalyzer") 
	<< "Is this PYTHIA?! This code doesn't know what to do for process ID = " << pid;
    }

  }
  
}
  
// ------------ method called once each job just before starting event loop  ------------
void 
EvtSelAnalyzer::beginJob()
{

  int numBins = 300;
  double xmax = 299.5;
  double xmax_STD = 299.5;
  double xmax_SPEC = 299.5;
  double xmax_AGR = 299.5;


  vtxGaussian_DATA = new TF1("vtxGaussian_DATA","[2]*TMath::Exp(-0.5*(x-[0])*(x-[0])/([1]*[1]))",-50,50);
  vtxGaussian_DATA->SetParameters(pvtxG_DATA_[0],pvtxG_DATA_[1],pvtxG_DATA_[2],pvtxG_DATA_[3]);
  vtxGaussian_MC = new TF1("vtxGaussian_MC","[2]*TMath::Exp(-0.5*(x-[0])*(x-[0])/([1]*[1]))",-50,50);
  vtxGaussian_MC->SetParameters(pvtxG_MC_[0],pvtxG_MC_[1],pvtxG_MC_[2],pvtxG_MC_[3]);


  hL1TechBits = f->make<TH1D>("hL1TechBits","L1 technical trigger bits before mask",64,-0.5,63.5);
  hL1AlgoBits = f->make<TH1D>("hL1AlgoBits","L1 algorithm trigger bits before mask",128,-0.5,127.5);
  hHPFracNtrk = f->make<TH2D>("hHPFracNtrk","High purity fraction vs. # of tracks; number of tracks; highPurity fraction",50,0,500,50,0,1);
  hHfTowers   = f->make<TH2D>("hHfTowers","Number of HF towers above threshold; positive side; negative side",80,-0.5,79.5,80,-0.5,79.5);
  hHLTPaths = f->make<TH1D>("hHLTPaths","HLT Paths",5,0,5);
  hHLTPaths->SetBit(TH1::kCanRebin);

  hRecMult = f->make<TH1D>("hRecMult","Charged mult. |#eta|<2.5)",numBins,-0.5,xmax);
  hRecMult_STD = f->make<TH1D>("hRecMult_STD","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
  hRecMult_SPEC = f->make<TH1D>("hRecMult_SPEC","Charged mult. |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
  hRecMult_AGR = f->make<TH1D>("hRecMult_AGR","Charged mult. |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);
  
  if(vtxWeight_) hRecMult_STD_W = f->make<TH1D>("hRecMult_STD_W","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);

  if(includeSelTrk_){
     hRecMult_sel = f->make<TH1D>("hRecMult_sel","Charged mult. |#eta|<2.5)",numBins,-0.5,xmax);
     hRecMult_STD_sel = f->make<TH1D>("hRecMult_STD_sel","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
     hRecMult_SPEC_sel = f->make<TH1D>("hRecMult_SPEC_sel","Charged mult. |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
     hRecMult_AGR_sel = f->make<TH1D>("hRecMult_AGR_sel","Charged mult. |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);
  }

  if(isGEN_) {
    hGenMultInel = f->make<TH1D>("hGenMultInel","Charged mult. (inel.) |#eta|<2.5)",numBins,-0.5,xmax);
    hGenMultNSD = f->make<TH1D>("hGenMultNSD","Charged mult. (NSD) |#eta|<2.5)",numBins,-0.5,xmax);
    hGenMultSD = f->make<TH1D>("hGenMultSD","Charged mult. (SD) |#eta|<2.5)",numBins,-0.5,xmax);
    hGenMultDD = f->make<TH1D>("hGenMultDD","Charged mult. (DD) |#eta|<2.5)",numBins,-0.5,xmax);
    hGenMultND = f->make<TH1D>("hGenMultND","Charged mult. (ND) |#eta|<2.5)",numBins,-0.5,xmax);

    hGenMultInel_STD = f->make<TH1D>("hGenMultInel_STD","Charged mult. (inel.) |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
    hGenMultNSD_STD = f->make<TH1D>("hGenMultNSD_STD","Charged mult. (NSD) |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
    hGenMultSD_STD = f->make<TH1D>("hGenMultSD_STD","Charged mult. (SD) |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
    hGenMultDD_STD = f->make<TH1D>("hGenMultDD_STD","Charged mult. (DD) |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
    hGenMultND_STD = f->make<TH1D>("hGenMultND_STD","Charged mult. (ND) |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);

    hGenMultInel_SPEC = f->make<TH1D>("hGenMultInel_SPEC","Charged mult. (inel.) |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
    hGenMultNSD_SPEC = f->make<TH1D>("hGenMultNSD_SPEC","Charged mult. (NSD) |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
    hGenMultSD_SPEC = f->make<TH1D>("hGenMultSD_SPEC","Charged mult. (SD) |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
    hGenMultDD_SPEC = f->make<TH1D>("hGenMultDD_SPEC","Charged mult. (DD) |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
    hGenMultND_SPEC = f->make<TH1D>("hGenMultND_SPEC","Charged mult. (ND) |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);

    hGenMultInel_AGR = f->make<TH1D>("hGenMultInel_AGR","Charged mult. (inel.) |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);
    hGenMultNSD_AGR = f->make<TH1D>("hGenMultNSD_AGR","Charged mult. (NSD) |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);
    hGenMultSD_AGR = f->make<TH1D>("hGenMultSD_AGR","Charged mult. (SD) |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);
    hGenMultDD_AGR = f->make<TH1D>("hGenMultDD_AGR","Charged mult. (DD) |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);
    hGenMultND_AGR = f->make<TH1D>("hGenMultND_AGR","Charged mult. (ND) |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);

    hGenRecMultNSD = f->make<TH1D>("hGenRecMultNSD","Charged mult. |#eta|<2.5)",numBins,-0.5,xmax);
    hGenRecMultNSD_STD = f->make<TH1D>("hGenRecMultNSD_STD","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
    hGenRecMultNSD_SPEC = f->make<TH1D>("hGenRecMultNSD_SPEC","Charged mult. |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
    hGenRecMultNSD_AGR = f->make<TH1D>("hGenRecMultNSD_AGR","Charged mult. |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);

    hGenRecMultInel = f->make<TH1D>("hGenRecMultInel","Charged mult. |#eta|<2.5)",numBins,-0.5,xmax);
    hGenRecMultInel_STD = f->make<TH1D>("hGenRecMultInel_STD","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
    hGenRecMultInel_SPEC = f->make<TH1D>("hGenRecMultInel_SPEC","Charged mult. |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
    hGenRecMultInel_AGR = f->make<TH1D>("hGenRecMultInel_AGR","Charged mult. |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);

    hGenRecMultNSD_STD_W = f->make<TH1D>("hGenRecMultNSD_STD_W","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);

    hGenRecMultSD = f->make<TH1D>("hGenRecMultSD","Charged mult. |#eta|<2.5)",numBins,-0.5,xmax);
    hGenRecMultSD_STD = f->make<TH1D>("hGenRecMultSD_STD","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);
    hGenRecMultSD_SPEC = f->make<TH1D>("hGenRecMultSD_SPEC","Charged mult. |#eta|<1.0 with min p_{T})",numBins,-0.5,xmax_SPEC);
    hGenRecMultSD_AGR = f->make<TH1D>("hGenRecMultSD_AGR","Charged mult. |#eta|<0.8 with min p_{T})",numBins,-0.5,xmax_AGR);

    if(isGEN_ && vtxWeight_) hGenRecMultSD_STD_W = f->make<TH1D>("hGenRecMultSD_STD_W","Charged mult. |#eta|<2.4 with min p_{T})",numBins,-0.5,xmax_STD);

    hGenToRecZeroBinNSD_STD = f->make<TH1D>("hGenToRecZeroBinNSD_STD","GEN mult. for REC mult. of zero",numBins,-0.5,xmax_STD);
    hGenToRec0123BinNSD_STD = f->make<TH1D>("hGenToRec0123BinNSD_STD","GEN mult. for REC mult. of zero",numBins,-0.5,xmax_STD);
    hGenVsRecMultNSD_STD = f->make<TH2F>("hGenVsRecMultNSD_STD","GEN vs REC;GEN mult.;REC mult.",numBins,-0.5,xmax_STD,
				      numBins,-0.5,xmax_STD);

    if(isGEN_ && vtxWeight_) hGenRecMultNSDvsZvtx_STD = f->make<TH2F>("hGenRecMultNSDvsZvtx_STD","Mult vs Zvtx;mult.;V_{z}",numBins,-0.5,xmax_STD,
							      80,-20,20);
    
  }

}

// ------------ method called once each job after completing the event loop  ------------
void 
EvtSelAnalyzer::endJob()
{

  hHLTPaths->LabelsDeflate();

}

//define this as a plug-in
DEFINE_FWK_MODULE(EvtSelAnalyzer);
