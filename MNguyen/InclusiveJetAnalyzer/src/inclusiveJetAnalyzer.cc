/*
  Based on the jet response analyzer
  Modified by Matt Nguyen, November 2010

*/

#include "MNguyen/InclusiveJetAnalyzer/interface/inclusiveJetAnalyzer.h"


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <Math/DistFunc.h>
#include "TMath.h"



#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"




using namespace std;
using namespace edm;
using namespace reco;


InclusiveJetAnalyzer::InclusiveJetAnalyzer(const edm::ParameterSet& iConfig) {
  

  jetTag_ = iConfig.getParameter<InputTag>("jetTag");
	vtxTag_ = iConfig.getUntrackedParameter<edm::InputTag>("vtxTag",edm::InputTag("hiSelectedVertex"));

  isMC_ = iConfig.getUntrackedParameter<bool>("isMC",false);

  if(isMC_){
    genjetTag_ = iConfig.getParameter<InputTag>("genjetTag");
    eventInfoTag_ = iConfig.getParameter<InputTag>("eventInfoTag");
  }
  verbose_ = iConfig.getUntrackedParameter<bool>("verbose",false);

  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality",false);
  useVtx_ = iConfig.getUntrackedParameter<bool>("useVtx",true);
  useJEC_ = iConfig.getUntrackedParameter<bool>("useJEC",true);

  if(!isMC_){
    L1gtReadout_ = iConfig.getParameter<edm::InputTag>("L1gtReadout");
    hltResName_ = iConfig.getUntrackedParameter<string>("hltTrgResults","TriggerResults::HLT");
    
    
    if (iConfig.exists("hltTrgNames"))
      hltTrgNames_ = iConfig.getUntrackedParameter<vector<string> >("hltTrgNames");
    
    if (iConfig.exists("hltProcNames"))
      hltProcNames_ = iConfig.getUntrackedParameter<vector<string> >("hltProcNames");
    else {
      hltProcNames_.push_back("FU");
      hltProcNames_.push_back("HLT");
    }
  }

  cout<<" jet collection : "<<jetTag_<<endl;
  if(isMC_)cout<<" genjet collection : "<<genjetTag_<<endl;


   
}



InclusiveJetAnalyzer::~InclusiveJetAnalyzer() { }



void 
InclusiveJetAnalyzer::beginRun(const edm::Run& run, 
			      const edm::EventSetup & es) {}

void 
InclusiveJetAnalyzer::beginJob() {

  centrality_ = 0;

  //string jetTagName = jetTag_.label()+"_tree"; 
  string jetTagTitle = jetTag_.label()+" Jet Analysis Tree"; 
  t = fs1->make<TTree>("t",jetTagTitle.c_str());


  //  TTree* t= new TTree("t","Jet Response Analyzer");
  t->Branch("run",&jets_.run,"run/I");
  t->Branch("evt",&jets_.evt,"evt/I");
  t->Branch("lumi",&jets_.lumi,"lumi/I");
  t->Branch("b",&jets_.b,"b/F");
  t->Branch("vx",&jets_.vx,"vx/F");
  t->Branch("vy",&jets_.vy,"vy/F");
  t->Branch("vz",&jets_.vz,"vz/F");
  t->Branch("hf",&jets_.hf,"hf/F");
  t->Branch("nref",&jets_.nref,"nref/I");
  t->Branch("bin",&jets_.bin,"bin/I");
  t->Branch("rawpt",jets_.rawpt,"rawpt[nref]/F");
  t->Branch("jtpt",jets_.jtpt,"jtpt[nref]/F");
  t->Branch("jteta",jets_.jteta,"jteta[nref]/F");
  t->Branch("jty",jets_.jty,"jty[nref]/F");
  t->Branch("jtphi",jets_.jtphi,"jtphi[nref]/F");

  if(isMC_){
    t->Branch("pthat",&jets_.pthat,"pthat/F");    

    // Only matched gen jets
    t->Branch("refpt",jets_.refpt,"refpt[nref]/F");
    t->Branch("refeta",jets_.refeta,"refeta[nref]/F");
    t->Branch("refy",jets_.refy,"refy[nref]/F");
    t->Branch("refphi",jets_.refphi,"refphi[nref]/F");
    t->Branch("refdphijt",jets_.refdphijt,"refdphijt[nref]/F");
    t->Branch("refdrjt",jets_.refdrjt,"refdrjt[nref]/F");
    // matched parton
    t->Branch("refparton_pt",jets_.refparton_pt,"refparton_pt[nref]/F");
    t->Branch("refparton_flavor",jets_.refparton_flavor,"refparton_flavor[nref]/F");

    // For all gen jets, matched or unmatched
    t->Branch("ngen",&jets_.ngen,"ngen/I");
    t->Branch("genmatchindex",jets_.genmatchindex,"genmatchindex[ngen]/I");
    t->Branch("genpt",jets_.genpt,"genpt[ngen]/F");
    t->Branch("geneta",jets_.geneta,"geneta[ngen]/F");
    t->Branch("geny",jets_.geny,"geny[ngen]/F");
    t->Branch("genphi",jets_.genphi,"genphi[ngen]/F");
    t->Branch("gendphijt",jets_.gendphijt,"gendphijt[ngen]/F");
    t->Branch("gendrjt",jets_.gendrjt,"gendrjt[ngen]/F");
  }
  
  if(!isMC_){
    t->Branch("nL1TBit",&jets_.nL1TBit,"nL1TBit/I");
    t->Branch("l1TBit",jets_.l1TBit,"l1TBit[nL1TBit]/O");

    t->Branch("nL1ABit",&jets_.nL1ABit,"nL1ABit/I");
    t->Branch("l1ABit",jets_.l1ABit,"l1ABit[nL1ABit]/O");

    t->Branch("nHLTBit",&jets_.nHLTBit,"nHLTBit/I");
    t->Branch("hltBit",jets_.hltBit,"hltBit[nHLTBit]/O");

  }

  TH1D::SetDefaultSumw2();
  
  
}


void 
InclusiveJetAnalyzer::analyze(const Event& iEvent, 
			     const EventSetup& iSetup) {
  
  int event = iEvent.id().event();
  int run = iEvent.id().run();
  int lumi = iEvent.id().luminosityBlock();
  
  jets_.run = run;
  jets_.evt = event;
  jets_.lumi = lumi;

  LogDebug("InclusiveJetAnalyzer")<<"START event: "<<event<<" in run "<<run<<endl;


 int bin = -1;
  double hf = 0.;
  double b = 999.;

  if(useCentrality_){
    //if(!isMC_){
      if(!centrality_) centrality_ = new CentralityProvider(iSetup);      
      centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
      //double c = centrality_->centralityValue();
      const reco::Centrality *cent = centrality_->raw();
      
      hf = cent->EtHFhitSum();

      bin = centrality_->getBin();
      b = centrality_->bMean();
      //}
      /*
    else{

      TFile * centFile = new TFile("/net/hidsk0001/d00/scratch/mnguyen/CMSSW_3_9_1_patch1/src/macros/Hydjet_CentTable.root");

      edm::Handle<reco::Centrality> cent;
      iEvent.getByLabel(edm::InputTag("hiCentrality"),cent);
      //cout<<" grabbed centrality "<<endl;
      CentralityBins::RunMap cmap = getCentralityFromFile(centFile, "makeCentralityTableTFile", "HFhitsHydjet_2760GeV", 1, 1);

      // Still getting cent from hits.  When tower tables become available, we need to switch
      hf = cent->EtHFhitSum();
      //cout<<" hf "<<hf<<endl;
      bin = cmap[run]->getBin(hf);
      b = cmap[run]->bMeanOfBin(bin);
    }
      */
  }
   



  //double jetPtMin = 35.;


   // loop the events
   
   jets_.bin = bin;
   jets_.hf = hf;
   

	 if (useVtx_) {
      edm::Handle<vector<reco::Vertex> >vertex;
      iEvent.getByLabel(vtxTag_, vertex);

      if(vertex->size()>0) {
        jets_.vx=vertex->begin()->x();
        jets_.vy=vertex->begin()->y();
        jets_.vz=vertex->begin()->z();
      }
	 }
   


   edm::Handle<pat::JetCollection> jets;
   iEvent.getByLabel(jetTag_, jets);

   
   // FILL JRA TREE

   jets_.b = b;
   jets_.nref = 0;
   
   if(!isMC_){
     fillL1Bits(iEvent);
     fillHLTBits(iEvent);
   }

   for(unsigned int j = 0 ; j < jets->size(); ++j){
     const pat::Jet& jet = (*jets)[j];
     
     //cout<<" jet pt "<<jet.pt()<<endl;
     //if(jet.pt() < jetPtMin) continue;
     if (useJEC_) jets_.rawpt[jets_.nref]=jet.correctedJet("Uncorrected").pt();
     jets_.jtpt[jets_.nref] = jet.pt();                            
     jets_.jteta[jets_.nref] = jet.eta();
     jets_.jtphi[jets_.nref] = jet.phi();
     jets_.jty[jets_.nref] = jet.eta();
     
	 
     if(jet.genJet()){
       jets_.refpt[jets_.nref] = jet.genJet()->pt();                            
       jets_.refeta[jets_.nref] = jet.genJet()->eta();
       jets_.refphi[jets_.nref] = jet.genJet()->phi();
       jets_.refy[jets_.nref] = jet.genJet()->eta();
       jets_.refdphijt[jets_.nref] = reco::deltaPhi(jet.phi(),jet.genJet()->phi());	
       jets_.refdrjt[jets_.nref] = reco::deltaR(jet.eta(),jet.phi(),jet.genJet()->eta(),jet.genJet()->phi());	       
     }            	
     else{
       jets_.refpt[jets_.nref] = -999.;
       jets_.refeta[jets_.nref] = -999.;
       jets_.refphi[jets_.nref] = -999.;
       jets_.refy[jets_.nref] = -999.;
       jets_.refdphijt[jets_.nref] = -999.;
       jets_.refdrjt[jets_.nref] = -999.;
     }

     // matched partons
     if (jet.genParton()) {
       jets_.refparton_pt[jets_.nref] = jet.genParton()->pt();
       jets_.refparton_flavor[jets_.nref] = jet.genParton()->pdgId();
     } else {
       jets_.refparton_pt[jets_.nref] = -999;
       jets_.refparton_flavor[jets_.nref] = -999;
     }

 
     jets_.nref++;
       
       
   }


   if(isMC_){

     edm::Handle<GenEventInfoProduct> hEventInfo;
     iEvent.getByLabel(eventInfoTag_,hEventInfo);
     //jets_.pthat = hEventInfo->binningValues()[0];

     // binning values and qscale appear to be equivalent, but binning values not always present
     jets_.pthat = hEventInfo->qScale();

     edm::Handle<vector<reco::GenJet> >genjets;
     iEvent.getByLabel(genjetTag_, genjets);

     jets_.ngen = 0;
     for(unsigned int igen = 0 ; igen < genjets->size(); ++igen){
       const reco::GenJet & genjet = (*genjets)[igen];
       
       float genjet_pt = genjet.pt();
       
       // threshold to reduce size of output in minbias PbPb
       if(genjet_pt>20.){

	 jets_.genpt[jets_.ngen] = genjet_pt;                            
	 jets_.geneta[jets_.ngen] = genjet.eta();
	 jets_.genphi[jets_.ngen] = genjet.phi();
	 jets_.geny[jets_.ngen] = genjet.eta();
	 
	 
	 // find matching patJet if there is one
	 
	 jets_.gendrjt[jets_.ngen] = -1.0;
	 jets_.genmatchindex[jets_.ngen] = -1;
	 
	 
	 for(unsigned int ijet = 0 ; ijet < jets->size(); ++ijet){
	   const pat::Jet& jet = (*jets)[ijet];
	   
	   if(jet.genJet()){
	     if(fabs(genjet.pt()-jet.genJet()->pt())<0.0001 &&
		fabs(genjet.eta()-jet.genJet()->eta())<0.0001 &&
		(fabs(genjet.phi()-jet.genJet()->phi())<0.0001 || fabs(genjet.phi()-jet.genJet()->phi() - 2.0*TMath::Pi()) < 0.0001 )){
	       
	       jets_.genmatchindex[jets_.ngen] = (int)ijet;
	       jets_.gendphijt[jets_.ngen] = reco::deltaPhi(jet.phi(),genjet.phi());	
	       jets_.gendrjt[jets_.ngen] = reco::deltaR(jet,genjet);	
	       
	       break;
	     }            		
	   }
	 }
	 jets_.ngen++;
       }
       
     }
     
   }
   t->Fill();
}


  

//--------------------------------------------------------------------------------------------------
void InclusiveJetAnalyzer::fillL1Bits(const edm::Event &iEvent)
{
  edm::Handle< L1GlobalTriggerReadoutRecord >  L1GlobalTrigger;

  iEvent.getByLabel(L1gtReadout_, L1GlobalTrigger); 
  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = L1GlobalTrigger->technicalTriggerWord();

  for (int i=0; i<64;i++)
    {
      jets_.l1TBit[i] = technicalTriggerWordBeforeMask.at(i);
    }
  jets_.nL1TBit = 64;

  int ntrigs = L1GlobalTrigger->decisionWord().size();
  jets_.nL1ABit = ntrigs;

  for (int i=0; i != ntrigs; i++) {
    bool accept = L1GlobalTrigger->decisionWord()[i];
    //jets_.l1ABit[i] = (accept == true)? 1:0;
    if(accept== true){
      jets_.l1ABit[i] = 1;
    }
    else{
      jets_.l1ABit[i] = 0;
    }
    
  }
}

//--------------------------------------------------------------------------------------------------
void InclusiveJetAnalyzer::fillHLTBits(const edm::Event &iEvent)
{
  // Fill HLT trigger bits.
  Handle<TriggerResults> triggerResultsHLT;
  getProduct(hltResName_, triggerResultsHLT, iEvent);

  const TriggerResults *hltResults = triggerResultsHLT.product();
  const TriggerNames & triggerNames = iEvent.triggerNames(*hltResults);

  jets_.nHLTBit = hltTrgNames_.size();

  for(size_t i=0;i<hltTrgNames_.size();i++){
   
    for(size_t j=0;j<triggerNames.size();++j) {
      
      if(triggerNames.triggerName(j) == hltTrgNames_[i]){
	
	//cout <<"hltTrgNames_(i) "<<hltTrgNames_[i]<<endl;
	//cout <<"triggerName(j) "<<triggerNames.triggerName(j)<<endl;
	//cout<<" result "<<triggerResultsHLT->accept(j)<<endl;
	jets_.hltBit[i] = triggerResultsHLT->accept(j);
      }
      
    }
  }
}

//--------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void InclusiveJetAnalyzer::getProduct(const std::string name, edm::Handle<TYPE> &prod,
					 const edm::Event &event) const
{
  // Try to access data collection from EDM file. We check if we really get just one
  // product with the given name. If not we throw an exception.

  event.getByLabel(edm::InputTag(name),prod);
  if (!prod.isValid()) 
    throw edm::Exception(edm::errors::Configuration, "InclusiveJetAnalyzer::GetProduct()\n")
      << "Collection with label '" << name << "' is not valid" <<  std::endl;
}

//--------------------------------------------------------------------------------------------------
template <typename TYPE>
inline bool InclusiveJetAnalyzer::getProductSafe(const std::string name, edm::Handle<TYPE> &prod,
					     const edm::Event &event) const
{
  // Try to safely access data collection from EDM file. We check if we really get just one
  // product with the given name. If not, we return false.

  if (name.size()==0)
    return false;

  try {
    event.getByLabel(edm::InputTag(name),prod);
    if (!prod.isValid()) 
      return false;
  } catch (...) {
    return false;
  }
  return true;
}



				     					
DEFINE_FWK_MODULE(InclusiveJetAnalyzer);
