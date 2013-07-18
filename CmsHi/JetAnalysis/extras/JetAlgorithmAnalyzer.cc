#ifndef RecoJets_JetProducers_plugins_JetAlgorithmAnalyzer_h
#define RecoJets_JetProducers_plugins_JetAlgorithmAnalyzer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/CodedException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoJets/JetProducers/interface/PileUpSubtractor.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ActiveAreaSpec.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/CMSIterativeConePlugin.hh"
#include "fastjet/ATLASConePlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <limits>
#include <cmath>
#include <boost/shared_ptr.hpp>


#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Math/interface/Vector3D.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "TNtuple.h"
#include "TH2D.h"
#include <vector>

static const int nSteps = 8;
static const double PI = 3.141592653589;

class JetAlgorithmAnalyzer : public edm::EDProducer{

public:
  //
  // construction/destruction
  //
  explicit JetAlgorithmAnalyzer(const edm::ParameterSet& iConfig);
  virtual ~JetAlgorithmAnalyzer();

protected:

  //
  // member functions
  //

   virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );
  virtual void runAlgorithm( edm::Event& iEvent, const edm::EventSetup& iSetup );
  virtual void output(  edm::Event & iEvent, edm::EventSetup const& iSetup );
  template< typename T >
  void writeBkgJets( edm::Event & iEvent, edm::EventSetup const& iSetup );

   void fillNtuple(int output, const  std::vector<fastjet::PseudoJet>& jets, int step);
   void fillTowerNtuple(const  std::vector<fastjet::PseudoJet>& jets, int step);
   void fillJetNtuple(const  std::vector<fastjet::PseudoJet>& jets, int step);
   void fillBkgNtuple(const PileUpSubtractor* subtractor, int step);

   // VirtualJetProducer Stuff 
   //
   // typedefs & structs
   //
   struct JetType {
      enum Type {
	 BasicJet,
	 GenJet,
	 CaloJet,
	 PFJet,
	 TrackJet,
	 LastJetType  // no real type, technical
      };
      static const char *names[];
      static Type byName(const std::string &name);
   };

   JetType::Type jetTypeE;

   inline bool makeCaloJet(const JetType::Type &fTag) {
      return fTag == JetType::CaloJet;
   }
   inline bool makePFJet(const JetType::Type &fTag) {
      return fTag == JetType::PFJet;
   }
   inline bool makeGenJet(const JetType::Type &fTag) {
      return fTag == JetType::GenJet;
   }
   inline bool makeTrackJet(const JetType::Type &fTag) {
      return fTag == JetType::TrackJet;
   }
   inline bool makeBasicJet(const JetType::Type &fTag) {
      return fTag == JetType::BasicJet;
   }

public:
   // typedefs
   typedef boost::shared_ptr<fastjet::ClusterSequence>        ClusterSequencePtr;
   typedef boost::shared_ptr<fastjet::JetDefinition::Plugin>  PluginPtr;
   typedef boost::shared_ptr<fastjet::JetDefinition>          JetDefPtr;
   typedef boost::shared_ptr<fastjet::ActiveAreaSpec>         ActiveAreaSpecPtr;
   typedef boost::shared_ptr<fastjet::RangeDefinition>        RangeDefPtr;
   std::string   jetType() const { return jetType_; }

 private:

  // trackjet clustering parameters
  bool useOnlyVertexTracks_;
  bool useOnlyOnePV_;
  float dzTrVtxMax_;

  double phi0_;
   int nFill_;
   float etaMax_;
   int iev_;
   bool avoidNegative_;
   bool doAnalysis_;
   bool doMC_;
   
   int centBin_;
   const CentralityBins* cbins_;
      TNtuple* ntTowers;
   TNtuple* ntJetTowers;
   TNtuple* ntTowersFromEvtPlane;
   TNtuple* ntJets;
   TNtuple* ntPU;
   TNtuple* ntRandom;
   TNtuple* ntuple;
  std::vector<TH2D*> hTowers;
  std::vector<TH2D*> hJetTowers;
  std::vector<TH2D*> hTowersFromEvtPlane;
  std::vector<TH2D*> hJets;
  std::vector<TH1D*> hPU;
  std::vector<TH1D*> hMean;
  std::vector<TH1D*> hSigma;
   std::vector<TH2D*> hRandom;
   const CaloGeometry *geo;
   edm::Service<TFileService> f;

   // VirtualJetProducer Stuff
 
protected:
   virtual void inputTowers();
   virtual bool isAnomalousTower(reco::CandidatePtr input);
   virtual void copyConstituents(const std::vector<fastjet::PseudoJet>&fjConstituents,
				 reco::Jet* jet);

   // This will run the actual algorithm. This method is pure virtual and
   // has no default.
  template< typename T >
  void writeJets( edm::Event & iEvent, edm::EventSetup const& iSetup );

   // This method copies the constituents from the fjConstituents method
   // to an output of CandidatePtr's.
   virtual std::vector<reco::CandidatePtr>
   getConstituents(const std::vector<fastjet::PseudoJet>&fjConstituents);


protected:
   std::string           moduleLabel_;               // label for this module
   edm::InputTag         src_;                       // input constituent source
   edm::InputTag         srcPVs_;                    // primary vertex source
   std::string           jetType_;                   // type of jet (Calo,PF,Basic,Gen)
   std::string           jetAlgorithm_;              // the jet algorithm to use
   double                rParam_;                    // the R parameter to use
   double                inputEtMin_;                // minimum et of input constituents
   double                inputEMin_;                 // minimum e of input constituents
   double                jetPtMin_;                  // minimum jet pt
   bool                  doPVCorrection_;            // correct to primary vertex?

   // for restricting inputs due to processing time
   bool                  restrictInputs_;            // restrict inputs to first "maxInputs" inputs.
   unsigned int          maxInputs_;                 // maximum number of inputs.

   // for fastjet jet area calculation
   bool                  doAreaFastjet_;             // calculate area w/ fastjet?
   // for fastjet rho calculation
   bool                  doRhoFastjet_;              // calculate rho w/ fastjet?

   // for pileup offset correction
   bool                  doPUOffsetCorr_;            // add the pileup calculation from offset correction?
   std::string           puSubtractorName_;
   // anomalous cell cuts
   unsigned int          maxBadEcalCells_;           // maximum number of bad ECAL cells
   unsigned int          maxRecoveredEcalCells_;     // maximum number of recovered ECAL cells
   unsigned int          maxProblematicEcalCells_;   // maximum number of problematic ECAL cells
   unsigned int          maxBadHcalCells_;           // maximum number of bad HCAL cells
   unsigned int          maxRecoveredHcalCells_;     // maximum number of recovered HCAL cells
   unsigned int          maxProblematicHcalCells_;   // maximum number of problematic HCAL cells
   
   std::vector<edm::Ptr<reco::Candidate> > inputs_;  // input candidates [View, PtrVector and CandColle
   reco::Particle::Point           vertex_;          // Primary vertex
   ClusterSequencePtr              fjClusterSeq_;    // fastjet cluster sequence
   JetDefPtr                       fjJetDefinition_; // fastjet jet definition
   PluginPtr                       fjPlugin_;        // fastjet plugin
   ActiveAreaSpecPtr               fjActiveArea_;    // fastjet active area definition
   RangeDefPtr                     fjRangeDef_;      // range definition
   std::vector<fastjet::PseudoJet> fjInputs_;        // fastjet inputs
   std::vector<fastjet::PseudoJet> fjJets_;          // fastjet jets
   
   std::string                     jetCollInstanceName_;       // instance name for output jet collection
   boost::shared_ptr<PileUpSubtractor>  subtractor_;
   
};


#endif
////////////////////////////////////////////////////////////////////////////////
//
// JetAlgorithmAnalyzer
// ------------------
//
//            04/21/2009 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////

//#include "RecoJets/JetProducers/plugins/JetAlgorithmAnalyzer.h"
//#include "JetAlgorithmAnalyzer.h"

#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/CodedException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "fastjet/SISConePlugin.hh"
#include "fastjet/CMSIterativeConePlugin.hh"
#include "fastjet/ATLASConePlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"

#include "CLHEP/Random/RandomEngine.h"

#include <iostream>
#include <memory>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace std;
using namespace edm;

static const double pi = 3.14159265358979;

CLHEP::HepRandomEngine* randomEngine;


////////////////////////////////////////////////////////////////////////////////
// construction / destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
JetAlgorithmAnalyzer::JetAlgorithmAnalyzer(const edm::ParameterSet& iConfig)
   : phi0_(0),
     nFill_(5),
     etaMax_(3),
     iev_(0),
     geo(0),
     cbins_(0),
     moduleLabel_   (iConfig.getParameter<string>       ("@module_label"))
   , src_           (iConfig.getParameter<edm::InputTag>("src"))
   , srcPVs_        (iConfig.getParameter<edm::InputTag>("srcPVs"))
   , jetType_       (iConfig.getParameter<string>       ("jetType"))
   , jetAlgorithm_  (iConfig.getParameter<string>       ("jetAlgorithm"))
   , rParam_        (iConfig.getParameter<double>       ("rParam"))
   , inputEtMin_    (iConfig.getParameter<double>       ("inputEtMin"))
   , inputEMin_     (iConfig.getParameter<double>       ("inputEMin"))
   , jetPtMin_      (iConfig.getParameter<double>       ("jetPtMin"))
   , doPVCorrection_(iConfig.getParameter<bool>         ("doPVCorrection"))
   , restrictInputs_(false)
   , maxInputs_(99999999)
   , doAreaFastjet_ (iConfig.getParameter<bool>         ("doAreaFastjet"))
   , doRhoFastjet_  (iConfig.getParameter<bool>         ("doRhoFastjet"))
   , doPUOffsetCorr_(iConfig.getParameter<bool>         ("doPUOffsetCorr"))
   , maxBadEcalCells_        (iConfig.getParameter<unsigned>("maxBadEcalCells"))
   , maxRecoveredEcalCells_  (iConfig.getParameter<unsigned>("maxRecoveredEcalCells"))
   , maxProblematicEcalCells_(iConfig.getParameter<unsigned>("maxProblematicEcalCells"))
   , maxBadHcalCells_        (iConfig.getParameter<unsigned>("maxBadHcalCells"))
   , maxRecoveredHcalCells_  (iConfig.getParameter<unsigned>("maxRecoveredHcalCells"))
   , maxProblematicHcalCells_(iConfig.getParameter<unsigned>("maxProblematicHcalCells"))
   , jetCollInstanceName_ ("")
{
   edm::Service<RandomNumberGenerator> rng;
   randomEngine = &(rng->getEngine());

   doAnalysis_  = iConfig.getUntrackedParameter<bool>("doAnalysis",true);
   doMC_  = iConfig.getUntrackedParameter<bool>("doMC",true);
   
   if(doAnalysis_) centBin_ = iConfig.getUntrackedParameter<int>("centrality",0);
   
   avoidNegative_  = iConfig.getParameter<bool>("avoidNegative");
  
  if ( iConfig.exists("UseOnlyVertexTracks") )
    useOnlyVertexTracks_ = iConfig.getParameter<bool>("UseOnlyVertexTracks");
  else 
    useOnlyVertexTracks_ = false;
  
  if ( iConfig.exists("UseOnlyOnePV") )
    useOnlyOnePV_        = iConfig.getParameter<bool>("UseOnlyOnePV");
  else
    useOnlyOnePV_ = false;

  if ( iConfig.exists("DrTrVtxMax") )
    dzTrVtxMax_          = iConfig.getParameter<double>("DzTrVtxMax");
  else
    dzTrVtxMax_ = false;

  if(doAnalysis_){
     ntTowers = f->make<TNtuple>("ntTowers","Algorithm Analysis Towers","eta:phi:et:step:event");
     ntTowersFromEvtPlane = f->make<TNtuple>("ntTowersFromEvtPlane","Algorithm Analysis Towers","eta:phi:et:step:event");
     ntJetTowers = f->make<TNtuple>("ntTowers","Algorithm Analysis Towers","eta:phi:et:step:event");

     ntJets = f->make<TNtuple>("ntJets","Algorithm Analysis Jets","eta:phi:et:step:event");
     ntPU = f->make<TNtuple>("ntPU","Algorithm Analysis Background","eta:mean:sigma:step:event");
     ntRandom = f->make<TNtuple>("ntRandom","Algorithm Analysis Background","eta:phi:et:pu:event");

     ntuple = f->make<TNtuple>("nt","debug","ieta:eta:iphi:phi:pt:em:had");

     for(int i = 0; i < nSteps; ++i){
       hTowers.push_back(f->make<TH2D>(Form("hTowers_step%d",i),"",200,-5.5,5.5,200,-0.02,6.3));
       hJets.push_back(f->make<TH2D>(Form("hJets_step%d",i),"",200,-5.5,5.5,200,-0.02,6.3));
       hTowersFromEvtPlane.push_back(f->make<TH2D>(Form("hTowersFromEvtPlane_step%d",i),"",200,-5.5,5.5,200,-0.02,6.3));

       hJetTowers.push_back(f->make<TH2D>(Form("hJetTowers_step%d",i),"",200,-5.5,5.5,200,-0.02,6.3));

       hPU.push_back(f->make<TH1D>(Form("hPU_step%d",i),"",200,-5.5,5.5));
       hMean.push_back(f->make<TH1D>(Form("hMean_step%d",i),"",200,-5.5,5.5));
       hSigma.push_back(f->make<TH1D>(Form("hSigma_step%d",i),"",200,-5.5,5.5));

     }

  }


  if (jetAlgorithm_=="IterativeCone") {
     fjPlugin_ = PluginPtr(new fastjet::CMSIterativeConePlugin(rParam_,1.0));
     fjJetDefinition_= JetDefPtr(new fastjet::JetDefinition(&*fjPlugin_));
  }
  else if (jetAlgorithm_=="Kt")
     fjJetDefinition_= JetDefPtr(new fastjet::JetDefinition(fastjet::kt_algorithm,rParam_));
  else if (jetAlgorithm_=="AntiKt")
     fjJetDefinition_= JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm,rParam_) );
  else
     throw cms::Exception("Invalid jetAlgorithm")
	<<"Jet algorithm for MyVirtualJetProducer is invalid, Abort!\n";

  jetTypeE=JetType::byName(jetType_);

  if ( iConfig.exists("jetCollInstanceName") ) {
     jetCollInstanceName_ = iConfig.getParameter<string>("jetCollInstanceName");
  }

  if ( doPUOffsetCorr_ ) {
     if ( jetTypeE != JetType::CaloJet ) {
	throw cms::Exception("InvalidInput") << "Can only offset correct jets of type CaloJet";
     }
     //     subtractor_ = boost::shared_ptr<PileUpSubtractor>(new PileUpSubtractor(iConfig));
     puSubtractorName_  =  iConfig.getParameter<string> ("subtractorName");
     if(puSubtractorName_.empty()){
	edm::LogWarning("VirtualJetProducer") << "Pile Up correction on; however, pile up type is not specified. Using default... \n";
	subtractor_ =  boost::shared_ptr<PileUpSubtractor>(new PileUpSubtractor(iConfig));
     }else{
	subtractor_ =  boost::shared_ptr<PileUpSubtractor>(PileUpSubtractorFactory::get()->create( puSubtractorName_, iConfig) );
     }
  }

  // do fasjet area / rho calcluation? => accept corresponding parameters
  if ( doAreaFastjet_ || doRhoFastjet_ ) {
     // default Ghost_EtaMax should be 5
     double ghostEtaMax = iConfig.getParameter<double>("Ghost_EtaMax");
     // default Active_Area_Repeats 1
     int    activeAreaRepeats = iConfig.getParameter<int> ("Active_Area_Repeats");
     // default GhostArea 0.01
     double ghostArea = iConfig.getParameter<double> ("GhostArea");
     fjActiveArea_ =  ActiveAreaSpecPtr(new fastjet::ActiveAreaSpec(ghostEtaMax,
								    activeAreaRepeats,
								    ghostArea));
     fjRangeDef_ = RangeDefPtr( new fastjet::RangeDefinition(ghostEtaMax) );
  }
  
  // restrict inputs to first "maxInputs" towers?
  if ( iConfig.exists("restrictInputs") ) {
     restrictInputs_ = iConfig.getParameter<bool>("restrictInputs");
     maxInputs_      = iConfig.getParameter<unsigned int>("maxInputs");
  }


  string alias=iConfig.getUntrackedParameter<string>("alias",moduleLabel_);

  // make the "produces" statements
  //  produces<reco::CaloJetCollection>();
  produces<reco::CaloJetCollection>("randomCones");
  produces<std::vector<bool> >("directions");
  
  if(doRhoFastjet_){
     produces<double>("rho");
     produces<double>("sigma");
  }
  
}


//______________________________________________________________________________
JetAlgorithmAnalyzer::~JetAlgorithmAnalyzer()
{
} 

void JetAlgorithmAnalyzer::fillNtuple(int output, const  std::vector<fastjet::PseudoJet>& jets, int step){
   if(!doAnalysis_) return;

   TNtuple* nt;
   TH2D* h;

   if(output == 1){
     nt = ntJets;
     h = hJets[step];
   }
   if(output == 0){
     nt = ntTowers;
     h = hTowers[step];
   }
   if(output == 2){
     nt = ntTowersFromEvtPlane;
     h = hTowersFromEvtPlane[step];
   }
   if(output == 3){
     nt = ntJetTowers;
     h = hJetTowers[step];
   }

   double totet = 0;
   int ntow = 0;
   int nj = jets.size();

   cout<<"step : "<<step<<endl;
   cout<<"Size of input : "<<nj<<endl;

   for(unsigned int i = 0; i < jets.size(); ++i){
     const fastjet::PseudoJet& jet = jets[i];

     double pt = jet.perp();

     if(output != 1){
       reco::CandidatePtr const & itow =  inputs_[ jet.user_index() ];
       pt =  itow->et();
     }
     
     double phi = jet.phi();
     if(output == 2){
       phi = phi - phi0_;
       if(phi < 0) phi += 2*PI;
       if(phi > 2*PI) phi -= 2*PI;
     }
     
     double eta = jet.eta();
     
     if(eta > 0 && eta < 0.08){
       //     if(fabs(eta) < 1.){ 
       totet += pt;
       ntow++;
     }

     nt->Fill(jet.eta(),phi,pt,step,iev_);
     h->Fill(jet.eta(),phi,pt);
   }
   
   cout<<"-----------------------------"<<endl;
   cout<<"STEP             = "<<step<<endl;   
   cout<<"Total ET         = "<<totet<<endl;
   cout<<"Towers counted   = "<<ntow<<endl;
   cout<<"Average tower ET = "<<totet/ntow<<endl;
   cout<<"-----------------------------"<<endl;     
   
}


void JetAlgorithmAnalyzer::fillTowerNtuple(const  std::vector<fastjet::PseudoJet>& jets, int step){
   fillNtuple(0,jets,step);
   fillNtuple(2,jets,step);
}
 
void JetAlgorithmAnalyzer::fillJetNtuple(const  std::vector<fastjet::PseudoJet>& jets, int step){
   fillNtuple(1,jets,step);

   for(unsigned int i = 0; i < jets.size(); ++i){
     const fastjet::PseudoJet& jet = jets[i];
     std::vector<fastjet::PseudoJet> fjConstituents = sorted_by_pt(fjClusterSeq_->constituents(jet));



     
     fillNtuple(3,fjConstituents,step);
   }
}

void JetAlgorithmAnalyzer::fillBkgNtuple(const PileUpSubtractor* subtractor, int step){
   if(!doAnalysis_) return;
   CaloTowerCollection col;
   for(int ieta = -29; ieta < 30; ++ieta){
     if(ieta == 0) continue;
       CaloTowerDetId id(ieta,1);
       const GlobalPoint& hitpoint = geo->getPosition(id);
       cout<<"iETA "<<ieta<<endl;
       double eta = hitpoint.eta();
       cout<<"eta "<<eta<<endl;
       math::PtEtaPhiMLorentzVector p4(1,eta,1.,1.);
       GlobalPoint pos(1,1,1);
       CaloTower c(id,1.,1.,1.,1,1, p4, pos,pos);      
       col.push_back(c);
       reco::CandidatePtr ptr(&col,col.size() - 1);
       double mean = subtractor->getMeanAtTower(ptr);  
       double sigma = subtractor->getSigmaAtTower(ptr);
       double pu = subtractor->getPileUpAtTower(ptr);
       ntPU->Fill(eta,mean,sigma,step,iev_);
       hPU[step]->Fill(eta,pu);
       hMean[step]->Fill(eta,mean);
       hSigma[step]->Fill(eta,sigma);
   }
}

void JetAlgorithmAnalyzer::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{

  //  if(doAnalysis_)   subtractor_->setDebug(1);
  //  else subtractor_->setDebug(0);

   if(!geo){
      edm::ESHandle<CaloGeometry> pGeo;
      iSetup.get<CaloGeometryRecord>().get(pGeo);
      geo = pGeo.product();
   }

   if(doAnalysis_ && centBin_ >= 0){
     edm::Handle<reco::Centrality> cent;
     iEvent.getByLabel(edm::InputTag("hiCentrality"),cent);
     if(!cbins_) cbins_ = getCentralityBinsFromDB(iSetup);
     int bin = cbins_->getBin(cent->EtHFhitSum());
     if(bin != centBin_) return;
   }

   LogDebug("VirtualJetProducer") << "Entered produce\n";
   //determine signal vertex
   vertex_=reco::Jet::Point(0,0,0);
   if (makeCaloJet(jetTypeE)&&doPVCorrection_) {
      LogDebug("VirtualJetProducer") << "Adding PV info\n";
      edm::Handle<reco::VertexCollection> pvCollection;
      iEvent.getByLabel(srcPVs_,pvCollection);
      if (pvCollection->size()>0) vertex_=pvCollection->begin()->position();
   }

   // For Pileup subtraction using offset correction:
   // set up geometry map
   if ( doPUOffsetCorr_ ) {
      subtractor_->setupGeometryMap(iEvent, iSetup);
   }

   // clear data
   LogDebug("VirtualJetProducer") << "Clear data\n";
   fjInputs_.clear();
   fjJets_.clear();
   inputs_.clear();  

   if(doAnalysis_ && doMC_){
     Handle<GenHIEvent> hi;
     iEvent.getByLabel("heavyIon",hi);
     phi0_ = hi->evtPlane();

     double b = hi->b();
     cout<<"===================================="<<endl;
     cout<<"IMPACT PARAMETER OF THE EVENT : "<<b<<endl;
     cout<<"===================================="<<endl;


   }

   // get inputs and convert them to the fastjet format (fastjet::PeudoJet)
   edm::Handle<reco::CandidateView> inputsHandle;
   iEvent.getByLabel(src_,inputsHandle);
   for (size_t i = 0; i < inputsHandle->size(); ++i) {
      inputs_.push_back(inputsHandle->ptrAt(i));
   }
   LogDebug("VirtualJetProducer") << "Got inputs\n";
  
   // Convert candidates to fastjet::PseudoJets.
   // Also correct to Primary Vertex. Will modify fjInputs_
   // and use inputs_
   fjInputs_.reserve(inputs_.size());
   inputTowers();
   LogDebug("VirtualJetProducer") << "Inputted towers\n";

   fillTowerNtuple(fjInputs_,0);
   fillBkgNtuple(subtractor_.get(),0);

   // For Pileup subtraction using offset correction:
   // Subtract pedestal. 

   if ( doPUOffsetCorr_ ) {
      subtractor_->reset(inputs_,fjInputs_,fjJets_);
      subtractor_->calculatePedestal(fjInputs_); 

      fillTowerNtuple(fjInputs_,1);
      fillBkgNtuple(subtractor_.get(),1);
      subtractor_->subtractPedestal(fjInputs_); 

      fillTowerNtuple(fjInputs_,2);
      fillBkgNtuple(subtractor_.get(),2);

      LogDebug("VirtualJetProducer") << "Subtracted pedestal\n";
   }

   // Run algorithm. Will modify fjJets_ and allocate fjClusterSeq_. 
   // This will use fjInputs_
   runAlgorithm( iEvent, iSetup );

   fillTowerNtuple(fjInputs_,3);
   fillBkgNtuple(subtractor_.get(),3);
   fillJetNtuple(fjJets_,3);

   if ( doPUOffsetCorr_ ) {
      subtractor_->setAlgorithm(fjClusterSeq_);
   }

   LogDebug("VirtualJetProducer") << "Ran algorithm\n";

   // For Pileup subtraction using offset correction:
   // Now we find jets and need to recalculate their energy,
   // mark towers participated in jet,
   // remove occupied towers from the list and recalculate mean and sigma
   // put the initial towers collection to the jet,   
   // and subtract from initial towers in jet recalculated mean and sigma of towers 
   if ( doPUOffsetCorr_ ) {
      vector<fastjet::PseudoJet> orphanInput;
      subtractor_->calculateOrphanInput(orphanInput);
      fillTowerNtuple(orphanInput,4);
      fillBkgNtuple(subtractor_.get(),4);
      fillJetNtuple(fjJets_,4);

      //only the id's of the orphan input are used, not their energy
      subtractor_->calculatePedestal(orphanInput);
      fillTowerNtuple(orphanInput,5);
      fillBkgNtuple(subtractor_.get(),5);
      fillJetNtuple(fjJets_,5);

      subtractor_->offsetCorrectJets();
      fillTowerNtuple(orphanInput,6);
      fillBkgNtuple(subtractor_.get(),6);
      fillJetNtuple(fjJets_,6);

   }
  
   // Write the output jets.
   // This will (by default) call the member function template
   // "writeJets", but can be overridden. 
   // this will use inputs_
   output( iEvent, iSetup );
   fillBkgNtuple(subtractor_.get(),7);
   fillJetNtuple(fjJets_,7);

   LogDebug("VirtualJetProducer") << "Wrote jets\n";

   ++iev_;

   doAnalysis_ = false;
   return;
}

void JetAlgorithmAnalyzer::output(edm::Event & iEvent, edm::EventSetup const& iSetup)
{
   // Write jets and constitutents. Will use fjJets_, inputs_                          
   // and fjClusterSeq_                                                                
   switch( jetTypeE ) {
   case JetType::CaloJet :
      writeJets<reco::CaloJet>( iEvent, iSetup);
      writeBkgJets<reco::CaloJet>( iEvent, iSetup);
      break;
   case JetType::PFJet :
     writeJets<reco::PFJet>( iEvent, iSetup);
      writeBkgJets<reco::PFJet>( iEvent, iSetup);
      break;
   case JetType::GenJet :
     writeJets<reco::GenJet>( iEvent, iSetup);
      writeBkgJets<reco::GenJet>( iEvent, iSetup);
      break;
   case JetType::TrackJet :
     writeJets<reco::TrackJet>( iEvent, iSetup);
      writeBkgJets<reco::TrackJet>( iEvent, iSetup);
      break;
   case JetType::BasicJet :
     writeJets<reco::BasicJet>( iEvent, iSetup);
      writeBkgJets<reco::BasicJet>( iEvent, iSetup);
      break;
   default:
      edm::LogError("InvalidInput") << " invalid jet type in VirtualJetProducer\n";
      break;
   };

}

template< typename T >
void JetAlgorithmAnalyzer::writeJets( edm::Event & iEvent, edm::EventSetup const& iSetup )
{
   // produce output jet collection

   using namespace reco;

   std::auto_ptr<std::vector<T> > jets(new std::vector<T>() );
   jets->reserve(fjJets_.size());

   for (unsigned int ijet=0;ijet<fjJets_.size();++ijet) {
      // allocate this jet
      T jet;
      // get the fastjet jet
      const fastjet::PseudoJet& fjJet = fjJets_[ijet];
      // get the constituents from fastjet
      std::vector<fastjet::PseudoJet> fjConstituents =
	 sorted_by_pt(fjClusterSeq_->constituents(fjJet));
      // convert them to CandidatePtr vector
      std::vector<CandidatePtr> constituents =
	 getConstituents(fjConstituents);

      // calcuate the jet area
      double jetArea=0.0;
      if ( doAreaFastjet_ ) {
	 fastjet::ClusterSequenceArea const * clusterSequenceWithArea =
	    dynamic_cast<fastjet::ClusterSequenceArea const *>(&*fjClusterSeq_);
	 jetArea = clusterSequenceWithArea->area(fjJet);
      }

      // write the specifics to the jet (simultaneously sets 4-vector, vertex).
      // These are overridden functions that will call the appropriate
      // specific allocator.
      writeSpecific(jet,
		    Particle::LorentzVector(fjJet.px(),
					    fjJet.py(),
					    fjJet.pz(),
					    fjJet.E()),
		    vertex_,
		    constituents, iSetup);

      jet.setJetArea (jetArea);

      if(doPUOffsetCorr_){
	 jet.setPileup(subtractor_->getPileUpEnergy(ijet));
      }else{
	 jet.setPileup (0.0);
      }

      // add to the list
      jets->push_back(jet);
   }

   // put the jets in the collection
   //   iEvent.put(jets);

   // calculate rho (median pT per unit area, for PU&UE subtraction down the line
   std::auto_ptr<double> rho(new double(0.0));
   std::auto_ptr<double> sigma(new double(0.0));

   if (doRhoFastjet_) {
      fastjet::ClusterSequenceArea const * clusterSequenceWithArea =
	 dynamic_cast<fastjet::ClusterSequenceArea const *> ( &*fjClusterSeq_ );
      double mean_area = 0;
      clusterSequenceWithArea->get_median_rho_and_sigma(*fjRangeDef_,false,*rho,*sigma,mean_area);
      iEvent.put(rho,"rho");
      iEvent.put(sigma,"sigma");
   }
}

template< typename T >
void JetAlgorithmAnalyzer::writeBkgJets( edm::Event & iEvent, edm::EventSetup const& iSetup )
{
   // produce output jet collection

   using namespace reco;

   std::vector<fastjet::PseudoJet> fjFakeJets_;
   std::vector<std::vector<reco::CandidatePtr> > constituents_;
   vector<double> phiRandom;
   vector<double> etaRandom;
   vector<double> et;
   vector<double> pileUp;
   std::auto_ptr<std::vector<bool> > directions(new std::vector<bool>());
   directions->reserve(nFill_);

   phiRandom.reserve(nFill_);
   etaRandom.reserve(nFill_);
   et.reserve(nFill_);
   pileUp.reserve(nFill_);

   fjFakeJets_.reserve(nFill_);
   constituents_.reserve(nFill_);

   constituents_.reserve(nFill_);
   for(int ijet = 0; ijet < nFill_; ++ijet){
      vector<reco::CandidatePtr> vec;
      constituents_.push_back(vec);
      directions->push_back(true);
   }

   for(int ijet = 0; ijet < nFill_; ++ijet){
      phiRandom[ijet] = 2*pi*randomEngine->flat() - pi;
      etaRandom[ijet] = 2*etaMax_*randomEngine->flat() - etaMax_;
      et[ijet] = 0;
      pileUp[ijet] = 0;
   }

   for (vector<fastjet::PseudoJet>::const_iterator input_object = fjInputs_.begin (),
	   fjInputsEnd = fjInputs_.end();
	input_object != fjInputsEnd; ++input_object) {

      const reco::CandidatePtr & tower=inputs_[input_object->user_index()];
      const CaloTower* ctc = dynamic_cast<const CaloTower*>(tower.get());
      int ieta = ctc->id().ieta();
      int iphi = ctc->id().iphi();
      CaloTowerDetId id(ieta,iphi);
      const GlobalPoint& hitpoint = geo->getPosition(id);
      double towEta = hitpoint.eta();
      double towPhi = hitpoint.phi();

      for(int ir = 0; ir < nFill_; ++ir){
	 if(reco::deltaR(towEta,towPhi,etaRandom[ir],phiRandom[ir]) > rParam_) continue;

	 constituents_[ir].push_back(tower);

	 double towet = tower->et();
	 double putow = subtractor_->getPileUpAtTower(tower);
	 double etadd = towet - putow; 
	 if(avoidNegative_ && etadd < 0.) etadd = 0;
	 et[ir] += etadd;
	 pileUp[ir] += towet - etadd;
      }
   }
   cout<<"Start filling jets"<<endl;

   for(int ir = 0; ir < nFill_; ++ir){

      if(doAnalysis_) ntRandom->Fill(etaRandom[ir],phiRandom[ir],et[ir],pileUp[ir],iev_);

      if(et[ir] < 0){
	//	 cout<<"Flipping vector"<<endl;
	 (*directions)[ir] = false;
	 et[ir] = -et[ir];
      }else{
	//         cout<<"Keep vector same"<<endl;
         (*directions)[ir] = true;
      }
      cout<<"Lorentz"<<endl;

      math::PtEtaPhiMLorentzVector p(et[ir],etaRandom[ir],phiRandom[ir],0);
      fastjet::PseudoJet jet(p.px(),p.py(),p.pz(),p.energy());
      fjFakeJets_.push_back(jet);
   }

   std::auto_ptr<std::vector<T> > jets(new std::vector<T>() );
   jets->reserve(fjFakeJets_.size());
      
   for (unsigned int ijet=0;ijet<fjFakeJets_.size();++ijet) {
      // allocate this jet
      T jet;
      // get the fastjet jet
      const fastjet::PseudoJet& fjJet = fjFakeJets_[ijet];

      // convert them to CandidatePtr vector
      std::vector<CandidatePtr> constituents =
	 constituents_[ijet];

      writeSpecific(jet,
		    Particle::LorentzVector(fjJet.px(),
					    fjJet.py(),
					    fjJet.pz(),
					    fjJet.E()),
		    vertex_,
		    constituents, iSetup);
      
      // calcuate the jet area
      double jetArea=0.0;
      jet.setJetArea (jetArea);
      if(doPUOffsetCorr_){
	 jet.setPileup(pileUp[ijet]);
      }else{
	 jet.setPileup (0.0);
      }

      // add to the list
      jets->push_back(jet);
   }
  
   // put the jets in the collection
   iEvent.put(jets,"randomCones");
   iEvent.put(directions,"directions");
}

//void JetAlgorithmAnalyzer::runAlgorithm( edm::Event & iEvent, edm::EventSetup const& iSetup, std::vector<fastjet::PseudoJet>& input )
void JetAlgorithmAnalyzer::runAlgorithm( edm::Event & iEvent, edm::EventSetup const& iSetup)
{
  if ( !doAreaFastjet_ && !doRhoFastjet_) {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs_, *fjJetDefinition_ ) );
  } else {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequenceArea( fjInputs_, *fjJetDefinition_ , *fjActiveArea_ ) );
  }
  fjJets_ = fastjet::sorted_by_pt(fjClusterSeq_->inclusive_jets(jetPtMin_));

}

void JetAlgorithmAnalyzer::inputTowers( )
{
   std::vector<edm::Ptr<reco::Candidate> >::const_iterator inBegin = inputs_.begin(),
      inEnd = inputs_.end(), i = inBegin;
   for (; i != inEnd; ++i ) {
      reco::CandidatePtr input = *i;
      if (isnan(input->pt()))           continue;
      if (input->et()    <inputEtMin_)  continue;
      if (input->energy()<inputEMin_)   continue;
      if (isAnomalousTower(input))      continue;

      //Check consistency of kinematics

      const CaloTower* ctc = dynamic_cast<const CaloTower*>(input.get());
      if(ctc){
	 int ieta = ctc->id().ieta();
	 int iphi = ctc->id().iphi();

	 ntuple->Fill(ieta, input->eta(), iphi, input->phi(),input->et(),ctc->emEt(),ctc->hadEt());



	 if(abs(ieta) < 5){

	    if(0){
	       math::RhoEtaPhiVector v(1.4,input->eta(),input->phi());
	       GlobalPoint point(v.x(),v.y(),v.z());
	       //      const DetId d = geo->getClosestCell(point);
	       //      HcalDetId hd(d);
	       HcalDetId hd(0);
	       if(hd.ieta() != ieta || hd.iphi() != iphi){
		  cout<<"Inconsistent kinematics!!!   ET = "<<input->pt()<<endl;
		  cout<<"ieta candidate : "<<ieta<<" ieta detid : "<<hd.ieta()<<endl;
		  cout<<"iphi candidate : "<<iphi<<" iphi detid : "<<hd.iphi()<<endl;
	       }

	    }

	    if(0){
	       HcalDetId det(HcalBarrel,ieta,iphi,1);

	       if(geo->present(det)){
		  double eta = geo->getPosition(det).eta();
		  double phi = geo->getPosition(det).phi();

		  if(input->eta() != eta || input->phi() != phi){
		     cout<<"Inconsistent kinematics!!!   ET = "<<input->pt()<<endl;
		     cout<<"eta candidate : "<<input->eta()<<" eta detid : "<<eta<<endl;
		     cout<<"phi candidate : "<<input->phi()<<" phi detid : "<<phi<<endl;
		  }
	       }else{
		  cout<<"DetId not present in the Calo Geometry : ieta = "<<ieta<<" iphi = "<<iphi<<endl;
	       }
	    }
	 }

      }

      if (makeCaloJet(jetTypeE)&&doPVCorrection_) {
	 const CaloTower* tower=dynamic_cast<const CaloTower*>(input.get());
	 math::PtEtaPhiMLorentzVector ct(tower->p4(vertex_));
	 fjInputs_.push_back(fastjet::PseudoJet(ct.px(),ct.py(),ct.pz(),ct.energy()));
      }
      else {
	 fjInputs_.push_back(fastjet::PseudoJet(input->px(),input->py(),input->pz(),
						input->energy()));
      }
      fjInputs_.back().set_user_index(i - inBegin);
   }
}

//______________________________________________________________________________
bool JetAlgorithmAnalyzer::isAnomalousTower(reco::CandidatePtr input)
{
   if (!makeCaloJet(jetTypeE)) return false;
   const CaloTower* tower=dynamic_cast<const CaloTower*>(input.get());
   if (0==tower) return false;
   if (tower->numBadEcalCells()        >maxBadEcalCells_        ||
       tower->numRecoveredEcalCells()  >maxRecoveredEcalCells_  ||
       tower->numProblematicEcalCells()>maxProblematicEcalCells_||
       tower->numBadHcalCells()        >maxBadHcalCells_        ||
       tower->numRecoveredHcalCells()  >maxRecoveredHcalCells_  ||
       tower->numProblematicHcalCells()>maxProblematicHcalCells_) return true;
   return false;
}
//______________________________________________________________________________
void JetAlgorithmAnalyzer::copyConstituents(const vector<fastjet::PseudoJet>& fjConstituents,
					    reco::Jet* jet)
{
   for (unsigned int i=0;i<fjConstituents.size();++i)
      jet->addDaughter(inputs_[fjConstituents[i].user_index()]);
}


//______________________________________________________________________________
vector<reco::CandidatePtr>
JetAlgorithmAnalyzer::getConstituents(const vector<fastjet::PseudoJet>&fjConstituents)
{
   vector<reco::CandidatePtr> result;
   for (unsigned int i=0;i<fjConstituents.size();i++) {
      int index = fjConstituents[i].user_index();
      reco::CandidatePtr candidate = inputs_[index];
      result.push_back(candidate);
   }
   return result;
}

JetAlgorithmAnalyzer::JetType::Type
JetAlgorithmAnalyzer::JetType::byName(const string &name)
{
   const char **pos = std::find(names, names + LastJetType, name);
   if (pos == names + LastJetType) {
      std::string errorMessage="Requested jetType not supported: "+name+"\n";
      throw cms::Exception("Configuration",errorMessage);
   }
   return (Type)(pos-names);
}

const char *JetAlgorithmAnalyzer::JetType::names[] = {
   "BasicJet","GenJet","CaloJet","PFJet","TrackJet"
};


////////////////////////////////////////////////////////////////////////////////
// define as cmssw plugin
////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(JetAlgorithmAnalyzer);

