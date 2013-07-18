#ifndef MNguyen_InclusiveJetAnalyzer_inclusiveJetAnalyzer_
#define MNguyen_InclusiveJetAnalyzer_inclusiveJetAnalyzer_

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"



#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

//

/**\class InclusiveJetAnalyzer 

\author Matt Nguyen
\date   November 2010
*/




class InclusiveJetAnalyzer : public edm::EDAnalyzer {
 public:

  explicit InclusiveJetAnalyzer(const edm::ParameterSet&);

  ~InclusiveJetAnalyzer();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

  virtual void beginJob();

  void fillL1Bits(const edm::Event &iEvent);

  void fillHLTBits(const edm::Event &iEvent);

  template <typename TYPE>
    void                          getProduct(const std::string name, edm::Handle<TYPE> &prod,
					     const edm::Event &event) const;    
  template <typename TYPE>
    bool                          getProductSafe(const std::string name, edm::Handle<TYPE> &prod,
						 const edm::Event &event) const;
  

 private:
  


  edm::InputTag   jetTag_, vtxTag_, genjetTag_, eventInfoTag_, L1gtReadout_; 


  /// verbose ?
  bool   verbose_;

  bool useCentrality_;
  bool useVtx_;
  bool useJEC_;
  bool isMC_;


  TTree *t;
  edm::Service<TFileService> fs1;

  CentralityProvider * centrality_;



  std::string                   hltResName_;         //HLT trigger results name
  std::vector<std::string>      hltProcNames_;       //HLT process name(s)
  std::vector<std::string>      hltTrgNames_;        //HLT trigger name(s)

  std::vector<int>              hltTrgBits_;         //HLT trigger bit(s)
  std::vector<bool>             hltTrgDeci_;         //HLT trigger descision(s)
  std::vector<std::string>      hltTrgUsedNames_;    //HLT used trigger name(s)
  std::string                   hltUsedResName_;     //used HLT trigger results name



  static const int MAXJETS = 50000;
  static const int MAXHLTBITS = 500000;


  struct JRA{
    
    int nref;
    int run;
    int evt;
    int lumi;
    int bin;
    float vx, vy, vz;
    float b;
    float hf;

    float rawpt[MAXJETS];
    float jtpt[MAXJETS];
    float jteta[MAXJETS];
    float jtphi[MAXJETS];
    float jty[MAXJETS];
    float refpt[MAXJETS];
    float refeta[MAXJETS];
    float refphi[MAXJETS];
    float refy[MAXJETS];
    float refdphijt[MAXJETS];
    float refdrjt[MAXJETS];
    float refparton_pt[MAXJETS];
    float refparton_flavor[MAXJETS];

    float pthat;
    int ngen;
    int genmatchindex[MAXJETS];
    float genpt[MAXJETS];
    float geneta[MAXJETS];
    float genphi[MAXJETS];
    float geny[MAXJETS];
    float gendphijt[MAXJETS];
    float gendrjt[MAXJETS];

    // hlt
    int nHLTBit;
    bool hltBit[MAXHLTBITS];

    // l1
    int nL1TBit;
    bool l1TBit[MAXHLTBITS];
    int nL1ABit;
    bool l1ABit[MAXHLTBITS];

  };

  JRA jets_;

};

#endif
