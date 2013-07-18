
#ifndef TrackSpectraAnalyzer_TrackSpectraAnalyzerHist
#define TrackSpectraAnalyzer_TrackSpectraAnalyzerHist

#include <vector>
namespace edm { class ParameterSet ; }
class TFile;
class TTree;
class TH2F;
class TH3F;


class TrackSpectraAnalyzerHist
{
 public:
   TrackSpectraAnalyzerHist(const edm::ParameterSet& pset);
   ~TrackSpectraAnalyzerHist();
   
   void prepareHist();
   void fillTrack2DHist(float eta, float pt, int type);
   void fillJet2DHist(float jeta, float jet, int tbit);
   void fillTrackJet3DHist(float eta, float pt, float jet, int tbit);
   void writeHist();
   
 private:
   TFile * outputFile;

   bool isGEN_;

   double  etaBinW_;
   double  etaMax_;    

   double ptBinW_;
   double ptMin_, ptMax_;

   double jetaBinW_;
   double jetaMax_;

   double jetBinW_;
   double jetMin_, jetMax_;


   std::vector<float> etaBins, ptBins;
   std::vector<float> jetaBins, jetBins;
   std::vector<float> multBins;


   TH2F* hTrkPtEta;
   TH2F* hGenTrkPtEta;

   TH2F* hJetEtEtaMB;
   TH2F* hJetEtEta6U;
   TH2F* hJetEtEta15U;
   TH2F* hJetEtEta30U;
   TH2F* hJetEtEta50U;

   TH3F* hTrkPtEtaMB;
   TH3F* hTrkPtEta6U;
   TH3F* hTrkPtEta15U;
   TH3F* hTrkPtEta30U;
   TH3F* hTrkPtEta50U; 


};


#endif
