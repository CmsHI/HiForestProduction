#include <string>
#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TH3.h"
#include "TObject.h"

class TrackCorrector 
{
  public:
    TrackCorrector();
    TrackCorrector( std::string fileName );
    void load(std::string dirName);
    double getWeight( double pT, double eta, double occ );
    double getEventWeight( int M);
    double getZeroMultFrac();
    double getEventWeightEPOS( int M);
    double getZeroMultFracEPOS();
    virtual ~TrackCorrector();

  private:

    static const double trigEff[30];
    static const double trigFak[30];
    static const double zeroMFraction;
    static const double trigEffEPOS[30];
    static const double trigFakEPOS[30];
    static const double zeroMFractionEPOS;

    TFile * table;
    TH3F * rfak;
    TH3F * rsec;
    TH3F * reff;
    TH3F * rmul;
};

const double TrackCorrector::trigEff[30] = {
        0.0, 0.0448119, 0.513363, 0.77141, 0.883644, 0.926159, 0.948583, 
        0.960119, 0.966817, 0.972921, 0.976898, 0.980945, 0.983142, 0.986, 0.988247, 
        0.989673, 0.991689, 0.993215, 0.993808, 0.994555, 0.996014, 0.996399, 0.996832, 
        0.996901, 0.997528, 0.997738, 0.998142, 0.998412, 0.99869, 0.998518 };

const double TrackCorrector::trigFak[30] = {
        0.0, 0.0150659, 0.0134671, 0.00927223, 0.00814763, 0.00572825, 0.00545534, 
        0.00500624, 0.00392851, 0.00290192, 0.0027658, 0.00256728, 0.00172565, 0.00151442, 0.00191661, 
        0.0013661, 0.000786178, 0.000747828, 0.000770812, 0.000412993, 0.000491937, 0.000288172, 0.000299312, 
        0.00021031, 0.000264633, 0.000192915, 0.000148894, 0.000151408, 0.000180004, 0.000156213};

const double TrackCorrector::zeroMFraction = 0.00142393;

const double TrackCorrector::trigEffEPOS[30] = {
        0.0, 0.0367376, 0.470091, 0.716075, 0.830268, 0.884619, 0.91243, 
        0.929862, 0.941618, 0.951745, 0.959496, 0.965881, 0.971707, 0.976154, 0.978857, 
        0.983329, 0.985951, 0.989692, 0.991847, 0.993275, 0.99445, 0.995134, 0.997041, 
        0.997154, 0.997692, 0.998125, 0.998405, 0.998711, 0.999165, 0.999453 };

const double TrackCorrector::trigFakEPOS[30] = {
        0.0, 0.0255556, 0.0299233, 0.021994, 0.0175246, 0.0131437, 0.0105594, 
        0.00910782, 0.00750184, 0.00538194, 0.00502911, 0.00393098, 0.00322919, 0.00251204, 0.00205294, 
        0.00182698, 0.00121582, 0.000822846, 0.000699091, 0.000763068, 0.000475838, 0.000478164, 0.000217806, 
        0.000221613, 0.000137669, 0.00011212, 2.8529e-05, 0.000143349, 5.76037e-05, 0.0};

const double TrackCorrector::zeroMFractionEPOS = 0.00398087;

TrackCorrector::TrackCorrector( std::string fileName ) 
{
 table = new TFile(fileName.c_str());
}

TrackCorrector::TrackCorrector()
{
  table = new TFile("trackCorrections_HIJING_DiJet120_v0.root");
}

void
TrackCorrector::load(std::string dirName)
{
  std::string cn[6] = {"heff","hrec","hsim","hfak","hmul","hsec"};
//  std::map<std::string,TH3F*> corr3D;

//  for(int i=0;i<6;i++) 
//  {
//    corr3D[cn[i]] = (TH3F*) table->Get(Form("%s/%s3D",dirName.c_str(),cn[i].c_str()));
//  }
  TH3F * corr3Deff = (TH3F*) table->Get(Form("%s/heff3D",dirName.c_str()));
  TH3F * corr3Dsim = (TH3F*) table->Get(Form("%s/hsim3D",dirName.c_str()));
  TH3F * corr3Drec = (TH3F*) table->Get(Form("%s/hrec3D",dirName.c_str()));
  TH3F * corr3Dfak = (TH3F*) table->Get(Form("%s/hfak3D",dirName.c_str()));
  TH3F * corr3Dsec = (TH3F*) table->Get(Form("%s/hsec3D",dirName.c_str()));
  TH3F * corr3Dmul = (TH3F*) table->Get(Form("%s/hmul3D",dirName.c_str()));

  reff = (TH3F*) corr3Deff->Clone("reff");
  rfak = (TH3F*) corr3Dfak->Clone("rfak");
  rmul = (TH3F*) corr3Dmul->Clone("rmul");
  rsec = (TH3F*) corr3Dsec->Clone("rsec");

  reff->Divide(corr3Deff,corr3Dsim,1,1,"B");
  rmul->Divide(corr3Dmul,corr3Dsim,1,1,"B");
  rfak->Divide(corr3Dfak,corr3Drec,1,1,"B");
  rsec->Divide(corr3Dsec,corr3Drec,1,1,"B");
}

  
TrackCorrector::~TrackCorrector()
{  
  delete table;  
}

double
TrackCorrector::getWeight(double pT, double eta, double occ ) 
{
  double eff = reff->GetBinContent(
                  reff->GetXaxis()->FindBin(eta),
                  reff->GetYaxis()->FindBin(pT),
                  reff->GetZaxis()->FindBin(occ) );
  if( eff >= 0.9999 || eff <= 0.0001) eff = 1;

  double fak = rfak->GetBinContent(
                  rfak->GetXaxis()->FindBin(eta),
                  rfak->GetYaxis()->FindBin(pT),
                  rfak->GetZaxis()->FindBin(occ) );
  if( fak >= 0.9999 || fak <= 0.0001) fak = 0;

  double sec = rsec->GetBinContent(
                  rsec->GetXaxis()->FindBin(eta),
                  rsec->GetYaxis()->FindBin(pT),
                  rsec->GetZaxis()->FindBin(occ) );
  if( sec >= 0.9999 || sec <= 0.0001) sec = 0;

  double mul = rmul->GetBinContent(
                  rmul->GetXaxis()->FindBin(eta),
                  rmul->GetYaxis()->FindBin(pT),
                  rmul->GetZaxis()->FindBin(occ) );
  if( mul >= 0.9999 || mul <= 0.0001) mul = 0;


  return (1. - fak ) * ( 1. - sec ) / eff  / (1. + mul );
}

double 
TrackCorrector::getEventWeight( int M )
{
  if( M<1 || M>29) return 1;
  return (1. - trigFak[M] ) / trigEff[M];
}

double 
TrackCorrector::getZeroMultFrac()
{
  return zeroMFraction;
}

double
TrackCorrector::getEventWeightEPOS( int M )
{
  if( M<1 || M>29) return 1;
  return (1. - trigFakEPOS[M] ) / trigEffEPOS[M];
}

double
TrackCorrector::getZeroMultFracEPOS()
{
  return zeroMFractionEPOS;
}

