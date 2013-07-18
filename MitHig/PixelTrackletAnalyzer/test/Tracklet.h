#include <vector>
#include "RecoHit.h"
using namespace std;

#define PI 3.141592653589

class Tracklet
{
   public:
      Tracklet(double eta1,double eta2,double phi1,double phi2);
      ~Tracklet(){};

      double eta1() { return eta1_; }
      double eta2() { return eta2_; }
      double phi1() { return phi1_; }
      double phi2() { return phi2_; }

      double deta() { return eta1_-eta2_; }
      double dphi();
      double dR();
      double dR2();
      
 
      void   setIt1(int i) { it1_ = i; }
      void   setIt2(int i) { it2_ = i; }
      int    getIt1() { return it1_; }
      int    getIt2() { return it2_; }

      void   setId(int i) { id_ = i; }
      void   setSId(int i) { sid_ = i; }
      void   setType(int i) { ptype_ = i; }
      void   setId1(int i) { id1_ = i; }
      void   setId2(int i) { id2_ = i; }
      int    getId() { return id_; }
      int    getSId() { return sid_; }
      int    getType() { return ptype_; }
      int    getId1() { return id1_; }
      int    getId2() { return id2_; }


   private:

      double eta1_;
      double eta2_;

      double phi1_;
      double phi2_;

      int it1_;   // first iterator
      int it2_;   // second iterator
      
      int sid_;   // signal exist?
      int ptype_;  // process type
      int id_;
      int id1_;
      int id2_;
      
};

Tracklet::Tracklet(double eta1, double eta2, double phi1, double phi2) {

   eta1_ = eta1;
   eta2_ = eta2;

   phi1_ = phi1;
   while (phi1_>2*PI) phi1_-=2*PI;

   phi2_ = phi2;
   while (phi2_>2*PI) phi2_-=2*PI;
}

double Tracklet::dphi() 
{
   double dphi=phi1_-phi2_;

   if (dphi>0){
      while (dphi>2*PI) dphi-=2*PI;
      if (dphi>PI) dphi=2*PI-dphi;
   } else {
      while (dphi<-2*PI) dphi+=2*PI;
      if (dphi<-PI) dphi=-2*PI-dphi;
   }


   return dphi; 
}

double Tracklet::dR() 
{
   return sqrt(dR2());
}

double Tracklet::dR2() 
{
   double dPhi=dphi();
   double dEta=deta();

   return dPhi*dPhi+dEta*dEta;
}

bool compareDeltaR(Tracklet a,Tracklet b) { return fabs(a.dR2())<fabs(b.dR2());}
bool compareDeltaEta(Tracklet a,Tracklet b) {return fabs(a.deta())<fabs(b.deta());}

vector<Tracklet> cleanTracklets(vector<Tracklet> input, int matchNumber, SelectionCriteria cuts);
vector<Tracklet> recoProtoTracklets(vector<RecoHit> firstLayerHits, vector<RecoHit> secondLayerHits);

vector<Tracklet> recoProtoTracklets(vector<RecoHit> hits1, vector<RecoHit> hits2)
{
  vector<Tracklet> protoTracklets;

  for (int i = 0; i < (int) hits1.size(); i++)
    {
      for (int j = 0; j < (int) hits2.size(); j++)
	{
	  Tracklet mytracklet(hits1[i].eta,hits2[j].eta,hits1[i].phi,hits2[j].phi);
	  mytracklet.setIt1(i);
	  mytracklet.setIt2(j);
	  protoTracklets.push_back(mytracklet);
	}
    }

  return protoTracklets;
}

vector<Tracklet> cleanTracklets(vector<Tracklet> input, int matchNumber,SelectionCriteria cuts)
{
  vector<Tracklet> output;

  if(cuts.useDeltaPhi_)
    sort( input.begin() , input.end() , compareDeltaR);
  else
    sort( input.begin() , input.end() , compareDeltaEta);

  if (cuts.verbose_) {
    for (unsigned int i = 0; i < input.size(); i++)
      {
	cout <<input[i].deta()<<" "<<input[i].getIt1()<<" "<<input[i].getIt2()<<endl;
      }
  }

  int used1[1000];
  int used2[1000];

  for (int i=0;i<1000;i++) {
    used1[i]=0;
    used2[i]=0;
  }
  if (cuts.verbose_) cout<<"Printing Hits"<<endl;

  for (unsigned int i = 0; i < input.size(); i++){


    if(cuts.useDeltaPhi_)
      if (cuts.verbose_) cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta R : "<<input[i].dR()<<endl;
      else
	if (cuts.verbose_) cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta Eta : "<<input[i].deta()<<endl;

    int i1=input[i].getIt1();
    int i2=input[i].getIt2();

    if (used1[i1]==0&&used2[i2]==matchNumber) {
      Tracklet tmp = input[i];
      output.push_back(tmp);
      used1[i1]++;
      if (cuts.checkSecondLayer_) used2[i2]++;
    }
    if (used1[i1]==0&&used2[i2]<matchNumber) {
      if (cuts.checkSecondLayer_) used2[i2]++;
    }
  }
  if (cuts.verbose_) {
    cout <<"Output:"<<endl;
    for (unsigned int i = 0; i < output.size(); i++)
      {
	cout <<output[i].deta()<<" "<<output[i].getIt1()<<" "<<output[i].getIt2()<<endl;
      }
  }

  return output;
}











