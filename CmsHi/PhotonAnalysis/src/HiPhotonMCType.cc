#include "CmsHi/PhotonAnalysis/interface/HiPhotonMCType.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

using namespace edm;
using namespace reco;

#define PI 3.141592653589793238462643383279502884197169399375105820974945

HiPhotonMCType::HiPhotonMCType ()
{
} 

HiPhotonMCType::EType HiPhotonMCType::getType(const reco::Candidate &part)
{
   using namespace edm;
   using namespace reco;

   if (part.pdgId()!=22) return kInvalid;  // Not a photon


   if (part.mother()==0) {
      if (part.status()==1) { 
         // Given particle is not a mother and has no mother.   
         return kInvalid;
      }
      return kLO;
   }

   int mId = part.mother()->pdgId();

   if (mId==kGamma) {
      return kDirect;
   }
   
   if (mId==kGluon || mId==kUq || mId==kDq || mId==kSq || mId==kCq 
                   || mId==kBq ||mId==kTq) {
      // It's a parton
      UInt_t nSiblings = part.mother()->numberOfDaughters();
      if (nSiblings>2) return kFrag;
      else return kBrems;      
   }

   // check strange "particles"
   if (mId==110 || mId==990 || mId>1000000)
      return kUnknown;
   
   switch (mId) {
      case 1103:
      case 2101: 
      case 2103:
      case 2203:
      case 3101:
      case 3103:
      case 3201:
      case 3203:
      case 3303:
      case 4101:
      case 4103:
      case 4201:
      case 4203:
      case 4301:
      case 4303:
      case 4403:
      case 5101:
      case 5103:
      case 5201:
      case 5203:
      case 5301:
      case 5303:
      case 5401:
      case 5403:
      case 5503:
         return kUnknown;
      default:
         break;
   }

   return kDecay;

}

