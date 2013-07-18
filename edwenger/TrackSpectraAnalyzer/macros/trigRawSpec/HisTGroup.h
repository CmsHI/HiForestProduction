#ifndef HisTGROUP_H
#define HisTGROUP_H

#include <map>
#include <cassert>
#include "TString.h"
#include "TFile.h"
#include <iostream>

template <typename TData>
class HisTGroup
{
  public:
    // methods
    HisTGroup(TString name,
	Int_t xn=-1, Double_t xmin=-1, Double_t xmax=-1,
	Int_t yn=-1, Double_t ymin=-1, Double_t ymax=-1);
    void Add(TData * h1, TString iname,Double_t sc=-1);
    void Add1D(TString iname);
    void Add2D(TString iname);
    void Add(TFile * inFile, TString hname,TString iname, Double_t sc=-1);
    TData * Average();
    TData * Divide(TString iname1, TString iname2);
    void Print() const;
    void Save();
    TData * Sum();

    // Accessors
    TData * H(TString iname) const {
      if (hm_.find(iname)==hm_.end()) {
	cout << name_ << ": " << iname << " not found" << endl;
	assert(false);
      }
      if (!hm_.find(iname)->second) cout << name_ << ": " << iname << " not found" << endl;
      assert(hm_.find(iname)->second);
      return hm_.find(iname)->second;
    }
    TData * R(TString iname) const {
      if (hr_.find(iname)==hr_.end()) {
	cout << name_ << ": " << iname << " not found" << endl;
	assert(false);
      }
      if (!hr_.find(iname)->second) cout << name_ << ": " << iname << " not found" << endl;
      assert(hr_.find(iname)->second);
      return hr_.find(iname)->second;
    }
    const TData * GetH(TString iname) const {
      if (hm_.find(iname)==hm_.end()) {
	cout << name_ << ": " << iname << " not found" << endl;
	assert(false);
      }
      if (!hm_.find(iname)->second) cout << name_ << ": " << iname << " not found" << endl;
      assert(hm_.find(iname)->second);
      return hm_.find(iname)->second;
    }
    const TData * GetR(TString iname) const {
      if (hr_.find(iname)==hr_.end()) {
	cout << name_ << ": " << iname << " not found" << endl;
	assert(false);
      }
      if (!hr_.find(iname)->second) cout << name_ << ": " << iname << " not found" << endl;
      assert(hr_.find(iname)->second);
      return hr_.find(iname)->second;
    }

    // data
    TString name_;
    std::map<TString,Double_t> scales_;

    // group properties
    Int_t xnbins_;
    Double_t xmin_;
    Double_t xmax_;
    Int_t ynbins_;
    Double_t ymin_;
    Double_t ymax_;

  protected:
    // group private members
    std::map<TString,TData*> hm_;
    std::map<TString,TData*> hr_;
};

// === Constructors ===
template <typename TData>
HisTGroup<TData>::HisTGroup(TString name,
    Int_t xn, Double_t xmin, Double_t xmax,
    Int_t yn, Double_t ymin, Double_t ymax) :
  name_(name),
  xnbins_(xn),
  xmin_(xmin),
  xmax_(xmax),
  ynbins_(yn),
  ymin_(ymin),
  ymax_(ymax)
{
  TData::SetDefaultSumw2();
}

// === Adding to Group ===
template <typename TData>
void HisTGroup<TData>::Add(TData * h1, TString iname, Double_t sc)
{
  if (!h1) cout << name_ << ": " << iname << " hist not found" << endl;
  assert(h1);
  scales_[iname] = sc;

  if (sc<0) {
    hm_[iname] = h1;
  }

  if (sc>0) {
    TString schname = TString("sc")+h1->GetName();
    TData * hsc = (TData*)h1->Clone(schname);
    hsc->Scale(sc);
    hm_[iname] = hsc;
  }
}

template <typename TData>
void HisTGroup<TData>::Add1D(TString iname)
{
  TString hname = TString("h")+name_+iname;
  TData * h1 = 0;
  if (xnbins_>0)
    h1 = new TData(hname,hname,xnbins_,xmin_,xmax_);
  assert(h1);
  Add(h1,iname);
}

template <typename TData>
void HisTGroup<TData>::Add2D(TString iname)
{
  TString hname = TString("h")+name_+iname;
  TData * h1 = 0;
  if (xnbins_>0 && ynbins_>0)
    h1 = new TData(hname,hname,xnbins_,xmin_,xmax_,ynbins_,ymin_,ymax_);
  assert(h1);
  Add(h1,iname);
}

template <typename TData>
void HisTGroup<TData>::Add(TFile * inFile, TString hname,TString iname, Double_t sc)
{
  assert(inFile);
  TData * h1;
  inFile->GetObject(hname.Data(),h1);
  Add(h1,iname,sc);
}

// === Relation: Average ===
template <typename TData>
TData * HisTGroup<TData>::Average()
{
  // check if average already taken
  if (hr_.find("Ave")!=hr_.end()) return hr_["Ave"];

  if (hr_.find("Sum")==hr_.end()) Sum();
  hr_["Ave"] = (TData*)hr_["Sum"]->Clone("h"+name_+"_Ave");
  hr_["Ave"]->Scale(1./hm_.size());
  return hr_["Ave"];
}

template <typename TData>
TData * HisTGroup<TData>::Divide(TString iname1, TString iname2)
{
  TString quoIName = iname1+"Div"+iname2;
  TString quoHName = "h"+name_+"_"+quoIName;
  hr_[quoIName] = (TData*)GetH(iname1)->Clone(quoHName);
  hr_[quoIName]->Divide(GetH(iname2));
  return hr_[quoIName];
}

// === Helper Functions ===
template <typename TData>
void HisTGroup<TData>::Print() const
{
  std::cout << "group: " << name_ << std::endl;
  std::cout << "xnbins: " << xnbins_ << " xmin: " << xmin_ << " xmax: " << xmax_ << std::endl;
  if (ynbins_>0)
    std::cout << "ynbins: " << ynbins_ << " ymin: " << ymin_ << " ymax: " << ymax_ << std::endl;
  cout << "hm size: " << hm_.size() << endl;
  cout << "hr size: " << hr_.size() << endl;
  if (hr_.find("Sum")!=hr_.end()) {
    std::cout << "hrSum: " << GetR("Sum")->GetName() << " has: " << GetR("Sum")->GetEntries() << std::endl;
  }
  if (hr_.find("Ave")!=hr_.end()) {
    std::cout << "hrAve: " << GetR("Ave")->GetName() << " has: " << GetR("Ave")->GetEntries() << std::endl;
  }
  for (typename std::map<TString, TData*>::const_iterator 
      iter=hr_.begin(); iter != hr_.end(); ++iter) {
    std::cout << "hr: " << iter->first << " " << iter->second << std::endl;
  }
}

template <typename TData>
void HisTGroup<TData>::Save()
{
  for (typename std::map<TString, TData*>::iterator 
      iter=hm_.begin(); iter != hm_.end(); ++iter) {
    iter->second->Write();
    //std::cout << "Save: " << iter->second->GetName() << std::endl;
  }
}

// === Relation: Sum ===
template <typename TData>
TData * HisTGroup<TData>::Sum()
{
  // check if Sum already taken
  if (hr_.find("Sum")!=hr_.end()) return hr_["Sum"];

  typename std::map<TString, TData*>::iterator iter; // typename keyword needed here, b/c comiple is in doubt whether TData is a type
  for (iter=hm_.begin(); iter != hm_.end(); ++iter) {
    if (iter==hm_.begin()) {
      hr_["Sum"] = (TData*)iter->second->Clone("h"+name_+"_Sum");
      assert(hr_["Sum"]);
      //std::cout << "first to add: " << iter->first << endl;
    } else {
      //std::cout << "add more: " << iter->first << endl;
      hr_["Sum"]->Add(iter->second);
    }
  }
  return hr_["Sum"];
}

#endif
