#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"

using namespace std;

#endif

#include "CommonParameters.h"

static int iPlot = -99;

int centMode = 1;
int entryMode = 0;

//---------------------------------------------------------------------
void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, 
			  const Int_t rows, const Float_t leftOffset=0.,
                          const Float_t bottomOffset=0., 
			  const Float_t leftMargin=0.2, 
			  const Float_t bottomMargin=0.2,
                          const Float_t edge=0.05);

void plotBalance(int cbin = 0,
		 TString infname = "file1.root",
		 TString refname = "file2.root",
		 TString mixname = "file3.root",
		 bool drawXLabel = false,
		 bool drawLeg = false);


void drawText(const char *text, float xp, float yp, int size = 22);
//void drawText(const char *text, float xp1, float yp1, float xp2, float yp2, int size);

//--------------------------------------------------------------
// drawPatch() is a crazy way of removing 0 in the second and third 
// pad which is partially shown due to no margin between the pads
// if anybody has a better way of doing it let me know! - Andre
//--------------------------------------------------------------
void drawPatch(float x1, float y1, float x2, float y2); 
//---------------------------------------------------------------------

void plotFigure(int iplot = 0){
  TString infname = "/d101/yetkin/analysis/d0204/ntuple_data_pPb_akPu3PF_forest71_20130204_01.root";
  TString refname = "/d101/yetkin/analysis/d0205/ntuple_data_PbPb_akPu3PF_forest71_20130205_01.root";
  TString mixname = "/d101/yetkin/analysis/d0204/hijing.root";

  infname = "ntuple_data_pPb_akPu3PF_forest71_20130204_01.root";
  refname = "ntuple_data_PbPb_akPu3PF_forest71_20130205_01.root";
  mixname = "hijing.root";

  string hfNames[] = {
    "30<E_{T}^{HF[#eta > 4]}<70",
    "20<E_{T}^{HF[#eta > 4]}<30",
    "15<E_{T}^{HF[#eta > 4]}<20",
    "10<E_{T}^{HF[#eta > 4]}<15",
    "5<E_{T}^{HF[#eta > 4]}<10",
    "0<E_{T}^{HF[#eta > 4]}<5"
  };


  string ntrkNames[] = {
    "180 #leq N_{trk}^{offline}",
    "150 #leq N_{trk}^{offline} < 180",
    "110 #leq N_{trk}^{offline} < 150",
    "90 #leq N_{trk}^{offline} < 110",
    "60 #leq N_{trk}^{offline} < 90",
    "N_{trk}^{offline} < 60"
  };

  string *binNames = hfNames;
  if(centMode == 0) binNames = ntrkNames;

  iPlot = iplot;

  TH1::SetDefaultSumw2();
	
   TCanvas *c1 = new TCanvas("c1","",1050,700);
   makeMultiPanelCanvas(c1,3,2,0.0,0.0,0.2,0.2,0.02);
   TLatex *jetf_PbPb;

   if(iPlot == 3 || iPlot == 8 || iPlot == 11 || iPlot == 12 || iPlot == 13 || iPlot == 14){
     for(int i = 0; i < 6; ++i){
       c1->cd(i+1)->SetLogy();
     }
   }


  c1->cd(1);

  for(int i = 0; i < 6; ++i){
    c1->cd(6-i);
    plotBalance(i,infname,refname,mixname,i==1,i==0);
    double y1 = 0.85;
    double x1 = 0.5;

    if(i == 2 || i == 5) x1 = 0.6;
    drawText(binNames[i].data(),x1,y1);
  }



  c1->cd(1);

  drawText("CMS Preliminary",0.22,0.9,18);   
  drawText(Form("pPb %s",LUM_pPb),0.22,0.8,18);
  if(iPlot == 0 || iPlot == 3) drawText(Form("PbPb %s",LUM_PbPb),0.22,0.7,18);


  c1->cd(2);
  drawText("anti-k_{T} (R=0.3)",0.02,0.9,18);
  drawText("PU subt. PF jets",0.02,0.8,18);
  drawText("|#eta|<2",0.02,0.7,18);

  c1->cd(3);
  drawText(Form("p_{T,1} > %d GeV/c",leadCut),0.02,0.9,18);
  drawText(Form("p_{T,2} > %d GeV/c",subleadCut),0.02,0.8,18);
  drawText("#Delta#phi_{12} > #frac{2}{3}#pi",0.02,0.7,18);


  const char* date = "20130127";

  string figures[150] = {"imbalance","etaDijet","",    "dphi","",
		      "","",        "", "pudifference", "ntrk",
		      "","pu1","pu2","trkMax1","trkMax2",
                      "","","","","",
		      "","eta1","eta2"
  };

  string formats[] = {"gif","pdf","eps","C"};

  for(int it = 0; it < 4; ++it){
    c1->Print(Form("./fig/dijet_%s_lead%d_sub%d_all_cent_%s%s.%s",figures[iPlot].data(),leadCut,subleadCut,date,subtract?"_subt":"",formats[it].data()));
  }  

}

void plotBalance(int cbin,
		 TString infname,
		 TString pythia,
		 TString mix,
		 bool drawXLabel,
		 bool drawLeg)
{

  bool showPbPb = 0;
  if(iPlot == 0 || iPlot == 3) showPbPb = 1;

  if(iPlot != 3) dijet = dijet&&deltaPhi;

  TCut centHF("");
  TCut centNtrk("");

  if(centMode == 0){
    cout<<"Obsolete centrality usage."<<endl;
    return;
    if(cbin==0) centHF = "bin>=20 && bin<26";
    if(cbin==1) centHF = "bin>=26 && bin<27";
    if(cbin==2) centHF = "bin>=27 && bin<28";
    if(cbin==3) centHF = "bin>=28 && bin<29";
    if(cbin==4) centHF = "bin>=29 && bin<31";
    if(cbin==5) centHF = "bin>=31";

    if(cbin==0) centNtrk = "ntrk >= 180";
    if(cbin==1) centNtrk = "150 <= ntrk && ntrk < 180";
    if(cbin==2) centNtrk = "110 <= ntrk && ntrk < 150";
    if(cbin==3) centNtrk = "90 <= ntrk && ntrk < 150";
    if(cbin==4) centNtrk = "60 <= ntrk && ntrk < 90";
    if(cbin==5) centNtrk = "ntrk < 60";
  }

  if(centMode == 1){
    centHF = hfPlusBins[cbin];
    centNtrk = centHF;
  }

  cout<<"plotting ntrk bin : "<<cbin<<endl;

  // open the data file

  if(iPlot > 100) entryMode = 1;

  string treeName = "ntdijet";
  string wName = "ntw";
  if(entryMode == 1){
    treeName = "ntjet";
    wName = "ntjetw";
  }

  if(entryMode == 2){
    treeName = "nttrk";
    wName = "nttrkw";
  }

  TFile *inf = new TFile(infname.Data());
  TTree *nt =(TTree*)inf->FindObjectAny(treeName.data());
  TTree *ntevt =(TTree*)inf->FindObjectAny("ntevt");
  if(entryMode < 1) nt->AddFriend(ntevt);

  // open the pythia (MC) file
  TFile *infReference = new TFile(pythia.Data());
  TTree *ntReference;
  TTree *ntevtReference;

  ntReference = (TTree*) infReference->FindObjectAny(treeName.data());
  ntevtReference = (TTree*) infReference->FindObjectAny("ntevt");
  if(entryMode < 1) ntReference->AddFriend(ntevtReference);
  
  // open the datamix file
  TFile *infMix = new TFile(mix.Data());
  TTree *ntMix =(TTree*)infMix->FindObjectAny(treeName.data());
  TTree *ntevtMix =(TTree*)infMix->FindObjectAny("ntevt");
  //  TFile *infW = new TFile("weights_hydjet.root");
  //  TTree *ntw =(TTree*)infW->FindObjectAny("ntw");

  TTree *ntw =(TTree*)infMix->FindObjectAny(wName.data());
  ntMix->AddFriend(ntw);
  if(entryMode < 1) ntMix->AddFriend(ntevtMix);

  int Nbin = 10;
  double max = 1.;
  double min = 0;

  if(iPlot == 1 || iPlot == 21 || iPlot == 22){
    Nbin = 20;
    max = 2;
    min = -2;
  };

  if(iPlot == 3){
    Nbin = 30;
    max = pi;
  };

  if(iPlot == 8){
    Nbin = 50;
    max = 25;
    min = -25;
  };

  if(iPlot == 9){
    Nbin = 50;
    max = 250;
  };

  if(iPlot == 11 || iPlot == 12){
    Nbin = 50;
    max = 25;
  };

  if(iPlot == 13 || iPlot == 14){
    Nbin = 75;
    max = 150;
  }

  if(iPlot == 101){
    Nbin = 20;
    max = 2;
    min = -2;
  }

  // projection histogram
  TH1D *h = new TH1D(Form("h",cbin),"",Nbin,min,max);
  TH1D *hReference = new TH1D(Form("hReference",cbin),"",Nbin,min,max);
  TH1D *hDataMix = new TH1D(Form("hDataMix",cbin),"",Nbin,min,max);

  TH1D *hB = new TH1D(Form("hB",cbin),"",Nbin,min,max);
  TH1D *hReferenceB = new TH1D(Form("hReferenceB",cbin),"",Nbin,min,max);
  TH1D *hDataMixB = new TH1D(Form("hDataMixB",cbin),"",Nbin,min,max);

  TH1D *hFull = new TH1D("hFull","",Nbin,min,max);
  TH1D *hReferenceFull = new TH1D("hReferenceFull","",Nbin,min,max);
  TH1D *hDataMixFull = new TH1D("hDataMixFull","",Nbin,min,max);

  TH1D* hNorm = new TH1D("hNorm","",1000,0,1000);
  TH1D* hNormReference = new TH1D("hNormReference","",1000,0,1000);
  TH1D* hNormDataMix = new TH1D("hNormDataMix","",1000,0,1000);

  hB->SetLineStyle(2);
  hReferenceB->SetLineStyle(2);
  hDataMixB->SetLineStyle(2);

  //  ntReference->SetAlias("pt1","et1");
  //  ntReference->SetAlias("pt2","et2");

  nt->SetAlias("pt1","jtpt1");
  nt->SetAlias("pt2","jtpt2");
  nt->SetAlias("eta1","jteta1");
  nt->SetAlias("eta2","jteta2");
  nt->SetAlias("phi1","jtphi1");
  nt->SetAlias("phi2","jtphi2");
  ntMix->SetAlias("pt1","jtpt1");
  ntMix->SetAlias("pt2","jtpt2");
  ntMix->SetAlias("eta1","jteta1");
  ntMix->SetAlias("eta2","jteta2");
  ntMix->SetAlias("phi1","jtphi1");
  ntMix->SetAlias("phi2","jtphi2");
  ntReference->SetAlias("pt1","jtpt1");
  ntReference->SetAlias("pt2","jtpt2");
  ntReference->SetAlias("eta1","jteta1");
  ntReference->SetAlias("eta2","jteta2");
  ntReference->SetAlias("phi1","jtphi1");
  ntReference->SetAlias("phi2","jtphi2");

  nt->SetAlias("cent","hfp");
  ntMix->SetAlias("cent","hfp");
  ntReference->SetAlias("cent","hfp");

  nt->SetAlias("dphi","acos(cos(phi1-phi2))");
  ntMix->SetAlias("dphi","acos(cos(phi1-phi2))");
  ntReference->SetAlias("dphi","acos(cos(phi1-phi2))");

  if(iPlot == 0){
    nt->SetAlias("var","pt2/pt1");
    ntReference->SetAlias("var","pt2/pt1");
    ntMix->SetAlias("var","pt2/pt1");
  }

  if(iPlot == 1){
    nt->SetAlias("var","(eta1+eta2)/2");
    ntReference->SetAlias("var","(eta1+eta2)/2");
    ntMix->SetAlias("var","(eta1+eta2)/2");
  }

  if(iPlot == 3){
    nt->SetAlias("var","acos(cos(phi1-phi2))");
    ntReference->SetAlias("var","acos(cos(phi1-phi2))");
    ntMix->SetAlias("var","acos(cos(phi1-phi2))");
  }

  if(iPlot == 8){
    nt->SetAlias("var","pu1-pu2");
    ntReference->SetAlias("var","pu1-pu2");
    ntMix->SetAlias("var","pu1-pu2");
  }

  if(iPlot == 11){
    nt->SetAlias("var","pu1");
    ntReference->SetAlias("var","pu1");
    ntMix->SetAlias("var","pu1");
  }

  if(iPlot == 12){
    nt->SetAlias("var","pu2");
    ntReference->SetAlias("var","pu2");
    ntMix->SetAlias("var","pu2");
  }

  if(iPlot == 13){
    nt->SetAlias("var","trkMax1");
    ntReference->SetAlias("var","trkMax1");
    ntMix->SetAlias("var","trkMax1");
  }

  if(iPlot == 14){
    nt->SetAlias("var","trkMax2");
    ntReference->SetAlias("var","trkMax2");
    ntMix->SetAlias("var","trkMax2");
  }

  if(iPlot == 9){
    nt->SetAlias("var","ntrk");
    ntReference->SetAlias("var","ntrk");
    ntMix->SetAlias("var","ntrk");
  }

  if(iPlot == 21){
    nt->SetAlias("var","eta1");
    ntReference->SetAlias("var","eta1");
    ntMix->SetAlias("var","eta1");
  }

  if(iPlot == 22){
    nt->SetAlias("var","eta2");
    ntReference->SetAlias("var","eta2");
    ntMix->SetAlias("var","eta2");
  }



  if(iPlot == 101){
    nt->SetAlias("var","eta");
    ntReference->SetAlias("var","eta");
    ntMix->SetAlias("var","eta");
  }

  if(iPlot == 51){
    nt->SetAlias("var","matchPt1/pt1");
    ntReference->SetAlias("var","matchPt1/pt1");
    ntMix->SetAlias("var","matchPt1/pt1");
  }

  if(iPlot == 52){
    nt->SetAlias("var","matchPt2/pt2");
    ntReference->SetAlias("var","matchPt2/pt2");
    ntMix->SetAlias("var","matchPt2/pt2");
  }


  nt->Draw("var>>hFull",dijet&&noise&&jetID&&centHF&&vtx); 
  nt->Draw("var>>hB",side&&noise&&jetID&&centHF&&vtx);
  nt->Draw("pt1>>hNorm",lead&&noise&&jetID&&centHF&&vtx);
   
  ntMix->Draw("var>>hDataMixFull",weight*(dijet&&jetID&&centHF&&vtx));
  ntMix->Draw("var>>hDataMixB",weight*(side&&jetID&&centHF&&vtx));
  ntMix->Draw("pt1>>hNormDataMix",weight*(lead&&jetID&&centHF&&vtx));
  
  ntReference->Draw("var>>hReferenceFull",dijet&&noise&&jetID&&centHF&&vtx);
  ntReference->Draw("var>>hReferenceB",side&&noise&&jetID&&centHF&&vtx);
  ntReference->Draw("pt1>>hNormReference",lead&&noise&&jetID&&centHF&&vtx);

  hDataMix->Add(hDataMixFull);
  h->Add(hFull);
  hReference->Add(hReferenceFull);

  if(subtract){
    hDataMix->Add(hDataMixB,-1);
    h->Add(hB,-1);
    hReference->Add(hReferenceB,-1);
  }

  hB->SetFillStyle(3005);
  hB->SetFillColor(15);

  // calculate the statistical error and normalize
  h->SetLineColor(dataColor);
  h->SetMarkerColor(dataColor);
  h->Sumw2();
  if(normLead){
    h->Scale(1./hNorm->Integral());
    hB->Scale(1./hNorm->Integral());
  }else{
    hB->Scale(1./h->Integral());
    h->Scale(1./h->Integral());
  }
  h->SetMarkerStyle(20);

  if(hReference->Integral() > 0){
    hReference->Scale(1./hNormReference->Integral());
  }

  hReference->SetLineColor(kBlue);
  hReference->SetFillColor(kAzure-8);
  hReference->SetFillStyle(3005);

  if(normLead){
    hDataMixB->Scale(1./hNormDataMix->Integral());
    hDataMix->Scale(1./hNormDataMix->Integral());
  }else{
    hDataMixB->Scale(1./hDataMix->Integral());
    hDataMix->Scale(1./hDataMix->Integral());
  }
  hDataMix->SetLineColor(mixColor);
  hDataMix->SetFillColor(mixColor);
  hDataMix->SetFillStyle(3004);
  
  hDataMix->SetMarkerSize(0);
  hDataMix->SetStats(0);

  hDataMix->GetXaxis()->SetLabelSize(22);
  hDataMix->GetXaxis()->SetLabelFont(43);
  hDataMix->GetXaxis()->SetTitleSize(28);
  hDataMix->GetXaxis()->SetTitleFont(43);
  hDataMix->GetXaxis()->SetTitleOffset(2.2);
  hDataMix->GetXaxis()->CenterTitle();
  
  
  hDataMix->GetYaxis()->SetLabelSize(22);
  hDataMix->GetYaxis()->SetLabelFont(43);
  hDataMix->GetYaxis()->SetTitleSize(28);
  hDataMix->GetYaxis()->SetTitleFont(43);
  hDataMix->GetYaxis()->SetTitleOffset(2.2);
  hDataMix->GetYaxis()->CenterTitle();
  
  if(drawXLabel) hDataMix->SetXTitle("A_{J} = (p_{T,1}-p_{T,2})/(p_{T,1}+p_{T,2})");  
  hDataMix->SetYTitle("Event Fraction");

  if(iPlot==1) hDataMix->SetXTitle("#eta_{dijet} = (#eta_{1}+#eta_{2})/2");
  if(iPlot==21) hDataMix->SetXTitle("#eta_{1}");
  if(iPlot==22) hDataMix->SetXTitle("#eta_{2}");

  if(iPlot==8) hDataMix->SetXTitle("p_{T,1}^{background}-p_{T,2}^{background}");
  if(iPlot==9) hDataMix->SetXTitle("N_{trk}^{offline}");

  if(iPlot==11) hDataMix->SetXTitle("p_{T,1}^{background}");
  if(iPlot==12) hDataMix->SetXTitle("p_{T,2}^{background}");

  if(iPlot==13) hDataMix->SetXTitle("p_{T} leading track in leading jet");
  if(iPlot==14) hDataMix->SetXTitle("p_{T} leading track in subleading jet");


  if(drawXLabel){
    if(iPlot==0) hDataMix->SetXTitle("p_{T,2}/p_{T,1}");
    if(iPlot==3) hDataMix->SetXTitle("#Delta #phi_{1,2}");
  }

  if(iPlot==0) hDataMix->SetMaximum(0.32);

  if(iPlot == 1 || iPlot==21 || iPlot==22){
    hDataMix->SetMaximum(0.3);
    hDataMix->SetMinimum(0);
  }


  if(iPlot==3){
    hDataMix->SetMaximum(2.52);
    hDataMix->SetMinimum(0.00004);
  }

  if(iPlot==8){
    hDataMix->SetMaximum(10);
    hDataMix->SetMinimum(0.00004);
  }

  if(iPlot==9){
    hDataMix->SetMaximum(0.15);
    hDataMix->SetMinimum(0.);
  }

  if(iPlot==11 || iPlot == 12){
    hDataMix->SetMaximum(10);
    hDataMix->SetMinimum(0.00004);
  }



  //hDataMix->GetXaxis()->SetNdivisions(905,true);
  hDataMix->GetYaxis()->SetNdivisions(505,true);

  hReference->SetMarkerColor(ppColor);
  hReference->SetLineColor(ppColor);
  hReference->SetMarkerStyle(25);

	hDataMix->Draw();//"hist");
	hDataMix->Draw("hist same");
	if(showPbPb) hReference->Draw("same");

	cout<<"PbPb ENTRIES : "<<endl;
	cout<<hReference->GetEntries()<<endl;

        cout<<"PbPb integral : "<<endl;
        cout<<hReference->Integral()<<endl;

	cout<<"pPb integral : "<<endl;
        cout<<h->Integral()<<endl;

	cout<<"Mix integral : "<<endl;
        cout<<hDataMix->Integral()<<endl;

	h->SetLineWidth(1);
	h->Draw("same");
	h->SetLineWidth(2);
	h->Draw("same");

	//	hDataMixB->Draw("same hist");
	//        hB->Draw("same hist");

	 cout<<" mean value of data "<<h->GetMean()<<endl;

  if(drawLeg){
    TLegend *t3=new TLegend(0.01,0.7,0.3,0.95);
    t3->AddEntry(h,"pPb #sqrt{s}=5.02 TeV","p");
    if(showPbPb) t3->AddEntry(hReference,"PbPb #sqrt{s}=2.76 TeV","p");
    t3->AddEntry(hDataMix,"PYTHIA+HIJING","lf");

    t3->SetFillColor(0);
    t3->SetBorderSize(0);
    t3->SetFillStyle(0);
    t3->SetTextFont(63);
    t3->SetTextSize(15);
    t3->Draw();
  }

}

void drawPatch(float x1, float y1, float x2, float y2){
  TLegend *t1=new TLegend(x1,y1,x2,y2);
  t1->SetFillColor(kWhite);
  t1->SetBorderSize(0);
  t1->SetFillStyle(1001);
  t1->Draw("");
}

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}





void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
      return;
   }
   canv->Clear();
   
   TPad* pad[columns][rows];

   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth = 
   (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
   (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
   (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}


