
#include <iostream>
#include <utility>
#include <TROOT.h>
#include <TStyle.h>

#include "TFile.h"
#include "TCanvas.h"

#include "TH1F.h"
#include "TH1D.h"

#include "TH2F.h"
#include "TH2D.h"

#include "TF1.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TChain.h"
#include "TGraph.h"

#include "TLatex.h"

Int_t color[12]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,kRed+2,1,1};


void printCanvases(TCanvas * MyCanvas, const char * name, int log=0, int doCMS=2)
{
   MyCanvas->cd();
   MyCanvas->SetLogx(log);
   printf("canvas name: %s\n",name);

   // add some text labels                                                                                                                        
   double ndcX = 0.25;
   double ndcY = 0.9;
   TLatex * tex;

   if (doCMS>0) {
      /*
      tex = new TLatex(0.5,0.5,"CMS Preliminary");
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      tex->SetNDC();
      tex->Draw();
      */
   }


   //MyCanvas->Print(Form("./fig/%s.eps",name));
   MyCanvas->Print(Form("./fig/%s.gif",name));
   MyCanvas->Print(Form("./fig/%s.pdf",name));

   /*
   tex = new TLatex(ndcX,ndcY,"(a)");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   tex->Draw();

   MyCanvas->Print(Form("%s_A.eps",name));
   MyCanvas->Print(Form("%s_A.gif",name));
   MyCanvas->Print(Form("%s_A.pdf",name));

   if (tex) tex->Delete();
   tex = new TLatex(ndcX,ndcY,"(b)");
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   tex->Draw();

   MyCanvas->Print(Form("%s_B.eps",name));
   MyCanvas->Print(Form("%s_B.gif",name));
   MyCanvas->Print(Form("%s_B.pdf",name));
   */

}


void th1Style1( TH1 *h=0, int mcol = 1, int marker = 20, int msize = 1, int lcol = 1, int lsize = 1)
{
   h->SetMarkerStyle(marker);
   h->SetMarkerColor(mcol);
   h->SetMarkerSize(msize);

   h->SetLineColor(lcol);
   h->SetLineWidth(lsize);
}

void th1Style1( TGraph *g=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1)
{
   g->SetMarkerStyle(marker);
   if(mcol>95 || mcol<14){
      g->SetMarkerColor(mcol);
   }else{
      g->SetMarkerColor(color[mcol-14]);
   }
   g->SetMarkerSize(msize);
   if(lcol>95 || lcol<14){
      g->SetLineColor(lcol);
   }else{
      g->SetLineColor(color[lcol-14]);
   }
   g->SetLineWidth(lsize);
   g->SetLineStyle(lstyle);

   if(dStyle==1) g->Draw("PZsame");
   if(dStyle==2) g->Draw("Xlsame");
   if(dStyle==3) g->Draw("histsame");
}

void th1Style2( TGraph *g=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1, bool drawIt=false)
{
  g->SetMarkerStyle(marker);
  if(mcol>95 || mcol<14){
    g->SetMarkerColor(mcol);
  }else{
    g->SetMarkerColor(color[mcol-14]);
  }
  g->SetMarkerSize(msize);
  if(lcol>95 || lcol<14){
    g->SetLineColor(lcol);
  }else{
    g->SetLineColor(color[lcol-14]);
  }
  g->SetLineWidth(lsize);
  g->SetLineStyle(lstyle);

  if(drawIt){
    if(dStyle==1) g->Draw("PZsame");
    if(dStyle==2) g->Draw("Xlsame");
    if(dStyle==3) g->Draw("histsame");
  }
}

void th1Style2( TF1 *g=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1, bool drawIt=false)
{
  /*
  g->SetMarkerStyle(marker);
  if(mcol>95 || mcol<14){
    g->SetMarkerColor(mcol);
  }else{
    g->SetMarkerColor(color[mcol-14]);
  }
  g->SetMarkerSize(msize);
  */

  if(lcol>95 || lcol<14){
    g->SetLineColor(lcol);
  }else{
    g->SetLineColor(color[lcol-14]);
  }
  g->SetLineWidth(lsize);
  g->SetLineStyle(lstyle);

  if(drawIt){
    if(dStyle==1) g->Draw("same");
    if(dStyle==2) g->Draw("same");
    if(dStyle==3) g->Draw("same");
  }
}


void th1Style1( TH1 *h=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1)
{
   h->SetMarkerStyle(marker);
   if(mcol>95 || mcol<14){
      h->SetMarkerColor(mcol);
   }else{
      h->SetMarkerColor(color[mcol-14]);
   }
   h->SetMarkerSize(msize);
   if(lcol>95 || lcol<14){
      h->SetLineColor(lcol);
   }else{
      h->SetLineColor(color[lcol-14]);
   }
   h->SetLineWidth(lsize);
   h->SetLineStyle(lstyle);

   if(dStyle==1) h->Draw("PZsame");
   if(dStyle==2) h->Draw("Csame");
   if(dStyle==3) h->Draw("histsame");
}



void th1Style2( TH1 *h=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1, bool drawstat=false)
{
  h->SetMarkerStyle(marker);
  if(mcol>95 || mcol<14){
    h->SetMarkerColor(mcol);
  }else{
    h->SetMarkerColor(color[mcol-14]);
  }
  h->SetMarkerSize(msize);
  if(lcol>95 || lcol<14){
    h->SetLineColor(lcol);
  }else{
    h->SetLineColor(color[lcol-14]);
  }
  h->SetLineWidth(lsize);
  h->SetLineStyle(lstyle);

  if(drawstat){
    if(dStyle==1) h->Draw("PZsames");
    if(dStyle==2) h->Draw("Csames");
    if(dStyle==3) h->Draw("histsames");
  }else{
    if(dStyle==1) h->Draw("PZsame");
    if(dStyle==2) h->Draw("Csame");
    if(dStyle==3) h->Draw("histsame");
  }
}


