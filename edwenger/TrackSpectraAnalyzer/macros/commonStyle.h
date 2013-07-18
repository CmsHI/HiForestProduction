#include "TLatex.h"

Int_t color[12]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,1,1,1};

void printCanvases(TCanvas * MyCanvas, const char * name, int log=0)
{
   MyCanvas->cd();
   MyCanvas->SetLogx(log);
   printf("canvas name: %s\n",name);

   MyCanvas->Print(Form("%s.C",name));
   MyCanvas->Print(Form("%s.gif",name));
   MyCanvas->Print(Form("%s.pdf",name));


}


void th1Style1( TH1 *h=0, int mcol = 1, int marker = 20, int msize = 1, int lcol = 1, int lsize = 1)
//void th1Style1( TH1 *h=0)
{
   h->SetMarkerStyle(marker);
   h->SetMarkerColor(mcol);
   h->SetMarkerSize(msize);

   h->SetLineColor(lcol);
   h->SetLineWidth(lsize);
}

/*
void th1Style1( TGraph *g=0, int mcol = 1, int marker = 20, int msize = 1, int lcol = 1, int lsize = 1, int dStyle=1)
//void th1Style1( TH1 *h=0)                                                                                                                
{
   g->SetMarkerStyle(marker);
   if(mcol>95 || mcol<5){
      g->SetMarkerColor(mcol);
   }else{
      g->SetMarkerColor(color[mcol-5]);
   }
   g->SetMarkerSize(msize);
   if(lcol>95 || lcol<5){
      g->SetLineColor(lcol);
   }else{
      g->SetLineColor(color[lcol-5]);
   }
   g->SetLineWidth(lsize);


   
   //g->GetXaxis()->SetRange(5,100); 

   if(dStyle==1) g->Draw("PZsame");
   if(dStyle==2) g->Draw("Xlsame");
   if(dStyle==3) g->Draw("histsame");
}
*/

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

void th1Style1( TH1 *g=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1)
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
   //if(dStyle==2) g->Draw("l");
   if(dStyle==3) g->Draw("histsame");
   if(dStyle==4) g->Draw("phistsame");

}





