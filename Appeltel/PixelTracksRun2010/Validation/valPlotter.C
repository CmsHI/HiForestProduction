{

  const int ETABINS = 6;
  const int CBINS = 9;

  char etaStrings[ETABINS][256] = {"0.5","0.8","1","1.5","2", "2.4"};

  char centStrings[CBINS][256] = {"5-10%", "10-15%", "15-20%", "20-30%", "30-40%",
                                  "40-50%", "50-60%", "60-70%", "70-80%"};

  char directory[256] = "DQMData/Tracking/Track/hiMergedPtSplit12_tpRecoAssocMerged12s/";


  double ptMax = 5.;

  TPaveText *pt[6][9][20];
  TText *text[6][9][20];

  TLine *l[6][9][4];


  gStyle->SetOptStat(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
//  gStyle->SetFrameBorderMode(0);
//  gStyle->SetFrameLineColor(0);
  gStyle->SetTitleColor(0);
//  gStyle->SetTitleBorderSize(0); 
  gROOT->ForceStyle();


  TCanvas * c[ETABINS][CBINS];
  TCanvas * e[CBINS];
  TFile * f[ETABINS][CBINS];

  
  TH1F * efficPt[ETABINS][CBINS];
  TH1F * effic[CBINS];
  TH1F * momRes[ETABINS][CBINS];
  TH1F * fakePt[ETABINS][CBINS];
  TH1F * fake[CBINS];

  for( int cbin = 0; cbin<CBINS; cbin++ )
  {


     for( int eta = 0; eta<ETABINS; eta++ )
     {

       c[eta][cbin] = new TCanvas(Form("plots_%d_%d",eta, cbin),
                                  Form("Track Validation for eta = %s at %s Centrality ",
                                       etaStrings[eta],centStrings[cbin]),
                                       1600.,1200.);

       c[eta][cbin]->SetFillColor(0);
       c[eta][cbin]->SetBorderSize(2);
       c[eta][cbin]->SetFrameLineWidth(2);

       c[eta][cbin]->Divide(2,2);

       f[eta][cbin] = TFile::Open(Form("validation_eta%s_centbin%d.root",etaStrings[eta],cbin+1));
//       f[eta][cbin] = TFile::Open("validation10_20_30.root");
      
       c[eta][cbin]->cd(1);
       efficPt[eta][cbin] = (TH1F*) f[eta][cbin]->Get(Form("%sefficPt",directory));
       efficPt[eta][cbin]->SetTitle("");
       efficPt[eta][cbin]->GetYaxis()->SetRangeUser(0,1);
       efficPt[eta][cbin]->Draw();
       pt[eta][cbin][1] = new TPaveText(0.00862069,0.925678,0.9658046,0.9997627,"blNDC");
     //  pt[eta][cbin][1]->SetName(Form("title%d_%d_%d",eta,cbin,1));
       pt[eta][cbin][1]->SetBorderSize(2);
       text[eta][cbin][1] = pt[eta][cbin][1]->AddText("Efficiency vs p_{T}");
       pt[eta][cbin][1]->SetFillColor(0);
       pt[eta][cbin][1]->Draw();
       l[eta][cbin][1] = new TLine(0,0.5,ptMax,0.5);
       l[eta][cbin][1]->SetLineColor(kRed);
       l[eta][cbin][1]->Draw();

       c[eta][cbin]->cd(2);
       momRes[eta][cbin] = (TH1F*) f[eta][cbin]->Get(Form("%sptres_vs_pt_Sigma",directory));
       momRes[eta][cbin]->SetTitle("");
       momRes[eta][cbin]->GetYaxis()->SetRangeUser(0,0.1);
       momRes[eta][cbin]->Draw();
       pt[eta][cbin][3] = new TPaveText(0.00862069,0.925678,0.9658046,0.9997627,"blNDC");
       pt[eta][cbin][3]->SetBorderSize(2);
       text[eta][cbin][3] = pt[eta][cbin][3]->AddText("#sigma(#deltap_{T}/p_{T}) vs p_{T}");
       pt[eta][cbin][3]->SetFillColor(0);
       pt[eta][cbin][3]->Draw();

       
       c[eta][cbin]->cd(3);
       fakePt[eta][cbin] = (TH1F*) f[eta][cbin]->Get(Form("%sfakeratePt",directory));
       fakePt[eta][cbin]->SetTitle("");
       fakePt[eta][cbin]->GetYaxis()->SetRangeUser(0,0.5);
       fakePt[eta][cbin]->Draw();
       pt[eta][cbin][0] = new TPaveText(0.00862069,0.925678,0.9658046,0.9997627,"blNDC");
       pt[eta][cbin][0]->SetBorderSize(2);
       text[eta][cbin][0] = pt[eta][cbin][0]->AddText("Fake Rate vs p_{T}");
       pt[eta][cbin][0]->SetFillColor(0);
       pt[eta][cbin][0]->Draw();
       l[eta][cbin][0] = new TLine(0,0.05,ptMax,0.05);
       l[eta][cbin][0]->SetLineColor(kRed);
       l[eta][cbin][0]->Draw();

       c[eta][cbin]->cd(4);

       pt[eta][cbin][6] = new TPaveText(0.1,0.7,0.9,0.9,"blNDC");
       pt[eta][cbin][6]->SetBorderSize(0);
       text[eta][cbin][6] = pt[eta][cbin][6]->AddText(
           Form("P_{T} dependent variables"));
       pt[eta][cbin][6]->SetFillColor(0);
       pt[eta][cbin][6]->Draw();

       pt[eta][cbin][7] = new TPaveText(0.1,0.58,0.9,0.65,"blNDC");
       pt[eta][cbin][7]->SetBorderSize(0);
       text[eta][cbin][7] = pt[eta][cbin][7]->AddText(
           Form("|#eta|<%s",
                etaStrings[eta]));
       pt[eta][cbin][7]->SetFillColor(0);
       pt[eta][cbin][7]->Draw();

       pt[eta][cbin][8] = new TPaveText(0.1,0.3,0.9,0.45,"blNDC");
       pt[eta][cbin][8]->SetBorderSize(0);
       text[eta][cbin][8] = pt[eta][cbin][8]->AddText(
           Form("Centrality Range %s",
                centStrings[cbin]));
       pt[eta][cbin][8]->SetFillColor(0);
       pt[eta][cbin][8]->Draw();

       pt[eta][cbin][9] = new TPaveText(0.94,0.02,0.97,0.08,"blNDC");
       pt[eta][cbin][9]->SetBorderSize(0);
       text[eta][cbin][9] = pt[eta][cbin][9]->AddText(Form("%d",
           (cbin)*(ETABINS+1) + eta+4));
       pt[eta][cbin][9]->SetFillColor(0);
       pt[eta][cbin][9]->Draw();

       c[eta][cbin]->SaveAs(Form("page_%d.ps",(cbin)*(ETABINS+1) + eta+4));

       if( eta == ETABINS-1 )
       {

         e[cbin] = new TCanvas(Form("plots_%d", cbin),
                                  Form("Track Validation at %s Centrality ",
                                       centStrings[cbin]),
                                       1600.,1200.);
         e[cbin]->SetFillColor(0);
         e[cbin]->SetBorderSize(2);
         e[cbin]->SetFrameLineWidth(2);
         e[cbin]->Divide(2,2);

         e[cbin]->cd(2);
         effic[cbin] = (TH1F*) f[eta][cbin]->Get(Form("%seffic",directory));
         effic[cbin]->SetTitle("");
         effic[cbin]->GetYaxis()->SetRangeUser(0,1);
         effic[cbin]->Draw();
         pt[eta][cbin][2] = new TPaveText(0.00862069,0.925678,0.9658046,0.9997627,"blNDC");
         pt[eta][cbin][2]->SetBorderSize(2);
         text[eta][cbin][2] = pt[eta][cbin][2]->AddText("Efficiency vs #eta");
         pt[eta][cbin][2]->SetFillColor(0);
         pt[eta][cbin][2]->Draw();
         l[eta][cbin][2] = new TLine(-2.5,0.5,2.5,0.5);
         l[eta][cbin][2]->SetLineColor(kRed);
         l[eta][cbin][2]->Draw();

         e[cbin]->cd(4);
         fake[cbin] = (TH1F*) f[eta][cbin]->Get(Form("%sfakerate",directory));
         fake[cbin]->SetTitle("");
         fake[cbin]->GetYaxis()->SetRangeUser(0,0.5);
         fake[cbin]->Draw();
         pt[eta][cbin][4] = new TPaveText(0.00862069,0.925678,0.9658046,0.9997627,"blNDC");
         pt[eta][cbin][4]->SetBorderSize(2);
         text[eta][cbin][4] = pt[eta][cbin][4]->AddText("Fake Rate vs #eta");
         pt[eta][cbin][4]->SetFillColor(0);
         pt[eta][cbin][4]->Draw();
         l[eta][cbin][3] = new TLine(-2.5,0.05,2.5,0.05);
         l[eta][cbin][3]->SetLineColor(kRed);
         l[eta][cbin][3]->Draw();

         e[cbin]->cd(3);
         fake[cbin] = (TH1F*) f[eta][cbin]->Get(Form("%sptres_vs_eta_Sigma",directory));
         fake[cbin]->SetTitle("");
         fake[cbin]->GetYaxis()->SetRangeUser(0,0.1);
         fake[cbin]->Draw();
         pt[eta][cbin][5] = new TPaveText(0.00862069,0.925678,0.9658046,0.9997627,"blNDC");
         pt[eta][cbin][5]->SetBorderSize(2);
         text[eta][cbin][5] = pt[eta][cbin][5]->AddText("#sigma(#delta p_{T}/p_{T}) vs #eta");
         pt[eta][cbin][5]->SetFillColor(0);
         pt[eta][cbin][5]->Draw();

       e[cbin]->cd(1);

       pt[eta][cbin][10] = new TPaveText(0.1,0.7,0.9,0.9,"blNDC");
       pt[eta][cbin][10]->SetBorderSize(0);
       text[eta][cbin][10] = pt[eta][cbin][10]->AddText(
           Form("Tracking Validation"));
       pt[eta][cbin][10]->SetFillColor(0);
       pt[eta][cbin][10]->Draw();

       pt[eta][cbin][11] = new TPaveText(0.1,0.55,0.9,0.65,"blNDC");
       pt[eta][cbin][11]->SetBorderSize(0);
       text[eta][cbin][11] = pt[eta][cbin][11]->AddText(
       Form("Centrality Range %s",
                                       centStrings[cbin]));
       pt[eta][cbin][11]->SetFillColor(0);
       pt[eta][cbin][11]->Draw();

       e[cbin]->cd(4);

       pt[eta][cbin][12] = new TPaveText(0.94,0.02,0.97,0.08,"blNDC");
       pt[eta][cbin][12]->SetBorderSize(0);
       text[eta][cbin][12] = pt[eta][cbin][12]->AddText(Form("%d",
           (ETABINS+1)*cbin + 3));
       pt[eta][cbin][12]->SetFillColor(0);
       pt[eta][cbin][12]->Draw();


       e[cbin]->SaveAs(Form("page_%d.ps",(ETABINS+1)*cbin + 3));

       }


     }
   }

}



