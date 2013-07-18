{

   TFile * f[3];
   TString htit[3];

   bool isMC = false;


   if( !isMC)
   {
     f[0] = new TFile("vertexAnalysis_prompt202792.root");
     f[1] = new TFile("vertexAnalysis_prompt210759.root");
     f[2] = new TFile("vertexAnalysis_prompt210759HighMult190.root");
     htit[0] = "Run 202792, SingleTrack";
     htit[1] = "Run 210759, SingleTrack";
     htit[2] = "Run 210759, HighMult190";

   }
   

   if( isMC)
   {
     f[0] = new TFile("vertexAnalysis_hijingMB.root");
     f[1] = new TFile("vertexAnalysis_hijing200trk.root");
     f[2] = new TFile("vertexAnalysis_hijingPU2.root");
     htit[0] = "HIJING Min-Bias";
     htit[1] = "HIJING 200trk trig.";
     htit[2] = "HIJING+HIJING Pileup";
   }

   bool doLog = true;
 
   // draw G-plus cut
 
   double gplusX[17] = {0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,2.0,3.0,5.0};
   double gplusY[17] = {30.,26.,21.,18.,16.,14.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.};

   TLine * gplusH[16];
   TLine * gplusV[16];

   TLine * veto = new TLine(0.05,0.,0.05,10.);
   veto->SetLineColor(kRed);
   veto->SetLineWidth(3);

   for( unsigned int i =0; i<16; i++)
   {
     gplusH[i] = new TLine( gplusX[i], gplusY[i], gplusX[i+1], gplusY[i]);
     gplusV[i] = new TLine( gplusX[i+1], gplusY[i], gplusX[i+1], gplusY[i+1]);
     gplusH[i]->SetLineColor(kRed);
     gplusH[i]->SetLineWidth(3);
     gplusV[i]->SetLineColor(kRed);
     gplusV[i]->SetLineWidth(3);


   }



   TCanvas c1("c1","Dz vs Ntrk",1100,500);
   TCanvas c2("c2","Dxy vs Dz",1100,500);

   c1.Divide(3,1);
   c2.Divide(3,1);


   TH1F * dzntrk[3];
   TH1F * dxydz[3];

   for(unsigned int i=0;i<3;i++)
   {
      dzntrk[i] = (TH1F *) f[i]->Get("vtxAna/assocVtxDzNtrk");
      dzntrk[i]->GetXaxis()->SetRangeUser(0,5);
      dzntrk[i]->SetTitle(htit[i]);
      dzntrk[i]->GetYaxis()->SetRangeUser(0,50);
      c1.cd(i+1);
      if( doLog) c1.GetPad(i+1)->SetLogz(1);
      dzntrk[i]->Draw("colz");
      for( int j=0;j<15;j++) gplusV[j]->Draw();
      for( int j=0;j<16;j++) gplusH[j]->Draw();

      dxydz[i] = (TH1F *) f[i]->Get("vtxAna/assocVtxDxyDz");
      dxydz[i]->SetTitle(htit[i]);
      dxydz[i]->GetXaxis()->SetRangeUser(0,0.5);
      dxydz[i]->GetYaxis()->SetRangeUser(0,10.);
      c2.cd(i+1);
      if( doLog) c2.GetPad(i+1)->SetLogz(1);
      dxydz[i]->Draw("colz");
      veto->Draw();

   }

} 
