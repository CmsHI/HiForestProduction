// later other measurement can be included


TGraphErrors* CMS_7TEV(Double_t scale){

   //7 TeV measurement at CMS                                                                                                                                                                                 
   TGraphErrors *gre = new TGraphErrors(33);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,0.15,9.94559);
   gre->SetPointError(0,0,0.4460492);
   gre->SetPoint(1,0.25,6.4274);
   gre->SetPointError(1,0,0.2867568);
   gre->SetPoint(2,0.35,3.69525);
   gre->SetPointError(2,0,0.1654095);
   gre->SetPoint(3,0.45,2.2442);
   gre->SetPointError(3,0,0.1009501);
   gre->SetPoint(4,0.55,1.40735);
   gre->SetPointError(4,0,0.06369867);
   gre->SetPoint(5,0.65,0.920301);
   gre->SetPointError(5,0,0.04193949);
   gre->SetPoint(6,0.75,0.616821);
   gre->SetPointError(6,0,0.0283265);
   gre->SetPoint(7,0.85,0.426732);
   gre->SetPointError(7,0,0.01975661);
   gre->SetPoint(8,0.95,0.299085);
   gre->SetPointError(8,0,0.01397616);
   gre->SetPoint(9,1.1,0.183225);
   gre->SetPointError(9,0,0.00845129);
   gre->SetPoint(10,1.3,0.0995915);
   gre->SetPointError(10,0,0.004670335);
   gre->SetPoint(11,1.5,0.056565);
   gre->SetPointError(11,0,0.002704398);
   gre->SetPoint(12,1.7,0.0337805);
   gre->SetPointError(12,0,0.001650384);
   gre->SetPoint(13,1.9,0.0208494);
   gre->SetPointError(13,0,0.001043639);
   gre->SetPoint(14,2.1,0.01323271);
   gre->SetPointError(14,0,0.0006804352);
   gre->SetPoint(15,2.3,0.00852133);
   gre->SetPointError(15,0,0.0004520424);
   gre->SetPoint(16,2.5,0.00555772);
   gre->SetPointError(16,0,0.0003051332);
   gre->SetPoint(17,2.7,0.00369416);
   gre->SetPointError(17,0,0.0002105437);
   gre->SetPoint(18,2.9,0.00275361);
   gre->SetPointError(18,0,0.0001617788);
   gre->SetPoint(19,3.1,0.001784498);
   gre->SetPointError(19,0,0.0001103706);
   gre->SetPoint(20,3.3,0.001251193);
   gre->SetPointError(20,0,8.105581e-05);
   gre->SetPoint(21,3.5,0.000932109);
   gre->SetPointError(21,0,6.296976e-05);
   gre->SetPoint(22,3.7,0.00072809);
   gre->SetPointError(22,0,5.106009e-05);
   gre->SetPoint(23,3.9,0.000483879);
   gre->SetPointError(23,0,3.643427e-05);
   gre->SetPoint(24,4.1,0.000397255);
   gre->SetPointError(24,0,3.101453e-05);
   gre->SetPoint(25,4.3,0.000281682);
   gre->SetPointError(25,0,2.426278e-05);
   gre->SetPoint(26,4.5,0.0001898557);
   gre->SetPointError(26,0,1.777722e-05);
   gre->SetPoint(27,4.7,0.0001404168);
   gre->SetPointError(27,0,1.377145e-05);
   gre->SetPoint(28,4.9,0.0001152827);
   gre->SetPointError(28,0,1.228788e-05);
   gre->SetPoint(29,5.1,8.94332e-05);
   gre->SetPointError(29,0,1.021356e-05);
   gre->SetPoint(30,5.3,6.56232e-05);
   gre->SetPointError(30,0,8.759586e-06);
   gre->SetPoint(31,5.5,5.48406e-05);
   gre->SetPointError(31,0,8.112416e-06);
   gre->SetPoint(32,5.7,4.50414e-05);
   gre->SetPointError(32,0,7.518979e-06);

   /*
   gre->SetPoint(33,5.9,3.60633e-05);
   gre->SetPointError(33,0,5.641449e-06);
   gre->SetPoint(34,6.25,2.018451e-05);
   gre->SetPointError(34,0,3.081411e-06);
   gre->SetPoint(35,6.75,1.966985e-05);
   gre->SetPointError(35,0,2.63818e-06);
   gre->SetPoint(36,7.25,1.572381e-05);
   gre->SetPointError(36,0,2.668063e-06);
   gre->SetPoint(37,7.75,8.84664e-06);
   gre->SetPointError(37,0,2.136684e-06);
   gre->SetPoint(38,8.25,1.264757e-05);
   gre->SetPointError(38,0,3.756436e-06);
   gre->SetPoint(39,8.75,6.31772e-06);
   gre->SetPointError(39,0,2.762383e-06);
   gre->SetPoint(40,9.25,5.99949e-06);
   gre->SetPointError(40,0,3.246975e-06);
   gre->SetPoint(41,9.25,5.99949e-06);
   gre->SetPointError(41,0,3.239439e-06);
   /*
   

   /*
   Double_t x[33], y[33];
   Double_t xerr[33], yerr[33];

   for(Int_t i=0; i<33;i++){
      gre->GetPoint(i,x[i],y[i]);
      gre->SetPoint(i,x[i],scale*y[i]);
      //gre->GetPointError(i,x[i],yerr[i]);
      xerr[i] = 0, yerr[i] = 0;
      gre->SetPointError(i,scale*xerr[i],scale*yerr[i]);
   }
   */
   return gre;
}

