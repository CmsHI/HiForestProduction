


void fitSlices(TH2* hCorr, TF1* func = 0, bool variableBins = 0){

  bool useFit = func != 0;

  int nBins = hCorr->GetNbinsX();

  TH1D* hMean;
  TH1D* hSigma;

  if(variableBins){
    hMean = new TH1D(Form("%s_1",hCorr->GetName()),"",hCorr->GetNbinsX(),hCorr->GetXaxis()->GetXbins()->GetArray());
    hSigma = new TH1D(Form("%s_2",hCorr->GetName()),"",hCorr->GetNbinsX(),hCorr->GetXaxis()->GetXbins()->GetArray());
  }else{
    hMean = new TH1D(Form("%s_1",hCorr->GetName()),"",nBins,hCorr->GetXaxis()->GetXmin(),hCorr->GetXaxis()->GetXmax());
    hSigma = new TH1D(Form("%s_2",hCorr->GetName()),"",nBins,hCorr->GetXaxis()->GetXmin(),hCorr->GetXaxis()->GetXmax());
  }

  for(int i = 1; i < nBins+1; ++i){
    int bin = nBins - i;
    TH1D* h = hCorr->ProjectionY(Form("%s_bin%d",hCorr->GetName(),bin),i,i);

    if(h->GetEntries() < 5) continue;
    if(useFit){

      func->SetParameter(0,h->GetMaximum());
      func->SetParameter(1,h->GetMean());
      func->SetParameter(2,h->GetRMS());

      h->Fit(func);

      hMean->SetBinContent(i,func->GetParameter(1));
      hMean->SetBinError(i,func->GetParError(1));
      hSigma->SetBinContent(i,func->GetParameter(2));
      hSigma->SetBinError(i,func->GetParError(2));
    }else{

      hMean->SetBinContent(i,h->GetMean());
      hMean->SetBinError(i,h->GetMeanError());
      hSigma->SetBinContent(i,h->GetRMS());
      hSigma->SetBinError(i,h->GetRMSError());

    }

  }
}






