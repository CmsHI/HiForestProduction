/*
 * =====================================================================================
 *
 *       Filename:  extract_pt.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/30/13 04:49:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

{
	TFile *f = new TFile("spectra.root");
	TH1D* hPt[40];
	for ( int i = 0; i < 40; i++ ) {
		cout << Form("spectra/hPt_%i",i) << endl;
		hPt[i] = (TH1D*)f->Get(Form("spectra/hPt_%i",i));
	}
	for ( int i = 0; i < 40; i++ ) {
		cout << hPt[i]->GetMean() << ", ";
		if ( !((i+1)%5) ) cout << endl;
	}
	cout << endl;
	for ( int i = 0; i < 40; i++ ) {
		cout << hPt[i]->GetRMS() * hPt[i]->GetRMS() << ", ";
		if ( !((i+1)%5) ) cout << endl;
	}
}
