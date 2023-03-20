/*
 * generateToyMC.C
 *
 *  Created on: 14 lis 2016
 *      Author: Leszek Kosarzewski
 */

#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include "style.C"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooHist.h"


using namespace std;
using namespace RooFit;


/////////////////////////////////////////////////////////////
void generateToyMC(TString fileName = "toyMC.root") {

	style();
	gStyle->SetOptFit(1);
	gStyle->SetStatY(0.9);
	gROOT->ForceStyle();

	Double_t hist_minX = 5.0;
	Double_t hist_maxX = 15.0;

	TFile *file = new TFile(fileName, "RECREATE");

	TH1D *histSB;
	TH1D *histBkg;

	TCanvas *cnv = new TCanvas("cnv", "cnv");

	// Create parameters, set initial values and ranges
	RooRealVar Bkg_A("Bkg_A","A", -0.5);
	RooRealVar Bkg_B("Bkg_B","B", -0.05);

	RooRealVar Sig_mu("Sig_mu","#mu", 8.0);
	RooRealVar Sig_sigma("Sig_sigma","#sigma", 0.3);
	RooRealVar Sig_mu_2("Sig_mu_2","#mu", 9.0);
	RooRealVar Sig_sigma_2("Sig_sigma_2","#sigma", 0.4);
	RooRealVar Sig_mu_3("Sig_mu_3","#mu", 12.0);
	RooRealVar Sig_sigma_3("Sig_sigma_3","#sigma", 0.2);

	RooRealVar Bkg_frac("Bkg_frac","Bkg_frac", 0.7);
	RooRealVar Bkg_frac2("Bkg_frac","Bkg_frac", 0.2);
	RooRealVar Sig_frac_1("Sig_frac_1","Sig_frac_1", 0.5);
	RooRealVar Sig_frac_2("Sig_frac_2","Sig_frac_2", 0.6);

	// Create variable, set range and it`s unit
	RooRealVar mass("m_{ee}", "m_{ee}", hist_minX, hist_maxX, "GeV/c^{2}");
	mass.setBins(100); // set number of bins for histogram


	// Create PDFs
	RooExponential *exp = new RooExponential("exp","Exp background", mass, Bkg_A);
	RooExponential *exp2 = new RooExponential("exp2","Exp background 2", mass, Bkg_B);
	RooGaussian *gauss = new RooGaussian("gauss", "Gaussian signal", mass, Sig_mu, Sig_sigma);
	RooGaussian *gauss2 = new RooGaussian("gauss2", "Gaussian signal 2", mass, Sig_mu_2, Sig_sigma_2);
	RooGaussian *gauss3 = new RooGaussian("gauss3", "Gaussian signal 3", mass, Sig_mu_3, Sig_sigma_3);

	// Add PDFs together with recursive fractions (forces all ratios to be in [0,1])
	//RooAddPdf *total = new RooAddPdf("total", "signal+background", RooArgList(*exp, *exp2, *gauss, *gauss2, *gauss3), RooArgList(Bkg_frac, Bkg_frac2, Sig_frac_1, Sig_frac_2), kTRUE);
	RooAddPdf *total = new RooAddPdf("total", "signal+background", RooArgList(*exp, *gauss, *gauss2, *gauss3), RooArgList(Bkg_frac, Sig_frac_1, Sig_frac_2), kTRUE);

	// Generate binned MC data - 100k events
	RooDataHist *RooHistSig = total->generateBinned(mass, 100000);
	RooDataHist *RooHistBkg = exp2->generateBinned(mass, 40000);

	RooHistSig->add(*RooHistBkg);

	// Export TH1D with 100 bins
	histSB = (TH1D *)RooHistSig->createHistogram("histSB", mass, 100);
	histBkg = (TH1D *)RooHistBkg->createHistogram("histBkg", mass, 100);

	histSB->SetTitle("Signal; m_{ee} [GeV/c^2]; #frac{dN}{dm_{ee}} [#frac{c^2}{GeV}]");
	histBkg->SetTitle("Background exp2; m_{ee} [GeV/c^2]; #frac{dN}{dm_{ee}} [#frac{c^2}{GeV}]");

	// Plot simulated data
	cnv->cd();
	RooPlot* frame = mass.frame();
	frame->SetTitle("Data");
	RooHistSig->plotOn(frame);
	RooHistBkg->plotOn(frame, MarkerColor(kRed));
	frame->Draw();
	//--------------------------------------------

	file->cd();

	histSB->Write("hSB");
	histBkg->Write("hBkg");

	cnv->cd()->SaveAs("toyMC.png");

	cout<<"Finished!"<<endl;


	//-----------------
}

