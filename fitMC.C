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
#include "TLatex.h"

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
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"


using namespace std;
using namespace RooFit;


/////////////////////////////////////////////////////////////
void fitMC(TString fileName = "toyMC.root") {

	style();
	gStyle->SetOptFit(1);
	gStyle->SetStatY(0.9);
	gROOT->ForceStyle();

	Double_t hist_minX = 5.0;
	Double_t hist_maxX = 15.0;

	TFile *file = new TFile(fileName, "READ");

	TH1D *histSB;
	TH1D *histBkg;

	histSB = (TH1D *)file->Get("hSB");
	histBkg = (TH1D *)file->Get("hBkg");

	RooRealVar mass("mass", "mass", 0.0, 15.0, "GeV/c^2");

	RooDataHist *rooHistSig = new RooDataHist("rooHistSig", "Signal+bkg. hist", RooArgList(mass), histSB);
	RooDataHist *rooHistBkg = new RooDataHist("rooHistBkg", "Background hist", RooArgList(mass), histBkg);



	TCanvas *cnv = new TCanvas("cnv", "cnv");

	// Create parameters, set initial values and ranges
	RooRealVar Bkg_A("Bkg_A","A", -0.5, -10.0, 10.0);
	RooRealVar Bkg_B("Bkg_B","B", -0.05, -10.0, 10.0);

	RooRealVar Sig_mu("Sig_mu","#mu", 8.0, -10.0, 20.0);
	RooRealVar Sig_sigma("Sig_sigma","#sigma", 0.3, -10.0, 20.0);
	RooRealVar Sig_mu_2("Sig_mu_2","#mu", 9.0, -10.0, 20.0);
	RooRealVar Sig_sigma_2("Sig_sigma_2","#sigma", 0.4, -10.0, 20.0);
	RooRealVar Sig_mu_3("Sig_mu_3","#mu", 12.0, -10.0, 20.0);
	RooRealVar Sig_sigma_3("Sig_sigma_3","#sigma", 0.2, -10.0, 20.0);

	RooRealVar Bkg_frac("Bkg_frac","Bkg_frac", 0.7, 0.0, 1.0);
	RooRealVar Bkg_frac2("Bkg_frac","Bkg_frac", 0.2, 0.0, 1.0);
	RooRealVar Sig_frac_1("Sig_frac_1","Sig_frac_1", 0.5, 0.0, 1.0);
	RooRealVar Sig_frac_2("Sig_frac_2","Sig_frac_2", 0.6, 0.0, 1.0);


	// Create PDFs
	RooExponential *exp = new RooExponential("exp","Exp background", mass, Bkg_A);
	RooExponential *exp2 = new RooExponential("exp2","Exp background 2", mass, Bkg_B);
	RooGaussian *gauss = new RooGaussian("gauss", "Gaussian signal", mass, Sig_mu, Sig_sigma);
	RooGaussian *gauss2 = new RooGaussian("gauss2", "Gaussian signal 2", mass, Sig_mu_2, Sig_sigma_2);
	RooGaussian *gauss3 = new RooGaussian("gauss3", "Gaussian signal 3", mass, Sig_mu_3, Sig_sigma_3);

	// Add PDFs together with recursive fractions (forces all ratios to be in [0,1])
	//RooAddPdf *total = new RooAddPdf("total", "signal+background", RooArgList(*exp, *exp2, *gauss, *gauss2, *gauss3), RooArgList(Bkg_frac, Bkg_frac2, Sig_frac_1, Sig_frac_2), kTRUE);
	RooAddPdf *total = new RooAddPdf("total", "signal+background", RooArgList(*exp, *exp2, *gauss, *gauss2, *gauss3), RooArgList(Bkg_frac, Bkg_frac2, Sig_frac_1, Sig_frac_2), kTRUE);

	RooFitResult *resultBkg = exp2->fitTo(*rooHistBkg, Save());

	RooAbsPdf *hessePDF = resultBkg->createHessePdf(RooArgSet(Bkg_B));

	RooFitResult *result = total->fitTo(*rooHistSig, Save(), ExternalConstraints(*hessePDF));


	// Plot simulated data
	cnv->cd();
	RooPlot* frame = mass.frame();
	frame->SetTitle("Data");
	rooHistBkg->plotOn(frame, MarkerColor(kRed));
	rooHistSig->plotOn(frame);
	total->plotOn(frame, VisualizeError(*result, 1));
	total->plotOn(frame);
	total->plotOn(frame, Components(*exp2), LineColor(kRed), LineStyle(kDashed), VisualizeError(*result, 1), FillColor(kRed), FillStyle(3001));
	total->plotOn(frame, Components(*exp2), LineColor(kRed), LineStyle(kDashed));
	total->plotOn(frame, Components(*gauss), LineColor(kTeal), LineStyle(kSolid), VisualizeError(*result, 1), FillColor(kTeal), FillStyle(3001));
	total->plotOn(frame, Components(*gauss), LineColor(kTeal), LineStyle(kSolid));
	total->plotOn(frame, Components(*gauss), LineColor(kTeal), LineStyle(kSolid), VisualizeError(*result, 1), FillColor(kTeal), FillStyle(3001));
	total->plotOn(frame, Components(RooArgSet(*gauss, *exp, *exp2)), LineColor(kOrange), LineStyle(kSolid));
	total->plotOn(frame);
	frame->Draw();
	//--------------------------------------------


	// Calculate chi2 after extracting the number of parameters - used to calculate Ndof
	Int_t nparam = total->getParameters(*rooHistSig)->selectByAttrib("Constant",kFALSE)->getSize();
	cout<<"nparam = "<<nparam<<endl;
	Double_t chi2ndf = frame->chiSquare(nparam);

	TLatex *latex = new TLatex();
	latex->SetTextSize(0.035);
	latex->SetNDC();
	latex->DrawLatex(0.25, 0.82, Form("#frac{#chi^{2}}{N_{dof}} = %.2f", chi2ndf));

	// pulls

	TCanvas *cnv_pull = new TCanvas("cnv_pull", "cnv_pull");
	cnv_pull->cd();


	RooHist* hpull;
	hpull = frame->pullHist("h_rooHistSig", "total_Norm[mass]");

	hpull->Draw();

	result->Print();
	result->correlationMatrix().Print();

	RooGenericPdf yield1("yield1", "yield from gauss fit 1", "(1.0-Bkg_frac)*(1-Bkg_frac2)*Sig_frac_1",RooArgSet(Bkg_frac, Bkg_frac2, Sig_frac_1));

	cout<<"yield1 = "<<yield1.getVal()<<endl;

	file->cd();

	histSB->Write("hSB");
	histBkg->Write("hBkg");

	cnv->cd()->SaveAs("toyMC.png");

	cout<<"Finished!"<<endl;


	//-----------------
}

