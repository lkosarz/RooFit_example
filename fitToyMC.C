/*
 * fitToyMC.C
 *
 *  Created on: 14 lis 2016
 *      Author: Leszek Kosarzewski
 */

#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TAttText.h"
#include "TLatex.h"

#include "style.C"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooHist.h"
#include "RooConstVar.h"


using namespace std;
using namespace RooFit;


/////////////////////////////////////////////////////////////
void fitToyMC(TString fileName = "toyMC.root") {

	style();
	gStyle->SetOptFit(1);
	gStyle->SetStatY(0.9);
	gROOT->ForceStyle();

	Double_t hist_minX = 5.0;
	Double_t hist_maxX = 15.0;

	Double_t xclude_m_lo = 11.2;
	Double_t xclude_m_hi = 12.8;

	Double_t chi2_m_lo = 7.2;
	Double_t chi2_m_hi = 13.0;

	TFile *file = new TFile(fileName, "READ");


	TH1D * histSB = (TH1D *)file->Get("hSB");
	TH1D * histBkg = (TH1D *)file->Get("hBkg");

	TCanvas *cnv = new TCanvas("cnv", "cnv", 1000, 800);

	// Create parameters, set initial values and ranges
	RooRealVar Bkg_A("Bkg_A","A", -0.59, -10.0, 10.0);
	RooRealVar Bkg_B("Bkg_B","B", -0.059, -10.0, 10.0);
	//RooRealVar *Bkg_B = new RooRealVar("Bkg_B","B", -0.259, -10.0, 10.0);
	RooRealVar Bkg_Ba("Bkg_Ba","Ba", -0.25, -10.0, 10.0);

	RooRealVar Sig_mu("#mu","#mu", 8.0, 7.5, 8.5);
	RooRealVar Sig_sigma("#sigma","#sigma", 0.3, 0.01, 1.0);
	RooRealVar Sig_mu_2("#mu_{2}","#mu2", 9.0, 8.5, 9.5);
	RooRealVar Sig_sigma_2("#sigma_{2}","#sigma2", 0.4, 0.01, 1.0);
	RooRealVar Sig_mu_3("#mu_{3}","#mu3", 12.0, 11.0, 13.0);
	RooRealVar Sig_sigma_3("#sigma_{3}","#sigma3", 0.2, 0.01, 1.0);

	RooRealVar Bkg_frac("Bkg_frac","Bkg_frac", 0.9, 0.0, 1.0);
	RooRealVar Bkg_frac2("Bkg_frac2","Bkg_frac2", 0.2, 0.0, 1.0);
	RooRealVar Sig_frac_1("Sig_frac_1","Sig_frac_1", 0.5, 0.0, 1.0);
	RooRealVar Sig_frac_2("Sig_frac_2","Sig_frac_2", 0.5, 0.0, 1.0);

	// Create variable, set range and it`s unit
	RooRealVar mass("m_{ee}", "m_{ee}", hist_minX, hist_maxX, "GeV/c^{2}");

	// Setup ranges for fitting or plotting
	mass.setRange("lo", hist_minX, xclude_m_lo);
	mass.setRange("hi", xclude_m_hi, hist_maxX);
	mass.setRange("chi2", chi2_m_lo, chi2_m_hi);	// range in which to calculate chi2/ndf
	mass.setRange("full", hist_minX, hist_maxX);

	mass.setRange("gauss", 6.6, 9.4);
	mass.setRange("gauss2", 7.6, 10.4);
	mass.setRange("gauss3", 10.6, 13.4);

	//-------------------

	//Sig_mu.setConstant(kTRUE);
	//Sig_sigma.setConstant(kTRUE);
/*	Sig_mu_3.setConstant(kTRUE);
	Sig_sigma_3.setConstant(kTRUE);
	Sig_frac_2.setVal(1.0);
	Sig_frac_2.setConstant(kTRUE);
	Sig_frac_1.setVal(1.0);
	Sig_frac_1.setConstant(kTRUE);
*/
	//Bkg_frac.setVal(1.0);
	//Bkg_frac.setConstant(kTRUE);
	//Bkg_A.setConstant(kTRUE);

	//Sig_mu.setConstant(kTRUE);
	//Sig_mu_2.setConstant(kTRUE);

	//-------------------

	// Create PDFs
	RooExponential *exp = new RooExponential("exp","Exp background", mass, Bkg_A);
	RooExponential *exp2 = new RooExponential("exp2","Exp background 2", mass, Bkg_B);
	RooGaussian *gauss = new RooGaussian("gauss", "Gaussian signal", mass, Sig_mu, Sig_sigma);
	RooGaussian *gauss2 = new RooGaussian("gauss2", "Gaussian signal 2", mass, Sig_mu_2, Sig_sigma_2);
	RooGaussian *gauss3 = new RooGaussian("gauss3", "Gaussian signal 3", mass, Sig_mu_3, Sig_sigma_3);

	// Add PDFs together with recursive fractions (forces all ratios to be in [0,1])
	RooAddPdf *total = new RooAddPdf("total", "signal+background", RooArgList(*exp, *exp2, *gauss, *gauss2, *gauss3), RooArgList(Bkg_frac, Bkg_frac2, Sig_frac_1, Sig_frac_2), kTRUE);
	//RooAddPdf *total = new RooAddPdf("total", "signal+background", RooArgList(*exp), RooArgList(Bkg_frac), kTRUE);

	// Make RooDataHist from ROOT TH1D
	RooDataHist *rooHistSig = new RooDataHist("rooHistSig", "Signal+bkg. hist", RooArgList(mass), histSB);
	RooDataHist *rooHistBkg = new RooDataHist("rooHistBkg", "Background hist", RooArgList(mass), histBkg);


	// Perform a initial background only likelihood fit
	RooFitResult *resultBkg = exp2->fitTo(*rooHistBkg, Range("full"),Save());
	RooAbsPdf *paramPDF = resultBkg->createHessePdf(RooArgSet(Bkg_B));
	//RooExponential *exp2_tmp = new RooExponential(*exp2);
	//RooExponential *exp2_tmp = (RooExponential *)exp2->Clone("exp2_tmp");
	//Bkg_B.setConstant(kTRUE);
	Bkg_Ba.setVal(Bkg_B.getVal());
	Bkg_Ba.setError(Bkg_B.getError());
	//RooGaussian *constr = new RooGaussian("constr","constr",Bkg_B,RooConst(Bkg_B.getVal()),RooConst(Bkg_B.getError()));
	RooGaussian *constr = new RooGaussian("constr","constr",Bkg_B,RooConst(0.25),RooConst(0.05));
	RooExponential *exp2_tmp = new RooExponential("exp2_tmp","Exp background 2 tmp", mass, Bkg_Ba);

	//RooProdPdf *total = new RooProdPdf("total", "signal+background", RooArgSet(*totalc, *paramPDF));

/*
	TCanvas *cnv_b = new TCanvas("cnv_b", "cnv_b", 500, 600);
	cnv_b->cd();
	RooPlot* frame_b = Bkg_B.frame();
	frame_b->SetTitle("B");
	paramPDF->plotOn(frame_b);
	frame_b->Draw();*/

	// Perform a likelihood fit
	//RooFitResult *result = total->fitTo(*rooHistSig, Range("full"), Save());
	RooFitResult *result = total->fitTo(*rooHistSig, Range("full"), ExternalConstraints(*paramPDF), Save()); // use parameter PDF from initial fit as external constraint
	//RooFitResult *result = total->fitTo(*rooHistSig, Range("full"), ExternalConstraints(*constr), Save()); // use parameter PDF from initial fit as external constraint
	//RooFitResult *result = total->fitTo(*rooHistSig, Range("full"), Constrain(Bkg_B), Save()); // use parameter PDF from initial fit as external constraint
	//RooFitResult *result = total->chi2FitTo(*rooHistSig, Range("lo,hi"),Save());

	TCanvas *cnv_bkg = new TCanvas("cnv_bkg", "cnv_bkg", 1000, 800);
	cnv_bkg->cd();
	RooPlot* frame_bkg = mass.frame();
	frame_bkg->SetTitle("Bkg data");
	rooHistBkg->plotOn(frame_bkg, MarkerColor(kRed));
	exp2_tmp->plotOn(frame_bkg, NormRange("full"), Range("full"), VisualizeError(*resultBkg, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kRed));
	exp2_tmp->plotOn(frame_bkg, NormRange("full"), Range("full"), LineColor(kRed), LineStyle(kDashed));
	frame_bkg->Draw();

	// kTRUE = use linear error propagation WARNING! fails in presence of strong correlations
	Bool_t linearErrors = kFALSE; // use MC sampling method

	// Plot components of the fit (blue line - total, red line - range used for chi2 calculation)
	cnv->cd();
	RooPlot* frame = mass.frame();
	frame->SetTitle("Data");
	rooHistBkg->plotOn(frame, MarkerColor(kRed));
	exp2_tmp->plotOn(frame, NormRange("full"), Range("full"), VisualizeError(*resultBkg, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kRed));
	exp2_tmp->plotOn(frame, NormRange("full"), Range("full"), LineColor(kRed), LineStyle(kDashed)); // draw red dashed line to show the initial fit (hidden beneath uncertainties for the gray line)
	rooHistSig->plotOn(frame);
	total->plotOn(frame, Components(RooArgSet(*exp, *exp2)), NormRange("full"), Range("full"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kBlue+3));
	total->plotOn(frame, Components(RooArgSet(*exp, *exp2)), NormRange("full"), Range("full"), LineColor(kBlue+3), LineStyle(kSolid));
	total->plotOn(frame, NormRange("full"), Range("full"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kBlue));
	total->plotOn(frame, NormRange("full"), Range("full"));
	//total->plotOn(frame, NormRange("lo,hi"), Range("lo,hi"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kBlue), LineColor(kBlue), LineStyle(kSolid));
	//total->plotOn(frame, NormRange("lo,hi"), Range("lo,hi"));
	//total->plotOn(frame, NormRange("lo,hi"), Range("lo,hi"));
	total->plotOn(frame, Components(*exp), NormRange("full"), Range("full"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kGreen));
	total->plotOn(frame, Components(*exp), NormRange("full"), Range("full"), LineColor(kGreen), LineStyle(kSolid));
	total->plotOn(frame, Components(*exp2), NormRange("full"), Range("full"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kGray));
	total->plotOn(frame, Components(*exp2), NormRange("full"), Range("full"), LineColor(kGray), LineStyle(kSolid));
	total->plotOn(frame, Components(*gauss), NormRange("full"), Range("gauss"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kTeal));
	total->plotOn(frame, Components(*gauss), NormRange("full"), Range("gauss"), LineColor(kTeal), LineStyle(kSolid));
	total->plotOn(frame, Components(*gauss2), NormRange("full"), Range("gauss2"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kOrange));
	total->plotOn(frame, Components(*gauss2), NormRange("full"), Range("gauss2"), LineColor(kOrange), LineStyle(kSolid));
	total->plotOn(frame, Components(*gauss3), NormRange("full"), Range("gauss3"), VisualizeError(*result, 1), DrawOption("LF"), LineWidth(2), FillStyle(3001), FillColor(kMagenta));
	total->plotOn(frame, Components(*gauss3), NormRange("full"), Range("gauss3"), LineColor(kMagenta), LineStyle(kSolid));
	total->plotOn(frame, LineColor(kRed), NormRange("full"), Range("chi2"));	// Plot it last for chi2/ndf calculation
	total->paramOn(frame,Layout(0.68, 1.0, 1.0), ShowConstants(kTRUE), Format("NE", AutoPrecision(2)));

	// Calculate chi2 after extracting the number of parameters - used to calculate Ndof
	Int_t nparam = total->getParameters(*rooHistSig)->selectByAttrib("Constant",kFALSE)->getSize();
	cout<<"nparam = "<<nparam<<endl;
	Double_t chi2ndf = frame->chiSquare(nparam);

	frame->getAttText()->SetTextSize(0.03);

	frame->Draw();

	TLatex *latex = new TLatex();
	latex->SetTextSize(0.035);
	latex->SetNDC();
	latex->DrawLatex(0.25, 0.82, Form("#frac{#chi^{2}}{N_{dof}} = %.2f", chi2ndf));


	// Kolmogorov-Smirnov test
	TH1* h_func = total->createHistogram("h_func", mass);
	TH1* h_data = rooHistSig->createHistogram("h_data", mass);

	//h_func->Scale(1.0*binWidth);

	h_data->SetMarkerStyle(21);
	h_data->SetMarkerColor(kRed);

	//h_data->Draw("same");
	//h_func->Draw("same");

	double KStest = h_data->KolmogorovTest(h_func);
	latex->DrawLatex(0.25, 0.75, Form("K-S test = %.2f", KStest)); // K-S test = 1 means very high probability of data coming from the distribution described by the model

	//--------

	// draw Pulls
	TCanvas* cnv_pull = new TCanvas("cnv_pull","cnv_pull",800,300);// 800 400
	TCanvas* cnv_pull_bkg = new TCanvas("cnv_pull_bkg","cnv_pull_bkg",800,300);// 800 400

	RooHist* hpull;
	RooHist* hpull_bkg;
	hpull = frame->pullHist("h_rooHistSig", "total_Norm[m_{ee}]_Range[full]_NormRange[full]");
	hpull_bkg = frame->pullHist("h_rooHistBkg", "total_Norm[m_{ee}]_Comp[exp2]_Range[full]_NormRange[full]");

	hpull->SetTitle("Pull signal+bkg;  m_{ee} [GeV/c^{2}]; Pull signal+bkg");
	hpull_bkg->SetTitle("Pull bkg;  m_{ee} [GeV/c^{2}]; Pull bkg");

	hpull->SetLineColor(kBlack);
	hpull->SetMarkerColor(kBlack);
	hpull_bkg->SetLineColor(kRed);
	hpull_bkg->SetMarkerColor(kRed);

	// fit lines
	TF1 *line_pull = new TF1("line_pull", "[0]*x+[1]", hist_minX, hist_maxX);
	TF1 *line_pull_bkg = new TF1("line_pull_bkg", "[0]*x+[1]", hist_minX, hist_maxX);

	hpull->Fit(line_pull, "0", "",  hist_minX, hist_maxX);
	hpull_bkg->Fit(line_pull_bkg, "0", "",  hist_minX, hist_maxX);

	cnv_pull->cd();
	hpull->Draw();
	line_pull->Draw("same");

	cnv_pull_bkg->cd();
	hpull_bkg->Draw();
	line_pull_bkg->Draw("same");

	// Print results
	result->Print();
	result->correlationMatrix().Print();


	// Fit is always normalized to the counts in the fit range
	Double_t yield_hist = rooHistSig->sum(kFALSE);
	RooRealVar yieldHist("yieldHist","yieldHist", 0.0, 0.0, 1e9);
	//yieldHist.setConstant(kTRUE);
	yieldHist.setVal(yield_hist);

	// Calculate ratios from fit with uncertainties including correlations
	RooGenericPdf yieldbkg("yieldbkg", "yield from background fit", "Bkg_frac",RooArgSet(Bkg_frac));
	RooGenericPdf yieldbkg2("yieldbkg2", "yield from background fit", "(1.0-Bkg_frac)*Bkg_frac2",RooArgSet(Bkg_frac, Bkg_frac2));
	RooGenericPdf yield1("yield1", "yield from gauss fit 1", "(1.0-Bkg_frac)*(1-Bkg_frac2)*Sig_frac_1",RooArgSet(Bkg_frac, Bkg_frac2, Sig_frac_1));
	RooGenericPdf yield2("yield2", "yield from gauss fit 2", "(1.0-Bkg_frac)*(1-Bkg_frac2)*(1.0-Sig_frac_1)*Sig_frac_2",RooArgSet(Bkg_frac, Bkg_frac2, Sig_frac_1, Sig_frac_2));
	RooGenericPdf yield3("yield3", "yield from gauss fit 3", "(1.0-Bkg_frac)*(1-Bkg_frac2)*(1.0-Sig_frac_1)*(1.0-Sig_frac_2)",RooArgSet(Bkg_frac, Bkg_frac2, Sig_frac_1, Sig_frac_2));

	// Propagate uncertainties from the covariance matrix from fit
	Double_t yieldbkg_err = yield_hist*yieldbkg.getPropagatedError(*result);
	Double_t yieldbkg2_err = yield_hist*yieldbkg2.getPropagatedError(*result);
	Double_t yield1_err = yield_hist*yield1.getPropagatedError(*result);
	Double_t yield2_err = yield_hist*yield2.getPropagatedError(*result);
	Double_t yield3_err = yield_hist*yield3.getPropagatedError(*result);

	// Since we are using standard (not extended) likelihood, the total yield (sig+bkg) is equal to total histogram counts
	cout<<endl;
	cout<<"Yield from hist = "<<yield_hist<<endl;
	cout<<"Yields form fit:"<<endl;
	cout<<"background = "<<yield_hist*yieldbkg.getVal()<<"+/-"<<yieldbkg_err<<endl;
	cout<<"background 2 = "<<yield_hist*yieldbkg2.getVal()<<"+/-"<<yieldbkg2_err<<endl;
	cout<<"yield1 = "<<yield_hist*yield1.getVal()<<"+/-"<<yield1_err<<endl;
	cout<<"yield2 = "<<yield_hist*yield2.getVal()<<"+/-"<<yield2_err<<endl;
	cout<<"yield3 = "<<yield_hist*yield3.getVal()<<"+/-"<<yield3_err<<endl;

	RooGenericPdf ratio("ratio", "ratio of gauss 2 to 1", "yield2/yield1",RooArgSet(yield1, yield2));

	Double_t ratio_err = yieldbkg.getPropagatedError(*result);

	cout<<endl;
	cout<<"Ratios:"<<nparam<<endl;
	cout<<"Gauss2(teal)/Gauss1(orange) = "<<ratio.getVal()<<"+/-"<<ratio_err<<endl;
	cout<<endl;


	// yields from bin counting:


	Double_t dataTotal_gauss = rooHistSig->sumEntries("1","gauss");

	RooAbsReal *integ_bkg = exp->createIntegral(mass,NormSet(mass),Range("gauss"));
	RooAbsReal *integ_bkg2 = exp2->createIntegral(mass,NormSet(mass),Range("gauss"));
	RooAbsReal *integ_gauss2 = gauss2->createIntegral(mass,NormSet(mass),Range("gauss2")); // gauss2 overlaps wit gauss

	Double_t bkgFrac = integ_bkg->getVal();
	Double_t bkgFrac2 = integ_bkg2->getVal();
	Double_t gaussFrac2 = integ_gauss2->getVal();

	Double_t numBkg = yield_hist*yieldbkg.getVal()*bkgFrac;
	Double_t numBkg2 = yield_hist*yieldbkg2.getVal()*bkgFrac2;
	Double_t numGauss2 = yield_hist*yield2.getVal()*gaussFrac2; // gauss2 contribution from fit to gauss bin counts

	Double_t counts_gauss = dataTotal_gauss - numBkg - numBkg2 - numGauss2;
	Double_t counts_gauss_err = counts_gauss*yield1_err/(yield_hist*yield1.getVal()); // error calculation on bin counts


	cout<<endl;
	cout<<"Yields form bin counts:"<<endl;
	cout<<"yield1 = "<<counts_gauss<<"+/-"<<counts_gauss_err<<endl;

	//--------------------------------------------

	cnv->cd()->SaveAs("fitMC.png");
	cnv_pull->cd()->SaveAs("fitMC_pull.png");
	cnv_pull_bkg->cd()->SaveAs("fitMC_pullBkg.png");

	cout<<"Finished!"<<endl;


	//-----------------
}

