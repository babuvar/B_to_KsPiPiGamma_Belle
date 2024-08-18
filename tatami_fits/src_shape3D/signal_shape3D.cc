#include "../src_aux/options.h"
#include <iostream>
#include <fstream>
#include <string>
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooSimPdfBuilder.h"
#include "RooTruthModel.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLegend.h"
#include "libRooTatami/RooTatamiHelper.h"
#include "libRooTatami/RooDtCPSignal.h"
#include "libRooTatami/RooDtLifetime.h"
#include <boost/lexical_cast.hpp>
#include "RooArgusBG.h"
#include "RooExponential.h"
#include "libRooTatami/RooDtBkg.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "TH2F.h"
#include "TStyle.h"


using namespace std;
using namespace RooFit;

bool fixedparameters = false; //Whether or not to fix parameters from Tristan's optimal values

// Main method
int main(int argc, char *argv[]) {

	//Deal with options
        Options *opt = new Options();
        if (!(opt->ParseOptions(argc, argv)))
        {
                std::exit(EXIT_FAILURE);
        } else
        {
                opt->PrintOptions();
        }

	//Get the input file (one file for now!)
	string filename = opt->GetFilenames().at(0);

	TChain * tree = new TChain("h1");
	tree->Add(filename.c_str());
	

	// Variable definitions
	RooRealVar mbc("mbc","m_{bc} [GeV]",5.20,5.30) ;
        RooRealVar de("de","#DeltaE [GeV]",-0.2,0.2) ;
        RooRealVar dt("dt","#DeltaT [ps]",-10.,10.) ;
        mbc.setRange("fitrange",5.20,5.30);
        de.setRange("fitrange",-0.2,0.2);
        dt.setRange("fitrange",-10.,10.);
        mbc.setRange("Signal",5.27,5.29);
        de.setRange("Signal",-0.1,0.05);
	RooCategory dalitzCategory("dalitzCategory","Split by charge") ;
                dalitzCategory.defineType("Plus",1) ;
                dalitzCategory.defineType("Minus",-1) ;
        RooRealVar tauB("\\tau_{B}", "tauB", 1.3, 2.0);
        RooRealVar w("w","mistag rate",0.,0.5) ;
        RooRealVar dw("dw","mistag rate",-0.5,0.5) ;
        RooCategory iflavor("iflavor","iflavor") ;
                iflavor.defineType("untagged",0);
                iflavor.defineType("positive",1);
                iflavor.defineType("negative",-1);
        RooRealVar A("A","A",-2.,2.) ;
        RooRealVar S("S","S",-2.,2.) ;

        RooCategory irbin("irbin", "irbin");
        for (int i = 0; i < 7; ++i){irbin.defineType(boost::lexical_cast<std::string > (i).c_str(), i);}

        RooCategory inexp("inexp", "inexp");
        inexp.defineType("7", 7);
        inexp.defineType("31", 31);	

	//Conditional variables
	RooRealVar CosThetaB("coscms", "CosThetaB", -1, 1);
        RooRealVar z_err_sig("vtxzerr", "z_err_sig", 0, 0.1);
        RooRealVar z_err_tag("tagzerr", "z_err_tag", 0, 0.05);
        RooRealVar chisq_tracks_sig("vtxchi2", "chisq_tracks_sig", 0, 600);
        RooRealVar chisq_tracks_tag("tagchi2", "chisq_tracks_tag", 0, 600);
        RooRealVar ecms("ecm", "ecms", 5.279, 5.4);
        RooRealVar dgf_tracks_sig("vtxndf", "dgf_tracks_sig", 0, 30);
        RooRealVar dgf_tracks_tag("tagndf", "dgf_tracks_tag", 0, 30);
        RooRealVar keeptagl_tag("taglepto", "keeptagl_tag", 0, 12);
        RooRealVar ntrk_sig("vtxntrk", "ntrk_sig", 0, 12);
        RooRealVar ntrk_tag("tagntrk", "ntrk_tag", 0, 12);
        RooRealVar issignal("issignal", "issignal", 0, 1);
        RooRealVar issgnevt("issgnevt", "issgnevt", 0, 1);

	const double delta_m = 0.507;
	const double tau_b = 1.5344;	

	// Declare conditional observables for fit and plot
	RooArgSet AllConditionalVariables(inexp, CosThetaB, ntrk_sig, ntrk_tag);
        AllConditionalVariables.add(z_err_sig);
        AllConditionalVariables.add(z_err_tag);
        AllConditionalVariables.add(chisq_tracks_sig);
        AllConditionalVariables.add(chisq_tracks_tag);
        AllConditionalVariables.add(dgf_tracks_sig);
        AllConditionalVariables.add(dgf_tracks_tag);
        AllConditionalVariables.add(keeptagl_tag);
        AllConditionalVariables.add(iflavor);
        AllConditionalVariables.add(w);
        AllConditionalVariables.add(dw);
        AllConditionalVariables.add(ecms);
        AllConditionalVariables.add(irbin);
        AllConditionalVariables.add(issignal);
        AllConditionalVariables.add(issgnevt);

        RooArgSet mainVariables(mbc, de, dt, dalitzCategory);
        mainVariables.add(AllConditionalVariables);

	// Import ttree
	RooDataSet * data = new RooDataSet("data", "data", mainVariables, Import(*tree));


//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//                                                                                              //
//  The signal is modelled as :                                                              //
//  Mbc : Crystalball                                                                           //
//  dE  :  Crystalball                                                                          //
//  dT  : Tatami resolution function                                                            //
//----------------------------------------------------------------------------------------------//

	//Continuum Mbc Model
	RooRealVar mbcmeanCB("\\mu_{sig}^{M_{bc}}","B^{#pm} mass", 5.27959,5.275,5.285);
        RooRealVar mbcsigmaCB("\\sigma_{sig}^{M_{bc}}","B^{#pm} width", 0.002881,5e-4,5e-3);
        RooRealVar mbcalphasig("\\alpha_{sig}^{M_{bc}}","B^{#pm} width", 1.16,0.,5.);
        RooRealVar mbcnsig("n_{sig}^{M_{bc}}","B^{#pm} width", 19,0,20);
        RooCBShape mbcCBall("signalMBC", "Crystal Ball shape", mbc, mbcmeanCB, mbcsigmaCB, mbcalphasig, mbcnsig);

	//Continuum dE Model
	RooRealVar demean("\\mu_{sig}^{\\Delta E}","B^{#pm} mass", -0.0067,-0.09,0.02) ;
        RooRealVar desigma("\\sigma_{sig}^{\\Delta E}","B^{#pm} width", 0.0366,0,0.05);
        RooRealVar dealpha("\\alpha_{sig}^{\\Delta E}","B^{#pm} width", 0.60, 0, 5);
        RooRealVar den("n_{sig}^{\\Delta E}","B^{#pm} width", 6.2,1, 50);
	RooCBShape deCBall("signalDE", "Crystal Ball shape", de, demean, desigma, dealpha, den);

	//Continuum dt Model
	RooDtCPSignal dtsignal("signalDT", "DtSignal", dt, S, A,
                        inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig,
                        dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag,
                        keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);


	//3D PDF
	RooProdPdf * signal3D =  new RooProdPdf("signal","signal",RooArgSet(mbcCBall, deCBall, dtsignal));

	//If parameters are fixed
	if(fixedparameters == true){
	//Mbc
	mbcmeanCB.setConstant(true);
	mbcsigmaCB.setConstant(true);
	mbcalphasig.setConstant(true);
	mbcnsig.setConstant(true);
	//dE
	demean.setConstant(true);
	desigma.setConstant(true);
	dealpha.setConstant(true);
	den.setConstant(true);
	}


//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//
	
	//Fitting
	RooFitResult * fit3D = signal3D->fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	

	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	data->plotOn(xframe_mbc);
	mbcCBall.plotOn(xframe_mbc, LineColor(kBlue), ProjWData(*data));
	//mbcCBall.paramOn(xframe_mbc, ProjWData(*data));	
	RooHist* pull_mbc = xframe_mbc->pullHist() ;
	RooPlot* f_pull_mbc = mbc.frame(Title("")) ;
	f_pull_mbc->addPlotable(pull_mbc,"BX0") ;
	f_pull_mbc->SetMaximum(4);
	f_pull_mbc->SetMinimum(-4);
	f_pull_mbc->GetYaxis()->SetLabelSize(0.11);
        f_pull_mbc->GetXaxis()->SetLabelSize(0.11);


	//dE
	RooPlot *xframe_de = de.frame(Bins(40));
	data->plotOn(xframe_de);
	deCBall.plotOn(xframe_de, LineColor(kBlue), ProjWData(*data));
	//deCBall.paramOn(xframe_de, ProjWData(*data));
	RooHist* pull_de = xframe_de->pullHist() ;
        RooPlot* f_pull_de = de.frame(Title("")) ;
        f_pull_de->addPlotable(pull_de,"BX0") ;
        f_pull_de->SetMaximum(4);
        f_pull_de->SetMinimum(-4);
	f_pull_de->GetYaxis()->SetLabelSize(0.11);
        f_pull_de->GetXaxis()->SetLabelSize(0.11);	

	//dt
	RooPlot *xframe_dt = dt.frame(Bins(40));
        data->plotOn(xframe_dt);
        dtsignal.plotOn(xframe_dt, LineColor(kBlue), ProjWData(*data));
        //dtsignal.paramOn(xframe_dt, ProjWData(*data));
        RooHist* pull_dt = xframe_dt->pullHist() ;
        RooPlot* f_pull_dt = dt.frame(Title("")) ;
        f_pull_dt->addPlotable(pull_dt,"BX0") ;
        f_pull_dt->SetMaximum(6);
        f_pull_dt->SetMinimum(-6);
        f_pull_dt->GetYaxis()->SetLabelSize(0.11);
        f_pull_dt->GetXaxis()->SetLabelSize(0.11);

	//Mbc-DE 2D	
	TH2F * mbc_de_datahist = (TH2F *)data->createHistogram("mbc,de", 60, 60);
        TH2F * mbc_de_pdfhist = (TH2F *)signal3D->createHistogram("mbc,de", 60, 60);
	mbc_de_pdfhist->Scale(mbc_de_datahist->GetEntries());
	TH2F * mbc_de_resid = new TH2F("residual","residual",60,5.2, 5.3, 60, -0.2, 0.2);
	mbc_de_resid->Add(mbc_de_datahist, 1.0);
	mbc_de_resid->Add(mbc_de_pdfhist, -1.0);
	TH2F * mbc_de_pull = new TH2F("pull","pull",60,5.2, 5.3, 60, -0.2, 0.2);
	for(int i = 0; i< 60; i++){
	for(int j= 0; j< 60; j++){
	float value = mbc_de_resid->GetBinContent(i,j);
	float uncert = sqrt(mbc_de_datahist->GetBinContent(i,j));
	mbc_de_pull->SetBinContent(i,j, (value/uncert) );
	}}



	//Canvases and Pads
	TCanvas* can = new TCanvas("c","c",1300,800) ; can->Divide(2,1);
	TPad *pad1 = new TPad("pad1", "pad 1",0.0,0.0,0.333,1.0,0);
	TPad *pad2 = new TPad("pad2", "pad 2",0.333,0.0,0.667,1.0,0);
	TPad *pad3 = new TPad("pad3", "pad 3",0.667,0.0,1.0,1.0,0);
	pad1->Draw(); 	pad2->Draw();   pad3->Draw();
	pad1->cd();
	TPad *pad11 = new TPad("pad11", "pad 11",0.0,0.2,1.0,1.0,0);
	TPad *pad12 = new TPad("pad12", "pad 12",0.0,0.0,1.0,0.2,0);
	pad11->Draw(); pad12->Draw();	
	pad2->cd();
	TPad *pad21 = new TPad("pad21", "pad 21",0.0,0.2,1.0,1.0,0);
        TPad *pad22 = new TPad("pad22", "pad 22",0.0,0.0,1.0,0.2,0);
	pad21->Draw(); pad22->Draw();
	pad3->cd();
	TPad *pad31 = new TPad("pad31", "pad 31",0.0,0.2,1.0,1.0,0);
        TPad *pad32 = new TPad("pad32", "pad 32",0.0,0.0,1.0,0.2,0);
        pad31->Draw(); pad32->Draw();

	//Draw
	pad11->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_mbc->Draw() ;
	pad12->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_mbc->Draw() ;
	pad21->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_de->Draw() ;
	pad22->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_de->Draw() ;
        pad31->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_dt->Draw() ;
        pad32->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_dt->Draw() ;
	can->SaveAs("plots/fix_shape/signal_shape3D.png");

	//2-D Plot
	gStyle->SetOptStat(0);
        TCanvas* can2 = new TCanvas("c2","c2",1300,800) ; can2->Divide(2,2);
        can2->cd(1); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); mbc_de_datahist->Draw("colz") ;
        can2->cd(2); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); mbc_de_pdfhist->Draw("colz") ;
	can2->cd(3); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); mbc_de_resid->Draw("colz") ;
	can2->cd(4); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); mbc_de_pull->Draw("colz") ;
        can2->SaveAs("plots/fix_shape/signal_shape3D_MbcDe.png");

	//Save fit results
	RooArgSet(fit3D->floatParsFinal()).writeToFile("results/fix_shape/signal_3D.txt");
	
	//Save fit_model in RooWorkspace root file
	RooWorkspace *rws = new RooWorkspace("signal", "signal");
	rws->import(mbcCBall); 
        rws->import(deCBall);
        rws->import(dtsignal);
	rws->writeToFile("pdf_models/fix_shape/signal_shape3D.root");

//----------------------------------------------------------------------------------------------//		

        return 0;
}
