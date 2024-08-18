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
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "libRooTatami/RooDtBkg.h"
#include "RooWorkspace.h"


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
//  The BBbar is modelled as :                                                              //
//  Mbc : Gaussian + Argus                                                                      //
//  dE  : Exponential                                                                           //
//  dT  : Tatami resolution function                                                            //
//----------------------------------------------------------------------------------------------//

	//BBbar Mbc Model
	//Argus
	RooRealVar argparBBbkg_c("argusBBbkg_c","argus shape parameter",-5,-140,2) ;
	RooRealVar argparBBbkg_m0("argusBBbkg_m0","argus shape parameter",5.29, 5.28, 5.30) ;
	RooArgusBG argusBBbkg("argusBBbkg","Argus PDF",mbc,argparBBbkg_m0,argparBBbkg_c);
	//Gaussian
	RooRealVar mbcmeanBBbkg("\\mu_{BBbkg}^{M_{bc}}","B^{#pm} mass",5.2797, 5.270,5.290) ;
        RooRealVar mbcsigmaBBbkg("\\sigma_{BBbkg}^{M_{bc}}","B^{#pm} width",0.0020,0.001,0.05) ;
        RooGaussian gaussBBbkg("gaussft","gaussian PDF",mbc,mbcmeanBBbkg,mbcsigmaBBbkg) ;
	//combined
	RooRealVar fmbcbbbkg("f_{BBbkg}^{M_{bc}}","#signal events",0.015, 0.,0.9) ;
	RooAddPdf sumMbcBBbkg("bbMBC","g+a2",RooArgList(gaussBBbkg,argusBBbkg),RooArgList(fmbcbbbkg)) ;

	//BBbar dE Model
	RooRealVar dealphaBBbkg("\\alpha_{BBbkg}^{\\Delta E}","dealphaBBbkg", -6.04, -10., 0.) ;
	RooExponential expoDeBBbkg("bbDE","background p.d.f.",de,dealphaBBbkg) ;

	//BBbar dt Model
	RooRealVar tau_BBbkg("\\tau_{BBbkg}", "tau_bkgBB", 1.4, 0.5, 1.8);
        RooRealVar mu_lgBB("\\mu_{l}^{BBbkg}", "mu_lBB", -0.62, -0.7, 0.5);
        RooRealVar mu_dgBB("\\mu_{d}^{BBbkg}", "mu_dBB", -0.04, -0.1, 0.1);
        RooRealVar f_deltgBB("f_{delt}^{BBbkg}", "f_deltBB", 0.15, 0.1,1);
        RooRealVar f_tailgBB("f_{tail}^{BBbkg}", "f_tailBB", 0.218, 0, 1);
        RooRealVar S_maingBB("S_{main}^{BBbkg}", "S_mainBB", 2.79, 0.4, 10);
        RooRealVar S_tailgBB("S_{tail}^{BBbkg}", "S_tailBB", 3.43, 1, 40);

	RooDtBkg DtBkgBB("bbDT", "DtBkgBB", dt,
                tau_BBbkg, mu_lgBB, mu_dgBB, f_deltgBB, f_tailgBB, S_maingBB, S_tailgBB,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);

	//3d PDF
	RooProdPdf * bb3D = new RooProdPdf("bb3D","bb3D",RooArgList(sumMbcBBbkg, expoDeBBbkg, DtBkgBB));

	//If parameters are fixed
	if(fixedparameters == true){
	//Mbc
	argparBBbkg_c.setConstant(true);
	argparBBbkg_m0.setConstant(true);
	mbcmeanBBbkg.setConstant(true);
	mbcsigmaBBbkg.setConstant(true);
	fmbcbbbkg.setConstant(true);
	//dE
	dealphaBBbkg.setConstant(true);
	//dt
	tau_BBbkg.setConstant(true);
	mu_lgBB.setConstant(true);
	mu_dgBB.setConstant(true);
	f_deltgBB.setConstant(true);
	f_tailgBB.setConstant(true);
	S_maingBB.setConstant(true);
	S_tailgBB.setConstant(true);
	}


//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//
	
	//Fitting
	RooFitResult * fit3D = bb3D->fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	

	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	data->plotOn(xframe_mbc);
	sumMbcBBbkg.plotOn(xframe_mbc, LineColor(kBlue), ProjWData(*data));
	//sumMbcBBbkg.paramOn(xframe_mbc, ProjWData(*data));	
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
	expoDeBBbkg.plotOn(xframe_de, LineColor(kBlue), ProjWData(*data));
	//expoDeBBbkg.paramOn(xframe_de, ProjWData(*data));
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
        DtBkgBB.plotOn(xframe_dt, LineColor(kBlue), ProjWData(*data));
        //DtBkgBB.paramOn(xframe_dt, ProjWData(*data));
        RooHist* pull_dt = xframe_dt->pullHist() ;
        RooPlot* f_pull_dt = dt.frame(Title("")) ;
        f_pull_dt->addPlotable(pull_dt,"BX0") ;
        f_pull_dt->SetMaximum(6);
        f_pull_dt->SetMinimum(-6);
        f_pull_dt->GetYaxis()->SetLabelSize(0.11);
        f_pull_dt->GetXaxis()->SetLabelSize(0.11);
	
	
	
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
	can->SaveAs("plots/fix_shape/bb_shape3D.png");

	//Save fit results
	RooArgSet(fit3D->floatParsFinal()).writeToFile("results/fix_shape/bb_3D.txt");
	
	//Save fit_model in RooWorkspace root file
	RooWorkspace *rws = new RooWorkspace("bb", "bb");
	rws->import(sumMbcBBbkg);
        rws->import(expoDeBBbkg);
        rws->import(DtBkgBB);
	rws->writeToFile("pdf_models/fix_shape/bb_shape3D.root");
	

//----------------------------------------------------------------------------------------------//		

        return 0;
}

