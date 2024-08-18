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
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "libRooTatami/RooDtBkg.h"
#include "RooWorkspace.h"
#include "RooGaussian.h"


using namespace std;
using namespace RooFit;

bool fixedparameters = false; //Whether or not to fix parameters from Tristan's optimal values


// Main method
int main(int argc, char *argv[]) {


	//Get the input file (one file for now!)
	TChain * tree = new TChain("h1");
	tree->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/rootfiles_for_shape/MCrare_bkg_finalized_x49.root");
	

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
        RooRealVar tauB("\\tau_{B}", "tauB", 0.8, 2.0);
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
	RooRealVar mcpdg("mcpdg", "mcpdg", -310000, 310000);

	const double delta_m = 0.507;
	//const double tau_b = 1.5344;	
	//const double tau_b = 1.5905;	

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
	AllConditionalVariables.add(mcpdg);

        RooArgSet mainVariables(mbc, de, dt, dalitzCategory);
        mainVariables.add(AllConditionalVariables);

	// Import ttree
	RooDataSet * data = new RooDataSet("data", "data", mainVariables, Import(*tree), Cut("abs(dt) < 10.0 && (abs(mcpdg) == 511)"));


//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//                                                                                              //
//  The BB missingFSP bkg is modelled as :                                                               //
//  Mbc : Argus + crystalball                                                                   //
//  dE  : Exponential                                                                           //
//  dT  : Tatami resolution function                                                            //
//----------------------------------------------------------------------------------------------//

	//MBC
	//CrystalBall
	RooRealVar bb_missFSP_MBC_mean("bb_missFSP_MBC_mean", "bb_missFSP_MBC_mean", 5.2803, 5.275, 5.285) ;
        RooRealVar bb_missFSP_MBC_sigma("bb_missFSP_MBC_sigma", "bb_missFSP_MBC_sigma", 0.0033531, 5e-4, 5e-3);
        RooRealVar bb_missFSP_MBC_alpha("bb_missFSP_MBC_alpha", "bb_missFSP_MBC_alpha", 0.51645, 0.,5.);
        RooRealVar bb_missFSP_MBC_n("bb_missFSP_MBC_n", "bb_missFSP_MBC_n", 14.000, 1, 1000);
        RooCBShape bb_missFSP_MBC("bb_missFSP_MBC", "bb_missFSP_MBC", mbc, bb_missFSP_MBC_mean, bb_missFSP_MBC_sigma, bb_missFSP_MBC_alpha, bb_missFSP_MBC_n);

	//dE
	//Exponential
        RooRealVar bb_missFSP_DE_alpha("bb_missFSP_DE_alpha", "bb_missFSP_DE_alpha", -3.95308, -20., 0.1);
        RooExponential bb_missFSP_DE_expo("bb_missFSP_DE_expo", "bb_missFSP_DE_expo", de, bb_missFSP_DE_alpha);
	//Gaussian
	RooRealVar bb_missFSP_DE_mean("bb_missFSP_DE_mean", "bb_missFSP_DE_mean", -0.1, -0.2, 0.05);
        RooRealVar bb_missFSP_DE_sigma("bb_missFSP_DE_sigma", "bb_missFSP_DE_sigma", 0.04, 0.001, 0.2);
        RooGaussian bb_missFSP_DE_gauss("bb_missFSP_DE_gauss", "bb_missFSP_DE_gauss", de, bb_missFSP_DE_mean, bb_missFSP_DE_sigma);
	//Combined
	RooRealVar bb_missFSP_DE_f("bb_missFSP_DE_f", "bb_missFSP_DE_f", 0.3, 0.01, 0.99);
        RooAddPdf bb_missFSP_DE("bb_missFSP_DE", "bb_missFSP_DE", RooArgList(bb_missFSP_DE_expo, bb_missFSP_DE_gauss), RooArgList(bb_missFSP_DE_f)) ;

	//dt
	RooRealVar bb_missFSP_DT_tau("bb_missFSP_DT_tau", "bb_missFSP_DT_tau", 1.40, 0.5,1.8);
        RooRealVar bb_missFSP_DT_mu_l("bb_missFSP_DT_mu_l", "bb_missFSP_DT_mu_l", -0.798748,-0.8,0.5);
        RooRealVar bb_missFSP_DT_mu_d("bb_missFSP_DT_mu_d", "bb_missFSP_DT_mu_d", -0.0776618, -0.1,0.1);
        RooRealVar bb_missFSP_DT_f_delt("bb_missFSP_DT_f_delt", "bb_missFSP_DT_f_delt", 0.150, 0.1,1);
        RooRealVar bb_missFSP_DT_f_tail("bb_missFSP_DT_f_tail", "bb_missFSP_DT_f_tail", 0.32953, 0, 1.);
        RooRealVar bb_missFSP_DT_S_main("bb_missFSP_DT_S_main", "bb_missFSP_DT_S_main", 2.4685, 0.4,10);
        RooRealVar bb_missFSP_DT_S_tail("bb_missFSP_DT_S_tail", "bb_missFSP_DT_S_tail", 3.3250, 1, 40);

        //RooDtBkg bb_missFSP_DT("bb_missFSP_DT", "bb_missFSP_DT", dt, bb_missFSP_DT_tau, bb_missFSP_DT_mu_l, bb_missFSP_DT_mu_d, bb_missFSP_DT_f_delt, bb_missFSP_DT_f_tail, bb_missFSP_DT_S_main, bb_missFSP_DT_S_tail, inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);

	
	const double tau_b = 1.4896;  
        RooDtCPSignal bb_missFSP_DT("bb_missFSP_DT", "bb_missFSP_DT", dt, S, A, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);
	        
        //RooDtLifetime bb_missFSP_DT("bb_missFSP_DT","bb_missFSP_DT", dt, tauB, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, true, true, 1.0);
	        //

	//If parameters are fixed
	if(fixedparameters == true){
        //Mbc
        bb_missFSP_MBC_mean.setConstant(true);
        bb_missFSP_MBC_sigma.setConstant(true);
        bb_missFSP_MBC_alpha.setConstant(true);
        bb_missFSP_MBC_n.setConstant(true);
        //dE
        bb_missFSP_DE_alpha.setConstant(true);
        bb_missFSP_DE_mean.setConstant(true);
        bb_missFSP_DE_sigma.setConstant(true);
        bb_missFSP_DE_f.setConstant(true);
        //dt
        bb_missFSP_DT_tau.setConstant(true);
        bb_missFSP_DT_mu_l.setConstant(true);
        bb_missFSP_DT_mu_d.setConstant(true);
        bb_missFSP_DT_f_delt.setConstant(true);
        bb_missFSP_DT_f_tail.setConstant(true);
        bb_missFSP_DT_S_main.setConstant(true);
        bb_missFSP_DT_S_tail.setConstant(true);
	}

//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//
	
	//Fitting
	RooFitResult * fitMbcbkg = bb_missFSP_MBC.fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	RooFitResult * fitDebkg = bb_missFSP_DE.fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	RooFitResult * fitDtbkg = bb_missFSP_DT.fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	

	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	data->plotOn(xframe_mbc);
	bb_missFSP_MBC.plotOn(xframe_mbc, LineColor(kBlue), ProjWData(*data));
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
	bb_missFSP_DE.plotOn(xframe_de, LineColor(kBlue), ProjWData(*data));
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
        bb_missFSP_DT.plotOn(xframe_dt, LineColor(kBlue), ProjWData(*data));
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
	can->SaveAs("plots/fix_shape_new/bb_missFSP.png");
	
	/*
	//Show only Mbc and Delta-E, which are actually fixed
	TCanvas* can = new TCanvas("c","c",1300,800) ; can->Divide(2,1);
        TPad *pad1 = new TPad("pad1", "pad 1",0.0,0.0,0.5,1.0,0);
        TPad *pad2 = new TPad("pad2", "pad 2",0.5,0.0,1.0,1.0,0);
        pad1->Draw();   pad2->Draw();
        pad1->cd();
        TPad *pad11 = new TPad("pad11", "pad 11",0.0,0.2,1.0,1.0,0);
        TPad *pad12 = new TPad("pad12", "pad 12",0.0,0.0,1.0,0.2,0);
        pad11->Draw(); pad12->Draw();
        pad2->cd();
        TPad *pad21 = new TPad("pad21", "pad 21",0.0,0.2,1.0,1.0,0);
        TPad *pad22 = new TPad("pad22", "pad 22",0.0,0.0,1.0,0.2,0);
        pad21->Draw(); pad22->Draw();

	//Draw
	pad11->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_mbc->Draw() ;
        pad12->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_mbc->Draw() ;
        pad21->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_de->Draw() ;
        pad22->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_de->Draw() ;
	can->SaveAs("plots/fix_shape_new/bb_missFSP.png");
	*/
	//Save fit results
	RooArgSet(fitMbcbkg->floatParsFinal()).writeToFile("results/fix_shape_new/bb_missFSP_Mbc.txt");
	RooArgSet(fitDebkg->floatParsFinal()).writeToFile("results/fix_shape_new/bb_missFSP_dE.txt");
	RooArgSet(fitDtbkg->floatParsFinal()).writeToFile("results/fix_shape_new/bb_missFSP_dt.txt");

	//Save fit_model in RooWorkspace root file
	RooWorkspace *rws = new RooWorkspace("bb_missFSP", "bb_missFSP");
	rws->import(bb_missFSP_MBC);
        rws->import(bb_missFSP_DE);
        rws->import(bb_missFSP_DT);
	rws->writeToFile("pdf_models/fix_shape_new/bb_missFSP_shape.root");

	
	//Generate randomized parameters from RooFitResult and save in a root-file 
        RooDataSet randPars_mbc("mbc","mbc",fitMbcbkg->floatParsFinal()) ;
        RooDataSet randPars_de("de", "de", fitDebkg->floatParsFinal()) ;
        RooDataSet randPars_dt("dt", "dt", fitDtbkg->floatParsFinal()) ;
        for (Int_t i=0 ; i<1000 ; i++) {
        randPars_mbc.add(fitMbcbkg->randomizePars()) ;
        randPars_de.add(fitDebkg->randomizePars()) ;
        randPars_dt.add(fitDtbkg->randomizePars()) ;}

        RooWorkspace *w_randPars = new RooWorkspace("rws_bb_missFSP_randPars", "rws_bb_missFSP_randPars");
        w_randPars->import(randPars_mbc);
        w_randPars->import(randPars_de);
        w_randPars->import(randPars_dt);
        w_randPars->writeToFile("RandomizedParameters/fix_shape_new/bb_missFSP_randPars.root");

	//Print final Correlation and Covariance matrices
	//MBC
	const TMatrixDSym& cor_mbc =fitMbcbkg->correlationMatrix();
        const TMatrixDSym& cov_mbc =fitMbcbkg->covarianceMatrix() ;
        cout<<endl<<"***correlationMBC***"<<endl;
        cor_mbc.Print();
        cout<<endl<<"***correlationMBC***"<<endl;
        cout<<endl<<"***covarianceMBC***"<<endl;
        cov_mbc.Print();
        cout<<endl<<"***covarianceMBC***"<<endl;
        //DE
        const TMatrixDSym& cor_de =fitDebkg->correlationMatrix();
        const TMatrixDSym& cov_de =fitDebkg->covarianceMatrix() ;
        cout<<endl<<"***correlationDE***"<<endl;
        cor_de.Print();
        cout<<endl<<"***correlationDE***"<<endl;
        cout<<endl<<"***covarianceDE***"<<endl;
        cov_de.Print();
        cout<<endl<<"***covarianceDE***"<<endl;
        //DT
        const TMatrixDSym& cor_dt =fitDtbkg->correlationMatrix();
        const TMatrixDSym& cov_dt =fitDtbkg->covarianceMatrix() ;
        cout<<endl<<"***correlationDT***"<<endl;
        cor_dt.Print();
        cout<<endl<<"***correlationDT***"<<endl;
        cout<<endl<<"***covarianceDT***"<<endl;
        cov_dt.Print();
        cout<<endl<<"***covarianceDT***"<<endl;



//----------------------------------------------------------------------------------------------//		

        return 0;
}

