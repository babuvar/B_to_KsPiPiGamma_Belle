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


using namespace std;
using namespace RooFit;

bool fixedparameters = false; //Whether or not to fix parameters from Tristan's optimal values

// Main method
int main(int argc, char *argv[]) {


	//Get the input file (one file for now!)
	TChain * tree = new TChain("h1");
	tree->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/rootfiles_for_shape/continuum_5streams.root");
	

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
	RooDataSet * data = new RooDataSet("data", "data", mainVariables, Import(*tree), Cut("abs(dt) < 10.0"));


//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//                                                                                              //
//  The continuum is modelled as :                                                              //
//  Mbc : Argus                                                                                 //
//  dE  : Exponential                                                                           //
//  dT  : Tatami resolution function                                                            //
//----------------------------------------------------------------------------------------------//

	//MBC
	RooRealVar  contMBC_c("contMBC_c", "contMBC_c", -7.39000, -200, -1);
        RooRealVar  contMBC_m0("contMBC_m0","contMBC_m0", 5.2890, 5.28, 5.30);
	RooArgusBG  contMBC("contMBC","contMBC", mbc, contMBC_m0, contMBC_c);

	//dE
	RooRealVar  contDE_al("contDE_al", "contDE_al", -0.915794,-10.,0.1) ;
        RooExponential contDE("contDE","contDE", de, contDE_al);

	//dt
	RooRealVar contDT_tau("contDT_tau", "contDT_tau", 1.4, 0.5, 1.8);
        RooRealVar contDT_mu_l("contDT_mu_l", "contDT_mu_l", -0.699099, -0.7, 0.5);
        RooRealVar contDT_mu_d("contDT_mu_d", "contDT_mu_d", -0.0126835, -0.1, 0.1);
        RooRealVar contDT_f_delt("contDT_f_delt", "contDT_f_delt", 0.110, 0.1, 1);
        RooRealVar contDT_f_tail("contDT_f_tail", "contDT_f_tail", 0.16327, 0, 1);
        RooRealVar contDT_S_main("contDT_S_main", "contDT_S_main", 1.7924, 0.4, 10);
        RooRealVar contDT_S_tail("contDT_S_tail" , "contDT_S_tail", 3.9538, 1, 40);

        RooDtBkg contDT("contDT", "contDT", dt,
                contDT_tau, contDT_mu_l, contDT_mu_d, contDT_f_delt, contDT_f_tail, contDT_S_main, contDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);



	//If parameters are fixed
	if(fixedparameters == true){
	//Mbc
	contMBC_c.setConstant(true);
	contMBC_m0.setConstant(true);
	//dE
	contDE_al.setConstant(true);
	//dt
	contDT_tau.setConstant(true);
	contDT_mu_l.setConstant(true);
	contDT_mu_d.setConstant(true);
	contDT_f_delt.setConstant(true);
	contDT_f_tail.setConstant(true);
	contDT_S_main.setConstant(true);
	contDT_S_tail.setConstant(true);
	}


//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//
	
	//Fitting
	RooFitResult * fitMbcbkg = contMBC.fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	RooFitResult * fitDebkg = contDE.fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	RooFitResult * fitDtbkg = contDT.fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	

	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	data->plotOn(xframe_mbc);
	contMBC.plotOn(xframe_mbc, LineColor(kBlue), ProjWData(*data));
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
	contDE.plotOn(xframe_de, LineColor(kBlue), ProjWData(*data));
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
        contDT.plotOn(xframe_dt, LineColor(kBlue), ProjWData(*data));
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
	can->SaveAs("plots/fix_shape_new/cont_shape.png");
	
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
	can->SaveAs("plots/fix_shape_new/cont_shape.png");
	*/

	//Save fit results
	RooArgSet(fitMbcbkg->floatParsFinal()).writeToFile("results/fix_shape_new/cont_Mbc.txt");
	RooArgSet(fitDebkg->floatParsFinal()).writeToFile("results/fix_shape_new/cont_dE.txt");
	RooArgSet(fitDtbkg->floatParsFinal()).writeToFile("results/fix_shape_new/cont_dt.txt");
	
	//Save fit_model in RooWorkspace root file
	RooWorkspace *rws = new RooWorkspace("continuum", "continuum");
	rws->import(contMBC); 
        rws->import(contDE);
        rws->import(contDT);
	rws->writeToFile("pdf_models/fix_shape_new/cont_shape.root");


	//Generate randomized parameters from RooFitResult and save in a root-file 
        RooDataSet randPars_mbc("mbc","mbc",fitMbcbkg->floatParsFinal()) ;
        RooDataSet randPars_de("de", "de", fitDebkg->floatParsFinal()) ;
        RooDataSet randPars_dt("dt", "dt", fitDtbkg->floatParsFinal()) ;
        for (Int_t i=0 ; i<1000 ; i++) {
        randPars_mbc.add(fitMbcbkg->randomizePars()) ;
        randPars_de.add(fitDebkg->randomizePars()) ;
        randPars_dt.add(fitDtbkg->randomizePars()) ;}

        RooWorkspace *w_randPars = new RooWorkspace("rws_cont_randPars", "rws_cont_randPars");
        w_randPars->import(randPars_mbc);
        w_randPars->import(randPars_de);
        w_randPars->import(randPars_dt);
        w_randPars->writeToFile("RandomizedParameters/fix_shape_new/cont_randPars.root");

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

