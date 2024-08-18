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
	TChain * tree = new TChain("variables");
	tree->Add("/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle2/bootstrapping/parent_samples/RecoDecaysMC_reduced_cols_wRMVA_BCS.root");	

	//Get ntuple ready
	double d_issignal, d_mcdt, d_dt, d_dterr, d_qr, d_delE;
	Int_t i_mctype;

	tree->SetBranchAddress("B0_isSignal", &d_issignal);
	tree->SetBranchAddress("B0_DeltaT", &d_dt);
	tree->SetBranchAddress("B0_DeltaTErr", &d_dterr);
	tree->SetBranchAddress("B0_mcDeltaT", &d_mcdt);
	tree->SetBranchAddress("B0_qr", &d_qr);
	tree->SetBranchAddress("MCtype", &i_mctype);
	tree->SetBranchAddress("B0_deltaE", &d_delE);

	//Define variables for fit
	RooCategory inexp("inexp", "inexp");
        inexp.defineType("7", 7);
        inexp.defineType("31", 31);
	RooRealVar dt("dt", "dt", -15, 15);
	RooRealVar dtErr("dtErr", "dtErr", 0.01, 2);
	RooRealVar dummy("dummy", "dummy", 0.01, 2);
	RooArgSet mainVariables(dt);
	RooArgSet AllConditionalVariables(dtErr, dummy, inexp);
	mainVariables.add(AllConditionalVariables);

	//Fill dataset
	RooDataSet * data = new RooDataSet("data", "data", mainVariables);
	for(int i = 0; i < tree->GetEntries(); i++){//tree-loop

	tree->GetEntry(i);
	if(i_mctype > 3){//select continuum
	if(abs(d_dt) < 15 && d_dterr > 0.01 && d_dterr < 2.0){
	if(d_delE > -0.5 && d_delE < 0.5){

	dummy.setVal(0.0);
	dtErr.setVal(d_dterr);
	dt.setVal(d_dt);
	inexp.setIndex(7);

	data->add(mainVariables);

	}
	}//abs(dt)
	}//select continuum

	}//tree-loop




//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//  dT  : Tatami resolution function                                                            //
//----------------------------------------------------------------------------------------------//


	//dt
	RooRealVar contDT_tau("contDT_tau", "contDT_tau", 0.0);//1.4, 0.5, 1.8);
        RooRealVar contDT_mu_l("contDT_mu_l", "contDT_mu_l", 0.0);//-0.699099, -0.7, 0.5);
        RooRealVar contDT_mu_d("contDT_mu_d", "contDT_mu_d", -0.0126835, -0.1, 0.1);
        RooRealVar contDT_f_delt("contDT_f_delt", "contDT_f_delt", 0.0);//0.110, 0.1, 1);
        RooRealVar contDT_f_tail("contDT_f_tail", "contDT_f_tail", 0.16327, 0, 1);
        RooRealVar contDT_S_main("contDT_S_main", "contDT_S_main", 0.0017924, 0.0004, 0.020);
        RooRealVar contDT_S_tail("contDT_S_tail" , "contDT_S_tail", 3.9538, 1, 40);

        RooDtBkg contDT("contDT", "contDT", dt, contDT_tau, contDT_mu_l, contDT_mu_d, contDT_f_delt, contDT_f_tail, contDT_S_main, contDT_S_tail, inexp, dummy, dummy, dummy, dtErr, true);

        //RooDtBkg contDT("contDT", "contDT", dt, contDT_tau, contDT_mu_l, contDT_mu_d, contDT_f_delt, contDT_f_tail, contDT_S_main, contDT_S_tail, inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);



	//If parameters are fixed
	if(fixedparameters == true){
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
	RooFitResult * fitDtbkg = contDT.fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	

	//Plotting

	//dt
	RooPlot *xframe_dt = dt.frame(Bins(60));
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
	
	TCanvas* can = new TCanvas("c","c",900,800) ; 
	TPad *pad1 = new TPad("pad1", "pad 1",0.0,0.2,1.0,1.0,0);
        TPad *pad2 = new TPad("pad2", "pad 2",0.0,0.0,1.0,0.2,0);
	pad1->Draw(); 	pad2->Draw();   
	//Draw
        pad1->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_dt->Draw() ;
        pad2->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_dt->Draw() ;
	can->SaveAs("plots/belle2_cont_dt.png");

	//Logplot
	pad1->cd();
	gPad->SetLogy();
	xframe_dt->SetMinimum(0.1);
	can->SaveAs("plots/belle2_cont_dt_log.png");

	//Save fit results
	RooArgSet(fitDtbkg->floatParsFinal()).writeToFile("results/belle2_cont_dt.txt");
	


//----------------------------------------------------------------------------------------------//		

        return 0;
}

