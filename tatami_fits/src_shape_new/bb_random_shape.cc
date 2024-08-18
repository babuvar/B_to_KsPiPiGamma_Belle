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
#include "RooAbsData.h"

using namespace std;
using namespace RooFit;

bool fixedparameters = false; //Whether or not to fix parameters from Tristan's optimal values

// Main method
int main(int argc, char *argv[]) {

	//Get the input file (one file for now!)
	TChain * tree = new TChain("h1");
	tree->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/rootfiles_for_shape/bb_5streams.root"); //BB bkg weight = 1/5
	int n_bb = tree->GetEntries();
	tree->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/rootfiles_for_shape/MCrare_bkg_finalized_x49_random.root"); //Rare MC weight = 1/49
	int n_rare = tree->GetEntries() - n_bb;

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
	RooRealVar evt_wt("evt_wt","evt_wt", 0.0, 100.0); //In order to scale bb and rare files appropriately

	const double delta_m = 0.507;
	//const double tau_b = 1.5344;	
	const double tau_b = 1.0841;
	
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

        RooArgSet mainVariables(mbc, de, dt, dalitzCategory, evt_wt);
        mainVariables.add(AllConditionalVariables);

	// Import ttree
	Float_t f_mbc, f_de, f_dt, f_dtgen, f_CosThetaB, f_ntrk_sig, f_ntrk_tag, f_z_err_sig, f_z_err_tag, f_chisq_tracks_sig, f_chisq_tracks_tag, f_dgf_tracks_sig, f_dgf_tracks_tag, f_keeptagl_tag, f_w, f_dw, f_ecms, f_issignal, f_issgnevt, f_genBFlavor, f_mc_type;
        Int_t i_inexp, i_iflavor, i_irbin, i_dalitzCategory;


        tree->SetBranchAddress("mbc",&f_mbc);
        tree->SetBranchAddress("de",&f_de);
        tree->SetBranchAddress("dt",&f_dt);
        tree->SetBranchAddress("dtgen",&f_dtgen);
        tree->SetBranchAddress("dalitzCategory",&i_dalitzCategory);
        tree->SetBranchAddress("inexp",&i_inexp);
        tree->SetBranchAddress("coscms",&f_CosThetaB);
        tree->SetBranchAddress("vtxntrk",&f_ntrk_sig);
        tree->SetBranchAddress("tagntrk",&f_ntrk_tag);
        tree->SetBranchAddress("vtxzerr",&f_z_err_sig);
        tree->SetBranchAddress("tagzerr",&f_z_err_tag);
        tree->SetBranchAddress("vtxchi2",&f_chisq_tracks_sig);
        tree->SetBranchAddress("tagchi2",&f_chisq_tracks_tag);
        tree->SetBranchAddress("vtxndf",&f_dgf_tracks_sig);
        tree->SetBranchAddress("tagndf",&f_dgf_tracks_tag);
        tree->SetBranchAddress("taglepto",&f_keeptagl_tag);
        tree->SetBranchAddress("iflavor",&i_iflavor);
        tree->SetBranchAddress("w",&f_w);
        tree->SetBranchAddress("dw",&f_dw);
        tree->SetBranchAddress("ecm",&f_ecms);
        tree->SetBranchAddress("irbin",&i_irbin);
        tree->SetBranchAddress("issignal",&f_issignal);
        tree->SetBranchAddress("issgnevt",&f_issgnevt);
        tree->SetBranchAddress("mc_type",&f_mc_type);
	tree->SetBranchAddress("dtgen",&f_dtgen);
        tree->SetBranchAddress("genbq",&f_genBFlavor);
	
	float ndat = 0.0, ndat_w = 0.0;

	//loop over tree
	RooDataSet * data = new RooDataSet("data", "data", mainVariables);	
	for(int i=0; i<tree->GetEntries(); i++){

	tree->GetEntry(i);

	if(abs(f_dt) < 10.0){
	mbc.setVal(f_mbc);
        de.setVal(f_de);
        dt.setVal(f_dt);
        dalitzCategory.setIndex(i_dalitzCategory);
        inexp.setIndex(i_inexp);
        CosThetaB.setVal(f_CosThetaB);
        ntrk_sig.setVal(f_ntrk_sig);
        ntrk_tag.setVal(f_ntrk_tag);
        z_err_sig.setVal(f_z_err_sig);
        z_err_tag.setVal(f_z_err_tag);
        chisq_tracks_sig.setVal(f_chisq_tracks_sig);
        chisq_tracks_tag.setVal(f_chisq_tracks_tag);
        dgf_tracks_sig.setVal(f_dgf_tracks_sig);
        dgf_tracks_tag.setVal(f_dgf_tracks_tag);
        keeptagl_tag.setVal(f_keeptagl_tag);
        iflavor.setIndex(i_iflavor);
	ecms.setVal(f_ecms);
        irbin.setIndex(i_irbin);
        issignal.setVal(f_issignal);
        issgnevt.setVal(f_issgnevt);

	//weight factors 1/5 and 1/49 are applied in order to get the right proportion of rare and bb bkg
	//The additional factor 26.65 is added so that the total statistics of the weighted sample
	// is equal to the total statistics of MCRAREx49 + BBbkgx5. This matters in the  calculation
	// of parameter uncertainties
	if(i < n_bb){evt_wt.setVal(26.65/5.0);}
	else{evt_wt.setVal(26.65/49.0);}
	

        data->add(mainVariables);
	ndat = ndat + 1.0; ndat_w = ndat_w + evt_wt.getVal();

        }//explicitly ensure fit range in Delta-t
        }

	cout<<"ndat = "<<ndat<<endl;
        cout<<"ndat_w = "<<ndat_w<<endl;
	
	//Weighted dataset
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,evt_wt.GetName()) ;

//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//                                                                                              //
//  The BBbar is modelled as :                                                              //
//  Mbc : Gaussian + Argus                                                                      //
//  dE  : Exponential                                                                           //
//  dT  : Tatami resolution function                                                            //
//----------------------------------------------------------------------------------------------//

	//MBC
	//Argus
	RooRealVar bbRandMBC_c("bbRandMBC_c","bbRandMBC_c", -70, -140, 20);
        RooRealVar bbRandMBC_m0("bbRandMBC_m0","bbRandMBC_m0", 5.29, 5.28, 5.30) ;
	RooArgusBG bbRandMBC("bbRandMBC","bbRandMBC", mbc, bbRandMBC_m0, bbRandMBC_c);

        //dE
        RooRealVar bbRandDE_al("bbRandDE_al", "bbRandDE_al", -6.01751, -10., 0.) ;
        RooExponential bbRandDE("bbRandDE", "bbRandDE", de, bbRandDE_al) ;

        //dt
        RooRealVar bbRandDT_tau("bbRandDT_tau", "bbRandDT_tau", 1.40, 0.5, 1.8);
        RooRealVar bbRandDT_mu_l("bbRandDT_mu_l", "bbRandDT_mu_l", -0.698928, -0.7, 0.5);
        RooRealVar bbRandDT_mu_d("bbRandDT_mu_d", "bbRandDT_mu_d", 0.0028751, -0.1, 0.1);
        RooRealVar bbRandDT_f_delt("bbRandDT_f_delt", "bbRandDT_f_delt", 0.150, 0.1,1);
        RooRealVar bbRandDT_f_tail("bbRandDT_f_tail", "bbRandDT_f_tail", 0.19111, 0, 1);
        RooRealVar bbRandDT_S_main("bbRandDT_S_main", "bbRandDT_S_main", 2.8550, 0.4, 10);
        RooRealVar bbRandDT_S_tail("bbRandDT_S_tail", "bbRandDT_S_tail", 3.4377, 1, 40);

        RooDtBkg bbRandDT("bbRandDT", "bbRandDT", dt, bbRandDT_tau, bbRandDT_mu_l, bbRandDT_mu_d, bbRandDT_f_delt, bbRandDT_f_tail, bbRandDT_S_main, bbRandDT_S_tail, inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);

	//RooDtCPSignal bbRandDT("bbRandDT", "bbRandDT", dt, S, A, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);

	//RooDtLifetime bbRandDT("bbRandDT","bbRandDT", dt, tauB, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, true, true, 1.0);

	//If parameters are fixed
	if(fixedparameters == true){
	//Mbc
	bbRandMBC_c.setConstant(true);
	bbRandMBC_m0.setConstant(true);
	//dE
	bbRandDE_al.setConstant(true);
	//dt
	bbRandDT_tau.setConstant(true);
	bbRandDT_mu_l.setConstant(true);
	bbRandDT_mu_d.setConstant(true);
	bbRandDT_f_delt.setConstant(true);
	bbRandDT_f_tail.setConstant(true);
	bbRandDT_S_main.setConstant(true);
	bbRandDT_S_tail.setConstant(true);
	}


//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//
	
	//bool w2err = false;
	bool w2err = true;

	//Fitting
	RooFitResult * fitMbcbkg = bbRandMBC.fitTo(wdata,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0), SumW2Error(w2err)) ;
	RooFitResult * fitDebkg = bbRandDE.fitTo(wdata,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0), SumW2Error(w2err)) ;
	RooFitResult * fitDtbkg = bbRandDT.fitTo(wdata,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0), SumW2Error(w2err)) ;
	

	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	wdata.plotOn(xframe_mbc);
	bbRandMBC.plotOn(xframe_mbc, LineColor(kBlue), ProjWData(wdata));
	RooHist* pull_mbc = xframe_mbc->pullHist() ;
	RooPlot* f_pull_mbc = mbc.frame(Title("")) ;
	f_pull_mbc->addPlotable(pull_mbc,"BX0") ;
	f_pull_mbc->SetMaximum(4);
	f_pull_mbc->SetMinimum(-4);
	f_pull_mbc->GetYaxis()->SetLabelSize(0.11);
        f_pull_mbc->GetXaxis()->SetLabelSize(0.11);


	//dE
	RooPlot *xframe_de = de.frame(Bins(40));
	wdata.plotOn(xframe_de);
	bbRandDE.plotOn(xframe_de, LineColor(kBlue), ProjWData(wdata));
	RooHist* pull_de = xframe_de->pullHist() ;
        RooPlot* f_pull_de = de.frame(Title("")) ;
        f_pull_de->addPlotable(pull_de,"BX0") ;
        f_pull_de->SetMaximum(4);
        f_pull_de->SetMinimum(-4);
	f_pull_de->GetYaxis()->SetLabelSize(0.11);
        f_pull_de->GetXaxis()->SetLabelSize(0.11);	

	//dt
	RooPlot *xframe_dt = dt.frame(Bins(40));
        wdata.plotOn(xframe_dt);
        bbRandDT.plotOn(xframe_dt, LineColor(kBlue), ProjWData(wdata));
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
	can->SaveAs("plots/fix_shape_new/bb_random_shape.png");
	
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
	can->SaveAs("plots/fix_shape_new/bb_random_shape.png");
	*/	

	//Save fit results
	RooArgSet(fitMbcbkg->floatParsFinal()).writeToFile("results/fix_shape_new/bb_random_Mbc.txt");
	RooArgSet(fitDebkg->floatParsFinal()).writeToFile("results/fix_shape_new/bb_random_dE.txt");
	RooArgSet(fitDtbkg->floatParsFinal()).writeToFile("results/fix_shape_new/bb_random_dt.txt");
	
	//Save fit_model in RooWorkspace root file
	RooWorkspace *rws = new RooWorkspace("bbRand", "bbRand");
	rws->import(bbRandMBC); 
        rws->import(bbRandDE);
        rws->import(bbRandDT);
	rws->writeToFile("pdf_models/fix_shape_new/bb_random_shape.root");


	//Generate randomized parameters from RooFitResult and save in a root-file 
	RooDataSet randPars_mbc("mbc","mbc",fitMbcbkg->floatParsFinal()) ;
	RooDataSet randPars_de("de", "de", fitDebkg->floatParsFinal()) ;
	RooDataSet randPars_dt("dt", "dt", fitDtbkg->floatParsFinal()) ;
	for (Int_t i=0 ; i<1000 ; i++) {
	randPars_mbc.add(fitMbcbkg->randomizePars()) ;
	randPars_de.add(fitDebkg->randomizePars()) ;
	randPars_dt.add(fitDtbkg->randomizePars()) ;}
	
	RooWorkspace *w_randPars = new RooWorkspace("rws_bb_randPars", "rws_bb_randPars");
	w_randPars->import(randPars_mbc);
        w_randPars->import(randPars_de);
        w_randPars->import(randPars_dt);
	w_randPars->writeToFile("RandomizedParameters/fix_shape_new/bb_randPars.root");


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

