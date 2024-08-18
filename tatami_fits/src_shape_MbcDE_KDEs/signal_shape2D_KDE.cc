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
#include "TObject.h"
#include "libRooTatami/RooTatamiHelper.h"
#include "libRooTatami/RooDtCPSignal.h"
#include "libRooTatami/RooDtLifetime.h"
#include <boost/lexical_cast.hpp>
#include "RooArgusBG.h"
#include "RooExponential.h"
#include "libRooTatami/RooDtBkg.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "TIterator.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"


using namespace std;
using namespace RooFit;


// Main method
int main(int argc, char *argv[]) {

	gROOT->SetBatch(kTRUE);

	//Deal with options
        Options *opt = new Options();
        if (!(opt->ParseOptions(argc, argv)))
        {
                std::exit(EXIT_FAILURE);
        } else
        {
                opt->PrintOptions();
        }


	TChain * tree_sig = new TChain("h1"); tree_sig->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles/for_shape/b_kspipigam_sig_500k_wRMVA_wB2MVA_sel_BCS_branches_Puresignal.root");

	TChain * tree_bb = new TChain("h1"); tree_bb->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles/for_shape/bb_5streams.root");

	TChain * tree_rare = new TChain("h1"); tree_rare->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles/for_shape/merged_MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches_shuffled_49PartsIn50_bkg.root");

TChain * tree_cont = new TChain("h1"); tree_cont->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles/for_shape/continuum_5streams.root");

TChain * tree_scf = new TChain("h1"); tree_scf->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles/for_shape/b_kspipigam_sig_500k_wRMVA_wB2MVA_sel_BCS_branches_SCF.root");

	// Variable definitions
        RooRealVar mbc("mbc","m_{bc} [GeV]",5.20,5.30) ;
        RooRealVar de("de","#DeltaE [GeV]",-0.2,0.2) ;

	RooArgSet mainVariables(mbc, de);	

	// Import ttree
	RooDataSet * data = new RooDataSet("data", "data", mainVariables, Import(*tree));

//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//----------------------------------------------------------------------------------------------//
	//Signal Mbc-dE Model
	RooNDKeysPdf mbcde_sig("mbcde_sig","mbcde_sig",RooArgSet(mbc,de),*data,"am") ;	

	//Save fit_model in RooWorkspace root file
	RooWorkspace *rws = new RooWorkspace("signal_KDE", "signal_KDE");
        rws->import(mbcde_sig);
	rws->writeToFile("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/pdf_models/fix_shape_MbcDE_KDE/signal_shape2D_KDE.root");


//----------------------------------------------------------------------------------------------//		

        return 0;
}

