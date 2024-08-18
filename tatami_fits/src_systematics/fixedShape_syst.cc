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
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooExponential.h"
#include "libRooTatami/RooDtBkg.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "TIterator.h"
#include "RooAddPdf.h"
#include "TF1.h"
#include "RooCFunction1Binding.h"
#include "RooFormulaVar.h"
//#include "RooCFunction3Binding.h"
#include "RooChebychev.h"
#include "RooBifurGauss.h"

using namespace std;
using namespace RooFit;

float randParVal(RooWorkspace * w, char* dimension, int toyIndex,char * parameter )
{
RooAbsData *randpars = w->data(dimension);
RooRealVar * par = (RooRealVar*) (randpars->get(toyIndex)->find(parameter) );

return par->getVal();
}


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

	//Get the input file (one file for now!)
	string filename = opt->GetFilenames().at(0);

	//Get the names of the output plot, result file, consolidated result directory and toy-fit index
	string canvasName = opt->GetPlotname();
	string resultName = opt->GetResultname();
	string ConResDir = opt->GetConsolidatedResultsDir();
	string component = opt->GetSMCSample();
	int toyindex = opt->GetToyIndex();

	TChain * tree = new TChain("h1");
	tree->Add(filename.c_str());

	//Set tree branch addressed to Float_t variables
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
	//const double tau_b = 1.525;

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



        tree->SetBranchAddress("issignal",&f_issignal);
        tree->SetBranchAddress("issgnevt",&f_issgnevt);
        tree->SetBranchAddress("mc_type",&f_mc_type);



	//Mistag rates 
	float dw_true[7][2], w_true[7][2];
	
	
	//Belle MC mistag calibration
	
        w_true[0][0] = 0.5; w_true[0][1] = 0.5; dw_true[0][0] = 0.0; dw_true[0][1] = 0.0;
        w_true[1][0] = 0.420827; w_true[1][1] = 0.412222; dw_true[1][0] = 0.0583019; dw_true[1][1] = 0.00408778;
        w_true[2][0] = 0.300296; w_true[2][1] = 0.307838; dw_true[2][0] = 0.00573998; dw_true[2][1] = 0.010326;
        w_true[3][0] = 0.219317; w_true[3][1] = 0.212765; dw_true[3][0] = -0.0392635; dw_true[3][1] = -0.00479522;
        w_true[4][0] = 0.154636; w_true[4][1] = 0.149933; dw_true[4][0] = 0.00474508; dw_true[4][1] = 0.00151989;
        w_true[5][0] = 0.0916131; w_true[5][1] = 0.0913264; dw_true[5][0] = -0.0118737; dw_true[5][1] = 0.0143633;
        w_true[6][0] = 0.0228891; w_true[6][1] = 0.0218754; dw_true[6][0] = -0.00585326; dw_true[6][1] = 0.00189979;
		

	// Import ttree
RooDataSet * data = new RooDataSet("data", "data", mainVariables);
        for(int i=0; i<tree->GetEntries(); i++){

        tree->GetEntry(i);

	if(i_irbin > 0){//omit 0th r bin
	if(abs(f_dt) < 10.0){//explicitly ensure fit range in Delta-t
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
	//correct avg-values from MC-truth
	int exp_bin = int(i_inexp/30);
        w.setVal(w_true[i_irbin][exp_bin]);
        dw.setVal(dw_true[i_irbin][exp_bin]);

        ecms.setVal(f_ecms);
        irbin.setIndex(i_irbin);
        issignal.setVal(f_issignal);
        issgnevt.setVal(f_issgnevt);

        data->add(mainVariables);
	}//explicitly ensure fit range in Delta-t
        }//omit 0th r bin	

        }//tree-reading loop


//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//----------------------------------------------------------------------------------------------//

	//Full model
	
//------------------------------------------------------------------------------------------------------
	
	RooWorkspace rws_fixed("rws_fixed", "rws_fixed");

	//Load rws to get fixed shape parameters
	//Import Mbc pdfs
 	rws_fixed.import("pdf_models/fix_shape_new/cont_shape.root:continuum:contMBC");
        rws_fixed.import("pdf_models/fix_shape_new/bb_missFSP_shape.root:bb_missFSP:bb_missFSP_MBC");
        rws_fixed.import("pdf_models/fix_shape_new/bb_random_shape.root:bbRand:bbRandMBC");
        rws_fixed.import("pdf_models/fix_shape_new/scf_shape.root:scf:scfMBC");
        rws_fixed.import("pdf_models/fix_shape_new/signal_pure_shape.root:signal:signalMBC");

        //Import dE pdfs
	rws_fixed.import("pdf_models/fix_shape_new/cont_shape.root:continuum:contDE");
        rws_fixed.import("pdf_models/fix_shape_new/bb_missFSP_shape.root:bb_missFSP:bb_missFSP_DE");
        rws_fixed.import("pdf_models/fix_shape_new/bb_random_shape.root:bbRand:bbRandDE");
        rws_fixed.import("pdf_models/fix_shape_new/scf_shape.root:scf:scfDE");
        rws_fixed.import("pdf_models/fix_shape_new/signal_pure_shape.root:signal:signalDE");

        //Import dt pdfs
        rws_fixed.import("pdf_models/fix_shape_new/cont_shape.root:continuum:contDT");
        rws_fixed.import("pdf_models/fix_shape_new/bb_missFSP_shape.root:bb_missFSP:bb_missFSP_DT");
        rws_fixed.import("pdf_models/fix_shape_new/bb_random_shape.root:bbRand:bbRandDT");
        rws_fixed.import("pdf_models/fix_shape_new/scf_shape.root:scf:scfDT");
        rws_fixed.import("pdf_models/fix_shape_new/signal_pure_shape.root:signal:signalDT");


//------------------------------------------------------------------------------------------------------

 	//Common Argus-threshold
	RooRealVar commonArgus_m0("commonArgus_m0", "commonArgus_m0", 5.29, 5.28, 5.30) ;

	//BB Random-bkg
	//MBC
	//Argus
	RooRealVar *bbRandMBC_c = rws_fixed.var("bbRandMBC_c"); bbRandMBC_c->setConstant(true);
	RooRealVar *bbRandMBC_m0 = rws_fixed.var("bbRandMBC_m0"); bbRandMBC_m0->setConstant(true);
        //RooArgusBG bbRandMBC("bbRandMBC","bbRandMBC", mbc, bbRandMBC_m0, *bbRandMBC_c);
        RooArgusBG bbRandMBC("bbRandMBC","bbRandMBC", mbc, commonArgus_m0, *bbRandMBC_c);

	//dE
	RooRealVar *bbRandDE_al = rws_fixed.var("bbRandDE_al"); bbRandDE_al->setConstant(true);
        RooExponential bbRandDE("bbRandDE", "bbRandDE", de, *bbRandDE_al) ;

	//dt
	RooRealVar *bbRandDT_tau = rws_fixed.var("bbRandDT_tau"); bbRandDT_tau->setConstant(true);
	RooRealVar *bbRandDT_mu_l = rws_fixed.var("bbRandDT_mu_l"); bbRandDT_mu_l->setConstant(true);
        RooRealVar *bbRandDT_mu_d = rws_fixed.var("bbRandDT_mu_d"); bbRandDT_mu_d->setConstant(true);
        RooRealVar *bbRandDT_f_delt = rws_fixed.var("bbRandDT_f_delt"); bbRandDT_f_delt->setConstant(true);
        RooRealVar *bbRandDT_f_tail = rws_fixed.var("bbRandDT_f_tail"); bbRandDT_f_tail->setConstant(true);
        RooRealVar *bbRandDT_S_main = rws_fixed.var("bbRandDT_S_main"); bbRandDT_S_main->setConstant(true);
        RooRealVar *bbRandDT_S_tail = rws_fixed.var("bbRandDT_S_tail"); bbRandDT_S_tail->setConstant(true);
	

        RooDtBkg bbRandDT("bbRandDT", "bbRandDT", dt,
                *bbRandDT_tau, *bbRandDT_mu_l, *bbRandDT_mu_d, *bbRandDT_f_delt, *bbRandDT_f_tail, *bbRandDT_S_main, *bbRandDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);

//------------------------------------------------------------------------------------------------------

	//Continuum
	//MBC
	RooRealVar *contMBC_c = rws_fixed.var("contMBC_c"); contMBC_c->setConstant(true);
	RooRealVar *contMBC_m0 = rws_fixed.var("contMBC_m0"); contMBC_m0->setConstant(true);
        //RooArgusBG  contMBC("contMBC","contMBC", mbc, *contMBC_m0, *contMBC_c);
	RooArgusBG  contMBC("contMBC","contMBC", mbc, commonArgus_m0, *contMBC_c);

	//dE
	RooRealVar *contDE_al = rws_fixed.var("contDE_al"); contDE_al->setConstant(true);
        RooExponential contDE("contDE","contDE", de, *contDE_al);

	//dt
	RooRealVar *contDT_tau = rws_fixed.var("contDT_tau"); contDT_tau->setConstant(true);
        RooRealVar *contDT_mu_l = rws_fixed.var("contDT_mu_l"); contDT_mu_l->setConstant(true);
        RooRealVar *contDT_mu_d = rws_fixed.var("contDT_mu_d"); contDT_mu_d->setConstant(true);
        RooRealVar *contDT_f_delt = rws_fixed.var("contDT_f_delt"); contDT_f_delt->setConstant(true);
        RooRealVar *contDT_f_tail = rws_fixed.var("contDT_f_tail"); contDT_f_tail->setConstant(true);
        RooRealVar *contDT_S_main = rws_fixed.var("contDT_S_main"); contDT_S_main->setConstant(true);
        RooRealVar *contDT_S_tail = rws_fixed.var("contDT_S_tail"); contDT_S_tail->setConstant(true);


        RooDtBkg contDT("contDT", "contDT", dt,
                *contDT_tau, *contDT_mu_l, *contDT_mu_d, *contDT_f_delt, *contDT_f_tail, *contDT_S_main, *contDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);

//------------------------------------------------------------------------------------------------------

	//BB Missing FSP-bkg
	//MBC
	//CrystalBall
	RooRealVar *bb_missFSP_MBC_mean = rws_fixed.var("bb_missFSP_MBC_mean"); bb_missFSP_MBC_mean->setConstant(true);
        RooRealVar *bb_missFSP_MBC_sigma = rws_fixed.var("bb_missFSP_MBC_sigma"); bb_missFSP_MBC_sigma->setConstant(true);
        RooRealVar *bb_missFSP_MBC_alpha = rws_fixed.var("bb_missFSP_MBC_alpha"); bb_missFSP_MBC_alpha->setConstant(true);
        RooRealVar *bb_missFSP_MBC_n = rws_fixed.var("bb_missFSP_MBC_n"); bb_missFSP_MBC_n->setConstant(true);
	RooCBShape bb_missFSP_MBC("bb_missFSP_MBC", "bb_missFSP_MBC", mbc, *bb_missFSP_MBC_mean, *bb_missFSP_MBC_sigma, *bb_missFSP_MBC_alpha, *bb_missFSP_MBC_n);

	//dE
	//Exponential
	RooRealVar *bb_missFSP_DE_alpha = rws_fixed.var("bb_missFSP_DE_alpha"); bb_missFSP_DE_alpha->setConstant(true);
	RooExponential bb_missFSP_DE_expo("bb_missFSP_DE_expo", "bb_missFSP_DE_expo", de, *bb_missFSP_DE_alpha);
        //Gaussian
        RooRealVar *bb_missFSP_DE_mean = rws_fixed.var("bb_missFSP_DE_mean"); bb_missFSP_DE_mean->setConstant(true);
        RooRealVar *bb_missFSP_DE_sigma = rws_fixed.var("bb_missFSP_DE_sigma"); bb_missFSP_DE_sigma->setConstant(true);
	RooGaussian bb_missFSP_DE_gauss("bb_missFSP_DE_gauss", "bb_missFSP_DE_gauss", de, *bb_missFSP_DE_mean, *bb_missFSP_DE_sigma);
        //Combined
        RooRealVar *bb_missFSP_DE_f = rws_fixed.var("bb_missFSP_DE_f"); bb_missFSP_DE_f->setConstant(true);
	RooAddPdf bb_missFSP_DE("bb_missFSP_DE", "bb_missFSP_DE", RooArgList(bb_missFSP_DE_expo, bb_missFSP_DE_gauss), RooArgList(*bb_missFSP_DE_f)) ;

        //dt
        /*
        RooRealVar *bb_missFSP_DT_tau = rws_fixed.var("bb_missFSP_DT_tau"); bb_missFSP_DT_tau->setConstant(true);
        RooRealVar *bb_missFSP_DT_mu_l = rws_fixed.var("bb_missFSP_DT_mu_l"); bb_missFSP_DT_mu_l->setConstant(true);
        RooRealVar *bb_missFSP_DT_mu_d = rws_fixed.var("bb_missFSP_DT_mu_d"); bb_missFSP_DT_mu_d->setConstant(true);
        RooRealVar *bb_missFSP_DT_f_delt = rws_fixed.var("bb_missFSP_DT_f_delt"); bb_missFSP_DT_f_delt->setConstant(true);
        RooRealVar *bb_missFSP_DT_f_tail = rws_fixed.var("bb_missFSP_DT_f_tail"); bb_missFSP_DT_f_tail->setConstant(true);
        RooRealVar *bb_missFSP_DT_S_main = rws_fixed.var("bb_missFSP_DT_S_main"); bb_missFSP_DT_S_main->setConstant(true);
        RooRealVar *bb_missFSP_DT_S_tail = rws_fixed.var("bb_missFSP_DT_S_tail"); bb_missFSP_DT_S_tail->setConstant(true);

        RooDtBkg bb_missFSP_DT("bb_missFSP_DT", "bb_missFSP_DT", dt,
                *bb_missFSP_DT_tau, *bb_missFSP_DT_mu_l, *bb_missFSP_DT_mu_d, *bb_missFSP_DT_f_delt, *bb_missFSP_DT_f_tail, *bb_missFSP_DT_S_main, *bb_missFSP_DT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	*/
	const double tau_bb_missFSP = 1.4896;
        RooDtCPSignal bb_missFSP_DT("bb_missFSP_DT", "bb_missFSP_DT", dt, RooConst(0.0), RooConst(0.0), inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_bb_missFSP, true, true, 1.0);

//------------------------------------------------------------------------------------------------------


	//Self cross feed
        //MBC
	//Gaussian
	RooRealVar *scfMBC_mean = rws_fixed.var("scfMBC_mean"); scfMBC_mean->setConstant(true);

        RooRealVar *scfMBC_sigmaL = rws_fixed.var("scfMBC_sigmaL"); scfMBC_sigmaL->setConstant(true);
	RooRealVar *scfMBC_sigmaR = rws_fixed.var("scfMBC_sigmaR"); scfMBC_sigmaR->setConstant(true);
	RooBifurGauss scfMBC_Gauss("scfMBC_Gauss", "scfMBC_Gauss", mbc, *scfMBC_mean, *scfMBC_sigmaL, *scfMBC_sigmaR) ;

        //Argus
        RooRealVar *scfMBC_c = rws_fixed.var("scfMBC_c"); scfMBC_c->setConstant(true);
        RooRealVar *scfMBC_m0 = rws_fixed.var("scfMBC_m0"); scfMBC_m0->setConstant(true);
        //RooArgusBG scfMBC_Argus("scfMBC_Argus", "scfMBC_Argus", mbc, *scfMBC_m0, *scfMBC_c) ;
	RooArgusBG scfMBC_Argus("scfMBC_Argus", "scfMBC_Argus", mbc, commonArgus_m0, *scfMBC_c) ;
        //combined
        RooRealVar *scfMBC_f = rws_fixed.var("scfMBC_f"); scfMBC_f->setConstant(true);
        RooAddPdf scfMBC("scfMBC", "scfMBC", RooArgList(scfMBC_Gauss, scfMBC_Argus), RooArgList(*scfMBC_f)) ;


        //dE
	//Chebychev
	RooRealVar *scfDE_p0 = rws_fixed.var("scfDE_p0"); scfDE_p0->setConstant(true);
	RooRealVar *scfDE_p1 = rws_fixed.var("scfDE_p1"); scfDE_p1->setConstant(true);
	RooChebychev scfDE("scfDE", "scfDE", de, RooArgList(*scfDE_p0, *scfDE_p1));

	//dt
	/*
	RooRealVar *scfDT_tau = rws_fixed.var("scfDT_tau"); scfDT_tau->setConstant(true);
	RooRealVar *scfDT_mu_l = rws_fixed.var("scfDT_mu_l"); scfDT_mu_l->setConstant(true);
        RooRealVar *scfDT_mu_d = rws_fixed.var("scfDT_mu_d"); scfDT_mu_d->setConstant(true);
        RooRealVar *scfDT_f_delt = rws_fixed.var("scfDT_f_delt"); scfDT_f_delt->setConstant(true);
        RooRealVar *scfDT_f_tail = rws_fixed.var("scfDT_f_tail"); scfDT_f_tail->setConstant(true);
        RooRealVar *scfDT_S_main = rws_fixed.var("scfDT_S_main"); scfDT_S_main->setConstant(true);
        RooRealVar *scfDT_S_tail = rws_fixed.var("scfDT_S_tail"); scfDT_S_tail->setConstant(true);

        RooDtBkg scfDT("scfDT", "scfDT", dt,
                *scfDT_tau, *scfDT_mu_l, *scfDT_mu_d, *scfDT_f_delt, *scfDT_f_tail, *scfDT_S_main, *scfDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	*/
	//SCF-pars
	const double tau_scf = 1.0906;
        RooRealVar A_scf("A_scf", "A_scf", 0.0, -2., 2.) ;
        RooRealVar S_scf("S_scf", "S_scf", 0.0, -2., 2.) ;
        //Signal-pars
        RooRealVar A("A", "A", -2., 2.) ;
        RooRealVar S("S", "S",  -2., 2.);

	RooDtCPSignal scfDT("scfDT", "scfDT", dt, RooConst(0.0), RooConst(0.0), inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);

        //RooDtCPSignal scfDT("scfDT", "scfDT", dt, S, A, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);

//------------------------------------------------------------------------------------------------------
//

        //Signal
        //MBC
        RooRealVar meanMBC_fudge("meanMBC_fudge", "meanMBC_fudge", 0.0005, -0.01, 0.01);
        RooRealVar sigmaMBC_fudge("sigmaMBC_fudge", "sigmaMBC_fudge", 1.0005, 0.5, 2.0);
        RooFormulaVar signalMBC_mean("signalMBC_mean", "5.2796+meanMBC_fudge", RooArgList(meanMBC_fudge));
        RooFormulaVar signalMBC_sigma("signalMBC_sigma", "0.0029357*sigmaMBC_fudge", RooArgList(sigmaMBC_fudge));
        RooRealVar *signalMBC_alpha = rws_fixed.var("signalMBC_alpha"); signalMBC_alpha->setConstant(true);
        RooRealVar *signalMBC_n = rws_fixed.var("signalMBC_n"); signalMBC_n->setConstant(true);
        RooCBShape signalMBC("signalMBC", "signalMBC", mbc, signalMBC_mean, signalMBC_sigma, *signalMBC_alpha, *signalMBC_n);

	//dE
        RooRealVar meanDE_fudge("meanDE_fudge", "meanDE_fudge", 0.0005, -0.01, 0.01);
        RooRealVar sigmaDE_fudge("sigmaDE_fudge", "sigmaDE_fudge", 1.0005,0.5,2.0);
        RooFormulaVar signalDE_mean("signalDE_mean", "-0.00689755+meanDE_fudge", RooArgList(meanDE_fudge));
        RooFormulaVar signalDE_sigma("signalDE_sigma", "0.037056*sigmaDE_fudge", RooArgList(sigmaDE_fudge));	
	RooRealVar *signalDE_alpha = rws_fixed.var("signalDE_alpha"); signalDE_alpha->setConstant(true);
	RooRealVar *signalDE_n = rws_fixed.var("signalDE_n"); signalDE_n->setConstant(true);
        RooCBShape signalDE("signalDE", "signalDE", de, signalDE_mean, signalDE_sigma, *signalDE_alpha, *signalDE_n);

        //dt
	//RooRealVar A("A", "A", -2., 2.) ;
        //RooRealVar S("S", "S", -2., 2.) ;
        RooDtCPSignal signalDT("signalDT", "signalDT", dt, S, A,
                        inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig,
                        dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag,
                        keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);



//------------------------------------------------------------------------------------------------------

 

	//Here we load the randomized parameters
	TFile *F_bbmissFSP = new TFile("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/RandomizedParameters/fix_shape_new/bb_missFSP_randPars.root");
	RooWorkspace *w_bbmissFSP = (RooWorkspace *)F_bbmissFSP->Get("rws_bb_missFSP_randPars");
	TFile *F_cont = new TFile("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/RandomizedParameters/fix_shape_new/cont_randPars.root");
	RooWorkspace *w_cont  = (RooWorkspace *)F_cont->Get("rws_cont_randPars");
	TFile *F_bbrand = new TFile("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/RandomizedParameters/fix_shape_new/bb_randPars.root");
	RooWorkspace *w_bbrand = (RooWorkspace *)F_bbrand->Get("rws_bb_randPars");
	TFile *F_scf = new TFile("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/RandomizedParameters/fix_shape_new/scf_randPars.root");
	RooWorkspace *w_scf = (RooWorkspace *)F_scf->Get("rws_scf_randPars");
	TFile *F_signal = new TFile("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/RandomizedParameters/fix_shape_new/signal_pure_randPars.root");
	RooWorkspace *w_signal = (RooWorkspace *)F_signal->Get("rws_signal_randPars");	

	//BB-random
	if(component == "bbrand" || component == "all"){
	bbRandMBC_c->setVal(randParVal(w_bbrand, "mbc", toyindex, "bbRandMBC_c")); 
	bbRandMBC_m0->setVal(randParVal(w_bbrand, "mbc", toyindex, "bbRandMBC_m0")); 
	bbRandDE_al->setVal(randParVal(w_bbrand, "de", toyindex, "bbRandDE_al"));
	bbRandDT_tau->setVal(randParVal(w_bbrand, "dt", toyindex, "bbRandDT_tau"));
	bbRandDT_mu_l->setVal(randParVal(w_bbrand, "dt", toyindex, "bbRandDT_mu_l"));
	bbRandDT_mu_d->setVal(randParVal(w_bbrand, "dt", toyindex, "bbRandDT_mu_d"));
	bbRandDT_f_delt->setVal(randParVal(w_bbrand, "dt", toyindex, "bbRandDT_f_delt"));
	bbRandDT_f_tail->setVal(randParVal(w_bbrand, "dt", toyindex, "bbRandDT_f_tail"));
	bbRandDT_S_main->setVal(randParVal(w_bbrand, "dt", toyindex, "bbRandDT_S_main"));
	bbRandDT_S_tail->setVal(randParVal(w_bbrand, "dt", toyindex, "bbRandDT_S_tail"));
	}

	//Continuum
	if(component == "cont" || component == "all"){
	contMBC_c->setVal(randParVal(w_cont, "mbc", toyindex, "contMBC_c"));
	contMBC_m0->setVal(randParVal(w_cont, "mbc", toyindex, "contMBC_m0"));
	contDE_al->setVal(randParVal(w_cont, "de", toyindex, "contDE_al"));
	contDT_tau->setVal(randParVal(w_cont, "dt", toyindex, "contDT_tau"));
	contDT_mu_l->setVal(randParVal(w_cont, "dt", toyindex, "contDT_mu_l"));
	contDT_mu_d->setVal(randParVal(w_cont, "dt", toyindex, "contDT_mu_d"));
	contDT_f_delt->setVal(randParVal(w_cont, "dt", toyindex, "contDT_f_delt"));
	contDT_f_tail->setVal(randParVal(w_cont, "dt", toyindex, "contDT_f_tail"));
	contDT_S_main->setVal(randParVal(w_cont, "dt", toyindex, "contDT_S_main"));
	contDT_S_tail->setVal(randParVal(w_cont, "dt", toyindex, "contDT_S_tail"));
	}

	//BB-missFSP
	if(component == "bbmissfsp" || component == "all"){
	bb_missFSP_MBC_mean->setVal(randParVal(w_bbmissFSP, "mbc", toyindex, "bb_missFSP_MBC_mean"));
	bb_missFSP_MBC_sigma->setVal(randParVal(w_bbmissFSP, "mbc", toyindex, "bb_missFSP_MBC_sigma"));
	bb_missFSP_MBC_alpha->setVal(randParVal(w_bbmissFSP, "mbc", toyindex, "bb_missFSP_MBC_alpha"));
	bb_missFSP_MBC_n->setVal(randParVal(w_bbmissFSP, "mbc", toyindex, "bb_missFSP_MBC_n"));
	bb_missFSP_DE_alpha->setVal(randParVal(w_bbmissFSP, "de", toyindex, "bb_missFSP_DE_alpha"));
	bb_missFSP_DE_mean->setVal(randParVal(w_bbmissFSP, "de", toyindex, "bb_missFSP_DE_mean"));
	bb_missFSP_DE_sigma->setVal(randParVal(w_bbmissFSP, "de", toyindex, "bb_missFSP_DE_sigma"));
	bb_missFSP_DE_f->setVal(randParVal(w_bbmissFSP, "de", toyindex, "bb_missFSP_DE_f"));
	//bb_missFSP_DT_tau->setVal(randParVal(w_bbmissFSP, "dt", toyindex, "bb_missFSP_DT_tau"));
	//bb_missFSP_DT_mu_l->setVal(randParVal(w_bbmissFSP, "dt", toyindex, "bb_missFSP_DT_mu_l"));
	//bb_missFSP_DT_mu_d->setVal(randParVal(w_bbmissFSP, "dt", toyindex, "bb_missFSP_DT_mu_d"));
	//bb_missFSP_DT_f_delt->setVal(randParVal(w_bbmissFSP, "dt", toyindex, "bb_missFSP_DT_f_delt"));	
  	//bb_missFSP_DT_f_tail->setVal(randParVal(w_bbmissFSP, "dt", toyindex, "bb_missFSP_DT_f_tail"));
	//bb_missFSP_DT_S_main->setVal(randParVal(w_bbmissFSP, "dt", toyindex, "bb_missFSP_DT_S_main"));
	//bb_missFSP_DT_S_tail->setVal(randParVal(w_bbmissFSP, "dt", toyindex, "bb_missFSP_DT_S_tail"));
	}

	//SCF
	if(component == "scf" || component == "all"){
	scfMBC_mean->setVal(randParVal(w_scf, "mbc", toyindex, "scfMBC_mean"));
	scfMBC_sigmaL->setVal(randParVal(w_scf, "mbc", toyindex, "scfMBC_sigmaL"));
	scfMBC_sigmaR->setVal(randParVal(w_scf, "mbc", toyindex, "scfMBC_sigmaR"));
	scfMBC_c->setVal(randParVal(w_scf, "mbc", toyindex, "scfMBC_c"));
	scfMBC_m0->setVal(randParVal(w_scf, "mbc", toyindex, "scfMBC_m0"));
	scfMBC_f->setVal(randParVal(w_scf, "mbc", toyindex, "scfMBC_f"));
	scfDE_p0->setVal(randParVal(w_scf, "de", toyindex, "scfDE_p0"));
	scfDE_p1->setVal(randParVal(w_scf, "de", toyindex, "scfDE_p1"));
	//scfDT_tau->setVal(randParVal(w_scf, "dt", toyindex, "scfDT_tau"));
	//scfDT_mu_l->setVal(randParVal(w_scf, "dt", toyindex, "scfDT_mu_l"));
	//scfDT_mu_d->setVal(randParVal(w_scf, "dt", toyindex, "scfDT_mu_d"));
	//scfDT_f_delt->setVal(randParVal(w_scf, "dt", toyindex, "scfDT_f_delt"));
	//scfDT_f_tail->setVal(randParVal(w_scf, "dt", toyindex, "scfDT_f_tail"));
	//scfDT_S_main->setVal(randParVal(w_scf, "dt", toyindex, "scfDT_S_main"));
	//scfDT_S_tail->setVal(randParVal(w_scf, "dt", toyindex, "scfDT_S_tail"));
	}

	//Signal
	if(component == "signal" || component == "all"){
	signalMBC_alpha->setVal(randParVal(w_signal, "mbc", toyindex, "signalMBC_alpha"));
        signalMBC_n->setVal(randParVal(w_signal, "mbc", toyindex, "signalMBC_n"));
        signalDE_alpha->setVal(randParVal(w_signal, "de", toyindex, "signalDE_alpha"));
        signalDE_n->setVal(randParVal(w_signal, "de", toyindex, "signalDE_n"));
	}
	
	cout<<"stage after random parameter fixing"<<endl;
//------------------------------------------------------------------------------------------------------	


 	//3D pdf's
	RooProdPdf * cont3D = new RooProdPdf("cont3D", "cont3D", RooArgList(contMBC, contDE, contDT));
        RooProdPdf * bbRand3D = new RooProdPdf("bbRand3D", "bbRand3D", RooArgList(bbRandMBC, bbRandDE, bbRandDT));
        RooProdPdf * bb_missFSP3D = new RooProdPdf("bb_missFSP3D", "bb_missFSP3D", RooArgList(bb_missFSP_MBC, bb_missFSP_DE, bb_missFSP_DT));
        RooProdPdf * scf3D = new RooProdPdf("scf3D", "scf3D", RooArgList(scfMBC, scfDE, scfDT));
        RooProdPdf * signal3D =  new RooProdPdf("signal3D","signal3D", RooArgSet(signalMBC, signalDE, signalDT));


	//OVERALL NSIG and NBKG floated
	//Excluding rbin-0
	//RooRealVar f_scf("f_scf","scf fraction", 0.1394);//, 0.0, 0.6);//fixed from MC truth info	
	//RooRealVar f_bbRand("f_bbRand","fraction of rare bkg", 0.0696);//fixed from MC truth info
	//RooRealVar f_bb_missFSP("f_bb_missFSP","fraction of bb bkg", 0.0077);//fixed from MC truth info
	//Including rbin-0
	RooRealVar f_scf("f_scf","scf fraction", 0.1416);//, 0.0, 0.6);//fixed from MC truth info       
        RooRealVar f_bbRand("f_bbRand","fraction of rare bkg", 0.0602);//fixed from MC truth info
        RooRealVar f_bb_missFSP("f_bb_missFSP","fraction of bb bkg", 0.0066);//fixed from MC truth info

	RooRealVar Nsigscf("N_{sig+scf}","Signal + SCF yield", 840, 300, 1500);
	RooRealVar Nbkgtot("N_{bkg}","Combined background yield",4650, 1500, 6000);
	


	//Combining P.D.F.'s
	RooAddPdf * fullSig3D = new RooAddPdf("fullSig3D","fullSig3D", RooArgList(*scf3D, *signal3D), RooArgList(f_scf));
	RooAddPdf * fullBkg3D = new RooAddPdf("fullBkg3D","fullBkg3D", RooArgList(*bb_missFSP3D, *bbRand3D, *cont3D), RooArgList(f_bb_missFSP, f_bbRand));
	RooAddPdf * combined3D = new RooAddPdf("combined3D","combined3D", RooArgList(*fullSig3D, *fullBkg3D), RooArgList(Nsigscf, Nbkgtot));

	


//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult * fullFit = combined3D->fitTo(*data, ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;

 





/*	

	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40), Title(""));
	data->plotOn(xframe_mbc);
	combined3D->plotOn(xframe_mbc,Components("signalMBC"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
	combined3D->plotOn(xframe_mbc, LineColor(kBlue), ProjWData(*data));
	//combined3D->paramOn(xframe_mbc, ProjWData(*data));	
	RooHist* pull_mbc = xframe_mbc->pullHist() ;
	RooPlot* f_pull_mbc = mbc.frame(Title("")) ;
	f_pull_mbc->addPlotable(pull_mbc,"BX0") ;
	f_pull_mbc->SetMaximum(4);
	f_pull_mbc->SetMinimum(-4);
	f_pull_mbc->GetYaxis()->SetLabelSize(0.11);
        f_pull_mbc->GetXaxis()->SetLabelSize(0.11);


	//dE
	RooPlot *xframe_de = de.frame(Bins(40), Title(""));
	data->plotOn(xframe_de);
	combined3D->plotOn(xframe_de, Components("signalDE"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
	combined3D->plotOn(xframe_de, LineColor(kBlue), ProjWData(*data));
	//combined3D->paramOn(xframe_de, ProjWData(*data));
	RooHist* pull_de = xframe_de->pullHist() ;
        RooPlot* f_pull_de = de.frame(Title("")) ;
        f_pull_de->addPlotable(pull_de,"BX0") ;
        f_pull_de->SetMaximum(4);
        f_pull_de->SetMinimum(-4);
	f_pull_de->GetYaxis()->SetLabelSize(0.11);
        f_pull_de->GetXaxis()->SetLabelSize(0.11);	

	//dt
	RooPlot *xframe_dt = dt.frame(Bins(40), Title(""));
        data->plotOn(xframe_dt);
	combined3D->plotOn(xframe_dt,Components("signalDT"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
        combined3D->plotOn(xframe_dt, LineColor(kBlue), ProjWData(*data));
        combined3D->paramOn(xframe_dt, ProjWData(*data), Layout(0.60));
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
        pad31->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_dt->Draw() ;// xframe_dt->Draw() ;
        pad32->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_dt->Draw() ;
	can->SaveAs(canvasName.c_str());
*/	

	//Save fit results
	
	//To individual fit file
	RooArgSet(fullFit->floatParsFinal()).writeToFile(resultName.c_str());

	//To consolidated parameter result files
	ofstream fout; string confilename, fileSuffix;
	fileSuffix = "_" + component + ".txt";
	//A
	confilename = "/A"; confilename = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
	fout<<toyindex<<"\t"<<A.getVal()<<"\t"<<A.getError()<<"\t"<<fullFit->status()<<endl ;
	fout.close();
	//S
	confilename = "/S"; confilename = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
	fout<<toyindex<<"\t"<<S.getVal()<<"\t"<<S.getError()<<"\t"<<fullFit->status()<<endl ;
	fout.close();

	//Nsigscf
	confilename = "/Nsigscf"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsigscf.getVal()<<"\t"<<Nsigscf.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();
	
	//Nbkgtot
	confilename = "/Nbkgtot"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nbkgtot.getVal()<<"\t"<<Nbkgtot.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();

	

//----------------------------------------------------------------------------------------------//		

        return 0;
}

