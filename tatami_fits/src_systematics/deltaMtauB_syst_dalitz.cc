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
#include "RooSimultaneous.h"
#include "RooChebychev.h"
#include "RooBifurGauss.h"

//#include "RooCFunction3Binding.h"


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

	//const double delta_m = 0.507;
	//const double tau_b = 1.5344;	
	//const double tau_b = 1.525;
	//Randomized Delta-m and tauB
	TFile *F_phypars = new TFile("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way/RandomizedParameters/deltaMtauB/randPhyPars.root");
        RooWorkspace *w_phypars = (RooWorkspace *)F_phypars->Get("rws_randPhyPars");
        const double delta_m = randParVal(w_phypars, "data", toyindex, "DeltaM");
        const double tau_b = randParVal(w_phypars, "data", toyindex, "tauB");

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
	//Naming convention : Dp for daltiz+ and Dm for Dalitz-
	// Dalitz region I-bar (S.Akar) corresponds to my Dalitz-Plus category
	// and Dalitz region I corresponds to my Dalitz-Minus category	
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
        //RooArgusBG bbMBC_argus("bbMBC_argus","bbMBC_argus", mbc, *bbMBC_m0, *bbMBC_c);
	RooArgusBG bbRandMBC_Dp("bbRandMBC_Dp","bbRandMBC_Dp", mbc, commonArgus_m0, *bbRandMBC_c);
        RooArgusBG bbRandMBC_Dm("bbRandMBC_Dm","bbRandMBC_Dm", mbc, commonArgus_m0, *bbRandMBC_c);

	

	//dE
        RooRealVar *bbRandDE_al = rws_fixed.var("bbRandDE_al"); bbRandDE_al->setConstant(true);
        RooExponential bbRandDE_Dp("bbRandDE_Dp", "bbRandDE_Dp", de, *bbRandDE_al);
	RooExponential bbRandDE_Dm("bbRandDE_Dm", "bbRandDE_Dm", de, *bbRandDE_al);

	//dt
	RooRealVar *bbRandDT_tau = rws_fixed.var("bbRandDT_tau"); bbRandDT_tau->setConstant(true);
        RooRealVar *bbRandDT_mu_l = rws_fixed.var("bbRandDT_mu_l"); bbRandDT_mu_l->setConstant(true);
        RooRealVar *bbRandDT_mu_d = rws_fixed.var("bbRandDT_mu_d"); bbRandDT_mu_d->setConstant(true);
        RooRealVar *bbRandDT_f_delt = rws_fixed.var("bbRandDT_f_delt"); bbRandDT_f_delt->setConstant(true);
        RooRealVar *bbRandDT_f_tail = rws_fixed.var("bbRandDT_f_tail"); bbRandDT_f_tail->setConstant(true);
        RooRealVar *bbRandDT_S_main = rws_fixed.var("bbRandDT_S_main"); bbRandDT_S_main->setConstant(true);
        RooRealVar *bbRandDT_S_tail = rws_fixed.var("bbRandDT_S_tail"); bbRandDT_S_tail->setConstant(true);


        RooDtBkg bbRandDT_Dp("bbRandDT_Dp", "bbRandDT_Dp", dt, *bbRandDT_tau, *bbRandDT_mu_l, *bbRandDT_mu_d, *bbRandDT_f_delt, *bbRandDT_f_tail, *bbRandDT_S_main, *bbRandDT_S_tail, inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	RooDtBkg bbRandDT_Dm("bbRandDT_Dm", "bbRandDT_Dm", dt, *bbRandDT_tau, *bbRandDT_mu_l, *bbRandDT_mu_d, *bbRandDT_f_delt, *bbRandDT_f_tail, *bbRandDT_S_main, *bbRandDT_S_tail, inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);

//------------------------------------------------------------------------------------------------------
	//Continuum
	//MBC
	RooRealVar *contMBC_c = rws_fixed.var("contMBC_c"); contMBC_c->setConstant(true);
        RooRealVar *contMBC_m0 = rws_fixed.var("contMBC_m0"); contMBC_m0->setConstant(true);
        //RooArgusBG  contMBC("contMBC","contMBC", mbc, *contMBC_m0, *contMBC_c);
	RooArgusBG  contMBC_Dp("contMBC_Dp","contMBC_Dp", mbc, commonArgus_m0, *contMBC_c);
	RooArgusBG  contMBC_Dm("contMBC_Dm","contMBC_Dm", mbc, commonArgus_m0, *contMBC_c);

	//dE
	RooRealVar *contDE_al = rws_fixed.var("contDE_al"); contDE_al->setConstant(true);
        RooExponential contDE_Dp("contDE_Dp","contDE_Dp", de, *contDE_al);
	RooExponential contDE_Dm("contDE_Dm","contDE_Dm", de, *contDE_al);

	//dt
	RooRealVar *contDT_tau = rws_fixed.var("contDT_tau"); contDT_tau->setConstant(true);
        RooRealVar *contDT_mu_l = rws_fixed.var("contDT_mu_l"); contDT_mu_l->setConstant(true);
        RooRealVar *contDT_mu_d = rws_fixed.var("contDT_mu_d"); contDT_mu_d->setConstant(true);
        RooRealVar *contDT_f_delt = rws_fixed.var("contDT_f_delt"); contDT_f_delt->setConstant(true);
        RooRealVar *contDT_f_tail = rws_fixed.var("contDT_f_tail"); contDT_f_tail->setConstant(true);
        RooRealVar *contDT_S_main = rws_fixed.var("contDT_S_main"); contDT_S_main->setConstant(true);
        RooRealVar *contDT_S_tail = rws_fixed.var("contDT_S_tail"); contDT_S_tail->setConstant(true);

        RooDtBkg contDT_Dp("contDT_Dp", "contDT_Dp", dt,
                *contDT_tau, *contDT_mu_l, *contDT_mu_d, *contDT_f_delt, *contDT_f_tail, *contDT_S_main, *contDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	RooDtBkg contDT_Dm("contDT_Dm", "contDT_Dm", dt,
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
        RooCBShape bb_missFSP_MBC_Dp("bb_missFSP_MBC_Dp", "bb_missFSP_MBC_Dp", mbc, *bb_missFSP_MBC_mean, *bb_missFSP_MBC_sigma, *bb_missFSP_MBC_alpha, *bb_missFSP_MBC_n);
	RooCBShape bb_missFSP_MBC_Dm("bb_missFSP_MBC_Dm", "bb_missFSP_MBC_Dm", mbc, *bb_missFSP_MBC_mean, *bb_missFSP_MBC_sigma, *bb_missFSP_MBC_alpha, *bb_missFSP_MBC_n);

	//dE
	//Exponential
	RooRealVar *bb_missFSP_DE_alpha = rws_fixed.var("bb_missFSP_DE_alpha"); bb_missFSP_DE_alpha->setConstant(true);
        RooExponential bb_missFSP_DE_expo_Dp("bb_missFSP_DE_expo_Dp", "bb_missFSP_DE_expo_Dp", de, *bb_missFSP_DE_alpha);
	RooExponential bb_missFSP_DE_expo_Dm("bb_missFSP_DE_expo_Dm", "bb_missFSP_DE_expo_Dm", de, *bb_missFSP_DE_alpha);
        //Gaussian
        RooRealVar *bb_missFSP_DE_mean = rws_fixed.var("bb_missFSP_DE_mean"); bb_missFSP_DE_mean->setConstant(true);
        RooRealVar *bb_missFSP_DE_sigma = rws_fixed.var("bb_missFSP_DE_sigma"); bb_missFSP_DE_sigma->setConstant(true);
        RooGaussian bb_missFSP_DE_gauss_Dp("bb_missFSP_DE_gauss_Dp", "bb_missFSP_DE_gauss_Dp", de, *bb_missFSP_DE_mean, *bb_missFSP_DE_sigma);
	RooGaussian bb_missFSP_DE_gauss_Dm("bb_missFSP_DE_gauss_Dm", "bb_missFSP_DE_gauss_Dm", de, *bb_missFSP_DE_mean, *bb_missFSP_DE_sigma);
        //Combined
        RooRealVar *bb_missFSP_DE_f = rws_fixed.var("bb_missFSP_DE_f"); bb_missFSP_DE_f->setConstant(true);
	RooAddPdf bb_missFSP_DE_Dp("bb_missFSP_DE_Dp", "bb_missFSP_DE_Dp", RooArgList(bb_missFSP_DE_expo_Dp, bb_missFSP_DE_gauss_Dp), RooArgList(*bb_missFSP_DE_f)) ;
	RooAddPdf bb_missFSP_DE_Dm("bb_missFSP_DE_Dm", "bb_missFSP_DE_Dm", RooArgList(bb_missFSP_DE_expo_Dm, bb_missFSP_DE_gauss_Dm), RooArgList(*bb_missFSP_DE_f)) ;

	
	//dt
	const double tau_bb_missFSP = 1.4896;
        RooDtCPSignal bb_missFSP_DT_Dp("bb_missFSP_DT_Dp", "bb_missFSP_DT_Dp", dt, RooConst(0.0), RooConst(0.0), inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_bb_missFSP, true, true, 1.0);
	RooDtCPSignal bb_missFSP_DT_Dm("bb_missFSP_DT_Dm", "bb_missFSP_DT_Dm", dt, RooConst(0.0), RooConst(0.0), inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_bb_missFSP, true, true, 1.0);

//------------------------------------------------------------------------------------------------------
	//Self cross feed
        //MBC
	//Gaussian
	RooRealVar *scfMBC_mean = rws_fixed.var("scfMBC_mean"); scfMBC_mean->setConstant(true);
        RooRealVar *scfMBC_sigmaL = rws_fixed.var("scfMBC_sigmaL"); scfMBC_sigmaL->setConstant(true);
	RooRealVar *scfMBC_sigmaR = rws_fixed.var("scfMBC_sigmaR"); scfMBC_sigmaR->setConstant(true);
        RooBifurGauss scfMBC_Gauss_Dp("scfMBC_Gauss_Dp", "scfMBC_Gauss_Dp", mbc, *scfMBC_mean, *scfMBC_sigmaL, *scfMBC_sigmaR) ;
	RooBifurGauss scfMBC_Gauss_Dm("scfMBC_Gauss_Dm", "scfMBC_Gauss_Dm", mbc, *scfMBC_mean, *scfMBC_sigmaL, *scfMBC_sigmaR) ;
        //Argus
        RooRealVar *scfMBC_c = rws_fixed.var("scfMBC_c"); scfMBC_c->setConstant(true);
        RooRealVar *scfMBC_m0 = rws_fixed.var("scfMBC_m0"); scfMBC_m0->setConstant(true);
        //RooArgusBG scfMBC_Argus("scfMBC_Argus", "scfMBC_Argus", mbc, *scfMBC_m0, *scfMBC_c) ;
	RooArgusBG scfMBC_Argus_Dp("scfMBC_Argus_Dp", "scfMBC_Argus_Dp", mbc, commonArgus_m0, *scfMBC_c) ;
	RooArgusBG scfMBC_Argus_Dm("scfMBC_Argus_Dm", "scfMBC_Argus_Dm", mbc, commonArgus_m0, *scfMBC_c) ;

        //combined
        RooRealVar *scfMBC_f = rws_fixed.var("scfMBC_f"); scfMBC_f->setConstant(true);
        RooAddPdf scfMBC_Dp("scfMBC_Dp", "scfMBC_Dp", RooArgList(scfMBC_Gauss_Dp, scfMBC_Argus_Dp), RooArgList(*scfMBC_f)) ;
	RooAddPdf scfMBC_Dm("scfMBC_Dm", "scfMBC_Dm", RooArgList(scfMBC_Gauss_Dm, scfMBC_Argus_Dm), RooArgList(*scfMBC_f)) ;


        //dE
	//Exponential
	RooRealVar *scfDE_p0 = rws_fixed.var("scfDE_p0"); scfDE_p0->setConstant(true);
        RooRealVar *scfDE_p1 = rws_fixed.var("scfDE_p1"); scfDE_p1->setConstant(true);
        RooChebychev scfDE_Dp("scfDE_Dp", "scfDE_Dp", de, RooArgList(*scfDE_p0, *scfDE_p1));
	RooChebychev scfDE_Dm("scfDE_Dm", "scfDE_Dm", de, RooArgList(*scfDE_p0, *scfDE_p1));

	//dt
	//SCF-pars
	const double tau_scf = 1.0906;
        RooRealVar A_scf("A_scf", "A_scf", 0.0, -2., 2.) ;
        RooRealVar S_scf("S_scf", "S_scf", 0.0, -2., 2.) ;
	//Signal-pars
        RooRealVar A("A", "A", -2., 2.) ;
        RooRealVar S_plus("S_plus", "S_plus",  -2., 2.);
        RooRealVar S_minus("S_minus", "S_minus",  -2., 2.);
        RooFormulaVar S_I("S_I", "(S_plus+S_minus)/2", RooArgList(S_plus, S_minus));
        RooFormulaVar S_Ib("S_Ib", "(S_plus-S_minus)/2", RooArgList(S_plus, S_minus));
        //RooRealVar S_I("S_I", "S_I",  -2., 2.);
        //RooRealVar S_Ib("S_Ib", "S_Ib",  -2., 2.);



        //RooDtCPSignal scfDT_Dp("scfDT_Dp", "scfDT_Dp", dt, S_Ib, A, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);
	//RooDtCPSignal scfDT_Dm("scfDT_Dm", "scfDT_Dm", dt, S_I, A, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);

       RooDtCPSignal scfDT_Dp("scfDT_Dp", "scfDT_Dp", dt, S_scf, A_scf, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);
        RooDtCPSignal scfDT_Dm("scfDT_Dm", "scfDT_Dm", dt, S_scf, A_scf, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);


	//RooDtCPSignal scfDT_Dp("scfDT_Dp", "scfDT_Dp", dt, RooConst(0.0), RooConst(0.0), inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);
	//RooDtCPSignal scfDT_Dm("scfDT_Dm", "scfDT_Dm", dt, RooConst(0.0), RooConst(0.0), inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_scf, true, true, 1.0);

//------------------------------------------------------------------------------------------------------
	
        //Signal
        //MBC
        RooRealVar meanMBC_fudge("meanMBC_fudge", "meanMBC_fudge", 0.0005, -0.01, 0.01);
        RooRealVar sigmaMBC_fudge("sigmaMBC_fudge", "sigmaMBC_fudge", 1.0005, 0.5, 2.0);
        RooFormulaVar signalMBC_mean("signalMBC_mean", "5.2796+meanMBC_fudge", RooArgList(meanMBC_fudge));
        RooFormulaVar signalMBC_sigma("signalMBC_sigma", "0.0029357*sigmaMBC_fudge", RooArgList(sigmaMBC_fudge));
	RooRealVar *signalMBC_alpha = rws_fixed.var("signalMBC_alpha"); signalMBC_alpha->setConstant(true);
        RooRealVar *signalMBC_n = rws_fixed.var("signalMBC_n"); signalMBC_n->setConstant(true);
        RooCBShape signalMBC_Dp("signalMBC_Dp", "signalMBC_Dp", mbc, signalMBC_mean, signalMBC_sigma, *signalMBC_alpha, *signalMBC_n);
	RooCBShape signalMBC_Dm("signalMBC_Dm", "signalMBC_Dm", mbc, signalMBC_mean, signalMBC_sigma, *signalMBC_alpha, *signalMBC_n);

	//dE
        RooRealVar meanDE_fudge("meanDE_fudge", "meanDE_fudge", 0.0005, -0.01, 0.01);
        RooRealVar sigmaDE_fudge("sigmaDE_fudge", "sigmaDE_fudge", 1.0005,0.5,2.0);
        RooFormulaVar signalDE_mean("signalDE_mean", "-0.00689755+meanDE_fudge", RooArgList(meanDE_fudge));
        RooFormulaVar signalDE_sigma("signalDE_sigma", "0.037056*sigmaDE_fudge", RooArgList(sigmaDE_fudge));	
	RooRealVar *signalDE_alpha = rws_fixed.var("signalDE_alpha"); signalDE_alpha->setConstant(true);
        RooRealVar *signalDE_n = rws_fixed.var("signalDE_n"); signalDE_n->setConstant(true);
        RooCBShape signalDE_Dp("signalDE_Dp", "signalDE_Dp", de, signalDE_mean, signalDE_sigma, *signalDE_alpha, *signalDE_n);
	RooCBShape signalDE_Dm("signalDE_Dm", "signalDE_Dm", de, signalDE_mean, signalDE_sigma, *signalDE_alpha, *signalDE_n);

        //dt
        RooDtCPSignal signalDT_Dp("signalDT_Dp", "signalDT_Dp", dt, S_Ib, A, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);
	RooDtCPSignal signalDT_Dm("signalDT_Dm", "signalDT_Dm", dt, S_I, A, inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig, dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag, keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);

//------------------------------------------------------------------------------------------------------	



	//3D pdf's
	//Dalitz+
	RooProdPdf * cont3D_Dp = new RooProdPdf("cont3D_Dp", "cont3D_Dp", RooArgList(contMBC_Dp, contDE_Dp, contDT_Dp));
	RooProdPdf * bbRand3D_Dp = new RooProdPdf("bbRand3D_Dp", "bbRand3D_Dp", RooArgList(bbRandMBC_Dp, bbRandDE_Dp, bbRandDT_Dp));
        RooProdPdf * bb_missFSP3D_Dp = new RooProdPdf("bb_missFSP3D_Dp", "bb_missFSP3D_Dp", RooArgList(bb_missFSP_MBC_Dp, bb_missFSP_DE_Dp, bb_missFSP_DT_Dp));
        RooProdPdf * scf3D_Dp = new RooProdPdf("scf3D_Dp", "scf3D_Dp", RooArgList(scfMBC_Dp, scfDE_Dp, scfDT_Dp));
        RooProdPdf * signal3D_Dp =  new RooProdPdf("signal3D_Dp","signal3D_Dp", RooArgSet(signalMBC_Dp, signalDE_Dp, signalDT_Dp));
	//Dalitz-
	RooProdPdf * cont3D_Dm = new RooProdPdf("cont3D_Dm", "cont3D_Dm", RooArgList(contMBC_Dm, contDE_Dm, contDT_Dm));
	RooProdPdf * bbRand3D_Dm = new RooProdPdf("bbRand3D_Dm", "bbRand3D_Dm", RooArgList(bbRandMBC_Dm, bbRandDE_Dm, bbRandDT_Dm));
        RooProdPdf * bb_missFSP3D_Dm = new RooProdPdf("bb_missFSP3D_Dm", "bb_missFSP3D_Dm", RooArgList(bb_missFSP_MBC_Dm, bb_missFSP_DE_Dm, bb_missFSP_DT_Dm));
        RooProdPdf * scf3D_Dm = new RooProdPdf("scf3D_Dm", "scf3D_Dm", RooArgList(scfMBC_Dm, scfDE_Dm, scfDT_Dm));
        RooProdPdf * signal3D_Dm =  new RooProdPdf("signal3D_Dm","signal3D_Dm", RooArgSet(signalMBC_Dm, signalDE_Dm, signalDT_Dm));

	//OVERALL NSIG and NBKG floated
	RooRealVar f_scf("f_scf", "scf fraction", 0.1416);//, 0.0, 0.6);//fixed from MC truth info	
	RooRealVar f_bbRand("f_bbRand","fraction of rare bkg", 0.0602);//fixed from MC truth info
        RooRealVar f_bb_missFSP("f_bb_missFSP","fraction of bb bkg", 0.0066);//fixed from MC truth info
	RooRealVar Nsigscf_Ib("Nsigscf_Ib", "Nsigscf_Ib", 840, 100, 1500);
	RooRealVar Nsigscf_I("Nsigscf_I", "Nsigscf_I", 840, 100, 1500);
	RooRealVar Nbkgtot_Ib("Nbkgtot_Ib", "Nbkgtot_Ib", 4650, 800, 6000);
	 RooRealVar Nbkgtot_I("Nbkgtot_I", "Nbkgtot_I", 4650, 800, 6000);


	//Combining P.D.F.'s
	//Dalitz+
	RooAddPdf * fullSig3D_Dp = new RooAddPdf("fullSig3D_Dp","fullSig3D_Dp", RooArgList(*scf3D_Dp, *signal3D_Dp), RooArgList(f_scf));
	RooAddPdf * fullBkg3D_Dp = new RooAddPdf("fullBkg3D_Dp","fullBkg3D_Dp", RooArgList(*bb_missFSP3D_Dp, *bbRand3D_Dp, *cont3D_Dp), RooArgList(f_bb_missFSP, f_bbRand));
	RooAddPdf * combined3D_Dp = new RooAddPdf("combined3D_Dp","combined3D_Dp", RooArgList(*fullSig3D_Dp, *fullBkg3D_Dp), RooArgList(Nsigscf_Ib, Nbkgtot_Ib));
	//Dalitz-
	RooAddPdf * fullSig3D_Dm = new RooAddPdf("fullSig3D_Dm","fullSig3D_Dm", RooArgList(*scf3D_Dm, *signal3D_Dm), RooArgList(f_scf));
	RooAddPdf * fullBkg3D_Dm = new RooAddPdf("fullBkg3D_Dm","fullBkg3D_Dm", RooArgList(*bb_missFSP3D_Dm, *bbRand3D_Dm, *cont3D_Dm), RooArgList(f_bb_missFSP, f_bbRand));
        RooAddPdf * combined3D_Dm = new RooAddPdf("combined3D_Dm","combined3D_Dm", RooArgList(*fullSig3D_Dm, *fullBkg3D_Dm), RooArgList(Nsigscf_I, Nbkgtot_I));
	

	//Combined simultaneous p.d.f.
	RooSimultaneous * combined3D = new RooSimultaneous("combined3D", "combined3D", dalitzCategory);	
	combined3D->addPdf(*combined3D_Dp, "Plus");
	combined3D->addPdf(*combined3D_Dm, "Minus");

//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult * fullFit = combined3D->fitTo(*data, ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;





	
/*
	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40), Title(""));
	data->plotOn(xframe_mbc, Cut("dalitzCategory==dalitzCategory::Plus"));
	combined3D->plotOn(xframe_mbc,Slice(dalitzCategory,"Plus"), Components("signalMBC_Dp"), LineColor(kRed), LineStyle(kDashed), ProjWData(dalitzCategory, *data));
	combined3D->plotOn(xframe_mbc, Slice(dalitzCategory,"Plus"), LineColor(kBlue), ProjWData(dalitzCategory, *data));
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
	data->plotOn(xframe_de, Cut("dalitzCategory==dalitzCategory::Plus"));
	combined3D->plotOn(xframe_de, Slice(dalitzCategory,"Plus"), Components("signalDE_Dp"), LineColor(kRed), LineStyle(kDashed), ProjWData(dalitzCategory, *data));
	combined3D->plotOn(xframe_de, Slice(dalitzCategory,"Plus"), LineColor(kBlue), ProjWData(dalitzCategory, *data));
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
        data->plotOn(xframe_dt, Cut("dalitzCategory==dalitzCategory::Plus"));
	combined3D->plotOn(xframe_dt, Slice(dalitzCategory,"Plus"), Components("signalDT_Dp"), LineColor(kRed), LineStyle(kDashed), ProjWData(dalitzCategory, *data));
        combined3D->plotOn(xframe_dt, Slice(dalitzCategory,"Plus"), LineColor(kBlue), ProjWData(dalitzCategory, *data));
        //combined3D->paramOn(xframe_dt, ProjWData(*data), Layout(0.60));
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

        //S_Minus
        confilename = "/S_minus"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<S_minus.getVal()<<"\t"<<S_minus.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();

	//S_Plus
	confilename = "/S_plus"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<S_plus.getVal()<<"\t"<<S_plus.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();
/*
        //S_I
        confilename = "/S_I"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<S_I.getVal()<<"\t"<<S_I.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();

        //S_Ib
        confilename = "/S_Ib"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<S_Ib.getVal()<<"\t"<<S_Ib.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();
*/
        //Nsigscf_Ib
        confilename = "/Nsigscf_Ib"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsigscf_Ib.getVal()<<"\t"<<Nsigscf_Ib.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();

	//Nsigscf_I
	confilename = "/Nsigscf_I"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsigscf_I.getVal()<<"\t"<<Nsigscf_I.getError()<<"\t"<<fullFit->status()<<endl ;
        fout.close();
	

	

//----------------------------------------------------------------------------------------------//		

        return 0;
}

