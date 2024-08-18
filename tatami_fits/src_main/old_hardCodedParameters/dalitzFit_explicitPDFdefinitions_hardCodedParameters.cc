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

//#include "RooCFunction3Binding.h"


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

	//Get the input file (one file for now!)
	string filename = opt->GetFilenames().at(0);

	//Get the names of the output plot, result file, consolidated result directory and toy-fit index
	string canvasName = opt->GetPlotname();
	string resultName = opt->GetResultname();
	string ConResDir = opt->GetConsolidatedResultsDir();
	string smc_sample = opt->GetSMCSample();
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
	//const double tau_b = 1.5344;	
	const double tau_b = 1.525;

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

	//Common Argus-threshold
	RooRealVar commonArgus_m0("commonArgus_m0", "commonArgus_m0", 5.29, 5.28, 5.30) ;

	//BB-bkg
	//MBC
	//Argus
	RooRealVar bbMBC_c("bbMBC_c","bbMBC_c", -8.11070); //-5,-140,2) ;
        RooRealVar bbMBC_m0("bbMBC_m0","bbMBC_m0", 5.2885); //5.29, 5.28, 5.30) ;
        //RooArgusBG bbMBC_argus("bbMBC_argus","bbMBC_argus", mbc, bbMBC_m0, bbMBC_c);
        RooArgusBG bbMBC_argus_Dp("bbMBC_argus_Dp","bbMBC_argus_Dp", mbc, commonArgus_m0, bbMBC_c);
	RooArgusBG bbMBC_argus_Dm("bbMBC_argus_Dm","bbMBC_argus_Dm", mbc, commonArgus_m0, bbMBC_c);
	//Gaussian
	RooRealVar bbMBC_mu("bbMBC_mu", "bbMBC_mu", 5.2749); //5.2797, 5.270,5.290) ;
        RooRealVar bbMBC_sig("bbMBC_sig", "bbMBC_sig", 0.0020); //0.0020, 0.001, 0.05) ;
        RooGaussian bbMBC_gauss_Dp("bbMBC_gauss_Dp","bbMBC_gauss_Dp", mbc, bbMBC_mu, bbMBC_sig) ;
	RooGaussian bbMBC_gauss_Dm("bbMBC_gauss_Dm","bbMBC_gauss_Dm", mbc, bbMBC_mu, bbMBC_sig) ;
        //combined
        RooRealVar bbMBC_f("bbMBC_f","bbMBC_f", 0.0); //, 0.0 ,0.9) ;
        RooAddPdf bbMBC_Dp("bbMBC_Dp","bbMBC_Dp", RooArgList(bbMBC_gauss_Dp, bbMBC_argus_Dp), RooArgList(bbMBC_f)) ;	
	RooAddPdf bbMBC_Dm("bbMBC_Dm","bbMBC_Dm", RooArgList(bbMBC_gauss_Dm, bbMBC_argus_Dm), RooArgList(bbMBC_f)) ;

	//dE
	RooRealVar bbDE_al("bbDE_al", "bbDE_al", -6.01751); //-6.04, -10., 0.) ;
        RooExponential bbDE_Dp("bbDE_Dp", "bbDE_Dp", de, bbDE_al) ;
	RooExponential bbDE_Dm("bbDE_Dm", "bbDE_Dm", de, bbDE_al) ;

	//dt
	RooRealVar bbDT_tau("bbDT_tau", "bbDT_tau", 1.40); //1.4, 0.5, 1.8);
        RooRealVar bbDT_mu_l("bbDT_mu_l", "bbDT_mu_l", -0.698928); //-0.62, -0.7, 0.5);
        RooRealVar bbDT_mu_d("bbDT_mu_d", "bbDT_mu_d", 0.0028751); //-0.04, -0.1, 0.1);
        RooRealVar bbDT_f_delt("bbDT_f_delt", "bbDT_f_delt", 0.150); //0.15, 0.1,1);
        RooRealVar bbDT_f_tail("bbDT_f_tail", "bbDT_f_tail", 0.19111); //0.218, 0, 1);
        RooRealVar bbDT_S_main("bbDT_S_main", "bbDT_S_main", 2.8550); //2.79, 0.4, 10);
        RooRealVar bbDT_S_tail("bbDT_S_tail", "bbDT_S_tail", 3.4377); //3.43, 1, 40);

        RooDtBkg bbDT_Dp("bbDT_Dp", "bbDT_Dp", dt,
                bbDT_tau, bbDT_mu_l, bbDT_mu_d, bbDT_f_delt, bbDT_f_tail, bbDT_S_main, bbDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	RooDtBkg bbDT_Dm("bbDT_Dm", "bbDT_Dm", dt,
                bbDT_tau, bbDT_mu_l, bbDT_mu_d, bbDT_f_delt, bbDT_f_tail, bbDT_S_main, bbDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);


//------------------------------------------------------------------------------------------------------
	//Continuum
	//MBC
	RooRealVar  contMBC_c("contMBC_c", "contMBC_c", -7.39000); //-9.8, -200, -1);
        RooRealVar  contMBC_m0("contMBC_m0","contMBC_m0", 5.2890);// 5.29, 5.28, 5.30);
        //RooArgusBG  contMBC("contMBC","contMBC", mbc, contMBC_m0, contMBC_c);
	RooArgusBG  contMBC_Dp("contMBC_Dp","contMBC_Dp", mbc, commonArgus_m0, contMBC_c);
	RooArgusBG  contMBC_Dm("contMBC_Dm","contMBC_Dm", mbc, commonArgus_m0, contMBC_c);

	//dE
	RooRealVar  contDE_al("contDE_al", "contDE_al", -0.915794); //-1.07,-10.,0.1) ;
        RooExponential contDE_Dp("contDE_Dp","contDE_Dp", de, contDE_al);
	RooExponential contDE_Dm("contDE_Dm","contDE_Dm", de, contDE_al);

	//dt
	RooRealVar contDT_tau("contDT_tau", "contDT_tau", 1.40); //1.4, 0.5, 1.8);
        RooRealVar contDT_mu_l("contDT_mu_l", "contDT_mu_l", -0.699099); //-0.61, -0.7, 0.5);
        RooRealVar contDT_mu_d("contDT_mu_d", "contDT_mu_d", -0.0126835); //-0.011, -0.1, 0.1);
        RooRealVar contDT_f_delt("contDT_f_delt", "contDT_f_delt", 0.110); //0.11, 0.1, 1);
        RooRealVar contDT_f_tail("contDT_f_tail", "contDT_f_tail", 0.16327); //0.1778, 0, 1);
        RooRealVar contDT_S_main("contDT_S_main", "contDT_S_main", 1.7924); //1.825, 0.4, 10);
        RooRealVar contDT_S_tail("contDT_S_tail" , "contDT_S_tail", 3.9538); //3.77, 1, 40);

        RooDtBkg contDT_Dp("contDT_Dp", "contDT_Dp", dt,
                contDT_tau, contDT_mu_l, contDT_mu_d, contDT_f_delt, contDT_f_tail, contDT_S_main, contDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	RooDtBkg contDT_Dm("contDT_Dm", "contDT_Dm", dt,
                contDT_tau, contDT_mu_l, contDT_mu_d, contDT_f_delt, contDT_f_tail, contDT_S_main, contDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);


//------------------------------------------------------------------------------------------------------
	//Rare-bkg
	//MBC
	//CrystalBall
	RooRealVar rareMBC_mean("rareMBC_mean", "rareMBC_mean", 5.2803); //5.2812, 5.275, 5.285) ;
        RooRealVar rareMBC_sigma("rareMBC_sigma", "rareMBC_sigma", 0.0033531); //0.00300, 5e-4, 5e-3);
        RooRealVar rareMBC_alpha("rareMBC_alpha", "rareMBC_alpha", 0.51645); //0.25, 0.,5.);
        RooRealVar rareMBC_n("rareMBC_n", "rareMBC_n", 20.000); //5, 0, 20);
        RooCBShape rareMBC_CB_Dp("rareMBC_CB_Dp", "rareMBC_CB_Dp", mbc, rareMBC_mean, rareMBC_sigma, rareMBC_alpha, rareMBC_n);
	RooCBShape rareMBC_CB_Dm("rareMBC_CB_Dm", "rareMBC_CB_Dm", mbc, rareMBC_mean, rareMBC_sigma, rareMBC_alpha, rareMBC_n);
        //Argus
        RooRealVar rareMBC_c("rareMBC_c", "rareMBC_c", -58.2059); //-44, -140, -1) ;
        RooRealVar rareMBC_m0("rareMBC_m0}", "rareMBC_m0", 5.2889); //5.29, 5.28, 5.30);
        //RooArgusBG rareMBC_Argus("rareMBC_Argus", "rareMBC_Argus", mbc, rareMBC_m0, rareMBC_c) ;
	RooArgusBG rareMBC_Argus_Dp("rareMBC_Argus_Dp", "rareMBC_Argus_Dp", mbc, commonArgus_m0, rareMBC_c) ;
	RooArgusBG rareMBC_Argus_Dm("rareMBC_Argus_Dm", "rareMBC_Argus_Dm", mbc, commonArgus_m0, rareMBC_c) ;
	//Combined
	RooRealVar rareMBC_f("rareMBC_f","rareMBC_f", 0.19513); //0.29, 0., 1.0) ;
        RooAddPdf rareMBC_Dp("rareMBC_Dp","rareMBC_Dp",RooArgList(rareMBC_CB_Dp, rareMBC_Argus_Dp), RooArgList(rareMBC_f)) ;	
	RooAddPdf rareMBC_Dm("rareMBC_Dm","rareMBC_Dm",RooArgList(rareMBC_CB_Dm, rareMBC_Argus_Dm), RooArgList(rareMBC_f)) ;

	//dE
	RooRealVar rareDE_alpha("rareDE_alpha", "rareDE_alpha", -3.95308); //-3.62, -10., 0.1);
        RooExponential rareDE_Dp("rareDE_Dp", "rareDE_Dp", de, rareDE_alpha);
	RooExponential rareDE_Dm("rareDE_Dm", "rareDE_Dm", de, rareDE_alpha);
	
	//dt
	RooRealVar rareDT_tau("rareDT_tau", "rareDT_tau", 1.40); //1.4, 0.5,1.8);
        RooRealVar rareDT_mu_l("rareDT_mu_l", "rareDT_mu_l", -0.798748); //-0.7,-0.8,0.5);
        RooRealVar rareDT_mu_d("rareDT_mu_d", "rareDT_mu_d", -0.0776618); //-0.079, -0.1,0.1);
        RooRealVar rareDT_f_delt("rareDT_f_delt", "rareDT_f_delt", 0.150); //0.15, 0.1,1);
        RooRealVar rareDT_f_tail("rareDT_f_tail", "rareDT_f_tail", 0.32953); //0.322, 0, 1.);
        RooRealVar rareDT_S_main("rareDT_S_main", "rareDT_S_main", 2.4685); //2.487, 0.4,10);
        RooRealVar rareDT_S_tail("rareDT_S_tail", "rareDT_S_tail", 3.3250); //3.31, 1, 40);


        RooDtBkg rareDT_Dp("rareDT_Dp", "rareDT_Dp", dt,
                rareDT_tau, rareDT_mu_l, rareDT_mu_d, rareDT_f_delt, rareDT_f_tail, rareDT_S_main, rareDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	RooDtBkg rareDT_Dm("rareDT_Dm", "rareDT_Dm", dt,
                rareDT_tau, rareDT_mu_l, rareDT_mu_d, rareDT_f_delt, rareDT_f_tail, rareDT_S_main, rareDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	

//------------------------------------------------------------------------------------------------------
	//Self cross feed
        //MBC
	//Gaussian
	RooRealVar scfMBC_mean("scfMBC_mean", "scfMBC_mean", 5.2793); //5.27910, 5.270, 5.290) ;
        RooRealVar scfMBC_sigma("scfMBC_sigma", "scfMBC_sigma", 0.0034532); //0.00361, 0.001, 0.05) ;
        RooGaussian scfMBC_Gauss_Dp("scfMBC_Gauss_Dp", "scfMBC_Gauss_Dp", mbc, scfMBC_mean, scfMBC_sigma) ;
	RooGaussian scfMBC_Gauss_Dm("scfMBC_Gauss_Dm", "scfMBC_Gauss_Dm", mbc, scfMBC_mean, scfMBC_sigma) ;
        //Argus
        RooRealVar scfMBC_c("scfMBC_c", "scfMBC_c", -110.955); //-96.9, -140, -1) ;
        RooRealVar scfMBC_m0("scfMBC_m0", "scfMBC_m0", 5.2887); //5.29, 5.28, 5.30);
        //RooArgusBG scfMBC_Argus("scfMBC_Argus", "scfMBC_Argus", mbc, scfMBC_m0, scfMBC_c) ;
	RooArgusBG scfMBC_Argus_Dp("scfMBC_Argus_Dp", "scfMBC_Argus_Dp", mbc, commonArgus_m0, scfMBC_c) ;
	RooArgusBG scfMBC_Argus_Dm("scfMBC_Argus_Dm", "scfMBC_Argus_Dm", mbc, commonArgus_m0, scfMBC_c) ;

        //combined
        RooRealVar scfMBC_f("scfMBC_f", "scfMBC_f", 0.20454); //0.217, 0., 0.9) ;
        RooAddPdf scfMBC_Dp("scfMBC_Dp", "scfMBC_Dp", RooArgList(scfMBC_Gauss_Dp, scfMBC_Argus_Dp), RooArgList(scfMBC_f)) ;
	RooAddPdf scfMBC_Dm("scfMBC_Dm", "scfMBC_Dm", RooArgList(scfMBC_Gauss_Dm, scfMBC_Argus_Dm), RooArgList(scfMBC_f)) ;


        //dE
	//Exponential
	RooRealVar scfDE_alpha("scfDE_alpha" , "scfDE_alpha", -3.32831); //-3.72, -10., 0.) ;
        RooExponential scfDE_Expo_Dp("scfDE_Expo_Dp", "scfDE_Expo_Dp", de, scfDE_alpha) ;
	RooExponential scfDE_Expo_Dm("scfDE_Expo_Dm", "scfDE_Expo_Dm", de, scfDE_alpha) ;
        //Gaussian
        RooRealVar scfDE_mean("scfDE_mean", "scfDE_mean", -0.0472227); //-0.038, -0.15, 0.05) ;
        RooRealVar scfDE_sigma("scfDE_sigma", "scfDE_sigma", 0.058366); //0.058, 0.001, 0.1) ;
        RooGaussian scfDE_Gauss_Dp("scfDE_Gauss_Dp","scfDE_Gauss_Dp", de, scfDE_mean, scfDE_sigma) ;
	RooGaussian scfDE_Gauss_Dm("scfDE_Gauss_Dm","scfDE_Gauss_Dm", de, scfDE_mean, scfDE_sigma) ;
        //combined
        RooRealVar scfDE_f("scfDE_f","scfDE_f", 0.78138); //0.793, 0., 0.99) ;
        RooAddPdf scfDE_Dp("scfDE_Dp","scfDE_Dp",RooArgList(scfDE_Expo_Dp, scfDE_Gauss_Dp), RooArgList(scfDE_f)) ;
	RooAddPdf scfDE_Dm("scfDE_Dm","scfDE_Dm",RooArgList(scfDE_Expo_Dm, scfDE_Gauss_Dm), RooArgList(scfDE_f)) ;

	//dt
	RooRealVar scfDT_tau("scfDT_tau", "scfDT_tau", 1.40); //1.4, 0.5, 1.8);
        RooRealVar scfDT_mu_l("scfDT_mu_l", "scfDT_mu_l", -0.690042); //-0.67, -0.7, 0.5);
        RooRealVar scfDT_mu_d("scfDT_mu_d", "scfDT_mu_d", -0.0838778); //-0.040, -0.1, 0.1);
        RooRealVar scfDT_f_delt("scfDT_f_delt", "scfDT_f_delt", 0.150); //0.15, 0.1, 1);
        RooRealVar scfDT_f_tail("scfDT_f_tail", "scfDT_f_tail", 0.34548); //0.319,0, 1);
        RooRealVar scfDT_S_main("scfDT_S_main", "scfDT_S_main", 2.8412); //2.84, 0.4, 10);
        RooRealVar scfDT_S_tail("scfDT_S_tail", "scfDT_S_tail", 3.0532); //3.11, 1, 40);

        RooDtBkg scfDT_Dp("scfDT_Dp", "scfDT_Dp", dt,
                scfDT_tau, scfDT_mu_l, scfDT_mu_d, scfDT_f_delt, scfDT_f_tail, scfDT_S_main, scfDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);
	RooDtBkg scfDT_Dm("scfDT_Dm", "scfDT_Dm", dt,
                scfDT_tau, scfDT_mu_l, scfDT_mu_d, scfDT_f_delt, scfDT_f_tail, scfDT_S_main, scfDT_S_tail,
                inexp, ntrk_sig, z_err_sig, ntrk_tag, z_err_tag, true);

//------------------------------------------------------------------------------------------------------
	
        //Signal
        //MBC
        RooRealVar meanMBC_fudge("meanMBC_fudge", "meanMBC_fudge", 0.0005, -0.01, 0.01);
        RooRealVar sigmaMBC_fudge("sigmaMBC_fudge", "sigmaMBC_fudge", 1.0005, 0.5, 2.0);
        RooFormulaVar signalMBC_mean("signalMBC_mean", "5.2796+meanMBC_fudge", RooArgList(meanMBC_fudge));
        RooFormulaVar signalMBC_sigma("signalMBC_sigma", "0.0029357*sigmaMBC_fudge", RooArgList(sigmaMBC_fudge));
        RooRealVar signalMBC_alpha("signalMBC_alpha", "signalMBC_alpha", 1.2037); //1.16, 0., 5.);
        RooRealVar signalMBC_n("signalMBC_n", "signalMBC_n", 20.0); //19, 0, 20);
        RooCBShape signalMBC_Dp("signalMBC_Dp", "signalMBC_Dp", mbc, signalMBC_mean, signalMBC_sigma, signalMBC_alpha, signalMBC_n);
	RooCBShape signalMBC_Dm("signalMBC_Dm", "signalMBC_Dm", mbc, signalMBC_mean, signalMBC_sigma, signalMBC_alpha, signalMBC_n);

	//dE
        RooRealVar meanDE_fudge("meanDE_fudge", "meanDE_fudge", 0.0005, -0.01, 0.01);
        RooRealVar sigmaDE_fudge("sigmaDE_fudge", "sigmaDE_fudge", 1.0005,0.5,2.0);
        RooFormulaVar signalDE_mean("signalDE_mean", "-0.00689755+meanDE_fudge", RooArgList(meanDE_fudge));
        RooFormulaVar signalDE_sigma("signalDE_sigma", "0.037056*sigmaDE_fudge", RooArgList(sigmaDE_fudge));	
        RooRealVar signalDE_alpha("signalDE_alpha", "signalDE_alpha", 0.59009); //0.60, 0, 5);
        RooRealVar signalDE_n("signalDE_n", "signalDE_n", 7.6855); //6.2,1, 50);
        RooCBShape signalDE_Dp("signalDE_Dp", "signalDE_Dp", de, signalDE_mean, signalDE_sigma, signalDE_alpha, signalDE_n);
	RooCBShape signalDE_Dm("signalDE_Dm", "signalDE_Dm", de, signalDE_mean, signalDE_sigma, signalDE_alpha, signalDE_n);

        //dt
	RooRealVar A("A", "A", -2., 2.) ;
	RooRealVar S_plus("S_plus", "S_plus",  -2., 2.);
	RooRealVar S_minus("S_minus", "S_minus",  -2., 2.);
        //RooRealVar S_Ib("S_Ib", "S_Ib", -2., 2.) ; //Dalitz - I-bar
	//RooRealVar S_I("S_I", "S_I", -2., 2.) ; //Dalitz - I
	RooFormulaVar S_I("S_I", "(S_plus+S_minus)/2", RooArgList(S_plus, S_minus));
	RooFormulaVar S_Ib("S_Ib", "(S_plus-S_minus)/2", RooArgList(S_plus, S_minus));

        RooDtCPSignal signalDT_Dp("signalDT_Dp", "signalDT_Dp", dt, S_Ib, A,
                        inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig,
                        dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag,
                        keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);
	RooDtCPSignal signalDT_Dm("signalDT_Dm", "signalDT_Dm", dt, S_I, A,
                        inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig,
                        dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag,
                        keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);

//------------------------------------------------------------------------------------------------------	


	//3D pdf's
	//Dalitz+
	RooProdPdf * cont3D_Dp = new RooProdPdf("cont3D_Dp", "cont3D_Dp", RooArgList(contMBC_Dp, contDE_Dp, contDT_Dp));
        RooProdPdf * bb3D_Dp = new RooProdPdf("bb3D_Dp", "bb3D_Dp", RooArgList(bbMBC_Dp, bbDE_Dp, bbDT_Dp));
        RooProdPdf * rare3D_Dp = new RooProdPdf("rare3D_Dp", "rare3D_Dp", RooArgList(rareMBC_Dp, rareDE_Dp, rareDT_Dp));
        RooProdPdf * scf3D_Dp = new RooProdPdf("scf3D_Dp", "scf3D_Dp", RooArgList(scfMBC_Dp, scfDE_Dp, scfDT_Dp));
        RooProdPdf * signal3D_Dp =  new RooProdPdf("signal3D_Dp","signal3D_Dp", RooArgSet(signalMBC_Dp, signalDE_Dp, signalDT_Dp));
	//Dalitz-
	RooProdPdf * cont3D_Dm = new RooProdPdf("cont3D_Dm", "cont3D_Dm", RooArgList(contMBC_Dm, contDE_Dm, contDT_Dm));
        RooProdPdf * bb3D_Dm = new RooProdPdf("bb3D_Dm", "bb3D_Dm", RooArgList(bbMBC_Dm, bbDE_Dm, bbDT_Dm));
        RooProdPdf * rare3D_Dm = new RooProdPdf("rare3D_Dm", "rare3D_Dm", RooArgList(rareMBC_Dm, rareDE_Dm, rareDT_Dm));
        RooProdPdf * scf3D_Dm = new RooProdPdf("scf3D_Dm", "scf3D_Dm", RooArgList(scfMBC_Dm, scfDE_Dm, scfDT_Dm));
        RooProdPdf * signal3D_Dm =  new RooProdPdf("signal3D_Dm","signal3D_Dm", RooArgSet(signalMBC_Dm, signalDE_Dm, signalDT_Dm));

	//OVERALL NSIG and NBKG floated
	RooRealVar f_scf("f_scf", "scf fraction", 0.3097);//, 0.0, 0.6);//fixed from MC truth info	
	RooRealVar f_rare("f_rare", "fraction of rare bkg", 0.0412);//fixed from MC truth info
	RooRealVar f_bb("f_bb", "fraction of bb bkg", 0.0356);//fixed from MC truth info
	RooRealVar Nsigscf_Ib("Nsigscf_Ib", "Nsigscf_Ib", 840, 300, 1500);
	RooRealVar Nsigscf_I("Nsigscf_I", "Nsigscf_I", 840, 300, 1500);
	RooRealVar Nbkgtot_Ib("Nbkgtot_Ib", "Nbkgtot_Ib", 4650, 1500, 6000);
	 RooRealVar Nbkgtot_I("Nbkgtot_I", "Nbkgtot_I", 4650, 1500, 6000);


	//Combining P.D.F.'s
	//Dalitz+
	RooAddPdf * fullSig3D_Dp = new RooAddPdf("fullSig3D_Dp","fullSig3D_Dp", RooArgList(*scf3D_Dp, *signal3D_Dp), RooArgList(f_scf));
	RooAddPdf * fullBkg3D_Dp = new RooAddPdf("fullBkg3D_Dp","fullBkg3D_Dp", RooArgList(*rare3D_Dp, *bb3D_Dp, *cont3D_Dp), RooArgList(f_rare, f_bb));
	RooAddPdf * combined3D_Dp = new RooAddPdf("combined3D_Dp","combined3D_Dp", RooArgList(*fullSig3D_Dp, *fullBkg3D_Dp), RooArgList(Nsigscf_Ib, Nbkgtot_Ib));
	//Dalitz-
	RooAddPdf * fullSig3D_Dm = new RooAddPdf("fullSig3D_Dm","fullSig3D_Dm", RooArgList(*scf3D_Dm, *signal3D_Dm), RooArgList(f_scf));
        RooAddPdf * fullBkg3D_Dm = new RooAddPdf("fullBkg3D_Dm","fullBkg3D_Dm", RooArgList(*rare3D_Dm, *bb3D_Dm, *cont3D_Dm), RooArgList(f_rare, f_bb));
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
/*
	//To consolidated parameter result files
	ofstream fout; string confilename, fileSuffix;
	fileSuffix = "_" + smc_sample + ".txt";
	//A
	confilename = "/A"; confilename = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
	fout<<toyindex<<"\t"<<A.getVal()<<"\t"<<A.getError()<<"\t"<<0.0<<endl ;
	fout.close();
	//S
	confilename = "/S"; confilename = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
	fout<<toyindex<<"\t"<<S.getVal()<<"\t"<<S.getError()<<"\t"<<0.0<<endl ;
	fout.close();

	//Nsigscf
	confilename = "/Nsigscf"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsigscf.getVal()<<"\t"<<Nsigscf.getError()<<"\t"<<0.0<<endl ;
        fout.close();
	
	//Nbkgtot
	confilename = "/Nbkgtot"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nbkgtot.getVal()<<"\t"<<Nbkgtot.getError()<<"\t"<<0.0<<endl ;
        fout.close();
*/
	

//----------------------------------------------------------------------------------------------//		

        return 0;
}

