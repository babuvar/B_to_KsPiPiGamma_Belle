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

	//Get the input file (one file for now!)
	string filename = opt->GetFilenames().at(0);

	//Get the names of the output plot, result file, consolidated result directory and toy-fit index
	string canvasName = opt->GetPlotname();
	string resultName = opt->GetResultname();
	string ConResDir = opt->GetConsolidatedResultsDir();
	string smc_sample = opt->GetSMCSample();
	int toyindex = opt->GetToyIndex();

	TChain * tree = new TChain("h1"); tree->Add(filename.c_str());
	TChain * tree_KDE = new TChain("h1");
	tree_KDE->Add("/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/rootfiles/for_shape/b_kspipigam_sig_500k_wRMVA_wB2MVA_sel_BCS_branches_Puresignal.root");
		

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
	RooDataSet * data_KDE = new RooDataSet("data_KDE", "data_KDE", mainVariables, Import(*tree_KDE));

//----------------------------------------------------------------------------------------------//
//                                True yield values                                             //
//----------------------------------------------------------------------------------------------//
	
	int true_Nsig = 0, true_Ncont = 0, true_Nrare = 0, true_Nscf = 0, true_Nbb = 0;
	int true_Nsigscf = 0, true_Nbkgtot = 0;	


	//loop over tree
	for(int i=0; i<tree->GetEntries(); i++){

	tree->GetEntry(i);
	/*
	if(f_issignal == 1){ true_Nsig++; } //signal entry
	if(f_issignal != 1 && f_issgnevt == 1){ true_Nscf++; } //SCF entry
	if(f_mc_type == 20 || f_mc_type == 40){ true_Ncont++; } //continuum entry : charm or uds
	if(f_mc_type == 10 || f_mc_type == 30){ true_Nbb++; } //bb  entry : charged or mixed
	if(f_mc_type == 0 && f_issgnevt != 1){ true_Nrare++; } //rare entry
	*/
	if(f_issgnevt== 1){ true_Nsigscf++; }
	else{ true_Nbkgtot++; }

	}


//----------------------------------------------------------------------------------------------//
////                                 Mistag  Rates                                                //
////----------------------------------------------------------------------------------------------//

	//Mistag rates
	float dw_true[7][2], w_true[7][2];

	//MC-Truth mistag rates
	w_true[0][0] = 0.481623; w_true[0][1] = 0.473388; dw_true[0][0] = -0.0317004; dw_true[0][1] = 0.0259078;
        w_true[1][0] = 0.418345; w_true[1][1] = 0.414445; dw_true[1][0] = 0.0386381; dw_true[1][1] = 0.00769985;
        w_true[2][0] = 0.302825; w_true[2][1] = 0.309781; dw_true[2][0] = 0.0111002; dw_true[2][1] = 0.0110814;
        w_true[3][0] = 0.215118; w_true[3][1] = 0.21229; dw_true[3][0] = -0.0199043; dw_true[3][1] = -0.00184689;
        w_true[4][0] = 0.153102; w_true[4][1] = 0.151663; dw_true[4][0] = -0.00262208; dw_true[4][1] = -0.00299853;
        w_true[5][0] = 0.0874591; w_true[5][1] = 0.0889312; dw_true[5][0] = 0.00666419; dw_true[5][1] = 0.00236291;
        w_true[6][0] = 0.0210031; w_true[6][1] = 0.0226292; dw_true[6][0] = -0.0030608; dw_true[6][1] = 0.000663429;

	/*
	//Belle MC mistag calibration
	w_true[0][0] = 0.5; w_true[0][1] = 0.5; dw_true[0][0] = 0.0; dw_true[0][1] = 0.0;
        w_true[1][0] = 0.420827; w_true[1][1] = 0.412222; dw_true[1][0] = 0.0583019; dw_true[1][1] = 0.00408778;
        w_true[2][0] = 0.300296; w_true[2][1] = 0.307838; dw_true[2][0] = 0.00573998; dw_true[2][1] = 0.010326;
        w_true[3][0] = 0.219317; w_true[3][1] = 0.212765; dw_true[3][0] = -0.0392635; dw_true[3][1] = -0.00479522;
        w_true[4][0] = 0.154636; w_true[4][1] = 0.149933; dw_true[4][0] = 0.00474508; dw_true[4][1] = 0.00151989;
        w_true[5][0] = 0.0916131; w_true[5][1] = 0.0913264; dw_true[5][0] = -0.0118737; dw_true[5][1] = 0.0143633;
        w_true[6][0] = 0.0228891; w_true[6][1] = 0.0218754; dw_true[6][0] = -0.00585326; dw_true[6][1] = 0.00189979;
	*/	

//----------------------------------------------------------------------------

	// Import ttree
	RooDataSet * data = new RooDataSet("data", "data", mainVariables);
        for(int i=0; i<tree->GetEntries(); i++){

        tree->GetEntry(i);

        if(i_irbin > 0){//omit 0th r bin
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
        //cout<<"i_irbin = "<<i_irbin<<"\t i_inexp = "<<i_inexp<<"\texp_bin = "<<exp_bin<<endl;
        dw.setVal(dw_true[i_irbin][exp_bin]);
        //w.setVal(f_w);
        //dw.setVal(f_dw);
        ecms.setVal(f_ecms);
        irbin.setIndex(i_irbin);
        issignal.setVal(f_issignal);
        issgnevt.setVal(f_issgnevt);

        data->add(mainVariables);
        }//omit 0th r bin

        }//tree-reading loop


//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//----------------------------------------------------------------------------------------------//

	//Make a workspace to import previously developed pdf's
	RooWorkspace rws("workspace", "workspace");


	//Import Mbc pdf's
	rws.import("pdf_models/fix_shape/cont_shape3D.root:continuum:contMBC"); 
	rws.import("pdf_models/fix_shape/bb_shape3D.root:bb:bbMBC");
        rws.import("pdf_models/fix_shape/rarebkg_shape3D.root:rare:rareMBC");
        rws.import("pdf_models/fix_shape/scf_shape3D.root:scf:scfMBC");

	//Import dE pdf's
	rws.import("pdf_models/fix_shape/cont_shape3D.root:continuum:contDE");
        rws.import("pdf_models/fix_shape/bb_shape3D.root:bb:bbDE");
        rws.import("pdf_models/fix_shape/rarebkg_shape3D.root:rare:rareDE");
        rws.import("pdf_models/fix_shape/scf_shape3D.root:scf:scfDE");

	//Import dt pdf's
        rws.import("pdf_models/fix_shape/cont_shape3D.root:continuum:contDT");
        rws.import("pdf_models/fix_shape/bb_shape3D.root:bb:bbDT");
        rws.import("pdf_models/fix_shape/rarebkg_shape3D.root:rare:rareDT");
        rws.import("pdf_models/fix_shape/scf_shape3D.root:scf:scfDT");
        rws.import("pdf_models/fix_shape/signal_shape3D.root:signal:signalDT");

	
	//Fixing all parameters of imported functions
	TIterator * iter = rws.components().createIterator();
	TObject * var = iter->Next();
	while (var){
	//Check that var-TObject is a parameter, and not a pdf, and if so, fix it's value
	if(strcmp(var->ClassName(),"RooRealVar") == 0 ){
	RooRealVar * par = (RooRealVar *) var;
	par->setConstant(true);
	}
	var = iter->Next();
	}

	//Unfix A and S	
	rws.var("A")->setConstant(false);
	rws.var("S")->setConstant(false);

	//Full model
	//MBC - dE
	RooNDKeysPdf * signalMBCDE = new RooNDKeysPdf("signalMBCDE","signalMBCDE",RooArgSet(mbc,de),*data_KDE,"am") ;

	//MBC
	RooAbsPdf * contMBC = rws.pdf("contMBC");	
	RooAbsPdf * bbMBC = rws.pdf("bbMBC");
	RooAbsPdf * rareMBC = rws.pdf("rareMBC");
	RooAbsPdf * scfMBC = rws.pdf("scfMBC");

	//dE
	RooAbsPdf * contDE = rws.pdf("contDE");
	RooAbsPdf * bbDE = rws.pdf("bbDE");
        RooAbsPdf * rareDE = rws.pdf("rareDE");
        RooAbsPdf * scfDE = rws.pdf("scfDE");

	//dt
        RooAbsPdf * contDT = rws.pdf("contDT");
        RooAbsPdf * bbDT = rws.pdf("bbDT");
        RooAbsPdf * rareDT = rws.pdf("rareDT");
        RooAbsPdf * scfDT = rws.pdf("scfDT");
        RooAbsPdf * signalDT = rws.pdf("signalDT");

	//3D pdf's
	RooProdPdf * cont3D = new RooProdPdf("bkg","bkg",RooArgList(*contMBC, *contDE, *contDT));
        RooProdPdf * bb3D = new RooProdPdf("bbbkg","bbbkg",RooArgList(*bbMBC, *bbDE, *bbDT));
        RooProdPdf * rare3D = new RooProdPdf("cfbkg","cfbkg",RooArgList(*rareMBC, *rareDE, *rareDT));
        RooProdPdf * scf3D = new RooProdPdf("rarebkg","rarebkg",RooArgList(*scfMBC, *scfDE, *scfDT));
        RooProdPdf * signal3D =  new RooProdPdf("signal","signal",RooArgSet(*signalMBCDE, *signalDT));

	/*
 	//ALL YIELDS ARE FLOATING
	//Component yields
	RooRealVar Nsig("N_{sig}","Signal yield", 840, 400, 1200);
        RooRealVar Ncont("N_{cont}","Continuum background yield",4650, 3500, 6000);
        RooRealVar Nrare("N_{rare}","Rare background yield",200, 0, 1000);
        RooRealVar Nscf("N_{scf}","Self cross feed yield",90, 0, 1000);
        RooRealVar Nbb("N_{bb}","B-Bbar background yield",300, 0, 800);

	//Full background model, and combined overall model	
	RooAddPdf * combined3D = new RooAddPdf("combined","combined", RooArgList(*signal3D, *cont3D, *bb3D, *rare3D, *scf3D), RooArgList(Nsig, Ncont, Nbb, Nrare, Nscf));
	*/

	//OVERALL NSIG and NBKG floated
	RooRealVar f_scf("f_scf","scf fraction", 0.3174, 0.0, 0.6);//fixed from MC truth info	
	RooRealVar f_rare("f_rare","fraction of rare bkg", 0.0363);//fixed from MC truth info
	RooRealVar f_bb("f_bb","fraction of bb bkg", 0.0307);//fixed from MC truth info
	RooRealVar Nsigscf("N_{sig/scf}","Signal + SCF yield", 840, 300, 1500);
	RooRealVar Nbkgtot("N_{allbkg}","Combined background yield",4650, 3500, 6000);
	//RooRealVar Nsigscf("N_{sig/scf}","Signal + SCF yield", true_Nsigscf);
	//RooRealVar Nbkgtot("N_{allbkg}","Combined background yield", true_Nbkgtot);

	//Combining P.D.F.'s
	RooAddPdf * fullSig3D = new RooAddPdf("fullSig3D","fullSig3D", RooArgList(*scf3D, *signal3D), RooArgList(f_scf));
	RooAddPdf * fullBkg3D = new RooAddPdf("fullBkg3D","fullBkg3D", RooArgList(*rare3D, *bb3D, *cont3D), RooArgList(f_rare, f_bb));
	RooAddPdf * combined3D = new RooAddPdf("combined","combined", RooArgList(*fullSig3D, *fullBkg3D), RooArgList(Nsigscf, Nbkgtot));

//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult * fullFit = combined3D->fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	
	/*
	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	data->plotOn(xframe_mbc);
	combined3D->plotOn(xframe_mbc,Components("signalMBCDE"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
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
	RooPlot *xframe_de = de.frame(Bins(40));
	data->plotOn(xframe_de);
	combined3D->plotOn(xframe_de, Components("signalMBCDE"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
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
	RooPlot *xframe_dt = dt.frame(Bins(40));
        data->plotOn(xframe_dt);
	combined3D->plotOn(xframe_dt,Components("signalDT"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
        combined3D->plotOn(xframe_dt, LineColor(kBlue), ProjWData(*data));
        //combined3D->paramOn(xframe_dt, ProjWData(*data));
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
	can->SaveAs(canvasName.c_str());
	*/

	//Save fit results
	
	//To individual fit file
	RooArgSet(fullFit->floatParsFinal()).writeToFile(resultName.c_str());

	//To consolidated parameter result files
	ofstream fout; string confilename, fileSuffix;
	fileSuffix = "_" + smc_sample + ".txt";
	//A
	confilename = "/A"; confilename = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
	fout<<toyindex<<"\t"<<rws.var("A")->getVal()<<"\t"<<rws.var("A")->getError()<<"\t"<<0.0<<endl ;
	fout.close();
	//S
	confilename = "/S"; confilename = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
	fout<<toyindex<<"\t"<<rws.var("S")->getVal()<<"\t"<<rws.var("S")->getError()<<"\t"<<0.0<<endl ;
	fout.close();
	/*
	//Nsig
	confilename = "/Nsig"; confilename  = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsig.getVal()<<"\t"<<Nsig.getError()<<"\t"<<true_Nsig<<endl ;
        fout.close();
	//Ncont
	confilename = "/Ncont"; confilename  = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Ncont.getVal()<<"\t"<<Ncont.getError()<<"\t"<<true_Ncont<<endl ;
        fout.close();
	//Nscf
	confilename = "/Nrare"; confilename  = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nrare.getVal()<<"\t"<<Nrare.getError()<<"\t"<<true_Nrare<<endl ;
        fout.close();
	//Nscf
	confilename = "/Nscf"; confilename  = ConResDir + confilename + fileSuffix;
	fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nscf.getVal()<<"\t"<<Nscf.getError()<<"\t"<<true_Nscf<<endl ;
        fout.close();
	//Nbb
	confilename = "/Nbb"; confilename  = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nbb.getVal()<<"\t"<<Nbb.getError()<<"\t"<<true_Nbb<<endl ;
        fout.close();
	

	//Nsigscf
	confilename = "/Nsigscf"; confilename  = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsigscf.getVal()<<"\t"<<Nsigscf.getError()<<"\t"<<true_Nsigscf<<endl ;
        fout.close();

	//Nbkgtot
	confilename = "/Nbkgtot"; confilename  = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nbkgtot.getVal()<<"\t"<<Nbkgtot.getError()<<"\t"<<true_Nbkgtot<<endl ;
        fout.close();
	*/
	

//----------------------------------------------------------------------------------------------//		

        return 0;
}

