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
	

	//const double delta_m = 0.507;
	//const double delta_m = 0.499; //test3
	//const double delta_m = 0.535; //test4
	const double delta_m = 0.835; //test5
	const double tau_b = 1.5344; //test1	
	//const double tau_b = 1.525; //test2

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

	//Corrected mistag rates
	float dw_true[7][2], w_true[7][2];
        w_true[0][0] = 0.484503; w_true[0][1] = 0.475516; dw_true[0][0] = -0.0317004; dw_true[0][1] = 0.0259078;
        w_true[1][0] = 0.418366; w_true[1][1] = 0.414454; dw_true[1][0] = 0.0386381; dw_true[1][1] = 0.00769985;
        w_true[2][0] = 0.302811; w_true[2][1] = 0.3098; dw_true[2][0] = 0.0111002; dw_true[2][1] = 0.0110814;
        w_true[3][0] = 0.214935; w_true[3][1] = 0.212282; dw_true[3][0] = -0.0199043; dw_true[3][1] = -0.00184689;
        w_true[4][0] = 0.153083; w_true[4][1] = 0.151648; dw_true[4][0] = -0.00262208; dw_true[4][1] = -0.00299853;
        w_true[5][0] = 0.0875172; w_true[5][1] = 0.0889206; dw_true[5][0] = 0.00666419; dw_true[5][1] = 0.00236291;
        w_true[6][0] = 0.0210111; w_true[6][1] = 0.0226293; dw_true[6][0] = -0.0030608; dw_true[6][1] = 0.000663429;

	// Read tree in a loop and add to a roodataset
	int true_Nsig = 0, true_Ncont = 0, true_Nrare = 0, true_Nscf = 0, true_Nbb = 0;
	 
	RooDataSet * data = new RooDataSet("data", "data", mainVariables);
	for(int i=0; i<tree->GetEntries(); i++){

	tree->GetEntry(i);


	//True yield values-start
	if(f_issignal == 1){ true_Nsig++; } //signal entry
        if(f_issignal != 1 && f_issgnevt == 1){ true_Nscf++; } //SCF entry
        if(f_mc_type == 20 || f_mc_type == 40){ true_Ncont++; } //continuum entry : charm or uds
        if(f_mc_type == 10 || f_mc_type == 30){ true_Nbb++; } //bb  entry : charged or mixed
        if(f_mc_type == 0){ true_Nrare++; } //rare entry
	//True yield values-end

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
	dw.setVal(-dw_true[i_irbin][exp_bin]);

	ecms.setVal(f_ecms);
	irbin.setIndex(i_irbin);
	issignal.setVal(f_issignal);
	issgnevt.setVal(f_issgnevt);

	if(f_issignal == 1){
	data->add(mainVariables);
	}

	}//tree-reading loop


//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//----------------------------------------------------------------------------------------------//




	//Import dt pdf's
	RooRealVar A("A","A",-2.,2.) ;
	RooRealVar S("S","S",-2.,2.) ;


	RooDtCPSignal * signalDT = new RooDtCPSignal("signalDT", "signalDT", dt, S, A,
                        inexp, CosThetaB, ecms, ntrk_sig, z_err_sig, chisq_tracks_sig,
                        dgf_tracks_sig, ntrk_tag, z_err_tag, chisq_tracks_tag, dgf_tracks_tag,
                        keeptagl_tag, iflavor, w, dw, delta_m, tau_b, true, true, 1.0);	


//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult * fullFit = signalDT->fitTo(*data,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0)) ;
	


	//dt
	RooPlot *xframe_dt = dt.frame(Bins(40));
        data->plotOn(xframe_dt);
	signalDT->plotOn(xframe_dt, LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
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
        pad31->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); xframe_dt->Draw() ;
        pad32->cd(); gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_dt->Draw() ;
	can->SaveAs(canvasName.c_str());
	


	//Save fit results
	
	//To individual fit file
	RooArgSet(fullFit->floatParsFinal()).writeToFile(resultName.c_str());

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


//----------------------------------------------------------------------------------------------//		

        return 0;
}






