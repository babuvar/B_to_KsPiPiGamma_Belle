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
	tree->SetBranchAddress("irbin",&i_irbin);
	tree->SetBranchAddress("issgnevt",&f_issgnevt);

	// Variable definitions
        RooRealVar mbc("mbc","m_{bc} [GeV]",5.20,5.30) ;
        RooRealVar de("de","#DeltaE [GeV]",-0.2,0.2) ;

	RooArgSet mainVariables(mbc, de);	

	// Import ttree
	int true_nSigSCF = 0, true_nBkgTot = 0;
RooDataSet * data = new RooDataSet("data", "data", mainVariables);        for(int i=0; i<tree->GetEntries(); i++){

        tree->GetEntry(i);


	if(i_irbin > 0){//omit 0th r bin

	if(f_issgnevt == 1){true_nSigSCF++;}
	else{true_nBkgTot++;}

        mbc.setVal(f_mbc);
        de.setVal(f_de);

        data->add(mainVariables);
        }//omit 0th r bin

        }//tree-reading loop


//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//----------------------------------------------------------------------------------------------//

	//Make a workspace to import previously developed pdf's
	RooWorkspace rws("workspace", "workspace");

	//Import Mbc pdf's
	rws.import("pdf_models/fix_shape/cont_shape.root:continuum:contMBC"); 
	rws.import("pdf_models/fix_shape/bb_shape.root:bb:bbMBC");
        rws.import("pdf_models/fix_shape/rarebkg_shape.root:rare:rareMBC");
        rws.import("pdf_models/fix_shape/scf_shape.root:scf:scfMBC");
        rws.import("pdf_models/fix_shape/signal_shape.root:signal:signalMBC");

	//Import dE pdf's
	rws.import("pdf_models/fix_shape/cont_shape.root:continuum:contDE");
        rws.import("pdf_models/fix_shape/bb_shape.root:bb:bbDE");
        rws.import("pdf_models/fix_shape/rarebkg_shape.root:rare:rareDE");
        rws.import("pdf_models/fix_shape/scf_shape.root:scf:scfDE");
        rws.import("pdf_models/fix_shape/signal_shape.root:signal:signalDE");

	
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


	//Full model
	//MBC
	RooAbsPdf * contMBC = rws.pdf("contMBC");	
	RooAbsPdf * bbMBC = rws.pdf("bbMBC");
	RooAbsPdf * rareMBC = rws.pdf("rareMBC");
	RooAbsPdf * scfMBC = rws.pdf("scfMBC");
	RooAbsPdf * signalMBC = rws.pdf("signalMBC");

	//dE
	RooAbsPdf * contDE = rws.pdf("contDE");
	RooAbsPdf * bbDE = rws.pdf("bbDE");
        RooAbsPdf * rareDE = rws.pdf("rareDE");
        RooAbsPdf * scfDE = rws.pdf("scfDE");
        RooAbsPdf * signalDE = rws.pdf("signalDE");


	//3D pdf's
	RooProdPdf * cont2D = new RooProdPdf("bkg","bkg",RooArgList(*contMBC, *contDE));
        RooProdPdf * bb2D = new RooProdPdf("bbbkg","bbbkg",RooArgList(*bbMBC, *bbDE));
        RooProdPdf * rare2D = new RooProdPdf("cfbkg","cfbkg",RooArgList(*rareMBC, *rareDE));
        RooProdPdf * scf2D = new RooProdPdf("rarebkg","rarebkg",RooArgList(*scfMBC, *scfDE));
        RooProdPdf * signal2D =  new RooProdPdf("signal","signal",RooArgSet(*signalMBC, *signalDE));

/*	
 	//ALL YIELDS ARE FLOATING
	//Component yields
	//RooRealVar Nsig("N_{sig}","Signal yield", 840, 400, 1200);
        //RooRealVar Ncont("N_{cont}","Continuum background yield",4650, 3500, 6000);
        //RooRealVar Nrare("N_{rare}","Rare background yield",200, 0, 1000);
        //RooRealVar Nscf("N_{scf}","Self cross feed yield",90, 0, 1000);
        //RooRealVar Nbb("N_{bb}","B-Bbar background yield",300, 0, 800);
        RooRealVar Nsig("N_{sig}","Signal yield", true_Nsig);
        RooRealVar Ncont("N_{cont}","Continuum background yield", true_Ncont);
        RooRealVar Nrare("N_{rare}","Rare background yield", true_Nrare);
        RooRealVar Nscf("N_{scf}","Self cross feed yield", true_Nscf);
        RooRealVar Nbb("N_{bb}","B-Bbar background yield", true_Nbb);
	//Full background model, and combined overall model	
	RooAddPdf * combined3D = new RooAddPdf("combined","combined", RooArgList(*signal3D, *cont3D, *bb3D, *rare3D, *scf3D), RooArgList(Nsig, Ncont, Nbb, Nrare, Nscf));
*/	

	//OVERALL NSIG and NBKG floated
	RooRealVar f_scf("f_scf","scf fraction", 0.3100);//, 0.0, 0.6);//fixed from MC truth info	
	RooRealVar f_rare("f_rare","fraction of rare bkg", 0.0414);//fixed from MC truth info
	RooRealVar f_bb("f_bb","fraction of bb bkg", 0.0359);//fixed from MC truth info
	RooRealVar Nsigscf("N_{sig/scf}","Signal + SCF yield", 840, 300, 1500);
	RooRealVar Nbkgtot("N_{allbkg}","Combined background yield",4650, 2500, 5000);
	//RooRealVar Nsigscf("N_{sig/scf}","Signal + SCF yield", true_Nsigscf);
	//RooRealVar Nbkgtot("N_{allbkg}","Combined background yield", true_Nbkgtot);

	//Combining P.D.F.'s
	RooAddPdf * fullSig2D = new RooAddPdf("fullSig2D","fullSig2D", RooArgList(*scf2D, *signal2D), RooArgList(f_scf));
	RooAddPdf * fullBkg2D = new RooAddPdf("fullBkg2D","fullBkg2D", RooArgList(*rare2D, *bb2D, *cont2D), RooArgList(f_rare, f_bb));
	RooAddPdf * combined2D = new RooAddPdf("combined","combined", RooArgList(*fullSig2D, *fullBkg2D), RooArgList(Nsigscf, Nbkgtot));

//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult * fullFit = combined2D->fitTo(*data, Save(), NumCPU(10, 0)) ;
	
/*
	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	data->plotOn(xframe_mbc);
	combined2D->plotOn(xframe_mbc,Components("signalMBC"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
	combined2D->plotOn(xframe_mbc, LineColor(kBlue), ProjWData(*data));
	//combined2D->paramOn(xframe_mbc, ProjWData(*data));	
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
	combined2D->plotOn(xframe_de, Components("signalDE"), LineColor(kRed), LineStyle(kDashed), ProjWData(*data));
	combined2D->plotOn(xframe_de, LineColor(kBlue), ProjWData(*data));
	//combined2D->paramOn(xframe_de, ProjWData(*data));
	RooHist* pull_de = xframe_de->pullHist() ;
        RooPlot* f_pull_de = de.frame(Title("")) ;
        f_pull_de->addPlotable(pull_de,"BX0") ;
        f_pull_de->SetMaximum(4);
        f_pull_de->SetMinimum(-4);
	f_pull_de->GetYaxis()->SetLabelSize(0.11);
        f_pull_de->GetXaxis()->SetLabelSize(0.11);	

	
	
	
	//Canvases and Pads
	TCanvas* can = new TCanvas("c","c",1300,800) ; can->Divide(2,1);
	TPad *pad1 = new TPad("pad1", "pad 1",0.0,0.0,0.5,1.0,0);
	TPad *pad2 = new TPad("pad2", "pad 2",0.5,0.0,1.0,1.0,0);
	pad1->Draw(); 	pad2->Draw();   
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
	can->SaveAs(canvasName.c_str());
*/	

	//Save fit results
	
	//To individual fit file
	RooArgSet(fullFit->floatParsFinal()).writeToFile(resultName.c_str());

	//To consolidated parameter result files
	ofstream fout; string confilename, fileSuffix;
	fileSuffix = "_" + smc_sample + ".txt";

	//Nsigscf
	confilename = "/Nsigscf"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsigscf.getVal()<<"\t"<<Nsigscf.getError()<<"\t"<<true_nSigSCF<<endl ;
        fout.close();
	
	//Nbkgtot
	confilename = "/Nbkgtot"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nbkgtot.getVal()<<"\t"<<Nbkgtot.getError()<<"\t"<<true_nBkgTot<<endl ;
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

