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

	//Get input A and S
	float A_input = opt->GetAcp_Ainput(); cout<<"A_input = "<<A_input<<endl;
	float S_input = opt->GetAcp_Sinput(); cout<<"S_input = "<<S_input<<endl;

	//Get the input file (one file for now!)
	string filename = opt->GetFilenames().at(0);

	//Get the names of the output plot, result file, consolidated result directory and toy-fit index
	string canvasName = opt->GetPlotname();
	string resultName = opt->GetResultname();
	string ConResDir = opt->GetConsolidatedResultsDir();
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
	
	//Event weight to induce an artificial asymmetry
	RooRealVar acp_weight("acp_weight","acp_weight", 0, 2);

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

        RooArgSet mainVariables(mbc, de, dt, dalitzCategory, acp_weight);
        mainVariables.add(AllConditionalVariables);

	// Read tree in a loop and add to a roodataset
	int true_Nsig = 0, true_Ncont = 0, true_Nrare = 0, true_Nscf = 0, true_Nbb = 0;
	 
	RooDataSet * data = new RooDataSet("data", "data", mainVariables);
	for(int i=0; i<tree->GetEntries(); i++){

	tree->GetEntry(i);

	//Remove buggy events; To be invertigated later, why this happens
 	if(f_dtgen == 0.0) continue;

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
	w.setVal(f_w);
	dw.setVal(f_dw);
	ecms.setVal(f_ecms);
	irbin.setIndex(i_irbin);
	issignal.setVal(f_issignal);
	issgnevt.setVal(f_issgnevt);

	//Calculate acp_weight
	if(f_issgnevt == 1){
	cout<<"----------------------------------------------"<<endl;
	float acp_w = 1.0 + (f_genBFlavor * ((A_input * cos(delta_m*f_dtgen)) + (S_input * sin(delta_m*f_dtgen))));
	acp_weight.setVal(acp_w);
	//cout<<"genBFlavor = "<<f_genBFlavor<<endl;
	//cout<<"A_input = "<<A_input<<endl;
	//cout<<"S_input = "<<S_input<<endl;
	//cout<<"delta_m = "<<delta_m<<endl;
	//cout<<"f_dtgen = "<<f_dtgen<<endl;
	//cout<<"acp_weight.getVal() = "<<acp_weight.getVal()<<endl;
	//cout<<"f_acp_weight = "<<acp_w<<endl;
	}
	else{
	acp_weight.setVal(1.0);
	}

	data->add(mainVariables);


	}//tree-reading loop

	// weighted dataset
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,acp_weight.GetName()) ;

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

	//Import dt pdf's
        rws.import("pdf_models/fix_shape/cont_shape.root:continuum:contDT");
        rws.import("pdf_models/fix_shape/bb_shape.root:bb:bbDT");
        rws.import("pdf_models/fix_shape/rarebkg_shape.root:rare:rareDT");
        rws.import("pdf_models/fix_shape/scf_shape.root:scf:scfDT");
        rws.import("pdf_models/fix_shape/signal_shape.root:signal:signalDT");

	
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
        RooProdPdf * signal3D =  new RooProdPdf("signal","signal",RooArgSet(*signalMBC, *signalDE, *signalDT));

	//Component yields
	RooRealVar Nsig("N_{sig}","Signal yield", 840, 400, 1200);
        RooRealVar Ncont("N_{cont}","Continuum background yield",4650, 3500, 6000);
        RooRealVar Nrare("N_{rare}","Rare background yield",200, 0, 1000);
        RooRealVar Nscf("N_{scf}","Self cross feed yield",90, 0, 1000);
        RooRealVar Nbb("N_{bb}","B-Bbar background yield",300, 0, 800);

	//Full background model, and combined overall model	
	RooAddPdf * combined3D = new RooAddPdf("combined","combined", RooArgList(*signal3D, *cont3D, *bb3D, *rare3D, *scf3D), RooArgList(Nsig, Ncont, Nbb, Nrare, Nscf));

//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult * fullFit = combined3D->fitTo(wdata,ConditionalObservables(AllConditionalVariables), Save(), NumCPU(10, 0), SumW2Error(kTRUE)) ;
	

	//Plotting
	//Mbc
	RooPlot *xframe_mbc = mbc.frame(Bins(40));
	wdata.plotOn(xframe_mbc);
	combined3D->plotOn(xframe_mbc,Components("signalMBC"), LineColor(kRed), LineStyle(kDashed), ProjWData(wdata));
	combined3D->plotOn(xframe_mbc, LineColor(kBlue), ProjWData(wdata));
	//combined3D->paramOn(xframe_mbc, ProjWData(wdata));	
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
	combined3D->plotOn(xframe_de, Components("signalDE"), LineColor(kRed), LineStyle(kDashed), ProjWData(wdata));
	combined3D->plotOn(xframe_de, LineColor(kBlue), ProjWData(wdata));
	//combined3D->paramOn(xframe_de, ProjWData(wdata));
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
	combined3D->plotOn(xframe_dt,Components("signalDT"), LineColor(kRed), LineStyle(kDashed), ProjWData(wdata));
        combined3D->plotOn(xframe_dt, LineColor(kBlue), ProjWData(wdata));
        //combined3D->paramOn(xframe_dt, ProjWData(wdata));
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
	


	//Save fit results
	
	//To individual fit file
	RooArgSet(fullFit->floatParsFinal()).writeToFile(resultName.c_str());

	//To consolidated parameter result files
	ofstream fout; string confilename;
	//A
	confilename = "/A.txt"; confilename = ConResDir + confilename;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<rws.var("A")->getVal()<<"\t"<<rws.var("A")->getError()<<"\t"<<A_input<<endl ;
        fout.close();
        //S
        confilename = "/S.txt"; confilename = ConResDir + confilename;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<rws.var("S")->getVal()<<"\t"<<rws.var("S")->getError()<<"\t"<<S_input<<endl ;
        fout.close();
        //Nsig
        confilename = "/Nsig.txt"; confilename  = ConResDir + confilename;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nsig.getVal()<<"\t"<<Nsig.getError()<<"\t"<<true_Nsig<<endl ;
        fout.close();
        //Ncont
        confilename = "/Ncont.txt"; confilename  = ConResDir + confilename;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Ncont.getVal()<<"\t"<<Ncont.getError()<<"\t"<<true_Ncont<<endl ;
        fout.close();
        //Nscf
        confilename = "/Nrare.txt"; confilename  = ConResDir + confilename;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nrare.getVal()<<"\t"<<Nrare.getError()<<"\t"<<true_Nrare<<endl ;
        fout.close();
        //Nscf
        confilename = "/Nscf.txt"; confilename  = ConResDir + confilename;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nscf.getVal()<<"\t"<<Nscf.getError()<<"\t"<<true_Nscf<<endl ;
        fout.close();
        //Nbb
        confilename = "/Nbb.txt"; confilename  = ConResDir + confilename;
        fout.open(confilename, ios::app);
        fout<<toyindex<<"\t"<<Nbb.getVal()<<"\t"<<Nbb.getError()<<"\t"<<true_Nbb<<endl ;
        fout.close();

//----------------------------------------------------------------------------------------------//		

        return 0;
}

