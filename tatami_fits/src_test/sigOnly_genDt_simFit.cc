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
#include "RooBCPGenDecay.h"
#include "TF1.h"
#include "RooCFunction1Binding.h"
//#include "RooCFunction3Binding.h"
#include "RooAbsReal.h"
#include "RooSimWSTool.h"

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
	

	tree->SetBranchAddress("issignal",&f_issignal);
        tree->SetBranchAddress("issgnevt",&f_issgnevt);
        tree->SetBranchAddress("dtgen",&f_dtgen);
        tree->SetBranchAddress("genbq",&f_genBFlavor);
	tree->SetBranchAddress("iflavor",&i_iflavor);
	tree->SetBranchAddress("irbin",&i_irbin);
	tree->SetBranchAddress("inexp",&i_inexp);

        RooRealVar dt("dt", "dt", -10.0, 10.0, "ps");

	//Flavor q-category
	RooCategory sample("sample","sample");
	sample.defineType("B0",+1);
	sample.defineType("B0b",-1);

	//Mistag w,dw-category
	RooCategory irbin("irbin", "irbin");
        for (int i = 1; i < 7; ++i){irbin.defineType(boost::lexical_cast<std::string > (i).c_str(), i);}

	//SVD category
	RooCategory inexp("inexp", "inexp");
        inexp.defineType("7",7);
        inexp.defineType("31",31);	

	RooArgSet Variables(dt);
	RooArgSet ConditionalVariables(sample, irbin, inexp);
	Variables.add(ConditionalVariables);
	RooDataSet * data = new RooDataSet("data", "data", Variables);


	//Get data from tree
        for(int i=0; i<tree->GetEntries(); i++){

        tree->GetEntry(i);

        if(f_issignal == 1 && abs(f_dtgen) < 10.0 && i_irbin > 0){

	dt.setVal(f_dtgen);
	inexp.setIndex(i_inexp);
	irbin.setIndex(i_irbin);
        sample.setIndex(int(f_genBFlavor)); //gen-flavor
	//sample.setIndex(int(i_iflavor)); //FT-flavor	

	data->add(Variables);

        }

        }//tree-reading loop

//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//----------------------------------------------------------------------------------------------//




	//Define PDF
  RooRealVar dm("#Delta m","dm",0.507);//                                      
  RooRealVar A("A","A", 0.0, -1.0, 1.0);//                                     
  RooRealVar S("S","S", 0.0, -1.0, 1.0);//                                     
  RooRealVar tau("tau","tau", 1.525);//, 1.3, 1.7);//                          
  RooRealVar w("w","w", 0.0);
  RooRealVar dw("dw","dw", 0.0);
  RooRealVar mu("mu","mu", 0.0);//, -1.0, 1.0);


  RooTruthModel deltaFunc("deltaFunc", "deltaFunc", dt);
  RooBCPGenDecay Model("Model", "Model", dt, sample, tau, dm, w, A, S, dw, mu, deltaFunc);

 RooWorkspace rws("rws", "rws");
 rws.import(RooArgSet(Model, inexp, irbin));
 //rws.import(RooArgSet(Model, inexp));
 RooSimWSTool sct(rws);
 RooSimultaneous *SimPdf = sct.build("SimPdf", "Model", SplitParam("w,dw", "irbin,inexp"));
 //RooSimultaneous *SimPdf = sct.build("SimPdf", "Model", SplitParam("S", "inexp"));
 RooAbsPdf *simPdf = (RooAbsPdf*)SimPdf;
 
/*
//Set Belle-FT Mistag values
//w
rws.var("w_{1;7}")->setVal(0.420827); rws.var("w_{1;31}")->setVal(0.412222);
rws.var("w_{2;7}")->setVal(0.300296); rws.var("w_{2;31}")->setVal(0.307838);
rws.var("w_{3;7}")->setVal(0.219317); rws.var("w_{3;31}")->setVal(0.212765);
rws.var("w_{4;7}")->setVal(0.154636); rws.var("w_{4;31}")->setVal(0.149933);
rws.var("w_{5;7}")->setVal(0.0916131); rws.var("w_{5;31}")->setVal(0.0913264);
rws.var("w_{6;7}")->setVal(0.0228891); rws.var("w_{6;31}")->setVal(0.0218754);
//dw
rws.var("dw_{1;7}")->setVal(0.0583019); rws.var("dw_{1;31}")->setVal(0.00408778);
rws.var("dw_{2;7}")->setVal(0.00573998); rws.var("dw_{2;31}")->setVal(0.010326);
rws.var("dw_{3;7}")->setVal(-0.0392635); rws.var("dw_{3;31}")->setVal(-0.00479522);
rws.var("dw_{4;7}")->setVal(0.00474508); rws.var("dw_{4;31}")->setVal(0.00151989);
rws.var("dw_{5;7}")->setVal(-0.0118737); rws.var("dw_{5;31}")->setVal(0.0143633);
rws.var("dw_{6;7}")->setVal(-0.00585326); rws.var("dw_{6;31}")->setVal(0.00189979);
*/


//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult* fitRes = simPdf->fitTo(*data, ConditionalObservables(ConditionalVariables), Save());

const float afix = 0.5;
TF1 *f1 = new TF1("myfunc","(1-[2])*(([0]*sin(0.507*x)) + ([1]*cos(0.507*x)))",-10, 10);

RooArgSet pars(S,A,mu);
//RooAbsReal *osc = bindFunction(f1, dt, pars);

	//dt PLOTING
	RooPlot *xframe_1 =dt.frame(Bins(50),Title("B^{0}"));
  data->plotOn(xframe_1,Cut("sample==sample::B0"), MarkerColor(kRed));
    data->plotOn(xframe_1,Cut("sample==sample::B0b"), MarkerColor(kBlue));
  simPdf->plotOn(xframe_1,Slice(sample,"B0"),ProjWData(sample,*data),LineColor(kRed));
    simPdf->plotOn(xframe_1,Slice(sample,"B0b"),ProjWData(sample,*data),LineColor(kBlue));
  simPdf->paramOn(xframe_1);



//   RooPlot *xframe_3 = dt.frame(Bins(50),Title("Asymmetry"));
//   combData->plotOn(xframe_3, Asymmetry(sample));
//   osc->plotOn(xframe_3);

TCanvas* can = new TCanvas("c","c", 900, 1200) ;
//  can->Divide(1,2) ;
	
  can->cd() ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
//  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;

	//Save fit results

	can->SaveAs(canvasName.c_str());

	//To individual fit file
	RooArgSet(fitRes->floatParsFinal()).writeToFile(resultName.c_str());	

	//To consolidated parameter result files
	ofstream fout; string confilename, fileSuffix;
        fileSuffix = "_" + smc_sample + ".txt";
	//A
	confilename = "/A"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        //fout<<toyindex<<"\t"<<A.getVal()<<"\t"<<A.getError()<<"\t"<<0.0<<endl ;
        fout<<toyindex<<"\t"<<rws.var("A")->getVal()<<"\t"<<rws.var("A")->getError()<<"\t"<<0.0<<endl ;
        fout.close();
        //S
        confilename = "/S"; confilename = ConResDir + confilename + fileSuffix;
        fout.open(confilename, ios::app);
        //fout<<toyindex<<"\t"<<S.getVal()<<"\t"<<S.getError()<<"\t"<<0.0<<endl ;
	fout<<toyindex<<"\t"<<rws.var("S")->getVal()<<"\t"<<rws.var("S")->getError()<<"\t"<<0.0<<endl ;
        fout.close();


//----------------------------------------------------------------------------------------------//		

        return 0;
}







