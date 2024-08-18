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

        RooRealVar dt("dt", "dt", -10.0, 10.0, "ps");
        RooDataSet * data_p = new RooDataSet("data_p", "data_p", RooArgSet(dt));
        RooDataSet * data_n = new RooDataSet("data_n", "data_n", RooArgSet(dt));

        for(int i=0; i<tree->GetEntries(); i++){

        tree->GetEntry(i);

        dt.setVal(f_dtgen);

        if(f_issignal == 1){

        if(f_genBFlavor == 1)   data_p->add(RooArgSet(dt));
        if(f_genBFlavor == -1)   data_n->add(RooArgSet(dt));

        }

        }//tree-reading loop

//----------------------------------------------------------------------------------------------//	
//                                  FIT MODELLING                                               //
//----------------------------------------------------------------------------------------------//




	//Import dt pdf's
  RooRealVar dm("#Delta m","dm",0.507);//                                      
  RooRealVar A("A","A", 0.0, -1.0, 1.0);//                                     
  RooRealVar S("S","S", 0.0, -1.0, 1.0);//                                     
  RooRealVar tau("tau","tau", 1.525);//, 1.3, 1.7);//                          
  RooRealVar avgMistag("avgMistag","avgMistag", 0.0);
  RooRealVar delMistag("delMistag","delMistag", 0.0);
  RooRealVar mu("mu","mu", 0.0);//, -1.0, 1.0);

//CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  RooCategory sample("sample","sample");
  sample.defineType("B0",+1);
  sample.defineType("B0b",-1);

  RooTruthModel deltaFunc("deltaFunc", "deltaFunc", dt);
  RooBCPGenDecay Model("Model", "Model", dt, sample, tau, dm, avgMistag, A, S, delMistag, mu, deltaFunc);


RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dt),Index(sample),Import("B0",*data_p),Import("B0b",*data_n)); 




//----------------------------------------------------------------------------------------------//
//                                 FITTING, PLOTS and RESULTS                                   //
//----------------------------------------------------------------------------------------------//

	
	//Fitting
	RooFitResult* fitRes = Model.fitTo(*combData,Save());

const float afix = 0.5;
TF1 *f1 = new TF1("myfunc","(1-[2])*(([0]*sin(0.507*x)) + ([1]*cos(0.507*x)))",-10, 10);

RooArgSet pars(S,A,mu);
//RooAbsReal *osc = bindFunction(f1, dt, pars);

	//dt PLOTING
	RooPlot *xframe_1 =dt.frame(Bins(50),Title("B^{0}"));
  combData->plotOn(xframe_1,Cut("sample==sample::B0"), MarkerColor(kRed));
    combData->plotOn(xframe_1,Cut("sample==sample::B0b"), MarkerColor(kBlue));
  Model.plotOn(xframe_1,Slice(sample,"B0"),ProjWData(sample,*combData),LineColor(kRed));
    Model.plotOn(xframe_1,Slice(sample,"B0b"),ProjWData(sample,*combData),LineColor(kBlue));
  Model.paramOn(xframe_1,data_p);



//   RooPlot *xframe_3 = dt.frame(Bins(50),Title("Asymmetry"));
//   combData->plotOn(xframe_3, Asymmetry(sample));
//   osc->plotOn(xframe_3);

TCanvas* can = new TCanvas("c","c", 900, 1200) ;
//  can->Divide(1,2) ;
	
  can->cd() ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
//  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;

	//Save fit results

	can->SaveAs(canvasName.c_str());

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







