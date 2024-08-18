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
#include "RooAbsData.h"
#include "TH1D.h"

using namespace std;
using namespace RooFit;

string convertToString(char* a)
{
    string s = a;
    return s;
}



RooDataSet * getRandomizedDataSet(char *parname, double mu, double sigL, double sigR, int numPoints){

RooRealVar par(parname, parname, -10.0, 10.0);
//Bifurcated Gaussians
RooRealVar mean("mean", "mean", 0.0, -10.0, 10.0); mean.setConstant(true);
RooRealVar sigmaL("sigmaL", "sigmaL", 1.0, -10.0, 10.0); sigmaL.setConstant(true);
RooRealVar sigmaR("sigmaR", "sigmaR", 1.0, -10.0, 10.0); sigmaR.setConstant(true);
mean.setVal(mu); sigmaL.setVal(sigL); sigmaR.setVal(sigR);
RooBifurGauss Gauss("Gauss", "Gauss", par, mean, sigmaL, sigmaR) ;
RooDataSet* data = Gauss.generate(par, numPoints, Name("data")) ;

return data;
}

void plot_RandomizedParameter(char *parname, double mu, double sigL, double sigR, int numPoints, RooDataSet* data){

int nbins = 40;
double uLim = mu + (4*sigR); double lLim = mu - (4*sigL); 

TH1D hist(parname, parname, nbins, lLim, uLim);

for(int i = 0; i< numPoints; i++){
hist.Fill( ((RooRealVar*) (data->get(i)->find(parname)) )->getVal() );
}

//TCanvas *can = new TCanvas(parname, parname, 800, 800);
TCanvas *can = new TCanvas("can", "can", 800, 800);
can->cd();
hist.Draw("E0");
string Parname = convertToString(parname);
string plotname = "RandomizedParameters/deltaMtauB/" + Parname +  ".png";

can->SaveAs(plotname.c_str());
can->Delete(); 
hist.Delete();

}


// Main method
int main() {

gROOT->SetBatch(kTRUE);

//Generate random DeltaM and tauB
//w
RooDataSet* data = getRandomizedDataSet("DeltaM", 0.507, 0.0019, 0.0019, 100000);
RooDataSet* data_tauB = getRandomizedDataSet("tauB", 1.534, 0.004, 0.004, 100000);

 
//Merge datasets
data->merge(data_tauB); 

RooWorkspace *w_randPhyPars = new RooWorkspace("rws_randPhyPars", "rws_randPhyPars");
w_randPhyPars->import(*data);
w_randPhyPars->writeToFile("RandomizedParameters/deltaMtauB/randPhyPars.root");


plot_RandomizedParameter("DeltaM", 0.507, 0.0019, 0.0019, 100000, data);
plot_RandomizedParameter("tauB", 1.534, 0.004, 0.004, 100000, data);

return 0;

}

















