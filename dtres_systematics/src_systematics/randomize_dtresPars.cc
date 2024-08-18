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
#include "TRandom.h"


using namespace std;
using namespace RooFit;

string convertToString(char* a)
{
    string s = a;
    return s;
}


void get_err(int dtres_item, int dtres_subitem, int expno, double& errL, double& errR, double&val0, string& parName){

double val1, val2;

val0 = get_dtres_parVal_b0(dtres_item, dtres_subitem, expno, 1, parName);
dtres_systematics_b0(dtres_item, 1.0, expno, 1);
val1 = get_dtres_parVal_b0(dtres_item, dtres_subitem, expno, 1, parName);
dtres_systematics_b0(dtres_item, -1.0, expno, 1);
val2 = get_dtres_parVal_b0(dtres_item, dtres_subitem, expno, 1, parName);
errR = val1 - val0;
errL = val1 - val2;
}


RooDataSet* RandomizedParameter(int dtres_item, int dtres_subitem, int numPtsData, int numPtsPlot, int expno, TRandom * rand){

string parName; double sigma, errL, errR, val0, parVal;
get_err(dtres_item, dtres_subitem, expno, errL, errR, val0, parName);
int nbins = 40;
double uLim = val0 + (4*errR);
double lLim = val0 - (4*errL);

parName = parName + "_exp" + to_string(expno);
TH1D hist(parName.c_str(), parName.c_str(), nbins, lLim, uLim);

int numPtsPlotR, numPtsPlotL;
if(errR != 0 && errL != 0){
numPtsPlotR = round((float(numPtsPlot) * errR) / (errR + errL));
numPtsPlotL = round((float(numPtsPlot) * errL) / (errR + errL));
}
else{
numPtsPlotR = numPtsPlot/2; numPtsPlotL = numPtsPlot/2;
}
//cout<<"numPtsPlotR = "<<numPtsPlotR<<"\t numPtsPlotL = "<<numPtsPlotL<<endl;

//for(int i = 0; i< numPtsPlotR; i++){
int countR = 0; while (countR < numPtsPlotR){
sigma = rand->Gaus(0,1); 
if(sigma > 0) {parVal = val0 +(sigma *errR); hist.Fill( parVal ); countR++;}
}
//for(int i = 0; i< numPtsPlotL; i++){
int countL = 0; while (countL < numPtsPlotL){
sigma = rand->Gaus(0,1);
if(sigma < 0) {parVal = val0 + (sigma *errL); hist.Fill( parVal ); countL++;}
}

//TCanvas *can = new TCanvas(parname, parname, 800, 800);
TCanvas *can = new TCanvas("can", "can", 800, 800);
can->cd();
hist.Draw("E0");
string plotname = "RandomizedParameters/dtresPars/" + parName + ".png";
can->SaveAs(plotname.c_str());
can->Delete(); hist.Delete();

//dataset
string parname = "dtres_item" + to_string(dtres_item)+ "_" + to_string(dtres_subitem) + "_exp" + to_string(expno);
RooRealVar par(parname.c_str(), parname.c_str(), -10.0, 10.0);
RooDataSet* data = new RooDataSet("data","data",RooArgSet(par)) ;
int numPtsDataR, numPtsDataL;
if(errR != 0 && errL != 0){
numPtsDataR = round((float(numPtsData)* errR) / (errR + errL));
numPtsDataL = round((float(numPtsData)* errL) / (errR + errL));
}
else{
numPtsDataR = numPtsData/2; numPtsDataL = numPtsData/2;
}
//cout<<"numPtsDataR  = "<<numPtsDataR<<"\t numPtsDataL = "<<numPtsDataL<<endl;
//for(int i = 0; i < numPtsDataR; i++){
countR = 0; while (countR < numPtsDataR){
sigma = rand->Gaus(0,1);
if(sigma > 0) {par.setVal(sigma); data->add(RooArgSet(par)); countR++;}
}
//for(int i = 0; i < numPtsDataL; i++){
countL = 0; while (countL < numPtsDataL){
sigma = rand->Gaus(0,1);
if(sigma <0) {par.setVal(sigma); data->add(RooArgSet(par)); countL++;}
}

return data;
}


// Main method
int main() {

gROOT->SetBatch(kTRUE);


TRandom *rand = new TRandom();

RooDataSet *data = RandomizedParameter(0, 0, 1000, 100000, 29, rand);
data->merge(RandomizedParameter(0, 0, 1000, 100000, 31, rand));

for(int i = 1; i <= 38; i++ ){
data->merge(RandomizedParameter(i, 0, 1000, 100000, 29, rand));
data->merge(RandomizedParameter(i, 0, 1000, 100000, 31, rand));
}
//Deal with items with subitems
//29
data->merge(RandomizedParameter(6, 1, 1000, 100000, 29, rand));
data->merge(RandomizedParameter(7, 1, 1000, 100000, 29, rand));
data->merge(RandomizedParameter(8, 1, 1000, 100000, 29, rand));
data->merge(RandomizedParameter(9, 1, 1000, 100000, 29, rand));
data->merge(RandomizedParameter(11, 1, 1000, 100000, 29, rand));
data->merge(RandomizedParameter(38, 1, 1000, 100000, 29, rand));
///31
data->merge(RandomizedParameter(6, 1, 1000, 100000, 31, rand));
data->merge(RandomizedParameter(7, 1, 1000, 100000, 31, rand));
data->merge(RandomizedParameter(8, 1, 1000, 100000, 31, rand));
data->merge(RandomizedParameter(9, 1, 1000, 100000, 31, rand));
data->merge(RandomizedParameter(11, 1, 1000, 100000, 31, rand));
data->merge(RandomizedParameter(38, 1, 1000, 100000, 31, rand));




RooWorkspace *w_rand_dtresPars = new RooWorkspace("rws_rand_dtresPars", "rws_rand_dtresPars");
w_rand_dtresPars->import(*data);
w_rand_dtresPars->writeToFile("RandomizedParameters/dtresPars/rand_DTresPars.root");


return 0;

}

















