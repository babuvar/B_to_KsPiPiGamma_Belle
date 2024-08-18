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

RooRealVar par(parname, parname, -2.0, 2.0);
//Bifurcated Gaussians
RooRealVar mean("mean", "mean", 0.0, -2.0, 2.0); mean.setConstant(true);
RooRealVar sigmaL("sigmaL", "sigmaL", 1.0, -2.0, 2.0); sigmaL.setConstant(true);
RooRealVar sigmaR("sigmaR", "sigmaR", 1.0, -2.0, 2.0); sigmaR.setConstant(true);
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
string plotname = "RandomizedParameters/WrongTags/" + Parname +  ".png";

can->SaveAs(plotname.c_str());
can->Delete(); 
hist.Delete();

}


// Main method
int main() {

gROOT->SetBatch(kTRUE);

//Generate random w-tags : Values taken from https://belle.kek.jp/group/indirectcp/cpfit/2010/wtag_mdlh_2010
//SVD1
//w
RooDataSet* data = getRandomizedDataSet("w10", 0.420827, 0.006001569, 0.007235697, 100000);
RooDataSet* data_w20 = getRandomizedDataSet("w20", 0.300296, 0.006430566, 0.007129388, 100000);
RooDataSet* data_w30 = getRandomizedDataSet("w30", 0.219317, 0.007693083, 0.007417778, 100000);             
RooDataSet* data_w40 = getRandomizedDataSet("w40", 0.154636, 0.006416449, 0.006885875, 100000);             
RooDataSet* data_w50 = getRandomizedDataSet("w50", 0.0916131, 0.008807757, 0.006761047, 100000); 
RooDataSet* data_w60 = getRandomizedDataSet("w60", 0.0228891, 0.004587614, 0.004336734, 100000); 
//dw
RooDataSet* data_dw10 = getRandomizedDataSet("dw10", 0.0583019, 0.009165492, 0.008923343, 100000);             
RooDataSet* data_dw20 = getRandomizedDataSet("dw20", 0.00573998, 0.009147374, 0.009180684, 100000);  
RooDataSet* data_dw30 = getRandomizedDataSet("dw30", -0.0392635, 0.009988773, 0.01035175, 100000);  
RooDataSet* data_dw40 = getRandomizedDataSet("dw40", 0.00474508, 0.008879500, 0.009020215, 100000);             
RooDataSet* data_dw50 = getRandomizedDataSet("dw50", -0.0118737, 0.009294976, 0.009308208, 100000);  
RooDataSet* data_dw60 = getRandomizedDataSet("dw60", -0.00585326, 0.005787974, 0.005693457, 100000);  
//SVD2
//w
RooDataSet* data_w11 = getRandomizedDataSet("w11", 0.412222, 0.003577812, 0.004152612, 100000);             
RooDataSet* data_w21 = getRandomizedDataSet("w21", 0.307838, 0.002803811, 0.003243236, 100000); 
RooDataSet* data_w31 = getRandomizedDataSet("w31", 0.212765, 0.003486607, 0.003721417, 100000);  
RooDataSet* data_w41 = getRandomizedDataSet("w41", 0.149933, 0.004241595, 0.003315138, 100000);  
RooDataSet* data_w51 = getRandomizedDataSet("w51", 0.0913264, 0.003696399, 0.003180302, 100000); 
RooDataSet* data_w61 = getRandomizedDataSet("w61", 0.0218754, 0.003077622, 0.002175087, 100000); 
//dw
RooDataSet* data_dw11 = getRandomizedDataSet("dw11", 0.00408778, 0.003927049, 0.003961548, 100000);  
RooDataSet* data_dw21 = getRandomizedDataSet("dw21", 0.010326, 0.003698619, 0.003543129, 100000);  
RooDataSet* data_dw31 = getRandomizedDataSet("dw31", -0.00479522, 0.004179366, 0.004129422, 100000);  
RooDataSet* data_dw41 = getRandomizedDataSet("dw41", 0.00151989, 0.004602366, 0.004169570, 100000);  
RooDataSet* data_dw51 = getRandomizedDataSet("dw51", 0.0143633, 0.003914627, 0.003998982, 100000);  
RooDataSet* data_dw61 = getRandomizedDataSet("dw61", 0.0018998, 0.002360543, 0.002433324, 100000);  

RooDataSet* data_dummy = getRandomizedDataSet("dummy", 0.2, 0.1, 0.1, 100000);

 
//Merge datasets
data->merge(data_w20); data->merge(data_w30); data->merge(data_w40); data->merge(data_w50); data->merge(data_w60);
data->merge(data_w11); data->merge(data_w21); data->merge(data_w31); data->merge(data_w41); data->merge(data_w51); data->merge(data_w61);
data->merge(data_dw10); data->merge(data_dw20); data->merge(data_dw30); data->merge(data_dw40); data->merge(data_dw50); data->merge(data_dw60);
data->merge(data_dw11); data->merge(data_dw21); data->merge(data_dw31); data->merge(data_dw41); data->merge(data_dw51); data->merge(data_dw61);

data->merge(data_dummy);

RooWorkspace *w_randWTags = new RooWorkspace("rws_randWTags", "rws_randWTags");
w_randWTags->import(*data);
w_randWTags->writeToFile("RandomizedParameters/WrongTags/randWTagMetrics.root");


plot_RandomizedParameter("w10", 0.420827, 0.006001569, 0.007235697, 100000, data);
plot_RandomizedParameter("w20", 0.300296, 0.006430566, 0.007129388, 100000, data);
plot_RandomizedParameter("w30", 0.219317, 0.007693083, 0.007417778, 100000, data);             
plot_RandomizedParameter("w40", 0.154636, 0.006416449, 0.006885875, 100000, data);             
plot_RandomizedParameter("w50", 0.0916131, 0.008807757, 0.006761047, 100000, data); 
plot_RandomizedParameter("w60", 0.0228891, 0.004587614, 0.004336734, 100000, data); 
plot_RandomizedParameter("dw10", 0.0583019, 0.009165492, 0.008923343, 100000, data);             
plot_RandomizedParameter("dw20", 0.00573998, 0.009147374, 0.009180684, 100000, data);  
plot_RandomizedParameter("dw30", -0.0392635, 0.009988773, 0.01035175, 100000, data);  
plot_RandomizedParameter("dw40", 0.00474508, 0.008879500, 0.009020215, 100000, data);             
plot_RandomizedParameter("dw50", -0.0118737, 0.009294976, 0.009308208, 100000, data);  
plot_RandomizedParameter("dw60", -0.00585326, 0.005787974, 0.005693457, 100000, data);  
plot_RandomizedParameter("w11", 0.412222, 0.003577812, 0.004152612, 100000, data);             
plot_RandomizedParameter("w21", 0.307838, 0.002803811, 0.003243236, 100000, data); 
plot_RandomizedParameter("w31", 0.212765, 0.003486607, 0.003721417, 100000, data);  
plot_RandomizedParameter("w41", 0.149933, 0.004241595, 0.003315138, 100000, data);  
plot_RandomizedParameter("w51", 0.0913264, 0.003696399, 0.003180302, 100000, data); 
plot_RandomizedParameter("w61", 0.0218754, 0.003077622, 0.002175087, 100000, data); 
plot_RandomizedParameter("dw11", 0.00408778, 0.003927049, 0.003961548, 100000, data);  
plot_RandomizedParameter("dw21", 0.010326, 0.003698619, 0.003543129, 100000, data);  
plot_RandomizedParameter("dw31", -0.00479522, 0.004179366, 0.004129422, 100000, data);  
plot_RandomizedParameter("dw41", 0.00151989, 0.004602366, 0.004169570, 100000, data);  
plot_RandomizedParameter("dw51", 0.0143633, 0.003914627, 0.003998982, 100000, data);  
plot_RandomizedParameter("dw61", 0.0018998, 0.002360543, 0.002433324, 100000, data);  

plot_RandomizedParameter("dummy", 0.2, 0.1, 0.1, 100000, data);

return 0;

}

















