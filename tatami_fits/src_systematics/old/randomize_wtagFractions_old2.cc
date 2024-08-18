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

using namespace std;
using namespace RooFit;


//RooDataSet * getRandomizedDataSet(RooRealVar par, double mu, double sigL, double sigR, int numPoints){
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


// Main method
int main() {

//Roorealvars for w-tag metrics
//SVD-1
RooRealVar w10("w10", "w10", -2.0, 2.0); RooRealVar dw10("dw10", "dw10", -2.0, 2.0);
RooRealVar w20("w20", "w20", -2.0, 2.0); RooRealVar dw20("dw20", "dw20", -2.0, 2.0);
RooRealVar w30("w30", "w30", -2.0, 2.0); RooRealVar dw30("dw30", "dw30", -2.0, 2.0);
RooRealVar w40("w40", "w40", -2.0, 2.0); RooRealVar dw40("dw40", "dw40", -2.0, 2.0);
RooRealVar w50("w50", "w50", -2.0, 2.0); RooRealVar dw50("dw50", "dw50", -2.0, 2.0);
RooRealVar w60("w60", "w60", -2.0, 2.0); RooRealVar dw60("dw60", "dw60", -2.0, 2.0);
//SVD-2
RooRealVar w11("w11", "w11", -2.0, 2.0); RooRealVar dw11("dw11", "dw11", -2.0, 2.0);
RooRealVar w21("w21", "w21", -2.0, 2.0); RooRealVar dw21("dw21", "dw21", -2.0, 2.0);
RooRealVar w31("w31", "w31", -2.0, 2.0); RooRealVar dw31("dw31", "dw31", -2.0, 2.0);
RooRealVar w41("w41", "w41", -2.0, 2.0); RooRealVar dw41("dw41", "dw41", -2.0, 2.0);
RooRealVar w51("w51", "w51", -2.0, 2.0); RooRealVar dw51("dw51", "dw51", -2.0, 2.0);
RooRealVar w61("w61", "w61", -2.0, 2.0); RooRealVar dw61("dw61", "dw61", -2.0, 2.0);


RooDataSet* data = getRandomizedDataSet("w10", 0.420827, 0.006001569, 0.007235697, 1000);
RooDataSet* data_w20 = getRandomizedDataSet("w20", 0.300296, 0.006430566, 0.007129388, 1000);
/*
//Since the errors provided on the w-tag metrics are asymmetric, we randomize them using
//an asymmetric gaussian.
RooRealVar mean("mean", "mean", 0.0, -2.0, 2.0); mean.setConstant(true);
RooRealVar sigmaL("sigmaL", "sigmaL", 1.0, -2.0, 2.0); sigmaL.setConstant(true); 
RooRealVar sigmaR("sigmaR", "sigmaR", 1.0, -2.0, 2.0); sigmaR.setConstant(true);

//Generate random w-tags : Values taken from https://belle.kek.jp/group/indirectcp/cpfit/2010/wtag_mdlh_2010
//SVD1
//w
mean.setVal(0.420827); sigmaL.setVal(0.006001569); sigmaR.setVal(0.007235697); 
RooBifurGauss Gaussw10("Gaussw10", "Gaussw10", w10, mean, sigmaL, sigmaR) ; RooDataSet* data = Gaussw10.generate(w10, 1000, Name("data")) ; 

mean.setVal(0.300296); sigmaL.setVal(0.006430566); sigmaR.setVal(0.007129388);        
RooBifurGauss Gaussw20("Gaussw20", "Gaussw20", w20, mean, sigmaL, sigmaR) ; RooDataSet* data_w20 = Gaussw20.generate(w20, 1000) ;
     
mean.setVal(0.219317); sigmaL.setVal(0.007693083); sigmaR.setVal(0.007417778);             
RooBifurGauss Gaussw30("Gaussw30", "Gaussw30", w30, mean, sigmaL, sigmaR) ; RooDataSet* data_w30 = Gaussw30.generate(w30, 1000) ;

mean.setVal(0.154636); sigmaL.setVal(0.006416449); sigmaR.setVal(0.006885875);             
RooBifurGauss Gaussw40("Gaussw40", "Gaussw40", w40, mean, sigmaL, sigmaR) ; RooDataSet* data_w40 = Gaussw40.generate(w40, 1000) ;

mean.setVal(0.0916131); sigmaL.setVal(0.008807757); sigmaR.setVal(0.006761047); 
RooBifurGauss Gaussw50("Gaussw50", "Gaussw50", w50, mean, sigmaL, sigmaR) ; RooDataSet* data_w50 = Gaussw50.generate(w50, 1000) ;

mean.setVal(0.0228891); sigmaL.setVal(0.004587614); sigmaR.setVal(0.004336734); 
RooBifurGauss Gaussw60("Gaussw60", "Gaussw60", w60, mean, sigmaL, sigmaR) ; RooDataSet* data_w60 = Gaussw60.generate(w60, 1000) ;

//dw
mean.setVal(0.0583019); sigmaL.setVal(0.009165492); sigmaR.setVal(0.008923343); RooDataSet* data_dw10 = Gauss.generate(dw10, 1000) ;            


mean.setVal(0.00573998); sigmaL.setVal(0.009147374); sigmaR.setVal(0.009180684); RooDataSet* data_dw20 = Gauss.generate(dw20, 1000) ;
mean.setVal(-0.0392635); sigmaL.setVal(0.009988773); sigmaR.setVal(0.001035175); RooDataSet* data_dw30 = Gauss.generate(dw30, 1000) ;
mean.setVal(0.00474508); sigmaL.setVal(0.008879500); sigmaR.setVal(0.009020215); RooDataSet* data_dw40 = Gauss.generate(dw40, 1000) ;            
mean.setVal(-0.0118737); sigmaL.setVal(0.009294976); sigmaR.setVal(0.009308208); RooDataSet* data_dw50 = Gauss.generate(dw50, 1000) ;
mean.setVal(-0.00585326); sigmaL.setVal(0.005787974); sigmaR.setVal(0.005693457); RooDataSet* data_dw60 = Gauss.generate(dw60, 1000) ;
//SVD2
//w
mean.setVal(0.412222); sigmaL.setVal(0.003577812); sigmaR.setVal(0.004152612); RooDataSet* data_w11 = Gauss.generate(w11, 1000) ;            
mean.setVal(0.307838); sigmaL.setVal(0.002803811); sigmaR.setVal(0.003243236); RooDataSet* data_w21 = Gauss.generate(w21, 1000) ;
mean.setVal(0.212765); sigmaL.setVal(0.003486607); sigmaR.setVal(0.003721417); RooDataSet* data_w31 = Gauss.generate(w31, 1000) ;
mean.setVal(0.149933); sigmaL.setVal(0.004241595); sigmaR.setVal(0.003315138); RooDataSet* data_w41 = Gauss.generate(w41, 1000) ;
mean.setVal(0.0913264); sigmaL.setVal(0.003696399); sigmaR.setVal(0.003180302); RooDataSet* data_w51 = Gauss.generate(w51, 1000) ;
mean.setVal(0.0218754); sigmaL.setVal(0.003077622); sigmaR.setVal(0.002175087); RooDataSet* data_w61 = Gauss.generate(w61, 1000) ;
//dw
mean.setVal(0.00408778); sigmaL.setVal(0.003927049); sigmaR.setVal(0.003961548); RooDataSet* data_dw11 = Gauss.generate(dw11, 1000) ;
mean.setVal(0.010326); sigmaL.setVal(0.003698619); sigmaR.setVal(0.003543129); RooDataSet* data_dw21 = Gauss.generate(dw21, 1000) ;
mean.setVal(-0.00479522); sigmaL.setVal(0.004179366); sigmaR.setVal(0.004129422); RooDataSet* data_dw31 = Gauss.generate(dw31, 1000) ;
mean.setVal(0.00151989); sigmaL.setVal(0.004602366); sigmaR.setVal(0.004169570); RooDataSet* data_dw41 = Gauss.generate(dw41, 1000) ;
mean.setVal(0.0143633); sigmaL.setVal(0.003914627); sigmaR.setVal(0.003998982); RooDataSet* data_dw51 = Gauss.generate(dw51, 1000) ;
mean.setVal(0.00189979); sigmaL.setVal(0.002360543); sigmaR.setVal(0.002433324); RooDataSet* data_dw61 = Gauss.generate(dw61, 1000) ;
*/

//for(int i = 1; i<= 50; i++){
//RooRealVar * par = (RooRealVar*) (data_w->get(i)->find("w10"));
//cout<<par->getVal()<<"\t"<<x.getVal()<<endl;
//}

//Merge datasets
data->merge(data_w20);
/*
data->merge(data_w30);
data->merge(data_w40);
data->merge(data_w50);
data->merge(data_w60);

data->merge(data_w11);
data->merge(data_w21);
data->merge(data_w31);
data->merge(data_w41);
data->merge(data_w51);
data->merge(data_w61);

data->merge(data_dw10);
data->merge(data_dw20);
data->merge(data_dw30);
data->merge(data_dw40);
data->merge(data_dw50);
data->merge(data_dw60);

data->merge(data_dw11);
data->merge(data_dw21);
data->merge(data_dw31);
data->merge(data_dw41);
data->merge(data_dw51);
data->merge(data_dw61);
*/



RooWorkspace *w_randWTags = new RooWorkspace("rws_randWTags", "rws_randWTags");
w_randWTags->import(*data);
w_randWTags->writeToFile("RandomizedParameters/WrongTags/randWTagMetrics.root");

for(int i = 1; i<= 50; i++){
RooRealVar * par1 = (RooRealVar*) (data->get(i)->find("w10"));
RooRealVar * par2 = (RooRealVar*) (data->get(i)->find("w20"));
cout<<par1->getVal()<<"\t"<<par2->getVal()<<endl;
}

RooPlot *xframe = w10.frame(Bins(40), Title(""));
data->plotOn(xframe);



TCanvas *can = new TCanvas("can", "can", 800, 800);
can->cd();
xframe->Draw() ; 
can->SaveAs("RandomizedParameters/WrongTags/w10.png");
        

return 0;

}

















