#include "RooDataSet.h"

using namespace RooFit;

  RooRealVar dt("dt","dt",-10.0,10.00,"ps");
  RooDataSet* data = new RooDataSet("data","data",RooArgSet(dt));

int genRooLifetimeFit(){


string filename = "S0p0_1M_correct_toFit.root";

string Sample = filename;
Sample.erase(filename.length()-5);
string imagename = "plots_life/" + Sample + ".png";
string resultname = "results_life/" +  Sample + ".txt";
const char* Imagename = imagename.c_str();
const char* Resultname = resultname.c_str();
const char* Filename = filename.c_str();

gStyle->SetOptStat(0);



TFile * f = new TFile(Filename);


TTree * tree = (TTree *) f->Get("h1");

float nEntries = 0;

 
float f_dtgen, f_genbq, f_issgnevt;

tree->SetBranchAddress("dtgen",&f_dtgen);
tree->SetBranchAddress("issgnevt",&f_issgnevt);




for(int i = 0 ; i < tree->GetEntries(); i++ )
{


tree->GetEntry(i);
dt.setVal(f_dtgen);



if(f_issgnevt == 1){
if(fabs(f_dtgen) < 10 ){
 nEntries++;
 data->add(RooArgSet(dt));
}


}//fabs(f_dtgen) condition



} 

 

//PDF 
 
  RooRealVar tau("tau","tau", 1.525, 1.2, 2.0);//

  
  
  RooTruthModel deltaFunc("deltaFunc", "deltaFunc", dt);


   RooDecay Model("Model","Model", dt, tau, deltaFunc, RooDecay::DoubleSided);  

 


 
  
  //Fit
   RooFitResult* fitRes = Model.fitTo(*data,Save());


  TCanvas* can = new TCanvas("c","c", 900, 1200) ;
  
    //dt PLOTING
  RooPlot *xframe_1 =dt.frame(Bins(50),Title("B^{0} lifetime fit"));
  data->plotOn(xframe_1);
  Model.plotOn(xframe_1 ,LineColor(kBlue));
  RooHist* pull_dt = xframe_1->pullHist() ;
  Model.paramOn(xframe_1);

  RooPlot* f_pull_dt = dt.frame(Title("")) ;
  f_pull_dt->addPlotable(pull_dt,"BX0") ;
  f_pull_dt->SetMaximum(6);
  f_pull_dt->SetMinimum(-6);
  f_pull_dt->GetYaxis()->SetLabelSize(0.11);
  f_pull_dt->GetXaxis()->SetLabelSize(0.11);
 
 TPad *pad1 = new TPad("pad1", "pad 1",0.0,0.2,1.0,1.0,0);
 TPad *pad2 = new TPad("pad2", "pad 2",0.0,0.0,1.0,0.2,0);
 pad1->Draw();   pad2->Draw();
 pad1->cd() ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ; pad1->SetLogy();
 pad2->cd() ; gPad->SetLeftMargin(0.15) ; gPad->SetBottomMargin(0.15); f_pull_dt->Draw() ; 

/*
  can->cd() ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
 can->SetLogy(); 
*/
can->SaveAs(Imagename);
RooArgSet(fitRes->floatParsFinal()).writeToFile(Resultname);

return 0;
 


}















