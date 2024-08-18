
using namespace RooFit;

  RooRealVar dt("dt","dt",-10.0,10.00,"ps");
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dt));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dt));

int genRooBCPFit(string filename){

string Sample = filename;

Sample.erase(filename.length()-5);
string imagename = "plots/" + Sample + ".png";
string resultname = "results/" +  Sample + ".txt";
const char* Imagename = imagename.c_str();
const char* Resultname = resultname.c_str();
const char* Filename = filename.c_str();

gStyle->SetOptStat(0);



TFile * f = new TFile(Filename);


TTree * tree = (TTree *) f->Get("h1");

float nEntries = 0;

 
float f_dtgen, f_genbq, f_issgnevt;

tree->SetBranchAddress("dtgen",&f_dtgen);
tree->SetBranchAddress("genbq",&f_genbq);
tree->SetBranchAddress("issgnevt",&f_issgnevt);




for(int i = 0 ; i < tree->GetEntries(); i++ )
{


tree->GetEntry(i);
dt.setVal(f_dtgen);



if(f_issgnevt == 1){
if(fabs(f_dtgen) < 10 ){
 nEntries++;
if(f_genbq == 1) data_p->add(RooArgSet(dt));
else if(f_genbq == -1) data_n->add(RooArgSet(dt));
}


}//fabs(f_dtgen) condition



} 

 

//PDF 
 
  RooRealVar dm("#Delta m","dm", 0.507);//
  RooRealVar A("A","A", 0.0, -1.0, 1.0);//
  RooFormulaVar C("C", "-1.0*A", RooArgList(A));
  RooRealVar S("S","S", 0.0, -1.0, 1.0);//
  RooRealVar tau("tau","tau", 1.525);//
  RooRealVar avgMistag("avgMistag","avgMistag", 0.0);
  RooRealVar delMistag("delMistag","delMistag", 0.0);
  RooRealVar mu("mu","mu", 0.0);//, -1.0, 1.0);

 
 
   //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("B0",+1);
  sample.defineType("B0b",-1);
  
  
  RooTruthModel deltaFunc("deltaFunc", "deltaFunc", dt);
  RooBCPGenDecay Model("Model", "Model", dt, sample, tau, dm, avgMistag, C, S, delMistag, mu, deltaFunc);
  
  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dt),Index(sample),Import("B0",*data_p),Import("B0b",*data_n));

 


 
  
  //Fit
   RooFitResult* fitRes = Model.fitTo(*combData,Save());


const float afix = 0.5;
TF1 *f1 = new TF1("myfunc","[2]+(([0]*sin(0.507*x)) + ([1]*cos(0.507*x)))",-10, 10);

RooArgSet pars(S,A,mu);
RooAbsReal *osc = bindFunction(f1, dt, pars);


  TCanvas* can = new TCanvas("c","c", 900, 1200) ;
  can->Divide(1,2) ;
  
    //dt PLOTING
  RooPlot *xframe_1 =dt.frame(Bins(50),Title("B^{0}"));
  combData->plotOn(xframe_1,Cut("sample==sample::B0"));
    combData->plotOn(xframe_1,Cut("sample==sample::B0b"));
  Model.plotOn(xframe_1,Slice(sample,"B0"),ProjWData(sample,*combData),LineColor(kRed));
    Model.plotOn(xframe_1,Slice(sample,"B0b"),ProjWData(sample,*combData),LineColor(kBlue));
  Model.paramOn(xframe_1);

 
 
   RooPlot *xframe_3 = dt.frame(Bins(50),Title("Asymmetry"));
   combData->plotOn(xframe_3, Asymmetry(sample));
   osc->plotOn(xframe_3);
 
 
 

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;
  
can->SaveAs(Imagename);
RooArgSet(fitRes->floatParsFinal()).writeToFile(Resultname);

return 0;
 


}















