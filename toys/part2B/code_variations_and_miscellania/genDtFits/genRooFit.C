
using namespace RooFit;

  RooRealVar dt("dt","dt",-10.0,10.00,"ps");
  RooDataSet* data_p=new RooDataSet("data_p","data_p",RooArgSet(dt));
  RooDataSet* data_n=new RooDataSet("data_n","data_n",RooArgSet(dt));

int genRooFit(){


gStyle->SetOptStat(0);




//TFile * f = new TFile("merged_SMC_S0p2_100k_wRMVA_wB2MVA_sel_BCS_branches.root");
//TFile * f = new TFile("merged_SMC_S0p4_100k_wRMVA_wB2MVA_sel_BCS_branches.root");
//TFile * f = new TFile("merged_SMC_S0p6_100k_wRMVA_wB2MVA_sel_BCS_branches.root");

//TFile * f = new TFile("merged_SMC_Sm0p2_100k_wRMVA_wB2MVA_sel_BCS_branches.root");
//TFile * f = new TFile("merged_SMC_Sm0p4_100k_wRMVA_wB2MVA_sel_BCS_branches.root");
TFile * f = new TFile("merged_SMC_Sm0p6_100k_wRMVA_wB2MVA_sel_BCS_branches.root");




TTree * tree = (TTree *) f->Get("h1");

float nEntries = 0;

 
float f_dtgen, f_genbq;

tree->SetBranchAddress("dtgen",&f_dtgen);
tree->SetBranchAddress("genbq",&f_genbq);





for(int i = 0 ; i < tree->GetEntries(); i++ )
{


tree->GetEntry(i);
dt.setVal(f_dtgen);


if(fabs(f_dtgen) < 10 ){
 nEntries++;
if(f_genbq == 1) data_p->add(RooArgSet(dt));
else if(f_genbq == -1) data_n->add(RooArgSet(dt));

}//fabs(f_dtgen) condition



} 

 

//PDF 
  RooRealVar N("N","N",nEntries); 
  RooRealVar dm("#Delta m","dm",0.507);//
  RooRealVar A("A","A", 0.0, -1.0, 1.0);//
  RooRealVar S("S","S", 0.0, -1.0, 1.0);//
  RooRealVar tau("tau","tau", 1.525);//, 1.3, 1.7);//

 //B0
 RooGenericPdf expo_p("expo_p","(@5/(4*@1))*exp(-abs(@0)/@1)*(1+ (@2*sin(@4*@0)) + (@3*cos(@4*@0)) )",RooArgList(dt, tau, S, A, dm, N));
 //B0-bar
  RooGenericPdf expo_n("expo_p","(@5/(4*@1))*exp(-abs(@0)/@1)*(1 - (@2*sin(@4*@0)) - (@3*cos(@4*@0)) )",RooArgList(dt, tau, S, A, dm, N));
 
 //Asymmetry function
 //RooGenericPdf osc("osc","(@1*sin(@3*@0)) + (@2*cos(@3*@0)) ", RooArgList(dt, S, A, dm));
 
   //CREATE INDEX CATEGORY AND JOIN  SAMPLEs
  //_____________________________________________________
  RooCategory sample("sample","sample");
  sample.defineType("B0",+1);
  sample.defineType("B0b",-1);

  
  
  RooDataSet* combData =new RooDataSet("combData","combData",RooArgSet(dt),Index(sample),Import("B0",*data_p),Import("B0b",*data_n));

 


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(expo_p,"B0");
  simPdf.addPdf(expo_n,"B0b");
  //_____________________________________________________
  
  
  //Fit
   RooFitResult* fitRes = simPdf.fitTo(*combData,Save());


const float afix = 0.5;
TF1 *f1 = new TF1("myfunc","([0]*sin(0.507*x)) + ([1]*cos(0.507*x))",-10, 10);

RooArgSet pars(S,A);
RooAbsReal *osc = bindFunction(f1, dt, pars);


  TCanvas* can = new TCanvas("c","c", 900, 1200) ;
  can->Divide(1,2) ;
  
    //dt PLOTING
  RooPlot *xframe_1 =dt.frame(Bins(50),Title("B^{0}"));
  combData->plotOn(xframe_1,Cut("sample==sample::B0"));
    combData->plotOn(xframe_1,Cut("sample==sample::B0b"));
  simPdf.plotOn(xframe_1,Slice(sample,"B0"),ProjWData(sample,*combData),LineColor(kRed));
    simPdf.plotOn(xframe_1,Slice(sample,"B0b"),ProjWData(sample,*combData),LineColor(kBlue));
  expo_p.paramOn(xframe_1,data_p);

 
 
   RooPlot *xframe_3 = dt.frame(Bins(50),Title("Asymmetry"));
   combData->plotOn(xframe_3, Asymmetry(sample));
   osc->plotOn(xframe_3);
   //simPdf.plotOn(xframe_3, Asymmetry(sample));
   //osc.plotOn(xframe_3);
 

  can->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe_1->GetYaxis()->SetTitleOffset(1.4) ; xframe_1->Draw() ;
  can->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe_3->GetYaxis()->SetTitleOffset(1.4) ; xframe_3->Draw() ;
  

return 0;



}















