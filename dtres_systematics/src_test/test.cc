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

/*
cout<<"First test compilation"<<endl;
print_dtres_params(29, 1);

for(int i=0; i<=40; i++){
dtres_systematics_b0(i, 1.5, 29, 1);
}

print_dtres_params(29, 1);

for(int i=0; i<=40; i++){
dtres_systematics_b0(i, 1.0, 29, 1);
}


print_dtres_params(29, 1);
*/

string parName;
for(int i = 2; i<=38; i++){
double parVal = get_dtres_parVal_b0(i, 0, 31, 1, parName);
cout<<"dtres_parVal["<<i<<"] ("<<parName<<") = "<<parVal<<endl;
}

        return 0;
}

















