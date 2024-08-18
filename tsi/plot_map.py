from ROOT import RooRealVar, RooGaussian, RooFit, RooDataSet, RooArgSet, kRed, kSolid, TCanvas, gROOT, TGraphErrors, TH2D, gStyle, gPad

import numpy as np
import pandas as pd
from array import array


gROOT.SetBatch(True)
gStyle.SetOptStat(0)

resfile = 'output.txt'
array = pd.read_csv(resfile, header = None, delimiter = '\t')
array = array.to_numpy()
S = array[:, 0];
S_syst = array[:, 1];	
A =  array[:, 2];
hist= TH2D('tsi', 'TSI : #delta S; Observed S; Observed A', 13, -0.65, 0.65, 7, -0.35, 0.35)
for i in range(S.shape[0]):
		hist.Fill(S[i], A[i], S_syst[i])

can = TCanvas('can', 'can', 800, 600)
can.cd()
hist.Draw('colz'); gPad.SetRightMargin(0.15) ;
can.SaveAs('tsi_map.png')


























