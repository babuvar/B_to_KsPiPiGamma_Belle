#Seems to behave a bit suspiciously. Go for the C++ version instead

from ROOT import RooRealVar, RooDataSet, gStyle, TFile, TTree, RooGenericPdf, RooCategory, RooSimultaneous, RooFitResult, RooPlot, RooArgSet, gROOT, RooArgList, RooFit, TF1, TCanvas, kRed, kBlue, gPad
import numpy as np


dt = RooRealVar("dt","dt",-10.0,10.00,"ps");
data_p = RooDataSet("data_p","data_p",RooArgSet(dt));
data_n = RooDataSet("data_n","data_n",RooArgSet(dt));

gStyle.SetOptStat(0)
gROOT.SetBatch(True)


rootdir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/'

fnames = ['merged_SMC_S0p2_100k_wRMVA_wB2MVA_sel_BCS_branches.root', 'merged_SMC_S0p4_100k_wRMVA_wB2MVA_sel_BCS_branches.root', 'merged_SMC_S0p6_100k_wRMVA_wB2MVA_sel_BCS_branches.root', 'merged_SMC_Sm0p2_100k_wRMVA_wB2MVA_sel_BCS_branches.root', 'merged_SMC_Sm0p4_100k_wRMVA_wB2MVA_sel_BCS_branches.root', 'merged_SMC_Sm0p6_100k_wRMVA_wB2MVA_sel_BCS_branches.root']




def read_and_genfit(fname, roodir):

	rootfile = '%s/%s'%(roodir, fname)

	f = TFile(rootfile)
	tree = f.Get('h1')

	#number of entried that goes into the datasets
	nEntries = 0


	for i in range(tree.GetEntries()):

		tree.GetEntry(i)
		dt.setVal(tree.dtgen)
		if np.fabs(tree.dtgen) < 10:
			nEntries = nEntries + 1
			if tree.genbq == 1:
				data_p.add(RooArgSet(dt))
			elif tree.genbq == -1:
				data_n.add(RooArgSet(dt)) 
			 
 

	#PDF 
	N = RooRealVar('N','N',nEntries)
	dm = RooRealVar('#Delta m','dm',0.507) # Evtgen value
	A = RooRealVar('A','A', 0.0, -1.0, 1.0)
	S = RooRealVar('S','S', 0.0, -1.0, 1.0)
	tau = RooRealVar('tau','tau', 1.525) # Evtgen value

	#B0
	expo_p = RooGenericPdf('expo_p','(@5/(4*@1))*exp(-abs(@0)/@1)*(1+ (@2*sin(@4*@0)) + (@3*cos(@4*@0)) )',RooArgList(dt, tau, S, A, dm, N))
	#B0-bar
	expo_n = RooGenericPdf('expo_n','(@5/(4*@1))*exp(-abs(@0)/@1)*(1 - (@2*sin(@4*@0)) - (@3*cos(@4*@0)) )',RooArgList(dt, tau, S, A, dm, N))
 
 
	#CREATE INDEX CATEGORY AND JOIN  SAMPLEs
	sample = RooCategory('sample','sample')
	sample.defineType('B0',+1)
	sample.defineType('B0b',-1)

  
  	#combined data
	combData = RooDataSet('combData','combData', RooArgSet(dt), RooFit.Index(sample), RooFit.Import('B0',data_p), RooFit.Import('B0b',data_n))

	#combined pdf
	simPdf = RooSimultaneous('simPdf','simPdf',sample)
	simPdf.addPdf(expo_p,'B0')
	simPdf.addPdf(expo_n,'B0b')
  
  
	#Fit
	fitRes = simPdf.fitTo(combData, RooFit.Save())

	#TF1 definition of asymmetry function
	f1 = TF1('myfunc','([0]*sin(0.507*x)) + ([1]*cos(0.507*x))',-10, 10)
	pars = RooArgSet(S,A)
	osc = RooFit.bindFunction(f1, dt, pars)


	#Plotting
	can = TCanvas('c','c', 900, 1200)
	can.Divide(1,2) 
  
	#dt PLOTING
	xframe_1 =dt.frame(RooFit.Bins(50), RooFit.Title('B^{0}'))
	combData.plotOn(xframe_1, RooFit.Cut('sample==sample::B0'))
	combData.plotOn(xframe_1, RooFit.Cut('sample==sample::B0b'))
	simPdf.plotOn(xframe_1, RooFit.Slice(sample,'B0'), RooFit.ProjWData(RooArgSet(sample), combData), RooFit.LineColor(kRed))
	simPdf.plotOn(xframe_1, RooFit.Slice(sample,'B0b'), RooFit.ProjWData(RooArgSet(sample), combData), RooFit.LineColor(kBlue))
	expo_p.paramOn(xframe_1)

 	#Asymmetry
	xframe_3 = dt.frame(RooFit.Bins(50), RooFit.Title('Asymmetry'))
	combData.plotOn(xframe_3, RooFit.Asymmetry(sample))
	osc.plotOn(xframe_3)
 

	can.cd(1) ; gPad.SetLeftMargin(0.15) ; xframe_1.GetYaxis().SetTitleOffset(1.4) ; xframe_1.Draw() 
	can.cd(2) ; gPad.SetLeftMargin(0.15) ; xframe_3.GetYaxis().SetTitleOffset(1.4) ; xframe_3.Draw() 
	
	#Save canvas
	can.SaveAs('%s.png'%fname[:len(fname) - 5])
	del can
	 

if __name__ == '__main__': 
	
	for fname in fnames:
		read_and_genfit(fname, rootdir)












