from ROOT import RooRealVar, RooGaussian, RooFit, RooDataSet, RooArgSet, kRed, kSolid, TCanvas, gROOT, TGraphErrors

import numpy as np
import pandas as pd
from array import array


gROOT.SetBatch(True)

results_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/consolidated_results'; out_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/ensemble_plots'



def plot_pulls(rfname, results_dir, Trueval, parname,  out_dir):

	
	resfile = '%s/%s.txt'%(results_dir, rfname)  

	array = pd.read_csv(resfile, header = None, delimiter = '\t')
	array = array.to_numpy()
	val = array[:, 1]; val_range = np.percentile(val, 95.0) - np.percentile(val, 5.0)
	val_min = np.percentile(val, 5.0) - (0.5 * val_range); val_max = np.percentile(val, 95.0) + (0.5 * val_range)
	e_val = array[:, 2]; err_range = np.percentile(e_val, 95.0) - np.percentile(e_val, 5.0)
	err_min = np.percentile(e_val, 5.0) - (0.5 * err_range); err_max = np.percentile(e_val, 95.0) + (0.5 * err_range)


	#'par' as in parameter
	#par = RooRealVar('par','', -6.0, 6.0,'') 
	par = RooRealVar('par','', -12.0, 12.0,'')
	value = RooRealVar('value', '', val_min, val_max, '')
	err = RooRealVar('err', '', err_min, err_max, '')
	data_par = RooDataSet('data_par','data_par',RooArgSet(par)) #pull dist.
	data_val = RooDataSet('data_val','data_val',RooArgSet(value)) #value dist.
	data_err = RooDataSet('data_err','data_err',RooArgSet(err)) #error dist.	

	for i in range(val.shape[0]):
		
		pull = (val[i] - Trueval)/e_val[i]
		par.setVal(pull); data_par.add(RooArgSet(par))
		value.setVal(val[i]); data_val.add(RooArgSet(value))
		err.setVal(e_val[i]); data_err.add(RooArgSet(err))


	#DEFINE PDF
	mean = RooRealVar("mean","#mu",0.0,-10.,10.0);
	sig = RooRealVar("sig","#sigma",1.0,0.3,2.2);
	Model = RooGaussian("model", "Gaussian pull", par, mean, sig);

	fitRes = Model.fitTo(data_par)

	#Plot pull
	data_frame = par.frame(RooFit.Bins(20), RooFit.Title('%s Fit-Pull'%parname))
	data_par.plotOn(data_frame)
	Model.plotOn(data_frame, RooFit.LineColor(kRed), RooFit.LineStyle(kSolid), RooFit.LineWidth(1));
	Model.paramOn(data_frame)

	#Mean residuals and mean errors
	fitVal_median = np.median(val)
	fitVal_median_residual = fitVal_median - Trueval
	fitVal_mad = 1.4826 * np.median(abs(val - fitVal_median))
	fitVal_median_error = fitVal_mad / np.sqrt(len(val))
	fitVal_mad_error = np.sqrt(2) * fitVal_mad / np.sqrt(len(val) - 1)


	#Plot value
	val_frame = value.frame(RooFit.Bins(20), RooFit.Title('%s Fit-estimate'%parname))
	data_val.plotOn(val_frame)

	#Plot_error
	err_frame = err.frame(RooFit.Bins(20), RooFit.Title('%s Fit-Error'%parname))
	data_err.plotOn(err_frame)

	can = TCanvas("can","",1200, 600); can.Divide(3,1)
	can.cd(1); data_frame.Draw()
	can.cd(2); val_frame.Draw()
	can.cd(3); err_frame.Draw()
	can.SaveAs('%s/%s.png'%(out_dir, rfname))

	print("median = ", fitVal_median)
	print("width = ", fitVal_mad)

	return 
	



if __name__ == '__main__':


	#Pull distributions for tauB
	Trueval = 1.525; rfname = 'tauB_S0p0'
	plot_pulls(rfname, results_dir, Trueval, 'tau_B', out_dir)



	


