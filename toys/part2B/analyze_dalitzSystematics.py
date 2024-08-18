from ROOT import RooRealVar, RooGaussian, RooFit, RooDataSet, RooArgSet, kRed, kSolid, TCanvas, gROOT, TGraphErrors

import numpy as np
import pandas as pd
from array import array


gROOT.SetBatch(True)



results_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/systematic_results'
out_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/systematic_plots'




#components = ['bb', 'cont', 'scf', 'rare', 'signal', 'all']
#components = ['bb', 'cont', 'rare', 'signal']
#components = ['bbrand', 'cont', 'scf', 'bbmissfsp', 'signal', 'all']
#components = ['dummy']
components = ['all']


def fit_distribution(component, results_dir, out_dir, parname):

	
	resfile = '%s/%s_%s.txt'%(results_dir, parname, component)  

	array = pd.read_csv(resfile, header = None, delimiter = '\t')
	array = array.to_numpy()
	val = array[:, 1]; val_range = np.percentile(val, 95.0) - np.percentile(val, 5.0)
	val_min = np.percentile(val, 5.0) - (2.5 * val_range); val_max = np.percentile(val, 95.0) + (2.5 * val_range)
	e_val = array[:, 2]; err_range = np.percentile(e_val, 95.0) - np.percentile(e_val, 5.0)
	val_distMean = np.mean(val)
	val_sd = np.std(val); val_sd = round(val_sd, 5)

	#'par' as in parameter
	value = RooRealVar('value', '', val_min, val_max, '')
	data_val = RooDataSet('data_val','data_val',RooArgSet(value)) #value dist.

	for i in range(val.shape[0]):
		value.setVal(val[i]); data_val.add(RooArgSet(value))


	#DEFINE PDF
	mean = RooRealVar("mean","#mu", val_distMean, val_distMean - (val_range/3), val_distMean + (val_range/3));
	sig = RooRealVar("sig","#sigma", val_range, 0.0, 2.5 * val_range);
	Model = RooGaussian("model", "Gaussian pull", value, mean, sig);

	fitRes = Model.fitTo(data_val)

	#Plot pull
	data_frame = value.frame(RooFit.Bins(20), RooFit.Title('%s estimated distribution'%parname))
	data_val.plotOn(data_frame)
	Model.plotOn(data_frame, RooFit.LineColor(kRed), RooFit.LineStyle(kSolid), RooFit.LineWidth(1));
	Model.paramOn(data_frame)
	data_frame.GetXaxis().SetTitle('Distribution std. dev. = %s'%val_sd)

	can = TCanvas("can","",1200, 600); 
	can.cd(); data_frame.Draw()
	can.SaveAs('%s/%s_%s.png'%(out_dir, parname, component))


	return mean.getVal(), mean.getError(), sig.getVal(), sig.getError()
	

if __name__ == '__main__':

	pars = ['S_plus', 'S_minus', 'A']
	#pars = ['S_I', 'S_Ib', 'A']	

	#Pull distributions 
	for component  in components:
		for par in pars:
			fit_distribution(component, results_dir, out_dir, par)








	


