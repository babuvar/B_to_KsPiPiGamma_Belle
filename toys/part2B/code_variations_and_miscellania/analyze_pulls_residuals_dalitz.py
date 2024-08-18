from ROOT import RooRealVar, RooGaussian, RooFit, RooDataSet, RooArgSet, kRed, kSolid, TCanvas, gROOT, TGraphErrors

import numpy as np
import pandas as pd
from array import array


gROOT.SetBatch(True)

results_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/consolidated_results'; out_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/ensemble_plots'



rfnames_Sm = ['S_minus_Sm0p6', 'S_minus_Sm0p4', 'S_minus_Sm0p2', 'S_minus_S0p0', 'S_minus_S0p2', 'S_minus_S0p4', 'S_minus_S0p6']
trueval_S = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
trueval_Sp = [-1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2]
trueval_Sm = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
rfnames_Sp = ['S_plus_Sm0p6', 'S_plus_Sm0p4', 'S_plus_Sm0p2', 'S_plus_S0p0', 'S_plus_S0p2', 'S_plus_S0p4', 'S_plus_S0p6']
trueval_A = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
rfnames_A = ['A_Sm0p6', 'A_Sm0p4', 'A_Sm0p2', 'A_S0p0', 'A_S0p2', 'A_S0p4', 'A_S0p6']

def get_fitValMedianWidth(rfname, results_dir, Trueval, isfixedtrueval, parname, input_S, out_dir):

	resfile = '%s/%s.txt'%(results_dir, rfname)

	array = pd.read_csv(resfile, header = None, delimiter = '\t')
	array = array.to_numpy()
	val = array[:, 1];
	fitVal_median = np.median(val)
	fitVal_median_residual = fitVal_median - Trueval
	fitVal_mad = 1.4826 * np.median(abs(val - fitVal_median))
	fitVal_median_error = fitVal_mad / np.sqrt(len(val))
	fitVal_mad_error = np.sqrt(2) * fitVal_mad / np.sqrt(len(val) - 1)
	
	return fitVal_median_residual, fitVal_median_error, fitVal_mad, fitVal_mad_error


def plot_pulls(rfname, results_dir, Trueval, isfixedtrueval, parname, input_S, out_dir):

	
	resfile = '%s/%s.txt'%(results_dir, rfname)  

	array = pd.read_csv(resfile, header = None, delimiter = '\t')
	array = array.to_numpy()
	val = array[:, 1]; val_range = np.percentile(val, 95.0) - np.percentile(val, 5.0)
	val_min = np.percentile(val, 5.0) - (0.5 * val_range); val_max = np.percentile(val, 95.0) + (0.5 * val_range)
	e_val = array[:, 2]; err_range = np.percentile(e_val, 95.0) - np.percentile(e_val, 5.0)
	err_min = np.percentile(e_val, 5.0) - (0.5 * err_range); err_max = np.percentile(e_val, 95.0) + (0.5 * err_range)
	if isfixedtrueval == False:
		trueval = array[:, 3]


	#'par' as in parameter
	#par = RooRealVar('par','', -6.0, 6.0,'') 
	par = RooRealVar('par','', -12.0, 12.0,'')
	value = RooRealVar('value', '', val_min, val_max, '')
	err = RooRealVar('err', '', err_min, err_max, '')
	data_par = RooDataSet('data_par','data_par',RooArgSet(par)) #pull dist.
	data_val = RooDataSet('data_val','data_val',RooArgSet(value)) #value dist.
	data_err = RooDataSet('data_err','data_err',RooArgSet(err)) #error dist.	

	for i in range(val.shape[0]):
		if isfixedtrueval == False:
			Trueval = trueval[i]
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

	#Plot value
	val_frame = value.frame(RooFit.Bins(20), RooFit.Title('%s Fit-estimate. Input S-value = %s'%(parname, input_S)))
	data_val.plotOn(val_frame)

	#Plot_error
	err_frame = err.frame(RooFit.Bins(20), RooFit.Title('%s Fit-Error'%parname))
	data_err.plotOn(err_frame)

	can = TCanvas("can","",1200, 600); can.Divide(3,1)
	can.cd(1); data_frame.Draw()
	can.cd(2); val_frame.Draw()
	can.cd(3); err_frame.Draw()
	can.SaveAs('%s/%s.png'%(out_dir, rfname))


	return mean.getVal(), mean.getError(), sig.getVal(), sig.getError()

	

def plot_linearity(n, metric, e_metric, inputS, parname, metricname, out_dir):


	e_inputS = [0.0]*len(metric)
	x = array( 'f', inputS)
	ex = array( 'f', e_inputS)
	y = array( 'f', metric)
	ey = array( 'f', e_metric)
	gr = TGraphErrors(n, x, y, ex,ey) 
	gr.SetTitle('%s %s vs input-S; Input S; Residual'%(parname, metricname)); gr.SetMarkerColor(4); gr.SetMarkerStyle(21)
	can = TCanvas("can","",800, 500); 
	can.cd(); gr.Draw("ALP")	
	can.SaveAs('%s/%s_pull_%s_linearity.png'%(out_dir, parname, metricname))
	can.SaveAs('%s/%s_pull_%s_linearity.root'%(out_dir, parname, metricname))

	return



if __name__ == '__main__':



	mean_res_Sp  = []; e_mean_Sp = []; width_Sp = []; e_width_Sp = []
	mean_res_Sm  = []; e_mean_Sm = []; width_Sm = []; e_width_Sm = []
	mean_res_A  = []; e_mean_A = []; width_A = []; e_width_A = []
	
	#Pull distributions for S_plus
	for i, (rfname, Trueval, S_val)  in enumerate(zip(rfnames_Sp, trueval_Sp, trueval_S)):
		mean_res, e_mean, width, e_width = get_fitValMedianWidth(rfname, results_dir, Trueval, True, 'S', S_val, out_dir)
		mean_res_Sp.append(mean_res); e_mean_Sp.append(e_mean)
		width_Sp.append(width); e_width_Sp.append(e_width)	
		
	#Pull distributions for S_minus
	for i, (rfname, Trueval, S_val)  in enumerate(zip(rfnames_Sm, trueval_Sm, trueval_S)):
		mean_res, e_mean, width, e_width = get_fitValMedianWidth(rfname, results_dir, Trueval, True, 'S', S_val, out_dir)
		mean_res_Sm.append(mean_res); e_mean_Sm.append(e_mean)
		width_Sm.append(width); e_width_Sm.append(e_width)

	#Pull distributions for A
	for i, (rfname, Trueval, S_val)  in enumerate(zip(rfnames_A, trueval_A, trueval_S)):
		mean_res, e_mean, width, e_width = get_fitValMedianWidth(rfname, results_dir, Trueval, True, 'A', S_val, out_dir)
		mean_res_A.append(mean_res); e_mean_A.append(e_mean)
		width_A.append(width); e_width_A.append(e_width)
	
	
	plot_linearity(7, mean_res_Sp, e_mean_Sp, trueval_S, 'S_plus', 'residual', out_dir)
	plot_linearity(7, width_Sp, e_width_Sp, trueval_S, 'S_plus', 'error', out_dir)
	plot_linearity(7, mean_res_Sm, e_mean_Sm, trueval_S, 'S_minus', 'residual', out_dir)
	plot_linearity(7, width_Sm, e_width_Sm, trueval_S, 'S_minus', 'error', out_dir)
	plot_linearity(7, mean_res_A, e_mean_A, trueval_S, 'A', 'residual', out_dir)
	plot_linearity(7, width_A, e_width_A, trueval_S, 'A', 'error', out_dir)
	





	


