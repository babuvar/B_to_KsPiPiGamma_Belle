################################################################################################################
##                    OPTIMIZE the ROOT-TMVA CS classifier variable                                           ##
################################################################################################################

from ROOT import TChain, TGraph, TCanvas
from ROOT import kBlack, kBlue, kRed
from array import array
from math import sqrt
from ROOT import gROOT
gROOT.SetBatch(True)
import os
from operator import add

################################################################################################################
##                                        Function definitions                                                ##
################################################################################################################


# outputs nsig[i], nbkg[i], nsig_tot and nbkg_tot, which will be used to obtain the various metrics of optimization
def get_stats(rootdir, fileNames, key, variables_used, opt_variable, sig_variable, selections, cut_nbins, cut_start, cut_width):

	print('starting stats count')

	#quantities of interest that will be output
	nsig = [0.0] * cut_nbins; nbkg = [0.0] * cut_nbins; nsig_tot = 0.0; nbkg_tot =0.0; ntot = [0.0] * cut_nbins # ntot will be used for data/mc bkg correction by looking at sideband

	fileList = [rootdir + '/' + s  for s in fileNames]
	

	chain = TChain(key)
	for rootFile in fileList:
		chain.Add(rootFile)

	branches = {}
	for branchName in variables_used:
		branches[branchName] = array('f', [-999])
		chain.SetBranchAddress(branchName, branches[branchName])

	for i in range(chain.GetEntries()):
		chain.GetEntry(i);
		selection_conditions = '' # selections before evaluating
		signal_condition = 'chain.%s == 1'%sig_variable # to determine whether signal or background
		
		for idx, selection in enumerate(selections):
			if idx < (len(selections) - 1):
				selection_conditions = selection_conditions + 'chain.%s and '%selection
			else:
				selection_conditions = selection_conditions + 'chain.%s'%selection

		if eval(selection_conditions):
			# count nsig_tot and nbkg_tot (dummy vals for data)
			if eval(signal_condition):
				nsig_tot = nsig_tot + 1
			else:
				nbkg_tot = nbkg_tot + 1
			# count cut-dependent nsig and nbkg
			for j in range(cut_nbins): 
				optCut_val = cut_start + (j*cut_width)
				optCut_condition = 'chain.%s > %s'%(opt_variable, optCut_val)
				if eval(optCut_condition):
					ntot[j] = ntot[j] + 1 # will be only used when looking at sideband-region
					#dummy values of nsig[j] and nbkg[j] for data
					if eval(signal_condition):
						nsig[j] = nsig[j] +1
					else:
						nbkg[j] = nbkg[j] +1

	print('Finished stats count')	
	return nsig_tot, nbkg_tot, nsig, nbkg, ntot

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# outputs metrics of interest, i.e., f.o.m., signal-efficiency, etc.
def get_metrics(component_dic, key, variables_used, opt_variable, sig_variable, selections, cut_nbins, cut_start, cut_width, doDataMC_bkg_correction):

	#initialize stats
	nsig_tot = 0.0; nbkg_tot = 0.0; 
	nsig = [0.0] * cut_nbins; nbkg = [0.0] * cut_nbins; ntot = [0.0] * cut_nbins


	#get basic signal/background stats from mc components within signal-window
	for component_name, component in component_dic['component_MC'].items():

		rootdir = component['dir']
		fileNames = component['files']
		weight = component['weight']

		nsig_tot_comp, nbkg_tot_comp, nsig_comp, nbkg_comp, ntot_comp = get_stats(rootdir, fileNames, key, variables_used, opt_variable, sig_variable, selections['signalwindow'], cut_nbins, cut_start, cut_width)		

		nsig_tot = nsig_tot + (nsig_tot_comp * weight)
		nbkg_tot = nbkg_tot + (nbkg_tot_comp * weight)
		nsig_comp = [i * weight for i in nsig_comp]; nsig = list( map(add, nsig, nsig_comp) )
		nbkg_comp = [i * weight for i in nbkg_comp]; nbkg = list( map(add, nbkg, nbkg_comp) )	
		ntot_comp = [i * weight for i in ntot_comp]; ntot = list( map(add, ntot, ntot_comp) )

	
	#look at sideband-region
	if doDataMC_bkg_correction == True :
		ntot_mcsb = [0.0] * cut_nbins; ntot_datasb = [0.0] * cut_nbins

		#look at mc-sideband
		for component_name, component in component_dic['component_MC'].items():
			rootdir = component['dir']
			fileNames = component['files']
			weight = component['weight']

			dummy1, dummy2, dummy3, dummy4, ntot_mcsb_comp = get_stats(rootdir, fileNames, key, variables_used, opt_variable, sig_variable, selections['sideband'], cut_nbins, cut_start, cut_width)

			ntot_mcsb_comp = [i * weight for i in ntot_mcsb_comp]; ntot_mcsb = list( map(add, ntot_mcsb, ntot_mcsb_comp) )

		#look at data-sideband
		for component_name, component in component_dic['component_Data'].items():
			rootdir = component['dir']
			fileNames = component['files']
			weight = component['weight']

			dummy1, dummy2, dummy3, dummy4, ntot_datasb_comp = get_stats(rootdir, fileNames, key, variables_used, opt_variable, sig_variable, selections['sideband'], cut_nbins, cut_start, cut_width)

			ntot_datasb_comp = [i * weight for i in ntot_datasb_comp]; ntot_datasb = list( map(add, ntot_datasb, ntot_datasb_comp) )


	
	#Evaluating metrics
	print('Starting to evaluate optimization metrics')
	fom = [0.0] * len(nsig); sigeff = [0.0] * len(nsig); cut = [0.0] * len(nsig); data_mc_bkg_ratio = [1.0] * len(nsig)
	
	#scaled nbkg : This scaling refers to the #data_mc_bkg_ratio, if applicable
	nbkg_scaled = [0] * len(nsig)
	fom_best = -999; i_best = 0
	
	for i in range(len(nsig)):
		#array of cuts
		cut[i] = cut_start + (i*cut_width)
		
		#calculate signal efficiency
		sigeff[i] = (nsig[i] * 100) / nsig_tot #in percentage
		
		#data_mc_bkg_ratio
		if doDataMC_bkg_correction == True :
			data_mc_bkg_ratio[i] = ntot_datasb[i] / ntot_mcsb[i]
		
		#scaling nbkg cut-wise
		nbkg_scaled[i] = nbkg[i] * data_mc_bkg_ratio[i]
		
		#calculate fom
		fom[i] = nsig[i] / sqrt( nsig[i] + nbkg_scaled[i])
		if fom[i] > fom_best:
			fom_best = fom[i]; i_best = i
		
	opt_cut = cut[i_best]
		
	print('Done evaluating optimization metrics')	
	return cut, sigeff, fom, opt_cut, data_mc_bkg_ratio

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	

def plot_metric(cut, metric, labels, opt_cut, opt_variable):

	print('Starting to plot...')
	cnv = TCanvas('cnv','cnv', 1200, 800) ; 

	numbins = len(metric); cut_array = array('f', cut); metric_array = array('f',  metric)
	
	labels['title'] = '%s'%labels['title']
	statement = 'optimal cut is %s > %.2f'%(opt_variable, opt_cut)
	labels['x'] = '%s (%s)'%(labels['x'], statement)
	
	gr = TGraph(numbins, cut_array, metric_array)
	gr.SetTitle(labels['title']); gr.GetXaxis().SetTitle(labels['x']); gr.GetYaxis().SetTitle(labels['y']);
	gr.SetMarkerStyle(20); gr.SetMarkerSize(0.7); gr.SetMarkerColor(kBlue);
	gr.GetXaxis().SetTitleSize(0.06); gr.GetXaxis().SetTitleOffset(0.8); gr.GetXaxis().CenterTitle();
	gr.GetYaxis().SetTitleSize(0.06); gr.GetYaxis().SetTitleOffset(0.8); gr.GetYaxis().CenterTitle();


	cnv.cd(); gr.Draw('APL')
	
	#make plots directory if it does not exist
	if not os.path.exists('plots_%s'%opt_variable):
		os.makedirs('plots_%s'%opt_variable)
	
	cnv.SaveAs('plots_%s/%s'%(opt_variable, labels['name']))
	print('Finished plotting')
	return

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def plot_metrics(cut, sigeff, fom, data_mc_bkg_ratio, labels, opt_cut, opt_variable, doDataMC_bkg_correction):

	print('Starting to plot...')
	
	plot_metric(cut, sigeff, labels['sigEff'], opt_cut, opt_variable)
	plot_metric(cut, fom, labels['fom'], opt_cut, opt_variable)
	if doDataMC_bkg_correction == True :
		plot_metric(cut, data_mc_bkg_ratio, labels['data_mc_bkg_ratio'], opt_cut, opt_variable)
	
	print('Finished plotting')
	return
	
	



