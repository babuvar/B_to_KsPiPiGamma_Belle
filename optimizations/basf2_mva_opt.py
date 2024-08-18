################################################################################################################
##                    OPTIMIZE the ROOT-TMVA CS classifier variable                                           ##
################################################################################################################

from opt_helper import get_metrics  
from opt_helper import plot_metrics

################################################################################################################
##                                        Configuration area                                                  ##
################################################################################################################
#Defining cut variation
cut_start = -0.49
cut_end = 0.31; 
cut_nbins = 40;
cut_width = (cut_end - cut_start) / cut_nbins


#parameters of optimization

#doDataMC_bkg_correction = True 
doDataMC_bkg_correction = False

#rectangular signal-window
selections_signalwindow = [
'mbc > 5.27',
'mbc < 5.29',
'de > -0.16',
'de < 0.06'
]

selections_sideband = [
'mbc > 5.20',
'mbc < 5.26',
'de > -0.16',
'de < 0.06'
]
selections_cleanup=[
'pi0veto < 0.2',
'etaveto < 0.25',
'g_e9e25 > 0.95',
'g_cos > -0.65',
'g_cos < 0.86',
'tagcl > 0',
'vtxcl > 0',
'm23 > 0.6',
'm23 < 0.9',
'x_p0_kid < 0.25',
'x_p1_kid < 0.25',
'x_p0_eid < 0.25',
'x_p1_eid < 0.25',
'x_m < 1.8',
#'( tagh < 50 || tagntrk == 1)',
'vtxh < 50'
]

selections = {
'signalwindow' : selections_signalwindow+selections_cleanup,
'sideband' : selections_sideband+selections_cleanup
}

#Allowing more that two components, each with a different scaling

rareSig_files = ['MCrare_wRMVA_wB2MVA_shuffled_49PartsIn50_sig.root']
rareSig_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/prepare_recFiles_for_tatFits/splitfiles_with_myb2XMLFile'

rareSig_dic = {
'files' : rareSig_files,
'dir' : rareSig_dir,
'weight' : (1.0/49.0)
}

rareBkg_files = ['MCrare_wRMVA_wB2MVA_shuffled_49PartsIn50_bkg.root']
rareBkg_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/prepare_recFiles_for_tatFits/splitfiles_with_myb2XMLFile'

rareBkg_dic = {
'files' : rareBkg_files,
'dir' : rareBkg_dir,
'weight' : (1.0/49.0)
}

genericBkg_files = [
'merged_GMC_st0_charged_wRMVA_wB2MVA.root', 'merged_GMC_st0_mixed_wRMVA_wB2MVA.root',
'merged_GMC_st0_charm_wRMVA_wB2MVA.root', 'merged_GMC_st0_uds_wRMVA_wB2MVA.root', 'merged_GMC1_charged_wRMVA_wB2MVA.root', 'merged_GMC1_mixed_wRMVA_wB2MVA.root', 'merged_GMC1_charm_wRMVA_wB2MVA.root', 'merged_GMC1_uds_wRMVA_wB2MVA.root', 'merged_GMC3_charged_wRMVA_wB2MVA.root', 'merged_GMC3_mixed_wRMVA_wB2MVA.root', 'merged_GMC3_charm_wRMVA_wB2MVA.root', 'merged_GMC3_uds_wRMVA_wB2MVA.root', 'merged_GMC4_charged_wRMVA_wB2MVA.root', 'merged_GMC4_mixed_wRMVA_wB2MVA.root', 'merged_GMC4_charm_wRMVA_wB2MVA.root', 'merged_GMC4_uds_wRMVA_wB2MVA.root', 'merged_GMC5_charged_wRMVA_wB2MVA.root', 'merged_GMC5_mixed_wRMVA_wB2MVA.root', 'merged_GMC5_charm_wRMVA_wB2MVA.root', 'merged_GMC5_uds_wRMVA_wB2MVA.root']
genericBkg_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/prepare_recFiles_for_tatFits/outFiles_with_myb2XMLFile'

genericBkg_dic = {
'files' : genericBkg_files,
'dir' : genericBkg_dir,
'weight' : (1.0/5.0)
}

component_MC = {
'rareSig' : rareSig_dic,
'rareBkg' : rareBkg_dic,
'genericBkg' : genericBkg_dic
}

data_files = ['merged_data_sideband_wRMVA_wB2MVA.root']
data_dir = '/home/belle/varghese/IPHC/B_Kspipigamma_belle/prepare_recFiles_for_tatFits/outFiles_with_myb2XMLFile'

data_dic = {
'files' : data_files,
'dir' : data_dir,
'weight' : 1.0
}


component_Data = {
'data' : data_dic
}

component_dic = {
'component_MC' : component_MC,
'component_Data' : component_Data
}




#rootfiles
key = 'h1' #name of tree in rootfiles

#general
variables_used = ['issignal', 'b2MVA', 'de', 'mbc', 'pi0veto', 'etaveto', 'g_e9e25', 'g_cos', 'tagcl', 'vtxcl', 'm23', 'x_p0_kid', 'x_p1_kid', 'x_p0_eid', 'x_p1_eid', 'x_m', 'tagh', 'tagntrk', 'vtxh'] # all variables used in program
opt_variable = 'b2MVA' # the variable to be optimized
sig_variable = 'issignal'# the variable that determines the difference between signal and background

#plot labels
labels_sigEff = {
'title' : 'Signal efficiency vs BASF2 CS-MVA cut',
'x' : 'cut value',
'y' : 'efficiency (%)',
'name' : 'sigEff.png'
}


labels_fom = {
'title' : 'Figure-of-merit vs BASF2 CS-MVA cut',
'x' : 'cut value',
'y' : 'f.o.m.',
'name' : 'fom.png'
}


labels_dmc = {
'title' : 'Data/MC background ratio vs BASF2 CS-MVA cut',
'x' : 'cut value',
'y' : 'ratio',
'name' : 'data_by_mc_ratio.png'
}



labels = {
'sigEff' : labels_sigEff,
'fom' : labels_fom,
'data_mc_bkg_ratio' : labels_dmc
}

	
	
##################################################################################################
# Main program
##################################################################################################

if __name__ == '__main__':	


	#get various metrics	
	cut, sigeff, fom, opt_cut, data_mc_bkg_ratio = get_metrics(component_dic, key, variables_used, opt_variable, sig_variable, selections, cut_nbins, cut_start, cut_width, doDataMC_bkg_correction)
	
	#plot calculate optimization metrics
	plot_metrics(cut, sigeff, fom, data_mc_bkg_ratio, labels, opt_cut, opt_variable, doDataMC_bkg_correction)



##################################################################################################




