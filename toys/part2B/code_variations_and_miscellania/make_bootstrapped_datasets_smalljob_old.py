#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to make toys dataframes to check fit systematics
This script creates toys datasets with SCF from the private-SMC.
Note that the SCF component looks different for MC-rare and private-SMC.
"""


import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import json
import subprocess
import shlex
import uproot
import sys

def shuffle(df):
	df = df.reindex(np.random.permutation(df.index))
	df.reset_index(drop=True, inplace=True)
	return df

def split_prop(var):
	if isinstance(var, str):
		return (var, {})
	else:
		return var

variables = ['nexp', 'nrun', 'eventid', 'ecm', 'ebeam', 'm', 'mbc', 'de', 'p', 'mcflags', 'ecms', 'pcms', 'costheta', 'coscms', 'mcpdg', 'b0pgen', 'b0bpgen', 'mcp', 'issignal', 'issgnevt', 'sigmc', 'sigk1ks', 'sigk1k0s', 'sigk1k2p', 'sigk1rho', 'sigk2ks', 'sigk2rho', 'sigk1kpp', 'sigk1omg', 'sigk1f0', 'issgnang', 'issgngam', 'ncands', 'ncand', 'a_ro_pi1', 'a_ro_pi2', 'a_ks_pi1', 'a_ks_pi2', 'a_pi_max', 'p1mc_px', 'p1mc_py', 'p1mc_pz', 'p1mc_e', 'p2mc_px', 'p2mc_py', 'p2mc_pz', 'p2mc_e', 'p3mc_px', 'p3mc_py', 'p3mc_pz', 'p3mc_e', 'p4mc_px', 'p4mc_py', 'p4mc_pz', 'p4mc_e', 'p1_px', 'p1_py', 'p1_pz', 'p1_e', 'p2_px', 'p2_py', 'p2_pz', 'p2_e', 'p3_px', 'p3_py', 'p3_pz', 'p3_e', 'p4_px', 'p4_py', 'p4_pz', 'p4_e', 'qr', 'w', 'dw', 'wtag', 'flavor', 'rbin', 'genbq', 'zrec', 'vtxcl', 'vtxh', 'vtxchi2', 'vtxndf', 'vtxntrk', 'vtxzerr', 'vtxyerr', 'vtxxerr', 'ztag', 'tagcl', 'tagh', 'tagchi2', 'tagchi2w', 'tagndf', 'tagndfwo', 'taglepto', 'tagzerr', 'tagyerr', 'tagxerr', 'tagntrk', 'm12', 'm23', 'm13', 'dt', 'dtgen', 'dtres', 'dtpull', 'dterr', 'trecgen', 'ttaggen', 'dz', 'dzgen', 'dzres', 'dzpull', 'dzerr', 'zrecgen', 'ztaggen', 'x_p', 'x_pcms', 'x_issign', 'x_m', 'x_ks_p', 'x_ks_m', 'x_ks_pcm', 'x_ks_cos', 'g_p', 'g_pcms', 'g_issign', 'g_cos', 'x_cos', 'g_e9e25', 'x_p0_p', 'x_p1_p', 'x_p0_cos', 'x_p1_cos', 'x_p0_kid', 'x_p1_kid', 'x_p0_eid', 'x_p1_eid', 'x_p0_mid', 'x_p1_mid', 'x_p0_cdc', 'x_p1_cdc', 'x_p0_svd', 'x_p1_svd', 'k0mm2', 'k0et', 'k0hso00', 'k0hso01', 'k0hso02', 'k0hso03', 'k0hso04', 'k0hso10', 'k0hso12', 'k0hso14', 'k0hso20', 'k0hso22', 'k0hso24', 'k0hoo0', 'k0hoo1', 'k0hoo2', 'k0hoo3', 'k0hoo4', 'costh', 'pi0veto', 'etaveto', 'rMVA', 'b2MVA', 'mc_type', 'inexp', 'irbin', 'iflavor', 'dalitzCategory']

gmc_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/genericMC_x6_wRMVA_wB2MVA_sel_BCS_branches.root'
rare_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches.root'

toy_dir = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/part2B/bootstrapped_datasets'

#deal with indices of the job
SMC_sample = sys.argv[3]
#smc_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/SMC_%s_1M_toFit_wRMVA_wB2MVA_sel_BCS_branches.root'%SMC_sample
smc_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/%s_1M_correct_toFit.root'%SMC_sample
number_toys = int(sys.argv[1])
script_index = int(sys.argv[2])
toy_index_lowlimit = script_index * number_toys
toy_index_highlimit = (script_index + 1) * number_toys
toy_indices = [*range(toy_index_lowlimit, toy_index_highlimit)]

#load dataframes
df_gmc = uproot.open(gmc_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
df_gmc = shuffle(df_gmc)
df_rare = uproot.open(rare_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
df_rare_bkg = df_rare.query('issgnevt != 1')
df_smc = uproot.open(smc_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')

# estimate true - yields
true_yields = {'gmc_yield' : len(df_gmc)/6.0,
               'rare_bkg_yield' : len(df_rare_bkg)/50.0,
               'smc_yield' : (len(df_rare) - len(df_rare_bkg)) / 50.0}

print(true_yields)





# poissonian distribution for toy-ensemble yields
dataset_gmc_yields = np.random.poisson(true_yields['gmc_yield'], number_toys)
dataset_rarebkg_yields = np.random.poisson(true_yields['rare_bkg_yield'], number_toys)
dataset_smc_yields = np.random.poisson(true_yields['smc_yield'], number_toys)


#generate each toy-dataset
for i, (gmc_yield, rarebkg_yield, smc_yield) in enumerate(zip(dataset_gmc_yields, dataset_rarebkg_yields, dataset_smc_yields)):

	df_toy_gmc = df_gmc.sample(gmc_yield, replace=True)
	df_toy_rarebkg = df_rare_bkg.sample(rarebkg_yield, replace=True) 
	df_toy_smc = df_smc.sample(smc_yield, replace=True)
	
	df_toy = pd.concat([df_toy_gmc, df_toy_rarebkg, df_toy_smc], ignore_index=True)
	df_toy = shuffle(df_toy)

	data = {key: df_toy[key].values for key in variables}
	rdf = r.RDF.MakeNumpyDataFrame(data)
	rdf.Snapshot('h1', '%s/toy_dataset_%s_%i.root'%(toy_dir, SMC_sample,toy_indices[i]))
	print('done ', i)










