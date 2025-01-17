#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to make toys dataframes to check fit systematics
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


#load dataframes
df_gmc = uproot.open(gmc_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
df_gmc = shuffle(df_gmc)
df_rare = uproot.open(rare_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
df_rare_bkg = df_rare.query('issgnevt != 1')


data_gmc = {key: df_gmc[key].values for key in variables}
rdf_gmc = r.RDF.MakeNumpyDataFrame(data_gmc)
rdf_gmc.Snapshot('h1', 'df_gmc.root')
data_rare_bkg = {key: df_rare_bkg[key].values for key in variables}
rdf_rare_bkg = r.RDF.MakeNumpyDataFrame(data_rare_bkg)
rdf_rare_bkg.Snapshot('h1', 'df_rare_bkg.root')











