import uproot

def split_prop(var):
	if isinstance(var, str):
		return (var, {})
	else:
		return var

variables = ['nexp', 'nrun', 'eventid', 'ecm', 'ebeam', 'm', 'mbc', 'de', 'p', 'mcflags', 'ecms', 'pcms', 'costheta', 'coscms', 'mcpdg', 'b0pgen', 'b0bpgen', 'mcp', 'issignal', 'issgnevt', 'sigmc', 'sigk1ks', 'sigk1k0s', 'sigk1k2p', 'sigk1rho', 'sigk2ks', 'sigk2rho', 'sigk1kpp', 'sigk1omg', 'sigk1f0', 'issgnang', 'issgngam', 'ncands', 'ncand', 'a_ro_pi1', 'a_ro_pi2', 'a_ks_pi1', 'a_ks_pi2', 'a_pi_max', 'p1mc_px', 'p1mc_py', 'p1mc_pz', 'p1mc_e', 'p2mc_px', 'p2mc_py', 'p2mc_pz', 'p2mc_e', 'p3mc_px', 'p3mc_py', 'p3mc_pz', 'p3mc_e', 'p4mc_px', 'p4mc_py', 'p4mc_pz', 'p4mc_e', 'p1_px', 'p1_py', 'p1_pz', 'p1_e', 'p2_px', 'p2_py', 'p2_pz', 'p2_e', 'p3_px', 'p3_py', 'p3_pz', 'p3_e', 'p4_px', 'p4_py', 'p4_pz', 'p4_e', 'qr', 'w', 'dw', 'wtag', 'flavor', 'rbin', 'genbq', 'zrec', 'vtxcl', 'vtxh', 'vtxchi2', 'vtxndf', 'vtxntrk', 'vtxzerr', 'vtxyerr', 'vtxxerr', 'ztag', 'tagcl', 'tagh', 'tagchi2', 'tagchi2w', 'tagndf', 'tagndfwo', 'taglepto', 'tagzerr', 'tagyerr', 'tagxerr', 'tagntrk', 'm12', 'm23', 'm13', 'dt', 'dtgen', 'dtres', 'dtpull', 'dterr', 'trecgen', 'ttaggen', 'dz', 'dzgen', 'dzres', 'dzpull', 'dzerr', 'zrecgen', 'ztaggen', 'x_p', 'x_pcms', 'x_issign', 'x_m', 'x_ks_p', 'x_ks_m', 'x_ks_pcm', 'x_ks_cos', 'g_p', 'g_pcms', 'g_issign', 'g_cos', 'x_cos', 'g_e9e25', 'x_p0_p', 'x_p1_p', 'x_p0_cos', 'x_p1_cos', 'x_p0_kid', 'x_p1_kid', 'x_p0_eid', 'x_p1_eid', 'x_p0_mid', 'x_p1_mid', 'x_p0_cdc', 'x_p1_cdc', 'x_p0_svd', 'x_p1_svd', 'k0mm2', 'k0et', 'k0hso00', 'k0hso01', 'k0hso02', 'k0hso03', 'k0hso04', 'k0hso10', 'k0hso12', 'k0hso14', 'k0hso20', 'k0hso22', 'k0hso24', 'k0hoo0', 'k0hoo1', 'k0hoo2', 'k0hoo3', 'k0hoo4', 'costh', 'pi0veto', 'etaveto', 'rMVA', 'b2MVA', 'mc_type', 'inexp', 'irbin', 'iflavor', 'dalitzCategory']



rare_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches.root'

gmc_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/genericMC_x6_wRMVA_wB2MVA_sel_BCS_branches.root'

smc_rootfile = '/gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/toys/parent_MCsamples/S0p0_1M_correct_wRMVA_wB2MVA_sel_BCS_branches.root'

#get rare and generic df's
df_rare = uproot.open(rare_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
df_gmc = uproot.open(gmc_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')
df_smc = uproot.open(smc_rootfile)['h1'].arrays([split_prop(i)[0] for i in variables], library='pd')


#df_sig = df_smc.query('(issignal == 1 | mcflags == 10) & (irbin > 0) & (abs(dt) < 10.0)')
#df_scf = df_smc.query('(issignal != 1 & mcflags != 10 ) & (irbin > 0) & (abs(dt) < 10.0)')
#df_cont = df_gmc.query('(mc_type == 20 | mc_type == 40) & (irbin > 0) & (abs(dt) < 10.0)')
#df_bb = df_gmc.query('(mc_type == 10 | mc_type == 30) & (irbin > 0) & (abs(dt) < 10.0)')
#df_rarebkg = df_rare.query('(issgnevt!= 1) & irbin > 0')
df_sig = df_smc.query('(issignal == 1 | mcflags == 10) & (dt < 10.0 & dt > -10.0)')
df_scf = df_smc.query('(issignal != 1 & mcflags != 10 ) & (dt < 10.0 & dt > -10.0)')
df_cont = df_gmc.query('(mc_type == 20 | mc_type == 40) & (dt < 10.0 & dt > -10.0)')
df_bb = df_gmc.query('(mc_type == 10 | mc_type == 30) & (dt < 10.0 & dt > -10.0)')
df_rarebkg = df_rare.query('(issgnevt!= 1) & (dt < 10.0 & dt > -10.0)')


df_rareRandom = df_rarebkg.query('mcpdg != 511 & mcpdg != -511')
df_rareMissFSP = df_rarebkg.query('mcpdg == 511 | mcpdg == -511')
df_raresig = df_rare.query('(issgnevt== 1) & irbin > 0')

n_sig = len(df_sig)
n_scf = len(df_scf)
n_rarebkg = len(df_rarebkg)/50.0
n_cont = len(df_cont)/6.0
n_BB = len(df_bb)/6.0

n_BBmissFSP = len(df_rareMissFSP)/50.0
n_BBRand = (len(df_rareRandom)/50.0) + (len(df_bb)/6.0)


f_scf = float(n_scf/(n_sig+n_scf))
f_BBmissFSP = float(n_BBmissFSP/(n_BBmissFSP+n_BBRand+n_cont))
f_BBRand = float(n_BBRand/(n_BBmissFSP+n_BBRand+n_cont))

print('signal (unnormalized) =', n_sig)
print('scf (unnormalized) =', n_scf)
print('rare(sig+scf) =',len(df_raresig)/50)
print('rare-bkg  (normalized) =', n_rarebkg)
print('continuum (normalized) =', n_cont)
print('BB (normalized) =', n_BB)
print('BBmissFSP (normalized) =', n_BBmissFSP)
print('BBRand (normalized) =', n_BBRand)
print('===========')
print('f_scf = ',f_scf)
print('f_BBmissFSP = ',f_BBmissFSP)
print('f_BBRand = ',f_BBRand)
print('===========')
