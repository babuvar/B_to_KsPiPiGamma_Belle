from ROOT import  TH2F, TChain, TCanvas, gROOT, gStyle

gROOT.SetBatch(True)
gStyle.SetOptStat(0)

def plot(h_list, name):

	can = TCanvas('can','',1600, 1600); can.Divide(2,2)
	can.cd(1); h_list[0].Draw('colz')
	can.cd(2); h_list[1].Draw('colz')
	can.cd(3); h_list[2].Draw('colz')
	can.SaveAs('component_correlations/%s.png'%name)

	return

chain = TChain('h1')
chain.Add('MCrare_x50_wRMVA_wB2MVA_sel_BCS_branches.root')
chain.Add('genericMC_x6_wRMVA_wB2MVA_sel_BCS_branches.root')

sig_list = []; scf_list = []; bb_list = []; cont_list = []; rare_list = [];
names = ['Signal', 'SCF', 'BB', 'Continuum', 'RareBkg']


#signal
sig_mbc_de = TH2F('sig_mbc_de', 'Signal M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); sig_list.append(sig_mbc_de)
sig_mbc_dt = TH2F('sig_mbc_dt', 'Signal M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); sig_list.append(sig_mbc_dt)
sig_de_dt = TH2F('sig_de_dt', 'Signal #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); sig_list.append(sig_de_dt)
#SCF
scf_mbc_de = TH2F('scf_mbc_de', 'Cross-feed M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); scf_list.append(scf_mbc_de)
scf_mbc_dt = TH2F('scf_mbc_dt', 'Cross-feed M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); scf_list.append(scf_mbc_dt)
scf_de_dt = TH2F('scf_de_dt', 'Cross-feed #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); scf_list.append(scf_de_dt)
#BB-random
bbrand_mbc_de = TH2F('bbrand_mbc_de', 'B #bar{B} random  M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); bb_list.append(bbrand_mbc_de)
bbrand_mbc_dt = TH2F('bbrand_mbc_dt', 'B #bar{B} random M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); bb_list.append(bbrand_mbc_dt)
bbrand_de_dt = TH2F('bbrand_de_dt', 'B #bar{B} random  #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); bb_list.append(bbrand_de_dt)
#Continuum
cont_mbc_de = TH2F('cont_mbc_de', 'Continuum M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); cont_list.append(cont_mbc_de)
cont_mbc_dt = TH2F('cont_mbc_dt', 'Continuum M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); cont_list.append(cont_mbc_dt)
cont_de_dt = TH2F('cont_de_dt', 'Continuum #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); cont_list.append(cont_de_dt)
#BB-missing FSP
bbmissFSP_mbc_de = TH2F('bbmissFSP_mbc_de', 'B #bar{B} missing-FSP M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); rare_list.append(bbmissFSP_mbc_de)
bbmissFSP_mbc_dt = TH2F('bbmissFSP_mbc_dt', 'B #bar{B} missing-FSP M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); rare_list.append(bbmissFSP_mbc_dt)
bbmissFSP_de_dt = TH2F('bbmissFSP_de_dt', 'B #bar{B} missing-FSP #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); rare_list.append(bbmissFSP_de_dt)

component_list = []; component_list.append(sig_list); component_list.append(scf_list); component_list.append(bb_list); component_list.append(cont_list); component_list.append(rare_list);



for i in range(chain.GetEntries()):

	chain.GetEntry(i)
	
	#check that things are in given range
	if chain.mbc > 5.2 and chain.mbc < 5.3 and chain.de > -0.2 and chain.de < 0.2 and chain.dt > -10.0 and chain.dt < 10.0:

		#appropriate weighting between genericMC-BB and rareMC samples
		weight = 1.0

		#Fill signal
		if chain.issignal == 1 or (chain.issgnevt == 1 and chain.mcflags == 10):
			weight = 1.0/50.0
			sig_mbc_de.Fill(chain.mbc, chain.de, weight)
			sig_mbc_dt.Fill(chain.mbc, chain.dt, weight)
			sig_de_dt.Fill(chain.de, chain.dt, weight)

		#Fill SCF
		if chain.issignal != 1 and chain.issgnevt == 1 and chain.mcflags != 10 :
			weight = 1.0/50.0
			scf_mbc_de.Fill(chain.mbc, chain.de, weight)
			scf_mbc_dt.Fill(chain.mbc, chain.dt, weight)
			scf_de_dt.Fill(chain.de, chain.dt, weight)

		#Fill BB
		if chain.mc_type == 10 or chain.mc_type == 30 or (chain.mc_type == 0 and chain.issgnevt != 1) :
			#appropriate weighting between genericMC-BB and rareMC samples
			weight = 1.0
			if chain.mc_type == 10 or chain.mc_type == 30:
				weight = 1.0/6.0
			if chain.mc_type == 0 and chain.issgnevt != 1:
				weight = 1.0/50.0
			#random
			if abs(chain.mcpdg) != 511:
				bbrand_mbc_de.Fill(chain.mbc, chain.de, weight)
				bbrand_mbc_dt.Fill(chain.mbc, chain.dt, weight)
				bbrand_de_dt.Fill(chain.de, chain.dt, weight)
			#missing FSP
			if abs(chain.mcpdg) == 511:
				bbmissFSP_mbc_de.Fill(chain.mbc, chain.de, weight)
				bbmissFSP_mbc_dt.Fill(chain.mbc, chain.dt, weight)
				bbmissFSP_de_dt.Fill(chain.de, chain.dt, weight)	
		#Fill continuum
		if chain.mc_type == 20 or chain.mc_type == 40:
			weight = 1.0/6.0
			cont_mbc_de.Fill(chain.mbc, chain.de, weight)
			cont_mbc_dt.Fill(chain.mbc, chain.dt, weight)
			cont_de_dt.Fill(chain.de, chain.dt, weight)


#correlation coeffs
corr_sig_mbc_de = sig_mbc_de.GetCorrelationFactor(); sig_mbc_de.SetTitle('Signal correlation = %.3f'%corr_sig_mbc_de)
corr_sig_mbc_dt = sig_mbc_dt.GetCorrelationFactor(); sig_mbc_dt.SetTitle('Signal correlation = %.3f'%corr_sig_mbc_dt)
corr_sig_de_dt = sig_de_dt.GetCorrelationFactor(); sig_de_dt.SetTitle('Signal correlation = %.3f'%corr_sig_de_dt)
corr_scf_mbc_de = scf_mbc_de.GetCorrelationFactor(); scf_mbc_de.SetTitle('Cross-feed correlation = %.3f'%corr_scf_mbc_de)
corr_scf_mbc_dt = scf_mbc_dt.GetCorrelationFactor(); scf_mbc_dt.SetTitle('Cross-feed correlation = %.3f'%corr_scf_mbc_dt)
corr_scf_de_dt = scf_de_dt.GetCorrelationFactor(); scf_de_dt.SetTitle('Cross-feed correlation = %.3f'%corr_scf_de_dt)
corr_bbrand_mbc_de = bbrand_mbc_de.GetCorrelationFactor(); bbrand_mbc_de.SetTitle('B #bar{B}-random correlation = %.3f'%corr_bbrand_mbc_de)
corr_bbrand_mbc_dt = bbrand_mbc_dt.GetCorrelationFactor(); bbrand_mbc_dt.SetTitle('B #bar{B}-random correlation = %.3f'%corr_bbrand_mbc_dt)
corr_bbrand_de_dt = bbrand_de_dt.GetCorrelationFactor(); bbrand_de_dt.SetTitle('B #bar{B}-random correlation = %.3f'%corr_bbrand_de_dt)
corr_cont_mbc_de = cont_mbc_de.GetCorrelationFactor(); cont_mbc_de.SetTitle('Continuum correlation = %.3f'%corr_cont_mbc_de)
corr_cont_mbc_dt = cont_mbc_dt.GetCorrelationFactor(); cont_mbc_dt.SetTitle('Continuum correlation = %.3f'%corr_cont_mbc_dt)
corr_cont_de_dt = cont_de_dt.GetCorrelationFactor(); cont_de_dt.SetTitle('Continuum correlation = %.3f'%corr_cont_de_dt)
corr_bbmissFSP_mbc_de = bbmissFSP_mbc_de.GetCorrelationFactor(); bbmissFSP_mbc_de.SetTitle('B #bar{B} missing-FSP correlation = %.3f'%corr_bbmissFSP_mbc_de)
corr_bbmissFSP_mbc_dt = bbmissFSP_mbc_dt.GetCorrelationFactor(); bbmissFSP_mbc_dt.SetTitle('B #bar{B} missing-FSP correlation = %.3f'%corr_bbmissFSP_mbc_dt)
corr_bbmissFSP_de_dt = bbmissFSP_de_dt.GetCorrelationFactor(); bbmissFSP_de_dt.SetTitle('B #bar{B} missing-FSP correlation = %.3f'%corr_bbmissFSP_de_dt)


#plot histograms
for i, (h_list, name) in enumerate(zip(component_list, names)):
	plot(h_list, name)












