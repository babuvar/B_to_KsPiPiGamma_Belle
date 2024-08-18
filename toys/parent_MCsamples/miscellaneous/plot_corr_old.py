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
#BB
bb_mbc_de = TH2F('bb_mbc_de', 'B #bar{B} M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); bb_list.append(bb_mbc_de)
bb_mbc_dt = TH2F('bb_mbc_dt', 'B #bar{B} M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); bb_list.append(bb_mbc_dt)
bb_de_dt = TH2F('bb_de_dt', 'B #bar{B} #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); bb_list.append(bb_de_dt)
#Continuum
cont_mbc_de = TH2F('cont_mbc_de', 'Continuum M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); cont_list.append(cont_mbc_de)
cont_mbc_dt = TH2F('cont_mbc_dt', 'Continuum M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); cont_list.append(cont_mbc_dt)
cont_de_dt = TH2F('cont_de_dt', 'Continuum #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); cont_list.append(cont_de_dt)
#rare bkg
rare_mbc_de = TH2F('rare_mbc_de', 'Rare M_{bc} - #Delta E; M_{bc} (GeV); #Delta E (GeV)', 50, 5.2, 5.3, 50, -0.2, 0.2); rare_list.append(rare_mbc_de)
rare_mbc_dt = TH2F('rare_mbc_dt', 'Rare M_{bc} - #Delta t; M_{bc} (GeV); #Delta t (ps)', 50, 5.2, 5.3, 50, -10.0, 10.0); rare_list.append(rare_mbc_dt)
rare_de_dt = TH2F('rare_de_dt', 'Rare #Delta E - #Delta t; #Delta E (GeV); #Delta t (ps)', 50, -0.2, 0.2, 50, -10.0, 10.0); rare_list.append(rare_de_dt)

component_list = []; component_list.append(sig_list); component_list.append(scf_list); component_list.append(bb_list); component_list.append(cont_list); component_list.append(rare_list);



for i in range(chain.GetEntries()):

	chain.GetEntry(i)
	
	#check that things are in given range
	if chain.mbc > 5.2 and chain.mbc < 5.3 and chain.de > -0.2 and chain.de < 0.2 and chain.dt > -10.0 and chain.dt < 10.0:

		#Fill signal
		if chain.issignal == 1:
			sig_mbc_de.Fill(chain.mbc, chain.de)
			sig_mbc_dt.Fill(chain.mbc, chain.dt)
			sig_de_dt.Fill(chain.de, chain.dt)

		#Fill SCF
		if chain.issignal != 1 and chain.issgnevt == 1:
			scf_mbc_de.Fill(chain.mbc, chain.de)
			scf_mbc_dt.Fill(chain.mbc, chain.dt)
			scf_de_dt.Fill(chain.de, chain.dt)

		#Fill BB
		if chain.mc_type == 10 or chain.mc_type == 30:
			bb_mbc_de.Fill(chain.mbc, chain.de)
			bb_mbc_dt.Fill(chain.mbc, chain.dt)
			bb_de_dt.Fill(chain.de, chain.dt)

		#Fill continuum
		if chain.mc_type == 20 or chain.mc_type == 40:
			cont_mbc_de.Fill(chain.mbc, chain.de)
			cont_mbc_dt.Fill(chain.mbc, chain.dt)
			cont_de_dt.Fill(chain.de, chain.dt)

		#Fill Rare bkg
		if chain.mc_type == 0 and chain.issgnevt != 1:
			rare_mbc_de.Fill(chain.mbc, chain.de)
			rare_mbc_dt.Fill(chain.mbc, chain.dt)
			rare_de_dt.Fill(chain.de, chain.dt)

#correlation coeffs
corr_sig_mbc_de = sig_mbc_de.GetCorrelationFactor(); sig_mbc_de.SetTitle('Signal correlation = %.3f'%corr_sig_mbc_de)
corr_sig_mbc_dt = sig_mbc_dt.GetCorrelationFactor(); sig_mbc_dt.SetTitle('Signal correlation = %.3f'%corr_sig_mbc_dt)
corr_sig_de_dt = sig_de_dt.GetCorrelationFactor(); sig_de_dt.SetTitle('Signal correlation = %.3f'%corr_sig_de_dt)
corr_scf_mbc_de = scf_mbc_de.GetCorrelationFactor(); scf_mbc_de.SetTitle('Cross-feed correlation = %.3f'%corr_scf_mbc_de)
corr_scf_mbc_dt = scf_mbc_dt.GetCorrelationFactor(); scf_mbc_dt.SetTitle('Cross-feed correlation = %.3f'%corr_scf_mbc_dt)
corr_scf_de_dt = scf_de_dt.GetCorrelationFactor(); scf_de_dt.SetTitle('Cross-feed correlation = %.3f'%corr_scf_de_dt)
corr_bb_mbc_de = bb_mbc_de.GetCorrelationFactor(); bb_mbc_de.SetTitle('B #bar{B} correlation = %.3f'%corr_bb_mbc_de)
corr_bb_mbc_dt = bb_mbc_dt.GetCorrelationFactor(); bb_mbc_dt.SetTitle('B #bar{B} correlation = %.3f'%corr_bb_mbc_dt)
corr_bb_de_dt = bb_de_dt.GetCorrelationFactor(); bb_de_dt.SetTitle('B #bar{B} correlation = %.3f'%corr_bb_de_dt)
corr_cont_mbc_de = cont_mbc_de.GetCorrelationFactor(); cont_mbc_de.SetTitle('Continuum correlation = %.3f'%corr_cont_mbc_de)
corr_cont_mbc_dt = cont_mbc_dt.GetCorrelationFactor(); cont_mbc_dt.SetTitle('Continuum correlation = %.3f'%corr_cont_mbc_dt)
corr_cont_de_dt = cont_de_dt.GetCorrelationFactor(); cont_de_dt.SetTitle('Continuum correlation = %.3f'%corr_cont_de_dt)
corr_rare_mbc_de = rare_mbc_de.GetCorrelationFactor(); rare_mbc_de.SetTitle('Rare-bkg correlation = %.3f'%corr_rare_mbc_de)
corr_rare_mbc_dt = rare_mbc_dt.GetCorrelationFactor(); rare_mbc_dt.SetTitle('Rare-bkg correlation = %.3f'%corr_rare_mbc_dt)
corr_rare_de_dt = rare_de_dt.GetCorrelationFactor(); rare_de_dt.SetTitle('Rare-bkg correlation = %.3f'%corr_rare_de_dt)


#plot histograms
for i, (h_list, name) in enumerate(zip(component_list, names)):
	plot(h_list, name)












