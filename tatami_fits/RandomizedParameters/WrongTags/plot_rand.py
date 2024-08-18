from ROOT import TFile, RooWorkspace, RooPlot, RooFit, TCanvas, RooRealVar

f_randWtags = TFile('randWTagMetrics.root')
rws_randWtags = f_randWtags.Get('rws_randWTags')
data_w10  = rws_randWtags.data('w10')


#Plot
#w10 = RooRealVar('w10', 'w10', -2.0, 2.0)
#frame = w10.frame(RooFit.Bins(40), RooFit.Title(""))
#data_w10.plotOn(xframe)

#can = TCanvas('can', 'can', 800, 800)
#can.cd()
#frame.Draw()
#can.SaveAs('plot.png')


