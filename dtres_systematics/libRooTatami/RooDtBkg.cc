/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

// Your description goes here... 

#include "Riostream.h" 

#include "RooDtBkg.h"
#include "RooAbsReal.h" 
#include "RooCategory.h"
#include "RooRealConstant.h"
#include <math.h> 
#include "TMath.h"

#include "tatami/tatami.h"

ClassImp(RooDtBkg)

RooDtBkg::RooDtBkg(const char *name, const char *title,
		RooAbsReal& _dt,
		RooAbsReal& _tau_bkg,
		RooAbsReal& _mu_l,
		RooAbsReal& _mu_d,
		RooAbsReal& _f_delt,
		RooAbsReal& _f_tail,
		RooAbsReal& _S_main,
		RooAbsReal& _S_tail,
		RooCategory& _expno,
		RooAbsReal& _rec_vtntrk,
		RooAbsReal& _rec_vterr,
		RooAbsReal& _asc_vtntrk,
		RooAbsReal& _asc_vterr,
		bool _mc,
		bool _addoutlier,
		double _alpha,
		double _beta) :
		RooAbsPdf(name, title),
		dt("dt", "dt", this, _dt),
		tau_bkg("tau_bkg", "tau_bkg", this, _tau_bkg),
		mu_l("mu_l", "mu_l", this, _mu_l),
		mu_d("mu_d", "mu_d", this, _mu_d),
		f_delt("f_delt", "f_delt", this, _f_delt),
		f_tail("f_tail", "f_tail", this, _f_tail),
		S_main("S_main", "S_main", this, _S_main),
		S_tail("S_tail", "S_tail", this, _S_tail),
		f_delt_mlt("f_delt_mlt", "f_delt_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		f_tail_mlt("f_tail_mlt", "f_tail_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		S_main_mlt("S_main_mlt", "S_main_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		S_tail_mlt("S_tail_mlt", "S_tail_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		expno("expno", "expno", this, _expno),
		rec_vtntrk("rec_vtntrk", "rec_vtntrk", this, _rec_vtntrk),
		rec_vterr("rec_vterr", "rec_vterr", this, _rec_vterr),
		asc_vtntrk("asc_vtntrk", "asc_vtntrk", this, _asc_vtntrk),
		asc_vterr("asc_vterr", "asc_vterr", this, _asc_vterr),
		mc(_mc),
		dt_ll(Belle::dt_resol_global::dt_llmt),
		dt_ul(Belle::dt_resol_global::dt_ulmt),
		addoutlier(_addoutlier),
		alpha(_alpha),
		beta(_beta) {
}

RooDtBkg::RooDtBkg(const char *name, const char *title,
		RooAbsReal& _dt,
		RooAbsReal& _tau_bkg,
		RooAbsReal& _mu_l,
		RooAbsReal& _mu_d,
		RooAbsReal& _f_delt,
		RooAbsReal& _f_tail,
		RooAbsReal& _S_main,
		RooAbsReal& _S_tail,
		RooAbsReal& _f_delt_mlt,
		RooAbsReal& _f_tail_mlt,
		RooAbsReal& _S_main_mlt,
		RooAbsReal& _S_tail_mlt,
		RooCategory& _expno,
		RooAbsReal& _rec_vtntrk,
		RooAbsReal& _rec_vterr,
		RooAbsReal& _asc_vtntrk,
		RooAbsReal& _asc_vterr,
		bool _mc,
		bool _addoutlier,
		double _alpha,
		double _beta) :
		RooAbsPdf(name, title),
		dt("dt", "dt", this, _dt),
		tau_bkg("tau_bkg", "tau_bkg", this, _tau_bkg),
		mu_l("mu_l", "mu_l", this, _mu_l),
		mu_d("mu_d", "mu_d", this, _mu_d),
		f_delt("f_delt", "f_delt", this, _f_delt),
		f_tail("f_tail", "f_tail", this, _f_tail),
		S_main("S_main", "S_main", this, _S_main),
		S_tail("S_tail", "S_tail", this, _S_tail),
		f_delt_mlt("f_delt_mlt", "f_delt_mlt", this, _f_delt_mlt),
		f_tail_mlt("f_tail_mlt", "f_tail_mlt", this, _f_tail_mlt),
		S_main_mlt("S_main_mlt", "S_main_mlt", this, _S_main_mlt),
		S_tail_mlt("S_tail_mlt", "S_tail_mlt", this, _S_tail_mlt),
		expno("expno", "expno", this, _expno),
		rec_vtntrk("rec_vtntrk", "rec_vtntrk", this, _rec_vtntrk),
		rec_vterr("rec_vterr", "rec_vterr", this, _rec_vterr),
		asc_vtntrk("asc_vtntrk", "asc_vtntrk", this, _asc_vtntrk),
		asc_vterr("asc_vterr", "asc_vterr", this, _asc_vterr),
		mc(_mc),
		dt_ll(Belle::dt_resol_global::dt_llmt),
		dt_ul(Belle::dt_resol_global::dt_ulmt),
		addoutlier(_addoutlier),
		alpha(_alpha),
		beta(_beta) {
}

RooDtBkg::RooDtBkg(const char *name, const char *title,
		RooAbsReal& _dt,
		RooAbsReal& _tau_bkg,
		RooAbsReal& _mu_l,
		RooAbsReal& _mu_d,
		RooAbsReal& _f_delt,
		RooAbsReal& _f_tail,
		RooAbsReal& _S_main,
		RooAbsReal& _S_tail,
		RooCategory& _expno,
		RooAbsReal& _rec_vtntrk,
		RooAbsReal& _rec_vterr,
		RooAbsReal& _asc_vtntrk,
		RooAbsReal& _asc_vterr,
		bool _mc,
		double _dt_ll,
		double _dt_ul,
		bool _addoutlier,
		double _alpha,
		double _beta) :
		RooAbsPdf(name, title),
		dt("dt", "dt", this, _dt),
		tau_bkg("tau_bkg", "tau_bkg", this, _tau_bkg),
		mu_l("mu_l", "mu_l", this, _mu_l),
		mu_d("mu_d", "mu_d", this, _mu_d),
		f_delt("f_delt", "f_delt", this, _f_delt),
		f_tail("f_tail", "f_tail", this, _f_tail),
		S_main("S_main", "S_main", this, _S_main),
		S_tail("S_tail", "S_tail", this, _S_tail),
		f_delt_mlt("f_delt_mlt", "f_delt_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		f_tail_mlt("f_tail_mlt", "f_tail_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		S_main_mlt("S_main_mlt", "S_main_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		S_tail_mlt("S_tail_mlt", "S_tail_mlt", this, (RooAbsReal&) RooRealConstant::value(-999)),
		expno("expno", "expno", this, _expno),
		rec_vtntrk("rec_vtntrk", "rec_vtntrk", this, _rec_vtntrk),
		rec_vterr("rec_vterr", "rec_vterr", this, _rec_vterr),
		asc_vtntrk("asc_vtntrk", "asc_vtntrk", this, _asc_vtntrk),
		asc_vterr("asc_vterr", "asc_vterr", this, _asc_vterr),
		mc(_mc),
		dt_ll(_dt_ll),
		dt_ul(_dt_ul),
		addoutlier(_addoutlier),
		alpha(_alpha),
		beta(_beta) {
}

RooDtBkg::RooDtBkg(const char *name, const char *title,
		RooAbsReal& _dt,
		RooAbsReal& _tau_bkg,
		RooAbsReal& _mu_l,
		RooAbsReal& _mu_d,
		RooAbsReal& _f_delt,
		RooAbsReal& _f_tail,
		RooAbsReal& _S_main,
		RooAbsReal& _S_tail,
		RooAbsReal& _f_delt_mlt,
		RooAbsReal& _f_tail_mlt,
		RooAbsReal& _S_main_mlt,
		RooAbsReal& _S_tail_mlt,
		RooCategory& _expno,
		RooAbsReal& _rec_vtntrk,
		RooAbsReal& _rec_vterr,
		RooAbsReal& _asc_vtntrk,
		RooAbsReal& _asc_vterr,
		bool _mc,
		double _dt_ll,
		double _dt_ul,
		bool _addoutlier,
		double _alpha,
		double _beta) :
		RooAbsPdf(name, title),
		dt("dt", "dt", this, _dt),
		tau_bkg("tau_bkg", "tau_bkg", this, _tau_bkg),
		mu_l("mu_l", "mu_l", this, _mu_l),
		mu_d("mu_d", "mu_d", this, _mu_d),
		f_delt("f_delt", "f_delt", this, _f_delt),
		f_tail("f_tail", "f_tail", this, _f_tail),
		S_main("S_main", "S_main", this, _S_main),
		S_tail("S_tail", "S_tail", this, _S_tail),
		f_delt_mlt("f_delt_mlt", "f_delt_mlt", this, _f_delt_mlt),
		f_tail_mlt("f_tail_mlt", "f_tail_mlt", this, _f_tail_mlt),
		S_main_mlt("S_main_mlt", "S_main_mlt", this, _S_main_mlt),
		S_tail_mlt("S_tail_mlt", "S_tail_mlt", this, _S_tail_mlt),
		expno("expno", "expno", this, _expno),
		rec_vtntrk("rec_vtntrk", "rec_vtntrk", this, _rec_vtntrk),
		rec_vterr("rec_vterr", "rec_vterr", this, _rec_vterr),
		asc_vtntrk("asc_vtntrk", "asc_vtntrk", this, _asc_vtntrk),
		asc_vterr("asc_vterr", "asc_vterr", this, _asc_vterr),
		mc(_mc),
		dt_ll(_dt_ll),
		dt_ul(_dt_ul),
		addoutlier(_addoutlier),
		alpha(_alpha),
		beta(_beta) {
}

RooDtBkg::RooDtBkg(const RooDtBkg& other, const char* name) :
		RooAbsPdf(other, name),
		dt("dt", this, other.dt),
		tau_bkg("tau_bkg", this, other.tau_bkg),
		mu_l("mu_l", this, other.mu_l),
		mu_d("mu_d", this, other.mu_d),
		f_delt("f_delt", this, other.f_delt),
		f_tail("f_tail", this, other.f_tail),
		S_main("S_main", this, other.S_main),
		S_tail("S_tail", this, other.S_tail),
		f_delt_mlt("f_delt_mlt", this, other.f_delt_mlt),
		f_tail_mlt("f_tail_mlt", this, other.f_tail_mlt),
		S_main_mlt("S_main_mlt", this, other.S_main_mlt),
		S_tail_mlt("S_tail_mlt", this, other.S_tail_mlt),
		expno("expno", this, other.expno),
		rec_vtntrk("rec_vtntrk", this, other.rec_vtntrk),
		rec_vterr("rec_vterr", this, other.rec_vterr),
		asc_vtntrk("asc_vtntrk", this, other.asc_vtntrk),
		asc_vterr("asc_vterr", this, other.asc_vterr),
		mc(other.mc),
		dt_ll(other.dt_ll),
		dt_ul(other.dt_ul),
		addoutlier(other.addoutlier),
		alpha(other.alpha),
		beta(other.beta) {
}

Double_t RooDtBkg::evaluate() const {
	double in_S_main, in_S_tail, in_f_delt, in_f_tail;
	if (S_main_mlt == -999) {
		in_S_main = S_main;
		in_S_tail = S_tail;
		in_f_delt = f_delt;
		in_f_tail = f_tail;
	} else {
		if (rec_vtntrk > 1 && asc_vtntrk > 1) {
			in_S_main = S_main_mlt;
			in_S_tail = S_tail_mlt;
			in_f_delt = f_delt_mlt;
			in_f_tail = f_tail_mlt;
		} else 	{
			in_S_main = S_main;
			in_S_tail = S_tail;
			in_f_delt = f_delt;
			in_f_tail = f_tail;
		}
	}

	const Belle::dtres_param_t* dtres_param = Belle::get_dtres_param(expno, (int) mc);

	double m_s_main = in_S_main * sqrt(rec_vterr * rec_vterr + asc_vterr * asc_vterr) * Belle::dt_resol_global::inv_bgc;
	double m_s_tail = in_S_tail * m_s_main;

	double prompt_part_main = Belle::gaussian(dt, mu_d, m_s_main);
	prompt_part_main /= Belle::norm_gaussian(dt_ll, dt_ul, mu_d, m_s_main);

	double prompt_part_tail = Belle::gaussian(dt, mu_d, m_s_tail);
	prompt_part_tail /= Belle::norm_gaussian(dt_ll, dt_ul, mu_d, m_s_tail);

	double prompt_part = (1 - in_f_tail) * prompt_part_main + in_f_tail * prompt_part_tail;

	double lifetime_part_main = Belle::Enp_conv_gauss(dt, tau_bkg, tau_bkg, mu_l, m_s_main);
	lifetime_part_main /= Belle::norm_Enp_conv_gauss(dt_ll, dt_ul, tau_bkg, tau_bkg, mu_l, m_s_main);

	double lifetime_part_tail = Belle::Enp_conv_gauss(dt, tau_bkg, tau_bkg, mu_l, m_s_tail);
	lifetime_part_tail /= Belle::norm_Enp_conv_gauss(dt_ll, dt_ul, tau_bkg, tau_bkg, mu_l, m_s_tail);

	double lifetime_part = (1 - in_f_tail) * lifetime_part_main + in_f_tail * lifetime_part_tail;

	//double pdf = in_f_delt * prompt_part + (1 - in_f_delt) * lifetime_part; //all other case we should use the original PDF
	double pdf = prompt_part; // we use double Gaussian only for J/psi pi0 as used in the previous analysis

	cout<<"fol_mul = "<<dtres_param->fol_mul<<endl;
        cout<<"fol_sgl = "<<dtres_param->fol_sgl<<endl;
        cout<<"sig_ol = "<<dtres_param->sig_ol<<endl;
        cout<<"alpha = "<<alpha<<endl;
        cout<<"beta = "<<beta<<endl;

	if (addoutlier)
		return Belle::AddOutlierWithBkg(expno, dt, 0, 0, pdf, rec_vtntrk, asc_vtntrk, dtres_param, 1 / alpha, 1 / alpha, dt_ll, dt_ul, alpha, beta);
	else
		return pdf;
}

Int_t RooDtBkg::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {
	if (matchArgs(allVars, analVars, dt)) return 1;
	return 0;
}

Double_t RooDtBkg::analyticalIntegral(Int_t code, const char* /*rangeName*/) const {
	assert(code == 1);
	return 1;
}

void RooDtBkg::set_alpha(double in_alpha) {
	alpha = in_alpha;
}
