/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROODTBKG
#define ROODTBKG

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooDtBkg : public RooAbsPdf {
public:

	RooDtBkg() : mc(0), dt_ll(0), dt_ul(0), addoutlier(0), alpha(0), beta(0) { };

	RooDtBkg(const char *name, const char *title,
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
			bool _addoutlier = true,
			double _alpha = 1.0,
			double _beta = 1.0);

	RooDtBkg(const char *name, const char *title,
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
			bool _addoutlier = true,
			double _alpha = 1.0,
			double _beta = 1.0);

	RooDtBkg(const char *name, const char *title,
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
			bool _addoutlier = true,
			double _alpha = 1.0,
			double _beta = 1.0);

	RooDtBkg(const char *name, const char *title,
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
			bool _addoutlier = true,
			double _alpha = 1.0,
			double _beta = 1.0);

	RooDtBkg(const RooDtBkg& other, const char* name = 0);

	virtual TObject* clone(const char* newname) const {
		return new RooDtBkg(*this, newname);
	}

	inline virtual ~RooDtBkg() {
	}

	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
	Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

	void set_alpha(double in_alpha);

protected:
	RooRealProxy dt;
	RooRealProxy tau_bkg;
	RooRealProxy mu_l;
	RooRealProxy mu_d;
	RooRealProxy f_delt;
	RooRealProxy f_tail;
	RooRealProxy S_main;
	RooRealProxy S_tail;
	RooRealProxy f_delt_mlt;
	RooRealProxy f_tail_mlt;
	RooRealProxy S_main_mlt;
	RooRealProxy S_tail_mlt;
	RooCategoryProxy expno;
	RooRealProxy rec_vtntrk;
	RooRealProxy rec_vterr;
	RooRealProxy asc_vtntrk;
	RooRealProxy asc_vterr;
	const bool mc;
	const double dt_ll;
	const double dt_ul;
	bool addoutlier;
	double alpha;
	double beta;


	Double_t evaluate() const;

private:

	ClassDef(RooDtBkg, 1) // Your description goes here...
};

#endif
