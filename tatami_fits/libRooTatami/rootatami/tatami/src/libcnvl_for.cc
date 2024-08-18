//
//    libcnvl_for.cc
//       ---- Interface to fortran
//    $Id: libcnvl_for.cc 9932 2006-11-12 14:26:53Z katayama $
//
//    $Log$
//    Revision 1.2  2004/10/19 07:09:55  kohji
//    1) SVD1 and SVD2 resolution parameters are updated for 2004 publication.
//    2) "expno" is added to add_outlier and dtres_systematics functions to treat
//    difference between SVD1 and SVD2 properly[cpfit_ml:0771].
//    3) New functions are added for cosh/sinh terms
//     (norm_)HAfRkRdetRnp_{fullrec,partial} for HAf=exp(-|t|Gamma)sinh(DGamma/2 t)
//     (norm_)HMfRkRdetRnp_{fullrec,partial} for HMf=exp(-|t|Gamma)cosh(DGamma/2 t)
//
//    Revision 1.1  2002/09/05 01:27:04  katayama
//    New package tatami from Nakadaira/Sumisawa san
//
//    Revision 1.4  2002/07/05 06:16:37  nakadair
//    Add interfece to class dt_resol_global.
//
//    Revision 1.3  2002/05/14 04:11:47  nakadair
//    Update for GCC3.
//
//    Revision 1.2  2002/05/09 07:08:35  nakadair
//    Add add_outlier(...).
//
//    Revision 1.1  2002/05/07 21:39:21  nakadair
//    Add convolution with resolution function and fortran interface.
//
//
#include "belle.h"
#include "tatami/libcnvl.h"
#include "tatami/dt_resolution.h"
#include "tatami/conv_coef.h"
#include "tatami/conv_coef3.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



__LINKAGE__

double diracdelta_(double* x){
  return DiracDelta(*x);
}

double gaussian_(double* x, double* m, double* s){
  return gaussian(*x, *m, *s);
}

double norm_gaussian_(double* ll, double* ul, double* m, double* s){
  return norm_gaussian(*ll, *ul, *m, *s);
}

double norm_gaussian_w_cutoff_(double* cutoff, double* m, double* s){
  return norm_gaussian_w_cutoff(*cutoff, *m, *s);
}

double ep_(double* t, double* tau){
  return Ep(*t, *tau);
}

double en_(double* t, double* tau){
  return En(*t, *tau);
}

double ef_(double* t, double* tau){
  return Ef(*t, *tau);
}

double enp_(double* t, double* tau_n, double* tau_p){
  return Enp(*t, *tau_n, *tau_p);
}

double nep_conv_gauss_(double* t, double* m, double* s){
  return nEp_conv_gauss(*t, *m, *s);
}

double nen_conv_gauss_(double* t, double* m, double* s){
  return nEn_conv_gauss(*t, *m, *s);
}

double ep_conv_gauss_(double* t, double* tau, double* m, double* s){
  return Ep_conv_gauss(*t, *tau, *m, *s);
}

double en_conv_gauss_(double* t, double* tau, double* m, double* s){
  return En_conv_gauss(*t, *tau, *m, *s);
}

double ef_conv_gauss_(double* t, double* tau, double* m, double* s){
  return Ef_conv_gauss(*t, *tau, *m, *s);
}

double enp_conv_gauss_(double* t, double* tau_n, double* tau_p,
		       double* m, double* s){
  return Enp_conv_gauss(*t, *tau_n, *tau_p, *m, *s);
}

double nmp_(double* t, double* xd){
  return nMp(*t, *xd);
}

double mp_(double* t, double* tau, double* dm){
  return Mp(*t, *tau, *dm);
}

double nmn_(double* t, double* xd){
  return nMn(*t, *xd);
}

double mn_(double* t, double* tau, double* dm){
  return Mn(*t, *tau, *dm);
}

double nmf_(double* t, double* xd){
  return nMf(*t, *xd);
}

double mf_(double* t, double* tau, double* dm){
  return Mf(*t, *tau, *dm);
}

double nap_(double* t, double* xd){
  return nAp(*t, *xd);
}

double ap_(double* t, double* tau, double* dm){
  return Ap(*t, *tau, *dm);
}

double nan_(double* t, double* xd){
  return nAn(*t, *xd);
}

double an_(double* t, double* tau, double* dm){
  return An(*t, *tau, *dm);
}

double naf_(double* t, double* xd){
  return nAf(*t, *xd);
}

double af_(double* t, double* tau, double* chi){
  return Af(*t, *tau, *chi);
}

double nmp_conv_gauss_(double* t, double* xd, double* m, double* s){
  return nMp_conv_gauss(*t, *xd, *m, *s);
}

double nmn_conv_gauss_(double* t, double* xd, double* m, double* s){
  return nMn_conv_gauss(*t, *xd, *m, *s);
}

double nap_conv_gauss_(double* t, double* xd, double* m, double* s){
  return nAp_conv_gauss(*t, *xd, *m, *s);
}

double nan_conv_gauss_(double* t, double* xd, double* m, double* s){
  return nAn_conv_gauss(*t, *xd, *m, *s);
}

double mp_conv_gauss_(double* t, double* tau, double* dm,
		      double* m, double* s){
  return Mp_conv_gauss(*t, *tau, *dm, *m, *s);
}

double mn_conv_gauss_(double* t, double* tau, double* dm,
		      double* m, double* s){
  return Mn_conv_gauss(*t, *tau, *dm, *m, *s);
}

double mf_conv_gauss_(double* t, double* tau, double* dm,
		      double* m, double* s){
  return Mf_conv_gauss(*t, *tau, *dm, *m, *s);
}

double ap_conv_gauss_(double* t, double* tau, double* dm,
		      double* m, double* s){
  return Ap_conv_gauss(*t, *tau, *dm, *m, *s);
}

double an_conv_gauss_(double* t, double* tau, double* dm,
		      double* m, double* s){
  return An_conv_gauss(*t, *tau, *dm, *m, *s);
}

double af_conv_gauss_(double* t, double* tau, double* dm,
		      double* m, double* s){
  return Af_conv_gauss(*t, *tau, *dm, *m, *s);
}

double coefees_(double* tau1, double* tau2){
  return CoefEEs(*tau1, *tau2);
}

double coefeeo_(double* tau1, double* tau2){
  return CoefEEo(*tau1, *tau2);
}

double coefeef_(double* tau1, double* tau2){
  return CoefEEf(*tau1, *tau2);
}

double sum_sigma_(double* s1, double* s2){
  return sum_sigma(*s1, *s2);
}

double sum_sigma3_(double* s1, double* s2, double* s3){
  return sum_sigma3(*s1, *s2, *s3);
}

double norm_nep_(double* ll, double* ul, double* o){
  return norm_nEp(*ll, *ul, *o);
}

double norm_nen_(double* ll, double* ul, double* o){
  return norm_nEn(*ll, *ul, *o);
}

double norm_nef_(double* ll, double* ul, double* o){
  return norm_nEf(*ll, *ul, *o);
}

double norm_ep_(double* ll, double* ul, double* tau, double* o){
  return norm_Ep(*ll, *ul, *tau, *o);
}

double norm_en_(double* ll, double* ul, double* tau, double* o){
  return norm_En(*ll, *ul, *tau, *o);
}

double norm_ef_(double* ll, double* ul, double* tau, double* o){
  return norm_Ef(*ll, *ul, *tau, *o);
}

double norm_ap_(double* ll, double* ul, double* tau, double* dm, double* o){
  return norm_Ap(*ll, *ul, *tau, *dm, *o);
}

double norm_mp_(double* ll, double* ul, double* tau, double* dm, double* o){
  return norm_Mp(*ll, *ul, *tau, *dm, *o);
}

double norm_an_(double* ll, double* ul, double* tau, double* dm, double* o){
  return norm_An(*ll, *ul, *tau, *dm, *o);
}

double norm_mn_(double* ll, double* ul, double* tau, double* dm, double* o){
  return norm_Mn(*ll, *ul, *tau, *dm, *o);
}

double norm_af_(double* ll, double* ul, double* tau, double* dm, double* o){
  return norm_Af(*ll, *ul, *tau, *dm, *o);
}

double norm_mf_(double* ll, double* ul, double* tau, double* dm, double* o){
  return norm_Mf(*ll, *ul, *tau, *dm, *o);
}

double norm_nep_conv_gauss_(double* ll, double* ul,
			    double* m, double* s, double* o){
  return norm_nEp_conv_gauss(*ll, *ul, *m, *s, *o);
}

double norm_nen_conv_gauss_(double* ll, double* ul,
			    double* m, double* s, double* o){
  return norm_nEn_conv_gauss(*ll, *ul, *m, *s, *o);
}

double norm_ep_conv_gauss_(double* ll, double* ul, double* tau
			   , double* m, double* s, double* o){
  return norm_Ep_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_en_conv_gauss_(double* ll, double* ul, double* tau,
			   double* m, double* s, double* o){
  return norm_En_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_nef_conv_gauss_(double* ll, double* ul,
			    double* m, double* s, double* o){
  return norm_nEf_conv_gauss(*ll, *ul, *m, *s, *o);
}

double norm_ef_conv_gauss_(double* ll, double* ul,
			   double* tau, double* m, double* s, double* o){
  return norm_Ef_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_enp_conv_gauss_(double* ll, double* ul,
			    double* tau_n, double* tau_p, double* m, double* s, double* o){
  return norm_Enp_conv_gauss(*ll, *ul, *tau_n, *tau_p, *m, *s, *o);
}

double norm_nan_conv_gauss_(double* ll, double* ul,
			    double* xd, double* m, double* s, double* o){
  return norm_nAn_conv_gauss(*ll, *ul, *xd, *m, *s, *o);
}

double norm_nap_conv_gauss_(double* ll, double* ul,
			    double* xd, double* m, double* s, double* o){
  return norm_nAp_conv_gauss(*ll, *ul, *xd, *m, *s, *o);
}

double norm_naf_conv_gauss_(double* ll, double* ul,
			    double* xd, double* m, double* s, double* o){
  return norm_nAf_conv_gauss(*ll, *ul, *xd, *m, *s, *o);
}

double norm_nmn_conv_gauss_(double* ll, double* ul,
			    double* xd, double* m, double* s, double* o){
  return norm_nMn_conv_gauss(*ll, *ul, *xd, *m, *s, *o);
}

double norm_nmp_conv_gauss_(double* ll, double* ul,
			    double* xd, double* m, double* s, double* o){
  return norm_nMp_conv_gauss(*ll, *ul, *xd, *m, *s, *o);
}

double norm_nmf_conv_gauss_(double* ll, double* ul,
			    double* xd, double* m, double* s, double* o){
  return norm_nMf_conv_gauss(*ll, *ul, *xd, *m, *s, *o);
}

double norm_an_conv_gauss_(double* ll, double* ul, double* tau, double* dm,
			   double* m, double* s, double* o){
  return norm_An_conv_gauss(*ll, *ul, *tau, *dm, *m, *s, *o);
}

double norm_ap_conv_gauss_(double* ll, double* ul, double* tau, double* dm,
			   double* m, double* s, double* o){
  return norm_Ap_conv_gauss(*ll, *ul, *tau, *dm, *m, *s, *o);
}

double norm_af_conv_gauss_(double* ll, double* ul, double* tau, double* dm,
			   double* m, double* s, double* o){
  return norm_Af_conv_gauss(*ll, *ul, *tau, *dm, *m, *s, *o);
}

double norm_mn_conv_gauss_(double* ll, double* ul, double* tau, double* dm,
			   double* m, double* s, double* o){
  return norm_Mn_conv_gauss(*ll, *ul, *tau, *dm, *m, *s, *o);
}

double norm_mp_conv_gauss_(double* ll, double* ul, double* tau, double* dm,
			   double* m, double* s, double* o){
  return norm_Mp_conv_gauss(*ll, *ul, *tau, *dm, *m, *s, *o);
}

double norm_mf_conv_gauss_(double* ll, double* ul, double* tau, double* dm,
			   double* m, double* s, double* o){
  return norm_Mf_conv_gauss(*ll, *ul, *tau, *dm, *m, *s, *o);
}

double ap_conv_ep_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Ap_conv_Ep(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double an_conv_en_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return An_conv_En(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mp_conv_ep_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mp_conv_Ep(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mn_conv_en_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mn_conv_En(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double ap_conv_en_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Ap_conv_En(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double an_conv_ep_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return An_conv_Ep(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mp_conv_en_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mp_conv_En(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mn_conv_ep_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mn_conv_Ep(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double ap_conv_ef_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Ap_conv_Ef(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double an_conv_ef_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return An_conv_Ef(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mp_conv_ef_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mp_conv_Ef(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mn_conv_ef_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mn_conv_Ef(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double af_conv_ep_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Af_conv_Ep(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double af_conv_en_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Af_conv_En(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mf_conv_ep_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mf_conv_Ep(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mf_conv_en_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mf_conv_En(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double af_conv_ef_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Af_conv_Ef(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double mf_conv_ef_(double* t, double* m1, double* tau1, double* dm,
		   double* m2, double* tau2){
  return Mf_conv_Ef(*t, *m1, *tau1, *dm, *m2, *tau2);
}

double ef_conv_enpgauss_(double* t, double* tau,
			 double* taunk, double* taupk,
			 double* mu, double* sigma, double* ll, double* ul){
  return Ef_conv_EnpGauss(*t, *tau, *taunk, *taupk, *mu, *sigma, *ll, *ul);
}

int ep_conv_en_coef_(double* fp, double* fn, double* tau_p, double* tau_n){
  Ep_conv_En_coef(*fp, *fn, *tau_p, *tau_n);
  return 0;
}

int en_conv_ep_coef_(double* fn, double* fp, double* tau_n, double* tau_p){
  En_conv_Ep_coef(*fn, *fp, *tau_n, *tau_p);
  return 0;
}

int ep_conv_ep_coef_(double* fp1, double* fp2,
		      double* tau_p1, double* tau_p2){
  Ep_conv_Ep_coef(*fp1, *fp2, *tau_p1, *tau_p2);
  return 0;
}

int en_conv_en_coef_(double* fn1, double* fn2,
		      double* tau_n1, double* tau_n2){
  En_conv_En_coef(*fn1, *fn2, *tau_n1, *tau_n2);
  return 0;
}

double ep_conv_en_(double* t, double* m_p, double* tau_p,
		   double* m_n, double* tau_n){
  return Ep_conv_En(*t, *m_p, *tau_p, *m_n, *tau_n);
}

double en_conv_ep_(double* t, double* m_n, double* tau_n,
		   double* m_p, double* tau_p){
  return En_conv_Ep(*t, *m_n, *tau_n, *m_p, *tau_p);
}

double xep_(double* t, double* tau){
  return xEp(*t, *tau);
}

double xen_(double* t, double* tau){
  return xEn(*t, *tau);
}

double xef_(double* t, double* tau){
  return xEf(*t, *tau);
}

double norm_xep_(double* ll, double* ul, double* tau, double* o){
  return norm_xEp(*ll, *ul, *tau, *o);
}

double norm_xen_(double* ll, double* ul, double* tau, double* o){
  return norm_xEn(*ll, *ul, *tau, *o);
}

double norm_xef_(double* ll, double* ul, double* tau, double* o){
  return norm_xEf(*ll, *ul, *tau, *o);
}

double xep_conv_gauss_(double* t, double* tau, double* m, double* s){
  return xEp_conv_gauss(*t, *tau, *m, *s);
}

double xen_conv_gauss_(double* t, double* tau, double* m, double* s){
  return xEn_conv_gauss(*t, *tau, *m, *s);
}

double xef_conv_gauss_(double* t, double* tau, double* m, double* s){
  return xEf_conv_gauss(*t, *tau, *m, *s);
}

double xep_conv_gauss_w_meanshift_(double* t, double* m1,
				   double* tau, double* m, double* s){
  return xEp_conv_gauss_w_meanshift(*t, *m1, *tau, *m, *s);
}

double xen_conv_gauss_w_meanshift_(double* t, double* m1,
				   double* tau, double* m, double* s){
  return xEn_conv_gauss_w_meanshift(*t, *m1, *tau, *m, *s);
}

double xef_conv_gauss_w_meanshift_(double* t, double* m1,
				   double* tau, double* m, double* s){
  return xEf_conv_gauss_w_meanshift(*t, *m1, *tau, *m, *s);
}

double xxep_(double* t, double* tau){
  return xxEp(*t, *tau);
}

double xxen_(double* t, double* tau){
  return xxEn(*t, *tau);
}

double xxef_(double* t, double* tau){
  return xxEf(*t, *tau);
}

double norm_xxep_(double* ll, double* ul, double* tau, double* o){
  return norm_xxEp(*ll, *ul, *tau, *o);
}

double norm_xxen_(double* ll, double* ul, double* tau, double* o){
  return norm_xxEn(*ll, *ul, *tau, *o);
}

double norm_xxef_(double* ll, double* ul, double* tau, double* o){
  return norm_xxEf(*ll, *ul, *tau, *o);
}

double xxep_conv_gauss_(double* t, double* tau, double* m, double* s){
  return xxEp_conv_gauss(*t, *tau, *m, *s);
}

double xxen_conv_gauss_(double* t, double* tau, double* m, double* s){
  return xxEn_conv_gauss(*t, *tau, *m, *s);
}

double xxef_conv_gauss_(double* t, double* tau, double* m, double* s){
  return xxEf_conv_gauss(*t, *tau, *m, *s);
}

double xxep_conv_gauss_w_meanshift_(double* t, double* m1, double* tau,
				    double* m, double* s){
  return xxEp_conv_gauss_w_meanshift(*t, *m1, *tau, *m, *s);
}

double xxen_conv_gauss_w_meanshift_(double* t, double* m1, double* tau,
				    double* m, double* s){
  return xxEn_conv_gauss_w_meanshift(*t, *m1, *tau, *m, *s);
}

double xxef_conv_gauss_w_meanshift_(double* t, double* m1, double* tau,
				    double* m, double* s){
  return xxEf_conv_gauss_w_meanshift(*t, *m1, *tau, *m, *s);
}

double ep_conv_ep_(double* t, double* m1, double* tau1,
		   double* m2, double* tau2){
  return Ep_conv_Ep(*t, *m1, *tau1, *m2, *tau2);
}

double en_conv_en_(double* t, double* m1, double* tau1,
		   double* m2, double* tau2){
  return En_conv_En(*t, *m1, *tau1, *m2, *tau2);
}

double ep_conv_ef_(double* t, double* m1, double* tau1,
		   double* m2, double* tau2){
  return Ep_conv_Ef(*t, *m1, *tau1, *m2, *tau2);
}

double en_conv_ef_(double* t, double* m1, double* tau1,
		   double* m2, double* tau2){
  return En_conv_Ef(*t, *m1, *tau1, *m2, *tau2);
}

double ef_conv_ep_(double* t, double* m1, double* tau1,
		   double* m2, double* tau2){
  return Ef_conv_Ep(*t, *m1, *tau1, *m2, *tau2);
}

double ef_conv_en_(double* t, double* m1, double* tau1,
		   double* m2, double* tau2){
  return Ef_conv_En(*t, *m1, *tau1, *m2, *tau2);
}

double ef_conv_ef_(double* t, double* m1, double* tau1,
		   double* m2, double* tau2){
  return Ef_conv_Ef(*t, *m1, *tau1, *m2, *tau2);
}

double xep_conv_ep_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEp_conv_Ep(*t, *m1, *tau1, *m2, *tau2);
}

double xen_conv_en_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEn_conv_En(*t, *m1, *tau1, *m2, *tau2);
}

double xep_conv_en_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEp_conv_En(*t, *m1, *tau1, *m2, *tau2);
}

double xen_conv_ep_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEn_conv_Ep(*t, *m1, *tau1, *m2, *tau2);
}

double xef_conv_ep_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEf_conv_Ep(*t, *m1, *tau1, *m2, *tau2);
}

double xef_conv_en_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEf_conv_En(*t, *m1, *tau1, *m2, *tau2);
}

double xef_conv_ef_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEf_conv_Ef(*t, *m1, *tau1, *m2, *tau2);
}

double xep_conv_ef_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEp_conv_Ef(*t, *m1, *tau1, *m2, *tau2);
}

double xen_conv_ef_(double* t, double* m1, double* tau1,
		    double* m2, double* tau2){
  return xEn_conv_Ef(*t, *m1, *tau1, *m2, *tau2);
}

double ep_conv_ep_ep_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mp, double* taup){
  return Ep_conv_Ep_Ep(*t, *m1, *tau1, *m2, *tau2, *mp, *taup);
}

double en_conv_en_ep_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mp, double* taup){
  return En_conv_En_Ep(*t, *m1, *tau1, *m2, *tau2, *mp, *taup);
}

double ep_conv_en_ep_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mp, double* taup){
  return Ep_conv_En_Ep(*t, *m1, *tau1, *m2, *tau2, *mp, *taup);
}

double en_conv_ep_ep_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mp, double* taup){
  return En_conv_Ep_Ep(*t, *m1, *tau1, *m2, *tau2, *mp, *taup);
}

double ef_conv_ep_ep_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mp, double* taup){
  return Ef_conv_Ep_Ep(*t, *m1, *tau1, *m2, *tau2, *mp, *taup);
}

double ef_conv_en_ep_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mp, double* taup){
  return Ef_conv_En_Ep(*t, *m1, *tau1, *m2, *tau2, *mp, *taup);
}

double ef_conv_ef_ep_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mp, double* taup){
  return Ef_conv_Ef_Ep(*t, *m1, *tau1, *m2, *tau2, *mp, *taup);
}

double ep_conv_ep_en_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mn, double* taun){
  return Ep_conv_Ep_En(*t, *m1, *tau1, *m2, *tau2, *mn, *taun);
}

double en_conv_en_en_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mn, double* taun){
  return En_conv_En_En(*t, *m1, *tau1, *m2, *tau2, *mn, *taun);
}

double ep_conv_en_en_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mn, double* taun){
  return Ep_conv_En_En(*t, *m1, *tau1, *m2, *tau2, *mn, *taun);
}

double en_conv_ep_en_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mn, double* taun){
  return En_conv_Ep_En(*t, *m1, *tau1, *m2, *tau2, *mn, *taun);
}

double ef_conv_ep_en_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mn, double* taun){
  return Ef_conv_Ep_En(*t, *m1, *tau1, *m2, *tau2, *mn, *taun);
}

double ef_conv_en_en_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mn, double* taun){
  return Ef_conv_En_En(*t, *m1, *tau1, *m2, *tau2, *mn, *taun);
}

double ef_conv_ef_en_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mn, double* taun){
  return Ef_conv_Ef_En(*t, *m1, *tau1, *m2, *tau2, *mn, *taun);
}

double ep_conv_ep_ef_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mf, double* tauf){
  return Ep_conv_Ep_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double en_conv_en_ef_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mf, double* tauf){
  return En_conv_En_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double ep_conv_en_ef_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mf, double* tauf){
  return Ep_conv_En_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double en_conv_ep_ef_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mf, double* tauf){
  return En_conv_Ep_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double ef_conv_ep_ef_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mf, double* tauf){
  return Ef_conv_Ep_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double ef_conv_en_ef_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mf, double* tauf){
  return Ef_conv_En_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double ef_conv_ef_ef_(double* t, double* m1, double* tau1,
		      double* m2, double* tau2, double* mf, double* tauf){
  return Ef_conv_Ef_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double xep_conv_ep_ef_(double* t, double* m1, double* tau1,
		       double* m2, double* tau2, double* mf, double* tauf){
  return xEp_conv_Ep_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double xen_conv_en_ef_(double* t, double* m1, double* tau1,
		       double* m2, double* tau2, double* mf, double* tauf){
  return xEn_conv_En_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double xep_conv_en_ef_(double* t, double* m1, double* tau1,
		       double* m2, double* tau2, double* mf, double* tauf){
  return xEp_conv_En_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double xen_conv_ep_ef_(double* t, double* m1, double* tau1,
		       double* m2, double* tau2, double* mf, double* tauf){
  return xEn_conv_Ep_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double xef_conv_ep_ef_(double* t, double* m1, double* tau1,
		       double* m2, double* tau2, double* mf, double* tauf){
  return xEf_conv_Ep_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double xef_conv_en_ef_(double* t, double* m1, double* tau1,
		       double* m2, double* tau2, double* mf, double* tauf){
  return xEf_conv_En_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double xef_conv_ef_ef_(double* t, double* m1, double* tau1,
		       double* m2, double* tau2, double* mf, double* tauf){
  return xEf_conv_Ef_Ef(*t, *m1, *tau1, *m2, *tau2, *mf, *tauf);
}

double ep_conv_ep_gauss_(double* t, double* m1, double* tau1,
			 double* m2, double* tau2, double* m, double* s){
  return Ep_conv_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double en_conv_en_gauss_(double* t, double* m1, double* tau1,
			 double* m2, double* tau2, double* m, double* s){
  return En_conv_En_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double ep_conv_en_gauss_(double* t, double* m1, double* tau1,
			 double* m2, double* tau2, double* m, double* s){
  return Ep_conv_En_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double en_conv_ep_gauss_(double* t, double* m1, double* tau1,
			 double* m2, double* tau2, double* m, double* s){
  return En_conv_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double ef_conv_ep_gauss_(double* t, double* m1, double* tau1,
			 double* m2, double* tau2, double* m, double* s){
  return Ef_conv_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double ef_conv_en_gauss_(double* t, double* m1, double* tau1,
			 double* m2, double* tau2, double* m, double* s){
  return Ef_conv_En_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double ef_conv_ef_gauss_(double* t, double* m1, double* tau1,
			 double* m2, double* tau2, double* m, double* s){
  return Ef_conv_Ef_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double xep_conv_ep_gauss_(double* t, double* m1, double* tau1,
			  double* m2, double* tau2, double* m, double* s){
  return xEp_conv_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double xen_conv_en_gauss_(double* t, double* m1, double* tau1,
			  double* m2, double* tau2, double* m, double* s){
  return xEn_conv_En_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double xep_conv_en_gauss_(double* t, double* m1, double* tau1,
			  double* m2, double* tau2, double* m, double* s){
  return xEp_conv_En_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double xen_conv_ep_gauss_(double* t, double* m1, double* tau1,
			  double* m2, double* tau2, double* m, double* s){
  return xEn_conv_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double xef_conv_ep_gauss_(double* t, double* m1, double* tau1,
			  double* m2, double* tau2, double* m, double* s){
  return xEf_conv_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double xef_conv_en_gauss_(double* t, double* m1, double* tau1,
			  double* m2, double* tau2, double* m, double* s){
  return xEf_conv_En_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double xef_conv_ef_gauss_(double* t, double* m1, double* tau1,
			  double* m2, double* tau2, double* m, double* s){
  return xEf_conv_Ef_gauss(*t, *m1, *tau1, *m2, *tau2, *m, *s);
}

double ep_conv_ep_ep_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mp, double* taup, double* m, double* s){
  return Ep_conv_Ep_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s);
}

double en_conv_en_ep_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mp, double* taup, double* m, double* s){
  return En_conv_En_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s);
}

double ep_conv_en_ep_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mp, double* taup, double* m, double* s){
  return Ep_conv_En_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s);
}

double en_conv_ep_ep_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2, double* mp,
			    double* taup, double* m, double* s){
  return En_conv_Ep_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s);
}

double ef_conv_ep_ep_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mp, double* taup, double* m, double* s){
  return Ef_conv_Ep_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s);
}

double ef_conv_en_ep_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mp, double* taup, double* m, double* s){
  return Ef_conv_En_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s);
}

double ef_conv_ef_ep_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mp, double* taup, double* m, double* s){
  return Ef_conv_Ef_Ep_gauss(*t, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s);
}

double ep_conv_ep_en_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mn, double* taun, double* m, double* s){
  return Ep_conv_Ep_En_gauss(*t, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s);
}

double en_conv_en_en_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mn, double* taun, double* m, double* s){
  return En_conv_En_En_gauss(*t, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s);
}

double ep_conv_en_en_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mn, double* taun, double* m, double* s){
  return Ep_conv_En_En_gauss(*t, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s);
}

double en_conv_ep_en_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mn, double* taun, double* m, double* s){
  return En_conv_Ep_En_gauss(*t, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s);
}

double ef_conv_ep_en_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mn, double* taun, double* m, double* s){
  return Ef_conv_Ep_En_gauss(*t, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s);
}

double ef_conv_en_en_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mn, double* taun, double* m, double* s){
  return Ef_conv_En_En_gauss(*t, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s);
}

double ef_conv_ef_en_gauss_(double* t, double* m1, double* tau1,
			    double* m2, double* tau2,
			    double* mn, double* taun, double* m, double* s){
  return Ef_conv_Ef_En_gauss(*t, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s);
}

double norm_xep_conv_gauss_(double* ll, double* ul, double* tau,
			    double* m, double* s, double* o){
  return norm_xEp_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_xen_conv_gauss_(double* ll, double* ul, double* tau,
			    double* m, double* s, double* o){
  return norm_xEn_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_xef_conv_gauss_(double* ll, double* ul, double* tau, 
			    double* m, double* s, double* o){
  return norm_xEf_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_xxep_conv_gauss_(double* ll, double* ul, double* tau,
			     double* m, double* s, double* o){
  return norm_xxEp_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_xxen_conv_gauss_(double* ll, double* ul, double* tau,
			     double* m, double* s, double* o){
  return norm_xxEn_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_xxef_conv_gauss_(double* ll, double* ul, double* tau,
			     double* m, double* s, double* o){
  return norm_xxEf_conv_gauss(*ll, *ul, *tau, *m, *s, *o);
}

double norm_ep_conv_ep_(double* ll, double* ul, double* m1, double* tau1,
			double* m2, double* tau2, double* o){
  return norm_Ep_conv_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_en_conv_en_(double* ll, double* ul, double* m1, double* tau1,
			double* m2, double* tau2, double* o){
  return norm_En_conv_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_ep_conv_en_(double* ll, double* ul, double* m_p, double* tau_p,
			double* m_n, double* tau_n, double* o){
  return norm_Ep_conv_En(*ll, *ul, *m_p, *tau_p, *m_n, *tau_n, *o);
}

double norm_en_conv_ep_(double* ll, double* ul, double* m1, double* tau1,
			double* m2, double* tau2, double* o){
  return norm_En_conv_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_ef_conv_ep_(double* ll, double* ul, double* m1, double* tau1,
			double* m2, double* tau2, double* o){
  return norm_Ef_conv_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_ef_conv_en_(double* ll, double* ul, double* m1, double* tau1,
			double* m2, double* tau2, double* o){
  return norm_Ef_conv_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_ef_conv_ef_(double* ll, double* ul, double* m1, double* tau1,
			double* m2, double* tau2, double* o){
  return norm_Ef_conv_Ef(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_xep_conv_ep_(double* ll, double* ul, double* m1, double* tau1,
			 double* m2, double* tau2, double* o){
  return norm_xEp_conv_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_xen_conv_en_(double* ll, double* ul, double* m1, double* tau1,
			 double* m2, double* tau2, double* o){
  return norm_xEn_conv_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_xep_conv_en_(double* ll, double* ul, double* m1, double* tau1,
			 double* m2, double* tau2, double* o){
  return norm_xEp_conv_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_xen_conv_ep_(double* ll, double* ul, double* m1, double* tau1,
			 double* m2, double* tau2, double* o){
  return norm_xEn_conv_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_xef_conv_ep_(double* ll, double* ul, double* m1, double* tau1,
			 double* m2, double* tau2, double* o){
  return norm_xEf_conv_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_xef_conv_en_(double* ll, double* ul, double* m1, double* tau1,
			 double* m2, double* tau2, double* o){
  return norm_xEf_conv_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_xef_conv_ef_(double* ll, double* ul, double* m1, double* tau1,
			 double* m2, double* tau2, double* o){
  return norm_xEf_conv_Ef(*ll, *ul, *m1, *tau1, *m2, *tau2, *o);
}

double norm_ep_conv_ep_ep_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mp, double* taup, double* o){
  return norm_Ep_conv_Ep_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *o);
}

double norm_en_conv_en_ep_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mp, double* taup, double* o){
  return norm_En_conv_En_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *o);
}

double norm_ep_conv_en_ep_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mp, double* taup, double* o){
  return norm_Ep_conv_En_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *o);
}

double norm_en_conv_ep_ep_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mp, double* taup, double* o){
  return norm_En_conv_Ep_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *o);
}

double norm_ef_conv_ep_ep_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mp, double* taup, double* o){
  return norm_Ef_conv_Ep_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *o);
}

double norm_ef_conv_en_ep_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mp, double* taup, double* o){
  return norm_Ef_conv_En_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *o);
}

double norm_ef_conv_ef_ep_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mp, double* taup, double* o){
  return norm_Ef_conv_Ef_Ep(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *o);
}

double norm_ep_conv_ep_en_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mn, double* taun, double* o){
  return norm_Ep_conv_Ep_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *o);
}

double norm_en_conv_en_en_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mn, double* taun, double* o){
  return norm_En_conv_En_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *o);
}

double norm_ep_conv_en_en_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mn, double* taun, double* o){
  return norm_Ep_conv_En_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *o);
}

double norm_en_conv_ep_en_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mn, double* taun, double* o){
  return norm_En_conv_Ep_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *o);
}

double norm_ef_conv_ep_en_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mn, double* taun, double* o){
  return norm_Ef_conv_Ep_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *o);
}

double norm_ef_conv_en_en_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mn, double* taun, double* oe){
  return norm_Ef_conv_En_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *oe);
}

double norm_ef_conv_ef_en_(double* ll, double* ul, double* m1, double* tau1,
			   double* m2, double* tau2,
			   double* mn, double* taun, double* o){
  return norm_Ef_conv_Ef_En(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *o);
}

double norm_ep_conv_ep_gauss_(double* ll, double* ul, double* m1, double* tau1,
			      double* m2, double* tau2,
			      double* m, double* s, double* o){
  return norm_Ep_conv_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_en_conv_en_gauss_(double* ll, double* ul, double* m1, double* tau1,
			      double* m2, double* tau2,
			      double* m, double* s, double* o){
  return norm_En_conv_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_ep_conv_en_gauss_(double* ll, double* ul, double* m1, double* tau1,
			      double* m2, double* tau2,
			      double* m, double* s, double* o){
  return norm_Ep_conv_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_en_conv_ep_gauss_(double* ll, double* ul, double* m1, double* tau1,
			      double* m2, double* tau2,
			      double* m, double* s, double* o){
  return norm_En_conv_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_ef_conv_ep_gauss_(double* ll, double* ul, double* m1, double* tau1,
			      double* m2, double* tau2,
			      double* m, double* s, double* o){
  return norm_Ef_conv_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_ef_conv_en_gauss_(double* ll, double* ul, double* m1, double* tau1,
			      double* m2, double* tau2,
			      double* m, double* s, double* o){
  return norm_Ef_conv_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_ef_conv_ef_gauss_(double* ll, double* ul, double* m1, double* tau1,
			      double* m2, double* tau2,
			      double* m, double* s, double* o){
  return norm_Ef_conv_Ef_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_xep_conv_ep_gauss_(double* ll, double* ul, double* m1, double* tau1,
			       double* m2, double* tau2,
			       double* m, double* s, double* o){
  return norm_xEp_conv_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_xen_conv_en_gauss_(double* ll, double* ul, double* m1, double* tau1,
			       double* m2, double* tau2,
			       double* m, double* s, double* o){
  return norm_xEn_conv_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_xep_conv_en_gauss_(double* ll, double* ul, double* m1, double* tau1,
			       double* m2, double* tau2,
			       double* m, double* s, double* o){
  return norm_xEp_conv_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_xen_conv_ep_gauss_(double* ll, double* ul, double* m1, double* tau1,
			       double* m2, double* tau2,
			       double* m, double* s, double* o){
  return norm_xEn_conv_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_xef_conv_ep_gauss_(double* ll, double* ul, double* m1, double* tau1,
			       double* m2, double* tau2,
			       double* m, double* s, double* o){
  return norm_xEf_conv_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_xef_conv_en_gauss_(double* ll, double* ul, double* m1, double* tau1,
			       double* m2, double* tau2,
			       double* m, double* s, double* o){
  return norm_xEf_conv_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_xef_conv_ef_gauss_(double* ll, double* ul, double* m1, double* tau1,
			       double* m2, double* tau2,
			       double* m, double* s, double* o){
  return norm_xEf_conv_Ef_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *m, *s, *o);
}

double norm_ep_conv_ep_ep_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mp, double* taup, double* m, double* s, double* o){
  return norm_Ep_conv_Ep_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s, *o);
}

double norm_en_conv_en_ep_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mp, double* taup, double* m, double* s, double* o){
  return norm_En_conv_En_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s, *o);
}

double norm_ep_conv_en_ep_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mp, double* taup, double* m, double* s, double* o){
  return norm_Ep_conv_En_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s, *o);
}

double norm_en_conv_ep_ep_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mp, double* taup, double* m, double* s, double* o){
  return norm_En_conv_Ep_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s, *o);
}

double norm_ef_conv_ep_ep_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mp, double* taup, double* m, double* s, double* o){
  return norm_Ef_conv_Ep_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s, *o);
}

double norm_ef_conv_en_ep_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mp, double* taup, double* m, double* s, double* o){
  return norm_Ef_conv_En_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s, *o);
}

double norm_ef_conv_ef_ep_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mp, double* taup, double* m, double* s, double* o){
  return norm_Ef_conv_Ef_Ep_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mp, *taup, *m, *s, *o);
}

double norm_ep_conv_ep_en_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mn, double* taun, double* m, double* s, double* o){
  return norm_Ep_conv_Ep_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s, *o);
}

double norm_en_conv_en_en_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mn, double* taun, double* m, double* s, double* o){
  return norm_En_conv_En_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s, *o);
}

double norm_ep_conv_en_en_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mn, double* taun, double* m, double* s, double* o){
  return norm_Ep_conv_En_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s, *o);
}

double norm_en_conv_ep_en_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mn, double* taun, double* m, double* s, double* o){
  return norm_En_conv_Ep_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s, *o);
}

double norm_ef_conv_ep_en_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mn, double* taun, double* m, double* s, double* o){
  return norm_Ef_conv_Ep_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s, *o);
}

double norm_ef_conv_en_en_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mn, double* taun, double* m, double* s, double* o){
  return norm_Ef_conv_En_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s, *o);
}

double norm_ef_conv_ef_en_gauss_(double* ll, double* ul,
				 double* m1, double* tau1, double* m2, double* tau2,
				 double* mn, double* taun, double* m, double* s, double* o){
  return norm_Ef_conv_Ef_En_gauss(*ll, *ul, *m1, *tau1, *m2, *tau2, *mn, *taun, *m, *s, *o);
}

double rk_fullrec_(double* x, double* dt, double* tau, double* ak, double* ck){
  return Rk_fullrec(*x, *dt, *tau, *ak, *ck);
}

double rk_partial_(double* x, double* dz){
  extern const dtres_param_t respar_;
  return Rk_partial(*x, *dz, &respar_);
}

double rnp_(double* x, int* flavor,
	    int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return Rnp(*x, *flavor, *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double rasc_(double* x,
	     int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return Rasc(*x, *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double rascrnp_(double* x, int* flavor,
		int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return RascRnp(*x, *flavor, *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double rdet_(double* x,
	     int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
	     int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return Rdet(*x, *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
	      *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double rdetrnp_(double* x, int* flavor,
		int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
		int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return RdetRnp(*x, *flavor, *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
		 *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double efrk_fullrec_(double* x, double* tau, double* ak, double* ck){
  return EfRk_fullrec(*x, *tau, *ak, *ck);
}

double afrk_fullrec_(double* x, double* tau, double* dm, double* ak, double* ck){
  return AfRk_fullrec(*x, *tau, *dm, *ak, *ck);
}

double mfrk_fullrec_(double* x, double* tau, double* dm, double* ak, double* ck){
  return MfRk_fullrec(*x, *tau, *dm, *ak, *ck);
}

double efrk_partial_(double* x, double* tau, double* dz){
  extern const dtres_param_t respar_;
  return EfRk_partial(*x, *tau, *dz, &respar_);
}

double afrk_partial_(double* x, double* tau, double* dm, double* dz){
  extern const dtres_param_t respar_;
  return AfRk_partial(*x, *tau, *dm, *dz, &respar_);
}

double mfrk_partial_(double* x, double* tau, double* dm, double* dz){
  extern const dtres_param_t respar_;
  return MfRk_partial(*x, *tau, *dm, *dz, &respar_);
}

double efrkrdet_fullrec_(double* x, double* tau, double* ak, double* ck,
			 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return EfRkRdet_fullrec(*x, *tau, *ak, *ck,
			  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double afrkrdet_fullrec_(double* x, double* tau, double* dm, double* ak, double* ck,
			 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return AfRkRdet_fullrec(*x, *tau, *dm, *ak, *ck,
			  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double mfrkrdet_fullrec_(double* x, double* tau, double* dm, double* ak, double* ck,
			 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return MfRkRdet_fullrec(*x, *tau, *dm, *ak, *ck,
			  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double efrkrdet_partial_(double* x, double* tau, double* dz,
			 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return EfRkRdet_partial(*x, *tau, *dz,
			  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double afrkrdet_partial_(double* x, double* tau, double* dm, double* dz,
			 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return AfRkRdet_partial(*x, *tau, *dm, *dz,
			  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double mfrkrdet_partial_(double* x, double* tau, double* dm, double* dz,
			 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return MfRkRdet_partial(*x, *tau, *dm, *dz,
			  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double efrkrdetrnp_fullrec_(double* x, int* flavor, double* tau, double* ak, double* ck,
			    int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			    int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return EfRkRdetRnp_fullrec(*x, *flavor, *tau, *ak, *ck,
			     *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			     *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double afrkrdetrnp_fullrec_(double* x, int* flavor, double* tau, double* dm, double* ak, double* ck,
			    int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			    int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return AfRkRdetRnp_fullrec(*x, *flavor, *tau, *dm, *ak, *ck,
			     *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			     *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double mfrkrdetrnp_fullrec_(double* x, int* flavor, double* tau, double* dm, double* ak, double* ck,
			    int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			    int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return MfRkRdetRnp_fullrec(*x, *flavor, *tau, *dm, *ak, *ck,
			     *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			     *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double efrkrdetrnp_partial_(double* x, int* flavor, double* tau, double* dz,
			    int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			    int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return EfRkRdetRnp_partial(*x, *flavor, *tau, *dz,
			     *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			     *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double afrkrdetrnp_partial_(double* x, int* flavor, double* tau, double* dm, double* dz,
			    int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			    int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return AfRkRdetRnp_partial(*x, *flavor, *tau, *dm, *dz,
			     *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			     *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double mfrkrdetrnp_partial_(double* x, int* flavor, double* tau, double* dm, double* dz,
			    int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			    int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return MfRkRdetRnp_partial(*x, *flavor, *tau, *dm, *dz,
			     *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			     *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_rk_partial_(double* ll, double* ul, double* dz){
  extern const dtres_param_t respar_;
  return norm_Rk_partial(*ll, *ul, *dz, &respar_);
}

double norm_rasc_(double* ll, double* ul,
		  int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_Rasc(*ll, *ul, *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_rrec_(double* ll, double* ul,
		  int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec){
  extern const dtres_param_t respar_;
  return norm_Rrec(*ll, *ul,
		   *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec, &respar_);
}

double norm_rnp_(double* ll, double* ul, int* flavor,
		 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_Rnp(*ll, *ul, *flavor, *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_rascrnp_(double* ll, double* ul, int* flavor,
		     int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_RascRnp(*ll, *ul, *flavor,
		      *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_rdet_(double* ll, double* ul,
		  int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
		  int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_Rdet(*ll, *ul, *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
		   *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_rdetrnp_(double* ll, double* ul, int* flavor,
		     int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
		     int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_RdetRnp(*ll, *ul, *flavor, *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
		      *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_efrk_fullrec_(double* ll, double* ul, double* tau, double* ak, double* ck){
  return norm_EfRk_fullrec(*ll, *ul, *tau, *ak, *ck);
}

double norm_afrk_fullrec_(double* ll, double* ul, double* tau, double* dm, double* ak, double* ck){
  return norm_AfRk_fullrec(*ll, *ul, *tau, *dm, *ak, *ck);
}

double norm_mfrk_fullrec_(double* ll, double* ul, double* tau, double* dm, double* ak, double* ck){
  return norm_MfRk_fullrec(*ll, *ul, *tau, *dm, *ak, *ck);
}

double norm_efrk_partial_(double* ll, double* ul, double* tau, double* dz){
  extern const dtres_param_t respar_;
  return norm_EfRk_partial(*ll, *ul, *tau, *dz, &respar_);
}

double norm_afrk_partial_(double* ll, double* ul, double* tau, double* dm, double* dz){
  extern const dtres_param_t respar_;
  return norm_AfRk_partial(*ll, *ul, *tau, *dm, *dz, &respar_);
}

double norm_mfrk_partial_(double* ll, double* ul, double* tau, double* dm, double* dz){
  extern const dtres_param_t respar_;
  return norm_MfRk_partial(*ll, *ul, *tau, *dm, *dz, &respar_);
}

double norm_efrkrdet_fullrec_(double* ll, double* ul, double* tau, double* ak, double* ck,
			      int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			      int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_EfRkRdet_fullrec(*ll, *ul, *tau, *ak, *ck,
			       *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			       *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_afrkrdet_fullrec_(double* ll, double* ul, double* tau, double* dm, double* ak, double* ck,
			      int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			      int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_AfRkRdet_fullrec(*ll, *ul, *tau, *dm, *ak, *ck,
			       *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			       *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_mfrkrdet_fullrec_(double* ll, double* ul, double* tau, double* dm, double* ak, double* ck,
			      int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			      int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_MfRkRdet_fullrec(*ll, *ul, *tau, *dm, *ak, *ck,
			       *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			       *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_efrkrdet_partial_(double* ll, double* ul, double* tau, double* dz,
			      int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			      int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_EfRkRdet_partial(*ll, *ul, *tau, *dz,
			       *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			       *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_afrkrdet_partial_(double* ll, double* ul, double* tau, double* dm, double* dz,
			      int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			      int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_AfRkRdet_partial(*ll, *ul, *tau, *dm, *dz,
			       *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			       *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_mfrkrdet_partial_(double* ll, double* ul, double* tau, double* dm, double* dz,
			      int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
			      int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc){
  extern const dtres_param_t respar_;
  return norm_MfRkRdet_partial(*ll, *ul, *tau, *dm, *dz,
			       *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
			       *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, &respar_);
}

double norm_efrkrdetrnp_fullrec_(double* ll, double* ul, int* flavor, double* tau, double* ak, double* ck,
				 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
				 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_EfRkRdetRnp_fullrec(*ll, *ul, *flavor, *tau, *ak, *ck,
				  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
				  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_afrkrdetrnp_fullrec_(double* ll, double* ul, int* flavor, double* tau, double* dm, double* ak, double* ck,
				 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
				 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_AfRkRdetRnp_fullrec(*ll, *ul, *flavor, *tau, *dm, *ak, *ck,
				  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
				  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_mfrkrdetrnp_fullrec_(double* ll, double* ul, int* flavor, double* tau, double* dm, double* ak, double* ck,
				 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
				 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_MfRkRdetRnp_fullrec(*ll, *ul, *flavor, *tau, *dm, *ak, *ck,
				  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
				  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_efrkrdetrnp_partial_(double* ll, double* ul, int* flavor, double* tau, double* dz,
				 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
				 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_EfRkRdetRnp_partial(*ll, *ul, *flavor, *tau, *dz,
				  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
				  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_afrkrdetrnp_partial_(double* ll, double* ul, int* flavor, double* tau, double* dm, double* dz,
				 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
				 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_AfRkRdetRnp_partial(*ll, *ul, *flavor, *tau, *dm, *dz,
				  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
				  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}

double norm_mfrkrdetrnp_partial_(double* ll, double* ul, int* flavor, double* tau, double* dm, double* dz,
				 int* ntrk_rec, double* sz_rec, double* chisq_z_rec, int* ndf_z_rec,
				 int* ntrk_asc, double* sz_asc, double* chisq_z_asc, int* ndf_z_asc, int* keeptagl){
  extern const dtres_param_t respar_;
  return norm_MfRkRdetRnp_partial(*ll, *ul, *flavor, *tau, *dm, *dz,
				  *ntrk_rec, *sz_rec, *chisq_z_rec, *ndf_z_rec,
				  *ntrk_asc, *sz_asc, *chisq_z_asc, *ndf_z_asc, *keeptagl, &respar_);
}




double add_outlier_(const int* expno,
		    const double* x,
		    const double* Lin,
		    const int* ntrk_rec, const int* ntrk_asc){
  extern const dtres_param_t respar_;
  return add_outlier(*expno, *x, *Lin, *ntrk_rec, *ntrk_asc, &respar_);
}

double add_norm_outlier_(const int* expno,
			 const double *ll, const double *ul,
			 const double *nLin,
			 const int *ntrk_rec, const int *ntrk_asc){
  extern const dtres_param_t respar_;
  return add_norm_outlier(*expno, *ll,  *ul, *nLin,
			  *ntrk_rec, *ntrk_asc, &respar_);
}

extern "C" void dt_resol_global_set_bgc(const double bgc);
extern "C" void dt_resol_global_set_bmass(const double x);
extern "C" void dt_resol_global_set_b0_mass(const double x);
extern "C" void dt_resol_global_set_bplus_mass(const double x);

int dt_resol_global_set_bgc_(const double* bgc){
  dt_resol_global_set_bgc(*bgc);
  return 0;
}
int dt_resol_global_set_bmass_(const double* x){
  dt_resol_global_set_bmass(*x);
  return 0;
}
int dt_resol_global_set_b0_mass_(const double* x){
  dt_resol_global_set_b0_mass(*x);
  return 0;
}
int dt_resol_global_set_bplus_mass_(const double* x){
  dt_resol_global_set_bplus_mass(*x);
  return 0;
}

__CLOSE_LINKAGE__
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
