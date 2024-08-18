//
//   libcnvl.cc
//     --- Description: Library for Convolurion (See belle note #307.)
//         Note that the definition of Ap, An, Mp and Mp is different 
//         from that of belle note #370 by factor tau.
//
//   $Id: libcnvl.cc 9944 2006-11-29 07:36:07Z katayama $
//
//   $Log$
//   Revision 1.2  2005/09/22 05:16:47  kohji
//   1. Fixed the out-of-range bug in test_tatami_toymc.cc
//          and norm_Mf() bug in libcnv.cc, pointed out by Sasa.
//   2. Added resolution parameters used for LP05,
//    dtres_param_svd2_lp05, dtres_param_posi_error_svd2_lp05,
//    	 dtres_param_nega_error_svd2_lp05.
//
//   Revision 1.1  2002/09/05 01:27:04  katayama
//   New package tatami from Nakadaira/Sumisawa san
//
//   Revision 1.11  2002/06/20 05:50:42  nakadair
//   Update norm_neg_sup(...), norm_nag_sup(...) and norm_nmg_sup(...)
//   to avoid the floation exception.
//
//   Revision 1.10  2002/06/19 05:06:31  nakadair
//   Bug fix in _IM(...). Change the threshold of change approximation
//   method in _I{A,M}(...), E{p,n}_conv_gauss(...).
//
//   Revision 1.9  2002/06/17 18:13:19  nakadair
//   Update _I{M,A}(...), n{M,A}{n,p}_conv_gauss(...) and
//    nE{p,n}_conv_gauss(...) to Suppress floating exceptions
//    for huge input.
//
//   Revision 1.8  2002/06/08 14:42:39  nakadair
//   Change the threshold of change approximation in nEp_conv_Gauss(...) and nEn_conv_Gauss(...).
//
//   Revision 1.7  2002/05/11 07:51:50  nakadair
//   Update for solaris and GCC3.
//
//   Revision 1.6  2002/05/10 13:21:25  nakadair
//   Change output conversion specification of double variable.
//   Change to include float.h.
//   Bug fix in EfRkRdetRnp_partial_sup(...), norm_EfRkRdetRnp_partial_sup(...)
//   and add_outlier(...).
//
//   Revision 1.5  2002/05/07 21:39:21  nakadair
//   Add convolution with resolution function and fortran interface.
//
//   Revision 1.4  2002/05/06 10:41:43  nakadair
//   Implement norm_xxE?_conv_gauss functions.
//   Bug fix in xxE{p,n}_conv_gauss, xEp_conv_En and xEn_conv_Ep.
//
//   Revision 1.3  2002/03/01 14:31:29  nakadair
//   update nAp_conv_gauss, nAn_conv_gauss, nMp_conv_gauss, nMn_conv_gauss.
//   Treatment of overflow is implemented.
//
//   Revision 1.2  2002/02/21 19:15:26  nakadair
//   Update belle-note No.
//
//   Revision 1.1.1.1  2002/02/21 18:56:07  nakadair
//   Initial import(Branch from lmfit)
//
//

const char* libcnvl_id = "$Id: libcnvl.cc 9944 2006-11-29 07:36:07Z katayama $";

#include "belle.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#ifdef __sun
#include <ieeefp.h>
#endif /* __sun */

#include <cfloat>
#include "tatami/libcnvl.h"
#include "tatami/cmlxerf.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


#ifndef M_SQRTPI
#define M_SQRTPI 1.77245385090551602729
#endif /* M_SQRTPI */

#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887242096981L
#endif /* M_SQRT2 */

#ifndef M_SQRT2PI
#define M_SQRT2PI 2.50662827463100050241
#endif /* M_SQRT2PI */

extern "C" {

typedef double (*nXi_t)(const double t, const double xd);

static double xXi_conv_gauss_by_int(nXi_t p_func,
				    const double t, const double xd,
				    const double mu, const double sigma);

double rewerf(double x, double y){
  return rewerf_(&x, &y);
}

double imwerf(double x, double y){
  return imwerf_(&x, &y);
}

double recexp(double x, double y){
  return recexp_(&x, &y);
}

double imcexp(double x, double y){
  return imcexp_(&x, &y);
}

double derfc(double x){
#ifdef USE_CERNLIB_ERFC
  return derfc_(&x); 
#else /* USE_CERNLIB_ERFC */
  return erfc(x);
#endif /* USE_CERNLIB_ERFC */
}

double DiracDelta(const double x){
  if(x==0.0){
    return FLT_MAX;
  }
  return 0.0;
}
  
double gaussian(const double x, const double m, const double s){
  if(!finite(s))
    return 0;
  if(s==0.0)
    return DiracDelta(x-m);
  double inv_s = 1.0/fabs(s);
  if(!finite(inv_s))
    return DiracDelta(x-m);
  const double inv_sqrt_2pi = 
    0.398942280401432702863218082711682654917240142822265625;
  const double dx = x-m;
  return inv_sqrt_2pi*inv_s*exp(-0.5*dx*dx*inv_s*inv_s);
}
  
double
norm_gaussian_w_cutoff(const double cutoff, const double m, const double s){
  double a_s = fabs(s);
  if(s==0.0) return 1;
  const double inv_sqrt2 = 0.707106781;
  double inv_s = 1.0/a_s;
  double x1 = (cutoff+m)*inv_sqrt2*inv_s;
  double x2 = (cutoff-m)*inv_sqrt2*inv_s;
  return 1.0-0.5*erfc(x1)-0.5*erfc(x2);
}

double norm_gaussian(const double ll, const double ul, const double m, const double s){
  double a_s = fabs(s);
  if(s==0.0) return 1;
  const double inv_sqrt2 = 0.707106781;
  double inv_s = 1.0/a_s;
  double x1 = (-ll+m)*inv_sqrt2*inv_s;
  double x2 = (ul-m)*inv_sqrt2*inv_s;
  return 1.0-0.5*erfc(x1)-0.5*erfc(x2);
}

double Ep(const double t, const double tau){
  const double at   = fabs(t);
  const double inv_atau = 1.0/fabs(tau);
  if (t<0.0){
    return 0.0;
  }
  return inv_atau*exp(-at*inv_atau);
}

double En(const double t, const double tau){
  const double at   = fabs(t);
  const double inv_atau = 1.0/fabs(tau);
  if (t>=0.0){
    return 0.0;
  }
  return inv_atau*exp(-at*inv_atau);
}

double Ef(const double t, const double tau){
  const double at   = fabs(t);
  const double inv_atau = 1.0/fabs(tau);
  return 0.5*inv_atau*exp(-at*inv_atau);
}

double Enp(const double t, const double tau_n, const double tau_p){
  const double at   = fabs(t);
  const double norm = 1.0/(fabs(tau_n)+fabs(tau_p));
  double inv_atau = (t>=0) ? 1.0/fabs(tau_p) : 1.0/fabs(tau_n) ;
  return norm*exp(-at*inv_atau);
}

/* Approximation of exp(x^2)*erfc(x) (x>>1) */
inline double approx_exp2erfc(const double x){
  //  static double approx_exp2erfc(const double x){
  const double inv_sqrt_pi
    = 0.56418958354775627928034964497783221304416656494140625;
  const double inv_x = 1.0/x;
  const double f = inv_sqrt_pi*(inv_x-0.5*inv_x*inv_x*inv_x);
  return f;
}

double nEp_conv_gauss(const double t, 
		      const double m, const double s){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  static const double Tc = DBL_MAX_10_EXP*log(10.0);
  
  if(s==0.0){
    return Ep(t-m, 1.0);
  }
  double inv_s = 1.0/fabs(s);
  double dt = -t+m;
  double ex = 0.5*s*s+dt;
  if(ex<Tc){
    double f = exp(ex)*erfc(inv_sqrt2*inv_s*(s*s+dt));
    f *= 0.5;
    return f;
  }
  const double gf = exp(-0.5*dt*dt*inv_s*inv_s);
  const double x  = inv_sqrt2*(s+dt*inv_s);
  return 0.5*approx_exp2erfc(x)*gf;
}

double nEn_conv_gauss(const double t, 
		      const double m, const double s){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  if(s==0.0){
    return En(t-m, 1.0);
  }
  static const double Tc = DBL_MAX_10_EXP*log(10.0);

  double inv_s = 1.0/fabs(s);
  double dt = t-m;
  double ex = 0.5*s*s+dt;
  if(ex<Tc){
    double f = exp(ex)*erfc(inv_sqrt2*inv_s*(s*s+dt));
    f *= 0.5;
    return f;
  }
  const double gf = exp(-0.5*dt*dt*inv_s*inv_s);
  const double x  = inv_sqrt2*(s+dt*inv_s);
  return 0.5*approx_exp2erfc(x)*gf;
}

double Ep_conv_gauss(const double t, const double tau,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double f = inv_atau*nEp_conv_gauss(nt, nm, ns);
  return f;
}

double En_conv_gauss(const double t, const double tau,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double f = inv_atau*nEn_conv_gauss(nt, nm, ns);
  return f;
}

double Ef_conv_gauss(const double t, const double tau,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  double f = nEn_conv_gauss(nt, nm, ns)+nEp_conv_gauss(nt, nm, ns);
  f *= 0.5*inv_atau;
  return f;
}

double Enp_conv_gauss(const double t,
		      const double tau_n, const double tau_p,
		      const double m, const double s){

  const double inv_atau_n = 1.0/fabs(tau_n);
  const double nt_n = t*inv_atau_n;
  const double nm_n = m*inv_atau_n;
  const double ns_n = fabs(s)*inv_atau_n;

  const double inv_atau_p = 1.0/fabs(tau_p);
  const double nt_p = t*inv_atau_p;
  const double nm_p = m*inv_atau_p;
  const double ns_p = fabs(s)*inv_atau_p;

  double f = nEn_conv_gauss(nt_n, nm_n, ns_n)+nEp_conv_gauss(nt_p, nm_p, ns_p);
  f *= 1.0/(fabs(tau_n)+fabs(tau_p));
  return f;
}

double nMp(const double t, const double xd){
  const double at   = fabs(t);
  if (t<0.0){
    return 0.0;
  }
  return exp(-at)*cos(xd*t);
}

double Mp(const double t,
	  const double tau, const double dm){
  if (t<0.0){
    return 0.0;
  }
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nMp(nt, xd);
}

double nMn(const double t, const double xd){
  if (t>=0.0){
    return 0.0;
  }
  const double at = fabs(t);
  return exp(-at)*cos(xd*t);
}

double Mn(const double t,
	  const double tau, const double dm){
  if (t>=0.0){
    return 0.0;
  }
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nMn(nt, xd);
}

double nMf(const double t, const double xd){
  const double at   = fabs(t);
  return exp(-at)*cos(xd*t);
}

double Mf(const double t,
	  const double tau, const double dm){
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nMf(nt, xd);
}

double nAp(const double t, const double xd){
  if (t<0.0){
    return 0.0;
  }
  const double at = fabs(t);
  return exp(-at)*sin(xd*t);
}

double Ap(const double t,
	  const double tau, const double dm){
  if (t<0.0){
    return 0.0;
  }
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nAp(nt, xd);
}

double nAn(const double t, const double xd){
  if (t>=0.0){
    return 0.0;
  }
  const double at = fabs(t);
  return exp(-at)*sin(xd*t);
}

double An(const double t,
	  const double tau, const double dm){
  if (t>=0.0){
    return 0.0;
  }
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nAn(nt, xd);
}

double nAf(const double t, const double xd){
  const double at   = fabs(t);
  return exp(-at)*sin(xd*t);
}

double Af(const double t,
	  const double tau, const double dm){
  const double nt = t/fabs(tau);
  const double xd = dm*fabs(tau);
  return nAf(nt, xd);
}

#ifndef IX_TREAT_NEGATIVE_GAMMA
#define IX_TREAT_NEGATIVE_GAMMA 1
#endif /* IX_TREAT_NEGATIVE_GAMMA */
  
static double _IM(const double x, const double m, const double s,
		  const double beta, const double gamma, const double b){
  const double sqrt_pi
    = 1.772453850905515881919427556567825376987457275390625;
  const double inv_sqrt_2pi = 
    0.398942280401432702863218082711682654917240142822265625;
#if defined(USE_CERNLIB_WERF)&&(IX_TREAT_NEGATIVE_GAMMA==0)
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = rewerf(inv_2sqrtbeta*b, inv_2sqrtbeta*gamma);
  f *= gaussian(x, m, s)*sqrt_pi*inv_2sqrtbeta;
  return f;
#else /* USE_CERNLIB_WERF */
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = rewerf(inv_2sqrtbeta*b, inv_2sqrtbeta*fabs(gamma));
  f *= gaussian(x, m, s)*sqrt_pi*inv_2sqrtbeta;
  if(/* gamma */ inv_2sqrtbeta*gamma >=0){
    return f;
  }
  f *= -1.0;
  const double inv_beta = 1.0/beta;
  const double dx = x-m;
  const double inv_s = 1.0/fabs(s);
  const double rex = 0.25*(gamma*gamma-b*b)*inv_beta-0.5*dx*dx*inv_s*inv_s;
  double imx = 0.5*fabs(gamma)*b*inv_beta;
  f += inv_sqrt_2pi*inv_s*sqrt(M_PI*inv_beta)*recexp(rex, imx);
  return f;
#endif
}

static double _IA(const double x, const double m, const double s,
		  const double beta, const double gamma, const double b){
  const double sqrt_pi
    = 1.772453850905515881919427556567825376987457275390625;
  const double inv_sqrt_2pi = 
    0.398942280401432702863218082711682654917240142822265625;
#if defined(USE_CERNLIB_WERF)&&(IX_TREAT_NEGATIVE_GAMMA==0)
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = imwerf(inv_2sqrtbeta*b, inv_2sqrtbeta*gamma);
  f *= gaussian(x, m, s)*sqrt_pi*inv_2sqrtbeta;
  return f;
#else /* USE_CERNLIB_WERF */
  const double inv_2sqrtbeta = 0.5/sqrt(beta);
  double f = imwerf(inv_2sqrtbeta*b, inv_2sqrtbeta*fabs(gamma));
  f *= gaussian(x, m, s)*sqrt_pi*inv_2sqrtbeta;
  if( /*  gamma */  inv_2sqrtbeta*gamma  >= 0){
    return f;
  }
  const double inv_beta = 1.0/beta;
  const double dx = x-m;
  const double inv_s = 1.0/fabs(s);
  const double rex = 0.25*(gamma*gamma-b*b)*inv_beta-0.5*dx*dx*inv_s*inv_s;
  const double imx = 0.5*fabs(gamma)*b*inv_beta;
  f += inv_sqrt_2pi*inv_s*sqrt(M_PI*inv_beta)*imcexp(rex, imx);
  return f;
#endif /* USE_CERNLIB_WERF */
}

double nMp_conv_gauss(const double t, const double xd,
		      const double m, const double s){
  if(!finite(s)){
    return 0.0;
  }
  if(s==0.0){
    return Mp(t-m, 1.0, xd);
  }
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)){
    return Mp(t-m, 1.0, xd);
  }
  const double inv_s2 = 1.0/(s*s);
  double f = _IM(t, m, s, 0.5*inv_s2, 1-(t-m)*inv_s2, xd);
  if(!finite(f)){
    return xXi_conv_gauss_by_int(nMp, t, xd, m, s);
  }
  return f;
}

double nMn_conv_gauss(const double t, const double xd,
		      const double m, const double s){
  if(!finite(s)){
    return 0.0;
  }
  if(s==0.0){
    return Mn(t-m, 1.0, xd);
  }
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)){
    return Mn(t-m, 1.0, xd);
  }
  const double inv_s2 = 1.0/(s*s);
  double f = _IM(t, m, s, 0.5*inv_s2, 1+(t-m)*inv_s2, xd);
  if(!finite(f)){
    return xXi_conv_gauss_by_int(nMn, t, xd, m, s);
  }
  return f;
}

double nAp_conv_gauss(const double t, const double xd,
		      const double m, const double s){
  if(!finite(s)){
    return 0.0;
  }
  if(s==0.0){
    return Ap(t-m, 1.0, xd);
  }
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)){
    return Ap(t-m, 1.0, xd);
  }
  const double inv_s2 = 1.0/s/s;
  double f = _IA(t, m, s, 0.5*inv_s2, 1-(t-m)*inv_s2, xd);
  if(!finite(f)){
    return xXi_conv_gauss_by_int(nAp, t, xd, m, s);
  }
  return f;
}

double nAn_conv_gauss(const double t, const double xd,
		      const double m, const double s){
  if(!finite(s)){
    return 0.0;
  }
  if(s==0.0){
    return An(t-m, 1.0, xd);
  }
  const double inv_s = 1.0/fabs(s);
  if(!finite(inv_s)){
    return An(t-m, 1.0, xd);
  }
  const double inv_s2 = 1.0/s/s;
  double f = -1.0*_IA(t, m, s, 0.5*inv_s2, 1+(t-m)*inv_s2, xd);
  if(!finite(f)){
    return xXi_conv_gauss_by_int(nAn, t, xd, m, s);
  }
  return f;
}

double Mp_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nMp_conv_gauss(nt, xd, nm,ns);
  return f;
}

double Mn_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nMn_conv_gauss(nt, xd, nm,ns);
  return f;
}

double Mf_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nMn_conv_gauss(nt, xd, nm,ns)+nMp_conv_gauss(nt, xd, nm,ns);
//    f *= 0.5;
  return f;
}

double Ap_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nAp_conv_gauss(nt, xd, nm,ns);
  return f;
}

double An_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nAn_conv_gauss(nt, xd, nm,ns);
  return f;
}

double Af_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double xd = dm*fabs(tau);
  double f = nAn_conv_gauss(nt, xd, nm, ns)+nAp_conv_gauss(nt, xd, nm, ns);
//    f *= 0.5;
  return f;
}

double sum_sigma(const double s1, const double s2){
  return sqrt(s1*s1+s2*s2);
}

double sum_sigma3(const double s1, const double s2, const double s3){
  return sqrt(s1*s1+s2*s2+s3*s3);
}

double CoefEEs(const double tau1, const double tau2){
  return tau1/(tau1-tau2);
}

double CoefEEo(const double tau1, const double tau2){
  return tau1/(tau1+tau2);
}

double CoefEEf(const double tau1, const double tau2){
  const double t12 = tau1*tau1;
  const double t22 = tau2*tau2;
  return t12/(t12-t22);
}

double norm_nEp(const double ll, const double ul, const double o){
  const double nul = (ul-o>=0.0) ? (ul-o) : 0.0;
  const double nll = (ll-o>=0.0) ? (ll-o) : 0.0;
  if(nul==nll){
    return 0;
  }
  double f = exp(-nll)-exp(-nul);
  return f;
}


double norm_nEn(const double ll, const double ul, const double o){
  const double nul = (ul-o<0.0) ? (ul-o) : 0.0;
  const double nll = (ll-o<0.0) ? (ll-o) : 0.0;
  if(nul==nll){
    return 0;
  }
  double f = exp(nul)-exp(nll);
  return f;
}

double norm_nEf(const double ll, const double ul, const double o){
  double f = norm_nEn(ll, ul, o)+norm_nEp(ll, ul, o);
  f *= 0.5;
  return f;
}

double norm_Ep(const double ll, const double ul,
		const double tau, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double no  = o*inv_atau;
  return norm_nEp(nll, nul, no);
}

double norm_En(const double ll, const double ul,
		const double tau, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double no  = o*inv_atau;
  return norm_nEn(nll, nul, no);
}

double norm_Ef(const double ll, const double ul,
		const double tau, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double no  = o*inv_atau;
  return norm_nEf(nll, nul, no);
}

//  static double norm_neg_sup(const double m, const double s){
inline static double norm_neg_sup(const double m, const double s){
  static const double Tc = DBL_MAX_10_EXP*log(10.0);

  const double as = fabs(s);
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;

  const double ex = 0.5*as*as-m;
  const double x = inv_sqrt2*(as-m/as);
  if(ex<Tc){
    double f = exp(ex);
    double g = erfc(x);
    if(g>0.0&&finite(f)){
      return f*g;
    }
  }
  double f = exp(-0.5*m*m/(s*s));
  double g = approx_exp2erfc(x);
  return f*g;
}


static double norm_nEp_conv_gauss_sub(const double ll, const double ul,
				      const double m, const double s, const double o){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double f1 = erfc(inv_sqrt2*(-ll+o+m)/s) - norm_neg_sup(ll-o-m,s);
  const double f2 = erfc(inv_sqrt2*(ul-o-m)/s)  + norm_neg_sup(ul-o-m,s);
  return 0.5*(f1+f2);
}

static double norm_nEn_conv_gauss_sub(const double ll, const double ul,
				      const double m, const double s, const double o){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double f1 = erfc(inv_sqrt2*(-ll+o+m)/s)  + norm_neg_sup(-ll+o+m,s);
  const double f2 = erfc(inv_sqrt2*(ul-o-m)/s) - norm_neg_sup(-ul+o+m,s);
  return 0.5*(f1+f2);
}

double norm_nEp_conv_gauss(const double ll, const double ul,
			   const double m, const double s, const double o){
  return 1.0-norm_nEp_conv_gauss_sub(ll, ul, m, s, o);
}

double norm_nEn_conv_gauss(const double ll, const double ul,
			   const double m, const double s, const double o){
  return 1.0-norm_nEn_conv_gauss_sub(ll, ul, m, s, o);
}

double norm_Ep_conv_gauss(const double ll, const double ul, const double tau,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  return 1.0-norm_nEp_conv_gauss_sub(nll, nul, nm, ns, no);
}

double norm_En_conv_gauss(const double ll, const double ul, const double tau,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  return 1.0-norm_nEn_conv_gauss_sub(nll, nul, nm, ns, no);
}

double norm_nEf_conv_gauss(const double ll, const double ul,
			   const double m, const double s, const double o){
  double f = norm_nEn_conv_gauss(ll, ul, m, s, o);
  f       += norm_nEp_conv_gauss(ll, ul, m, s, o);
  return 0.5*f;
}

double norm_Ef_conv_gauss(const double ll, const double ul, const double tau,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = (o==0.0) ? 0.0 : o*inv_atau;
  return 1.0-0.5*(norm_nEp_conv_gauss_sub(nll, nul, nm, ns, no)
		  +norm_nEn_conv_gauss_sub(nll, nul, nm, ns, no));
//    if(ll==-ul&&m==0.0&&o==0.0){
//      const double inv_sqrt2
//        = 0.707106781186547461715008466853760182857513427734375;
//      const double f = erfc(inv_sqrt2*(ul)/s)
//        + 0.5*(norm_neg_sup(ul,s)- norm_neg_sup(-ul,s));
//      return f;
//    }
#if 0
  double f = norm_En_conv_gauss(ll, ul, tau, m, s, o);
  f       += norm_Ep_conv_gauss(ll, ul, tau, m, s, o);
  return 0.5*f;
#endif
}


double norm_Enp_conv_gauss(const double ll, const double ul, 
			   const double tau_n, const double tau_p,
			   const double m, const double s, const double o){
  const double inv_atau_n = 1.0/fabs(tau_n);
  const double nll_n = ll*inv_atau_n;
  const double nul_n = ul*inv_atau_n;
  const double nm_n  = m*inv_atau_n;
  const double ns_n  = s*inv_atau_n;
  const double no_n  = o*inv_atau_n;

  const double inv_atau_p = 1.0/fabs(tau_p);
  const double nll_p = ll*inv_atau_p;
  const double nul_p = ul*inv_atau_p;
  const double nm_p  = m*inv_atau_p;
  const double ns_p  = s*inv_atau_p;
  const double no_p  = o*inv_atau_p;

  const double fn = fabs(tau_n)/(fabs(tau_n)+fabs(tau_p));
  const double fp = fabs(tau_p)/(fabs(tau_n)+fabs(tau_p));
  double f = fn*norm_nEn_conv_gauss_sub(nll_n, nul_n, nm_n, ns_n, no_n);
  f       += fp*norm_nEp_conv_gauss_sub(nll_p, nul_p, nm_p, ns_p, no_p);
//    double f = norm_nEn_conv_gauss_sub(nll_n, nul_n, nm_n, ns_n, no_n);
//    f       += norm_nEp_conv_gauss_sub(nll_p, nul_p, nm_p, ns_p, no_p);
//    f       /= (fabs(tau_n)+fabs(tau_p));
  return 1.0-f;
}



//  static double norm_nag_sup(const double m, const double s, const double xd){
inline double norm_nag_sup(const double m, const double s, const double xd){
  const double as = fabs(s);
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double inv_as = 1.0/as;
  
  double rex = inv_sqrt2*xd*as;
  double imx = inv_sqrt2*(as+m*inv_as);
  double f=0;
  if(imx>=0){
    double f = imwerf(rex, imx)+xd*rewerf(rex, imx);
    f *= -exp(-0.5*m*m*inv_as*inv_as);
    return f;
  }

  double rexx = -rex*rex+imx*imx-0.5*m*m*inv_as*inv_as;
  double imxx = -2.0*rex*imx;
  
  f = exp(-0.5*m*m*inv_as*inv_as)*(imwerf(-rex, -imx)+xd*rewerf(-rex, -imx));
  f += -2.0*(imcexp(rexx, imxx)+xd*recexp(rexx, imxx));
  return f;
}

inline double norm_nmg_sup(const double m, const double s, const double xd){
//  static double norm_nmg_sup(const double m, const double s, const double xd){
  const double as = fabs(s);
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double inv_as = 1.0/as;

  double rex = inv_sqrt2*xd*as;
  double imx = inv_sqrt2*(as+m*inv_as);
  double f=0;
  if(imx>=0){
    f = rewerf(rex, imx)-xd*imwerf(rex, imx);
    f *= -exp(-0.5*m*m*inv_as*inv_as);
    return f;
  }
  double rexx = -rex*rex+imx*imx-0.5*m*m*inv_as*inv_as;
  double imxx = -2.0*rex*imx;
  
  f = exp(-0.5*m*m*inv_as*inv_as)*(rewerf(-rex, -imx)-xd*imwerf(-rex, -imx));
  f += -2.0*(recexp(rexx, imxx)-xd*imcexp(rexx, imxx));
  return f;
}
  
double norm_nAn_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s, const double o){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double inv_as = 1.0/fabs(s);
  
  double f = -xd;
  f += 0.5*(xd*erfc(inv_sqrt2*inv_as*(-ll+o+m))-norm_nag_sup(ll-o-m, s, xd));
  f += 0.5*(xd*erfc(inv_sqrt2*inv_as*(ul-o-m)) +norm_nag_sup(ul-o-m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double norm_nAp_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s, const double o){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double inv_as = 1.0/fabs(s);
  
  double f = xd;
  f += -0.5*(xd*erfc(inv_sqrt2*inv_as*(-ll+o+m))+norm_nag_sup(-ll+o+m, s, xd));
  f += -0.5*(xd*erfc(inv_sqrt2*inv_as*(ul-o-m)) -norm_nag_sup(-ul+o+m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double norm_nAf_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s, const double o){
  double f = norm_nAn_conv_gauss(ll, ul, xd, m, s, o);
  f += norm_nAp_conv_gauss(ll, ul, xd, m, s, o);
  return f;
}

double norm_nMn_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s, const double o){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double inv_as = 1.0/fabs(s);
  
  double f = 1.0;
  f += -0.5*(erfc(inv_sqrt2*inv_as*(-ll+o+m))-norm_nmg_sup(ll-o-m, s, xd));
  f += -0.5*(erfc(inv_sqrt2*inv_as*(ul-o-m)) +norm_nmg_sup(ul-o-m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double norm_nMp_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s, const double o){
  const double inv_sqrt2
    = 0.707106781186547461715008466853760182857513427734375;
  const double inv_as = 1.0/fabs(s);
  
  double f = 1.0;
  f += -0.5*(erfc(inv_sqrt2*inv_as*(-ll+o+m))+norm_nmg_sup(-ll+o+m, s, xd));
  f += -0.5*(erfc(inv_sqrt2*inv_as*(ul-o-m)) -norm_nmg_sup(-ul+o+m, s, xd));
  f *= 1.0/(1.0+xd*xd);
  return f;
}

double norm_nMf_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s, const double o){
  double f = norm_nMn_conv_gauss(ll, ul, xd, m, s, o);
  f += norm_nMp_conv_gauss(ll, ul, xd, m, s, o);
  return f;
}

double norm_An_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nAn_conv_gauss(nll, nul, xd, nm, ns, no);
}

double norm_Ap_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nAp_conv_gauss(nll, nul, xd, nm, ns, no);
}

double norm_Af_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nAf_conv_gauss(nll, nul, xd, nm, ns, no);
}

double norm_Mn_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nMn_conv_gauss(nll, nul, xd, nm, ns, no);
}

double norm_Mp_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nMp_conv_gauss(nll, nul, xd, nm, ns, no);
}

double norm_Mf_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nll = ll*inv_atau;
  const double nul = ul*inv_atau;
  const double nm  = m*inv_atau;
  const double ns  = s*inv_atau;
  const double no  = o*inv_atau;
  const double xd  = dm*fabs(tau);
  return tau*norm_nMf_conv_gauss(nll, nul, xd, nm, ns, no);
}

double Ap_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  double f = A*inv_at2*(-inv_tau*Ap(nt, tau1, dm)-dm*Mp(nt, tau1, dm));
  f += A*dm*Ep(nt, tau2);
  return f;
}

double An_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = -inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  double f = -A*inv_at2*(-inv_tau*An(nt, tau1, dm)-dm*Mn(nt, tau1, dm));
  f += -A*dm*En(nt, tau2);
  return f;
}

double Mp_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  double f = A*inv_at2*(-inv_tau*Mp(nt, tau1, dm)+dm*Ap(nt, tau1, dm));
  f += A*inv_tau*Ep(nt, tau2);
  return f;
}

double Mn_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  double f = A*inv_at2*(-inv_tau*Mn(nt, tau1, dm)-dm*An(nt, tau1, dm));
  f += A*inv_tau*En(nt, tau2);
  return f;
}

double Ap_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  if ( nt > 0.0 ){
    double f = A*inv_at2*(inv_tau*Ap(nt, tau1, dm)+dm*Mp(nt, tau1, dm));
    return f;
  }
  double f = A*dm*En(nt, tau2);
  return f;
}

double An_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  if ( nt < 0.0 ){
    double f = A*inv_at2*(+inv_tau*An(nt, tau1, dm)-dm*Mn(nt, tau1, dm));
    return f;
  }
  double f = -A*dm*Ep(nt, tau2);
  return f;
}

double Mp_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  if ( nt > 0.0 ){
    double f = A*inv_at2*(inv_tau*Mp(nt, tau1, dm)-dm*Ap(nt, tau1, dm));
    return f;
  }
  double f = A*inv_tau*En(nt, tau2);
  return f;
}

double Mn_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  const double inv_at1 = 1.0/fabs(tau1);
  const double inv_at2 = 1.0/fabs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  const double nt = t-m1-m2;
  if ( nt < 0.0 ){
    double f = A*inv_at2*(inv_tau*Mn(nt, tau1, dm)+dm*An(nt, tau1, dm));
    return f;
  }
  double f = A*inv_tau*Ep(nt, tau2);
  return f;
}

double Ap_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Ap_conv_Ep(t, m1, tau1, dm, m2, tau2)
    +  Ap_conv_En(t, m1, tau1, dm, m2, tau2);
  f *= 0.5;
  return f;
}

double An_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = An_conv_Ep(t, m1, tau1, dm, m2, tau2)
    +  An_conv_En(t, m1, tau1, dm, m2, tau2);
  f *= 0.5;
  return f;
}

double Mp_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Mp_conv_Ep(t, m1, tau1, dm, m2, tau2)
    +  Mp_conv_En(t, m1, tau1, dm, m2, tau2);
  f *= 0.5;
  return f;
}

double Mn_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Mn_conv_Ep(t, m1, tau1, dm, m2, tau2)
    +  Mn_conv_En(t, m1, tau1, dm, m2, tau2);
  f *= 0.5;
  return f;
}

double Af_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Ap_conv_Ep(t, m1, tau1, dm, m2, tau2)
    +  An_conv_Ep(t, m1, tau1, dm, m2, tau2);
  return f;
}

double Af_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Ap_conv_En(t, m1, tau1, dm, m2, tau2)
    +  An_conv_En(t, m1, tau1, dm, m2, tau2);
  return f;
}

double Mf_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Mp_conv_Ep(t, m1, tau1, dm, m2, tau2)
    +  Mn_conv_Ep(t, m1, tau1, dm, m2, tau2);
  return f;
}

double Mf_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Mp_conv_En(t, m1, tau1, dm, m2, tau2)
    +  Mn_conv_En(t, m1, tau1, dm, m2, tau2);
  return f;
}

double Af_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Ap_conv_Ef(t, m1, tau1, dm, m2, tau2)
    +  An_conv_Ef(t, m1, tau1, dm, m2, tau2);
  return f;
}

double Mf_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2){
  double f = Mp_conv_Ef(t, m1, tau1, dm, m2, tau2)
    +  Mn_conv_Ef(t, m1, tau1, dm, m2, tau2);
  return f;
}

double Ef_conv_EnpGauss(const double t, const double tau,
			const double taunk, const double taupk,
			const double mu, const double sigma,
			const double ll, const double ul){
  double Lp  = Ep_conv_gauss(t, tau, mu, sigma);
  double nLp = norm_Ep_conv_gauss(ll, ul, tau, mu, sigma);
  if(nLp>0.0){
    Lp /= nLp;
  }

  double Ln  = En_conv_gauss(t, tau, mu, sigma);
  double nLn = norm_En_conv_gauss(ll, ul, tau, mu, sigma);
  if(nLn>0.0){
    Ln /= nLn;
  }

  double Lpk  = Ep_conv_gauss(t, taupk, mu, sigma);
  double nLpk = norm_Ep_conv_gauss(ll, ul, taupk, mu, sigma);
  if(nLpk>0.0){
    Lpk /= nLpk;
  }

  double Lnk  = En_conv_gauss(t, taunk, mu, sigma);
  double nLnk = norm_En_conv_gauss(ll, ul, taunk, mu, sigma);
  if(nLnk>0.0){
    Lnk /= nLnk;
  }
  
  const double inv_taupnk = 1.0/(taupk+taunk);

  const double inv_taumtaupk = 1.0/(tau-taupk);
  const double inv_tauptaupk = 1.0/(tau+taupk);

  const double inv_taumtaunk = 1.0/(tau-taunk);
  const double inv_tauptaunk = 1.0/(tau+taunk);
  
  const double fLp = 0.5*tau*inv_taupnk*(taupk*inv_taumtaupk
					 +taunk*inv_tauptaunk);

  const double fLn = 0.5*tau*inv_taupnk*(taupk*inv_tauptaupk
					 +taunk*inv_taumtaunk);

  const double fLpk = -taupk*taupk*taupk*inv_taupnk*inv_taumtaupk*inv_tauptaupk;
  const double fLnk = -taunk*taunk*taunk*inv_taupnk*inv_taumtaunk*inv_tauptaunk;

  const double Li = fLp*Lp + fLn*Ln + fLpk*Lpk + fLnk*Lnk;
  return Li;
}

double xEp(const double t, const double tau){
  const double inv_atau = 1.0/fabs(tau);
  return t*inv_atau*Ep(t, tau);
}

double xEn(const double t, const double tau){
  const double inv_atau = 1.0/fabs(tau);
  return -t*inv_atau*En(t, tau);
}

double xEf(const double t, const double tau){
  const double inv_atau = 1.0/fabs(tau);
  return fabs(t)*inv_atau*Ef(t, tau);
}


double norm_xEp(const double ll, const double ul,
		const double tau, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nul = (ul-o>=0.0) ? (ul-o)*inv_atau : 0.0;
  const double nll = (ll-o>=0.0) ? (ll-o)*inv_atau : 0.0;
  if(nul==nll){
    return 0;
  }
  const double vu = exp(-nul)*(-1.0-nul);
  const double vl = exp(-nll)*(-1.0-nll);
  return (vu-vl);
}

double norm_xEn(const double ll, const double ul,
		const double tau, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nul = (ul-o<0.0) ? (ul-o)*inv_atau : 0.0;
  const double nll = (ll-o<0.0) ? (ll-o)*inv_atau : 0.0;
  if(nul==nll){
    return 0;
  }
  const double vu = exp(nul)*(1.0-nul);
  const double vl = exp(nll)*(1.0-nll);
  return (vu-vl);
}

double norm_xEf(const double ll, const double ul,
		const double tau, const double o){
  double f = norm_xEn(ll, ul, tau, o)+norm_xEp(ll, ul, tau, o);
  f *= 0.5;
  return f;
}

double xEp_conv_gauss(const double t, const double tau,
		      const double m, const double s){
  if(s==0.0){
    return xEp(t-m, tau);
  }
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  const double f
    = (-ns2+nt-nm)*Ep_conv_gauss(t, tau, m, s) + ns2*gaussian(t, m, s);
  return f;
}

double xEn_conv_gauss(const double t, const double tau,
		      const double m, const double s){
  if(s==0.0){
    return xEn(t-m, tau);
  }
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  const double f
    = (-ns2-nt+nm)*En_conv_gauss(t, tau, m, s) + ns2*gaussian(t, m, s);
  return f;
}

double xEf_conv_gauss(const double t, const double tau,
		      const double m, const double s){
  if(s==0.0){
    return xEn(t-m, tau);
  }
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  const double f
    = 0.5*(-ns2+nt-nm)*Ep_conv_gauss(t, tau, m, s)
    + 0.5*(-ns2-nt+nm)*En_conv_gauss(t, tau, m, s)
    + ns2*gaussian(t, m, s);
  return f;
}
  
double xEp_conv_gauss_w_meanshift(const double t,
				  const double m1, const double tau,
				  const double m, const double s){
  return xEp_conv_gauss(t-m1, tau, m, s);
}

double xEn_conv_gauss_w_meanshift(const double t,
				  const double m1, const double tau,
				  const double m, const double s){
  return xEn_conv_gauss(t-m1, tau, m, s);
}

double xEf_conv_gauss_w_meanshift(const double t,
				  const double m1, const double tau,
				  const double m, const double s){
  return xEf_conv_gauss(t-m1, tau, m, s);
}

double xxEp(const double t, const double tau){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  return 0.5*nt*nt*Ep(t, tau);
}

double xxEn(const double t, const double tau){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  return 0.5*nt*nt*En(t, tau);
}

double xxEf(const double t, const double tau){
  const double inv_atau = 1.0/fabs(tau);
  const double nt = t*inv_atau;
  return 0.5*nt*nt*Ef(t, tau);
}

double norm_xxEp(const double ll, const double ul,
		const double tau, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nul = (ul-o>=0.0) ? (ul-o)*inv_atau : 0.0;
  const double nll = (ll-o>=0.0) ? (ll-o)*inv_atau : 0.0;
  if(nul==nll){
    return 0;
  }
  const double vu = -0.5*exp(-nul)*(1.0+2.0*nul+nul*nul);
  const double vl = -0.5*exp(-nll)*(1.0+2.0*nll+nll*nll);
  return (vu-vl);
}

double norm_xxEn(const double ll, const double ul,
		const double tau, const double o){
  const double inv_atau = 1.0/fabs(tau);
  const double nul = (ul-o<0.0) ? (ul-o)*inv_atau : 0.0;
  const double nll = (ll-o<0.0) ? (ll-o)*inv_atau : 0.0;
  if(nul==nll){
    return 0;
  }
  const double vu = 0.5*exp(nul)*(1.0-2.0*nul+nul*nul);
  const double vl = 0.5*exp(nll)*(1.0-2.0*nll+nll*nll);
  return (vu-vl);
}

double norm_xxEf(const double ll, const double ul,
		const double tau, const double o){
  double f = norm_xxEn(ll, ul, tau, o)+norm_xxEp(ll, ul, tau, o);
  f *= 0.5;
  return f;
}

double xxEp_conv_gauss(const double t, const double tau,
		       const double m, const double s){
  if(s==0.0){
    return xxEp(t-m, tau);
  }
  const double inv_atau = 1.0/fabs(tau);
//   const double nt = t*inv_atau;
//   const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double z = t-m-s*ns;
  const double z_s = z/s;
  const double ns2 = ns*ns;
  
  double f = 0.5*z*inv_atau*ns2*gaussian(t, m, s);
  f += 0.5*(1.0+z_s*z_s)*ns2*Ep_conv_gauss(t, tau, m, s);
  return f;
}

double xxEn_conv_gauss(const double t, const double tau,
		       const double m, const double s){
  if(s==0.0){
    return xxEn(t-m, tau);
  }
  const double inv_atau = 1.0/fabs(tau);
//   const double nt = t*inv_atau;
//   const double nm = m*inv_atau;
  const double ns = fabs(s)*inv_atau;
  const double z = -t+m-s*ns;
  const double z_s = z/s;
  const double ns2 = ns*ns;

  double f = 0.5*z*inv_atau*ns2*gaussian(t, m, s);
  f += 0.5*(1.0+z_s*z_s)*ns2*En_conv_gauss(t, tau, m, s);
  return f;
}

double xxEf_conv_gauss(const double t, const double tau,
		       const double m, const double s){
  return 0.5*(xxEp_conv_gauss(t, tau, m, s)+xxEn_conv_gauss(t, tau, m, s));
}

double xxEp_conv_gauss_w_meanshift(const double t,
				  const double m1, const double tau,
				  const double m, const double s){
  return xxEp_conv_gauss(t-m1, tau, m, s);
}

double xxEn_conv_gauss_w_meanshift(const double t,
				   const double m1, const double tau,
				   const double m, const double s){
  return xxEn_conv_gauss(t-m1, tau, m, s);
}

double xxEf_conv_gauss_w_meanshift(const double t,
				   const double m1, const double tau,
				   const double m, const double s){
  return xxEf_conv_gauss(t-m1, tau, m, s);
}
  
void Ep_conv_En_coef(double& fp, double& fn,
		     const double tau_p, const double tau_n){
  const double inv_tausum = 1.0/(tau_p+tau_n);
  fp = tau_p*inv_tausum;
  fn = tau_n*inv_tausum;
  return;
}

void En_conv_Ep_coef(double& fn, double& fp,
			    const double tau_n, const double tau_p){
  Ep_conv_En_coef(fp, fn, tau_p, tau_n);
  return;
}



void Ep_conv_Ep_coef(double& fp1, double& fp2,
		     const double tau_p1, const double tau_p2){
  const double inv_tausub = 1.0/(tau_p1-tau_p2);
  fp1 =  tau_p1*inv_tausub;
  fp2 = -tau_p2*inv_tausub;
  return;
}

void En_conv_En_coef(double& fn1, double& fn2,
		     const double tau_n1, const double tau_n2){
  const double inv_tausub = 1.0/(tau_n1-tau_n2);
  fn1 =  tau_n1*inv_tausub;
  fn2 = -tau_n2*inv_tausub;
  return;
}

double Ep_conv_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp(t-m1-m2, tau1);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return f1*Ep(t-m1-m2, tau1)+f2*Ep(t-m1-m2, tau2);
}

double En_conv_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn(t-m1-m2, tau1);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return f1*En(t-m1-m2, tau1)+f2*En(t-m1-m2, tau2);
}

double Ep_conv_En(const double t,
		  const double m_p, const double tau_p,
		  const double m_n, const double tau_n){
  double fp=0.0, fn=0.0;
  En_conv_Ep_coef(fp, fn, tau_p, tau_n);
  return fp*Ep(t-m_p-m_n, tau_p)+fn*En(t-m_p-m_n, tau_n);
}


double En_conv_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  return Ep_conv_En(t, m2, tau2, m1, tau1);
}

double Ef_conv_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  return 0.5*(En_conv_Ep(t, m1, tau1, m2, tau2)+Ep_conv_Ep(t, m1, tau1, m2, tau2));
}

double Ef_conv_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  return 0.5*(En_conv_En(t, m1, tau1, m2, tau2)
	      +Ep_conv_En(t, m1, tau1, m2, tau2));
}

double Ep_conv_Ef(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  return 0.5*(Ep_conv_En(t, m1, tau1, m2, tau2)
	      +Ep_conv_Ep(t, m1, tau1, m2, tau2));
}

double En_conv_Ef(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  return 0.5*(En_conv_En(t, m1, tau1, m2, tau2)
	      +En_conv_Ep(t, m1, tau1, m2, tau2));
}

double Ef_conv_Ef(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  return 0.5*(Ef_conv_En(t, m1, tau1, m2, tau2)
	      +Ef_conv_Ep(t, m1, tau1, m2, tau2));
}

//
double xEp_conv_Ep(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xxEp(t-m1-m2, tau1);
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*xEp(t-m1-m2, tau1);
  f += tau2*inv_tausub2*(-tau1*Ep(t-m1-m2, tau1)+tau2*Ep(t-m1-m2, tau2));
  return f;
}

double xEn_conv_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xxEn(t-m1-m2, tau1);
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*xEn(t-m1-m2, tau1);
  f += tau2*inv_tausub2*(-tau1*En(t-m1-m2, tau1)+tau2*En(t-m1-m2, tau2));
  return f;
}

double xEp_conv_En(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*xEp(t-m1-m2, tau1);
  f += tau2*inv_tausum2*(tau1*Ep(t-m1-m2, tau1)+tau2*En(t-m1-m2, tau2));
  return f;
}

double xEn_conv_Ep(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*xEn(t-m1-m2, tau1);
  f += tau2*inv_tausum2*(tau1*En(t-m1-m2, tau1)+tau2*Ep(t-m1-m2, tau2));
  return f;
}

double xEf_conv_Ep(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  return 0.5*(xEn_conv_Ep(t, m1, tau1, m2, tau2)
	      +xEp_conv_Ep(t, m1, tau1, m2, tau2));
}

double xEf_conv_En(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  return 0.5*(xEn_conv_En(t, m1, tau1, m2, tau2)
	      +xEp_conv_En(t, m1, tau1, m2, tau2));
}

double xEf_conv_Ef(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  return 0.5*(xEf_conv_En(t, m1, tau1, m2, tau2)
	      +xEf_conv_Ep(t, m1, tau1, m2, tau2));
}

double xEp_conv_Ef(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  return 0.5*(xEn_conv_En(t, m1, tau1, m2, tau2)+
	      xEn_conv_Ep(t, m1, tau1, m2, tau2));
}

double xEn_conv_Ef(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2){
  return 0.5*(xEn_conv_En(t, m1, tau1, m2, tau2)+
	      xEn_conv_Ep(t, m1, tau1, m2, tau2));
}

// 

double Ep_conv_Ep_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp_conv_Ep(t, m1+m2, tau1, mp, taup);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return f1*Ep_conv_Ep(t, m1+m2, tau1, mp, taup)+f2*Ep_conv_Ep(t, m1+m2, tau2, mp, taup);
}

double En_conv_En_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2,
		  const double mp, const double taup){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn_conv_Ep(t, m1+m2, tau1, mp, taup);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return f1*En_conv_Ep(t, m1+m2, tau1, mp, taup)+f2*En_conv_Ep(t, m1+m2, tau2, mp, taup);
}


double Ep_conv_En_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2,
		  const double mp, const double taup){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return fp*Ep_conv_Ep(t, m1+m2, tau1, mp, taup)+fn*En_conv_Ep(t, m1+m2, tau2, mp, taup);
}

double En_conv_Ep_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup){
  return Ep_conv_En_Ep(t, m2, tau2, m1, tau1, mp, taup);
}

double Ef_conv_Ep_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup){
  return 0.5*(En_conv_Ep_Ep(t, m1, tau1, m2, tau2, mp, taup)
	      +Ep_conv_Ep_Ep(t, m1, tau1, m2, tau2, mp, taup));
}

double Ef_conv_En_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup){
  return 0.5*(En_conv_En_Ep(t, m1, tau1, m2, tau2, mp, taup)
	      +Ep_conv_En_Ep(t, m1, tau1, m2, tau2, mp, taup));
}

double Ef_conv_Ef_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup){
  return 0.5*(Ef_conv_En_Ep(t, m1, tau1, m2, tau2, mp, taup)
	      +Ef_conv_Ep_Ep(t, m1, tau1, m2, tau2, mp, taup));
}

double Ep_conv_Ep_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp_conv_En(t, m1+m2, tau1, mn, taun);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return f1*Ep_conv_En(t, m1+m2, tau1, mn, taun)+f2*Ep_conv_En(t, m1+m2, tau2, mn, taun);
}

double En_conv_En_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2,
		  const double mn, const double taun){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn_conv_En(t, m1+m2, tau1, mn, taun);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return f1*En_conv_En(t, m1+m2, tau1, mn, taun)+f2*En_conv_En(t, m1+m2, tau2, mn, taun);
}


double Ep_conv_En_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2,
		  const double mn, const double taun){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return fp*Ep_conv_En(t, m1+m2, tau1, mn, taun)+fn*En_conv_En(t, m1+m2, tau2, mn, taun);
}

double En_conv_Ep_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun){
  return Ep_conv_En_En(t, m2, tau2, m1, tau1, mn, taun);
}

double Ef_conv_Ep_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun){
  return 0.5*(En_conv_Ep_En(t, m1, tau1, m2, tau2, mn, taun)
	      +Ep_conv_Ep_En(t, m1, tau1, m2, tau2, mn, taun));
}

double Ef_conv_En_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun){
  return 0.5*(En_conv_En_En(t, m1, tau1, m2, tau2, mn, taun)
	      +Ep_conv_En_En(t, m1, tau1, m2, tau2, mn, taun));
}

double Ef_conv_Ef_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun){
  return 0.5*(Ef_conv_En_En(t, m1, tau1, m2, tau2, mn, taun)
	      +Ef_conv_Ep_En(t, m1, tau1, m2, tau2, mn, taun));
}

// 
double Ep_conv_Ep_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf){
  return 0.5*(Ep_conv_Ep_Ep(t, m1, tau1, m2, tau2, mf, tauf)
	      +Ep_conv_Ep_En(t, m1, tau1, m2, tau2, mf, tauf));
}

double En_conv_En_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf){
  return 0.5*(En_conv_En_Ep(t, m1, tau1, m2, tau2, mf, tauf)
	      +En_conv_En_En(t, m1, tau1, m2, tau2, mf, tauf));
}


double Ep_conv_En_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf){
  return 0.5*(Ep_conv_En_Ep(t, m1, tau1, m2, tau2, mf, tauf)
	      +Ep_conv_En_En(t, m1, tau1, m2, tau2, mf, tauf));
}

double En_conv_Ep_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf){
  return 0.5*(En_conv_Ep_Ep(t, m1, tau1, m2, tau2, mf, tauf)
	      +En_conv_Ep_En(t, m1, tau1, m2, tau2, mf, tauf));
}

double Ef_conv_Ep_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf){
  return 0.5*(Ef_conv_Ep_Ep(t, m1, tau1, m2, tau2, mf, tauf)
	      +Ef_conv_Ep_En(t, m1, tau1, m2, tau2, mf, tauf));
}

double Ef_conv_En_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf){
  return 0.5*(Ef_conv_En_Ep(t, m1, tau1, m2, tau2, mf, tauf)
	      +Ef_conv_En_En(t, m1, tau1, m2, tau2, mf, tauf));
}

double Ef_conv_Ef_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf){
  return 0.5*(Ef_conv_Ef_Ep(t, m1, tau1, m2, tau2, mf, tauf)
	      +Ef_conv_Ef_En(t, m1, tau1, m2, tau2, mf, tauf));
}

double xEp_conv_Ep_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf){
//    if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
//      return xxEp_conv_Ef(t, m1+m2, tau1, mf, tauf);
//    }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*xEp_conv_Ef(t, m1+m2, tau1, mf, tauf);
  f += tau2*inv_tausub2*(-tau1*Ep_conv_Ef(t, m1+m2, tau1, mf, tauf)
			 +tau2*Ep_conv_Ef(t, m1+m2, tau2, mf, tauf));
  return f;
}

double xEn_conv_En_Ef(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double mf, const double tauf){
//    if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
//      return xxEn_conv_Ef(t, m1+m2, tau1, mf, tauf);
//    }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*xEn_conv_Ef(t, m1+m2, tau1, mf, tauf);
  f += tau2*inv_tausub2*(-tau1*En_conv_Ef(t, m1+m2, tau1, mf, tauf)
			 +tau2*En_conv_Ef(t, m1+m2, tau2, mf, tauf));
  return f;
}

double xEp_conv_En_Ef(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double mf, const double tauf){
//    if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
//      return xxEp_conv_Ef(t, m1+m2, tau1, mf, tauf);
//    }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*xEp_conv_Ef(t, m1+m2, tau1, mf, tauf);
  f += tau2*inv_tausum2*(tau1*Ep_conv_Ef(t, m1+m2, tau1, mf, tauf)
			 +tau2*En_conv_Ef(t, m1+m2, tau2, mf, tauf));
  return f;
}

double xEn_conv_Ep_Ef(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double mf, const double tauf){
//    if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
//      return xxEp_conv_Ef(t, m1+m2, tau1, mf, tauf);
//    }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*xEn_conv_Ef(t, m1+m2, tau1, mf, tauf);
  f += tau2*inv_tausum2*(tau1*En_conv_Ef(t, m1+m2, tau1, mf, tauf)
			 +tau2*Ep_conv_Ef(t, m1+m2, tau2, mf, tauf));
  return f;
}

double xEf_conv_Ep_Ef(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double mf, const double tauf){
  return 0.5*(xEn_conv_Ep_Ef(t, m1, tau1, m2, tau2, mf, tauf)
	      +xEp_conv_Ep_Ef(t, m1, tau1, m2, tau2, mf, tauf));
}

double xEf_conv_En_Ef(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double mf, const double tauf){
  return 0.5*(xEn_conv_En_Ef(t, m1, tau1, m2, tau2, mf, tauf)
	      +xEp_conv_En_Ef(t, m1, tau1, m2, tau2, mf, tauf));
}

double xEf_conv_Ef_Ef(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double mf, const double tauf){
  return 0.5*(xEf_conv_En_Ef(t, m1, tau1, m2, tau2, mf, tauf)
	      +xEf_conv_Ep_Ef(t, m1, tau1, m2, tau2, mf, tauf));
}


//////////////////////////////////////////////////////////////////

double Ep_conv_Ep_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp_conv_gauss(t-m1-m2, tau1, m, s);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return f1*Ep_conv_gauss(t-m1-m2, tau1, m, s)+f2*Ep_conv_gauss(t-m1-m2, tau2, m, s);
}

double En_conv_En_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn_conv_gauss(t-m1-m2, tau1, m, s);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return f1*En_conv_gauss(t-m1-m2, tau1, m, s)+f2*En_conv_gauss(t-m1-m2, tau2, m, s);
}


double Ep_conv_En_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return fp*Ep_conv_gauss(t-m1-m2, tau1, m, s)+fn*En_conv_gauss(t-m1-m2, tau2, m, s);
}

double En_conv_Ep_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s){
  return Ep_conv_En_gauss(t, m2, tau2, m1, tau1, m, s);
}

double Ef_conv_Ep_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s){
  return 0.5*(En_conv_Ep_gauss(t, m1, tau1, m2, tau2, m, s)
	      +Ep_conv_Ep_gauss(t, m1, tau1, m2, tau2, m, s));
}

double Ef_conv_En_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s){
  return 0.5*(En_conv_En_gauss(t, m1, tau1, m2, tau2, m, s)
	      +Ep_conv_En_gauss(t, m1, tau1, m2, tau2, m, s));
}

double Ef_conv_Ef_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s){
  return 0.5*(Ef_conv_En_gauss(t, m1, tau1, m2, tau2, m, s)
	      +Ef_conv_Ep_gauss(t, m1, tau1, m2, tau2, m, s));
}
//
double xEp_conv_Ep_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xxEp_conv_gauss(t-m1-m2, tau1, m, s);
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*xEp_conv_gauss(t-m1-m2, tau1, m, s);
  f += tau2*inv_tausub2*(-tau1*Ep_conv_gauss(t-m1-m2, tau1, m, s)
			 +tau2*Ep_conv_gauss(t-m1-m2, tau2, m, s));
  return f;
}

double xEn_conv_En_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xxEn_conv_gauss(t-m1-m2, tau1, m, s);
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*xEn_conv_gauss(t-m1-m2, tau1, m, s);
  f += tau2*inv_tausub2*(-tau1*En_conv_gauss(t-m1-m2, tau1, m, s)
			 +tau2*En_conv_gauss(t-m1-m2, tau2, m, s));
  return f;
}

double xEp_conv_En_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xxEp_conv_gauss(t-m1-m2, tau1, m, s);
  }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*xEp_conv_gauss(t-m1-m2, tau1, m, s);
  f += tau2*inv_tausum2*(tau1*Ep_conv_gauss(t-m1-m2, tau1, m, s)
			 +tau2*En_conv_gauss(t-m1-m2, tau2, m, s));
  return f;
}

double xEn_conv_Ep_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xxEp_conv_gauss(t-m1-m2, tau1, m, s);
  }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*xEn_conv_gauss(t-m1-m2, tau1, m, s);
  f += tau2*inv_tausum2*(tau1*En_conv_gauss(t-m1-m2, tau1, m, s)
			 +tau2*Ep_conv_gauss(t-m1-m2, tau2, m, s));
  return f;
}

double xEf_conv_Ep_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s){
  return 0.5*(xEn_conv_Ep_gauss(t, m1, tau1, m2, tau2, m, s)
	      +xEp_conv_Ep_gauss(t, m1, tau1, m2, tau2, m, s));
}

double xEf_conv_En_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s){
  return 0.5*(xEn_conv_En_gauss(t, m1, tau1, m2, tau2, m, s)
	      +xEp_conv_En_gauss(t, m1, tau1, m2, tau2, m, s));
}

double xEf_conv_Ef_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s){
  return 0.5*(xEf_conv_En_gauss(t, m1, tau1, m2, tau2, m, s)
	      +xEf_conv_Ep_gauss(t, m1, tau1, m2, tau2, m, s));
}

// 

double Ep_conv_Ep_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp_conv_Ep_gauss(t, m1+m2, tau1, mp, taup, m, s);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*Ep_conv_Ep_gauss(t, m1+m2, tau1, mp, taup, m, s)
	  +f2*Ep_conv_Ep_gauss(t, m1+m2, tau2, mp, taup, m, s));
}

double En_conv_En_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn_conv_Ep_gauss(t, m1+m2, tau1, mp, taup, m, s);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*En_conv_Ep_gauss(t, m1+m2, tau1, mp, taup, m, s)
	  +f2*En_conv_Ep_gauss(t, m1+m2, tau2, mp, taup, m, s));
}


double Ep_conv_En_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*Ep_conv_Ep_gauss(t, m1+m2, tau1, mp, taup, m, s)
	  +fn*En_conv_Ep_gauss(t, m1+m2, tau2, mp, taup, m, s));
}

double En_conv_Ep_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s){
  return Ep_conv_En_Ep_gauss(t, m2, tau2, m1, tau1, mp, taup, m, s);
}

double Ef_conv_Ep_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s){
  return 0.5*(En_conv_Ep_Ep_gauss(t, m1, tau1, m2, tau2, mp, taup, m, s)
	      +Ep_conv_Ep_Ep_gauss(t, m1, tau1, m2, tau2, mp, taup, m, s));
}

double Ef_conv_En_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s){
  return 0.5*(En_conv_En_Ep_gauss(t, m1, tau1, m2, tau2, mp, taup, m, s)
	      +Ep_conv_En_Ep_gauss(t, m1, tau1, m2, tau2, mp, taup, m, s));
}

double Ef_conv_Ef_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s){
  return 0.5*(Ef_conv_En_Ep_gauss(t, m1, tau1, m2, tau2, mp, taup, m, s)
	      +Ef_conv_Ep_Ep_gauss(t, m1, tau1, m2, tau2, mp, taup, m, s));
}

double Ep_conv_Ep_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp_conv_En_gauss(t, m1+m2, tau1, mn, taun, m, s);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*Ep_conv_En_gauss(t, m1+m2, tau1, mn, taun, m, s)
	  +f2*Ep_conv_En_gauss(t, m1+m2, tau2, mn, taun, m, s));
}

double En_conv_En_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn_conv_En_gauss(t, m1+m2, tau1, mn, taun, m, s);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*En_conv_En_gauss(t, m1+m2, tau1, mn, taun, m, s)
	  +f2*En_conv_En_gauss(t, m1+m2, tau2, mn, taun, m, s));
}


double Ep_conv_En_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*Ep_conv_En_gauss(t, m1+m2, tau1, mn, taun, m, s)
	  +fn*En_conv_En_gauss(t, m1+m2, tau2, mn, taun, m, s));
}

double En_conv_Ep_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s){
  return Ep_conv_En_En_gauss(t, m2, tau2, m1, tau1, mn, taun, m, s);
}

double Ef_conv_Ep_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s){
  return 0.5*(En_conv_Ep_En_gauss(t, m1, tau1, m2, tau2, mn, taun, m, s)
	      +Ep_conv_Ep_En_gauss(t, m1, tau1, m2, tau2, mn, taun, m, s));
}

double Ef_conv_En_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s){
  return 0.5*(En_conv_En_En_gauss(t, m1, tau1, m2, tau2, mn, taun, m, s)
	      +Ep_conv_En_En_gauss(t, m1, tau1, m2, tau2, mn, taun, m, s));
}

double Ef_conv_Ef_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s){
  return 0.5*(Ef_conv_En_En_gauss(t, m1, tau1, m2, tau2, mn, taun, m, s)
	      +Ef_conv_Ep_En_gauss(t, m1, tau1, m2, tau2, mn, taun, m, s));
}

// 

double Ep_conv_Ep_Ep_Ef(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double mp, const double taup,
			const double mf, const double tauf){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp_conv_Ep_Ef(t, m1+m2, tau1, mp, taup, mf, tauf);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*Ep_conv_Ep_Ef(t, m1+m2, tau1, mp, taup, mf, tauf)
	  +f2*Ep_conv_Ep_Ef(t, m1+m2, tau2, mp, taup, mf, tauf));
}

double En_conv_En_Ep_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double mf, const double tauf){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn_conv_Ep_Ef(t, m1+m2, tau1, mp, taup, mf, tauf);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*En_conv_Ep_Ef(t, m1+m2, tau1, mp, taup, mf, tauf)
	  +f2*En_conv_Ep_Ef(t, m1+m2, tau2, mp, taup, mf, tauf));
}


double Ep_conv_En_Ep_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double mf, const double tauf){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*Ep_conv_Ep_Ef(t, m1+m2, tau1, mp, taup, mf, tauf)
	  +fn*En_conv_Ep_Ef(t, m1+m2, tau2, mp, taup, mf, tauf));
}

double En_conv_Ep_Ep_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double mf, const double tauf){
  return Ep_conv_En_Ep_Ef(t, m2, tau2, m1, tau1, mp, taup, mf, tauf);
}

double Ef_conv_Ep_Ep_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double mf, const double tauf){
  return 0.5*(En_conv_Ep_Ep_Ef(t, m1, tau1, m2, tau2, mp, taup, mf, tauf)
	      +Ep_conv_Ep_Ep_Ef(t, m1, tau1, m2, tau2, mp, taup, mf, tauf));
}

double Ef_conv_En_Ep_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double mf, const double tauf){
  return 0.5*(En_conv_En_Ep_Ef(t, m1, tau1, m2, tau2, mp, taup, mf, tauf)
	      +Ep_conv_En_Ep_Ef(t, m1, tau1, m2, tau2, mp, taup, mf, tauf));
}

double Ef_conv_Ef_Ep_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double mf, const double tauf){
  return 0.5*(Ef_conv_En_Ep_Ef(t, m1, tau1, m2, tau2, mp, taup, mf, tauf)
	      +Ef_conv_Ep_Ep_Ef(t, m1, tau1, m2, tau2, mp, taup, mf, tauf));
}

double Ep_conv_Ep_En_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double mf, const double tauf){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEp_conv_En_Ef(t, m1+m2, tau1, mn, taun, mf, tauf);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*Ep_conv_En_Ef(t, m1+m2, tau1, mn, taun, mf, tauf)
	  +f2*Ep_conv_En_Ef(t, m1+m2, tau2, mn, taun, mf, tauf));
}

double En_conv_En_En_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double mf, const double tauf){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return xEn_conv_En_Ef(t, m1+m2, tau1, mn, taun, mf, tauf);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*En_conv_En_Ef(t, m1+m2, tau1, mn, taun, mf, tauf)
	  +f2*En_conv_En_Ef(t, m1+m2, tau2, mn, taun, mf, tauf));
}


double Ep_conv_En_En_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double mf, const double tauf){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*Ep_conv_En_Ef(t, m1+m2, tau1, mn, taun, mf, tauf)
	  +fn*En_conv_En_Ef(t, m1+m2, tau2, mn, taun, mf, tauf));
}

double En_conv_Ep_En_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double mf, const double tauf){
  return Ep_conv_En_En_Ef(t, m2, tau2, m1, tau1, mn, taun, mf, tauf);
}

double Ef_conv_Ep_En_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double mf, const double tauf){
  return 0.5*(En_conv_Ep_En_Ef(t, m1, tau1, m2, tau2, mn, taun, mf, tauf)
	      +Ep_conv_Ep_En_Ef(t, m1, tau1, m2, tau2, mn, taun, mf, tauf));
}

double Ef_conv_En_En_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double mf, const double tauf){
  return 0.5*(En_conv_En_En_Ef(t, m1, tau1, m2, tau2, mn, taun, mf, tauf)
	      +Ep_conv_En_En_Ef(t, m1, tau1, m2, tau2, mn, taun, mf, tauf));
}

double Ef_conv_Ef_En_Ef(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double mf, const double tauf){
  return 0.5*(Ef_conv_En_En_Ef(t, m1, tau1, m2, tau2, mn, taun, mf, tauf)
	      +Ef_conv_Ep_En_Ef(t, m1, tau1, m2, tau2, mn, taun, mf, tauf));
}

// ---------------------------------------------------------------------------
double int_polyexp2(const double ll, const double ul,
		    const double alpha, const double beta, const double gamma,
		    const double a){

  const double sqrt_pi
    = 1.772453850905515881919427556567825376987457275390625;
  const double inv_2a = 0.5/a;
  const double sqrt_a = sqrt(a);

  double f = -inv_2a*(alpha*ul+beta)*exp(-a*ul*ul);
  f       +=  inv_2a*(alpha*ll+beta)*exp(-a*ll*ll);
  f += (alpha*inv_2a+gamma)*sqrt_pi/sqrt_a*(1.0-0.5*erfc(sqrt_a*ul)-0.5*erfc(sqrt_a*ll));
  return f;
}

double int_polyexp_erfc(const double ll, const double ul,
			const double alpha, const double beta, const double gamma,
			const double a){
  const double inv_a = 1.0/a;


  double f = 0.0;

  const double a1 = -2.0*alpha*inv_a + beta;
  const double a0 = 2.0*alpha*alpha*inv_a*inv_a + beta*inv_a + a*gamma;

  f +=  inv_a*exp(a*ul)*(alpha*ul*ul+a1*ul+a0)*erfc(ul);
  f += -inv_a*exp(a*ll)*(alpha*ll*ll+a1*ll+a0)*erfc(ll);
  
  const double aa1 = a1 + a*alpha;
  const double aa0 = a0 + 0.5*a*a1 + 0.25*a*a*alpha;

  f += exp(0.25*a*a)*int_polyexp2(ll-0.5*a, ul-0.5*a,
				  alpha, aa1, aa0, 1.0);
  return f;
}

double norm_xEp_conv_gauss(const double ll, const double ul, const double tau,
			   const double m, const double s, const double o){
  const double sqrt2
    = 1.4142135623730951454746218587388284504413604736328125;
  const double inv_sqrt2
    = 0.7071067811865474617150084668537601828575134277343750;

  if(s==0.0){
    return norm_xEp(ll-m, ul-m, tau, o);
  }
  const double inv_atau = 1.0/fabs(tau);
  const double inv_s = 1.0/fabs(s);

  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  const double nll = (-(ll-m-o)-ns*s)*inv_sqrt2*inv_s;
  const double nul = (-(ul-m-o)-ns*s)*inv_sqrt2*inv_s;
  
  const double f
    = ns2*exp(1.5*ns2)*int_polyexp_erfc(nll, nul, 0.0, 1.0, 0.0, ns*sqrt2)
    + ns2*norm_gaussian(ll, ul, m+o, s);
  return f;
}

double norm_xEn_conv_gauss(const double ll, const double ul, const double tau,
			   const double m, const double s, const double o){
  const double sqrt2
    = 1.4142135623730951454746218587388284504413604736328125;
  const double inv_sqrt2
    = 0.7071067811865474617150084668537601828575134277343750;

  if(s==0.0){
    return norm_xEn(ll-m, ul-m, tau, o);
  }
  const double inv_atau = 1.0/fabs(tau);
  const double inv_s = 1.0/fabs(s);

  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;

  const double nll = ((ll-m-o)-ns*s)*inv_sqrt2*inv_s;
  const double nul = ((ul-m-o)-ns*s)*inv_sqrt2*inv_s;
  
  const double f
    = ns2*exp(1.5*ns2)*int_polyexp_erfc(nll, nul, 0.0, 1.0, 0.0, ns*sqrt2)
    + ns2*norm_gaussian(ll, ul, m+o, s);
  return f;
}

double norm_xEf_conv_gauss(const double ll, const double ul, const double tau,
			   const double m, const double s, const double o){
  return 0.5*(norm_xEp_conv_gauss(ll, ul, tau, m, s, o)
	      +norm_xEn_conv_gauss(ll, ul, tau, m, s, o));
}

double norm_xxEp_conv_gauss(const double ll, const double ul, const double tau,
			    const double m, const double s, const double o){
  if(s==0.0){
    return norm_xxEp(ll-m, ul-m, tau, o);
  }
  const double inv_2sqrt2
    = 0.1994711402007163514316090413558413274586200714111328125;
  const double sqrt2
    = 1.4142135623730951454746218587388284504413604736328125;
  const double inv_sqrt2
    = 0.7071067811865474617150084668537601828575134277343750;

  const double inv_atau = 1.0/fabs(tau);
  const double inv_s = 1.0/fabs(s);

  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;
  const double ns3 = ns*ns2;

  double f = inv_2sqrt2*ns2*int_polyexp2(ll+m+o, ul+m+o,
					 0.0, inv_atau, -ns2, -0.5/(s*s));
  
  const double nll = (-(ll-m-o)-ns*s)*inv_sqrt2*inv_s;
  const double nul = (-(ul-m-o)-ns*s)*inv_sqrt2*inv_s;

  f += (inv_sqrt2*ns3*exp(1.5*ns2)
	*int_polyexp_erfc(nll, nul, 1.0, 0.0, 0.5, ns*sqrt2));

  return f;
}

double norm_xxEn_conv_gauss(const double ll, const double ul, const double tau,
			    const double m, const double s, const double o){
  if(s==0.0){
    return norm_xxEn(ll-m, ul-m, tau, o);
  }
  const double inv_2sqrt2
    = 0.1994711402007163514316090413558413274586200714111328125;
  const double sqrt2
    = 1.4142135623730951454746218587388284504413604736328125;
  const double inv_sqrt2
    = 0.7071067811865474617150084668537601828575134277343750;

  const double inv_atau = 1.0/fabs(tau);
  const double inv_s = 1.0/fabs(s);

  const double ns = fabs(s)*inv_atau;
  const double ns2 = ns*ns;
  const double ns3 = ns*ns2;

  double f = inv_2sqrt2*ns2*int_polyexp2(ll+m+o, ul+m+o,
					 0.0, -inv_atau, -ns2, -0.5/(s*s));
  
  const double nll = (-(ll-m-o)-ns*s)*inv_sqrt2*inv_s;
  const double nul = (-(ul-m-o)-ns*s)*inv_sqrt2*inv_s;

  f += (inv_sqrt2*ns3*exp(1.5*ns2)
	*int_polyexp_erfc(nll, nul, 1.0, 0.0, 0.5, ns*sqrt2));


  return f;
}

double norm_xxEf_conv_gauss(const double ll, const double ul, const double tau,
			    const double m, const double s, const double o){
  return 0.5*(norm_xxEp_conv_gauss(ll, ul, tau, m, s, o)
	      +norm_xxEn_conv_gauss(ll, ul, tau, m, s, o));
}


double norm_Ep_conv_Ep(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEp(ll-m1-m2, ul-m1-m2, tau1, o);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*norm_Ep(ll-m1-m2, ul-m1-m2, tau1, o)
	  +f2*norm_Ep(ll-m1-m2, ul-m1-m2, tau2, o));
}

double norm_En_conv_En(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEn(ll-m1-m2, ul-m1-m2, tau1, o);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*norm_En(ll-m1-m2, ul-m1-m2, tau1, o)
	  +f2*norm_En(ll-m1-m2, ul-m1-m2, tau2, o));
}

double norm_Ep_conv_En(const double ll, const double ul,
		       const double m_p, const double tau_p,
		       const double m_n, const double tau_n, const double o){
  double fp=0.0, fn=0.0;
  En_conv_Ep_coef(fp, fn, tau_p, tau_n);
  return (fp*norm_Ep(ll-m_p-m_n, ul-m_p-m_n, tau_p, o)
	  +fn*norm_En(ll-m_p-m_n, ul-m_p-m_n, tau_n, o));
}


double norm_En_conv_Ep(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2, const double o){
  return norm_Ep_conv_En(ll, ul, m2, tau2, m1, tau1, o);
}

double norm_Ef_conv_Ep(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2, const double o){
  return 0.5*(norm_En_conv_Ep(ll, ul, m1, tau1, m2, tau2, o)
	      +norm_Ep_conv_Ep(ll, ul, m1, tau1, m2, tau2, o));
}

double norm_Ef_conv_En(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2, const double o){
  return 0.5*(norm_En_conv_En(ll, ul, m1, tau1, m2, tau2, o)
	      +norm_Ep_conv_En(ll, ul, m1, tau1, m2, tau2, o));
}

double norm_Ef_conv_Ef(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2, const double o){
  return 0.5*(norm_Ef_conv_En(ll, ul, m1, tau1, m2, tau2, o)
	      +norm_Ef_conv_Ep(ll, ul, m1, tau1, m2, tau2, o));
}
//
double norm_xEp_conv_Ep(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xxEp(ll-m1-m2, ul-m1-m2, tau1, o);
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*norm_xEp(ll-m1-m2, ul-m1-m2, tau1, o);
  f += tau2*inv_tausub2*(-tau1*norm_Ep(ll-m1-m2, ul-m1-m2, tau1, o)
			 +tau2*norm_Ep(ll-m1-m2, ul-m1-m2, tau2, o));
  return f;
}

double norm_xEn_conv_En(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xxEn(ll-m1-m2, ul-m1-m2, tau1);
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*norm_xEn(ll-m1-m2, ul-m1-m2, tau1, o);
  f += tau2*inv_tausub2*(-tau1*norm_En(ll-m1-m2, ul-m1-m2, tau1, o)
			 +tau2*norm_En(ll-m1-m2, ul-m1-m2, tau2, o));
  return f;
}

double norm_xEp_conv_En(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xxEp(ll-m1-m2, ul-m1-m2, tau1, o);
  }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*norm_xEp(ll-m1-m2, ul-m1-m2, tau1, o);
  f += tau2*inv_tausum2*(tau1*norm_Ep(ll-m1-m2, ul-m1-m2, tau1, o)
			 +tau2*norm_En(ll-m1-m2, ul-m1-m2, tau2, o));
  return f;
}

double norm_xEn_conv_Ep(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xxEp(ll-m1-m2, ul-m1-m2, tau1, o);
  }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*norm_xEn(ll-m1-m2, ul-m1-m2, tau1, o);
  f += tau2*inv_tausum2*(tau1*norm_En(ll-m1-m2, ul-m1-m2, tau1, o)
			 +tau2*norm_Ep(ll-m1-m2, ul-m1-m2, tau2, o));
  return f;
}

double norm_xEf_conv_Ep(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2, const double o){
  return 0.5*(norm_xEn_conv_Ep(ll, ul, m1, tau1, m2, tau2, o)+
	      norm_xEp_conv_Ep(ll, ul, m1, tau1, m2, tau2, o));
}

double norm_xEf_conv_En(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2, const double o){
  return 0.5*(norm_xEn_conv_En(ll, ul, m1, tau1, m2, tau2, o)
	      +norm_xEp_conv_En(ll, ul, m1, tau1, m2, tau2, o));
}

double norm_xEf_conv_Ef(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2, const double o){
  return 0.5*(norm_xEf_conv_En(ll, ul, m1, tau1, m2, tau2, o)
	      +norm_xEf_conv_Ep(ll, ul, m1, tau1, m2, tau2, o));
}

// 

double norm_Ep_conv_Ep_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEp_conv_Ep(ll, ul, m1+m2, tau1, mp, taup, o);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*norm_Ep_conv_Ep(ll, ul, m1+m2, tau1, mp, taup, o)
	  +f2*norm_Ep_conv_Ep(ll, ul, m1+m2, tau2, mp, taup, o));
}

double norm_En_conv_En_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEn_conv_Ep(ll, ul, m1+m2, tau1, mp, taup, o);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*norm_En_conv_Ep(ll, ul, m1+m2, tau1, mp, taup, o)
	  +f2*norm_En_conv_Ep(ll, ul, m1+m2, tau2, mp, taup, o));
}


double norm_Ep_conv_En_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup, const double o){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*norm_Ep_conv_Ep(ll, ul, m1+m2, tau1, mp, taup, o)
	  +fn*norm_En_conv_Ep(ll, ul, m1+m2, tau2, mp, taup, o));
}

double norm_En_conv_Ep_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup, const double o){
  return norm_Ep_conv_En_Ep(ll, ul, m2, tau2, m1, tau1, mp, taup, o);
}

double norm_Ef_conv_Ep_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup, const double o){
  return 0.5*(norm_En_conv_Ep_Ep(ll, ul, m1, tau1, m2, tau2, mp, taup, o)
	      +norm_Ep_conv_Ep_Ep(ll, ul, m1, tau1, m2, tau2, mp, taup, o));
}

double norm_Ef_conv_En_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup, const double o){
  return 0.5*(norm_En_conv_En_Ep(ll, ul, m1, tau1, m2, tau2, mp, taup, o)
	      +norm_Ep_conv_En_Ep(ll, ul, m1, tau1, m2, tau2, mp, taup, o));
}

double norm_Ef_conv_Ef_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup, const double o){
  return 0.5*(norm_Ef_conv_En_Ep(ll, ul, m1, tau1, m2, tau2, mp, taup, o)
	      +norm_Ef_conv_Ep_Ep(ll, ul, m1, tau1, m2, tau2, mp, taup, o));
}

double norm_Ep_conv_Ep_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEp_conv_En(ll, ul, m1+m2, tau1, mn, taun, o);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*norm_Ep_conv_En(ll, ul, m1+m2, tau1, mn, taun, o)
	  +f2*norm_Ep_conv_En(ll, ul, m1+m2, tau2, mn, taun, o));
}

double norm_En_conv_En_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEn_conv_En(ll, ul, m1+m2, tau1, mn, taun, o);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*norm_En_conv_En(ll, ul, m1+m2, tau1, mn, taun, o)
	  +f2*norm_En_conv_En(ll, ul, m1+m2, tau2, mn, taun, o));
}


double norm_Ep_conv_En_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double o){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*norm_Ep_conv_En(ll, ul, m1+m2, tau1, mn, taun, o)
	  +fn*norm_En_conv_En(ll, ul, m1+m2, tau2, mn, taun, o));
}

double norm_En_conv_Ep_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double o){
  return norm_Ep_conv_En_En(ll, ul, m2, tau2, m1, tau1, mn, taun, o);
}

double norm_Ef_conv_Ep_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double o){
  return 0.5*(norm_En_conv_Ep_En(ll, ul, m1, tau1, m2, tau2, mn, taun, o)
	      +norm_Ep_conv_Ep_En(ll, ul, m1, tau1, m2, tau2, mn, taun, o));
}

double norm_Ef_conv_En_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double o){
  return 0.5*(norm_En_conv_En_En(ll, ul, m1, tau1, m2, tau2, mn, taun, o)
	      +norm_Ep_conv_En_En(ll, ul, m1, tau1, m2, tau2, mn, taun, o));
}

double norm_Ef_conv_Ef_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double o){
  return 0.5*(norm_Ef_conv_En_En(ll, ul, m1, tau1, m2, tau2, mn, taun, o)
	      +norm_Ef_conv_Ep_En(ll, ul, m1, tau1, m2, tau2, mn, taun, o));
}

//////////////////////////////////////////////////////////////////

double norm_Ep_conv_Ep_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEp_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*norm_Ep_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o)
	  +f2*norm_Ep_conv_gauss(ll-m1-m2, ul-m1-m2, tau2, m, s, o));
}

double norm_En_conv_En_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEn_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*norm_En_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o)
	  +f2*norm_En_conv_gauss(ll-m1-m2, ul-m1-m2, tau2, m, s, o));
}


double norm_Ep_conv_En_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s, const double o){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*norm_Ep_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o)
	  +fn*norm_En_conv_gauss(ll-m1-m2, ul-m1-m2, tau2, m, s, o));
}

double norm_En_conv_Ep_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s, const double o){
  return norm_Ep_conv_En_gauss(ll, ul, m2, tau2, m1, tau1, m, s, o);
}

double norm_Ef_conv_Ep_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s, const double o){
  return 0.5*(norm_En_conv_Ep_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o)
	      +norm_Ep_conv_Ep_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o));
}

double norm_Ef_conv_En_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s, const double o){
  return 0.5*(norm_En_conv_En_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o)
	      +norm_Ep_conv_En_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o));
}

double norm_Ef_conv_Ef_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s, const double o){
  return 0.5*(norm_Ef_conv_En_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o)
	      +norm_Ef_conv_Ep_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o));
}
//
double norm_xEp_conv_Ep_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    /* This case is not implemented yet */
    dout(Debugout::ERR,"libcnvl") << "[norm_xEp_conv_Ep_gauss] |tau1|==|tau2| case is not supported." << std::endl;
    abort();
    /* return norm_xxEp_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o); */
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*norm_xEp_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o);
  f += tau2*inv_tausub2*(-tau1*norm_Ep_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o)
			 +tau2*norm_Ep_conv_gauss(ll-m1-m2, ul-m1-m2, tau2, m, s, o));
  return f;
}

double norm_xEn_conv_En_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    /* This case is not implemented yet */
    dout(Debugout::ERR,"libcnvl") << "[norm_xEn_conv_En_gauss] |tau1|==|tau2| case is not supported." << std::endl;
    abort();
    /* return norm_xxEn_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o); */
  }
  const double inv_tausub = 1.0/(tau1-tau2);
  const double inv_tausub2 = inv_tausub*inv_tausub;

  double f = tau1*inv_tausub*norm_xEn_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o);
  f += tau2*inv_tausub2*(-tau1*norm_En_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o)
			 +tau2*norm_En_conv_gauss(ll-m1-m2, ul-m1-m2, tau2, m, s, o));
  return f;
}

double norm_xEp_conv_En_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    /* This case is not implemented yet */
    dout(Debugout::ERR,"libcnvl") << "[norm_xEp_conv_En_gauss] |tau1|==|tau2| case is not supported." << std::endl;
    abort();
    /* return norm_xxEp_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o); */
  }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*norm_xEp_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o);
  f += tau2*inv_tausum2*(tau1*norm_Ep_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o)
			 +tau2*norm_En_conv_gauss(ll-m1-m2, ul-m1-m2, tau2, m, s, o));
  return f;
}

double norm_xEn_conv_Ep_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    /* This case is not implemented yet */
    dout(Debugout::ERR,"libcnvl") << "[norm_xEn_conv_Ep_gauss] |tau1|==|tau2| case is not supported." << std::endl;
    abort();
    /* return norm_xxEp_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o); */
  }
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;

  double f = tau1*inv_tausum*norm_xEn_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o);
  f += tau2*inv_tausum2*(tau1*norm_En_conv_gauss(ll-m1-m2, ul-m1-m2, tau1, m, s, o)
			 +tau2*norm_Ep_conv_gauss(ll-m1-m2, ul-m1-m2, tau2, m, s, o));
  return f;
}

double norm_xEf_conv_Ep_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s, const double o){
  return 0.5*(norm_xEn_conv_Ep_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o)
	      +norm_xEp_conv_Ep_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o));
}

double norm_xEf_conv_En_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s, const double o){
  return 0.5*(norm_xEn_conv_En_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o)
	      +norm_xEp_conv_En_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o));
}

double norm_xEf_conv_Ef_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s, const double o){
  return 0.5*(norm_xEf_conv_En_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o)
	      +norm_xEf_conv_Ep_gauss(ll, ul, m1, tau1, m2, tau2, m, s, o));
}

// 

double norm_Ep_conv_Ep_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEp_conv_Ep_gauss(ll, ul, m1+m2, tau1, mp, taup, m, s, o);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*norm_Ep_conv_Ep_gauss(ll, ul, m1+m2, tau1, mp, taup, m, s, o)
	  +f2*norm_Ep_conv_Ep_gauss(ll, ul, m1+m2, tau2, mp, taup, m, s, o));
}

double norm_En_conv_En_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEn_conv_Ep_gauss(ll, ul, m1+m2, tau1, mp, taup, m, s, o);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*norm_En_conv_Ep_gauss(ll, ul, m1+m2, tau1, mp, taup, m, s, o)
	  +f2*norm_En_conv_Ep_gauss(ll, ul, m1+m2, tau2, mp, taup, m, s, o));
}


double norm_Ep_conv_En_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s, const double o){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*norm_Ep_conv_Ep_gauss(ll, ul, m1+m2, tau1, mp, taup, m, s, o)
	  +fn*norm_En_conv_Ep_gauss(ll, ul, m1+m2, tau2, mp, taup, m, s, o));
}

double norm_En_conv_Ep_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s, const double o){
  return norm_Ep_conv_En_Ep_gauss(ll, ul, m2, tau2, m1, tau1, mp, taup, m, s, o);
}

double norm_Ef_conv_Ep_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s, const double o){
  return 0.5*(norm_En_conv_Ep_Ep_gauss(ll, ul, m1, tau1, m2, tau2, mp, taup, m, s, o)
	      +norm_Ep_conv_Ep_Ep_gauss(ll, ul, m1, tau1, m2, tau2, mp, taup, m, s, o));
}

double norm_Ef_conv_En_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s, const double o){
  return 0.5*(norm_En_conv_En_Ep_gauss(ll, ul, m1, tau1, m2, tau2, mp, taup, m, s, o)
	      +norm_Ep_conv_En_Ep_gauss(ll, ul, m1, tau1, m2, tau2, mp, taup, m, s, o));
}

double norm_Ef_conv_Ef_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s, const double o){
  return 0.5*(norm_Ef_conv_En_Ep_gauss(ll, ul, m1, tau1, m2, tau2, mp, taup, m, s, o)
	      +norm_Ef_conv_Ep_Ep_gauss(ll, ul, m1, tau1, m2, tau2, mp, taup, m, s, o));
}

double norm_Ep_conv_Ep_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEp_conv_En_gauss(ll, ul, m1+m2, tau1, mn, taun, m, s, o);
  }
  double f1=0.0, f2=0.0;
  Ep_conv_Ep_coef(f1, f2, tau1, tau2);
  return (f1*norm_Ep_conv_En_gauss(ll, ul, m1+m2, tau1, mn, taun, m, s, o)
	  +f2*norm_Ep_conv_En_gauss(ll, ul, m1+m2, tau2, mn, taun, m, s, o));
}

double norm_En_conv_En_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s, const double o){
  if(tau1==tau2||tau1==-tau2||!finite(1.0/(tau1-tau2))){
    return norm_xEn_conv_En_gauss(ll, ul, m1+m2, tau1, mn, taun, m, s, o);
  }
  double f1=0.0, f2=0.0;
  En_conv_En_coef(f1, f2, tau1, tau2);
  return (f1*norm_En_conv_En_gauss(ll, ul, m1+m2, tau1, mn, taun, m, s, o)
	  +f2*norm_En_conv_En_gauss(ll, ul, m1+m2, tau2, mn, taun, m, s, o));
}


double norm_Ep_conv_En_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s, const double o){
  double fp=0.0, fn=0.0;
  Ep_conv_En_coef(fp, fn, tau1, tau2);
  return (fp*norm_Ep_conv_En_gauss(ll, ul, m1+m2, tau1, mn, taun, m, s, o)
	  +fn*norm_En_conv_En_gauss(ll, ul, m1+m2, tau2, mn, taun, m, s, o));
}

double norm_En_conv_Ep_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s, const double o){
  return norm_Ep_conv_En_En_gauss(ll, ul, m2, tau2, m1, tau1, mn, taun, m, s, o);
}

double norm_Ef_conv_Ep_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s, const double o){
  return 0.5*(norm_En_conv_Ep_En_gauss(ll, ul, m1, tau1, m2, tau2, mn, taun, m, s, o)
	      +norm_Ep_conv_Ep_En_gauss(ll, ul, m1, tau1, m2, tau2, mn, taun, m, s, o));
}

double norm_Ef_conv_En_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s, const double o){
  return 0.5*(norm_En_conv_En_En_gauss(ll, ul, m1, tau1, m2, tau2, mn, taun, m, s, o)
	      +norm_Ep_conv_En_En_gauss(ll, ul, m1, tau1, m2, tau2, mn, taun, m, s, o));
}

double norm_Ef_conv_Ef_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s, const double o){
  return 0.5*(norm_Ef_conv_En_En_gauss(ll, ul, m1, tau1, m2, tau2, mn, taun, m, s, o)
	      +norm_Ef_conv_Ep_En_gauss(ll, ul, m1, tau1, m2, tau2, mn, taun, m, s, o));
}

static double xXi_conv_gauss_by_int(nXi_t p_func,
				    const double t, const double xd,
				    const double mu, const double sigma){
  double area=0.0;
  if(sigma==0.0){
    return (*p_func)(t, xd);
  }
  /* integral (-20sigma -- +20 sigma)*/
  const double t_ini = -20.0*sigma+mu;
  const double t_end = +20.0*sigma+mu;
  const int ndiv = 400;
  const double t_del = (t_end-t_ini)/ndiv;

  for(int i=0;i<ndiv;++i){
    double t_run = t_ini+t_del*i;
    const double err = gaussian(t_run, t-mu, sigma);
    area += err*(*p_func)(t_run, xd);
  }
  area -= 0.5*gaussian(t_ini, t-mu, sigma)*(*p_func)(t_ini, xd);
  area -= 0.5*gaussian(t_end, t-mu, sigma)*(*p_func)(t_end, xd);
  area *= t_del;
  return area;
}



static double norm_Ax_sup(const double x,
			  const double tau, const double dm){

  const double dmt   = dm*x;
  const double dmtau = dm*tau;
  double f = -sin(dmt)-dmtau*cos(dmt);
  f *= tau*exp(-x/tau)/(1+dmtau*dmtau);
  return f;
}

static double norm_Mx_sup(const double x,
			  const double tau, const double dm){

  const double dmt   = dm*x;
  const double dmtau = dm*tau;
  double f = -cos(dmt)-dmtau*sin(dmt);
  f *= tau*exp(-x/tau)/(1+dmtau*dmtau);
  return f;
}

double norm_Ap(const double ll, const double ul,
	       const double tau, const double dm, const double o){
  double nll = (ll<ul) ? ll : ul;
  double nul = (ll<ul) ? ul : ll;
  nll -= o;
  nul -= o;
  if(nul<=0.0){
    return 0.0;
  }  
  double f = norm_Ax_sup(nul, tau, dm);
  f -= norm_Ax_sup((nll>0.0) ? nll : 0.0, tau, dm);
  return (ll<ul) ? f : -f;
}

double norm_Mp(const double ll, const double ul,
	       const double tau, const double dm, const double o){
  double nll = (ll<ul) ? ll : ul;
  double nul = (ll<ul) ? ul : ll;
  nll -= o;
  nul -= o;
  if(nul<=0.0){
    return 0.0;
  }  
  double f = norm_Mx_sup(nul, tau, dm);
  f -= norm_Mx_sup((nll>0.0) ? nll : 0.0, tau, dm);
  return (ll<ul) ? f : -f;
}

double norm_An(const double ll, const double ul,
	       const double tau, const double dm, const double o){
  double nll = (ll<ul) ? ll : ul;
  double nul = (ll<ul) ? ul : ll;
  nll -= o;
  nul -= o;
  if(nll>0.0){
    return 0.0;
  }  
  double f = norm_Ax_sup((nul<0.0) ? -nul : 0.0, tau, dm);
  f -= norm_Ax_sup(-nll, tau, dm);
  return (ll<ul) ? f : -f;
}

double norm_Mn(const double ll, const double ul,
	       const double tau, const double dm, const double o){
  double nll = (ll<ul) ? ll : ul;
  double nul = (ll<ul) ? ul : ll;
  nll -= o;
  nul -= o;
  if(nll>0.0){
    return 0.0;
  }  
  double f = norm_Mx_sup(-nll, tau, dm);
  f -= norm_Mx_sup((nul<0.0) ? -nul : 0.0, tau, dm);
  return (ll<ul) ? f : -f;
}

double norm_Af(const double ll, const double ul,
	       const double tau, const double dm, const double o){
  double f = norm_An(ll+o, ul+o, tau, dm);
  f += norm_Ap(ll+o, ul+o, tau, dm);
//   f *= 0.5;
  return f;
}

double norm_Mf(const double ll, const double ul,
	       const double tau, const double dm, const double o){
  double f = norm_Mn(ll+o, ul+o, tau, dm);
  f += norm_Mp(ll+o, ul+o, tau, dm);
//   f *= 0.5;
  return f;
}


} /* end of extern "C" */
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
