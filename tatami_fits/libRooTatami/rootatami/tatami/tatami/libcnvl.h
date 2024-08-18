//
//   libcnvl.h
//     --- Description: Library for Convolurion (See belle note #307.)
//         Note that the definition of Ap, An, Mp and Mp is different 
//         from that of belle note #370 by factor tau.
//
//   $Id: libcnvl.h 9932 2006-11-12 14:26:53Z katayama $
//
//   $Log$
//   Revision 1.1  2002/09/05 01:27:07  katayama
//   New package tatami from Nakadaira/Sumisawa san
//
//   Revision 1.5  2002/05/09 06:25:12  nakadair
//   Change indent.(No significant change.)
//
//   Revision 1.4  2002/05/07 21:38:23  nakadair
//   Add convolution with Resolution function and fortran interface.
//
//   Revision 1.3  2002/05/06 10:38:42  nakadair
//   Implement norm_xxE?_conv_gauss functions.
//
//   Revision 1.2  2002/02/21 19:16:13  nakadair
//   Update belle-note No.
//
//   Revision 1.1.1.1  2002/02/21 18:56:07  nakadair
//   Initial import(Branch from lmfit)
//
//

#ifndef __LIBCNVL_H__
#define __LIBCNVL_H__

extern const char* libcnvl_id;

#include "belle.h"
#include "tatami/cmlxerf.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//  #define USE_CERNLIB_ERFC
#define USE_CERNLIB_WERF

#ifdef __cplusplus
#define OFFSET_DEFAULT_ARGUMENT = 0.0
extern "C" {
#else /* __cplusplus */
#define OFFSET_DEFAULT_ARGUMENT
#endif /* __cplusplus */


double rewerf(double x, double y);
double imwerf(double x, double y);
double recexp(double x, double y);
double imcexp(double x, double y);

#ifdef USE_CERNLIB_ERFC
extern double erfc_(double*);
extern double derfc_(double*);
#endif /* USE_CERNLIB_ERFC */

double derfc(double x);

double DiracDelta(const double x);

double gaussian(const double x, const double m, const double s);

double norm_gaussian(const double ll, const double ul, const double m, const double s);

double norm_gaussian_w_cutoff(const double cutoff, const double m, const double s);

double Ep(const double t, const double tau);

double En(const double t, const double tau);

double Ef(const double t, const double tau);

double Enp(const double t, const double tau_n, const double tau_p);

double nEp_conv_gauss(const double t, 
		      const double m, const double s);

double nEn_conv_gauss(const double t, 
		      const double m, const double s);

double Ep_conv_gauss(const double t, const double tau,
		     const double m, const double s);

double En_conv_gauss(const double t, const double tau,
		     const double m, const double s);

double Ef_conv_gauss(const double t, const double tau,
		     const double m, const double s);
double Enp_conv_gauss(const double t,
		      const double tau_n, const double tau_p,
		      const double m, const double s);

double nMp(const double t, const double xd);

double Mp(const double t,
	  const double tau, const double dm);

double nMn(const double t, const double xd);

double Mn(const double t,
	  const double tau, const double dm);

double nMf(const double t, const double xd);

double Mf(const double t,
	  const double tau, const double dm);

double nAp(const double t, const double xd);

double Ap(const double t,
	  const double tau, const double dm);

double nAn(const double t, const double xd);

double An(const double t,
	  const double tau, const double dm);

double nAf(const double t, const double xd);

double Af(const double t,
	  const double tau, const double chi);

double nMp_conv_gauss(const double t, const double xd,
		      const double m, const double s);

double nMn_conv_gauss(const double t, const double xd,
		      const double m, const double s);

double nAp_conv_gauss(const double t, const double xd,
		      const double m, const double s);

double nAn_conv_gauss(const double t, const double xd,
		      const double m, const double s);

double Mp_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s);

double Mn_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s);

double Mf_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s);

double Ap_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s);

double An_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s);

double Af_conv_gauss(const double t,
		     const double tau, const double dm,
		     const double m, const double s);

double CoefEEs(const double tau1, const double tau2);

double CoefEEo(const double tau1, const double tau2);

double CoefEEf(const double tau1, const double tau2);

double sum_sigma(const double s1, const double s2);
double sum_sigma3(const double s1, const double s2, const double s3);

double norm_nEp(const double ll, const double ul,
		const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nEn(const double ll, const double ul,
		const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nEf(const double ll, const double ul,
		const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ep(const double ll, const double ul,
		const double tau,
	       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En(const double ll, const double ul,
		const double tau,
	       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef(const double ll, const double ul,
	       const double tau,
	       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ap(const double ll, const double ul,
	       const double tau, const double dm,
	       const double o OFFSET_DEFAULT_ARGUMENT);
double norm_Mp(const double ll, const double ul,
	       const double tau, const double dm,
	       const double o OFFSET_DEFAULT_ARGUMENT);
double norm_An(const double ll, const double ul,
	       const double tau, const double dm,
	       const double o OFFSET_DEFAULT_ARGUMENT);
double norm_Mn(const double ll, const double ul,
	       const double tau, const double dm,
	       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Af(const double ll, const double ul,
	       const double tau, const double dm,
	       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Mf(const double ll, const double ul,
	       const double tau, const double dm,
	       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nEp_conv_gauss(const double ll, const double ul,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nEn_conv_gauss(const double ll, const double ul,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ep_conv_gauss(const double ll, const double ul, const double tau,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_gauss(const double ll, const double ul, const double tau,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nEf_conv_gauss(const double ll, const double ul,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_gauss(const double ll, const double ul, const double tau,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Enp_conv_gauss(const double ll, const double ul, 
			   const double tau_n, const double tau_p,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nAn_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nAp_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nAf_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nMn_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nMp_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_nMf_conv_gauss(const double ll, const double ul, const double xd,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_An_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ap_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Af_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Mn_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Mp_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Mf_conv_gauss(const double ll, const double ul, 
			  const double tau, const double dm,
			  const double m, const double s,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double Ap_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double An_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mp_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mn_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Ap_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double An_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mp_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mn_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Ap_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double An_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mp_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mn_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Af_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Af_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mf_conv_Ep(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mf_conv_En(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Af_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Mf_conv_Ef(const double t, const double m1,
		  const double tau1, const double dm,
		  const double m2, const double tau2);

double Ef_conv_EnpGauss(const double t, const double tau,
			const double taunk, const double taupk,
			const double mu, const double sigma,
			const double ll, const double ul);

void Ep_conv_En_coef(double& fp, double& fn,
		     const double tau_p, const double tau_n);

void En_conv_Ep_coef(double& fn, double& fp,
		     const double tau_n, const double tau_p);

void Ep_conv_Ep_coef(double& fp1, double& fp2,
		     const double tau_p1, const double tau_p2);

void En_conv_En_coef(double& fn1, double& fn2,
		     const double tau_n1, const double tau_n2);

double Ep_conv_En(const double t,
		  const double m_p, const double tau_p,
		  const double m_n, const double tau_n);

double En_conv_Ep(const double t,
		  const double m_n, const double tau_n,
		  const double m_p, const double tau_p);

double xEp(const double t, const double tau);
double xEn(const double t, const double tau);
double xEf(const double t, const double tau);

double norm_xEp(const double ll, const double ul,
		const double tau,
		const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEn(const double ll, const double ul,
		const double tau,
		const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf(const double ll, const double ul,
		const double tau,
		const double o OFFSET_DEFAULT_ARGUMENT);

double xEp_conv_gauss(const double t, const double tau,
		      const double m, const double s);
double xEn_conv_gauss(const double t, const double tau,
		      const double m, const double s);

double xEf_conv_gauss(const double t, const double tau,
		      const double m, const double s);

double xEp_conv_gauss_w_meanshift(const double t,
				  const double m1, const double tau,
				  const double m, const double s);

double xEn_conv_gauss_w_meanshift(const double t,
				  const double m1, const double tau,
				  const double m, const double s);

double xEf_conv_gauss_w_meanshift(const double t,
				  const double m1, const double tau,
				  const double m, const double s);

double xxEp(const double t, const double tau);
double xxEn(const double t, const double tau);
double xxEf(const double t, const double tau);

double norm_xxEp(const double ll, const double ul,
		const double tau,
		 const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xxEn(const double ll, const double ul,
		const double tau,
		 const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xxEf(const double ll, const double ul,
		const double tau,
		 const double o OFFSET_DEFAULT_ARGUMENT);


double xxEp_conv_gauss(const double t, const double tau,
		       const double m, const double s);
double xxEn_conv_gauss(const double t, const double tau,
		       const double m, const double s);

double xxEf_conv_gauss(const double t, const double tau,
		       const double m, const double s);

double xxEp_conv_gauss_w_meanshift(const double t,
				   const double m1, const double tau,
				   const double m, const double s);

double xxEn_conv_gauss_w_meanshift(const double t,
				   const double m1, const double tau,
				   const double m, const double s);

double xxEf_conv_gauss_w_meanshift(const double t,
				   const double m1, const double tau,
				   const double m, const double s);

double Ep_conv_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double En_conv_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double Ep_conv_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double En_conv_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double Ep_conv_Ef(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double En_conv_Ef(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double Ef_conv_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double Ef_conv_En(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);

double Ef_conv_Ef(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2);
//
double xEp_conv_Ep(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

double xEn_conv_En(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

double xEp_conv_En(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

double xEn_conv_Ep(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

double xEf_conv_Ep(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

double xEf_conv_En(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);
double xEf_conv_Ef(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

double xEp_conv_Ef(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

double xEn_conv_Ef(const double t,
		   const double m1, const double tau1,
		   const double m2, const double tau2);

//

double Ep_conv_Ep_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2,
		  const double mp, const double taup);

double En_conv_En_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2,
		  const double mp, const double taup);

double Ep_conv_En_Ep(const double t,
		  const double m1, const double tau1,
		  const double m2, const double tau2,
		  const double mp, const double taup);

double En_conv_Ep_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup);

double Ef_conv_Ep_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup);

double Ef_conv_En_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup);

double Ef_conv_Ef_Ep(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mp, const double taup);

double Ep_conv_Ep_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun);

double En_conv_En_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun);

double Ep_conv_En_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun);

double En_conv_Ep_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun);

double Ef_conv_Ep_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun);

double Ef_conv_En_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun);

double Ef_conv_Ef_En(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mn, const double taun);

// 
double Ep_conv_Ep_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf);

double En_conv_En_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf);

double Ep_conv_En_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf);

double En_conv_Ep_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf);

double Ef_conv_Ep_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf);

double Ef_conv_En_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf);

double Ef_conv_Ef_Ef(const double t,
		     const double m1, const double tau1,
		     const double m2, const double tau2,
		     const double mf, const double tauf);

double xEp_conv_Ep_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf);

double xEn_conv_En_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf);

double xEp_conv_En_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf);

double xEn_conv_Ep_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf);

double xEf_conv_Ep_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf);

double xEf_conv_En_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf);

double xEf_conv_Ef_Ef(const double t,
		      const double m1, const double tau1,
		      const double m2, const double tau2,
		      const double mf, const double tauf);

double Ep_conv_Ep_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s);

double En_conv_En_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s);

double Ep_conv_En_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s);

double En_conv_Ep_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s);

double Ef_conv_Ep_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s);

double Ef_conv_En_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s);

double Ef_conv_Ef_gauss(const double t,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double m, const double s);
//
double xEp_conv_Ep_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s);

double xEn_conv_En_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s);

double xEp_conv_En_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s);

double xEn_conv_Ep_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s);

double xEf_conv_Ep_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s);

double xEf_conv_En_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s);

double xEf_conv_Ef_gauss(const double t,
			 const double m1, const double tau1,
			 const double m2, const double tau2,
			 const double m, const double s);

// 

double Ep_conv_Ep_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s);

double En_conv_En_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s);


double Ep_conv_En_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s);

double En_conv_Ep_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s);

double Ef_conv_Ep_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s);

double Ef_conv_En_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s);

double Ef_conv_Ef_Ep_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mp, const double taup,
			   const double m, const double s);

double Ep_conv_Ep_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s);

double En_conv_En_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s);


double Ep_conv_En_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s);

double En_conv_Ep_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s);

double Ef_conv_Ep_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s);

double Ef_conv_En_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s);

double Ef_conv_Ef_En_gauss(const double t,
			   const double m1, const double tau1,
			   const double m2, const double tau2,
			   const double mn, const double taun,
			   const double m, const double s);
// ----------------------------------------------------------------------

double norm_xEp_conv_gauss(const double ll, const double ul, const double tau,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEn_conv_gauss(const double ll, const double ul, const double tau,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf_conv_gauss(const double ll, const double ul, const double tau,
			   const double m, const double s,
			   const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xxEp_conv_gauss(const double ll, const double ul, const double tau,
			    const double m, const double s,
			    const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xxEn_conv_gauss(const double ll, const double ul, const double tau,
			    const double m, const double s,
			    const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xxEf_conv_gauss(const double ll, const double ul, const double tau,
			    const double m, const double s,
			    const double o OFFSET_DEFAULT_ARGUMENT);


double norm_Ep_conv_Ep(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2,
		       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_En(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2,
		       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ep_conv_En(const double ll, const double ul,
		       const double m_p, const double tau_p,
		       const double m_n, const double tau_n,
		       const double o OFFSET_DEFAULT_ARGUMENT);


double norm_En_conv_Ep(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2,
		       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ep(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2,
		       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_En(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2,
		       const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ef(const double ll, const double ul,
		       const double m1, const double tau1,
		       const double m2, const double tau2,
		       const double o OFFSET_DEFAULT_ARGUMENT);
//
double norm_xEp_conv_Ep(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEn_conv_En(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEp_conv_En(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEn_conv_Ep(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf_conv_Ep(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf_conv_En(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf_conv_Ef(const double ll, const double ul,
			const double m1, const double tau1,
			const double m2, const double tau2,
			const double o OFFSET_DEFAULT_ARGUMENT);

// 

double norm_Ep_conv_Ep_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_En_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup,
			  const double o OFFSET_DEFAULT_ARGUMENT);


double norm_Ep_conv_En_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_Ep_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ep_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_En_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ef_Ep(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mp, const double taup,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ep_conv_Ep_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_En_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun,
			  const double o OFFSET_DEFAULT_ARGUMENT);


double norm_Ep_conv_En_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_Ep_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ep_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun,
			  const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_En_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun, const double oe);

double norm_Ef_conv_Ef_En(const double ll, const double ul,
			  const double m1, const double tau1,
			  const double m2, const double tau2,
			  const double mn, const double taun,
			  const double o OFFSET_DEFAULT_ARGUMENT);

//////////////////////////////////////////////////////////////////

double norm_Ep_conv_Ep_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s,
			     const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_En_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s,
			     const double o OFFSET_DEFAULT_ARGUMENT);


double norm_Ep_conv_En_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s,
			     const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_Ep_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s,
			     const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ep_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s,
			     const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_En_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s,
			     const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ef_gauss(const double ll, const double ul,
			     const double m1, const double tau1,
			     const double m2, const double tau2,
			     const double m, const double s,
			     const double o OFFSET_DEFAULT_ARGUMENT);
//
double norm_xEp_conv_Ep_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s,
			      const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEn_conv_En_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s,
			      const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEp_conv_En_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s,
			      const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEn_conv_Ep_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s,
			      const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf_conv_Ep_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s,
			      const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf_conv_En_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s,
			      const double o OFFSET_DEFAULT_ARGUMENT);

double norm_xEf_conv_Ef_gauss(const double ll, const double ul,
			      const double m1, const double tau1,
			      const double m2, const double tau2,
			      const double m, const double s,
			      const double o OFFSET_DEFAULT_ARGUMENT);

// 

double norm_Ep_conv_Ep_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_En_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);


double norm_Ep_conv_En_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_Ep_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ep_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_En_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ef_Ep_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mp, const double taup,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ep_conv_Ep_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_En_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);


double norm_Ep_conv_En_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_En_conv_Ep_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ep_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_En_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

double norm_Ef_conv_Ef_En_gauss(const double ll, const double ul,
				const double m1, const double tau1,
				const double m2, const double tau2,
				const double mn, const double taun,
				const double m, const double s,
				const double o OFFSET_DEFAULT_ARGUMENT);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __LIBCNVL_H__ */

