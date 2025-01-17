//
//   dt_resolution.cc
//     --- Description: Library for Resolution Function of dt
//   $Id: dt_resolution.cc 11208 2011-03-01 14:53:58Z hitoshi $
//
//   $Log$
//   Revision 1.5  2004/10/19 07:09:55  kohji
//   1) SVD1 and SVD2 resolution parameters are updated for 2004 publication.
//   2) "expno" is added to add_outlier and dtres_systematics functions to treat
//   difference between SVD1 and SVD2 properly[cpfit_ml:0771].
//   3) New functions are added for cosh/sinh terms
//    (norm_)HAfRkRdetRnp_{fullrec,partial} for HAf=exp(-|t|Gamma)sinh(DGamma/2 t)
//    (norm_)HMfRkRdetRnp_{fullrec,partial} for HMf=exp(-|t|Gamma)cosh(DGamma/2 t)
//
//   Revision 1.4  2003/10/01 01:24:23  katayama
//   compile with CC
//
//   Revision 1.3  2003/08/06 06:43:32  katayama
//   New version from Nakadaira san
//
//   Revision 1.2  2002/12/10 12:27:49  sumisawa
//   Bug fix in (norm_)AfRkRdet(...), (norm_)MfRkRdet(...) and modify test_tatami_jpsiks.cc to speed up
//
//   Revision 1.1  2002/09/05 01:27:04  katayama
//   New package tatami from Nakadaira/Sumisawa san
//
//   Revision 1.19  2002/07/12 09:56:54  nakadair
//   Add new functions: EfRdet(...), AfRdet(...), MfRdet(...),
//   norm_EfRdet(...), norm_AfRdet(...) and norm_MfRdet(...).
//
//   Revision 1.18  2002/07/12 09:18:19  nakadair
//   Bug fix in EfRkRdetRnp_partial(...), AfRkRdetRnp_partial(...),
//   MfRkRdetRnp_partial(...), norm_EfRkRdetRnp_partial(...),
//   norm_AfRkRdetRnp_partial(...) and norm_MfRkRdetRnp_partial(...).
//
//   Revision 1.17  2002/07/11 09:41:03  nakadair
//   Move parameters constants to dt_resolution_const.{h,cc}.
//   Add new function Add(_)Outlier(_)With(_)Bkg(...).
//   Support the systematics study based on new resolution parameter
//   determination scheme.
//
//   Revision 1.16  2002/07/10 06:50:01  nakadair
//   Enable systematic study on delta-t cut-off.
//
//   Revision 1.15  2002/07/06 21:38:00  nakadair
//   Bug fix on treatment of single track case in RascRnp(...) and
//   norm_RascRnp(...).  Correct wrong comment in Rrec(...).
//   Change the constrait to Smain in Rasc_param(...) and Rrec_param(...)
//
//   Revision 1.14  2002/07/06 05:47:39  nakadair
//   Update constant dt_resol_gloval::inv_bgc.
//
//   Revision 1.13  2002/07/05 06:36:22  nakadair
//   Update for GCC3.
//
//   Revision 1.12  2002/07/05 06:15:29  nakadair
//   Add Calc_AkCk_Mass(...), Calc_AkCk_Flavor(...), Calc_AkCk(...), and Add(_)Outlier(...).
//   Add new members and C interfece to class dt_resol_global.
//
//   Revision 1.11  2002/07/02 23:58:22  nakadair
//   Add utility functions for systematics study.
//   Update Rk_parameters.
//
//   Revision 1.10  2002/06/28 03:49:53  nakadair
//   Update constant in dtres_param_default and dtres_param_default_mc.
//
//   Revision 1.9  2002/06/27 13:29:49  nakadair
//   Change Rk_partial_param(...), Rnp_param(...), swap_rnp_param(...),
//   Rasc_param(...) and Rrec_param(...) from inline to global function.
//   Add class dt_resol_global.
//
//   Revision 1.8  2002/06/17 16:17:00  nakadair
//   Bug fix in swap_rnp_param(...).
//
//   Revision 1.7  2002/06/17 10:37:28  nakadair
//   Bug fix in RdetRnp(...), XfRkRdetRnp_xxx(...),
//   norm_RdetRnp(...) and  XfRkRdetRnp_xxx(...).
//
//   Revision 1.6  2002/06/12 15:18:05  nakadair
//   Bug fix in norm_?fRkRdetRnp_{partial,full_sup}(...). Correct the coefficients of exponential terms.
//
//   Revision 1.5  2002/06/08 14:55:56  nakadair
//   Bug fix in Rnp(...) and so on. (See ChangeLog.)
//
//   Revision 1.4  2002/05/10 13:21:24  nakadair
//   Change output conversion specification of double variable.
//   Change to include float.h.
//   Bug fix in EfRkRdetRnp_partial_sup(...), norm_EfRkRdetRnp_partial_sup(...)
//   and add_outlier(...).
//
//   Revision 1.3  2002/05/09 07:08:34  nakadair
//   Add add_outlier(...).
//
//   Revision 1.2  2002/05/09 06:23:55  nakadair
//   Change void calc_akck(...) from inline to global function.
//
//   Revision 1.1  2002/05/07 21:39:20  nakadair
//   Add convolution with resolution function and fortran interface.
//
//

#include "belle.h"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include "tatami/libcnvl.h"
#include "tatami/dt_resolution.h"
#include "tatami/conv_coef.h"
#include "tatami/conv_coef3.h"

#include <iostream>
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
#ifndef __POSTIVEPARAM_MIN__
#define __POSTIVEPARAM_MIN__ FLT_MIN
/* #define __POSTIVEPARAM_MIN__ 1.0e-15 */
#endif /* __POSTIVEPARAM_MIN__ */

#ifndef __FORCE_PARAM_CONSTRAINT__
#define __FORCE_PARAM_CONSTRAINT__ 1
#endif /* __FORCE_PARAM_CONSTRAINT__ */

inline void constraint(double* x, const double ll, const double ul){
  if(*x<ll){
    *x = ll;
  }else if(*x>ul){
    *x = ul;
  }
  return;
}

inline void constraint(double* x, const double ll){
  if(*x<ll){
    *x = ll;
  }
  return;
}

__LINKAGE__

int cpdefpar_(int* expmc){
  const int siz = sizeof(dtres_param_t);
  extern dtres_param_t respar_;
  if(*expmc==1){
    std::memcpy(&respar_, &dtres_param_default, siz);
    return 0;
  }else if (*expmc==2){
    std::memcpy(&respar_, &dtres_param_default_mc, siz);
    return 0;
  }
  return -1;
}


void calc_akck(const double cos_theta_b,
	       const double beta, const double Eb_cms, const double pb_cms,
	       double* ak, double* ck){
  const double& inv_mb = dt_resol_global::inv_mb; /* 1.0 / 5.279 */
  *ak = Eb_cms*inv_mb;
  *ck = pb_cms*cos_theta_b/beta*inv_mb;
}

/* Rk(t) */

double Rk_fullrec(const double x,
		  const double dt,
		  const double tau,
		  const double ak, const double ck){
  const double adt = fabs(dt);
  //  if(x==0.0){
  if(ck==0.0){
    return DiracDelta(x-(ak-1)*dt);
  }
  const double xx = x-(ak-1)*dt-ck*adt;
  return (ck>0.0) ? Ep(xx, ck*tau) : En(xx, -ck*tau);
}

void Rk_partial_param(const double dz,
		      const dtres_param_t * const param,
		      double* tau_k_1_p, double* tau_k_2_p,
		      double* tau_k_1_n, double* tau_k_2_n,
		      double* sigma_k, double* fk2){
  const double& inv_bgc = dt_resol_global::inv_bgc;;
  const double adt = fabs(dz*inv_bgc);
  *tau_k_1_p = (param->tau_k_1_p)[0] + (param->tau_k_1_p)[1]*adt;
  *tau_k_2_p = (param->tau_k_2_p)[0] + (param->tau_k_2_p)[1]*adt;
  *tau_k_1_n = (param->tau_k_1_n)[0] + (param->tau_k_1_n)[1]*adt;
  *tau_k_2_n = (param->tau_k_2_n)[0] + (param->tau_k_2_n)[1]*adt;
  *sigma_k   = (param->sigma_k)[0]   + (param->sigma_k)[1]*adt;
  *fk2       = param->fk2;
#if defined(__FORCE_PARAM_CONSTRAINT__)&&__FORCE_PARAM_CONSTRAINT__
  constraint(tau_k_1_p, __POSTIVEPARAM_MIN__);
  constraint(tau_k_2_p, __POSTIVEPARAM_MIN__);
  constraint(tau_k_1_n, __POSTIVEPARAM_MIN__);
  constraint(tau_k_2_n, __POSTIVEPARAM_MIN__);
  constraint(sigma_k,   __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  return;
}

double Rk_partial(const double x, const double dz,
		  const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  const double Li_1 = Enp_conv_gauss(x, tau_k_1_n, tau_k_1_p, 0.0, sigma_k);
  const double Li_2 = Enp_conv_gauss(x, tau_k_2_n, tau_k_2_p, 0.0, sigma_k);
  const double Li = (1.0-fk2)*Li_1 + fk2*Li_2;
  return Li;
}

inline void calc_vtxparam_asc(const int ntrk_asc, const double sz_asc,
			      const double chisq_z_asc, const int ndf_z_asc,
			      double*  xi_asc, double* st_asc){
  const double& inv_bgc = dt_resol_global::inv_bgc;;
  *xi_asc = (ndf_z_asc > 0 ? chisq_z_asc/ndf_z_asc : 1);
  *st_asc = sz_asc*inv_bgc;
  return;
}

inline void calc_vtxparam_rec(const int ntrk_rec, const double sz_rec,
			      const double chisq_z_rec, const int ndf_z_rec,
			      double*  xi_rec, double* st_rec){
  const double& inv_bgc = dt_resol_global::inv_bgc;;
  *xi_rec = (ndf_z_rec > 0 ? chisq_z_rec/ndf_z_rec : 1);
  *st_rec = sz_rec*inv_bgc;
  return;
}

void Rasc_param(const int ntrk_asc,
		const double xi_asc, const double st_asc,
		const dtres_param_t * const param,
		double* ftail_asc,
		double* mu_main_asc, double* Smain_asc,
		double* mu_tail_asc, double* Stail_asc){
  if(ntrk_asc>1){ /*  multiple track vertex  */
    *ftail_asc = (param->ftl_asc_mlt)[0] + (param->ftl_asc_mlt)[1] * xi_asc;
    *Smain_asc = ((param->Sasc)[0] + (param->Sasc)[1] * xi_asc ) * st_asc;
    *Stail_asc = param->Stl_asc_mlt*(*Smain_asc);
    constraint(ftail_asc, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }else{ /*  single track vertex */
    *ftail_asc = param->ftl_asc;
    *Smain_asc = param->Smn_asc * st_asc;
    *Stail_asc = param->Stl_asc * st_asc;
    constraint(ftail_asc, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_asc, __POSTIVEPARAM_MIN__);
    constraint(Stail_asc, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
  *mu_main_asc = 0.0;
  *mu_tail_asc = 0.0;
  return;
}

double Rasc(const double x,
	    const int ntrk_asc, const double sz_asc,
	    const double chisq_z_asc, const int ndf_z_asc,
	    const dtres_param_t * const param){

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double Li_mn = gaussian(x, mu_main_asc, Smain_asc);
  if(ftail_asc>0.0){ /* Mainly single track case */
    const double Li_tl = gaussian(x, mu_tail_asc, Stail_asc);
    const double Li    = (1.0-ftail_asc)*Li_mn + ftail_asc*Li_tl;
    return Li;
  }
  return Li_mn;
}

void Rrec_param(const int ntrk_rec,
		const double xi_rec, const double st_rec,
		const dtres_param_t * const param,
		double* ftail_rec,
		double* mu_main_rec, double* Smain_rec,
		double* mu_tail_rec, double* Stail_rec){
  if(ntrk_rec>1){ /*  multiple track vertex  */
    *ftail_rec = (param->ftl_rec_mlt)[0] + (param->ftl_rec_mlt)[1] * xi_rec;
    *Smain_rec = ((param->Srec)[0] + (param->Srec)[1] * xi_rec ) * st_rec;
    *Stail_rec = param->Stl_rec_mlt*(*Smain_rec);
    constraint(ftail_rec, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_rec, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }else{ /*  single track vertex */
    *ftail_rec = param->ftl_rec;
    *Smain_rec = param->Smn_rec * st_rec;
    *Stail_rec = param->Stl_rec * st_rec;
    constraint(ftail_rec, 0.0, 1.0);
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(Smain_rec, __POSTIVEPARAM_MIN__);
    constraint(Stail_rec, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
  *mu_main_rec = 0.0;
  *mu_tail_rec = 0.0;
  return;
}

double Rrec(const double x,
	    const int ntrk_rec, const double sz_rec,
	    const double chisq_z_rec, const int ndf_z_rec,
	    const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  const double Li_mn = gaussian(x, mu_main_rec, Smain_rec);
  if(ftail_rec>0.0){ /* Mainly single track case */
    const double Li_tl = gaussian(x, mu_tail_rec, Stail_rec);
    const double Li    = (1.0-ftail_rec)*Li_mn + ftail_rec*Li_tl;
    return Li;
  }
  return Li_mn;
}

void Rnp_param(const int flavor, const int ntrk_asc, const int keeptagl,
	       const double Smain_asc, const double Stail_asc,
	       const dtres_param_t * const param,
	       double* fd, double* fp,
	       double* tau_np_p, double* tau_np_n,
	       double* tau_np_p_tl, double* tau_np_n_tl){

  const int flvidx = (flavor==0) ? 0 : 1;
  const int keeptlidx = (keeptagl==0)? 0 : 1;
  if(ntrk_asc>1){ /*  multiple track vertex  */
    *fd = (param->fd_np_mlt)[flvidx][keeptlidx];
    *fp = (param->fp_np_mlt)[flvidx];
    *tau_np_p = (param->tau_np_p_mlt)[flvidx][0] + (param->tau_np_p_mlt)[flvidx][1]*Smain_asc;
    *tau_np_n = (param->tau_np_n_mlt)[flvidx][0] + (param->tau_np_n_mlt)[flvidx][1]*Smain_asc;
    *tau_np_p_tl = 0.0;
    *tau_np_n_tl = 0.0;
  }else{ /*  single track vertex */
    *fd = (param->fd_np_sgl)[flvidx][keeptlidx];
    *fp = (param->fp_np_sgl)[flvidx];
    *tau_np_p    = (param->tau_np_p_sgl)[flvidx][0] + (param->tau_np_p_sgl)[flvidx][1]*Smain_asc;
    *tau_np_n    = (param->tau_np_n_sgl)[flvidx][0] + (param->tau_np_n_sgl)[flvidx][1]*Smain_asc;
    *tau_np_p_tl = (param->tau_np_p_sgl)[flvidx][0] + (param->tau_np_p_sgl)[flvidx][1]*Stail_asc;
    *tau_np_n_tl = (param->tau_np_n_sgl)[flvidx][0] + (param->tau_np_n_sgl)[flvidx][1]*Stail_asc;
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
    constraint(tau_np_p_tl, __POSTIVEPARAM_MIN__);
    constraint(tau_np_n_tl, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  constraint(tau_np_p,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_n,  __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  constraint(fp, 0.0, 1.0);
  constraint(fd, 0.0, 1.0);
  return;
}

void Rnp_param_03(const int flavor, const int ntrk_asc, const int keeptagl,
		  const double xi_asc, const double st_asc,
		  const double Smain_asc, const double Stail_asc,
		  const dtres_param_t * const param,
		  double* fd, double* fp,
		  double* tau_np_p, double* tau_np_n,
		  double* tau_np_p_tl, double* tau_np_n_tl){

  const int flvidx = (flavor==0) ? 0 : 1;
  const int keeptlidx = (keeptagl==0)? 0 : 1;
  if(ntrk_asc>1){ /*  multiple track vertex  */
    const double Smain_asc_in
      = (param->Snp<=DTRES_EXTERAM_THRE) ? Smain_asc
      : (1.0+param->Snp*xi_asc)*st_asc;
    *fd = (param->fd_np_mlt)[flvidx][keeptlidx];
    *fp = (param->fp_np_mlt)[flvidx];
    *tau_np_p = param->Snp_global*((param->tau_np_p_mlt)[flvidx][0] + (param->tau_np_p_mlt)[flvidx][1]*Smain_asc_in);
    *tau_np_n = param->Snp_global*((param->tau_np_n_mlt)[flvidx][0] + (param->tau_np_n_mlt)[flvidx][1]*Smain_asc_in);
    const double Stail_asc_in
      = (param->Snp<=DTRES_EXTERAM_THRE) ? Stail_asc
      : (1.0+param->Snp*xi_asc)*st_asc;
    *tau_np_p_tl = param->Snp_global*((param->tau_np_p_mlt)[flvidx][0] + (param->tau_np_p_mlt)[flvidx][1]*Stail_asc_in);
    *tau_np_n_tl = param->Snp_global*((param->tau_np_n_mlt)[flvidx][0] + (param->tau_np_n_mlt)[flvidx][1]*Stail_asc_in);
  }else{ /*  single track vertex */
    const double Smain_asc_in
      = (param->Snp<=DTRES_EXTERAM_THRE) ? Smain_asc : st_asc;
    const double Stail_asc_in
      = (param->Snp<=DTRES_EXTERAM_THRE) ? Stail_asc : st_asc;
    *fd = (param->fd_np_sgl)[flvidx][keeptlidx];
    *fp = (param->fp_np_sgl)[flvidx];
    *tau_np_p    = param->Snp_global*((param->tau_np_p_sgl)[flvidx][0] + (param->tau_np_p_sgl)[flvidx][1]*Smain_asc_in);
    *tau_np_n    = param->Snp_global*((param->tau_np_n_sgl)[flvidx][0] + (param->tau_np_n_sgl)[flvidx][1]*Smain_asc_in);
    *tau_np_p_tl = param->Snp_global*((param->tau_np_p_sgl)[flvidx][0] + (param->tau_np_p_sgl)[flvidx][1]*Stail_asc_in);
    *tau_np_n_tl = param->Snp_global*((param->tau_np_n_sgl)[flvidx][0] + (param->tau_np_n_sgl)[flvidx][1]*Stail_asc_in);
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  constraint(tau_np_p,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_n,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_p_tl, __POSTIVEPARAM_MIN__);
  constraint(tau_np_n_tl, __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  constraint(fp, 0.0, 1.0);
  constraint(fd, 0.0, 1.0);
  return;
}

void Rnp_param_10(const int flavor, const int ntrk_asc, const int keeptagl,
		  const double xi_asc, const double st_asc,
		  const dtres_param_t * const param,
		  double* fd, double* fp,
		  double* tau_np_p, double* tau_np_n,
		  double* tau_np_p_tl, double* tau_np_n_tl){

  const int flvidx = (flavor==0) ? 0 : 1;
  const int keeptlidx = (keeptagl==0)? 0 : 1;

  if(ntrk_asc>1){ /*  multiple track vertex  */
    double st_asc_forf = (st_asc > dt_resol_global::rnp_kink_st? dt_resol_global::rnp_kink_st: st_asc);
    double xi_asc_forf = (xi_asc > dt_resol_global::rnp_kink_xi? dt_resol_global::rnp_kink_xi: xi_asc);

    double f_neg = (param->fn_np_mlt)[flvidx];
    double f_delta = (param->fd_np_mlt)[flvidx][keeptlidx] + (param->fd_np_st_mlt)[flvidx]*st_asc_forf
      +(param->fd_np_xi_mlt)[flvidx]*xi_asc_forf + (param->fd_np_stxi_mlt)[flvidx]*st_asc_forf*xi_asc_forf;
    
    constraint(&f_delta, 0.0, 1.0);
    constraint(&f_neg, 0.0, 1.0);

    *fd = (1. - f_neg) * f_delta;
    *fp = (*fd < 1. ? (1. - f_neg) * (1. - f_delta) / (1. - *fd): 0.5);
    *tau_np_p = param->Snp_global*((param->tau_np_p_mlt)[flvidx][0] + (param->tau_np_p_mlt)[flvidx][1]*st_asc +
				   (param->tau_np_p_xi_mlt)[flvidx]*xi_asc + (param->tau_np_p_stxi_mlt)[flvidx]*st_asc*xi_asc);
    *tau_np_n = param->Snp_global*((param->tau_np_n_mlt)[flvidx][0] + (param->tau_np_n_mlt)[flvidx][1]*st_asc +
				   (param->tau_np_n_xi_mlt)[flvidx]*xi_asc + (param->tau_np_n_stxi_mlt)[flvidx]*st_asc*xi_asc);
  }else{ /*  single track vertex */
    *fd = (param->fd_np_sgl)[flvidx][keeptlidx];
    *fp = (param->fp_np_sgl)[flvidx];
    //double f_delta = (param->fd_np_sgl)[flvidx][keeptlidx];
    //double f_neg = 1. - (param->fp_np_sgl)[flvidx];    
    //constraint(&f_delta, 0.0, 1.0);
    //constraint(&f_neg, 0.0, 1.0);
    //*fd = (1. - f_neg) * f_delta;
    //*fp = (*fd < 1. ? (1. - f_neg) * (1. - f_delta) / (1. - *fd): 0.5);
    *tau_np_p    = param->Snp_global*((param->tau_np_p_sgl)[flvidx][0] + (param->tau_np_p_sgl)[flvidx][1]*st_asc);
    *tau_np_n    = param->Snp_global*((param->tau_np_n_sgl)[flvidx][0] + (param->tau_np_n_sgl)[flvidx][1]*st_asc);
  }
#if defined(__FORCE_PARAM_CONSTRAINT__) &&  __FORCE_PARAM_CONSTRAINT__
  constraint(tau_np_p,  __POSTIVEPARAM_MIN__);
  constraint(tau_np_n,  __POSTIVEPARAM_MIN__);
#endif /* __FORCE_PARAM_CONSTRAINT__ */
  constraint(fp, 0.0, 1.0);
  constraint(fd, 0.0, 1.0);

  *tau_np_p_tl = *tau_np_p;
  *tau_np_n_tl = *tau_np_n;
  return;
}

void swap_rnp_param(double * fp,
		    double* tau_np_p, double* tau_np_n,
		    double* tau_np_p_tl, double* tau_np_n_tl){
  *fp = 1.0 - *fp;
  const double taunp_n_tmp    = *tau_np_n;
  *tau_np_n = *tau_np_p;
  *tau_np_p = taunp_n_tmp;
  const double taunp_n_tl_tmp = *tau_np_n_tl;
  *tau_np_n_tl = *tau_np_p_tl;
  *tau_np_p_tl = taunp_n_tl_tmp;
  return;
}


double Rnp(const double x,
	   const int flavor,
	   const int ntrk_asc, const double sz_asc,
	   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
	   const dtres_param_t * const param){

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;

  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double Li_d = DiracDelta(x);
  const double Li_p = Ep(x, tau_np_p);
  const double Li_n = En(x, tau_np_n);
  const double Li_e = fp*Li_p + (1.0 - fp) * Li_n;
  const double Li_m = fd*Li_d + (1.0 - fd) * Li_e;
  if(ntrk_asc>1){ /* multiple track case */
    return Li_m;
  }
  const double Li_p_tl = Ep(x, tau_np_p_tl);
  const double Li_n_tl = En(x, tau_np_n_tl);
  const double Li_e_tl = fp*Li_p_tl + (1.0 - fp) * Li_n_tl;
  const double Li_tl   = fd*Li_d + (1.0 - fd) * Li_e_tl;
  double Li  = (1.0 - ftail_asc) * Li_m + ftail_asc * Li_tl;
  return Li;
}

double RascRnp(const double x,
	       const int flavor,
	       const int ntrk_asc, const double sz_asc,
	       const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
	       const dtres_param_t * const param){

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;

  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double Li_md = gaussian(x, mu_main_asc, Smain_asc);
  const double Li_mp = Ep_conv_gauss(x, tau_np_p, mu_main_asc, Smain_asc);
  const double Li_mn = En_conv_gauss(x, tau_np_n, mu_main_asc, Smain_asc);
  const double Li_me = fp*Li_mp + (1.0 - fp) * Li_mn;
  const double Li_mt = fd*Li_md + (1.0 - fd) * Li_me;
  if(ftail_asc==0.0){ /* Mainly multiple track case */
    return Li_mt;
  }
  const double Li_td = gaussian(x, mu_tail_asc, Stail_asc);
  const double Li_tp = Ep_conv_gauss(x, tau_np_p_tl, mu_tail_asc, Stail_asc);
  const double Li_tn = En_conv_gauss(x, tau_np_n_tl, mu_tail_asc, Stail_asc);
  const double Li_te = fp*Li_tp + (1.0 - fp) * Li_tn;
  const double Li_tt = fd*Li_td + (1.0 - fd) * Li_te;
  double Li  = (1.0 - ftail_asc) * Li_mt + ftail_asc * Li_tt;
  return Li;
}

double Rdet(const double x,
	    const int ntrk_rec, const double sz_rec,
	    const double chisq_z_rec, const int ndf_z_rec,
	    const int ntrk_asc, const double sz_asc,
	    const double chisq_z_asc, const int ndf_z_asc,
	    const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = gaussian(x, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = gaussian(x, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = gaussian(x, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = gaussian(x, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = gaussian(x, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double RdetRnp_sup(const double x, 
			  const double fd, const double fp,
			  const double tau_np_p, const double tau_np_n,
			  const double mu_det, const double sigma_det){
  const double Li 
    = fd*gaussian(x, mu_det, sigma_det)
    + (1.0-fd)*(fp*Ep_conv_gauss(x, tau_np_p, mu_det, sigma_det)
		+(1.0-fp)*En_conv_gauss(x, tau_np_n, mu_det, sigma_det));
  return Li;
}

double RdetRnp(const double x,
	       const int flavor,
	       const int ntrk_rec, const double sz_rec,
	       const double chisq_z_rec, const int ndf_z_rec,
	       const int ntrk_asc, const double sz_asc,
	       const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
	       const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);
  
  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);
  
  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = RdetRnp_sup(x, fd, fp, tau_np_p, tau_np_n, mu_mm, sigma_mm);

  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = RdetRnp_sup(x, fd, fp, tau_np_p_tl, tau_np_n_tl, mu_mt, sigma_mt);

    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = RdetRnp_sup(x, fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);

      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = RdetRnp_sup(x, fd, fp, tau_np_p_tl, tau_np_n_tl, mu_tt, sigma_tt);

      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = RdetRnp_sup(x, fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double EfRk_fullrec(const double x, const double tau,
		    const double ak, const double ck){
  if(ak==0.0){
    if(debugout("INFO")) std::printf("[EfRk_fullrec] ak==% e where ak should not be zero.", ak);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double Li
    = (x<0.0) ? (0.5*(1-r_ckak)*En(x,tau*(ak-ck))): (0.5*(1+r_ckak)*Ep(x, tau*(ak+ck)));
  return Li;
}


double AfRk_fullrec(const double x,
		    const double tau, const double dm,
		    const double ak, const double ck){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  double Li = 0.0;
  if(x<0.0){
    const double ndm  = dm/(ak-ck);
    const double ntau = tau*(ak-ck);
    Li =  inv_ak*fact*(An(x, ntau, ndm)-ndmtau*Mn(x, ntau, ndm));
  }else{
    const double ndm  = dm/(ak+ck);
    const double ntau = tau*(ak+ck);
    Li =  inv_ak*fact*(Ap(x, ntau, ndm)-ndmtau*Mp(x, ntau, ndm));
  }
  return Li;
}

double MfRk_fullrec(const double x,
		    const double tau, const double dm,
		    const double ak, const double ck){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[MfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  double Li = 0.0;
  if(x<0.0){
    const double ndm  = dm/(ak-ck);
    const double ntau = tau*(ak-ck);
    Li =  inv_ak*fact*(Mn(x, ntau, ndm)+ndmtau*An(x, ntau, ndm));
  }else{
    const double ndm  = dm/(ak+ck);
    const double ntau = tau*(ak+ck);
    Li =  inv_ak*fact*(Mp(x, ntau, ndm)+ndmtau*Ap(x, ntau, ndm));
  }
#if 0
  const double cktau = ck*tau;
  if(ck>0.0){
    Li += fact*ndmtau*cktau*Ep(x, cktau);
  }else if (ck<0.0){
    Li += fact*ndmtau*cktau*En(x, -cktau);
  }
#endif
  return Li;
}

double EfRk_partial(const double x, const double tau,
		    const double dz, const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_EfEnp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fxEn, &fxEp,
		 tau, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_EfEnp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fxEn, &fxEp,
		 tau, tau_k_2_n, tau_k_2_p, fk2);

  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  En_conv_gauss(x, tau,       0.0, sigma_k);
  if(fEp!=0.0)    Li += fEp    *  Ep_conv_gauss(x, tau,       0.0, sigma_k);
  if(fEn_k1!=0.0) Li += fEn_k1 *  En_conv_gauss(x, tau_k_1_n, 0.0, sigma_k);
  if(fEp_k1!=0.0) Li += fEp_k1 *  Ep_conv_gauss(x, tau_k_1_p, 0.0, sigma_k);
  if(fEn_k2!=0.0) Li += fEn_k2 *  En_conv_gauss(x, tau_k_2_n, 0.0, sigma_k);
  if(fEp_k2!=0.0) Li += fEp_k2 *  Ep_conv_gauss(x, tau_k_2_p, 0.0, sigma_k);
  if(fxEn!=0.0)   Li += fxEn   * xEn_conv_gauss(x, tau,       0.0, sigma_k);
  if(fxEp!=0.0)   Li += fxEp   * xEp_conv_gauss(x, tau,       0.0, sigma_k);
  return Li;
}

double AfRk_partial(const double x, const double tau, const double dm,
		    const double dz, const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);
  
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * An_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fAp!=0.0)    Li += fAp    * Ap_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fMn!=0.0)    Li += fMn    * Mn_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fMp!=0.0)    Li += fMp    * Mp_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fEn_k1!=0.0) Li += fEn_k1 * En_conv_gauss(x, tau_k_1_n, 0.0, sigma_k);
  if(fEp_k1!=0.0) Li += fEp_k1 * Ep_conv_gauss(x, tau_k_1_p, 0.0, sigma_k);
  if(fEn_k2!=0.0) Li += fEn_k2 * En_conv_gauss(x, tau_k_2_n, 0.0, sigma_k);
  if(fEp_k2!=0.0) Li += fEp_k2 * Ep_conv_gauss(x, tau_k_2_p, 0.0, sigma_k);
  return Li;
}

double MfRk_partial(const double x, const double tau, const double dm,
		    const double dz, const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);
  
  double Li = 0.0;
  if(fMn!=0.0)    Li += fMn    * Mn_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fMp!=0.0)    Li += fMp    * Mp_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fAn!=0.0)    Li += fAn    * An_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fAp!=0.0)    Li += fAp    * Ap_conv_gauss(x, tau, dm,   0.0, sigma_k);
  if(fEn_k1!=0.0) Li += fEn_k1 * En_conv_gauss(x, tau_k_1_n, 0.0, sigma_k);
  if(fEp_k1!=0.0) Li += fEp_k1 * Ep_conv_gauss(x, tau_k_1_p, 0.0, sigma_k);
  if(fEn_k2!=0.0) Li += fEn_k2 * En_conv_gauss(x, tau_k_2_n, 0.0, sigma_k);
  if(fEp_k2!=0.0) Li += fEp_k2 * Ep_conv_gauss(x, tau_k_2_p, 0.0, sigma_k);
  return Li;
}

inline double EfRkRdet_full_sup(const double x,
				const double fact_n, const double ntau_n,
				const double fact_p, const double ntau_p,
				const double mu_det, const double sigma_det){
  const double Li
    = fact_n*En_conv_gauss(x, ntau_n, mu_det, sigma_det)
    + fact_p*Ep_conv_gauss(x, ntau_p, mu_det, sigma_det);
  return Li;
}

double EfRkRdet_fullrec(const double x, const double tau,
			const double ak, const double ck,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param){
  if(ak==0.0){
    if(debugout("INFO")) std::printf("[EfRk_fullrec] ak==% e where ak should not be zero.", ak);
    return -DBL_MAX;
  }

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double r_ckak = ck/ak;
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdet_full_sup(x,
					    fact_n, ntau_n, fact_p, ntau_p,
					    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdet_full_sup(x,
					      fact_n, ntau_n, fact_p, ntau_p,
					      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdet_full_sup(x,
						fact_n, ntau_n, fact_p, ntau_p,
						mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdet_full_sup(x,
						fact_n, ntau_n, fact_p, ntau_p,
						mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EfRkRdet_full_sup(x,
					      fact_n, ntau_n, fact_p, ntau_p,
					      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double AfRkRdet_full_sup(const double x,
				const double fact_n, const double ntau_n, const double ndm_n,
				const double fact_p, const double ntau_p, const double ndm_p, const double ndmtau,
				const double mu, const double sigma){
  const double Li
    = fact_n*(An_conv_gauss(x, ntau_n, ndm_n, mu, sigma)-ndmtau*Mn_conv_gauss(x, ntau_n, ndm_n, mu, sigma))
    + fact_p*(Ap_conv_gauss(x, ntau_p, ndm_p, mu, sigma)-ndmtau*Mp_conv_gauss(x, ntau_p, ndm_p, mu, sigma));
  return Li;
}

double AfRkRdet_fullrec(const double x,
			const double tau, const double dm,
			const double ak, const double ck,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = AfRkRdet_full_sup(x,
					    fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
					    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = AfRkRdet_full_sup(x,
					      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
					      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = AfRkRdet_full_sup(x,
						fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = AfRkRdet_full_sup(x,
						fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = AfRkRdet_full_sup(x,
					      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
					      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double MfRkRdet_full_sup(const double x,
				const double fact_n, const double ntau_n, const double ndm_n,
				const double fact_p, const double ntau_p, const double ndm_p,
				const double ndmtau, const double cktau,
				const double mu, const double sigma){
  double Li
    = fact_n*(Mn_conv_gauss(x, ntau_n, ndm_n, mu, sigma)+ndmtau*An_conv_gauss(x, ntau_n, ndm_n, mu, sigma))
    + fact_p*(Mp_conv_gauss(x, ntau_p, ndm_p, mu, sigma)+ndmtau*Ap_conv_gauss(x, ntau_p, ndm_p, mu, sigma));
#if 0
  if (cktau<0.0){
    Li += fact_n*ndmtau*cktau*En(x, -cktau);
  }else if(cktau>0.0){
    Li += fact_p*ndmtau*cktau*Ep(x, cktau);
  }
#endif
  return Li;
}

double MfRkRdet_fullrec(const double x,
			const double tau, const double dm,
			const double ak, const double ck,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[MfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = MfRkRdet_full_sup(x,
					    fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
					    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = MfRkRdet_full_sup(x,
					      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
					      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = MfRkRdet_full_sup(x,
						fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = MfRkRdet_full_sup(x,
						fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = MfRkRdet_full_sup(x,
					      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
					      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double EfRkRdet_partial_sup(const double x, const double tau,
				   const double tau_k_1_p, const double tau_k_2_p, 
				   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				   const double fEn, const double fEp, const double fxEn, const double fxEp,
				   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEn!=0.0)    Li += fEn    *  En_conv_gauss(x, tau,       mu_det, nsigma);
  if(fEp!=0.0)    Li += fEp    *  Ep_conv_gauss(x, tau,       mu_det, nsigma);
  if(fEn_k1!=0.0) Li += fEn_k1 *  En_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0) Li += fEp_k1 *  Ep_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0) Li += fEn_k2 *  En_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0) Li += fEp_k2 *  Ep_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fxEn!=0.0)   Li += fxEn   * xEn_conv_gauss(x, tau,       mu_det, nsigma);
  if(fxEp!=0.0)   Li += fxEp   * xEp_conv_gauss(x, tau,       mu_det, nsigma);
  return Li;
}

double EfRkRdet_partial(const double x, const double tau,
			const double dz,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_EfEnp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fxEn, &fxEp,
		 tau, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_EfEnp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fxEn, &fxEp,
		 tau, tau_k_2_n, tau_k_2_p, fk2);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdet_partial_sup(x, tau,
					       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
					       fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdet_partial_sup(x, tau,
						 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						 fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdet_partial_sup(x, tau,
						   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						   fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdet_partial_sup(x, tau,
						   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						   fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EfRkRdet_partial_sup(x, tau,
						 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						 fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double AfRkRdet_partial_sup(const double x, const double tau, const double dm,
				   const double tau_k_1_p, const double tau_k_2_p, 
				   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				   const double fAn, const double fAp, const double fMn, const double fMp,
				   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fAn!=0.0)    Li += fAn    * An_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fAp!=0.0)    Li += fAp    * Ap_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fMn!=0.0)    Li += fMn    * Mn_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fMp!=0.0)    Li += fMp    * Mp_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fEn_k1!=0.0) Li += fEn_k1 * En_conv_gauss(x, tau_k_1_n, mu_det,  nsigma);
  if(fEp_k1!=0.0) Li += fEp_k1 * Ep_conv_gauss(x, tau_k_1_p, mu_det,  nsigma);
  if(fEn_k2!=0.0) Li += fEn_k2 * En_conv_gauss(x, tau_k_2_n, mu_det,  nsigma);
  if(fEp_k2!=0.0) Li += fEp_k2 * Ep_conv_gauss(x, tau_k_2_p, mu_det,  nsigma);
  return Li;
}

double AfRkRdet_partial(const double x, const double tau, const double dm,
			const double dz,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);
 

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);
  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = AfRkRdet_partial_sup(x, tau, dm,
					       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
					       fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = AfRkRdet_partial_sup(x, tau, dm,
						 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						 fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = AfRkRdet_partial_sup(x, tau, dm,
						   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						   fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = AfRkRdet_partial_sup(x, tau, dm,
						   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						   fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = AfRkRdet_partial_sup(x, tau, dm,
						 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						 fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;

}

inline double MfRkRdet_partial_sup(const double x, const double tau, const double dm,
				   const double tau_k_1_p, const double tau_k_2_p, 
				   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				   const double fMn, const double fMp, const double fAn, const double fAp,
				   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fMn!=0.0)    Li += fMn    * Mn_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fMp!=0.0)    Li += fMp    * Mp_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fAn!=0.0)    Li += fAn    * An_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fAp!=0.0)    Li += fAp    * Ap_conv_gauss(x, tau, dm,   mu_det,  nsigma);
  if(fEn_k1!=0.0) Li += fEn_k1 * En_conv_gauss(x, tau_k_1_n, mu_det,  nsigma);
  if(fEp_k1!=0.0) Li += fEp_k1 * Ep_conv_gauss(x, tau_k_1_p, mu_det,  nsigma);
  if(fEn_k2!=0.0) Li += fEn_k2 * En_conv_gauss(x, tau_k_2_n, mu_det,  nsigma);
  if(fEp_k2!=0.0) Li += fEp_k2 * Ep_conv_gauss(x, tau_k_2_p, mu_det,  nsigma);
  return Li;
}

double MfRkRdet_partial(const double x, const double tau, const double dm,
			const double dz,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = MfRkRdet_partial_sup(x, tau, dm,
					       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
					       fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = MfRkRdet_partial_sup(x, tau, dm,
						 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						 fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = MfRkRdet_partial_sup(x, tau, dm,
						   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						   fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = MfRkRdet_partial_sup(x, tau, dm,
						   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						   fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = MfRkRdet_partial_sup(x, tau, dm,
						 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						 fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
  
}
//////
inline double EfRkRdetRnp_full_sup(const double x,
				   const double fact_n, const double ntau_n,
				   const double fact_p, const double ntau_p,
				   const double fd, const double fp, const double tau_np_p, const double tau_np_n,
				   const double mu_det, const double sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  double fEn = fact_n*fd;
  double fEp = fact_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0, fxEn = 0.0, fxEp = 0.0;
  add_EnEn_coef(&fEn, &fEn_np, &fxEn, ntau_n, tau_np_n, fact_n*nfn);
  add_EnEp_coef(&fEn, &fEp_np,        ntau_n, tau_np_p, fact_n*nfp);
  add_EpEn_coef(&fEp, &fEn_np,        ntau_p, tau_np_n, fact_p*nfn);
  add_EpEp_coef(&fEp, &fEp_np, &fxEp, ntau_p, tau_np_p, fact_p*nfp);
  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  En_conv_gauss(x, ntau_n,   mu_det, sigma_det);
  if(fEp!=0.0)    Li += fEp    *  Ep_conv_gauss(x, ntau_p,   mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  En_conv_gauss(x, tau_np_n, mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  Ep_conv_gauss(x, tau_np_p, mu_det, sigma_det);
  if(fxEn!=0.0)   Li += fxEn   * xEn_conv_gauss(x, ntau_n,   mu_det, sigma_det);
  if(fxEp!=0.0)   Li += fxEp   * xEp_conv_gauss(x, ntau_p,   mu_det, sigma_det);
  return Li;
}

double EfRkRdetRnp_fullrec(const double x, const int flavor,  const double tau,
			   const double ak, const double ck,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param){
  if(ak==0.0){
    if(debugout("INFO")) std::printf("[EfRk_fullrec] ak==% e where ak should not be zero.", ak);
    return -DBL_MAX;
  }

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double r_ckak = ck/ak;
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdetRnp_full_sup(x,
					       fact_n, ntau_n, fact_p, ntau_p,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdetRnp_full_sup(x,
						 fact_n, ntau_n, fact_p, ntau_p,
						 fd, fp, tau_np_p_tl, tau_np_n_tl,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdetRnp_full_sup(x,
						   fact_n, ntau_n, fact_p, ntau_p,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdetRnp_full_sup(x,
						   fact_n, ntau_n, fact_p, ntau_p,
						   fd, fp, tau_np_p_tl, tau_np_n_tl,
						   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EfRkRdetRnp_full_sup(x,
						 fact_n, ntau_n, fact_p, ntau_p,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double AfRkRdetRnp_full_sup(const double x,
				   const double fact_n, const double ntau_n, const double ndm_n,
				   const double fact_p, const double ntau_p, const double ndm_p, const double ndmtau,
				   const double fd, const double fp, const double tau_np_p, const double tau_np_n,
				   const double mu_det, const double sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_mn_n = -fact_n*ndmtau;
  const double w_mn_p = -fact_p*ndmtau;
  double fAn = fact_n*fd;
  double fAp = fact_p*fd;
  double fMn = w_mn_n*fd;
  double fMp = w_mn_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(&fAn, &fMn, &fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_AnEp_coef(&fAn, &fMn, &fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_MnEn_coef(&fMn, &fAn, &fEn_np, ntau_n, ndm_n, tau_np_n, w_mn_n*nfn);
  add_MnEp_coef(&fMn, &fAn, &fEp_np, ntau_n, ndm_n, tau_np_p, w_mn_n*nfp);
  add_ApEn_coef(&fAp, &fMp, &fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_ApEp_coef(&fAp, &fMp, &fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);
  add_MpEn_coef(&fMp, &fAp, &fEn_np, ntau_p, ndm_p, tau_np_n, w_mn_p*nfn);
  add_MpEp_coef(&fMp, &fAp, &fEp_np, ntau_p, ndm_p, tau_np_p, w_mn_p*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * An_conv_gauss(x, ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    * Ap_conv_gauss(x, ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    * Mn_conv_gauss(x, ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    * Mp_conv_gauss(x, ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np * En_conv_gauss(x, tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np * Ep_conv_gauss(x, tau_np_p,      mu_det, sigma_det);
  return Li;
}

double AfRkRdetRnp_fullrec(const double x, const int flavor, 
			   const double tau, const double dm,
			   const double ak, const double ck,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = AfRkRdetRnp_full_sup(x,
					       fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = AfRkRdetRnp_full_sup(x,
						 fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						 fd, fp, tau_np_p_tl, tau_np_n_tl,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = AfRkRdetRnp_full_sup(x,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = AfRkRdetRnp_full_sup(x,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						   fd, fp, tau_np_p_tl, tau_np_n_tl,
						   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = AfRkRdetRnp_full_sup(x,
						 fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double MfRkRdetRnp_full_sup(const double x,
				   const double fact_n, const double ntau_n, const double ndm_n,
				   const double fact_p, const double ntau_p, const double ndm_p,
				   const double ndmtau, const double cktau,
				   const double fd, const double fp, const double tau_np_p, const double tau_np_n,
				   const double mu_det, const double sigma_det){

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_an_n = fact_n*ndmtau;
  const double w_an_p = fact_p*ndmtau;
  double fMn = fact_n*fd;
  double fMp = fact_p*fd;
  double fAn = w_an_n*fd;
  double fAp = w_an_p*fd;
  double fEn_k = 0.0, fEp_k = 0.0, fxEn_k = 0.0, fxEp_k = 0.0;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(&fAn, &fMn, &fEn_np, ntau_n, ndm_n, tau_np_n, w_an_n*nfn);
  add_AnEp_coef(&fAn, &fMn, &fEp_np, ntau_n, ndm_n, tau_np_p, w_an_n*nfp);
  add_MnEn_coef(&fMn, &fAn, &fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_MnEp_coef(&fMn, &fAn, &fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_ApEn_coef(&fAp, &fMp, &fEn_np, ntau_p, ndm_p, tau_np_n, w_an_p*nfn);
  add_ApEp_coef(&fAp, &fMp, &fEp_np, ntau_p, ndm_p, tau_np_p, w_an_p*nfp);
  add_MpEn_coef(&fMp, &fAp, &fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_MpEp_coef(&fMp, &fAp, &fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);
#if 0
  if(cktau<0.0){
    const double w_ex_n = fact_n*ndmtau*cktau;
    add_EnEn_coef(&fEn_k, &fEn_np, &fxEn_k, -cktau, tau_np_n, w_ex_n*nfn);
    add_EnEp_coef(&fEn_k, &fEp_np,          -cktau, tau_np_p, w_ex_n*nfp);
  }else if(cktau>0.0){
    const double w_ex_p = fact_p*ndmtau*cktau;
    add_EpEn_coef(&fEp_k, &fEn_np,          cktau, tau_np_n, w_ex_p*nfn);
    add_EpEp_coef(&fEp_k, &fEp_np, &fxEp_k, cktau, tau_np_p, w_ex_p*nfp);
  }
#endif  
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    *  An_conv_gauss(x,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    *  Ap_conv_gauss(x,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    *  Mn_conv_gauss(x,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    *  Mp_conv_gauss(x,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_k!=0.0)  Li += fEn_k  *  En_conv_gauss(x, -cktau,         mu_det, sigma_det);
  if(fEp_k!=0.0)  Li += fEp_k  *  Ep_conv_gauss(x,  cktau,         mu_det, sigma_det);
  if(fxEn_k!=0.0) Li += fxEn_k * xEn_conv_gauss(x, -cktau,         mu_det, sigma_det);
  if(fxEp_k!=0.0) Li += fxEp_k * xEp_conv_gauss(x,  cktau,         mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  En_conv_gauss(x,  tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  Ep_conv_gauss(x,  tau_np_p,      mu_det, sigma_det);
  return Li;
}

double MfRkRdetRnp_fullrec(const double x, const int flavor, 
			   const double tau, const double dm,
			   const double ak, const double ck,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[MfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = MfRkRdetRnp_full_sup(x,
					       fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = MfRkRdetRnp_full_sup(x,
						 fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						 fd, fp, tau_np_p_tl, tau_np_n_tl,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = MfRkRdetRnp_full_sup(x,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = MfRkRdetRnp_full_sup(x,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						   fd, fp, tau_np_p_tl, tau_np_n_tl,
						   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = MfRkRdetRnp_full_sup(x,
						 fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double EfRkRdetRnp_partial_sup(const double x, const double tau,
				      const double tau_k_1_p, const double tau_k_2_p, 
				      const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				      const double tau_np_p, const double tau_np_n,
				      const double fEn, const double fEp,
				      const double fxEn, const double fxEp, const double fxxEn, const double fxxEp,
				      const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				      const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
				      const double fEn_np, const double fEp_np, 
				      const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEn!=0.0)     Li += fEn     *   En_conv_gauss(x, tau,       mu_det, nsigma);
  if(fEp!=0.0)     Li += fEp     *   Ep_conv_gauss(x, tau,       mu_det, nsigma);
  if(fxEn!=0.0)    Li += fxEn    *  xEn_conv_gauss(x, tau,       mu_det, nsigma);
  if(fxEp!=0.0)    Li += fxEp    *  xEp_conv_gauss(x, tau,       mu_det, nsigma);
  if(fxxEn!=0.0)   Li += fxxEn   * xxEn_conv_gauss(x, tau,       mu_det, nsigma);
  if(fxxEp!=0.0)   Li += fxxEp   * xxEp_conv_gauss(x, tau,       mu_det, nsigma);
  if(fEn_k1!=0.0)  Li += fEn_k1  *   En_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)  Li += fEp_k1  *   Ep_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)  Li += fEn_k2  *   En_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)  Li += fEp_k2  *   Ep_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0) Li += fxEn_k1 *  xEn_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0) Li += fxEp_k1 *  xEp_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0) Li += fxEn_k2 *  xEn_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0) Li += fxEp_k2 *  xEp_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)  Li += fEn_np  *   En_conv_gauss(x, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)  Li += fEp_np  *   Ep_conv_gauss(x, tau_np_p,  mu_det, nsigma);
  return Li;
}

double EfRkRdetRnp_partial(const double x, const int flavor,  const double tau,
			   const double dz,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;
  
//  dout(Debugout::INFO,"dt_resolution")<<"===== EfRkRdetRnp_partial ====="<<std::endl;
//  dout(Debugout::INFO,"dt_resolution")<<"x: "<<x<<" flavor: "<<flavor<<" tau: "<<tau<<" dz: "<<dz<<std::endl;
//  dout(Debugout::INFO,"dt_resolution")<<"ntrk_rec: "<<ntrk_rec<<" sz_rec: "<<sz_rec<<" chisq_z_rec: "<<chisq_z_rec<<" ndf_z_rec: "<<ndf_z_rec<<std::endl;
//  dout(Debugout::INFO,"dt_resolution")<<"ntrk_asc: "<<ntrk_asc<<" sz_asc: "<<sz_asc<<" chisq_z_asc: "<<chisq_z_asc<<" ndf_z_asc: "<<ndf_z_asc<<std::endl;
//  dout(Debugout::INFO,"dt_resolution")<<"fd: "<<fd<<" fp: "<<fp<<" tau_np_p: "<<tau_np_p<<" tau_np_n: "
//    <<tau_np_n<<" tau_np_p_tl: "<<tau_np_p_tl<<" tau_np_n_tl: "<<tau_np_n_tl<<std::endl;
	  
  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0, fxxEn=0.0, fxxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;

  add_EfEnp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fxEn, &fxEp,
		 tau, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_EfEnp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fxEn, &fxEp,
		 tau, tau_k_2_n, tau_k_2_p, fk2*fd);
  
  add_EfEnpEn_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fEn_np,
		   &fxEn, &fxEp, &fxEn_k1, &fxxEn,
		   tau, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_EfEnpEn_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fEn_np,
		   &fxEn, &fxEp, &fxEn_k2, &fxxEn,
		   tau, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_EfEnpEp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fEp_np,
		   &fxEn, &fxEp, &fxEp_k1, &fxxEp,
		   tau, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_EfEnpEp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fEp_np,
		   &fxEn, &fxEp, &fxEp_k2, &fxxEp,
		   tau, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdetRnp_partial_sup(x, tau,
						  tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						  tau_np_p, tau_np_n,
						  fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						  fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						  fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						  fEn_np, fEp_np, 
						  mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fEn_tl=0.0, fEp_tl=0.0, fxEn_tl=0.0, fxEp_tl=0.0, fxxEn_tl=0.0, fxxEp_tl=0.0; 
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
    
    add_EfEnp_coef(&fEn_tl, &fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fxEn_tl, &fxEp_tl,
		   tau, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_EfEnp_coef(&fEn_tl, &fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fxEn_tl, &fxEp_tl,
		   tau, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_EfEnpEn_coef(&fEn_tl, &fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEn_k1_tl, &fxxEn_tl,
		     tau, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_EfEnpEn_coef(&fEn_tl, &fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEn_k2_tl, &fxxEn_tl,
		     tau, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_EfEnpEp_coef(&fEn_tl, &fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEp_k1_tl, &fxxEp_tl,
		     tau, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_EfEnpEp_coef(&fEn_tl, &fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEp_k2_tl, &fxxEp_tl,
		     tau, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
      
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdetRnp_partial_sup(x, tau,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p_tl, tau_np_n_tl,
						    fEn_tl, fEp_tl, fxEn_tl, fxEp_tl, fxxEn_tl, fxxEp_tl,
						    fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
						    fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl,
						    fEn_np_tl, fEp_np_tl, 
						    mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/

      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdetRnp_partial_sup(x, tau,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p,  tau_np_n, 
						      fEn,  fEp,  fxEn,  fxEp,  fxxEn,  fxxEp, 
						      fEn_k1,  fEp_k1,  fEn_k2,  fEp_k2, 
						      fxEn_k1,  fxEp_k1,  fxEn_k2,  fxEp_k2, 
						      fEn_np,  fEp_np,  
						      mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdetRnp_partial_sup(x, tau,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p_tl,  tau_np_n_tl, 
						      fEn_tl,  fEp_tl,  fxEn_tl,  fxEp_tl,  fxxEn_tl,  fxxEp_tl, 
						      fEn_k1_tl,  fEp_k1_tl,  fEn_k2_tl,  fEp_k2_tl, 
						      fxEn_k1_tl,  fxEp_k1_tl,  fxEn_k2_tl,  fxEp_k2_tl, 
						      fEn_np_tl,  fEp_np_tl,  
						      mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EfRkRdetRnp_partial_sup(x, tau,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p, tau_np_n,
						    fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						    fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						    fEn_np, fEp_np, 
						    mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double AfRkRdetRnp_partial_sup(const double x, const double tau, const double dm,
				      const double tau_k_1_p, const double tau_k_2_p, 
				      const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				      const double tau_np_p, const double tau_np_n,
				      const double fAn, const double fAp, const double fMn, const double fMp,
				      const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				      const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2, 
				      const double fEn_np, const double fEp_np,
				      const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fAn!=0.0)      Li += fAn     *  An_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fAp!=0.0)      Li += fAp     *  Ap_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fMn!=0.0)      Li += fMn     *  Mn_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fMp!=0.0)      Li += fMp     *  Mp_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fEn_k1!=0.0)   Li += fEn_k1  *  En_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)   Li += fEp_k1  *  Ep_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)   Li += fEn_k2  *  En_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)   Li += fEp_k2  *  Ep_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0)  Li += fxEn_k1 * xEn_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0)  Li += fxEp_k1 * xEp_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0)  Li += fxEn_k2 * xEn_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0)  Li += fxEp_k2 * xEp_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)   Li += fEn_np  *  En_conv_gauss(x, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)   Li += fEp_np  *  Ep_conv_gauss(x, tau_np_p,  mu_det, nsigma);
  return Li;
}

double AfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dm,
			   const double dz,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param){
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);
 
  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;
  
  //dout(Debugout::INFO,"dt_resolution")<<"AfRkRdetRnp_partial()>> fk2: "<<fk2<<" fk1: "<<fk1<<std::endl;
  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1,
  		 tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
  // tau, dm, tau_k_1_n, tau_k_1_p, fk2*fd); // BUGBUGBUG	 

  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);

  add_AfEnpEn_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1, &fEn_np, &fxEn_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_AfEnpEn_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2, &fEn_np, &fxEn_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_AfEnpEp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1, &fEp_np, &fxEp_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_AfEnpEp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2, &fEp_np, &fxEp_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = AfRkRdetRnp_partial_sup(x, tau, dm,
						  tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						  tau_np_p, tau_np_n,
						  fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						  fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						  mu_mm, sigma_mm);
  if(ftail_asc>0.0){

    double fAn_tl=0.0, fAp_tl=0.0, fMn_tl=0.0, fMp_tl=0.0;
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
  
    add_AfEnp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k1_tl, &fEp_k1_tl,
		   tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
    //tau, dm, tau_k_1_n, tau_k_1_p, fk2*fd);//	BUGBUGBUG  

    add_AfEnp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k2_tl, &fEp_k2_tl,
		   tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_AfEnpEn_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl, &fxEn_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_AfEnpEn_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl, &fxEn_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_AfEnpEp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl, &fxEp_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_AfEnpEp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl, &fxEp_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);

    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = AfRkRdetRnp_partial_sup(x, tau, dm,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p_tl, tau_np_n_tl,
						    fAn_tl, fAp_tl, fMn_tl, fMp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
						    fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
						    mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = AfRkRdetRnp_partial_sup(x, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p, tau_np_n,
						      fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						      mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = AfRkRdetRnp_partial_sup(x, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p_tl, tau_np_n_tl,
						      fAn_tl, fAp_tl, fMn_tl, fMp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
						      fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
						      mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = AfRkRdetRnp_partial_sup(x, tau, dm,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p, tau_np_n,
						    fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						    mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;

}

inline double MfRkRdetRnp_partial_sup(const double x, const double tau, const double dm,
				      const double tau_k_1_p, const double tau_k_2_p, 
				      const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				      const double tau_np_p, const double tau_np_n,
				      const double fMn, const double fMp, const double fAn, const double fAp,
				      const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				      const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
				      const double fEn_np, const double fEp_np,
				      const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fMn!=0.0)      Li += fMn     *  Mn_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fMp!=0.0)      Li += fMp     *  Mp_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fAn!=0.0)      Li += fAn     *  An_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fAp!=0.0)      Li += fAp     *  Ap_conv_gauss(x, tau,   dm, mu_det, nsigma);
  if(fEn_k1!=0.0)   Li += fEn_k1  *  En_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)   Li += fEp_k1  *  Ep_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)   Li += fEn_k2  *  En_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)   Li += fEp_k2  *  Ep_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0)  Li += fxEn_k1 * xEn_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0)  Li += fxEp_k1 * xEp_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0)  Li += fxEn_k2 * xEn_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0)  Li += fxEp_k2 * xEp_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)   Li += fEn_np  *  En_conv_gauss(x, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)   Li += fEp_np  *  Ep_conv_gauss(x, tau_np_p,  mu_det, nsigma);
  return Li;
}

double MfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dm,
			   const double dz,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;
  
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);

  add_MfEnpEn_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1, &fEn_np, &fxEn_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_MfEnpEn_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2, &fEn_np, &fxEn_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_MfEnpEp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1, &fEp_np, &fxEp_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_MfEnpEp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2, &fEp_np, &fxEp_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = MfRkRdetRnp_partial_sup(x, tau, dm,
						  tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						  tau_np_p, tau_np_n,
						  fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						  fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						  mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fAn_tl=0.0, fAp_tl=0.0, fMn_tl=0.0, fMp_tl=0.0;
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
  
    add_MfEnp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k1_tl, &fEp_k1_tl,
		   tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_MfEnp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k2_tl, &fEp_k2_tl,
		   tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_MfEnpEn_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl, &fxEn_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_MfEnpEn_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl, &fxEn_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_MfEnpEp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl, &fxEp_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_MfEnpEp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl, &fxEp_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = MfRkRdetRnp_partial_sup(x, tau, dm,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p_tl, tau_np_n_tl,
						    fMn_tl, fMp_tl, fAn_tl, fAp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
						    fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
						    mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = MfRkRdetRnp_partial_sup(x, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p, tau_np_n,
						      fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						      mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = MfRkRdetRnp_partial_sup(x, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p_tl, tau_np_n_tl,
						      fMn_tl, fMp_tl, fAn_tl, fAp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
						      fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
						      mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = MfRkRdetRnp_partial_sup(x, tau, dm,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p, tau_np_n,
						    fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						    mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
  
}

double norm_Rk_partial(const double ll, const double ul, const double dz,
		       const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double Li_1 = norm_Enp_conv_gauss(ll, ul, tau_k_1_n, tau_k_1_p, 0.0, sigma_k);
  const double Li_2 = norm_Enp_conv_gauss(ll, ul, tau_k_2_n, tau_k_2_p, 0.0, sigma_k);
  const double Li = (1.0-fk2)*Li_1 + fk2*Li_2;
  return Li;
}

double norm_Rasc(const double ll, const double ul,
		 const int ntrk_asc, const double sz_asc,
		 const double chisq_z_asc, const int ndf_z_asc,
		 const dtres_param_t * const param){

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double Li_mn = norm_gaussian(ll, ul, mu_main_asc, Smain_asc);
  if(ftail_asc>0.0){ /* Mainly multiple track case */
    const double Li_tl = norm_gaussian(ll, ul, mu_tail_asc, Stail_asc);
    const double Li    = (1.0-ftail_asc)*Li_mn + ftail_asc*Li_tl;
    return Li;
  }
  return Li_mn;
}

double norm_Rrec(const double ll, const double ul,
		 const int ntrk_rec, const double sz_rec,
		 const double chisq_z_rec, const int ndf_z_rec,
		 const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  const double Li_mn = norm_gaussian(ll, ul, mu_main_rec, Smain_rec);
  if(ftail_rec>0.0){ /* Mainly multiple track case */
    const double Li_tl = norm_gaussian(ll, ul, mu_tail_rec, Stail_rec);
    const double Li    = (1.0-ftail_rec)*Li_mn + ftail_rec*Li_tl;
    return Li;
  }
  return Li_mn;
}

double norm_Rnp(const double ll, const double ul,
		const int flavor,
		const int ntrk_asc, const double sz_asc,
		const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
		const dtres_param_t * const param){

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;

  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double Li_d = (ll*ul<0.0) ? 1.0 : 0.0;
  const double Li_p = norm_Ep(ll, ul, tau_np_p);
  const double Li_n = norm_En(ll, ul, tau_np_n);
  const double Li_e = fp*Li_p + (1.0 - fp) * Li_n;
  const double Li_m = fd*Li_d + (1.0 - fd) * Li_e;
  if(ntrk_asc>1){ /* multiple track case */
    return Li_m;
  }
  const double Li_p_tl = norm_Ep(ll, ul, tau_np_p_tl);
  const double Li_n_tl = norm_En(ll, ul, tau_np_n_tl);
  const double Li_e_tl = fp*Li_p_tl + (1.0 - fp) * Li_n_tl;
  const double Li_tl   = fd*Li_d + (1.0 - fd) * Li_e_tl;
  double Li  = (1.0 - ftail_asc) * Li_m + ftail_asc * Li_tl;
  return Li;
}

double norm_RascRnp(const double ll, const double ul,
		    const int flavor,
		    const int ntrk_asc, const double sz_asc,
		    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
		    const dtres_param_t * const param){

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;

  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double Li_md = norm_gaussian(ll, ul, mu_main_asc, Smain_asc);
  const double Li_mp = norm_Ep_conv_gauss(ll, ul, tau_np_p, mu_main_asc, Smain_asc);
  const double Li_mn = norm_En_conv_gauss(ll, ul, tau_np_n, mu_main_asc, Smain_asc);
  const double Li_me = fp*Li_mp + (1.0 - fp) * Li_mn;
  const double Li_mt = fd*Li_md + (1.0 - fd) * Li_me;
  if(ftail_asc==0.0){ /* Mainly multiple track case */
    return Li_mt;
  }
  const double Li_td = norm_gaussian(ll, ul, mu_tail_asc, Stail_asc);
  const double Li_tp = norm_Ep_conv_gauss(ll, ul, tau_np_p_tl, mu_tail_asc, Stail_asc);
  const double Li_tn = norm_En_conv_gauss(ll, ul, tau_np_n_tl, mu_tail_asc, Stail_asc);
  const double Li_te = fp*Li_tp + (1.0 - fp) * Li_tn;
  const double Li_tt = fd*Li_td + (1.0 - fd) * Li_te;
  double Li  = (1.0 - ftail_asc) * Li_mt + ftail_asc * Li_tt;
  return Li;
}

double norm_Rdet(const double ll, const double ul,
		 const int ntrk_rec, const double sz_rec,
		 const double chisq_z_rec, const int ndf_z_rec,
		 const int ntrk_asc, const double sz_asc,
		 const double chisq_z_asc, const int ndf_z_asc,
		 const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_gaussian(ll, ul, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_gaussian(ll, ul, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_gaussian(ll, ul, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_gaussian(ll, ul, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_gaussian(ll, ul, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_RdetRnp_sup(const double ll, const double ul, 
			       const double fd, const double fp,
			       const double tau_np_p, const double tau_np_n,
			       const double mu_det, const double sigma_det){
  const double Li 
    = fd*norm_gaussian(ll, ul, mu_det, sigma_det)
    + (1.0-fd)*(fp*norm_Ep_conv_gauss(ll, ul, tau_np_p, mu_det, sigma_det)
		+(1.0-fp)*norm_En_conv_gauss(ll, ul, tau_np_n, mu_det, sigma_det));
  return Li;
}

double norm_RdetRnp(const double ll, const double ul,
		    const int flavor,
		    const int ntrk_rec, const double sz_rec,
		    const double chisq_z_rec, const int ndf_z_rec,
		    const int ntrk_asc, const double sz_asc,
		    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
		    const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);
  
  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);
  
  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_RdetRnp_sup(ll, ul, fd, fp, tau_np_p, tau_np_n, mu_mm, sigma_mm);

  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_RdetRnp_sup(ll, ul, fd, fp, tau_np_p_tl, tau_np_n_tl, mu_mt, sigma_mt);

    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_RdetRnp_sup(ll, ul, fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);

      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_RdetRnp_sup(ll, ul, fd, fp, tau_np_p_tl, tau_np_n_tl, mu_tt, sigma_tt);

      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_RdetRnp_sup(ll, ul, fd, fp, tau_np_p, tau_np_n, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double norm_EfRk_fullrec(const double ll, const double ul, const double tau,
			 const double ak, const double ck){
  if(ak==0.0){
    if(debugout("INFO")) std::printf("[EfRk_fullrec] ak==% e where ak should not be zero.", ak);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  double Li = 0.0;
  const double nll = (ll<ul) ? ll : ul;
  const double nul = (ll<ul) ? ul : ll;

  if(nll<0.0){
    if(nul<0.0){
      Li = 0.5*(1-r_ckak)*norm_En(nll, nul,tau*(ak-ck));
    }else{
      Li = 0.5*(1-r_ckak)*norm_En(nll, 0.0,tau*(ak-ck));
      Li += 0.5*(1+r_ckak)*norm_Ep(0.0, nul, tau*(ak+ck));
    }
  }else{
    Li = 0.5*(1+r_ckak)*norm_Ep(nll, nul, tau*(ak+ck));
  }
  return (ll<ul) ? Li : -Li;
}


double norm_AfRk_fullrec(const double ll, const double ul,
			 const double tau, const double dm,
			 const double ak, const double ck){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  double Li = 0.0;

  const double nll = (ll<ul) ? ll : ul;
  const double nul = (ll<ul) ? ul : ll;

  if(nll<0.0){
    const double ndm_n  = dm/(ak-ck);
    const double ntau_n = tau*(ak-ck);
    if(nul<0.0){
      Li = inv_ak*fact*(norm_An(nll, nul, ntau_n, ndm_n)-ndmtau*norm_Mn(nll, nul, ntau_n, ndm_n));
    }else{
      const double ndm_p  = dm/(ak+ck);
      const double ntau_p = tau*(ak+ck);
      Li =  inv_ak*fact*(norm_An(nll, 0.0, ntau_n, ndm_n)-ndmtau*norm_Mn(nll, 0.0, ntau_n, ndm_n));
      Li += inv_ak*fact*(norm_Ap(0.0, nul, ntau_p, ndm_p)-ndmtau*norm_Mp(0.0, nul, ntau_p, ndm_p));
    }
  }else{
    const double ndm_p  = dm/(ak+ck);
    const double ntau_p = tau*(ak+ck);
    Li = inv_ak*fact*(norm_Ap(nll, nul, ntau_p, ndm_p)-ndmtau*norm_Mp(nll, nul, ntau_p, ndm_p));
  }  
  return (ll<ul) ? Li : -Li;
}

double norm_MfRk_fullrec(const double ll, const double ul,
			 const double tau, const double dm,
			 const double ak, const double ck){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[MfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double fact = 1.0/(1.0+ndmtau*ndmtau);
  double Li = 0.0;
  const double nll = (ll<ul) ? ll : ul;
  const double nul = (ll<ul) ? ul : ll;

  if(nll<0.0){
    const double ndm_n  = dm/(ak-ck);
    const double ntau_n = tau*(ak-ck);
    if(nul<0.0){
      Li = inv_ak*fact*(norm_Mn(nll, nul, ntau_n, ndm_n)+ndmtau*norm_An(nll, nul, ntau_n, ndm_n));
    }else{
      const double ndm_p  = dm/(ak+ck);
      const double ntau_p = tau*(ak+ck);
      Li =  inv_ak*fact*(norm_Mn(nll, 0.0, ntau_n, ndm_n)+ndmtau*norm_An(nll, 0.0, ntau_n, ndm_n));
      Li += inv_ak*fact*(norm_Mp(0.0, nul, ntau_p, ndm_p)+ndmtau*norm_Ap(0.0, nul, ntau_p, ndm_p));
    }
  }else{
    const double ndm_p  = dm/(ak+ck);
    const double ntau_p = tau*(ak+ck);
    Li = inv_ak*fact*(norm_Mp(nll, nul, ntau_p, ndm_p)+ndmtau*norm_Ap(nll, nul, ntau_p, ndm_p));
  }
#if 0
  const double cktau = ck*tau;
  if(ck>0.0){
    Li += 2.0*fact*ndmtau*cktau*norm_Ep(nll, nul, cktau);
  }else if (ck<0.0){
    Li += 2.0*fact*ndmtau*cktau*norm_En(nll, nul, -cktau);
  }
#endif
  return (ll<ul) ? Li : -Li;
}

double norm_EfRk_partial(const double ll, const double ul, const double tau,
			 const double dz, const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_EfEnp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fxEn, &fxEp,
		 tau, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_EfEnp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fxEn, &fxEp,
		 tau, tau_k_2_n, tau_k_2_p, fk2);

  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  norm_En_conv_gauss(ll, ul, tau,       0.0, sigma_k);
  if(fEp!=0.0)    Li += fEp    *  norm_Ep_conv_gauss(ll, ul, tau,       0.0, sigma_k);
  if(fEn_k1!=0.0) Li += fEn_k1 *  norm_En_conv_gauss(ll, ul, tau_k_1_n, 0.0, sigma_k);
  if(fEp_k1!=0.0) Li += fEp_k1 *  norm_Ep_conv_gauss(ll, ul, tau_k_1_p, 0.0, sigma_k);
  if(fEn_k2!=0.0) Li += fEn_k2 *  norm_En_conv_gauss(ll, ul, tau_k_2_n, 0.0, sigma_k);
  if(fEp_k2!=0.0) Li += fEp_k2 *  norm_Ep_conv_gauss(ll, ul, tau_k_2_p, 0.0, sigma_k);
  if(fxEn!=0.0)   Li += fxEn   * norm_xEn_conv_gauss(ll, ul, tau,       0.0, sigma_k);
  if(fxEp!=0.0)   Li += fxEp   * norm_xEp_conv_gauss(ll, ul, tau,       0.0, sigma_k);
  return Li;
}

double norm_AfRk_partial(const double ll, const double ul, const double tau, const double dm,
			 const double dz, const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);
  
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * norm_An_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fAp!=0.0)    Li += fAp    * norm_Ap_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fMn!=0.0)    Li += fMn    * norm_Mn_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fMp!=0.0)    Li += fMp    * norm_Mp_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fEn_k1!=0.0) Li += fEn_k1 * norm_En_conv_gauss(ll, ul, tau_k_1_n, 0.0, sigma_k);
  if(fEp_k1!=0.0) Li += fEp_k1 * norm_Ep_conv_gauss(ll, ul, tau_k_1_p, 0.0, sigma_k);
  if(fEn_k2!=0.0) Li += fEn_k2 * norm_En_conv_gauss(ll, ul, tau_k_2_n, 0.0, sigma_k);
  if(fEp_k2!=0.0) Li += fEp_k2 * norm_Ep_conv_gauss(ll, ul, tau_k_2_p, 0.0, sigma_k);
  return Li;
}

double norm_MfRk_partial(const double ll, const double ul, const double tau, const double dm,
			 const double dz, const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);
  
  double Li = 0.0;
  if(fMn!=0.0)    Li += fMn    * norm_Mn_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fMp!=0.0)    Li += fMp    * norm_Mp_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fAn!=0.0)    Li += fAn    * norm_An_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fAp!=0.0)    Li += fAp    * norm_Ap_conv_gauss(ll, ul, tau, dm,   0.0, sigma_k);
  if(fEn_k1!=0.0) Li += fEn_k1 * norm_En_conv_gauss(ll, ul, tau_k_1_n, 0.0, sigma_k);
  if(fEp_k1!=0.0) Li += fEp_k1 * norm_Ep_conv_gauss(ll, ul, tau_k_1_p, 0.0, sigma_k);
  if(fEn_k2!=0.0) Li += fEn_k2 * norm_En_conv_gauss(ll, ul, tau_k_2_n, 0.0, sigma_k);
  if(fEp_k2!=0.0) Li += fEp_k2 * norm_Ep_conv_gauss(ll, ul, tau_k_2_p, 0.0, sigma_k);
  return Li;
}

inline double norm_EfRkRdet_full_sup(const double ll, const double ul,
				     const double fact_n, const double ntau_n,
				     const double fact_p, const double ntau_p,
				     const double mu_det, const double sigma_det){
  const double Li
    = fact_n*norm_En_conv_gauss(ll, ul, ntau_n, mu_det, sigma_det)
    + fact_p*norm_Ep_conv_gauss(ll, ul, ntau_p, mu_det, sigma_det);
  return Li;
}

double norm_EfRkRdet_fullrec(const double ll, const double ul, const double tau,
			     const double ak, const double ck,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param){
  if(ak==0.0){
    if(debugout("INFO")) std::printf("[EfRk_fullrec] ak==% e where ak should not be zero.", ak);
    return -DBL_MAX;
  }

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double r_ckak = ck/ak;
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdet_full_sup(ll, ul,
						 fact_n, ntau_n, fact_p, ntau_p,
						 mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdet_full_sup(ll, ul,
						   fact_n, ntau_n, fact_p, ntau_p,
						   mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EfRkRdet_full_sup(ll, ul,
						     fact_n, ntau_n, fact_p, ntau_p,
						     mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EfRkRdet_full_sup(ll, ul,
						     fact_n, ntau_n, fact_p, ntau_p,
						     mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EfRkRdet_full_sup(ll, ul,
						   fact_n, ntau_n, fact_p, ntau_p,
						   mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_AfRkRdet_full_sup(const double ll, const double ul,
				     const double fact_n, const double ntau_n, const double ndm_n,
				     const double fact_p, const double ntau_p, const double ndm_p, const double ndmtau,
				     const double mu, const double sigma){
  const double Li
    = fact_n*(norm_An_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma)-ndmtau*norm_Mn_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma))
    + fact_p*(norm_Ap_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma)-ndmtau*norm_Mp_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma));
  return Li;
}

double norm_AfRkRdet_fullrec(const double ll, const double ul,
			     const double tau, const double dm,
			     const double ak, const double ck,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_AfRkRdet_full_sup(ll, ul,
						 fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						 mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_AfRkRdet_full_sup(ll, ul,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						   mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_AfRkRdet_full_sup(ll, ul,
						     fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						     mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_AfRkRdet_full_sup(ll, ul,
						     fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						     mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_AfRkRdet_full_sup(ll, ul,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						   mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_MfRkRdet_full_sup(const double ll, const double ul,
				     const double fact_n, const double ntau_n, const double ndm_n,
				     const double fact_p, const double ntau_p, const double ndm_p,
				     const double ndmtau, const double cktau,
				     const double mu, const double sigma){
  double Li
    = fact_n*(norm_Mn_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma)+ndmtau*norm_An_conv_gauss(ll, ul, ntau_n, ndm_n, mu, sigma))
    + fact_p*(norm_Mp_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma)+ndmtau*norm_Ap_conv_gauss(ll, ul, ntau_p, ndm_p, mu, sigma));
#if 0
  if (cktau<0.0){
    Li += fact_n*ndmtau*cktau*norm_En(ll, ul, -cktau);
  }else if(cktau>0.0){
    Li += fact_p*ndmtau*cktau*norm_Ep(ll, ul, cktau);
  }
#endif
  return Li;
}

double norm_MfRkRdet_fullrec(const double ll, const double ul,
			     const double tau, const double dm,
			     const double ak, const double ck,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[MfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_MfRkRdet_full_sup(ll, ul,
						 fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						 mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_MfRkRdet_full_sup(ll, ul,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						   mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_MfRkRdet_full_sup(ll, ul,
						     fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						     mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_MfRkRdet_full_sup(ll, ul,
						     fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						     mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_MfRkRdet_full_sup(ll, ul,
						   fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						   mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_EfRkRdet_partial_sup(const double ll, const double ul, const double tau,
					const double tau_k_1_p, const double tau_k_2_p, 
					const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					const double fEn, const double fEp, const double fxEn, const double fxEp,
					const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEn!=0.0)    Li += fEn    *  norm_En_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fEp!=0.0)    Li += fEp    *  norm_Ep_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fEn_k1!=0.0) Li += fEn_k1 *  norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0) Li += fEp_k1 *  norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0) Li += fEn_k2 *  norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0) Li += fEp_k2 *  norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fxEn!=0.0)   Li += fxEn   * norm_xEn_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fxEp!=0.0)   Li += fxEp   * norm_xEp_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  return Li;
}

double norm_EfRkRdet_partial(const double ll, const double ul, const double tau,
			     const double dz,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_EfEnp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fxEn, &fxEp,
		 tau, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_EfEnp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fxEn, &fxEp,
		 tau, tau_k_2_n, tau_k_2_p, fk2);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdet_partial_sup(ll, ul, tau,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdet_partial_sup(ll, ul, tau,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EfRkRdet_partial_sup(ll, ul, tau,
							tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EfRkRdet_partial_sup(ll, ul, tau,
							tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EfRkRdet_partial_sup(ll, ul, tau,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      fEn, fEp, fxEn, fxEp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_AfRkRdet_partial_sup(const double ll, const double ul, const double tau, const double dm,
					const double tau_k_1_p, const double tau_k_2_p, 
					const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					const double fAn, const double fAp, const double fMn, const double fMp,
					const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fAn!=0.0)    Li += fAn    * norm_An_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fAp!=0.0)    Li += fAp    * norm_Ap_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fMn!=0.0)    Li += fMn    * norm_Mn_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fMp!=0.0)    Li += fMp    * norm_Mp_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fEn_k1!=0.0) Li += fEn_k1 * norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det,  nsigma);
  if(fEp_k1!=0.0) Li += fEp_k1 * norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det,  nsigma);
  if(fEn_k2!=0.0) Li += fEn_k2 * norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det,  nsigma);
  if(fEp_k2!=0.0) Li += fEp_k2 * norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det,  nsigma);
  return Li;
}

double norm_AfRkRdet_partial(const double ll, const double ul, const double tau, const double dm,
			     const double dz,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);
 

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;

  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);
  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_AfRkRdet_partial_sup(ll, ul, tau, dm,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_AfRkRdet_partial_sup(ll, ul, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_AfRkRdet_partial_sup(ll, ul, tau, dm,
							tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_AfRkRdet_partial_sup(ll, ul, tau, dm,
							tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_AfRkRdet_partial_sup(ll, ul, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;

}

inline double norm_MfRkRdet_partial_sup(const double ll, const double ul, const double tau, const double dm,
					const double tau_k_1_p, const double tau_k_2_p, 
					const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					const double fMn, const double fMp, const double fAn, const double fAp,
					const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fMn!=0.0)    Li += fMn    * norm_Mn_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fMp!=0.0)    Li += fMp    * norm_Mp_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fAn!=0.0)    Li += fAn    * norm_An_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fAp!=0.0)    Li += fAp    * norm_Ap_conv_gauss(ll, ul, tau, dm,   mu_det,  nsigma);
  if(fEn_k1!=0.0) Li += fEn_k1 * norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det,  nsigma);
  if(fEp_k1!=0.0) Li += fEp_k1 * norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det,  nsigma);
  if(fEn_k2!=0.0) Li += fEn_k2 * norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det,  nsigma);
  if(fEp_k2!=0.0) Li += fEp_k2 * norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det,  nsigma);
  return Li;
}

double norm_MfRkRdet_partial(const double ll, const double ul, const double tau, const double dm,
			     const double dz,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param){

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, 1.0-fk2);
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_MfRkRdet_partial_sup(ll, ul, tau, dm,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_MfRkRdet_partial_sup(ll, ul, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_MfRkRdet_partial_sup(ll, ul, tau, dm,
							tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_MfRkRdet_partial_sup(ll, ul, tau, dm,
							tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_MfRkRdet_partial_sup(ll, ul, tau, dm,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
  
}
//////
inline double norm_EfRkRdetRnp_full_sup(const double ll, const double ul,
					const double fact_n, const double ntau_n,
					const double fact_p, const double ntau_p,
					const double fd, const double fp, const double tau_np_p, const double tau_np_n,
					const double mu_det, const double sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  double fEn = fact_n*fd;
  double fEp = fact_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0, fxEn = 0.0, fxEp = 0.0;
  add_EnEn_coef(&fEn, &fEn_np, &fxEn, ntau_n, tau_np_n, fact_n*nfn);
  add_EnEp_coef(&fEn, &fEp_np,        ntau_n, tau_np_p, fact_n*nfp);
  add_EpEn_coef(&fEp, &fEn_np,        ntau_p, tau_np_n, fact_p*nfn);
  add_EpEp_coef(&fEp, &fEp_np, &fxEp, ntau_p, tau_np_p, fact_p*nfp);
  double Li = 0.0;
  if(fEn!=0.0)    Li += fEn    *  norm_En_conv_gauss(ll, ul, ntau_n,   mu_det, sigma_det);
  if(fEp!=0.0)    Li += fEp    *  norm_Ep_conv_gauss(ll, ul, ntau_p,   mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  norm_En_conv_gauss(ll, ul, tau_np_n, mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  norm_Ep_conv_gauss(ll, ul, tau_np_p, mu_det, sigma_det);
  if(fxEn!=0.0)   Li += fxEn   * norm_xEn_conv_gauss(ll, ul, ntau_n,   mu_det, sigma_det);
  if(fxEp!=0.0)   Li += fxEp   * norm_xEp_conv_gauss(ll, ul, ntau_p,   mu_det, sigma_det);
  return Li;
}

double norm_EfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor,  const double tau,
				const double ak, const double ck,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param){
  if(ak==0.0){
    if(debugout("INFO")) std::printf("[EfRk_fullrec] ak==% e where ak should not be zero.", ak);
    return -DBL_MAX;
  }

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double r_ckak = ck/ak;
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = 0.5*(1-r_ckak);
  const double fact_p = 0.5*(1+r_ckak);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdetRnp_full_sup(ll, ul,
						    fact_n, ntau_n, fact_p, ntau_p,
						    fd, fp, tau_np_p, tau_np_n,
						    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdetRnp_full_sup(ll, ul,
						      fact_n, ntau_n, fact_p, ntau_p,
						      fd, fp, tau_np_p_tl, tau_np_n_tl,
						      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EfRkRdetRnp_full_sup(ll, ul,
							fact_n, ntau_n, fact_p, ntau_p,
							fd, fp, tau_np_p, tau_np_n,
							mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EfRkRdetRnp_full_sup(ll, ul,
							fact_n, ntau_n, fact_p, ntau_p,
							fd, fp, tau_np_p_tl, tau_np_n_tl,
							mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EfRkRdetRnp_full_sup(ll, ul,
						      fact_n, ntau_n, fact_p, ntau_p,
						      fd, fp, tau_np_p, tau_np_n,
						      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_AfRkRdetRnp_full_sup(const double ll, const double ul,
					const double fact_n, const double ntau_n, const double ndm_n,
					const double fact_p, const double ntau_p, const double ndm_p, const double ndmtau,
					const double fd, const double fp, const double tau_np_p, const double tau_np_n,
					const double mu_det, const double sigma_det){
  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_mn_n = -fact_n*ndmtau;
  const double w_mn_p = -fact_p*ndmtau;
  double fAn = fact_n*fd;
  double fAp = fact_p*fd;
  double fMn = w_mn_n*fd;
  double fMp = w_mn_p*fd;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(&fAn, &fMn, &fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_AnEp_coef(&fAn, &fMn, &fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_MnEn_coef(&fMn, &fAn, &fEn_np, ntau_n, ndm_n, tau_np_n, w_mn_n*nfn);
  add_MnEp_coef(&fMn, &fAn, &fEp_np, ntau_n, ndm_n, tau_np_p, w_mn_n*nfp);
  add_ApEn_coef(&fAp, &fMp, &fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_ApEp_coef(&fAp, &fMp, &fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);
  add_MpEn_coef(&fMp, &fAp, &fEn_np, ntau_p, ndm_p, tau_np_n, w_mn_p*nfn);
  add_MpEp_coef(&fMp, &fAp, &fEp_np, ntau_p, ndm_p, tau_np_p, w_mn_p*nfp);
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    * norm_An_conv_gauss(ll, ul, ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    * norm_Ap_conv_gauss(ll, ul, ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    * norm_Mn_conv_gauss(ll, ul, ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    * norm_Mp_conv_gauss(ll, ul, ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np * norm_En_conv_gauss(ll, ul, tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np * norm_Ep_conv_gauss(ll, ul, tau_np_p,      mu_det, sigma_det);
  return Li;
}

double norm_AfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor, 
				const double tau, const double dm,
				const double ak, const double ck,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);

  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_AfRkRdetRnp_full_sup(ll, ul,
						    fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						    fd, fp, tau_np_p, tau_np_n,
						    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_AfRkRdetRnp_full_sup(ll, ul,
						      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						      fd, fp, tau_np_p_tl, tau_np_n_tl,
						      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_AfRkRdetRnp_full_sup(ll, ul,
							fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
							fd, fp, tau_np_p, tau_np_n,
							mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_AfRkRdetRnp_full_sup(ll, ul,
							fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
							fd, fp, tau_np_p_tl, tau_np_n_tl,
							mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_AfRkRdetRnp_full_sup(ll, ul,
						      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau,
						      fd, fp, tau_np_p, tau_np_n,
						      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_MfRkRdetRnp_full_sup(const double ll, const double ul,
					const double fact_n, const double ntau_n, const double ndm_n,
					const double fact_p, const double ntau_p, const double ndm_p,
					const double ndmtau, const double cktau,
					const double fd, const double fp, const double tau_np_p, const double tau_np_n,
					const double mu_det, const double sigma_det){

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  const double w_an_n = fact_n*ndmtau;
  const double w_an_p = fact_p*ndmtau;
  double fMn = fact_n*fd;
  double fMp = fact_p*fd;
  double fAn = w_an_n*fd;
  double fAp = w_an_p*fd;
  double fEn_k = 0.0, fEp_k = 0.0, fxEn_k = 0.0, fxEp_k = 0.0;
  double fEn_np = 0.0, fEp_np = 0.0;
  add_AnEn_coef(&fAn, &fMn, &fEn_np, ntau_n, ndm_n, tau_np_n, w_an_n*nfn);
  add_AnEp_coef(&fAn, &fMn, &fEp_np, ntau_n, ndm_n, tau_np_p, w_an_n*nfp);
  add_MnEn_coef(&fMn, &fAn, &fEn_np, ntau_n, ndm_n, tau_np_n, fact_n*nfn);
  add_MnEp_coef(&fMn, &fAn, &fEp_np, ntau_n, ndm_n, tau_np_p, fact_n*nfp);
  add_ApEn_coef(&fAp, &fMp, &fEn_np, ntau_p, ndm_p, tau_np_n, w_an_p*nfn);
  add_ApEp_coef(&fAp, &fMp, &fEp_np, ntau_p, ndm_p, tau_np_p, w_an_p*nfp);
  add_MpEn_coef(&fMp, &fAp, &fEn_np, ntau_p, ndm_p, tau_np_n, fact_p*nfn);
  add_MpEp_coef(&fMp, &fAp, &fEp_np, ntau_p, ndm_p, tau_np_p, fact_p*nfp);
#if 0
  if(cktau<0.0){
    const double w_ex_n = fact_n*ndmtau*cktau;
    add_EnEn_coef(&fEn_k, &fEn_np, &fxEn_k, -cktau, tau_np_n, w_ex_n*nfn);
    add_EnEp_coef(&fEn_k, &fEp_np,          -cktau, tau_np_p, w_ex_n*nfp);
  }else if(cktau>0.0){
    const double w_ex_p = fact_p*ndmtau*cktau;
    add_EpEn_coef(&fEp_k, &fEn_np,          cktau, tau_np_n, w_ex_p*nfn);
    add_EpEp_coef(&fEp_k, &fEp_np, &fxEp_k, cktau, tau_np_p, w_ex_p*nfp);
  }
#endif  
  double Li = 0.0;
  if(fAn!=0.0)    Li += fAn    *  norm_An_conv_gauss(ll, ul,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fAp!=0.0)    Li += fAp    *  norm_Ap_conv_gauss(ll, ul,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fMn!=0.0)    Li += fMn    *  norm_Mn_conv_gauss(ll, ul,  ntau_n, ndm_n, mu_det, sigma_det);
  if(fMp!=0.0)    Li += fMp    *  norm_Mp_conv_gauss(ll, ul,  ntau_p, ndm_p, mu_det, sigma_det);
  if(fEn_k!=0.0)  Li += fEn_k  *  norm_En_conv_gauss(ll, ul, -cktau,         mu_det, sigma_det);
  if(fEp_k!=0.0)  Li += fEp_k  *  norm_Ep_conv_gauss(ll, ul,  cktau,         mu_det, sigma_det);
  if(fxEn_k!=0.0) Li += fxEn_k * norm_xEn_conv_gauss(ll, ul, -cktau,         mu_det, sigma_det);
  if(fxEp_k!=0.0) Li += fxEp_k * norm_xEp_conv_gauss(ll, ul,  cktau,         mu_det, sigma_det);
  if(fEn_np!=0.0) Li += fEn_np *  norm_En_conv_gauss(ll, ul,  tau_np_n,      mu_det, sigma_det);
  if(fEp_np!=0.0) Li += fEp_np *  norm_Ep_conv_gauss(ll, ul,  tau_np_p,      mu_det, sigma_det);
  return Li;
}

double norm_MfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor, 
				const double tau, const double dm,
				const double ak, const double ck,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[MfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }
  const double r_ckak = ck/ak;
  const double inv_ak = 1.0/ak;
  const double ndmtau = r_ckak*dm*tau;
  const double _fact = 1.0/(1.0+ndmtau*ndmtau);
  const double ndm_n  = dm/(ak-ck);
  const double ndm_p  = dm/(ak+ck);
  const double ntau_n = tau*(ak-ck);
  const double ntau_p = tau*(ak+ck);
  const double fact_n = inv_ak*_fact;
  const double fact_p = inv_ak*_fact;
  const double cktau = ck*tau;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_MfRkRdetRnp_full_sup(ll, ul,
						    fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						    fd, fp, tau_np_p, tau_np_n,
						    mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_MfRkRdetRnp_full_sup(ll, ul,
						      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						      fd, fp, tau_np_p_tl, tau_np_n_tl,
						      mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_MfRkRdetRnp_full_sup(ll, ul,
							fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
							fd, fp, tau_np_p, tau_np_n,
							mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_MfRkRdetRnp_full_sup(ll, ul,
							fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
							fd, fp, tau_np_p_tl, tau_np_n_tl,
							mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_MfRkRdetRnp_full_sup(ll, ul,
						      fact_n, ntau_n, ndm_n, fact_p, ntau_p, ndm_p, ndmtau, cktau,
						      fd, fp, tau_np_p, tau_np_n,
						      mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_EfRkRdetRnp_partial_sup(const double ll, const double ul, const double tau,
					   const double tau_k_1_p, const double tau_k_2_p, 
					   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					   const double tau_np_p, const double tau_np_n,
					   const double fEn, const double fEp,
					   const double fxEn, const double fxEp, const double fxxEn, const double fxxEp,
					   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					   const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
					   const double fEn_np, const double fEp_np, 
					   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEn!=0.0)     Li += fEn     *   norm_En_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fEp!=0.0)     Li += fEp     *   norm_Ep_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fxEn!=0.0)    Li += fxEn    *  norm_xEn_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fxEp!=0.0)    Li += fxEp    *  norm_xEp_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fxxEn!=0.0)   Li += fxxEn   * norm_xxEn_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fxxEp!=0.0)   Li += fxxEp   * norm_xxEp_conv_gauss(ll, ul, tau,       mu_det, nsigma);
  if(fEn_k1!=0.0)  Li += fEn_k1  *   norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)  Li += fEp_k1  *   norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)  Li += fEn_k2  *   norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)  Li += fEp_k2  *   norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0) Li += fxEn_k1 *  norm_xEn_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0) Li += fxEp_k1 *  norm_xEp_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0) Li += fxEn_k2 *  norm_xEn_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0) Li += fxEp_k2 *  norm_xEp_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)  Li += fEn_np  *   norm_En_conv_gauss(ll, ul, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)  Li += fEp_np  *   norm_Ep_conv_gauss(ll, ul, tau_np_p,  mu_det, nsigma);
  return Li;
}

double norm_EfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;
  
  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0, fxxEn=0.0, fxxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;

  add_EfEnp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fxEn, &fxEp,
		 tau, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_EfEnp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fxEn, &fxEp,
		 tau, tau_k_2_n, tau_k_2_p, fk2*fd);
  
  add_EfEnpEn_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fEn_np,
		   &fxEn, &fxEp, &fxEn_k1, &fxxEn,
		   tau, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_EfEnpEn_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fEn_np,
		   &fxEn, &fxEp, &fxEn_k2, &fxxEn,
		   tau, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_EfEnpEp_coef(&fEn, &fEp, &fEn_k1, &fEp_k1, &fEp_np,
		   &fxEn, &fxEp, &fxEp_k1, &fxxEp,
		   tau, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_EfEnpEp_coef(&fEn, &fEp, &fEn_k2, &fEp_k2, &fEp_np,
		   &fxEn, &fxEp, &fxEp_k2, &fxxEp,
		   tau, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdetRnp_partial_sup(ll, ul, tau,
						       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						       tau_np_p, tau_np_n,
						       fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						       fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						       fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						       fEn_np, fEp_np, 
						       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fEn_tl=0.0, fEp_tl=0.0, fxEn_tl=0.0, fxEp_tl=0.0, fxxEn_tl=0.0, fxxEp_tl=0.0; 
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
    
    add_EfEnp_coef(&fEn_tl, &fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fxEn_tl, &fxEp_tl,
		   tau, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_EfEnp_coef(&fEn_tl, &fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fxEn_tl, &fxEp_tl,
		   tau, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_EfEnpEn_coef(&fEn_tl, &fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEn_k1_tl, &fxxEn_tl,
		     tau, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_EfEnpEn_coef(&fEn_tl, &fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEn_k2_tl, &fxxEn_tl,
		     tau, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_EfEnpEp_coef(&fEn_tl, &fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEp_k1_tl, &fxxEp_tl,
		     tau, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_EfEnpEp_coef(&fEn_tl, &fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_tl, &fxEp_k2_tl, &fxxEp_tl,
		     tau, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdetRnp_partial_sup(ll, ul, tau,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p_tl, tau_np_n_tl,
							 fEn_tl, fEp_tl, fxEn_tl, fxEp_tl, fxxEn_tl, fxxEp_tl,
							 fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
							 fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl,
							 fEn_np_tl, fEp_np_tl, 
							 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/

      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EfRkRdetRnp_partial_sup(ll, ul, tau,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p,  tau_np_n, 
							   fEn,  fEp,  fxEn,  fxEp,  fxxEn,  fxxEp, 
							   fEn_k1,  fEp_k1,  fEn_k2,  fEp_k2, 
							   fxEn_k1,  fxEp_k1,  fxEn_k2,  fxEp_k2, 
							   fEn_np,  fEp_np,  
							   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EfRkRdetRnp_partial_sup(ll, ul, tau,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p_tl,  tau_np_n_tl, 
							   fEn_tl,  fEp_tl,  fxEn_tl,  fxEp_tl,  fxxEn_tl,  fxxEp_tl, 
							   fEn_k1_tl,  fEp_k1_tl,  fEn_k2_tl,  fEp_k2_tl, 
							   fxEn_k1_tl,  fxEp_k1_tl,  fxEn_k2_tl,  fxEp_k2_tl, 
							   fEn_np_tl,  fEp_np_tl,  
							   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EfRkRdetRnp_partial_sup(ll, ul, tau,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p, tau_np_n,
							 fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
							 fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							 fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
							 fEn_np, fEp_np, 
							 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double norm_AfRkRdetRnp_partial_sup(const double ll, const double ul, const double tau, const double dm,
					   const double tau_k_1_p, const double tau_k_2_p, 
					   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					   const double tau_np_p, const double tau_np_n,
					   const double fAn, const double fAp, const double fMn, const double fMp,
					   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					   const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2, 
					   const double fEn_np, const double fEp_np,
					   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fAn!=0.0)      Li += fAn     *  norm_An_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fAp!=0.0)      Li += fAp     *  norm_Ap_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fMn!=0.0)      Li += fMn     *  norm_Mn_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fMp!=0.0)      Li += fMp     *  norm_Mp_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fEn_k1!=0.0)   Li += fEn_k1  *  norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)   Li += fEp_k1  *  norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)   Li += fEn_k2  *  norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)   Li += fEp_k2  *  norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0)  Li += fxEn_k1 * norm_xEn_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0)  Li += fxEp_k1 * norm_xEp_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0)  Li += fxEn_k2 * norm_xEn_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0)  Li += fxEp_k2 * norm_xEp_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)   Li += fEn_np  *  norm_En_conv_gauss(ll, ul, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)   Li += fEp_np  *  norm_Ep_conv_gauss(ll, ul, tau_np_p,  mu_det, nsigma);
  return Li;
}

double norm_AfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau, const double dm,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);
 
  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;
  
  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;
  
  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
  //tau, dm, tau_k_1_n, tau_k_1_p, fk2*fd); // BUGBUGBUG
  add_AfEnp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);

  add_AfEnpEn_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1, &fEn_np, &fxEn_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_AfEnpEn_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2, &fEn_np, &fxEn_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_AfEnpEp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k1, &fEp_k1, &fEp_np, &fxEp_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_AfEnpEp_coef(&fAn, &fAp, &fMn, &fMp, &fEn_k2, &fEp_k2, &fEp_np, &fxEp_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_AfRkRdetRnp_partial_sup(ll, ul, tau, dm,
						       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						       tau_np_p, tau_np_n,
						       fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						       fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						       mu_mm, sigma_mm);
  if(ftail_asc>0.0){

    double fAn_tl=0.0, fAp_tl=0.0, fMn_tl=0.0, fMp_tl=0.0;
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
  
    add_AfEnp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k1_tl, &fEp_k1_tl,
		   tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
    //tau, dm, tau_k_1_n, tau_k_1_p, fk2*fd); // BUGBUGBUG
    add_AfEnp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k2_tl, &fEp_k2_tl,
		   tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_AfEnpEn_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl, &fxEn_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_AfEnpEn_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl, &fxEn_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_AfEnpEp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl, &fxEp_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_AfEnpEp_coef(&fAn_tl, &fAp_tl, &fMn_tl, &fMp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl, &fxEp_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);

    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_AfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p_tl, tau_np_n_tl,
							 fAn_tl, fAp_tl, fMn_tl, fMp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
							 fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
							 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_AfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p, tau_np_n,
							   fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							   fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
							   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_AfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p_tl, tau_np_n_tl,
							   fAn_tl, fAp_tl, fMn_tl, fMp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
							   fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
							   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_AfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p, tau_np_n,
							 fAn, fAp, fMn, fMp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							 fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
							 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;

}

inline double norm_MfRkRdetRnp_partial_sup(const double ll, const double ul, const double tau, const double dm,
					   const double tau_k_1_p, const double tau_k_2_p, 
					   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					   const double tau_np_p, const double tau_np_n,
					   const double fMn, const double fMp, const double fAn, const double fAp,
					   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					   const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
					   const double fEn_np, const double fEp_np,
					   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fMn!=0.0)      Li += fMn     *  norm_Mn_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fMp!=0.0)      Li += fMp     *  norm_Mp_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fAn!=0.0)      Li += fAn     *  norm_An_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fAp!=0.0)      Li += fAp     *  norm_Ap_conv_gauss(ll, ul, tau,   dm, mu_det, nsigma);
  if(fEn_k1!=0.0)   Li += fEn_k1  *  norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)   Li += fEp_k1  *  norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)   Li += fEn_k2  *  norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)   Li += fEp_k2  *  norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0)  Li += fxEn_k1 * norm_xEn_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0)  Li += fxEp_k1 * norm_xEp_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0)  Li += fxEn_k2 * norm_xEn_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0)  Li += fxEp_k2 * norm_xEp_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)   Li += fEn_np  *  norm_En_conv_gauss(ll, ul, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)   Li += fEp_np  *  norm_Ep_conv_gauss(ll, ul, tau_np_p,  mu_det, nsigma);
  return Li;
}

double norm_MfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau, const double dm,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;

  double fAn=0.0, fAp=0.0, fMn=0.0, fMp=0.0;
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;
  
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1,
		 tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_MfEnp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2,
		 tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);

  add_MfEnpEn_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1, &fEn_np, &fxEn_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_MfEnpEn_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2, &fEn_np, &fxEn_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_MfEnpEp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k1, &fEp_k1, &fEp_np, &fxEp_k1,
		   tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_MfEnpEp_coef(&fMn, &fMp, &fAn, &fAp, &fEn_k2, &fEp_k2, &fEp_np, &fxEp_k2,
		   tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_MfRkRdetRnp_partial_sup(ll, ul, tau, dm,
						       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						       tau_np_p, tau_np_n,
						       fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						       fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
						       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fAn_tl=0.0, fAp_tl=0.0, fMn_tl=0.0, fMp_tl=0.0;
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
  
    add_MfEnp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k1_tl, &fEp_k1_tl,
		   tau, dm, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_MfEnp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k2_tl, &fEp_k2_tl,
		   tau, dm, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_MfEnpEn_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl, &fxEn_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_MfEnpEn_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl, &fxEn_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_MfEnpEp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl, &fxEp_k1_tl,
		     tau, dm, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_MfEnpEp_coef(&fMn_tl, &fMp_tl, &fAn_tl, &fAp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl, &fxEp_k2_tl,
		     tau, dm, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_MfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p_tl, tau_np_n_tl,
							 fMn_tl, fMp_tl, fAn_tl, fAp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
							 fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
							 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_MfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p, tau_np_n,
							   fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							   fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
							   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_MfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p_tl, tau_np_n_tl,
							   fMn_tl, fMp_tl, fAn_tl, fAp_tl, fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
							   fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl, fEn_np_tl, fEp_np_tl,
							   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_MfRkRdetRnp_partial_sup(ll, ul, tau, dm,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p, tau_np_n,
							 fMn, fMp, fAn, fAp, fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							 fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2, fEn_np, fEp_np,
							 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
  
}

double EfRdet(const double x, const double tau,
	      const int ntrk_rec, const double sz_rec,
	      const double chisq_z_rec, const int ndf_z_rec,
	      const int ntrk_asc, const double sz_asc,
	      const double chisq_z_asc, const int ndf_z_asc,
	      const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = Ef_conv_gauss(x, tau,  mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = Ef_conv_gauss(x, tau,  mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = Ef_conv_gauss(x, tau,  mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = Ef_conv_gauss(x, tau,  mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = Ef_conv_gauss(x, tau,  mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double AfRdet(const double x, const double tau, const double dm,
	      const int ntrk_rec, const double sz_rec,
	      const double chisq_z_rec, const int ndf_z_rec,
	      const int ntrk_asc, const double sz_asc,
	      const double chisq_z_asc, const int ndf_z_asc,
	      const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = Af_conv_gauss(x, tau, dm, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = Af_conv_gauss(x, tau, dm, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = Af_conv_gauss(x, tau, dm, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = Af_conv_gauss(x, tau, dm, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = Af_conv_gauss(x, tau, dm, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double MfRdet(const double x, const double tau, const double dm,
	      const int ntrk_rec, const double sz_rec,
	      const double chisq_z_rec, const int ndf_z_rec,
	      const int ntrk_asc, const double sz_asc,
	      const double chisq_z_asc, const int ndf_z_asc,
	      const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = Mf_conv_gauss(x, tau, dm, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = Mf_conv_gauss(x, tau, dm, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = Mf_conv_gauss(x, tau, dm, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = Mf_conv_gauss(x, tau, dm, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = Mf_conv_gauss(x, tau, dm, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double norm_EfRdet(const double ll, const double ul, const double tau,
		   const int ntrk_rec, const double sz_rec,
		   const double chisq_z_rec, const int ndf_z_rec,
		   const int ntrk_asc, const double sz_asc,
		   const double chisq_z_asc, const int ndf_z_asc,
		   const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_Ef_conv_gauss(ll, ul, tau, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_Ef_conv_gauss(ll, ul, tau, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_Ef_conv_gauss(ll, ul, tau, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_Ef_conv_gauss(ll, ul, tau, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_Ef_conv_gauss(ll, ul, tau, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double norm_AfRdet(const double ll, const double ul,
		   const double tau, const double dm,
		   const int ntrk_rec, const double sz_rec,
		   const double chisq_z_rec, const int ndf_z_rec,
		   const int ntrk_asc, const double sz_asc,
		   const double chisq_z_asc, const int ndf_z_asc,
		   const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_Af_conv_gauss(ll, ul, tau, dm, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_Af_conv_gauss(ll, ul, tau, dm, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_Af_conv_gauss(ll, ul, tau, dm, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_Af_conv_gauss(ll, ul, tau, dm, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_Af_conv_gauss(ll, ul, tau, dm, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double norm_MfRdet(const double ll, const double ul,
		   const double tau, const double dm,
		   const int ntrk_rec, const double sz_rec,
		   const double chisq_z_rec, const int ndf_z_rec,
		   const int ntrk_asc, const double sz_asc,
		   const double chisq_z_asc, const int ndf_z_asc,
		   const dtres_param_t * const param){

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_Mf_conv_gauss(ll, ul, tau, dm, mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_Mf_conv_gauss(ll, ul, tau, dm, mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_Mf_conv_gauss(ll, ul, tau, dm, mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_Mf_conv_gauss(ll, ul, tau, dm, mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_Mf_conv_gauss(ll, ul, tau, dm, mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}



static int is_dtres_SVD1_by2005(const dtres_param_t * const param){
  if(param==&dtres_param_2002moriond || param==&dtres_param_2002moriond_mc)
    return dt_resol_global::moriond2002;
  else if(param==&dtres_param_2002ichep60 || param==&dtres_param_2002ichep60_mc)
    return dt_resol_global::ichep2002_60;
  else if(param==&dtres_param_2002ichep80 || param==&dtres_param_2002ichep80_mc)
    return dt_resol_global::ichep2002_80;
  else if(param==&dtres_param_2003lp03 || param==&dtres_param_2003lp03_mc)
    return dt_resol_global::lp2003_150;
  else if(param==&dtres_param_svd1_ichep04 || param==&dtres_param_svd1_ichep04_mc)
    return dt_resol_global::svd1_ichep04;
  else
    if(param->ftl_rec_mlt[0]==0. && param->ftl_rec_mlt[1]==0. &&
       param->ftl_asc_mlt[0]==0. && param->ftl_asc_mlt[1]==0. &&
       param->ftl_asc==0. &&
       param->fd_np_sgl[0][0] == param->fd_np_sgl[0][1] &&
       param->fd_np_mlt[0][0] == param->fd_np_mlt[0][1])
      return 9;
    else
      return 0;
}

double add_outlier(const int expno,
		   const double x,
		   const double Lin,
		   const int ntrk_rec, const int ntrk_asc,
		   const dtres_param_t * const param){
  const double fol
    = ((expno>29||!is_dtres_SVD1_by2005(param)||ntrk_rec>1)
       &&ntrk_asc>1) ? param->fol_mul : param->fol_sgl;
  const double m = 0.0;
  const double Li = (1.0-fol)*Lin + fol*gaussian(x, m, param->sig_ol);
  return Li;
}

double add_norm_outlier(const int expno,
			const double ll, const double ul,
			const double nLin,
			const int ntrk_rec, const int ntrk_asc,
			const dtres_param_t * const param){
  const double fol
    = ((expno>29||!is_dtres_SVD1_by2005(param)||ntrk_rec>1)
       &&ntrk_asc>1) ? param->fol_mul : param->fol_sgl;
  const double m = 0.0;
  const double Li
    = (1.0-fol)*nLin + fol*norm_gaussian(ll, ul, m, param->sig_ol);
  return Li;
}



void dtres_systematics(const dtres_syst_item_t item,
		       const double sigma,
		       const int expno,
		       dtres_param_t* param,
		       const dtres_param_t* const posi_err,
		       const dtres_param_t* const nega_err,
		       const double mc_param_factor, const int oldtype){
  const double fdscale
    = (oldtype||param==&dtres_param_2002moriond) ? mc_param_factor : 1.0;

  const double snpscale
    = (expno>29)? mc_param_factor : 1.0;
  
  const dtres_param_t& error = (sigma>=0.0) ? *posi_err : *nega_err;
  switch(item){
  case dtres_Rdet_Srec0:
    (param->Srec)[0] += sigma*fabs((error.Srec)[0]);
    break;
  case dtres_Rdet_Srec1:
    (param->Srec)[1] += sigma*fabs((error.Srec)[1]);
    break;
  case dtres_Rdet_Sasc0:
    (param->Sasc)[0] += sigma*fabs((error.Sasc)[0]);
    break;
  case dtres_Rdet_Sasc1:
    (param->Sasc)[1] += sigma*fabs((error.Sasc)[1]);
    break;
  case dtres_Rnp_Snp   :
    if(param->Snp>DTRES_EXTERAM_THRE){
      param->Snp += snpscale*fabs(error.Snp);
    }
    break;
  case dtres_Rnp_Snp_global:
//     if(param->Snp>DTRES_EXTERAM_THRE){
      param->Snp_global += sigma*fabs(error.Snp_global);
//     }
    break;
  case dtres_Rdet_Smn_onetrk:
    param->Smn_rec += sigma*fabs(error.Smn_rec);
    param->Smn_asc += sigma*fabs(error.Smn_asc);
    break;
  case dtres_Rdet_Stl_onetrk:
    param->Stl_rec += sigma*fabs(error.Stl_rec);
    param->Stl_asc += sigma*fabs(error.Stl_asc);
    break;
  case dtres_Rdet_ftl_onetrk:
    param->ftl_rec += sigma*fabs(error.ftl_rec);
    param->ftl_asc += sigma*fabs(error.ftl_asc);
    break;
    //
  case dtres_Rdet_rec_mlt_ftail:
    (param->ftl_rec_mlt)[0] += sigma*fabs((error.ftl_rec_mlt)[0]);
    (param->ftl_rec_mlt)[1] += sigma*fabs((error.ftl_rec_mlt)[1]);
    break;
  case dtres_Rdet_rec_mlt_stail:
    param->Stl_rec_mlt += sigma*fabs(error.Stl_rec_mlt);
    break;
  case dtres_Rdet_asc_mlt_ftail:
    (param->ftl_asc_mlt)[0] += sigma*fabs((error.ftl_asc_mlt)[0]);
    (param->ftl_asc_mlt)[1] += sigma*fabs((error.ftl_asc_mlt)[1]);
    break;
  case dtres_Rdet_asc_mlt_stail:
    param->Stl_asc_mlt += sigma*fabs(error.Stl_asc_mlt);
    break;
  case dtres_Rnp_fd_np_sgl0_Bzero:
    (param->fd_np_sgl)[0][0] += sigma*fdscale*fabs((error.fd_np_sgl)[0][0]);
    break;
  case dtres_Rnp_fd_np_sgl1_Bzero:
    (param->fd_np_sgl)[0][1] += sigma*fdscale*fabs((error.fd_np_sgl)[0][1]);
    break;
  case dtres_Rnp_fp_np_sgl_Bzero:
    (param->fp_np_sgl)[0] += sigma*mc_param_factor*fabs((error.fp_np_sgl)[0]);
    break;
  case dtres_Rnp_tau_np_p_sgl0_Bzero:
    (param->tau_np_p_sgl)[0][0] += sigma*mc_param_factor*fabs((error.tau_np_p_sgl)[0][0]);
    break;
  case dtres_Rnp_tau_np_p_sgl1_Bzero:
    (param->tau_np_p_sgl)[0][1] += sigma*mc_param_factor*fabs((error.tau_np_p_sgl)[0][1]);
    break;
  case dtres_Rnp_tau_np_n_sgl0_Bzero:
    (param->tau_np_n_sgl)[0][0] += sigma*mc_param_factor*fabs((error.tau_np_n_sgl)[0][0]);
    break;
  case dtres_Rnp_tau_np_n_sgl1_Bzero:
    (param->tau_np_n_sgl)[0][1] += sigma*mc_param_factor*fabs((error.tau_np_n_sgl)[0][1]);
    break;
  case dtres_Rnp_fd_np_mlt0_Bzero:
    (param->fd_np_mlt)[0][0] += sigma*fdscale*fabs((error.fd_np_mlt)[0][0]);
    break;
  case dtres_Rnp_fd_np_mlt1_Bzero:
    (param->fd_np_mlt)[0][1] += sigma*fdscale*fabs((error.fd_np_mlt)[0][1]);
    break;
  case dtres_Rnp_fd_np_st_mlt_Bzero:
    (param->fd_np_st_mlt)[0] += sigma*mc_param_factor*fabs((error.fd_np_st_mlt)[0]);
    break;
  case dtres_Rnp_fd_np_xi_mlt_Bzero:
    (param->fd_np_xi_mlt)[0] += sigma*mc_param_factor*fabs((error.fd_np_xi_mlt)[0]);
    break;
  case dtres_Rnp_fd_np_stxi_mlt_Bzero:
    (param->fd_np_stxi_mlt)[0] += sigma*mc_param_factor*fabs((error.fd_np_stxi_mlt)[0]);
    break;
  case dtres_Rnp_fp_np_mlt_Bzero:
    (param->fp_np_mlt)[0] += sigma*mc_param_factor*fabs((error.fp_np_mlt)[0]);
    break;
  case dtres_Rnp_fn_np_mlt_Bzero:
    (param->fn_np_mlt)[0] += sigma*mc_param_factor*fabs((error.fn_np_mlt)[0]);
    break;
  case dtres_Rnp_tau_np_p_mlt0_Bzero:
    (param->tau_np_p_mlt)[0][0] += sigma*mc_param_factor*fabs((error.tau_np_p_mlt)[0][0]);
    break;
  case dtres_Rnp_tau_np_p_mlt1_Bzero:
    (param->tau_np_p_mlt)[0][1] += sigma*mc_param_factor*fabs((error.tau_np_p_mlt)[0][1]);
    break;
  case dtres_Rnp_tau_np_p_xi_mlt_Bzero:
    (param->tau_np_p_xi_mlt)[0] += sigma*mc_param_factor*fabs((error.tau_np_p_xi_mlt)[0]);
    break;
  case dtres_Rnp_tau_np_p_stxi_mlt_Bzero:
    (param->tau_np_p_stxi_mlt)[0] += sigma*mc_param_factor*fabs((error.tau_np_p_stxi_mlt)[0]);
    break;
  case dtres_Rnp_tau_np_n_mlt0_Bzero:
    (param->tau_np_n_mlt)[0][0] += sigma*mc_param_factor*fabs((error.tau_np_n_mlt)[0][0]);
    break;
  case dtres_Rnp_tau_np_n_mlt1_Bzero:
    (param->tau_np_n_mlt)[0][1] += sigma*mc_param_factor*fabs((error.tau_np_n_mlt)[0][1]);
    break;
  case dtres_Rnp_tau_np_n_xi_mlt_Bzero:
    (param->tau_np_n_xi_mlt)[0] += sigma*mc_param_factor*fabs((error.tau_np_n_xi_mlt)[0]);
    break;
  case dtres_Rnp_tau_np_n_stxi_mlt_Bzero:
    (param->tau_np_n_stxi_mlt)[0] += sigma*mc_param_factor*fabs((error.tau_np_n_stxi_mlt)[0]);
    break;
  case dtres_Rnp_fd_np_sgl0_Bplus:
    (param->fd_np_sgl)[1][0] += sigma*fdscale*fabs((error.fd_np_sgl)[1][0]);
    break;
  case dtres_Rnp_fd_np_sgl1_Bplus:
    (param->fd_np_sgl)[1][1] += sigma*fdscale*fabs((error.fd_np_sgl)[1][1]);
    break;
  case dtres_Rnp_fp_np_sgl_Bplus:
    (param->fp_np_sgl)[1] += sigma*mc_param_factor*fabs((error.fp_np_sgl)[1]);
    break;
  case dtres_Rnp_tau_np_p_sgl0_Bplus:
    (param->tau_np_p_sgl)[1][0] += sigma*mc_param_factor*fabs((error.tau_np_p_sgl)[1][0]);
    break;
  case dtres_Rnp_tau_np_p_sgl1_Bplus:
    (param->tau_np_p_sgl)[1][1] += sigma*mc_param_factor*fabs((error.tau_np_p_sgl)[1][1]);
    break;
  case dtres_Rnp_tau_np_n_sgl0_Bplus:
    (param->tau_np_n_sgl)[1][0] += sigma*mc_param_factor*fabs((error.tau_np_n_sgl)[1][0]);
    break;
  case dtres_Rnp_tau_np_n_sgl1_Bplus:
    (param->tau_np_n_sgl)[1][1] += sigma*mc_param_factor*fabs((error.tau_np_n_sgl)[1][1]);
    break;
  case dtres_Rnp_fd_np_mlt0_Bplus:
    (param->fd_np_mlt)[1][0] += sigma*fdscale*fabs((error.fd_np_mlt)[1][0]);
    break;
  case dtres_Rnp_fd_np_mlt1_Bplus:
    (param->fd_np_mlt)[1][1] += sigma*fdscale*fabs((error.fd_np_mlt)[1][1]);
    break;
  case dtres_Rnp_fd_np_st_mlt_Bplus:
    (param->fd_np_st_mlt)[1] += sigma*mc_param_factor*fabs((error.fd_np_st_mlt)[1]);
    break;
  case dtres_Rnp_fd_np_xi_mlt_Bplus:
    (param->fd_np_xi_mlt)[1] += sigma*mc_param_factor*fabs((error.fd_np_xi_mlt)[1]);
    break;
  case dtres_Rnp_fd_np_stxi_mlt_Bplus:
    (param->fd_np_stxi_mlt)[1] += sigma*mc_param_factor*fabs((error.fd_np_stxi_mlt)[1]);
    break;
  case dtres_Rnp_fp_np_mlt_Bplus:
    (param->fp_np_mlt)[1] += sigma*mc_param_factor*fabs((error.fp_np_mlt)[1]);
    break;
  case dtres_Rnp_fn_np_mlt_Bplus:
    (param->fn_np_mlt)[1] += sigma*mc_param_factor*fabs((error.fn_np_mlt)[1]);
    break;
  case dtres_Rnp_tau_np_p_mlt0_Bplus:
    (param->tau_np_p_mlt)[1][0] += sigma*mc_param_factor*fabs((error.tau_np_p_mlt)[1][0]);
    break;
  case dtres_Rnp_tau_np_p_mlt1_Bplus:
    (param->tau_np_p_mlt)[1][1] += sigma*mc_param_factor*fabs((error.tau_np_p_mlt)[1][1]);
    break;
  case dtres_Rnp_tau_np_p_xi_mlt_Bplus:
    (param->tau_np_p_xi_mlt)[1] += sigma*mc_param_factor*fabs((error.tau_np_p_xi_mlt)[1]);
    break;
  case dtres_Rnp_tau_np_p_stxi_mlt_Bplus:
    (param->tau_np_p_stxi_mlt)[1] += sigma*mc_param_factor*fabs((error.tau_np_p_stxi_mlt)[1]);
    break;
  case dtres_Rnp_tau_np_n_mlt0_Bplus:
    (param->tau_np_n_mlt)[1][0] += sigma*mc_param_factor*fabs((error.tau_np_n_mlt)[1][0]);
    break;
  case dtres_Rnp_tau_np_n_mlt1_Bplus:
    (param->tau_np_n_mlt)[1][1] += sigma*mc_param_factor*fabs((error.tau_np_n_mlt)[1][1]);
    break;
  case dtres_Rnp_tau_np_n_xi_mlt_Bplus:
    (param->tau_np_n_xi_mlt)[1] += sigma*mc_param_factor*fabs((error.tau_np_n_xi_mlt)[1]);
    break;
  case dtres_Rnp_tau_np_n_stxi_mlt_Bplus:
    (param->tau_np_n_stxi_mlt)[1] += sigma*mc_param_factor*fabs((error.tau_np_n_stxi_mlt)[1]);
    break;
  case dtres_Rk_tau_k_1_p0:
    (param->tau_k_1_p)[0] += sigma*mc_param_factor*fabs((error.tau_k_1_p)[0]);
    break;
  case dtres_Rk_tau_k_1_p1:
    (param->tau_k_1_p)[1] += sigma*mc_param_factor*fabs((error.tau_k_1_p)[1]);
    break;
  case dtres_Rk_tau_k_2_p0:
    (param->tau_k_2_p)[0] += sigma*mc_param_factor*fabs((error.tau_k_2_p)[0]);
    break;
  case dtres_Rk_tau_k_2_p1:
    (param->tau_k_2_p)[1] += sigma*mc_param_factor*fabs((error.tau_k_2_p)[1]);
    break;
  case dtres_Rk_tau_k_1_n0:
    (param->tau_k_1_n)[0] += sigma*mc_param_factor*fabs((error.tau_k_1_n)[0]);
    break;
  case dtres_Rk_tau_k_1_n1:
    (param->tau_k_1_n)[1] += sigma*mc_param_factor*fabs((error.tau_k_1_n)[1]);
    break;
  case dtres_Rk_tau_k_2_n0:
    (param->tau_k_2_n)[0] += sigma*mc_param_factor*fabs((error.tau_k_2_n)[0]);
    break;
  case dtres_Rk_tau_k_2_n1:
    (param->tau_k_2_n)[1] += sigma*mc_param_factor*fabs((error.tau_k_2_n)[1]);
    break;
  case dtres_Rk_sigma_k0:
    (param->sigma_k)[0] += sigma*mc_param_factor*fabs((error.sigma_k)[0]);
    break;
  case dtres_Rk_sigma_k1:
    (param->sigma_k)[1] += sigma*mc_param_factor*fabs((error.sigma_k)[1]);
    break;
  case dtres_Rk_fk2:
    param->fk2 += sigma*mc_param_factor*fabs(error.fk2);
    break;
  case dtres_Rol_sig_ol:
    param->sig_ol += sigma*fabs(error.sig_ol);
    break;
  case dtres_Rol_fol_sgl:
    param->fol_sgl += sigma*fabs(error.fol_sgl);
    break;
  case dtres_Rol_fol_mul:
    param->fol_mul += sigma*fabs(error.fol_mul);
    break;
  case dtres_Dt_cutoff:
    dt_resol_global::dt_llmt -= sigma*fabs(dt_resol_global::dt_lmt_var);
    dt_resol_global::dt_ulmt += sigma*fabs(dt_resol_global::dt_lmt_var);
    break;
  default:
    return;
  }
}

const dtres_syst_item_t dtres_items_fullrecon[] = {
  dtres_Rdet_Srec0, dtres_Rdet_Srec1,
  dtres_Rdet_Sasc0, dtres_Rdet_Sasc1, dtres_Rnp_Snp   , dtres_Rnp_Snp_global,
  dtres_Rdet_Smn_onetrk, dtres_Rdet_Stl_onetrk, dtres_Rdet_ftl_onetrk,
  dtres_Rdet_rec_mlt_ftail, dtres_Rdet_rec_mlt_stail,
  dtres_Rdet_asc_mlt_ftail, dtres_Rdet_asc_mlt_stail,
  dtres_Rnp_fd_np_sgl0_Bzero, dtres_Rnp_fd_np_sgl1_Bzero, dtres_Rnp_fp_np_sgl_Bzero,
  dtres_Rnp_tau_np_p_sgl0_Bzero, dtres_Rnp_tau_np_p_sgl1_Bzero,
  dtres_Rnp_tau_np_n_sgl0_Bzero, dtres_Rnp_tau_np_n_sgl1_Bzero,
  dtres_Rnp_fd_np_mlt0_Bzero, dtres_Rnp_fd_np_mlt1_Bzero,
  dtres_Rnp_fd_np_xi_mlt_Bzero, dtres_Rnp_fd_np_stxi_mlt_Bzero,
  dtres_Rnp_fp_np_mlt_Bzero, dtres_Rnp_fn_np_mlt_Bzero,
  dtres_Rnp_tau_np_p_mlt0_Bzero, dtres_Rnp_tau_np_p_mlt1_Bzero,
  dtres_Rnp_tau_np_p_xi_mlt_Bzero, dtres_Rnp_tau_np_p_stxi_mlt_Bzero,
  dtres_Rnp_tau_np_n_mlt0_Bzero, dtres_Rnp_tau_np_n_mlt1_Bzero,
  dtres_Rnp_tau_np_n_xi_mlt_Bzero, dtres_Rnp_tau_np_n_stxi_mlt_Bzero,
  dtres_Rnp_fd_np_sgl0_Bplus, dtres_Rnp_fd_np_sgl1_Bplus, dtres_Rnp_fp_np_sgl_Bplus,
  dtres_Rnp_tau_np_p_sgl0_Bplus, dtres_Rnp_tau_np_p_sgl1_Bplus,
  dtres_Rnp_tau_np_n_sgl0_Bplus, dtres_Rnp_tau_np_n_sgl1_Bplus,
  dtres_Rnp_fd_np_mlt0_Bplus, dtres_Rnp_fd_np_mlt1_Bplus,
  dtres_Rnp_fd_np_xi_mlt_Bplus, dtres_Rnp_fd_np_stxi_mlt_Bplus,
  dtres_Rnp_fp_np_mlt_Bplus, dtres_Rnp_fn_np_mlt_Bplus,
  dtres_Rnp_tau_np_p_mlt0_Bplus, dtres_Rnp_tau_np_p_mlt1_Bplus,
  dtres_Rnp_tau_np_p_xi_mlt_Bplus, dtres_Rnp_tau_np_p_stxi_mlt_Bplus,
  dtres_Rnp_tau_np_n_mlt0_Bplus, dtres_Rnp_tau_np_n_mlt1_Bplus,
  dtres_Rnp_tau_np_n_xi_mlt_Bplus, dtres_Rnp_tau_np_n_stxi_mlt_Bplus,
  dtres_Rol_sig_ol, dtres_Rol_fol_sgl, dtres_Rol_fol_mul, dtres_Dt_cutoff,
};

const int num_dtres_items_fullrecon = sizeof(dtres_items_fullrecon)/sizeof(dtres_syst_item_t);

const dtres_syst_item_t dtres_items_fullrecon_b0[] = {
  dtres_Rdet_Srec0, dtres_Rdet_Srec1,
  dtres_Rdet_Sasc0, dtres_Rdet_Sasc1, dtres_Rnp_Snp   , dtres_Rnp_Snp_global,
  dtres_Rdet_Smn_onetrk, dtres_Rdet_Stl_onetrk, dtres_Rdet_ftl_onetrk,
  dtres_Rdet_rec_mlt_ftail, dtres_Rdet_rec_mlt_stail,
  dtres_Rdet_asc_mlt_ftail, dtres_Rdet_asc_mlt_stail,
  dtres_Rnp_fd_np_sgl0_Bzero, dtres_Rnp_fd_np_sgl1_Bzero, dtres_Rnp_fp_np_sgl_Bzero,
  dtres_Rnp_tau_np_p_sgl0_Bzero, dtres_Rnp_tau_np_p_sgl1_Bzero,
  dtres_Rnp_tau_np_n_sgl0_Bzero, dtres_Rnp_tau_np_n_sgl1_Bzero,
  dtres_Rnp_fd_np_mlt0_Bzero, dtres_Rnp_fd_np_mlt1_Bzero,
  dtres_Rnp_fd_np_st_mlt_Bzero, dtres_Rnp_fd_np_xi_mlt_Bzero, dtres_Rnp_fd_np_stxi_mlt_Bzero,
  dtres_Rnp_fp_np_mlt_Bzero, dtres_Rnp_fn_np_mlt_Bzero,
  dtres_Rnp_tau_np_p_mlt0_Bzero, dtres_Rnp_tau_np_p_mlt1_Bzero,
  dtres_Rnp_tau_np_p_xi_mlt_Bzero, dtres_Rnp_tau_np_p_stxi_mlt_Bzero,
  dtres_Rnp_tau_np_n_mlt0_Bzero, dtres_Rnp_tau_np_n_mlt1_Bzero,
  dtres_Rnp_tau_np_n_xi_mlt_Bzero, dtres_Rnp_tau_np_n_stxi_mlt_Bzero,
  dtres_Rol_sig_ol, dtres_Rol_fol_sgl, dtres_Rol_fol_mul, dtres_Dt_cutoff,
};

const int num_dtres_items_fullrecon_b0 = sizeof(dtres_items_fullrecon_b0)/sizeof(dtres_syst_item_t);

const dtres_syst_item_t dtres_items_fullrecon_bplus[] = {
  dtres_Rdet_Srec0, dtres_Rdet_Srec1,
  dtres_Rdet_Sasc0, dtres_Rdet_Sasc1, dtres_Rnp_Snp   , dtres_Rnp_Snp_global,
  dtres_Rdet_Smn_onetrk, dtres_Rdet_Stl_onetrk, dtres_Rdet_ftl_onetrk,
  dtres_Rdet_rec_mlt_ftail, dtres_Rdet_rec_mlt_stail,
  dtres_Rdet_asc_mlt_ftail, dtres_Rdet_asc_mlt_stail,
  dtres_Rnp_tau_np_p_sgl0_Bplus, dtres_Rnp_tau_np_p_sgl1_Bplus,
  dtres_Rnp_tau_np_n_sgl0_Bplus, dtres_Rnp_tau_np_n_sgl1_Bplus,
  dtres_Rnp_fd_np_mlt0_Bplus, dtres_Rnp_fd_np_mlt1_Bplus,
  dtres_Rnp_fd_np_st_mlt_Bplus, dtres_Rnp_fd_np_xi_mlt_Bplus, dtres_Rnp_fd_np_stxi_mlt_Bplus,
  dtres_Rnp_fp_np_mlt_Bplus, dtres_Rnp_fn_np_mlt_Bplus,
  dtres_Rnp_tau_np_p_mlt0_Bplus, dtres_Rnp_tau_np_p_mlt1_Bplus,
  dtres_Rnp_tau_np_p_xi_mlt_Bplus, dtres_Rnp_tau_np_p_stxi_mlt_Bplus,
  dtres_Rnp_tau_np_n_mlt0_Bplus, dtres_Rnp_tau_np_n_mlt1_Bplus,
  dtres_Rnp_tau_np_n_xi_mlt_Bplus, dtres_Rnp_tau_np_n_stxi_mlt_Bplus,
  dtres_Rol_sig_ol, dtres_Rol_fol_sgl, dtres_Rol_fol_mul, dtres_Dt_cutoff,
};

const int num_dtres_items_fullrecon_bplus = sizeof(dtres_items_fullrecon_bplus)/sizeof(dtres_syst_item_t);

const dtres_syst_item_t dtres_items_partial[] = {
  dtres_Rdet_Srec0, dtres_Rdet_Srec1,
  dtres_Rdet_Sasc0, dtres_Rdet_Sasc1, dtres_Rnp_Snp   , dtres_Rnp_Snp_global,
  dtres_Rdet_Smn_onetrk, dtres_Rdet_Stl_onetrk, dtres_Rdet_ftl_onetrk,
  dtres_Rdet_rec_mlt_ftail, dtres_Rdet_rec_mlt_stail,
  dtres_Rdet_asc_mlt_ftail, dtres_Rdet_asc_mlt_stail,
  dtres_Rnp_fd_np_sgl0_Bzero, dtres_Rnp_fd_np_sgl1_Bzero, dtres_Rnp_fp_np_sgl_Bzero,
  dtres_Rnp_tau_np_p_sgl0_Bzero, dtres_Rnp_tau_np_p_sgl1_Bzero,
  dtres_Rnp_tau_np_n_sgl0_Bzero, dtres_Rnp_tau_np_n_sgl1_Bzero,
  dtres_Rnp_fd_np_mlt0_Bzero, dtres_Rnp_fd_np_mlt1_Bzero,
  dtres_Rnp_fd_np_st_mlt_Bzero, dtres_Rnp_fd_np_xi_mlt_Bzero, dtres_Rnp_fd_np_stxi_mlt_Bzero,
  dtres_Rnp_fp_np_mlt_Bzero, dtres_Rnp_fn_np_mlt_Bzero,
  dtres_Rnp_tau_np_p_mlt0_Bzero, dtres_Rnp_tau_np_p_mlt1_Bzero,
  dtres_Rnp_tau_np_p_xi_mlt_Bzero, dtres_Rnp_tau_np_p_stxi_mlt_Bzero,
  dtres_Rnp_tau_np_n_mlt0_Bzero, dtres_Rnp_tau_np_n_mlt1_Bzero,
  dtres_Rnp_tau_np_n_xi_mlt_Bzero, dtres_Rnp_tau_np_n_stxi_mlt_Bzero,
  dtres_Rnp_fd_np_sgl0_Bplus, dtres_Rnp_fd_np_sgl1_Bplus, dtres_Rnp_fp_np_sgl_Bplus,
  dtres_Rnp_tau_np_p_sgl0_Bplus, dtres_Rnp_tau_np_p_sgl1_Bplus,
  dtres_Rnp_tau_np_n_sgl0_Bplus, dtres_Rnp_tau_np_n_sgl1_Bplus,
  dtres_Rnp_fd_np_mlt0_Bplus, dtres_Rnp_fd_np_mlt1_Bplus,
  dtres_Rnp_fd_np_st_mlt_Bplus, dtres_Rnp_fd_np_xi_mlt_Bplus, dtres_Rnp_fd_np_stxi_mlt_Bplus,
  dtres_Rnp_fp_np_mlt_Bplus, dtres_Rnp_fn_np_mlt_Bplus,
  dtres_Rnp_tau_np_p_mlt0_Bplus, dtres_Rnp_tau_np_p_mlt1_Bplus,
  dtres_Rnp_tau_np_p_xi_mlt_Bplus, dtres_Rnp_tau_np_p_stxi_mlt_Bplus,
  dtres_Rnp_tau_np_n_mlt0_Bplus, dtres_Rnp_tau_np_n_mlt1_Bplus,
  dtres_Rnp_tau_np_n_xi_mlt_Bplus, dtres_Rnp_tau_np_n_stxi_mlt_Bplus,
  dtres_Rk_tau_k_1_p0, dtres_Rk_tau_k_1_p1,
  dtres_Rk_tau_k_2_p0, dtres_Rk_tau_k_2_p1,
  dtres_Rk_tau_k_1_n0, dtres_Rk_tau_k_1_n1,
  dtres_Rk_tau_k_2_n0, dtres_Rk_tau_k_2_n1,
  dtres_Rk_sigma_k0, dtres_Rk_sigma_k1, dtres_Rk_fk2, 
  dtres_Rol_sig_ol, dtres_Rol_fol_sgl, dtres_Rol_fol_mul, dtres_Dt_cutoff,
};

const int num_dtres_items_partial = sizeof(dtres_items_partial)/sizeof(dtres_syst_item_t);

const dtres_syst_item_t dtres_items_partial_b0[] = {
  dtres_Rdet_Srec0, dtres_Rdet_Srec1,
  dtres_Rdet_Sasc0, dtres_Rdet_Sasc1, dtres_Rnp_Snp   , dtres_Rnp_Snp_global,
  dtres_Rdet_Smn_onetrk, dtres_Rdet_Stl_onetrk, dtres_Rdet_ftl_onetrk,
  dtres_Rdet_rec_mlt_ftail, dtres_Rdet_rec_mlt_stail,
  dtres_Rdet_asc_mlt_ftail, dtres_Rdet_asc_mlt_stail,
  dtres_Rnp_fd_np_sgl0_Bzero, dtres_Rnp_fd_np_sgl1_Bzero, dtres_Rnp_fp_np_sgl_Bzero,
  dtres_Rnp_tau_np_p_sgl0_Bzero, dtres_Rnp_tau_np_p_sgl1_Bzero,
  dtres_Rnp_tau_np_n_sgl0_Bzero, dtres_Rnp_tau_np_n_sgl1_Bzero,
  dtres_Rnp_fd_np_mlt0_Bzero, dtres_Rnp_fd_np_mlt1_Bzero,
  dtres_Rnp_fd_np_st_mlt_Bzero, dtres_Rnp_fd_np_xi_mlt_Bzero, dtres_Rnp_fd_np_stxi_mlt_Bzero,
  dtres_Rnp_fp_np_mlt_Bzero, dtres_Rnp_fn_np_mlt_Bzero,
  dtres_Rnp_tau_np_p_mlt0_Bzero, dtres_Rnp_tau_np_p_mlt1_Bzero,
  dtres_Rnp_tau_np_p_xi_mlt_Bzero, dtres_Rnp_tau_np_p_stxi_mlt_Bzero,
  dtres_Rnp_tau_np_n_mlt0_Bzero, dtres_Rnp_tau_np_n_mlt1_Bzero,
  dtres_Rnp_tau_np_n_xi_mlt_Bzero, dtres_Rnp_tau_np_n_stxi_mlt_Bzero,
  dtres_Rk_tau_k_1_p0, dtres_Rk_tau_k_1_p1,
  dtres_Rk_tau_k_2_p0, dtres_Rk_tau_k_2_p1,
  dtres_Rk_tau_k_1_n0, dtres_Rk_tau_k_1_n1,
  dtres_Rk_tau_k_2_n0, dtres_Rk_tau_k_2_n1,
  dtres_Rk_sigma_k0, dtres_Rk_sigma_k1, dtres_Rk_fk2, 
  dtres_Rol_sig_ol, dtres_Rol_fol_sgl, dtres_Rol_fol_mul, dtres_Dt_cutoff,
};

const int num_dtres_items_partial_b0 = sizeof(dtres_items_partial_b0)/sizeof(dtres_syst_item_t);

const dtres_syst_item_t dtres_items_partial_bplus[] = {
  dtres_Rdet_Srec0, dtres_Rdet_Srec1,
  dtres_Rdet_Sasc0, dtres_Rdet_Sasc1, dtres_Rnp_Snp   , dtres_Rnp_Snp_global,
  dtres_Rdet_Smn_onetrk, dtres_Rdet_Stl_onetrk, dtres_Rdet_ftl_onetrk,
  dtres_Rdet_rec_mlt_ftail, dtres_Rdet_rec_mlt_stail,
  dtres_Rdet_asc_mlt_ftail, dtres_Rdet_asc_mlt_stail,
  dtres_Rnp_fd_np_sgl0_Bplus, dtres_Rnp_fd_np_sgl1_Bplus, dtres_Rnp_fp_np_sgl_Bplus,
  dtres_Rnp_tau_np_p_sgl0_Bplus, dtres_Rnp_tau_np_p_sgl1_Bplus,
  dtres_Rnp_tau_np_n_sgl0_Bplus, dtres_Rnp_tau_np_n_sgl1_Bplus,
  dtres_Rnp_fd_np_mlt0_Bplus, dtres_Rnp_fd_np_mlt1_Bplus,
  dtres_Rnp_fd_np_st_mlt_Bplus, dtres_Rnp_fd_np_xi_mlt_Bplus, dtres_Rnp_fd_np_stxi_mlt_Bplus,
  dtres_Rnp_fp_np_mlt_Bplus, dtres_Rnp_fn_np_mlt_Bplus,
  dtres_Rnp_tau_np_p_mlt0_Bplus, dtres_Rnp_tau_np_p_mlt1_Bplus,
  dtres_Rnp_tau_np_p_xi_mlt_Bplus, dtres_Rnp_tau_np_p_stxi_mlt_Bplus,
  dtres_Rnp_tau_np_n_mlt0_Bplus, dtres_Rnp_tau_np_n_mlt1_Bplus,
  dtres_Rnp_tau_np_n_xi_mlt_Bplus, dtres_Rnp_tau_np_n_stxi_mlt_Bplus,
  dtres_Rk_tau_k_1_p0, dtres_Rk_tau_k_1_p1,
  dtres_Rk_tau_k_2_p0, dtres_Rk_tau_k_2_p1,
  dtres_Rk_tau_k_1_n0, dtres_Rk_tau_k_1_n1,
  dtres_Rk_tau_k_2_n0, dtres_Rk_tau_k_2_n1,
  dtres_Rk_sigma_k0, dtres_Rk_sigma_k1, dtres_Rk_fk2, 
  dtres_Rol_sig_ol, dtres_Rol_fol_sgl, dtres_Rol_fol_mul, dtres_Dt_cutoff,
};

const int num_dtres_items_partial_bplus = sizeof(dtres_items_partial_bplus)/sizeof(dtres_syst_item_t);

void Calc_AkCk_Mass(const double cos_theta_b, const double Eb_cms,
		    double* ak, double* ck, const double mb, const double beta){
  *ak = Eb_cms/mb;
  const double pb_cms = sqrt(Eb_cms*Eb_cms-mb*mb);
  *ck = pb_cms*cos_theta_b/(beta*mb);
  return;
}

void Calc_AkCk_Flavor(const double cos_theta_b, const double Eb_cms,
		      double* ak, double* ck, const int b_flavor, const double beta){
  const double mb = (b_flavor==0) ? dt_resol_global::mbzero : dt_resol_global::mbplus;
  Calc_AkCk_Mass(cos_theta_b, Eb_cms, ak, ck, mb, beta);
  return;
}

double Add_Outlier(const int expno,
		   const double x,
		   const double Lin,
		   const int ntrk_rec, const int ntrk_asc,
		   const dtres_param_t * const param,
		   const double nLi,
		   const double llmt,
		   const double ulmt,
		   const double alpha){
  const double fol
    = ((expno>29||!is_dtres_SVD1_by2005(param)||ntrk_rec>1)
       &&ntrk_asc>1) ? param->fol_mul : param->fol_sgl;
  const double m = 0.0;
  const double Lol = gaussian(x, m, param->sig_ol);
  const double nLol = norm_gaussian(llmt, ulmt, m, param->sig_ol);
  const double Li = (1.0-fol)*Lin/nLi + fol*alpha*Lol/nLol;
  return Li;
}

double Add_Outlier_with_Bkg(const int expno,
			    const double x, const double fsig,
			    const double Lsig, const double Lbkg,
			    const int ntrk_rec, const int ntrk_asc,
			    const dtres_param_t * const param,
			    const double nLsig, const double nLbkg,
			    const double llmt, const double ulmt,
			    const double alpha, const double beta){
  const double fol
    = ((expno>29||!is_dtres_SVD1_by2005(param)||ntrk_rec>1)
       &&ntrk_asc>1) ? param->fol_mul : param->fol_sgl;
  const double m = 0.0;
  const double Lol = gaussian(x, m, param->sig_ol);
  const double nLol = norm_gaussian(llmt, ulmt, m, param->sig_ol);
  const double Li = (1.0-fol)*(fsig*Lsig/nLsig+(1.0-fsig)*beta*Lbkg/nLbkg)+fol*alpha*Lol/nLol;
  return Li;
}

void dtres_systematics_expno(const dtres_syst_item_t item,
			     const double sigma,
			     const int expno, const int mc, const int version,
			     const double mc_param_factor, const int oldtype){
  dtres_param_t* param = get_dtres_param(expno, mc, version);
  const dtres_param_t* const posi_err = get_dtres_param_posi_error(expno, version);
  const dtres_param_t* const nega_err = get_dtres_param_nega_error(expno, version);
  dtres_systematics(item, sigma, expno, param,
		    posi_err, nega_err, mc_param_factor, oldtype);
 return; 
}

/// hyperbolic sine and cosine term
double HAfRkRdetRnp_fullrec(const double x, const int flavor, 
			    const double tau, const double dg,
			    const double ak, const double ck,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }

  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  
  const double r_ckak = ck/ak;
  const double ndgtau = r_ckak*ndg*tau;
  const double _fact  = 1.0/(1.0-ndgtau*ndgtau);

  const double ntau_n_0 = 1./gamma_p*(ak-ck);
  const double ntau_p_0 = 1./gamma_m*(ak+ck);

  const double ntau_n_1 = 1./gamma_m*(ak-ck);
  const double ntau_p_1 = 1./gamma_p*(ak+ck);

  const double fact_n_0 = 1./ak*(ntau_n_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_p_0 = 1./ak*(ntau_p_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_n_1 = 1./ak*(ntau_n_1)*_fact*(1.0+ndgtau)/2.0;
  const double fact_p_1 = 1./ak*(ntau_p_1)*_fact*(1.0+ndgtau)/2.0;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdetRnp_full_sup(x,
					       fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm)
                         -EfRkRdetRnp_full_sup(x,
					       fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm);
/*  dout(Debugout::INFO,"dt_resolution")<<"HAfRkRdetRnp_fullrec: "
    <<"p1: "<<EfRkRdetRnp_full_sup(x,
				   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
				   fd, fp, tau_np_p, tau_np_n,
				   mu_mm, sigma_mm)
    <<" p2: "<<EfRkRdetRnp_full_sup(x,
				   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
				   fd, fp, tau_np_p, tau_np_n,
				   mu_mm, sigma_mm)
    <<" Li_mm: "<<Li_mm<<std::endl;*/
    
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdetRnp_full_sup(x,
						 fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt)
                           -EfRkRdetRnp_full_sup(x,
						 fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdetRnp_full_sup(x,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm)
                             -EfRkRdetRnp_full_sup(x,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdetRnp_full_sup(x,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt)
                             -EfRkRdetRnp_full_sup(x,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt);
	const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    =  EfRkRdetRnp_full_sup(x,
						  fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						  fd, fp, tau_np_p, tau_np_n,
						  mu_tm, sigma_tm)
                            -EfRkRdetRnp_full_sup(x,
						  fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						  fd, fp, tau_np_p, tau_np_n,
						  mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double HMfRkRdetRnp_fullrec(const double x, const int flavor, 
			    const double tau, const double dg,
			    const double ak, const double ck,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }

  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  
  const double r_ckak = ck/ak;
  const double ndgtau = r_ckak*ndg*tau;
  const double _fact  = 1.0/(1.0-ndgtau*ndgtau);

  const double ntau_n_0 = 1./gamma_p*(ak-ck);
  const double ntau_p_0 = 1./gamma_m*(ak+ck);

  const double ntau_n_1 = 1./gamma_m*(ak-ck);
  const double ntau_p_1 = 1./gamma_p*(ak+ck);

  const double fact_n_0 = 1./ak*(ntau_n_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_p_0 = 1./ak*(ntau_p_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_n_1 = 1./ak*(ntau_n_1)*_fact*(1.0+ndgtau)/2.0;
  const double fact_p_1 = 1./ak*(ntau_p_1)*_fact*(1.0+ndgtau)/2.0;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EfRkRdetRnp_full_sup(x,
					       fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm)
                         +EfRkRdetRnp_full_sup(x,
					       fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EfRkRdetRnp_full_sup(x,
						 fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt)
                           +EfRkRdetRnp_full_sup(x,
						 fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EfRkRdetRnp_full_sup(x,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm)
                             +EfRkRdetRnp_full_sup(x,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EfRkRdetRnp_full_sup(x,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt)
                             +EfRkRdetRnp_full_sup(x,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt);
	const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EfRkRdetRnp_full_sup(x,
						 fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_tm, sigma_tm)
                           +EfRkRdetRnp_full_sup(x,
						 fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}


/// hyperbolic sine and cosine term
double norm_HAfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor, 
				 const double tau, const double dg,
				 const double ak, const double ck,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }

  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  
  const double r_ckak = ck/ak;
  const double ndgtau = r_ckak*ndg*tau;
  const double _fact  = 1.0/(1.0-ndgtau*ndgtau);

  const double ntau_n_0 = 1./gamma_p*(ak-ck);
  const double ntau_p_0 = 1./gamma_m*(ak+ck);

  const double ntau_n_1 = 1./gamma_m*(ak-ck);
  const double ntau_p_1 = 1./gamma_p*(ak+ck);

  const double fact_n_0 = 1./ak*(ntau_n_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_p_0 = 1./ak*(ntau_p_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_n_1 = 1./ak*(ntau_n_1)*_fact*(1.0+ndgtau)/2.0;
  const double fact_p_1 = 1./ak*(ntau_p_1)*_fact*(1.0+ndgtau)/2.0;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdetRnp_full_sup(ll,ul,
					       fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm)
                         -norm_EfRkRdetRnp_full_sup(ll,ul,
					       fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdetRnp_full_sup(ll,ul,
						 fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt)
                           -norm_EfRkRdetRnp_full_sup(ll,ul,
						 fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm)
                             -norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt)
                             -norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt);
	const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    =  norm_EfRkRdetRnp_full_sup(ll,ul,
						  fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						  fd, fp, tau_np_p, tau_np_n,
						  mu_tm, sigma_tm)
                            -norm_EfRkRdetRnp_full_sup(ll,ul,
						  fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						  fd, fp, tau_np_p, tau_np_n,
						  mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double norm_HMfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor, 
				 const double tau, const double dg,
				 const double ak, const double ck,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param){
  if(ak==0.0||ak==ck){
    if(debugout("INFO")) std::printf("[AfRk_fullrec] Invalid ak==% e, ck==% e,"
	   " where they should be ak!=0.0&&ak!=ck.", ak, ck);
    return -DBL_MAX;
  }

  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  
  const double r_ckak = ck/ak;
  const double ndgtau = r_ckak*ndg*tau;
  const double _fact  = 1.0/(1.0-ndgtau*ndgtau);

  const double ntau_n_0 = 1./gamma_p*(ak-ck);
  const double ntau_p_0 = 1./gamma_m*(ak+ck);

  const double ntau_n_1 = 1./gamma_m*(ak-ck);
  const double ntau_p_1 = 1./gamma_p*(ak+ck);

  const double fact_n_0 = 1./ak*(ntau_n_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_p_0 = 1./ak*(ntau_p_0)*_fact*(1.0-ndgtau)/2.0;
  const double fact_n_1 = 1./ak*(ntau_n_1)*_fact*(1.0+ndgtau)/2.0;
  const double fact_p_1 = 1./ak*(ntau_p_1)*_fact*(1.0+ndgtau)/2.0;

  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EfRkRdetRnp_full_sup(ll,ul,
					       fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm)
                         +norm_EfRkRdetRnp_full_sup(ll,ul,
					       fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
					       fd, fp, tau_np_p, tau_np_n,
					       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EfRkRdetRnp_full_sup(ll,ul,
						 fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt)
                           +norm_EfRkRdetRnp_full_sup(ll,ul,
						 fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/
      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm)
                             +norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt)
                             +norm_EfRkRdetRnp_full_sup(ll,ul,
						   fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						   fd, fp, tau_np_p, tau_np_n,
						   mu_tt, sigma_tt);
	const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EfRkRdetRnp_full_sup(ll,ul,
						 fact_n_0, ntau_n_0, fact_p_0, ntau_p_0,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_tm, sigma_tm)
                           +norm_EfRkRdetRnp_full_sup(ll,ul,
						 fact_n_1, ntau_n_1, fact_p_1, ntau_p_1,
						 fd, fp, tau_np_p, tau_np_n,
						 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

inline double EnRkRdetRnp_partial_sup(const double x, const double tau_n,
				      const double tau_k_1_p, const double tau_k_2_p,
				      const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				      const double tau_np_p, const double tau_np_n,
				      const double fEn, const double fEp,
				      const double fxEn, const double fxEp, const double fxxEn, const double fxxEp,
				      const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				      const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
				      const double fEn_np, const double fEp_np, 
				      const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEn!=0.0)     Li += fEn     *   En_conv_gauss(x, tau_n,     mu_det, nsigma);
  if(fxEn!=0.0)    Li += fxEn    *  xEn_conv_gauss(x, tau_n,     mu_det, nsigma);
  if(fxxEn!=0.0)   Li += fxxEn   * xxEn_conv_gauss(x, tau_n,     mu_det, nsigma);
  if(fEn_k1!=0.0)  Li += fEn_k1  *   En_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)  Li += fEp_k1  *   Ep_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)  Li += fEn_k2  *   En_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)  Li += fEp_k2  *   Ep_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0) Li += fxEn_k1 *  xEn_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0) Li += fxEp_k1 *  xEp_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0) Li += fxEn_k2 *  xEn_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0) Li += fxEp_k2 *  xEp_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)  Li += fEn_np  *   En_conv_gauss(x, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)  Li += fEp_np  *   Ep_conv_gauss(x, tau_np_p,  mu_det, nsigma);

  if(fEp!=0.0||fxEp!=0.0||fxxEp!=0.0){
    dout(Debugout::ERR,"dt_resolution")<<"tatami::EnRkRdetRnp_partial_sup()>> TATAMI INTERNAL ERROR !! Ep fraction!=0 !!"<<std::endl;
    dout(Debugout::ERR,"dt_resolution")<<"tatami::EnRkRdetRnp_partial_sup()>> fEp: "<<fEp<<" fxEp: "<<fxEp<<" fxxEp: "<<fxxEp<<std::endl;
    exit(1);
  }
  return Li;
}

inline double EpRkRdetRnp_partial_sup(const double x, const double tau_p,
				      const double tau_k_1_p, const double tau_k_2_p, 
				      const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
				      const double tau_np_p, const double tau_np_n,
				      const double fEn, const double fEp,
				      const double fxEn, const double fxEp, const double fxxEn, const double fxxEp,
				      const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
				      const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
				      const double fEn_np, const double fEp_np, 
				      const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEp!=0.0)     Li += fEp     *   Ep_conv_gauss(x, tau_p,     mu_det, nsigma);
  if(fxEp!=0.0)    Li += fxEp    *  xEp_conv_gauss(x, tau_p,     mu_det, nsigma);
  if(fxxEp!=0.0)   Li += fxxEp   * xxEp_conv_gauss(x, tau_p,     mu_det, nsigma);
  if(fEn_k1!=0.0)  Li += fEn_k1  *   En_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)  Li += fEp_k1  *   Ep_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)  Li += fEn_k2  *   En_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)  Li += fEp_k2  *   Ep_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0) Li += fxEn_k1 *  xEn_conv_gauss(x, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0) Li += fxEp_k1 *  xEp_conv_gauss(x, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0) Li += fxEn_k2 *  xEn_conv_gauss(x, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0) Li += fxEp_k2 *  xEp_conv_gauss(x, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)  Li += fEn_np  *   En_conv_gauss(x, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)  Li += fEp_np  *   Ep_conv_gauss(x, tau_np_p,  mu_det, nsigma);

  if(fEn!=0.0||fxEn!=0.0||fxxEn!=0.0){
    dout(Debugout::ERR,"dt_resolution")<<"tatami::EpRkRdetRnp_partial_sup()>> TATAMI INTERNAL ERROR !! En fraction!=0 !!"<<std::endl;
    dout(Debugout::ERR,"dt_resolution")<<"tatami::EpRkRdetRnp_partial_sup()>> fEn: "<<fEn<<" fxEn: "<<fxEn<<" fxxEn: "<<fxxEn<<std::endl;
    exit(1);
  }
  return Li;
}

double EnRkRdetRnp_partial_sup0(const double x, const int flavor,  const double tau_n,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param)
{
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;
  
  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0, fxxEn=0.0, fxxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;

  //void add_EnEnp_coef(double* fEn1, double* fEn2, double* fEp2, double* fxEn1,
  //			   const double tau1, const double tau2_n, const double tau2_p,

  add_EnEnp_coef(&fEn, &fEn_k1, &fEp_k1, &fxEn,
		 tau_n, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_EnEnp_coef(&fEn, &fEn_k2, &fEp_k2, &fxEn,
		 tau_n, tau_k_2_n, tau_k_2_p, fk2*fd);
  
  // inline void add_EnEnpEn_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3,
  //			     double* fxEn1, double* fxEn2, double* fxxEn1,

  add_EnEnpEn_coef(&fEn, &fEn_k1, &fEp_k1, &fEn_np,
		   &fxEn, &fxEn_k1, &fxxEn,
		   tau_n, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_EnEnpEn_coef(&fEn, &fEn_k2, &fEp_k2, &fEn_np,
		   &fxEn, &fxEn_k2, &fxxEn,
		   tau_n, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  // void add_EnEnpEp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEp3,
  //			     double* fxEn1, double* fxEp2,
  add_EnEnpEp_coef(&fEn, &fEn_k1, &fEp_k1, &fEp_np,
		   &fxEn, &fxEp_k1,
		   tau_n, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_EnEnpEp_coef(&fEn, &fEn_k2, &fEp_k2, &fEp_np,
		   &fxEn, &fxEp_k2,
		   tau_n, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EnRkRdetRnp_partial_sup(x, tau_n,
						  tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						  tau_np_p, tau_np_n,
						  fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						  fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						  fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						  fEn_np, fEp_np, 
						  mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fEn_tl=0.0, fEp_tl=0.0, fxEn_tl=0.0, fxEp_tl=0.0, fxxEn_tl=0.0, fxxEp_tl=0.0; 
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;

    //void add_EnEnp_coef(double* fEn1, double* fEn2, double* fEp2, double* fxEn1,
    //			   const double tau1, const double tau2_n, const double tau2_p,
    add_EnEnp_coef(&fEn_tl, &fEn_k1_tl, &fEp_k1_tl, &fxEn_tl, 
		   tau_n, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_EnEnp_coef(&fEn_tl, &fEn_k2_tl, &fEp_k2_tl, &fxEn_tl, 
		   tau_n, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    // inline void add_EnEnpEn_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3,
    //			     double* fxEn1, double* fxEn2, double* fxxEn1,
    add_EnEnpEn_coef(&fEn_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEn_k1_tl, &fxxEn_tl,
		     tau_n, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_EnEnpEn_coef(&fEn_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEn_k2_tl, &fxxEn_tl,
		     tau_n, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    // void add_EnEnpEp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEp3,
    //			     double* fxEn1, double* fxEp2,
    add_EnEnpEp_coef(&fEn_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_k1_tl,
		     tau_n, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_EnEnpEp_coef(&fEn_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_k2_tl,
		     tau_n, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
      
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EnRkRdetRnp_partial_sup(x, tau_n,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p_tl, tau_np_n_tl,
						    fEn_tl, fEp_tl, fxEn_tl, fxEp_tl, fxxEn_tl, fxxEp_tl,
						    fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
						    fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl,
						    fEn_np_tl, fEp_np_tl, 
						    mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/

      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EnRkRdetRnp_partial_sup(x, tau_n,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p,  tau_np_n, 
						      fEn,  fEp,  fxEn,  fxEp,  fxxEn,  fxxEp, 
						      fEn_k1,  fEp_k1,  fEn_k2,  fEp_k2, 
						      fxEn_k1,  fxEp_k1,  fxEn_k2,  fxEp_k2, 
						      fEn_np,  fEp_np,  
						      mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EnRkRdetRnp_partial_sup(x, tau_n,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p_tl,  tau_np_n_tl, 
						      fEn_tl,  fEp_tl,  fxEn_tl,  fxEp_tl,  fxxEn_tl,  fxxEp_tl, 
						      fEn_k1_tl,  fEp_k1_tl,  fEn_k2_tl,  fEp_k2_tl, 
						      fxEn_k1_tl,  fxEp_k1_tl,  fxEn_k2_tl,  fxEp_k2_tl, 
						      fEn_np_tl,  fEp_np_tl,  
						      mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EnRkRdetRnp_partial_sup(x, tau_n,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p, tau_np_n,
						    fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						    fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						    fEn_np, fEp_np, 
						    mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double EpRkRdetRnp_partial_sup0(const double x, const int flavor,  const double tau_p,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param)
{
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;
  
  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0, fxxEn=0.0, fxxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;

  add_EpEnp_coef(&fEp, &fEn_k1, &fEp_k1, &fxEp,
		 tau_p, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_EpEnp_coef(&fEp, &fEn_k2, &fEp_k2, &fxEp,
		 tau_p, tau_k_2_n, tau_k_2_p, fk2*fd);
  
  add_EpEnpEn_coef(&fEp, &fEn_k1, &fEp_k1, &fEn_np,
		   &fxEp, &fxEn_k1,
		   tau_p, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_EpEnpEn_coef(&fEp, &fEn_k2, &fEp_k2, &fEn_np,
		   &fxEp, &fxEn_k2,
		   tau_p, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_EpEnpEp_coef(&fEp, &fEn_k1, &fEp_k1, &fEp_np,
		   &fxEp, &fxEp_k1, &fxxEp,
		   tau_p, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_EpEnpEp_coef(&fEp, &fEn_k2, &fEp_k2, &fEp_np,
		   &fxEp, &fxEp_k2, &fxxEp,
		   tau_p, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = EpRkRdetRnp_partial_sup(x, tau_p,
						  tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						  tau_np_p, tau_np_n,
						  fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						  fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						  fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						  fEn_np, fEp_np, 
						  mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fEn_tl=0.0, fEp_tl=0.0, fxEn_tl=0.0, fxEp_tl=0.0, fxxEn_tl=0.0, fxxEp_tl=0.0; 
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
    
    add_EpEnp_coef(&fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fxEp_tl,
		   tau_p, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_EpEnp_coef(&fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fxEp_tl,
		   tau_p, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_EpEnpEn_coef(&fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl,
		     &fxEp_tl, &fxEn_k1_tl,
		     tau_p, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_EpEnpEn_coef(&fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl,
		     &fxEp_tl, &fxEn_k2_tl,
		     tau_p, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_EpEnpEp_coef(&fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl,
		     &fxEp_tl, &fxEp_k1_tl, &fxxEp_tl,
		     tau_p, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_EpEnpEp_coef(&fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl,
		     &fxEp_tl, &fxEp_k2_tl, &fxxEp_tl,
		     tau_p, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
      
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = EpRkRdetRnp_partial_sup(x, tau_p,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p_tl, tau_np_n_tl,
						    fEn_tl, fEp_tl, fxEn_tl, fxEp_tl, fxxEn_tl, fxxEp_tl,
						    fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
						    fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl,
						    fEn_np_tl, fEp_np_tl, 
						    mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/

      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = EpRkRdetRnp_partial_sup(x, tau_p,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p,  tau_np_n, 
						      fEn,  fEp,  fxEn,  fxEp,  fxxEn,  fxxEp, 
						      fEn_k1,  fEp_k1,  fEn_k2,  fEp_k2, 
						      fxEn_k1,  fxEp_k1,  fxEn_k2,  fxEp_k2, 
						      fEn_np,  fEp_np,  
						      mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = EpRkRdetRnp_partial_sup(x, tau_p,
						      tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						      tau_np_p_tl,  tau_np_n_tl, 
						      fEn_tl,  fEp_tl,  fxEn_tl,  fxEp_tl,  fxxEn_tl,  fxxEp_tl, 
						      fEn_k1_tl,  fEp_k1_tl,  fEn_k2_tl,  fEp_k2_tl, 
						      fxEn_k1_tl,  fxEp_k1_tl,  fxEn_k2_tl,  fxEp_k2_tl, 
						      fEn_np_tl,  fEp_np_tl,  
						      mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = EpRkRdetRnp_partial_sup(x, tau_p,
						    tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						    tau_np_p, tau_np_n,
						    fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						    fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						    fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						    fEn_np, fEp_np, 
						    mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double HMfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dg,
			    const double dz,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param)
{
  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  const double tau_n_0 = 1./gamma_p;
  const double tau_p_0 = 1./gamma_m;

  const double tau_n_1 = 1./gamma_m;
  const double tau_p_1 = 1./gamma_p;

  const double L_n_0=EnRkRdetRnp_partial_sup0(x, flavor, tau_n_0,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  const double L_p_0=EpRkRdetRnp_partial_sup0(x, flavor, tau_p_0,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  const double L_n_1=EnRkRdetRnp_partial_sup0(x, flavor, tau_n_1,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  const double L_p_1=EpRkRdetRnp_partial_sup0(x, flavor, tau_p_1,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  // cancel their normalization factors
  const double L_0 = tau_n_0 * L_n_0 + tau_p_0 * L_p_0;
  const double L_1 = tau_n_1 * L_n_1 + tau_p_1 * L_p_1;

  return (L_0+L_1)/2.;
}

double HAfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dg,
			    const double dz,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param)
{
  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  const double tau_n_0 = 1./gamma_p;
  const double tau_p_0 = 1./gamma_m;

  const double tau_n_1 = 1./gamma_m;
  const double tau_p_1 = 1./gamma_p;

  const double L_n_0=EnRkRdetRnp_partial_sup0(x, flavor, tau_n_0,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  const double L_p_0=EpRkRdetRnp_partial_sup0(x, flavor, tau_p_0,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  const double L_n_1=EnRkRdetRnp_partial_sup0(x, flavor, tau_n_1,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  const double L_p_1=EpRkRdetRnp_partial_sup0(x, flavor, tau_p_1,
					      dz,
					      ntrk_rec, sz_rec,
					      chisq_z_rec, ndf_z_rec,
					      ntrk_asc, sz_asc,
					      chisq_z_asc, ndf_z_asc, keeptagl,
					      param);
  // cancel their normalization factors
  const double L_0 = tau_n_0 * L_n_0 + tau_p_0 * L_p_0;
  const double L_1 = tau_n_1 * L_n_1 + tau_p_1 * L_p_1;

  return (L_0-L_1)/2.;
}

  
inline double norm_EnRkRdetRnp_partial_sup(const double ll, const double ul, const double tau_n,
					   const double tau_k_1_p, const double tau_k_2_p, 
					   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					   const double tau_np_p, const double tau_np_n,
					   const double fEn, const double fEp,
					   const double fxEn, const double fxEp, const double fxxEn, const double fxxEp,
					   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					   const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
					   const double fEn_np, const double fEp_np, 
					   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEn!=0.0)     Li += fEn     *   norm_En_conv_gauss(ll, ul, tau_n,     mu_det, nsigma);
  if(fxEn!=0.0)    Li += fxEn    *  norm_xEn_conv_gauss(ll, ul, tau_n,     mu_det, nsigma);
  if(fxxEn!=0.0)   Li += fxxEn   * norm_xxEn_conv_gauss(ll, ul, tau_n,     mu_det, nsigma);
  if(fEn_k1!=0.0)  Li += fEn_k1  *   norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)  Li += fEp_k1  *   norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)  Li += fEn_k2  *   norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)  Li += fEp_k2  *   norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0) Li += fxEn_k1 *  norm_xEn_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0) Li += fxEp_k1 *  norm_xEp_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0) Li += fxEn_k2 *  norm_xEn_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0) Li += fxEp_k2 *  norm_xEp_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)  Li += fEn_np  *   norm_En_conv_gauss(ll, ul, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)  Li += fEp_np  *   norm_Ep_conv_gauss(ll, ul, tau_np_p,  mu_det, nsigma);

  if(fEp!=0.0||fxEp!=0.0||fxxEp!=0.0){
    dout(Debugout::ERR,"dt_resolution")<<"tatami::norm_EnRkRdetRnp_partial_sup()>> TATAMI INTERNAL ERROR !! Ep fraction!=0 !!"<<std::endl;
    dout(Debugout::ERR,"dt_resolution")<<"tatami::norm_EnRkRdetRnp_partial_sup()>> fEp: "<<fEp<<" fxEp: "<<fxEp<<" fxxEp: "<<fxxEp<<std::endl;
    exit(1);
  }
  return Li;
}

inline double norm_EpRkRdetRnp_partial_sup(const double ll, const double ul, const double tau_p,
					   const double tau_k_1_p, const double tau_k_2_p, 
					   const double tau_k_1_n, const double tau_k_2_n, const double sigma_k,
					   const double tau_np_p, const double tau_np_n,
					   const double fEn, const double fEp,
					   const double fxEn, const double fxEp, const double fxxEn, const double fxxEp,
					   const double fEn_k1, const double fEp_k1, const double fEn_k2, const double fEp_k2,
					   const double fxEn_k1, const double fxEp_k1, const double fxEn_k2, const double fxEp_k2,
					   const double fEn_np, const double fEp_np, 
					   const double mu_det, const double sigma_det){
  double Li = 0.0;
  const double nsigma = sum_sigma(sigma_k, sigma_det);
  if(fEp!=0.0)     Li += fEp     *   norm_Ep_conv_gauss(ll, ul, tau_p,     mu_det, nsigma);
  if(fxEp!=0.0)    Li += fxEp    *  norm_xEp_conv_gauss(ll, ul, tau_p,     mu_det, nsigma);
  if(fxxEp!=0.0)   Li += fxxEp   * norm_xxEp_conv_gauss(ll, ul, tau_p,     mu_det, nsigma);
  if(fEn_k1!=0.0)  Li += fEn_k1  *   norm_En_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fEp_k1!=0.0)  Li += fEp_k1  *   norm_Ep_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fEn_k2!=0.0)  Li += fEn_k2  *   norm_En_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fEp_k2!=0.0)  Li += fEp_k2  *   norm_Ep_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fxEn_k1!=0.0) Li += fxEn_k1 *  norm_xEn_conv_gauss(ll, ul, tau_k_1_n, mu_det, nsigma);
  if(fxEp_k1!=0.0) Li += fxEp_k1 *  norm_xEp_conv_gauss(ll, ul, tau_k_1_p, mu_det, nsigma);
  if(fxEn_k2!=0.0) Li += fxEn_k2 *  norm_xEn_conv_gauss(ll, ul, tau_k_2_n, mu_det, nsigma);
  if(fxEp_k2!=0.0) Li += fxEp_k2 *  norm_xEp_conv_gauss(ll, ul, tau_k_2_p, mu_det, nsigma);
  if(fEn_np!=0.0)  Li += fEn_np  *   norm_En_conv_gauss(ll, ul, tau_np_n,  mu_det, nsigma);
  if(fEp_np!=0.0)  Li += fEp_np  *   norm_Ep_conv_gauss(ll, ul, tau_np_p,  mu_det, nsigma);

  if(fEn!=0.0||fxEn!=0.0||fxxEn!=0.0){
    dout(Debugout::ERR,"dt_resolution")<<"tatami::norm_EpRkRdetRnp_partial_sup()>> TATAMI INTERNAL ERROR !! En fraction!=0 !!"<<std::endl;
    dout(Debugout::ERR,"dt_resolution")<<"tatami::norm_EpRkRdetRnp_partial_sup()>> fEn: "<<fEn<<" fxEn: "<<fxEn<<" fxxEn: "<<fxxEn<<std::endl;
    exit(1);
  }

  return Li;
}

double norm_EnRkRdetRnp_partial_sup0(const double ll, const double ul, const int flavor,  const double tau_n,
				     const double dz,
				     const int ntrk_rec, const double sz_rec,
				     const double chisq_z_rec, const int ndf_z_rec,
				     const int ntrk_asc, const double sz_asc,
				     const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				     const dtres_param_t * const param)
{
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;
  
  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0, fxxEn=0.0, fxxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;

  //void add_EnEnp_coef(double* fEn1, double* fEn2, double* fEp2, double* fxEn1,
  //			   const double tau1, const double tau2_n, const double tau2_p,

  add_EnEnp_coef(&fEn, &fEn_k1, &fEp_k1, &fxEn,
		 tau_n, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_EnEnp_coef(&fEn, &fEn_k2, &fEp_k2, &fxEn,
		 tau_n, tau_k_2_n, tau_k_2_p, fk2*fd);
  
  // inline void add_EnEnpEn_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3,
  //			     double* fxEn1, double* fxEn2, double* fxxEn1,

  add_EnEnpEn_coef(&fEn, &fEn_k1, &fEp_k1, &fEn_np,
		   &fxEn, &fxEn_k1, &fxxEn,
		   tau_n, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_EnEnpEn_coef(&fEn, &fEn_k2, &fEp_k2, &fEn_np,
		   &fxEn, &fxEn_k2, &fxxEn,
		   tau_n, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  // void add_EnEnpEp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEp3,
  //			     double* fxEn1, double* fxEp2,
  add_EnEnpEp_coef(&fEn, &fEn_k1, &fEp_k1, &fEp_np,
		   &fxEn, &fxEp_k1,
		   tau_n, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_EnEnpEp_coef(&fEn, &fEn_k2, &fEp_k2, &fEp_np,
		   &fxEn, &fxEp_k2,
		   tau_n, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EnRkRdetRnp_partial_sup(ll, ul, tau_n,
						       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						       tau_np_p, tau_np_n,
						       fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						       fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						       fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						       fEn_np, fEp_np, 
						       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fEn_tl=0.0, fEp_tl=0.0, fxEn_tl=0.0, fxEp_tl=0.0, fxxEn_tl=0.0, fxxEp_tl=0.0; 
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;

    //void add_EnEnp_coef(double* fEn1, double* fEn2, double* fEp2, double* fxEn1,
    //			   const double tau1, const double tau2_n, const double tau2_p,
    add_EnEnp_coef(&fEn_tl, &fEn_k1_tl, &fEp_k1_tl, &fxEn_tl, 
		   tau_n, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_EnEnp_coef(&fEn_tl, &fEn_k2_tl, &fEp_k2_tl, &fxEn_tl, 
		   tau_n, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    // inline void add_EnEnpEn_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3,
    //			     double* fxEn1, double* fxEn2, double* fxxEn1,
    add_EnEnpEn_coef(&fEn_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEn_k1_tl, &fxxEn_tl,
		     tau_n, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_EnEnpEn_coef(&fEn_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl,
		     &fxEn_tl, &fxEn_k2_tl, &fxxEn_tl,
		     tau_n, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    // void add_EnEnpEp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEp3,
    //			     double* fxEn1, double* fxEp2,
    add_EnEnpEp_coef(&fEn_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_k1_tl,
		     tau_n, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_EnEnpEp_coef(&fEn_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl,
		     &fxEn_tl, &fxEp_k2_tl,
		     tau_n, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
      
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EnRkRdetRnp_partial_sup(ll, ul, tau_n,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p_tl, tau_np_n_tl,
							 fEn_tl, fEp_tl, fxEn_tl, fxEp_tl, fxxEn_tl, fxxEp_tl,
							 fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
							 fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl,
							 fEn_np_tl, fEp_np_tl, 
							 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/

      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EnRkRdetRnp_partial_sup(ll, ul, tau_n,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p,  tau_np_n, 
							   fEn,  fEp,  fxEn,  fxEp,  fxxEn,  fxxEp, 
							   fEn_k1,  fEp_k1,  fEn_k2,  fEp_k2, 
							   fxEn_k1,  fxEp_k1,  fxEn_k2,  fxEp_k2, 
							   fEn_np,  fEp_np,  
							   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EnRkRdetRnp_partial_sup(ll, ul, tau_n,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p_tl,  tau_np_n_tl, 
							   fEn_tl,  fEp_tl,  fxEn_tl,  fxEp_tl,  fxxEn_tl,  fxxEp_tl, 
							   fEn_k1_tl,  fEp_k1_tl,  fEn_k2_tl,  fEp_k2_tl, 
							   fxEn_k1_tl,  fxEp_k1_tl,  fxEn_k2_tl,  fxEp_k2_tl, 
							   fEn_np_tl,  fEp_np_tl,  
						      mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EnRkRdetRnp_partial_sup(ll, ul, tau_n,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p, tau_np_n,
							 fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
							 fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							 fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
							 fEn_np, fEp_np, 
							 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double norm_EpRkRdetRnp_partial_sup0(const double ll, const double ul, const int flavor,  const double tau_p,
				     const double dz,
				     const int ntrk_rec, const double sz_rec,
				     const double chisq_z_rec, const int ndf_z_rec,
				     const int ntrk_asc, const double sz_asc,
				     const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				     const dtres_param_t * const param)
{
  double xi_rec, st_rec, ftail_rec, mu_main_rec, Smain_rec, mu_tail_rec, Stail_rec;
  calc_vtxparam_rec(ntrk_rec, sz_rec, chisq_z_rec, ndf_z_rec, &xi_rec, &st_rec);
  Rrec_param(ntrk_rec, xi_rec, st_rec, param,
	     &ftail_rec, &mu_main_rec, &Smain_rec, &mu_tail_rec, &Stail_rec);

  double xi_asc, st_asc, ftail_asc, mu_main_asc, Smain_asc, mu_tail_asc, Stail_asc;
  calc_vtxparam_asc(ntrk_asc, sz_asc, chisq_z_asc, ndf_z_asc, &xi_asc, &st_asc);
  Rasc_param(ntrk_asc, xi_asc, st_asc, param,
	     &ftail_asc, &mu_main_asc, &Smain_asc, &mu_tail_asc, &Stail_asc);

  double fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl;
  Rnp_param(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc, Stail_asc, param,
	    &fd, &fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  swap_rnp_param(&fp, &tau_np_p, &tau_np_n, &tau_np_p_tl, &tau_np_n_tl);

  const double fe = 1.0 - fd;
  const double nfn = (1.0 - fp)*fe;
  const double nfp = fp*fe;

  double tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k, fk2;
  Rk_partial_param(dz, param,
		   &tau_k_1_p, &tau_k_2_p, &tau_k_1_n, &tau_k_2_n, &sigma_k, &fk2);
  const double fk1 = 1.0 - fk2;
  
  double fEn=0.0, fEp=0.0, fxEn=0.0, fxEp=0.0, fxxEn=0.0, fxxEp=0.0; 
  double fEn_k1=0.0, fEp_k1=0.0, fEn_k2=0.0, fEp_k2=0.0;
  double fxEn_k1=0.0, fxEp_k1=0.0, fxEn_k2=0.0, fxEp_k2=0.0;
  double fEn_np=0.0, fEp_np=0.0;

  add_EpEnp_coef(&fEp, &fEn_k1, &fEp_k1, &fxEp,
		 tau_p, tau_k_1_n, tau_k_1_p, fk1*fd);
  add_EpEnp_coef(&fEp, &fEn_k2, &fEp_k2, &fxEp,
		 tau_p, tau_k_2_n, tau_k_2_p, fk2*fd);
  
  add_EpEnpEn_coef(&fEp, &fEn_k1, &fEp_k1, &fEn_np,
		   &fxEp, &fxEn_k1,
		   tau_p, tau_k_1_n, tau_k_1_p, tau_np_n, fk1*nfn);
  add_EpEnpEn_coef(&fEp, &fEn_k2, &fEp_k2, &fEn_np,
		   &fxEp, &fxEn_k2,
		   tau_p, tau_k_2_n, tau_k_2_p, tau_np_n, fk2*nfn);
  add_EpEnpEp_coef(&fEp, &fEn_k1, &fEp_k1, &fEp_np,
		   &fxEp, &fxEp_k1, &fxxEp,
		   tau_p, tau_k_1_n, tau_k_1_p, tau_np_p, fk1*nfp);
  add_EpEnpEp_coef(&fEp, &fEn_k2, &fEp_k2, &fEp_np,
		   &fxEp, &fxEp_k2, &fxxEp,
		   tau_p, tau_k_2_n, tau_k_2_p, tau_np_p, fk2*nfp);

  
  const double mu_mm    = mu_main_rec+mu_main_asc;
  const double sigma_mm = sum_sigma(Smain_rec, Smain_asc);
  const double Li_mm    = norm_EpRkRdetRnp_partial_sup(ll, ul, tau_p,
						       tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
						       tau_np_p, tau_np_n,
						       fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
						       fEn_k1, fEp_k1, fEn_k2, fEp_k2,
						       fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
						       fEn_np, fEp_np, 
						       mu_mm, sigma_mm);
  if(ftail_asc>0.0){
    double fEn_tl=0.0, fEp_tl=0.0, fxEn_tl=0.0, fxEp_tl=0.0, fxxEn_tl=0.0, fxxEp_tl=0.0; 
    double fEn_k1_tl=0.0, fEp_k1_tl=0.0, fEn_k2_tl=0.0, fEp_k2_tl=0.0;
    double fxEn_k1_tl=0.0, fxEp_k1_tl=0.0, fxEn_k2_tl=0.0, fxEp_k2_tl=0.0;
    double fEn_np_tl=0.0, fEp_np_tl=0.0;
    
    add_EpEnp_coef(&fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fxEp_tl,
		   tau_p, tau_k_1_n, tau_k_1_p, fk1*fd);
    add_EpEnp_coef(&fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fxEp_tl,
		   tau_p, tau_k_2_n, tau_k_2_p, fk2*fd);
    
    add_EpEnpEn_coef(&fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEn_np_tl,
		     &fxEp_tl, &fxEn_k1_tl,
		     tau_p, tau_k_1_n, tau_k_1_p, tau_np_n_tl, fk1*nfn);
    add_EpEnpEn_coef(&fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEn_np_tl,
		     &fxEp_tl, &fxEn_k2_tl,
		     tau_p, tau_k_2_n, tau_k_2_p, tau_np_n_tl, fk2*nfn);
    add_EpEnpEp_coef(&fEp_tl, &fEn_k1_tl, &fEp_k1_tl, &fEp_np_tl,
		     &fxEp_tl, &fxEp_k1_tl, &fxxEp_tl,
		     tau_p, tau_k_1_n, tau_k_1_p, tau_np_p_tl, fk1*nfp);
    add_EpEnpEp_coef(&fEp_tl, &fEn_k2_tl, &fEp_k2_tl, &fEp_np_tl,
		     &fxEp_tl, &fxEp_k2_tl, &fxxEp_tl,
		     tau_p, tau_k_2_n, tau_k_2_p, tau_np_p_tl, fk2*nfp);
      
    
    const double mu_mt    = mu_main_rec+mu_tail_asc;
    const double sigma_mt = sum_sigma(Smain_rec, Stail_asc);
    const double Li_mt    = norm_EpRkRdetRnp_partial_sup(ll, ul, tau_p,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p_tl, tau_np_n_tl,
							 fEn_tl, fEp_tl, fxEn_tl, fxEp_tl, fxxEn_tl, fxxEp_tl,
							 fEn_k1_tl, fEp_k1_tl, fEn_k2_tl, fEp_k2_tl,
							 fxEn_k1_tl, fxEp_k1_tl, fxEn_k2_tl, fxEp_k2_tl,
							 fEn_np_tl, fEp_np_tl, 
							 mu_mt, sigma_mt);
    if(ftail_rec>0.0){ /* single track track (Asc) && single track vertex (rec)*/

      const double mu_tm    = mu_tail_rec+mu_main_asc;
      const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
      const double Li_tm    = norm_EpRkRdetRnp_partial_sup(ll, ul, tau_p,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p,  tau_np_n, 
							   fEn,  fEp,  fxEn,  fxEp,  fxxEn,  fxxEp, 
							   fEn_k1,  fEp_k1,  fEn_k2,  fEp_k2, 
							   fxEn_k1,  fxEp_k1,  fxEn_k2,  fxEp_k2, 
							   fEn_np,  fEp_np,  
							   mu_tm, sigma_tm);
      const double mu_tt    = mu_tail_rec+mu_tail_asc;
      const double sigma_tt = sum_sigma(Stail_rec, Stail_asc);
      const double Li_tt    = norm_EpRkRdetRnp_partial_sup(ll, ul, tau_p,
							   tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							   tau_np_p_tl,  tau_np_n_tl, 
							   fEn_tl,  fEp_tl,  fxEn_tl,  fxEp_tl,  fxxEn_tl,  fxxEp_tl, 
							   fEn_k1_tl,  fEp_k1_tl,  fEn_k2_tl,  fEp_k2_tl, 
							   fxEn_k1_tl,  fxEp_k1_tl,  fxEn_k2_tl,  fxEp_k2_tl, 
							   fEn_np_tl,  fEp_np_tl,  
							   mu_tt, sigma_tt);
      const double Li = ((1.0-ftail_rec)*(1.0-ftail_asc)*Li_mm + ftail_rec*(1.0-ftail_asc)*Li_tm
			 +(1.0-ftail_rec)*ftail_asc*Li_mt + ftail_rec*ftail_asc*Li_tt);
      return Li;
    }else{ /* single track track (Asc) && multiple track vertex (rec)*/
      const double Li = (1.0-ftail_asc)*Li_mm + ftail_asc*Li_mt;
      return Li;
    }
  }else if(ftail_rec>0.0){ /* multiple track track (Asc) && single track vertex (rec)*/
    const double mu_tm    = mu_tail_rec+mu_main_asc;
    const double sigma_tm = sum_sigma(Stail_rec, Smain_asc);
    const double Li_tm    = norm_EpRkRdetRnp_partial_sup(ll, ul, tau_p,
							 tau_k_1_p, tau_k_2_p, tau_k_1_n, tau_k_2_n, sigma_k,
							 tau_np_p, tau_np_n,
							 fEn, fEp, fxEn, fxEp, fxxEn, fxxEp,
							 fEn_k1, fEp_k1, fEn_k2, fEp_k2,
							 fxEn_k1, fxEp_k1, fxEn_k2, fxEp_k2,
							 fEn_np, fEp_np, 
							 mu_tm, sigma_tm);
    const double Li = (1.0-ftail_rec)*Li_mm + ftail_rec*Li_tm;
    return Li;
  }
  /* multiple track track (Asc) && multiple track vertex (rec)*/
  return Li_mm;
}

double norm_HMfRkRdetRnp_partial(const double ll, const double ul, const int flavor,
				 const double tau, const double dg,
				 const double dz,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param)
{
  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  const double tau_n_0 = 1./gamma_p;
  const double tau_p_0 = 1./gamma_m;

  const double tau_n_1 = 1./gamma_m;
  const double tau_p_1 = 1./gamma_p;

  const double L_n_0=norm_EnRkRdetRnp_partial_sup0(ll, ul, flavor, tau_n_0,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  const double L_p_0=norm_EpRkRdetRnp_partial_sup0(ll, ul, flavor, tau_p_0,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  const double L_n_1=norm_EnRkRdetRnp_partial_sup0(ll, ul, flavor, tau_n_1,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  const double L_p_1=norm_EpRkRdetRnp_partial_sup0(ll, ul, flavor, tau_p_1,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  // cancel their normalization factors
  const double L_0 = tau_n_0 * L_n_0 + tau_p_0 * L_p_0;
  const double L_1 = tau_n_1 * L_n_1 + tau_p_1 * L_p_1;

  return (L_0+L_1)/2.;
}

double norm_HAfRkRdetRnp_partial(const double ll, const double ul, const int flavor,
				 const double tau, const double dg,
				 const double dz,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param)
{
  const double gamma=1./tau;
  const double ndg=dg*0.5;
  const double gamma_p=gamma+ndg;
  const double gamma_m=gamma-ndg;
  const double tau_n_0 = 1./gamma_p;
  const double tau_p_0 = 1./gamma_m;

  const double tau_n_1 = 1./gamma_m;
  const double tau_p_1 = 1./gamma_p;

  const double L_n_0=norm_EnRkRdetRnp_partial_sup0(ll, ul, flavor, tau_n_0,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  const double L_p_0=norm_EpRkRdetRnp_partial_sup0(ll, ul, flavor, tau_p_0,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  const double L_n_1=norm_EnRkRdetRnp_partial_sup0(ll, ul, flavor, tau_n_1,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  const double L_p_1=norm_EpRkRdetRnp_partial_sup0(ll, ul, flavor, tau_p_1,
						   dz,
						   ntrk_rec, sz_rec,
						   chisq_z_rec, ndf_z_rec,
						   ntrk_asc, sz_asc,
						   chisq_z_asc, ndf_z_asc, keeptagl,
						   param);
  // cancel their normalization factors
  const double L_0 = tau_n_0 * L_n_0 + tau_p_0 * L_p_0;
  const double L_1 = tau_n_1 * L_n_1 + tau_p_1 * L_p_1;

  return (L_0-L_1)/2.;
}

__CLOSE_LINKAGE__

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
