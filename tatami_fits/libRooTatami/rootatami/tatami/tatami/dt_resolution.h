//
//  dt_resolution.h
//     --- Description: Library for Resolution Function of dt
//   $Id: dt_resolution.h 11208 2011-03-01 14:53:58Z hitoshi $
//
//   $Log$
//   Revision 1.3  2004/10/19 07:09:56  kohji
//   1) SVD1 and SVD2 resolution parameters are updated for 2004 publication.
//   2) "expno" is added to add_outlier and dtres_systematics functions to treat
//   difference between SVD1 and SVD2 properly[cpfit_ml:0771].
//   3) New functions are added for cosh/sinh terms
//    (norm_)HAfRkRdetRnp_{fullrec,partial} for HAf=exp(-|t|Gamma)sinh(DGamma/2 t)
//    (norm_)HMfRkRdetRnp_{fullrec,partial} for HMf=exp(-|t|Gamma)cosh(DGamma/2 t)
//
//   Revision 1.2  2003/08/06 06:43:33  katayama
//   New version from Nakadaira san
//
//   Revision 1.1  2002/09/05 01:27:06  katayama
//   New package tatami from Nakadaira/Sumisawa san
//
//   Revision 1.10  2002/07/12 09:56:37  nakadair
//   Add new functions: EfRdet(...), AfRdet(...), MfRdet(...),
//   norm_EfRdet(...), norm_AfRdet(...) and norm_MfRdet(...).
//
//   Revision 1.9  2002/07/11 09:32:30  nakadair
//   Move parameters constants to dt_resolution_const.{h,cc}.
//   Add new function Add(_)Outlier(_)With(_)Bkg(...).
//
//   Revision 1.8  2002/07/10 06:50:24  nakadair
//   Enable systematic study on delta-t cut-off.
//
//   Revision 1.7  2002/07/05 06:16:11  nakadair
//   Add Calc_AkCk_Mass(...), Calc_AkCk_Flavor(...), Calc_AkCk(...), and Add(_)Outlier(...).
//   Add new members and C interfece to class dt_resol_global.
//
//   Revision 1.6  2002/07/02 23:56:34  nakadair
//   Add utility functions for systematics study.
//   Update Rk_parameters.
//
//   Revision 1.5  2002/06/27 13:29:08  nakadair
//   Change Rk_partial_param(...), Rnp_param(...), swap_rnp_param(...),
//   Rasc_param(...) and Rrec_param(...) from inline to global function.
//   Add class dt_resol_global.
//
//   Revision 1.4  2002/06/11 09:36:30  nakadair
//   Add declaration of Rdet(...)
//
//   Revision 1.3  2002/05/09 07:08:23  nakadair
//   Add add_outlier(...).
//
//   Revision 1.2  2002/05/09 06:23:11  nakadair
//   Change void calc_akck(...) from inline to global function.
//
//   Revision 1.1  2002/05/07 21:38:23  nakadair
//   Add convolution with Resolution function and fortran interface.
//
//

#ifndef __DT_RESOLUTION_H__
#define __DT_RESOLUTION_H__

#ifdef __cplusplus
#define __LINKAGE__  extern "C" {
#define __CLOSE_LINKAGE__       }
#else /* __cplusplus */
#define __LINKAGE__  
#define __CLOSE_LINKAGE__
#endif /* __cplusplus */

#include "belle.h"
#include "tatami/dt_resolution_const.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

__LINKAGE__ 

typedef enum dtres_syst_item_t{
  none                              =0,
  dtres_Rdet_Srec0                  ,
  dtres_Rdet_Srec1                  ,
  dtres_Rdet_Sasc0                  ,
  dtres_Rdet_Sasc1                  ,
  dtres_Rnp_Snp                     ,
  dtres_Rnp_Snp_global              ,
  dtres_Rdet_Smn_onetrk             ,
  dtres_Rdet_Stl_onetrk             ,
  dtres_Rdet_ftl_onetrk             ,
  dtres_Rdet_rec_mlt_ftail          ,
  dtres_Rdet_rec_mlt_stail          ,
  dtres_Rdet_asc_mlt_ftail          ,
  dtres_Rdet_asc_mlt_stail          ,
  dtres_Rnp_fd_np_sgl0_Bzero        ,
  dtres_Rnp_fd_np_sgl1_Bzero        ,
  dtres_Rnp_fp_np_sgl_Bzero         ,
  dtres_Rnp_tau_np_p_sgl0_Bzero     ,
  dtres_Rnp_tau_np_p_sgl1_Bzero     ,
  dtres_Rnp_tau_np_n_sgl0_Bzero     ,
  dtres_Rnp_tau_np_n_sgl1_Bzero     ,
  dtres_Rnp_fd_np_mlt0_Bzero        ,
  dtres_Rnp_fd_np_mlt1_Bzero        ,
  dtres_Rnp_fd_np_st_mlt_Bzero      ,
  dtres_Rnp_fd_np_xi_mlt_Bzero      ,
  dtres_Rnp_fd_np_stxi_mlt_Bzero    ,
  dtres_Rnp_fp_np_mlt_Bzero         ,
  dtres_Rnp_fn_np_mlt_Bzero         ,
  dtres_Rnp_tau_np_p_mlt0_Bzero     ,
  dtres_Rnp_tau_np_p_mlt1_Bzero     ,
  dtres_Rnp_tau_np_p_xi_mlt_Bzero   ,
  dtres_Rnp_tau_np_p_stxi_mlt_Bzero ,
  dtres_Rnp_tau_np_n_mlt0_Bzero     ,
  dtres_Rnp_tau_np_n_mlt1_Bzero     ,
  dtres_Rnp_tau_np_n_xi_mlt_Bzero   ,
  dtres_Rnp_tau_np_n_stxi_mlt_Bzero ,
  dtres_Rnp_fd_np_sgl0_Bplus        ,
  dtres_Rnp_fd_np_sgl1_Bplus        ,
  dtres_Rnp_fp_np_sgl_Bplus         ,
  dtres_Rnp_tau_np_p_sgl0_Bplus     ,
  dtres_Rnp_tau_np_p_sgl1_Bplus     ,
  dtres_Rnp_tau_np_n_sgl0_Bplus     ,
  dtres_Rnp_tau_np_n_sgl1_Bplus     ,
  dtres_Rnp_fd_np_mlt0_Bplus        ,
  dtres_Rnp_fd_np_mlt1_Bplus        ,
  dtres_Rnp_fd_np_st_mlt_Bplus      ,
  dtres_Rnp_fd_np_xi_mlt_Bplus      ,
  dtres_Rnp_fd_np_stxi_mlt_Bplus    ,
  dtres_Rnp_fp_np_mlt_Bplus         ,
  dtres_Rnp_fn_np_mlt_Bplus         ,
  dtres_Rnp_tau_np_p_mlt0_Bplus     ,
  dtres_Rnp_tau_np_p_mlt1_Bplus     ,
  dtres_Rnp_tau_np_p_xi_mlt_Bplus   ,
  dtres_Rnp_tau_np_p_stxi_mlt_Bplus ,
  dtres_Rnp_tau_np_n_mlt0_Bplus     ,
  dtres_Rnp_tau_np_n_mlt1_Bplus     ,
  dtres_Rnp_tau_np_n_xi_mlt_Bplus   ,
  dtres_Rnp_tau_np_n_stxi_mlt_Bplus ,
  dtres_Rk_tau_k_1_p0               ,
  dtres_Rk_tau_k_1_p1               ,
  dtres_Rk_tau_k_2_p0               ,
  dtres_Rk_tau_k_2_p1               ,
  dtres_Rk_tau_k_1_n0               ,
  dtres_Rk_tau_k_1_n1               ,
  dtres_Rk_tau_k_2_n0               ,
  dtres_Rk_tau_k_2_n1               ,
  dtres_Rk_sigma_k0                 ,
  dtres_Rk_sigma_k1                 ,
  dtres_Rk_fk2                      ,
  dtres_Rol_sig_ol                  ,
  dtres_Rol_fol_sgl                 ,
  dtres_Rol_fol_mul                 ,
  dtres_Dt_cutoff
};

extern const dtres_syst_item_t dtres_items_fullrecon[];
extern const int num_dtres_items_fullrecon;

extern const dtres_syst_item_t dtres_items_fullrecon_b0[];
extern const int num_dtres_items_fullrecon_b0;

extern const dtres_syst_item_t dtres_items_fullrecon_bplus[];
extern const int num_dtres_items_fullrecon_bplus;

extern const dtres_syst_item_t dtres_items_partial[];
extern const int num_dtres_items_partial;

extern const dtres_syst_item_t dtres_items_partial_b0[];
extern const int num_dtres_items_partial_b0;

extern const dtres_syst_item_t dtres_items_partial_bplus[];
extern const int num_dtres_items_partial_bplus;

#ifdef __cplusplus
void dtres_systematics(const dtres_syst_item_t item, const double sigma,
		       const int expno,
		       dtres_param_t* param,
		       const dtres_param_t* const posi_err,
		       const dtres_param_t* const nega_err,
		       const double mc_param_factor = 2.0, const int oldtype=0);
#else /* __cplusplus */
void dtres_systematics(const dtres_syst_item_t item, const double sigma,
		       const int expno,
		       dtres_param_t* param,
		       const dtres_param_t* const posi_err,
		       const dtres_param_t* const nega_err,
		       const double mc_param_factor, const int oldtype);
#endif /* __cplusplus */

void dtres_systematics_expno(const dtres_syst_item_t item,
			     const double sigma,
			     const int expno, const int mc, const int version,
			     const double mc_param_factor, const int oldtype);
int cpdefpar_(int* expmc);

void calc_akck(const double cos_theta_b,
	       const double beta, const double Eb_cms, const double pb_cms,
	       double* ak, double* ck);

double Rk_fullrec(const double x,
		  const double dt, const double tau,
		  const double ak, const double ck);

double Rk_partial(const double x, const double dz,
		  const dtres_param_t * const param);

void Rk_partial_param(const double dz,
		      const dtres_param_t * const param,
		      double* tau_k_1_p, double* tau_k_2_p,
		      double* tau_k_1_n, double* tau_k_2_n,
		      double* sigma_k, double* fk2);

void Rnp_param(const int flavor, const int ntrk_asc, const int keeptagl,
	       const double Smain_asc, const double Stail_asc,
	       const dtres_param_t * const param,
	       double* fd, double* fp,
	       double* tau_np_p, double* tau_np_n,
	       double* tau_np_p_tl, double* tau_np_n_tl);

void Rnp_param_03(const int flavor, const int ntrk_asc, const int keeptagl,
		  const double xi_asc, const double st_asc,
		  const double Smain_asc, const double Stail_asc,
		  const dtres_param_t * const param,
		  double* fd, double* fp,
		  double* tau_np_p, double* tau_np_n,
		  double* tau_np_p_tl, double* tau_np_n_tl);

void Rnp_param_10(const int flavor, const int ntrk_asc, const int keeptagl,
		  const double xi_asc, const double st_asc,
		  const dtres_param_t * const param,
		  double* fd, double* fp,
		  double* tau_np_p, double* tau_np_n,
		  double* tau_np_p_tl, double* tau_np_n_tl);

void swap_rnp_param(double * fp,
		    double* tau_np_p, double* tau_np_n,
		    double* tau_np_p_tl, double* tau_np_n_tl);

double Rnp(const double x, 
	   const int flavor,
	   const int ntrk_asc, const double sz_asc,
	   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
	   const dtres_param_t * const param);

void Rasc_param(const int ntrk_asc,
		const double xi_asc, const double st_asc,
		const dtres_param_t * const param,
		double* ftail_asc,
		double* mu_main_asc, double* Smain_asc,
		double* mu_tail_asc, double* Stail_asc);

double Rasc(const double x, 
	    const int ntrk_asc, const double sz_asc,
	    const double chisq_z_asc, const int ndf_z_asc,
	    const dtres_param_t * const param);

double RascRnp(const double x,
	       const int flavor,
	       const int ntrk_asc, const double sz_asc,
	       const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
	       const dtres_param_t * const param);

void Rrec_param(const int ntrk_rec,
		const double xi_rec, const double st_rec,
		const dtres_param_t * const param,
		double* ftail_rec,
		double* mu_main_rec, double* Smain_rec,
		double* mu_tail_rec, double* Stail_rec);

double Rrec(const double x,
	    const int ntrk_rec, const double sz_rec,
	    const double chisq_z_rec, const int ndf_z_rec,
	    const dtres_param_t * const param);

double Rdet(const double x,
	    const int ntrk_rec, const double sz_rec,
	    const double chisq_z_rec, const int ndf_z_rec,
	    const int ntrk_asc, const double sz_asc,
	    const double chisq_z_asc, const int ndf_z_asc,
	    const dtres_param_t * const param);

double RdetRnp(const double x,
	       const int flavor,
	       const int ntrk_rec, const double sz_rec,
	       const double chisq_z_rec, const int ndf_z_rec,
	       const int ntrk_asc, const double sz_asc,
	       const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
	       const dtres_param_t * const param);

double EfRk_fullrec(const double x, const double tau,
		    const double ak, const double ck);

double AfRk_fullrec(const double x,
		    const double tau, const double dm,
		    const double ak, const double ck);

double MfRk_fullrec(const double x,
		    const double tau, const double dm,
		    const double ak, const double ck);

double EfRk_partial(const double x, const double tau,
		    const double dz, const dtres_param_t * const param);

double AfRk_partial(const double x, const double tau, const double dm,
		    const double dz, const dtres_param_t * const param);

double MfRk_partial(const double x, const double tau, const double dm,
		    const double dz, const dtres_param_t * const param);

double EfRkRdet_fullrec(const double x, const double tau,
			const double ak, const double ck,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param);

double AfRkRdet_fullrec(const double x,
			const double tau, const double dm,
			const double ak, const double ck,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param);

double MfRkRdet_fullrec(const double x,
			const double tau, const double dm,
			const double ak, const double ck,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param);

double EfRkRdet_partial(const double x, const double tau,
			const double dz,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param);

double AfRkRdet_partial(const double x, const double tau, const double dm,
			const double dz,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param);

double MfRkRdet_partial(const double x, const double tau, const double dm,
			const double dz,
			const int ntrk_rec, const double sz_rec,
			const double chisq_z_rec, const int ndf_z_rec,
			const int ntrk_asc, const double sz_asc,
			const double chisq_z_asc, const int ndf_z_asc,
			const dtres_param_t * const param);

double EfRkRdetRnp_fullrec(const double x, const int flavor,  const double tau,
			   const double ak, const double ck,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param);

double AfRkRdetRnp_fullrec(const double x, const int flavor, 
			   const double tau, const double dm,
			   const double ak, const double ck,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param);

double MfRkRdetRnp_fullrec(const double x, const int flavor, 
			   const double tau, const double dm,
			   const double ak, const double ck,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param);

double EfRkRdetRnp_partial(const double x, const int flavor,  const double tau,
			   const double dz,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param);

double AfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dm,
			   const double dz,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param);

double MfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dm,
			   const double dz,
			   const int ntrk_rec, const double sz_rec,
			   const double chisq_z_rec, const int ndf_z_rec,
			   const int ntrk_asc, const double sz_asc,
			   const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			   const dtres_param_t * const param);

double norm_Rk_partial(const double ll, const double ul, const double dz,
		       const dtres_param_t * const param);

double norm_Rasc(const double ll, const double ul,
		 const int ntrk_asc, const double sz_asc,
		 const double chisq_z_asc, const int ndf_z_asc,
		 const dtres_param_t * const param);

double norm_Rrec(const double ll, const double ul,
		 const int ntrk_rec, const double sz_rec,
		 const double chisq_z_rec, const int ndf_z_rec,
		 const dtres_param_t * const param);

double norm_Rnp(const double ll, const double ul,
		const int flavor,
		const int ntrk_asc, const double sz_asc,
		const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
		const dtres_param_t * const param);

double norm_RascRnp(const double ll, const double ul,
		    const int flavor,
		    const int ntrk_asc, const double sz_asc,
		    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
		    const dtres_param_t * const param);

double norm_Rdet(const double ll, const double ul,
		 const int ntrk_rec, const double sz_rec,
		 const double chisq_z_rec, const int ndf_z_rec,
		 const int ntrk_asc, const double sz_asc,
		 const double chisq_z_asc, const int ndf_z_asc,
		 const dtres_param_t * const param);

double norm_RdetRnp(const double ll, const double ul,
		    const int flavor,
		    const int ntrk_rec, const double sz_rec,
		    const double chisq_z_rec, const int ndf_z_rec,
		    const int ntrk_asc, const double sz_asc,
		    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
		    const dtres_param_t * const param);

double norm_EfRk_fullrec(const double ll, const double ul, const double tau,
			 const double ak, const double ck);

double norm_AfRk_fullrec(const double ll, const double ul,
			 const double tau, const double dm,
			 const double ak, const double ck);

double norm_MfRk_fullrec(const double ll, const double ul,
			 const double tau, const double dm,
			 const double ak, const double ck);

double norm_EfRk_partial(const double ll, const double ul, const double tau,
			 const double dz, const dtres_param_t * const param);

double norm_AfRk_partial(const double ll, const double ul, const double tau, const double dm,
			 const double dz, const dtres_param_t * const param);

double norm_MfRk_partial(const double ll, const double ul, const double tau, const double dm,
			 const double dz, const dtres_param_t * const param);

double norm_EfRkRdet_fullrec(const double ll, const double ul, const double tau,
			     const double ak, const double ck,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param);

double norm_AfRkRdet_fullrec(const double ll, const double ul,
			     const double tau, const double dm,
			     const double ak, const double ck,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param);

double norm_MfRkRdet_fullrec(const double ll, const double ul,
			     const double tau, const double dm,
			     const double ak, const double ck,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param);

double norm_EfRkRdet_partial(const double ll, const double ul, const double tau,
			     const double dz,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param);

double norm_AfRkRdet_partial(const double ll, const double ul, const double tau, const double dm,
			     const double dz,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param);

double norm_MfRkRdet_partial(const double ll, const double ul, const double tau, const double dm,
			     const double dz,
			     const int ntrk_rec, const double sz_rec,
			     const double chisq_z_rec, const int ndf_z_rec,
			     const int ntrk_asc, const double sz_asc,
			     const double chisq_z_asc, const int ndf_z_asc,
			     const dtres_param_t * const param);

double norm_EfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor,  const double tau,
				const double ak, const double ck,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param);

double norm_AfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor, 
				const double tau, const double dm,
				const double ak, const double ck,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param);

double norm_MfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor, 
				const double tau, const double dm,
				const double ak, const double ck,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param);

double norm_EfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param);

double norm_AfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau, const double dm,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param);

double norm_MfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau, const double dm,
				const double dz,
				const int ntrk_rec, const double sz_rec,
				const double chisq_z_rec, const int ndf_z_rec,
				const int ntrk_asc, const double sz_asc,
				const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				const dtres_param_t * const param);

double EfRdet(const double x, const double tau,
              const int ntrk_rec, const double sz_rec,
              const double chisq_z_rec, const int ndf_z_rec,
              const int ntrk_asc, const double sz_asc,
              const double chisq_z_asc, const int ndf_z_asc,
              const dtres_param_t * const param);

double AfRdet(const double x, const double tau, const double dm,
              const int ntrk_rec, const double sz_rec,
              const double chisq_z_rec, const int ndf_z_rec,
              const int ntrk_asc, const double sz_asc,
              const double chisq_z_asc, const int ndf_z_asc,
              const dtres_param_t * const param);

double MfRdet(const double x, const double tau, const double dm,
              const int ntrk_rec, const double sz_rec,
              const double chisq_z_rec, const int ndf_z_rec,
              const int ntrk_asc, const double sz_asc,
              const double chisq_z_asc, const int ndf_z_asc,
              const dtres_param_t * const param);

double norm_EfRdet(const double ll, const double ul, const double tau,
                   const int ntrk_rec, const double sz_rec,
                   const double chisq_z_rec, const int ndf_z_rec,
                   const int ntrk_asc, const double sz_asc,
                   const double chisq_z_asc, const int ndf_z_asc,
                   const dtres_param_t * const param);

double norm_AfRdet(const double ll, const double ul,
		   const double tau, const double dm,
                   const int ntrk_rec, const double sz_rec,
                   const double chisq_z_rec, const int ndf_z_rec,
                   const int ntrk_asc, const double sz_asc,
                   const double chisq_z_asc, const int ndf_z_asc,
                   const dtres_param_t * const param);

double norm_MfRdet(const double ll, const double ul,
		   const double tau, const double dm,
                   const int ntrk_rec, const double sz_rec,
                   const double chisq_z_rec, const int ndf_z_rec,
                   const int ntrk_asc, const double sz_asc,
                   const double chisq_z_asc, const int ndf_z_asc,
                   const dtres_param_t * const param);

double add_outlier(const int expno,
		   const double x,
		   const double Lin,
		   const int ntrk_rec, const int ntrk_asc,
		   const dtres_param_t * const param);

double add_norm_outlier(const int expno,
			const double ll, const double ul,
			const double nLin,
			const int ntrk_rec, const int ntrk_asc,
			const dtres_param_t * const param);

double Add_Outlier(const int expno,
		   const double x,
		   const double Lin,
		   const int ntrk_rec, const int ntrk_asc,
		   const dtres_param_t * const param,
		   const double nLi,
		   const double llmt,
		   const double ulmt,
		   const double alpha);

double Add_Outlier_with_Bkg(const int expno,
			    const double x, const double fsig,
			    const double Lsig, const double Lbkg,
			    const int ntrk_rec, const int ntrk_asc,
			    const dtres_param_t * const param,
			    const double nLsig, const double nLbkg,
			    const double llmt, const double ulmt,
			    const double alpha, const double beta);

void Calc_AkCk_Mass(const double cos_theta_b, const double Eb_cms,
		    double* ak, double* ck, const double mb, const double beta);

void Calc_AkCk_Flavor(const double cos_theta_b, const double Eb_cms,
		      double* ak, double* ck, const int b_flavor, const double beta);

double HAfRkRdetRnp_fullrec(const double x, const int flavor, 
			    const double tau, const double dg,
			    const double ak, const double ck,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param);

double HMfRkRdetRnp_fullrec(const double x, const int flavor, 
			    const double tau, const double dg,
			    const double ak, const double ck,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param);

double norm_HMfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor,
				 const double tau, const double dg,
				 const double ak, const double ck,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param);

double norm_HAfRkRdetRnp_fullrec(const double ll, const double ul, const int flavor,
				 const double tau, const double dg,
				 const double ak, const double ck,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param);

double HAfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dg,
			    const double dz,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param);

double HMfRkRdetRnp_partial(const double x, const int flavor,  const double tau, const double dg,
			    const double dz,
			    const int ntrk_rec, const double sz_rec,
			    const double chisq_z_rec, const int ndf_z_rec,
			    const int ntrk_asc, const double sz_asc,
			    const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
			    const dtres_param_t * const param);

double norm_HAfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau, const double dg,
				 const double dz,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param);

double norm_HMfRkRdetRnp_partial(const double ll, const double ul, const int flavor,  const double tau, const double dg,
				 const double dz,
				 const int ntrk_rec, const double sz_rec,
				 const double chisq_z_rec, const int ndf_z_rec,
				 const int ntrk_asc, const double sz_asc,
				 const double chisq_z_asc, const int ndf_z_asc, const int keeptagl,
				 const dtres_param_t * const param);

__CLOSE_LINKAGE__

#ifdef __cplusplus

inline void CalcAkCk(const double cos_theta_b, const double Eb_cms,
		     double* ak, double* ck, 
		     const double m_B,
		     const double beta = dt_resol_global::beta){
  Calc_AkCk_Mass(cos_theta_b, Eb_cms, ak, ck,m_B, beta);
  return;
}

inline void CalcAkCk(const double cos_theta_b, const double Eb_cms,
		     double* ak, double* ck,
		     const int b_flavor = 0, /* B0:=0, B+ := 1 */
		     const double beta = dt_resol_global::beta){
  Calc_AkCk_Flavor(cos_theta_b, Eb_cms, ak, ck, b_flavor, beta);
  return;
}

/* AddOutlier 

   returns (1-fol)*Li/nLi + fol*alpha*Lol/nLol;

   where nLol = \int^{ulmt}_{llmt} Lol(x) dx;

   optional input: nLi = \int^{ulmt}_{llmt} Lin(x) dx;

   !!! Be carefull not to divide Li by nLi two times !!!

*/

inline double AddOutlier(const int expno,
			 const double x,
			 const double Lin,
			 const int ntrk_rec, const int ntrk_asc,
			 const dtres_param_t * const param,
			 const double nLi = 1.0,
			 const double llmt = dt_resol_global::dt_llmt,
			 const double ulmt = dt_resol_global::dt_ulmt,
			 const double alpha = 1.0){
  return Add_Outlier(expno, x, Lin, ntrk_rec, ntrk_asc,
		     param, nLi, llmt, ulmt, alpha);
}

/* AddOutlier

   returns (1-fol)*(fsig*Lsig/nLsig+(1-fsig)*beta*Lbkg/nLbkg) + fol*alpha*Lol/nLol;

   where nLol = \int^{ulmt}_{llmt} Lol(x) dx;

   optional input: nLsig = \int^{ulmt}_{llmt} Lsig(x) dx,
                   nLbkg = \int^{ulmt}_{llmt} Lbkg(x) dx;

   !!! Be carefull not to divide Lsig by nLsig two times !!!
   !!! Be carefull not to divide Lbkg by nLbkg two times !!!

*/

inline double AddOutlierWithBkg(const int expno,
				const double x, const double fsig,
				const double Lsig, const double Lbkg,
				const int ntrk_rec, const int ntrk_asc,
				const dtres_param_t * const param,
				const double nLsig = 1.0, const double nLbkg = 1.0,
				const double llmt = dt_resol_global::dt_llmt,
				const double ulmt = dt_resol_global::dt_ulmt,
				const double alpha = 1.0, const double beta = 1.0){
  return Add_Outlier_with_Bkg(expno, x, fsig,Lsig, Lbkg,
			      ntrk_rec, ntrk_asc, param,
			      nLsig, nLbkg, llmt, ulmt, alpha, beta);
}

#endif /* __cplusplus */

#ifdef __cplusplus
/* Overload New Rnp_function */
inline void Rnp_param(const int flavor, const int ntrk_asc, const int keeptagl,
		      const double xi_asc, const double st_asc,
		      const double Smain_asc, const double Stail_asc,
		      const dtres_param_t * const param,
		      double* fd, double* fp,
		      double* tau_np_p, double* tau_np_n,
		      double* tau_np_p_tl, double* tau_np_n_tl){
  //if((param->fn_np_mlt)[0]>DTRES_EXTERAM_THRE){
  if((param->fn_np_mlt)[0] >= -0.00001){
    Rnp_param_10(flavor, ntrk_asc, keeptagl, xi_asc, st_asc,
		 param, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);
  }else{
    Rnp_param_03(flavor, ntrk_asc, keeptagl, xi_asc, st_asc, Smain_asc,  Stail_asc,
		 param, fd, fp, tau_np_p, tau_np_n, tau_np_p_tl, tau_np_n_tl);
  }
  
  return;
};

inline void dtres_systematics(const dtres_syst_item_t item,
			      const double sigma,
			      const int expno, const int mc = 0, const int version = 0,
			      const double mc_param_factor = 2.0, const int oldtype=0){
  dtres_systematics_expno(item, sigma,
			  expno, mc, version, mc_param_factor, oldtype);
  return;
};

#endif /* __cplusplus */


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __DT_RESOLUTION_H__ */

