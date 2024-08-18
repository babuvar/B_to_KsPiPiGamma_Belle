//               
//    dt_resolution_const.h   
//       ---- 
//    $Id: dt_resolution_const.h 11208 2011-03-01 14:53:58Z hitoshi $ 
//                
//    $Log$
//    Revision 1.5  2005/09/22 05:16:48  kohji
//    1. Fixed the out-of-range bug in test_tatami_toymc.cc
//           and norm_Mf() bug in libcnv.cc, pointed out by Sasa.
//    2. Added resolution parameters used for LP05,
//     dtres_param_svd2_lp05, dtres_param_posi_error_svd2_lp05,
//     	 dtres_param_nega_error_svd2_lp05.
//
//    Revision 1.4  2004/10/19 07:09:56  kohji
//    1) SVD1 and SVD2 resolution parameters are updated for 2004 publication.
//    2) "expno" is added to add_outlier and dtres_systematics functions to treat
//    difference between SVD1 and SVD2 properly[cpfit_ml:0771].
//    3) New functions are added for cosh/sinh terms
//     (norm_)HAfRkRdetRnp_{fullrec,partial} for HAf=exp(-|t|Gamma)sinh(DGamma/2 t)
//     (norm_)HMfRkRdetRnp_{fullrec,partial} for HMf=exp(-|t|Gamma)cosh(DGamma/2 t)
//
//    Revision 1.3  2003/10/01 01:24:23  katayama
//    compile with CC
//
//    Revision 1.2  2003/08/06 06:43:33  katayama
//    New version from Nakadaira san
//
//    Revision 1.1  2002/09/05 01:27:06  katayama
//    New package tatami from Nakadaira/Sumisawa san
//
//    Revision 1.2  2002/07/12 18:01:19  nakadair
//    Implement resolution parameter for 2002 ICHEP analysis
//    based on the lifetime fit using 80/fb data.
//
//    Revision 1.1  2002/07/12 09:20:18  nakadair
//    Iniaial import.
// 
//                

#ifndef __DT_RESOLUTION_CONST_H__
#define __DT_RESOLUTION_CONST_H__

#include "belle.h"
#include <cfloat>

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
#ifndef DTRES_EXTERAM_THRE
#define DTRES_EXTERAM_THRE -FLT_MAX
//#define DTRES_EXTERAM_THRE -999
#endif /* DTRES_EXTERAM_THRE */

#ifndef DTRES_EXTERAM_VAL
#define DTRES_EXTERAM_VAL (2.0*DTRES_EXTERAM_THRE)
#endif /* DTRES_EXTERAM_VAL */

#ifdef __cplusplus
#define __LINKAGE__  extern "C" {
#define __CLOSE_LINKAGE__       }
#else /* __cplusplus */
#define __LINKAGE__  
#define __CLOSE_LINKAGE__
#endif /* __cplusplus */

#ifdef __cplusplus

class dt_resol_global{
public:
  static double inv_bgc;    /* 1.0/(beta*gamma*c) */
  static double beta;       /* beta */ 
  static double mb;         /* mass_B */
  static double inv_mb;     /* 1.0/mass_B */
  static double mbzero;     /* (m_B0) */
  static double mbplus;     /* (m_B+) */
  static double dt_llmt;    /* lower limit of dt[ps]. (default: -70.0ps)*/ 
  static double dt_ulmt;    /* upper limit of dt[ps]. (default: +70.0ps)*/ 
  static double dt_lmt_var; /* Value for systmatic study (default: 30.0ps)*/ 
  static double rnp_kink_xi; /* for Rnp from 2010*/
  static double rnp_kink_st; /* for Rnp from 2010*/

public:
  static void set_bgc(const double bgc){
    inv_bgc = 1.0/bgc;
  }
  static void set_beta(const double x){
    beta = x;
  }
  static void set_bmass(const double x){
    inv_mb = 1.0/x;
  }
  static void set_b0_mass(const double x){
    mbzero     = x;
  }
  static void set_bplus_mass(const double x){
    mbplus     = x;
  }
public:
  enum {
    moriond2002   = 1,
    ichep2002_60  = 2,
    ichep2002_80  = 3,
    lp2003_150    = 4,
    svd1_ichep04  = 5,
    svd1_ichep06  = 6,
    svd1_2010mdlh = 7,
    svd1_2010nn   = 8,
    svd2_ichep04  = 10,
    svd2_lp05     = 11,
    svd2_ichep06  = 12,
    svd2_2010mdlh = 13,
    svd2_2010nn   = 14
  };
};

#else /* __cplusplus */

void dt_resol_global_set_bgc(const double bgc);
void dt_resol_global_set_bmass(const double x);
void dt_resol_global_set_b0_mass(const double x);
void dt_resol_global_set_bplus_mass(const double x);

extern const int dtres_version_moriond2002;
extern const int dtres_version_ichep2002_60;
extern const int dtres_version_ichep2002_80;
extern const int dtres_version_lp2003_150;
extern const int dtres_version_svd1_ichep04;
extern const int dtres_version_svd1_ichep06;
extern const int dtres_version_svd2_ichep04_v0;
extern const int dtres_version_svd2_ichep04;
extern const int dtres_version_svd1_2010mdlh;
extern const int dtres_version_svd1_2010nn1;
extern const int dtres_version_svd2_lp05;
extern const int dtres_version_svd2_ichep06;
extern const int dtres_version_svd2_2010mdlh;
extern const int dtres_version_svd2_2010nn;

#endif /* __cplusplus */

__LINKAGE__ 

typedef struct {
  /* R_det */
  /*  multiple track vertex */
  double Srec[2];
  double ftl_rec_mlt[2];
  double Stl_rec_mlt;
  double Sasc[2];
  double Snp, Snp_global;  /* for Rnp */
  double ftl_asc_mlt[2];
  double Stl_asc_mlt;
  /*  single track vertex */
  double Smn_rec;
  double Stl_rec;
  double ftl_rec;
  double Smn_asc;
  double Stl_asc;
  double ftl_asc;
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  double fd_np_sgl[2][2];
  double fp_np_sgl[2];
  double tau_np_p_sgl[2][2];
  double tau_np_n_sgl[2][2];
  /*  multiple track vertex  */
  double fd_np_mlt[2][2];
  double fd_np_st_mlt[2];
  double fd_np_xi_mlt[2];
  double fd_np_stxi_mlt[2];
  double fp_np_mlt[2];
  double fn_np_mlt[2];
  double tau_np_p_mlt[2][2];
  double tau_np_p_xi_mlt[2];
  double tau_np_p_stxi_mlt[2];
  double tau_np_n_mlt[2][2];
  double tau_np_n_xi_mlt[2];
  double tau_np_n_stxi_mlt[2];
  /* Rk parameters (For partial reconstracted event) */
  double tau_k_1_p[2];
  double tau_k_2_p[2];
  double tau_k_1_n[2];
  double tau_k_2_n[2];
  double sigma_k[2];
  double fk2;
  /* outlier parameters */
  double sig_ol;
  double fol_sgl;
  double fol_mul;
} dtres_param_t;


/*  dtres_param_t* get_dtres_param()                             */
/*   ... return the pointer to the default resolutionparameter   */
/*  argument:                                                    */
/*            expno   ... experimental number                    */
/*            mc      ... 0 for real data, = 1 for Monte Carlo   */
/*            version ... 0 for latest data set                  */

#ifdef __cplusplus
dtres_param_t* get_dtres_param(const int expno, const int mc = 0, const int version = 0);
dtres_param_t* get_dtres_param_posi_error(const int expno, const int version = 0);
dtres_param_t* get_dtres_param_nega_error(const int expno, const int version = 0);
#else  /* __cplusplus */
dtres_param_t* get_dtres_param(const int expno, const int mc, const int version);
dtres_param_t* get_dtres_param_posi_error(const int expno, const int version);
dtres_param_t* get_dtres_param_nega_error(const int expno, const int version);
#endif  /* __cplusplus */

/* Resolution Parameters for 2002 Moriond. */
extern dtres_param_t dtres_param_2002moriond;
extern dtres_param_t dtres_param_2002moriond_mc;
extern dtres_param_t dtres_param_posi_error_2002moriond;
extern dtres_param_t dtres_param_nega_error_2002moriond;

#if 1
/* Resolution Parameters for 2002 ICHEP.(from 60[fb^-1]) */
extern dtres_param_t dtres_param_2002ichep60;
extern dtres_param_t dtres_param_2002ichep60_mc;
extern dtres_param_t dtres_param_posi_error_2002ichep60;
extern dtres_param_t dtres_param_nega_error_2002ichep60;
#endif

/* Resolution Parameters for 2002 ICHEP.(from 80[fb^-1]) */
extern dtres_param_t dtres_param_2002ichep80;
extern dtres_param_t dtres_param_2002ichep80_mc;
extern dtres_param_t dtres_param_posi_error_2002ichep80;
extern dtres_param_t dtres_param_nega_error_2002ichep80;

/* Resolution Parameters for 2002 ICHEP.(from 80[fb^-1]) */
extern dtres_param_t dtres_param_2003lp03;
extern dtres_param_t dtres_param_2003lp03_mc;
extern dtres_param_t dtres_param_posi_error_2003lp03;
extern dtres_param_t dtres_param_nega_error_2003lp03;

/* SVD1 Resolution Parameters for ICHEP2004.(from 150[fb^-1]) (LP03 debugged) */
extern dtres_param_t dtres_param_svd1_ichep04;
extern dtres_param_t dtres_param_svd1_ichep04_mc;
extern dtres_param_t dtres_param_posi_error_svd1_ichep04;
extern dtres_param_t dtres_param_nega_error_svd1_ichep04;

/* SVD2 Resolution Parameters for ICHEP2004.(from 120[fb^-1]) (D*lnu duplicated events bug) */
extern dtres_param_t dtres_param_svd2_ichep04_v0;
extern dtres_param_t dtres_param_posi_error_svd2_ichep04_v0;
extern dtres_param_t dtres_param_nega_error_svd2_ichep04_v0;

/* SVD2 Resolution Parameters for 2004 publication.(from 120[fb^-1]) */
extern dtres_param_t dtres_param_svd2_ichep04;
extern dtres_param_t dtres_param_svd2_ichep04_mc;
extern dtres_param_t dtres_param_posi_error_svd2_ichep04;
extern dtres_param_t dtres_param_nega_error_svd2_ichep04;

extern dtres_param_t dtres_param_svd2_lp05;
extern dtres_param_t dtres_param_posi_error_svd2_lp05;
extern dtres_param_t dtres_param_nega_error_svd2_lp05;

extern dtres_param_t dtres_param_center_svd1_ichep06_mc;
extern dtres_param_t dtres_param_center_svd1_ichep06_data;
extern dtres_param_t dtres_param_posi_error_svd1_ichep06_data;
extern dtres_param_t dtres_param_nega_error_svd1_ichep06_data;

extern dtres_param_t dtres_param_center_svd2_ichep06_mc;
extern dtres_param_t dtres_param_center_svd2_ichep06_data;
extern dtres_param_t dtres_param_posi_error_svd2_ichep06_data;
extern dtres_param_t dtres_param_nega_error_svd2_ichep06_data;

extern dtres_param_t dtres_param_center_svd1_2010mdlh_mc;
extern dtres_param_t dtres_param_center_svd1_2010mdlh_data;
extern dtres_param_t dtres_param_posi_error_svd1_2010mdlh_data;
extern dtres_param_t dtres_param_nega_error_svd1_2010mdlh_data;

extern dtres_param_t dtres_param_center_svd2_2010mdlh_mc;
extern dtres_param_t dtres_param_center_svd2_2010mdlh_data;
extern dtres_param_t dtres_param_posi_error_svd2_2010mdlh_data;
extern dtres_param_t dtres_param_nega_error_svd2_2010mdlh_data;

extern dtres_param_t dtres_param_center_svd1_2010nn_mc;
extern dtres_param_t dtres_param_center_svd1_2010nn_data;
extern dtres_param_t dtres_param_center_svd2_2010nn_mc;
extern dtres_param_t dtres_param_center_svd2_2010nn_data;

/* Default Resolution Parameters */
/* (== parameter set for SVD1 ICHEP04 150[fb^-1] (LP03 debugged) */
/* (== parameter set for SVD2 2004 publication 120[fb^-1] ) */

extern dtres_param_t dtres_param_default;
extern dtres_param_t dtres_param_default_mc;
extern dtres_param_t dtres_param_posi_error;
extern dtres_param_t dtres_param_nega_error;

extern dtres_param_t dtres_param_default_svd2;
extern dtres_param_t dtres_param_default_svd2_mc;
extern dtres_param_t dtres_param_posi_error_svd2;
extern dtres_param_t dtres_param_nega_error_svd2;

__CLOSE_LINKAGE__

#endif /* __DT_RESOLUTION_CONST_H__ */

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
