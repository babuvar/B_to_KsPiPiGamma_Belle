//
//   dt_resolution_const.cc
//     --- Description: Definitions of conststant for Resolution Function Library.
//   $Id: dt_resolution_const.cc 11330 2011-07-25 11:29:28Z hitoshi $
//
//   $Log$
//   Revision 1.4  2005/09/22 05:16:47  kohji
//   1. Fixed the out-of-range bug in test_tatami_toymc.cc
//          and norm_Mf() bug in libcnv.cc, pointed out by Sasa.
//   2. Added resolution parameters used for LP05,
//    dtres_param_svd2_lp05, dtres_param_posi_error_svd2_lp05,
//    	 dtres_param_nega_error_svd2_lp05.
//
//   Revision 1.3  2004/10/19 07:09:55  kohji
//   1) SVD1 and SVD2 resolution parameters are updated for 2004 publication.
//   2) "expno" is added to add_outlier and dtres_systematics functions to treat
//   difference between SVD1 and SVD2 properly[cpfit_ml:0771].
//   3) New functions are added for cosh/sinh terms
//    (norm_)HAfRkRdetRnp_{fullrec,partial} for HAf=exp(-|t|Gamma)sinh(DGamma/2 t)
//    (norm_)HMfRkRdetRnp_{fullrec,partial} for HMf=exp(-|t|Gamma)cosh(DGamma/2 t)
//
//   Revision 1.2  2003/08/06 06:43:32  katayama
//   New version from Nakadaira san
//
//   Revision 1.1  2002/09/05 01:27:04  katayama
//   New package tatami from Nakadaira/Sumisawa san
//
//   Revision 1.2  2002/07/12 18:00:25  nakadair
//   Implement resolution parameter for 2002 ICHEP analysis
//   based on the lifetime fit using 80/fb data.
//
//   Revision 1.1  2002/07/12 09:19:15  nakadair
//   Initial import.
//
//

#include "belle.h"
#include "tatami/dt_resolution_const.h"



#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
/* 1.0/(beta*gamma*c) */
double dt_resol_global::inv_bgc = 78.48566945838871754705; /* (0.425*0.03)^{-1} */
double dt_resol_global::beta    = 0.39114064034485169394; /* 0.425/sqrt(1.0+0.425^2) */ 
double dt_resol_global::inv_mb  = 0.18941546387847103837; /* (5.2794)^{-1} */;
double dt_resol_global::mbzero      = 5.2794;     /* (m_B0) */;
double dt_resol_global::mbplus      = 5.2790;     /* (m_B+) */;
double dt_resol_global::dt_llmt     = -70.0;  /* lower limit of dt[ps].*/ 
double dt_resol_global::dt_ulmt     =  70.0;  /* upper limit of dt[ps].*/ 
double dt_resol_global::dt_lmt_var  =  30.0;  /* for syetmatic study */
double dt_resol_global::rnp_kink_xi = 3.5; /* for Rnp from 2010*/
double dt_resol_global::rnp_kink_st = 0.75; /* for Rnp from 2010*/
  //double dt_resol_global::rnp_kink_xi = 5; /* for Rnp from 2010*/
  //double dt_resol_global::rnp_kink_st = 0.8; /* for Rnp from 2010*/

const int dtres_version_moriond2002   = (int)dt_resol_global::moriond2002;
const int dtres_version_ichep2002_60  = (int)dt_resol_global::ichep2002_60;
const int dtres_version_ichep2002_80  = (int)dt_resol_global::ichep2002_80;
const int dtres_version_lp2003_150    = (int)dt_resol_global::lp2003_150;
const int dtres_version_svd1_ichep04  = (int)dt_resol_global::svd1_ichep04;
const int dtres_version_svd1_ichep06  = (int)dt_resol_global::svd1_ichep06;
const int dtres_version_svd1_2010mdlh = (int)dt_resol_global::svd1_2010mdlh;
const int dtres_version_svd2_ichep04  = (int)dt_resol_global::svd2_ichep04;
const int dtres_version_svd2_lp05     = (int)dt_resol_global::svd2_lp05;
const int dtres_version_svd2_ichep06  = (int)dt_resol_global::svd2_ichep06;
const int dtres_version_svd2_2010mdlh = (int)dt_resol_global::svd2_2010mdlh;
const int dtres_version_svd2_2010nn   = (int)dt_resol_global::svd2_2010nn;

__LINKAGE__

dtres_param_t* get_dtres_param(const int expno, const int mc,
			       const int version){
  if(expno<30){
    switch(version){
     case dt_resol_global::moriond2002:
      return (mc==0) ? &dtres_param_2002moriond : &dtres_param_2002moriond_mc;
     case dt_resol_global::ichep2002_60:
      return (mc==0) ? &dtres_param_2002ichep60 : &dtres_param_2002ichep60_mc;
     case dt_resol_global::ichep2002_80:
      return (mc==0) ? &dtres_param_2002ichep80 : &dtres_param_2002ichep80_mc;
     case dt_resol_global::lp2003_150:
      return (mc==0) ? &dtres_param_2003lp03    : &dtres_param_2003lp03_mc;
     case dt_resol_global::svd1_ichep04:
      return (mc==0) ? &dtres_param_svd1_ichep04: &dtres_param_svd1_ichep04_mc;
     case dt_resol_global::svd1_ichep06:
       return (mc==0) ? &dtres_param_center_svd1_ichep06_data: &dtres_param_center_svd1_ichep06_mc;
     case dt_resol_global::svd1_2010mdlh:
     default:
       return (mc==0) ? &dtres_param_center_svd1_2010mdlh_data: &dtres_param_center_svd1_2010mdlh_mc;
     case dt_resol_global::svd1_2010nn:
       return (mc==0) ? &dtres_param_center_svd1_2010nn_data: &dtres_param_center_svd1_2010nn_mc;
    }
  }
  else {
    if(version==dt_resol_global::svd2_ichep04){
      return (mc==0) ? &dtres_param_svd2_ichep04  : &dtres_param_svd2_ichep04_mc;
    }else if(version==dt_resol_global::svd2_lp05){
      return (mc==0) ? &dtres_param_svd2_lp05  : &dtres_param_svd2_ichep04_mc;
    }else if(version==dt_resol_global::svd2_ichep06){
      return (mc==0) ? &dtres_param_center_svd2_ichep06_data  : &dtres_param_center_svd2_ichep06_mc;
    }else if(version==dt_resol_global::svd2_2010nn){
      return (mc==0) ? &dtres_param_center_svd2_2010nn_data  : &dtres_param_center_svd2_2010nn_mc;
    }else{
      return (mc==0) ? &dtres_param_center_svd2_2010mdlh_data  : &dtres_param_center_svd2_2010mdlh_mc;
    }
  }
  /* return (mc==0) ? &dtres_param_default : &dtres_param_default_mc; */
}

dtres_param_t* get_dtres_param_posi_error(const int expno, const int version){
  if(expno<30){
    switch(version){
     case dt_resol_global::moriond2002:
      return &dtres_param_posi_error_2002moriond;
     case dt_resol_global::ichep2002_60:
      return &dtres_param_posi_error_2002ichep60;
     case dt_resol_global::ichep2002_80:
      return &dtres_param_posi_error_2002ichep80;
     case dt_resol_global::lp2003_150:
      return &dtres_param_posi_error_2003lp03;
     case dt_resol_global::svd1_ichep04:
      return &dtres_param_posi_error_svd1_ichep04;
     case dt_resol_global::svd1_2010mdlh:
     default:
      return &dtres_param_posi_error_svd1_2010mdlh_data;
    }
  }
  else {
    if(version==dt_resol_global::svd2_ichep04){
      return &dtres_param_posi_error_svd2_ichep04;
    }else if(version==dt_resol_global::svd2_lp05){
      return &dtres_param_posi_error_svd2_lp05;
    }else if(version==dt_resol_global::svd2_ichep06){
      return &dtres_param_posi_error_svd2_ichep06_data;
    } else {
      return &dtres_param_posi_error_svd2_2010mdlh_data;
    }
  }
  /* return &dtres_param_posi_error; */
}

dtres_param_t* get_dtres_param_nega_error(const int expno, const int version){
  if(expno<30){
    switch(version){
     case dt_resol_global::moriond2002:
      return &dtres_param_nega_error_2002moriond;
     case dt_resol_global::ichep2002_60:
      return &dtres_param_nega_error_2002ichep60;
     case dt_resol_global::ichep2002_80:
      return &dtres_param_nega_error_2002ichep80;
     case dt_resol_global::lp2003_150:
      return &dtres_param_nega_error_2003lp03;
     case dt_resol_global::svd1_ichep04:
      return &dtres_param_nega_error_svd1_ichep04;
     case dt_resol_global::svd1_2010mdlh:
     default:
      return &dtres_param_nega_error_svd1_2010mdlh_data;
    }
  }
  else {
    if(version==dt_resol_global::svd2_ichep04){
      return &dtres_param_nega_error_svd2_ichep04;
    }else if(version==dt_resol_global::svd2_lp05){
      return &dtres_param_nega_error_svd2_lp05;
    }else if(version==dt_resol_global::svd2_ichep06){
      return &dtres_param_nega_error_svd2_ichep06_data;
    } else {
      return &dtres_param_nega_error_svd2_2010mdlh_data;
    }
  }
  /* return &dtres_nega_nega_error; */
}


dtres_param_t dtres_param_2002moriond = {
  {0.80867, 0.15440}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.75304, 0.064044},/* Sasc[2] */
  DTRES_EXTERAM_VAL, 1.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.64685, /* Smn_rec */
  3.00, /* Stl_rec */
  0.082707, /* ftl_rec */
  0.64685, /* Smn_asc */
  3.00, /* Stl_asc */
  0.082707, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.78658, 0.78658}, {0.76303, 0.76303}}, /* fd_np_sgl[2][2] */
  {0.78965, 0.75683}, /* fp_np_sgl[2] */
  {{ 0.10835, 1.3212}, {-0.0018994, 1.1132}}, /* tau_np_p_sgl[2][2] */
  {{-0.28124, 1.5828}, {-0.37471,   1.5477}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.67617, 0.67617}, {0.65020, 0.65020}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.95471, 0.96322}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{-0.0098103, 0.92739}, { 0.036614, 0.67378}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{-0.19398,   1.9896},  {-0.26871,  2.0698}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.17875,  0.015811},/* tau_k_1_p[2] */
  {0.066674, 0.031771},/* tau_k_2_p[2] */
  {0.17353,  0.015699},/* tau_k_1_n[2] */
  {0.072707, 0.028984},/* tau_k_2_n[2] */
  {-0.026323, 0.063469},/* sigma_k[2] */
  0.46008,/* fk2 */
  /* outlier parameters */
  36.185,/* sig_ol */
  0.030574,/* fol_sgl */
  5.83107e-4 /* fol_mul */
};

dtres_param_t dtres_param_2002moriond_mc = {
  {0.99029, 0.053705}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.97238, 0.028433},/* Sasc[2] */
  DTRES_EXTERAM_VAL, 1.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  1.0522, /* Smn_rec */
  4.43616, /* Stl_rec */
  0.037956, /* ftl_rec */
  1.0522, /* Smn_asc */
  4.43616, /* Stl_asc */
  0.037956, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.78658, 0.78658}, {0.76303, 0.76303}}, /* fd_np_sgl[2][2] */
  {0.78965, 0.75683}, /* fp_np_sgl[2] */
  {{ 0.10835, 1.3212}, {-0.0018994, 1.1132}}, /* tau_np_p_sgl[2][2] */
  {{-0.28124, 1.5828}, {-0.37471,   1.5477}}, /* tau_np_n_sgl[2][2] */
    /*  multiple track vertex  */
  {{0.67617, 0.67617}, {0.65020, 0.65020}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.95471, 0.96322}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{-0.0098103, 0.92739}, { 0.036614, 0.67378}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{-0.19398,   1.9896},  {-0.26871,  2.0698}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.17875,  0.015811},/* tau_k_1_p[2] */
  {0.066674, 0.031771},/* tau_k_2_p[2] */
  {0.17353,  0.015699},/* tau_k_1_n[2] */
  {0.072707, 0.028984},/* tau_k_2_n[2] */
  {-0.026323, 0.063469},/* sigma_k[2] */
  0.46008,/* fk2 */
  /* outlier parameters */
  43.625,/* sig_ol */
  0.0021983,/* fol_sgl */
  9.8589e-5 /* fol_mul */
};

dtres_param_t dtres_param_posi_error_2002moriond = {
  {0.14626, 0.12578E-01}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.64024E-01, 0.49585E-02},/* Sasc[2] */
  0.0, 0.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.74397E-01, /* Smn_rec */
  2.2309, /* Stl_rec */
  0.82990E-01, /* ftl_rec */
  0.74397E-01, /* Smn_asc */
  2.2309, /* Stl_asc */
  0.82990E-01, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.10403E-01, 0.10403E-01}, {0.16956E-01, 0.16956E-01}}, /* fd_np_sgl[2][2] */
  {0.20453E-01, 0.25419E-01}, /* fp_np_sgl[2] */
  {{ 0.67649E-01, 0.99228E-01}, {0.65600E-01, 0.98701E-01}}, /* tau_np_p_sgl[2][2] */
  {{0.12974, 0.21296}, {0.11079,   0.20675}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.68679E-02, 0.68679E-02}, {0.10054E-01, 0.10054E-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {0.35046E-02, 0.38156E-02}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{0.11298E-01, 0.24578E-01}, { 0.11952E-01, 0.25257E-01}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{0.77661E-01,   0.18201},  {0.98542E-01,  0.23474}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.72523E-02, 0.22225E-02},/* tau_k_1_p[2] */
  {0.51597E-02, 0.29020E-02},/* tau_k_2_p[2] */
  {0.64781E-02, 0.23897E-02},/* tau_k_1_n[2] */
  {0.51559E-02, 0.36122E-02},/* tau_k_2_n[2] */
  {0.31491E-02, 0.17754E-02},/* sigma_k[2] */
  0.50539E-01,/* fk2 */
  /* outlier parameters */
  4.9605,/* sig_ol */
  0.35811E-02,/* fol_sgl */
  0.30159E-03 /* fol_mul */
};

dtres_param_t dtres_param_nega_error_2002moriond = {
  {-0.14999, -0.12750E-01}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {-0.64692E-01, -0.48536E-02},/* Sasc[2] */
  0.0, 0.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  -0.83469E-01, /* Smn_rec */
  -0.98867, /* Stl_rec */
  -0.45125E-01, /* ftl_rec */
  -0.83469E-01, /* Smn_asc */
  -0.98867, /* Stl_asc */
  -0.45125E-01, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{-0.10980E-01, -0.10980E-01}, {-0.18225E-01, -0.18225E-01}}, /* fd_np_sgl[2][2] */
  {-0.21369E-01, -0.26411E-01}, /* fp_np_sgl[2] */
  {{-0.66530E-01, -0.93746E-01}, {-0.64910E-01, -0.92146E-01}}, /* tau_np_p_sgl[2][2] */
  {{-0.14684, -0.18442}, {-0.12184,   -0.18242}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-0.69734E-02, -0.69734E-02}, {-0.10266E-01, -0.10266E-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-0.37070E-02, -0.41127E-02}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-0.11271E-01, -0.24147E-01}, { -0.11900E-01, -0.24716E-01}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-0.77206E-01,   -0.16940},  {-0.99131E-01,  -0.21300}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {-0.64564E-02, -0.24297E-02},/* tau_k_1_p[2] */
  {-0.52206E-02, -0.29522E-02},/* tau_k_2_p[2] */
  {-0.57495E-02, -0.27099E-02},/* tau_k_1_n[2] */
  {-0.53495E-02, -0.35106E-02},/* tau_k_2_n[2] */
  {-0.33081E-02, -0.17945E-02},/* sigma_k[2] */
   -0.48524E-01,/* fk2 */
  /* outlier parameters */
  -3.5189,/* sig_ol */
  -0.35889E-02,/* fol_sgl */
  -0.22493E-03 /* fol_mul */
};

dtres_param_t dtres_param_2002ichep60 = {
  {0.97923, 0.96833E-01}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.72911, 0.46476E-01},/* Sasc[2] */
  DTRES_EXTERAM_VAL, 1.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.90311, /* Smn_rec */
  1.8587, /* Stl_rec */
  0.13233, /* ftl_rec */
  0.90311, /* Smn_asc */
  1.8587, /* Stl_asc */
  0.13233, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.74727, 0.74727}, {0.73961, 0.73961}}, /* fd_np_sgl[2][2] */
  {0.86839, 0.81504}, /* fp_np_sgl[2] */
  {{0.14381, 1.5446},  {0.051255, 0.97426}}, /* tau_np_p_sgl[2][2] */
  {{0.032259, 1.8943}, {-0.21700, 1.2160}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.48961, 0.48961}, {0.41154, 0.41154}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.95860, 0.97100}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{-0.048612, 1.0079}, { 0.015022, 0.73559}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{-0.13793,  2.0569}, {-0.069551, 1.9270}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.17875,  0.015811},/* tau_k_1_p[2] */
  {0.066674, 0.031771},/* tau_k_2_p[2] */
  {0.17353,  0.015699},/* tau_k_1_n[2] */
  {0.072707, 0.028984},/* tau_k_2_n[2] */
  {-0.026323, 0.063469},/* sigma_k[2] */
  0.46008,/* fk2 */
  /* outlier parameters */
  40.886,/* sig_ol */
  0.25712E-01,/* fol_sgl */
  0.17456E-03 /* fol_mul */
};

dtres_param_t dtres_param_2002ichep60_mc = {
  {1.0133, 0.75507E-01}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.83553, 0.38561E-01},/* Sasc[2] */
  DTRES_EXTERAM_VAL, 1.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.95974, /* Smn_rec */
  3.5308,  /* Stl_rec */
  0.48866E-01, /* ftl_rec */
  0.95974, /* Smn_asc */
  3.5308, /* Stl_asc */
  0.48866E-01, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.81759, 0.81759}, {0.70884, 0.70884}}, /* fd_np_sgl[2][2] */
  {0.86839, 0.81504}, /* fp_np_sgl[2] */
  {{0.14381, 1.5446},  {0.051255, 0.97426}}, /* tau_np_p_sgl[2][2] */
  {{0.032259, 1.8943}, {-0.21700, 1.2160}},  /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.57727, 0.57727}, {0.55096, 0.55096}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.95860, 0.97100}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{-0.048612, 1.0079}, { 0.015022, 0.73559}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{-0.13793,  2.0569}, {-0.069551, 1.9270}},  /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.17875,  0.015811},  /* tau_k_1_p[2] */
  {0.066674, 0.031771},  /* tau_k_2_p[2] */
  {0.17353,  0.015699},  /* tau_k_1_n[2] */
  {0.072707, 0.028984},  /* tau_k_2_n[2] */
  {-0.026323, 0.063469}, /* sigma_k[2] */
  0.46008,/* fk2 */
  /* outlier parameters */
  39.823,/* sig_ol */
  0.17230E-01,/* fol_sgl */
  0.49039E-10 /* fol_mul */
};

dtres_param_t dtres_param_posi_error_2002ichep60 = {
  {0.12787,     0.90112E-02}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.53337E-01, 0.27372E-02},/* Sasc[2] */
  0.0, 0.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.10080, /* Smn_rec */
  1.1067, /* Stl_rec */
  0.24478, /* ftl_rec */
  0.10080, /* Smn_asc */
  1.1067, /* Stl_asc */
  0.24478, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.40158E-01, 0.40158E-01}, {0.62633E-01, 0.62633E-01}}, /* fd_np_sgl[2][2] */
  {0.24546E-01, 0.36790E-01}, /* fp_np_sgl[2] */
  {{0.88796E-01, 0.11905}, {0.56364E-01, 0.87761E-01}}, /* tau_np_p_sgl[2][2] */
  {{0.28369, 0.66028},     {0.11292,     0.36751}},    /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.50109E-01, 0.50109E-01}, {0.51791E-01, 0.51791E-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {0.29303E-02, 0.25720E-02}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{0.98059E-02, 0.22704E-01}, {0.88301E-02, 0.19921E-01}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{0.83093E-01, 0.18771},     {0.99094E-01,  0.23128}},    /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.72523E-02, 0.22225E-02}, /* tau_k_1_p[2] */
  {0.51597E-02, 0.29020E-02}, /* tau_k_2_p[2] */
  {0.64781E-02, 0.23897E-02}, /* tau_k_1_n[2] */
  {0.51559E-02, 0.36122E-02}, /* tau_k_2_n[2] */
  {0.31491E-02, 0.17754E-02}, /* sigma_k[2] */
  0.50539E-01,/* fk2 */
  /* outlier parameters */
  5.1933,/* sig_ol */
  0.22170E-02,/* fol_sgl */
  0.12741E-03 /* fol_mul */
};

dtres_param_t dtres_param_nega_error_2002ichep60 = {
  {-0.13695,     -0.87390E-02},/* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {-0.52915E-01, -0.26821E-02},/* Sasc[2] */
  0.0, 0.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  -0.67718, /* Smn_rec */
  -0.78133, /* Stl_rec */
  -0.11069, /* ftl_rec */
  -0.67718, /* Smn_asc */
  -0.78133, /* Stl_asc */
  -0.11069, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{-0.41958E-01, -0.41958E-01}, {-0.65599E-01, -0.65599E-01}}, /* fd_np_sgl[2][2] */
  {-0.27165E-01, -0.39926E-01}, /* fp_np_sgl[2] */
  {{-0.84942E-01, -0.11450}, {-0.53637E-01, -0.81113E-01}}, /* tau_np_p_sgl[2][2] */
  {{-0.31388,     -0.48587}, {-0.13901,     -0.27466}},     /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-0.52270E-01, -0.52270E-01}, {-0.53664E-01, -0.53664E-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-0.30836E-02, -0.27445E-02}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-0.97936E-02, -0.22315E-01}, {-0.87970E-02, -0.19591E-01}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-0.82486E-01, -0.17592},     {-0.99536E-01, -0.21200}},     /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {-0.64564E-02, -0.24297E-02},/* tau_k_1_p[2] */
  {-0.52206E-02, -0.29522E-02},/* tau_k_2_p[2] */
  {-0.57495E-02, -0.27099E-02},/* tau_k_1_n[2] */
  {-0.53495E-02, -0.35106E-02},/* tau_k_2_n[2] */
  {-0.33081E-02, -0.17945E-02},/* sigma_k[2] */
   -0.48524E-01,/* fk2 */
  /* outlier parameters */
  -3.8473,/* sig_ol */
  -0.20967E-02,/* fol_sgl */
  -0.89694E-04 /* fol_mul */
};

dtres_param_t dtres_param_2002ichep80 = {
  /* R_det */
  /*  multiple track vertex */
  {0.98748, 0.94156E-01}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.77752, 0.43639E-01},/* Sasc[2] */
  DTRES_EXTERAM_VAL, 1.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.97202, /* Smn_rec */
  0.97202, /* Stl_rec */
  0.0000, /* ftl_rec */
  0.97202, /* Smn_asc */
  0.97202, /* Stl_asc */
  0.0000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.70145, 0.70145}, {0.76352, 0.76352}}, /* fd_np_sgl[2][2] */
  {0.78483, 0.80424}, /* fp_np_sgl[2] */
  {{0.37369, 1.4719}, {0.21383, 1.2342}}, /* tau_np_p_sgl[2][2] */
  {{0.19275, 1.8114}, {-0.14846, 2.0861}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.55463, 0.55463}, {0.43972, 0.43972}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.95860, 0.97100}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{-0.48612E-01, 1.0079}, {0.15022E-01, 0.73559}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{-0.13793, 2.0569}, {-0.69551E-01, 1.9270}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.17875, 0.015811}, /* tau_k_1_p[2] */
  {0.066674, 0.031771}, /* tau_k_2_p[2] */
  {0.17353, 0.015699}, /* tau_k_1_n[2] */
  {0.072707, 0.028984}, /* tau_k_2_n[2] */
  {-0.026323, 0.063469}, /* sigma_k[2] */
  0.46008, /* fk2 */
  /* outlier parameters */
  41.990, /* sig_ol */
  0.26866E-01, /* fol_sgl */
  0.16477E-03, /* fol_mul */
};

dtres_param_t dtres_param_2002ichep80_mc = {
  /* R_det */
  /*  multiple track vertex */
  {0.98522, 0.75790E-01}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.84204, 0.38563E-01},/* Sasc[2] */
  DTRES_EXTERAM_VAL, 1.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.99845, /* Smn_rec */
  0.99845, /* Stl_rec */
  0.00000, /* ftl_rec */
  0.99845, /* Smn_asc */
  0.99845, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.77570, 0.77570}, {0.76243, 0.76243}}, /* fd_np_sgl[2][2] */
  {0.78483, 0.80424}, /* fp_np_sgl[2] */
  {{0.37369, 1.4719}, {0.21383, 1.2342}}, /* tau_np_p_sgl[2][2] */
  {{0.19275, 1.8114}, {-0.14846, 2.0861}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.57931, 0.57931}, {0.55346, 0.55346}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.95860, 0.97100}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{-0.48612E-01, 1.0079}, {0.15022E-01, 0.73559}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{-0.13793, 2.0569}, {-0.69551E-01, 1.9270}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.17875, 0.015811}, /* tau_k_1_p[2] */
  {0.066674, 0.031771}, /* tau_k_2_p[2] */
  {0.17353, 0.015699}, /* tau_k_1_n[2] */
  {0.072707, 0.028984}, /* tau_k_2_n[2] */
  {-0.026323, 0.063469}, /* sigma_k[2] */
  0.46008, /* fk2 */
  /* outlier parameters */
  37.705, /* sig_ol */
  0.18623E-01, /* fol_sgl */
  0.43408E-10, /* fol_mul */
};

dtres_param_t dtres_param_posi_error_2002ichep80 = {
  /* R_det */
  /*  multiple track vertex */
  {0.11676, 0.77193E-02}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.48048E-01, 0.23884E-02},/* Sasc[2] */
  0.0, 0.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.45240E-01, /* Smn_rec */
  0.00000, /* Stl_rec */
  0.00000, /* ftl_rec */
  0.45240E-01, /* Smn_asc */
  0.00000, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.39384E-01, 0.39384E-01}, {0.42221E-01, 0.42221E-01}}, /* fd_np_sgl[2][2] */
  {0.12991E-01, 0.16112E-01}, /* fp_np_sgl[2] */
  {{0.83706E-01, 0.10534}, {0.70416E-01, 0.10438}}, /* tau_np_p_sgl[2][2] */
  {{0.14627, 0.22107}, {0.19333, 0.31817}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.41498E-01, 0.41498E-01}, {0.44734E-01, 0.44734E-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {0.29303E-02, 0.25720E-02}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{0.98059E-02, 0.22704E-01}, {0.88301E-02, 0.19921E-01}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{0.83093E-01, 0.18771}, {0.99094E-01, 0.23128}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.72523E-02, 0.22225E-02}, /* tau_k_1_p[2] */
  {0.51597E-02, 0.29020E-02}, /* tau_k_2_p[2] */
  {0.64781E-02, 0.23897E-02}, /* tau_k_1_n[2] */
  {0.51559E-02, 0.36122E-02}, /* tau_k_2_n[2] */
  {0.31491E-02, 0.17754E-02}, /* sigma_k[2] */
  0.50539E-01, /* fk2 */
  /* outlier parameters */
  4.6291, /* sig_ol */
  0.19126E-02, /* fol_sgl */
  0.11321E-03, /* fol_mul */
};

dtres_param_t dtres_param_nega_error_2002ichep80 = {
  /* R_det */
  /*  multiple track vertex */
  {-0.12401, -0.75469E-02}, /* Srec[2] */ 
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {-0.47923E-01, -0.23480E-02},/* Sasc[2] */
  0.0, 0.0, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  -0.45434E-01, /* Smn_rec */
  0.00000, /* Stl_rec */
  0.00000, /* ftl_rec */
  -0.45434E-01, /* Smn_asc */
  0.00000, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{-0.41675E-01, -0.41675E-01}, {-0.45037E-01, -0.45037E-01}}, /* fd_np_sgl[2][2] */
  {-0.13423E-01, -0.16872E-01}, /* fp_np_sgl[2] */
  {{-0.82084E-01, -0.10161}, {-0.68864E-01, -0.99490E-01}}, /* tau_np_p_sgl[2][2] */
  {{-0.14070, -0.20537}, {-0.18898, -0.28569}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-0.43046E-01, -0.43046E-01}, {-0.46243E-01, -0.46243E-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-0.30836E-02, -0.27445E-02}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-0.97936E-02, -0.22315E-01}, {-0.87970E-02, -0.19591E-01}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-0.82486E-01, -0.17592}, {-0.99536E-01, -0.21200}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {-0.64564E-02, -0.24297E-02}, /* tau_k_1_p[2] */
  {-0.52206E-02, -0.29522E-02}, /* tau_k_2_p[2] */
  {-0.57495E-02, -0.27099E-02}, /* tau_k_1_n[2] */
  {-0.53495E-02, -0.35106E-02}, /* tau_k_2_n[2] */
  {-0.33081E-02, -0.17945E-02}, /* sigma_k[2] */
  -0.48524E-01, /* fk2 */
  /* outlier parameters */
  -3.5299, /* sig_ol */
  -0.18268E-02, /* fol_sgl */
  -0.81875E-04, /* fol_mul */
};

dtres_param_t dtres_param_2003lp03 = {
  /* R_det */
  /*  multiple track vertex */
  {0.970058, 0.086367}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.787972, 0.0347955},/* Sasc[2] */
  0.0586515, 0.924236, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.924087, /* Smn_rec */
  0.00000, /* Stl_rec */
  0.00000, /* ftl_rec */
  0.924087, /* Smn_asc */
  0.00000, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.687531, 0.687531}, {0.707721, 0.707721}}, /* fd_np_sgl[2][2] */
  {0.790758, 0.803726}, /* fp_np_sgl[2] */
  {{0.452823, 1.32027}, {0.300821, 1.13042}}, /* tau_np_p_sgl[2][2] */
  {{0.295111, 1.47565}, {0.0159646, 1.75911}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.55449, 0.55449}, {0.489829, 0.489829}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.929445, 0.935395}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.0557029, 0.724763}, {0.0518227, 0.629883}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.0253055, 1.28048}, {-0.0692348, 1.34919}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.800358e-01, 1.284213e-02}, /* tau_k_1_p[2] */
  {7.490356e-02, 3.666411e-02}, /* tau_k_2_p[2] */
  {1.803169e-01, 1.195082e-02}, /* tau_k_1_n[2] */
  {7.546680e-02, 3.480336e-02}, /* tau_k_2_n[2] */
  {-1.822232e-02, 4.673038e-02}, /* sigma_k[2] */
  5.756437e-01, /* fk2 */
  /* outlier parameters */
  38.7854, /* sig_ol */
  0.0312508, /* fol_sgl */
  0.000206659, /* fol_mul */
};

dtres_param_t dtres_param_2003lp03_mc = {
  /* R_det */
  /*  multiple track vertex */
  {1.02785, 0.0733288}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.922899, 0.0236024},/* Sasc[2] */
  0.060349, 1.04406, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.999832, /* Smn_rec */
  0.00000, /* Stl_rec */
  0.00000, /* ftl_rec */
  0.999832, /* Smn_asc */
  0.00000, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.7725, 0.7725}, {0.787286, 0.787286}}, /* fd_np_sgl[2][2] */
  {0.790758, 0.803726}, /* fp_np_sgl[2] */
  {{0.452823, 1.32027}, {0.300821, 1.13042}}, /* tau_np_p_sgl[2][2] */
  {{0.295111, 1.47565}, {0.0159646, 1.75911}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.619196, 0.619196}, {0.661041, 0.661041}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.929445, 0.935395}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.0557029, 0.724763}, {0.0518227, 0.629883}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.0253055, 1.28048}, {-0.0692348, 1.34919}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.800358e-01, 1.284213e-02}, /* tau_k_1_p[2] */
  {7.490356e-02, 3.666411e-02}, /* tau_k_2_p[2] */
  {1.803169e-01, 1.195082e-02}, /* tau_k_1_n[2] */
  {7.546680e-02, 3.480336e-02}, /* tau_k_2_n[2] */
  {-1.822232e-02, 4.673038e-02}, /* sigma_k[2] */
  5.756437e-01, /* fk2 */
  /* outlier parameters */
  36.2549, /* sig_ol */
  0.0189902, /* fol_sgl */
  6.4485e-05, /* fol_mul */
};

dtres_param_t dtres_param_posi_error_2003lp03 = {
  /* R_det */
  /*  multiple track vertex */
  {0.0533451, 0.00320566}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.031584, 0.00163838},/* Sasc[2] */
  0.00461318, 0.0438257, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.0226567, /* Smn_rec */
  0.00000, /* Stl_rec */
  0.00000, /* ftl_rec */
  0.0226567, /* Smn_asc */
  0.00000, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.0211739, 0.0211739}, {0.0299131, 0.0299131}}, /* fd_np_sgl[2][2] */
  {0.00908806, 0.0105695}, /* fp_np_sgl[2] */
  {{0.0463331, 0.0643792}, {0.0398841, 0.0605884}}, /* tau_np_p_sgl[2][2] */
  {{0.0857982, 0.141208}, {0.0816327, 0.152409}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.0231622, 0.0231622}, {0.032111, 0.032111}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {0.00250954, 0.00264214}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{0.00633286, 0.0145159}, {0.00571942, 0.0134357}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{0.0316072, 0.063646}, {0.0288697, 0.0648656}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {4.010734e-03, 1.291296e-03}, /* tau_k_1_p[2] */
  {2.212691e-03, 1.006885e-03}, /* tau_k_2_p[2] */
  {4.037715e-03, 1.455432e-03}, /* tau_k_1_n[2] */
  {2.215337e-03, 1.136337e-03}, /* tau_k_2_n[2] */
  {1.712383e-03, 9.324146e-04}, /* sigma_k[2] */
  2.576726e-02, /* fk2 */
  /* outlier parameters */
  1.64817, /* sig_ol */
  0.00108591, /* fol_sgl */
  6.05913e-05, /* fol_mul */
};

dtres_param_t dtres_param_nega_error_2003lp03 = {
  /* R_det */
  /*  multiple track vertex */
  {-0.0544444, -0.0031762}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {-0.032151, -0.00163149},/* Sasc[2] */
  -0.00438922, -0.042621, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  -0.0228761, /* Smn_rec */
  0.00000, /* Stl_rec */
  0.00000, /* ftl_rec */
  -0.0228761, /* Smn_asc */
  0.00000, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{-0.022269, -0.022269}, {-0.0313874, -0.0313874}}, /* fd_np_sgl[2][2] */
  {-0.00938177, -0.0110095}, /* fp_np_sgl[2] */
  {{-0.0459104, -0.062931}, {-0.0392732, -0.0590955}}, /* tau_np_p_sgl[2][2] */
  {{-0.0848693, -0.133996}, {-0.0805739, -0.144376}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-0.0243739, -0.0243739}, {-0.0336518, -0.0336518}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-0.00258786, -0.0027415}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-0.00633083, -0.0143718}, {-0.00571115, -0.0132913}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-0.0316049, -0.0614472}, {-0.029099, -0.0623312}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {-3.786448e-03, -1.347454e-03}, /* tau_k_1_p[2] */
  {-2.263418e-03, -9.866002e-04}, /* tau_k_2_p[2] */
  {-3.813974e-03, -1.520782e-03}, /* tau_k_1_n[2] */
  {-2.264278e-03, -1.117868e-03}, /* tau_k_2_n[2] */
  {-1.703007e-03, -9.464530e-04}, /* sigma_k[2] */
  -2.595315e-02, /* fk2 */
  /* outlier parameters */
  -1.5114, /* sig_ol */
  -0.00105297, /* fol_sgl */
  -5.19202e-05, /* fol_mul */
};

dtres_param_t dtres_param_svd1_ichep04 = {
  /* R_det */
  /*  multiple track vertex */
  {0.970106, 0.0863648}, /* Srec[2] */
  {0., 0.0}, /* ftl_rec_mlt[2] */
  0.,  /* Stl_rec_mlt */
  //  {0.799592, 0.0349734,0.0602796,0.929777}, /* Sasc[4] */
  {0.799592, 0.0349734},/* Sasc[2] */
  0.0602796, 0.92977, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.923743, /* Smn_rec */
  0., /* Stl_rec */
  0., /* ftl_rec */
  0.923743, /* Smn_asc */
  0., /* Stl_asc */
  0., /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.691997, 0.691997}, {0.71177, 0.71177}}, /* fd_np_sgl[2][2] */
  {0.790876, 0.803436}, /* fp_np_sgl[2] */
  {{0.457941,1.32583}, {0.301517,1.13921}}, /* tau_np_p_sgl[2][2] */
  {{0.292904,1.49776}, {0.0127143,1.77989}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.556886, 0.556886}, {0.491796, 0.491796}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.939597, 0.945704}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.0498904,0.704206}, {0.0474975,0.609198}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.0439374, 1.28731}, {-0.0512509,1.35101}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.800358e-01, 1.284213e-02}, /* tau_k_1_p[2] */
  {7.490356e-02, 3.666411e-02}, /* tau_k_2_p[2] */
  {1.803169e-01, 1.195082e-02}, /* tau_k_1_n[2] */
  {7.546680e-02, 3.480336e-02}, /* tau_k_2_n[2] */
  {-1.822232e-02, 4.673038e-02}, /* sigma_k[2] */
  5.756437e-01, /* fk2 */
  /* outlier parameters */
  38.8948, /* sig_ol */
  0.0311556, /* fol_sgl */
  0.000206552, /* fol_mul */
};

dtres_param_t dtres_param_svd1_ichep04_mc = {
  /* R_det */
  /*  multiple track vertex */
  {1.03596, 0.0733278}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
  0.0000, /* Stl_rec_mlt */
  {0.901069, 0.0282578},/* Sasc[2] */
  0.0710005, 1.00563, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.997662, /* Smn_rec */
  0.00000, /* Stl_rec */
  0.00000, /* ftl_rec */
  0.997662, /* Smn_asc */
  0.00000, /* Stl_asc */
  0.00000, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.76524, 0.76524}, {0.780648, 0.780648}}, /* fd_np_sgl[2][2] */
  {0.790876, 0.803436}, /* fp_np_sgl[2] */
  {{0.457941, 1.32583}, {0.301517, 1.13921}}, /* tau_np_p_sgl[2][2] */
  {{0.292904, 1.49776}, {0.0127143, 1.77989}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.610748, 0.610748}, {0.634781, 0.634781}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.939597, 0.945704}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.0498904, 0.704206}, {0.0474975, 0.609198}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.0439374, 1.28731}, {-0.0512509, 1.35101}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.800358e-01, 1.284213e-02}, /* tau_k_1_p[2] */
  {7.490356e-02, 3.666411e-02}, /* tau_k_2_p[2] */
  {1.803169e-01, 1.195082e-02}, /* tau_k_1_n[2] */
  {7.546680e-02, 3.480336e-02}, /* tau_k_2_n[2] */
  {-1.822232e-02, 4.673038e-02}, /* sigma_k[2] */
  5.756437e-01, /* fk2 */
  /* outlier parameters */
  36.084, /* sig_ol */
  0.0190383, /* fol_sgl */
  6.17239e-05 /* fol_mul */
};

dtres_param_t dtres_param_posi_error_svd1_ichep04 = {
  /* R_det */
  /*  multiple track vertex */
    {+0.142344713, +0.011705448}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
    0.0000, /* Stl_rec_mlt */
  {+0.183069076, +0.001844738},/* Sasc[2] */
  +0.010900729, +0.165881307, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
    0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
    +0.218238054, /* Smn_rec */
    0.00000, /* Stl_rec */
    0.00000, /* ftl_rec */
    +0.218238054, /* Smn_asc */
    0.00000, /* Stl_asc */
    0.00000, /* ftl_asc */
    /* Rnp parameters */
    /*  single track vertex */
    /* index ... 0 <--> B0_B0B */
    /*           1 <--> B+_B-  */
  {{+0.035120518, +0.035120518}, {+0.061178322, +0.061178322}}, /* fd_np_sgl[2][2] */
    {+0.00907758, +0.0105837}, /* fp_np_sgl[2] */
    {{+0.0466642, +0.0644879}, {+0.0401322, +0.060925}}, /* tau_np_p_sgl[2][2] */
    {{+0.0865975, +0.141832}, {+0.0822161, +0.152761}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{+0.040280979, +0.040280979}, {+0.049871221, +0.049871221}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
    {+0.00235392, +0.00249127}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
    {{+0.00600101, +0.0139987}, {+0.00541068, +0.0129163}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
    {{+0.0348137, +0.0707513}, {+0.0317484, +0.0721001}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
    {4.010734E-03, 1.291296E-03}, /* tau_k_1_p[2] */
    {2.212691E-03, 1.006885E-03}, /* tau_k_2_p[2] */
    {4.037715E-03, 1.455432E-03}, /* tau_k_1_n[2] */
    {2.215337E-03, 1.136337E-03}, /* tau_k_2_n[2] */
    {1.712383E-03, 9.324146E-04}, /* sigma_k[2] */
    2.576726E-02, /* fk2 */
    /* outlier parameters */
    +1.758414661, /* sig_ol */
    +0.00304396, /* fol_sgl */
    +0.000122048, /* fol_mul */
};

dtres_param_t dtres_param_nega_error_svd1_ichep04 = {
  /* R_det */
  /*  multiple track vertex */
    {-0.069003508, -0.010950977}, /* Srec[2] */
  {0.0000, 0.0}, /* ftl_rec_mlt[2] */
    0.0000, /* Stl_rec_mlt */
  {-0.037272377, -0.006118318},/* Sasc[2] */
  -0.015982148, -0.079279803, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
    0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
    -0.031568739, /* Smn_rec */
    0.00000, /* Stl_rec */
    0.00000, /* ftl_rec */
    -0.031568739, /* Smn_asc */
    0.00000, /* Stl_asc */
    0.00000, /* ftl_asc */
    /* Rnp parameters */
    /*  single track vertex */
    /* index ... 0 <--> B0_B0B */
    /*           1 <--> B+_B-  */
  {{-0.03654336, -0.03654336}, {-0.05110292, -0.05110292}}, /* fd_np_sgl[2][2] */
    {-0.00927331, -0.0108754}, /* fp_np_sgl[2] */
    {{-0.0460243, -0.0629668}, {-0.0394386, -0.0592732}}, /* tau_np_p_sgl[2][2] */
    {{-0.084934, -0.134368}, {-0.0804162, -0.144111}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-0.040693549, -0.040693549}, {-0.050384555, -0.050384555}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
    {-0.00241696, -0.00257006}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
    {{-0.00599047, -0.0138339}, {-0.00540709, -0.0127585}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
    {{-0.0345698, -0.0682507}, {-0.0316477, -0.0692694}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
    {-3.786448e-03, -1.347454e-03}, /* tau_k_1_p[2] */
    {-2.263418e-03, -9.866002e-04}, /* tau_k_2_p[2] */
    {-3.813974e-03, -1.520782e-03}, /* tau_k_1_n[2] */
    {-2.264278e-03, -1.117868e-03}, /* tau_k_2_n[2] */
    {-1.703007e-03, -9.464530e-04}, /* sigma_k[2] */
    -2.595315e-02, /* fk2 */    
    /* outlier parameters */
    -12.50853513, /* sig_ol */
    -0.006035018, /* fol_sgl */
    -8.10812E-05, /* fol_mul */
};

dtres_param_t dtres_param_svd2_ichep04_v0 = {
  /* R_det */
  /*  multiple track vertex */
  {0.687154, 0.129697}, /* Srec[2] */
  {0.101263, 0.0}, /* ftl_rec_mlt[2] */
  4.78078,  /* Stl_rec_mlt */
  {0.705497, 0.0496698},/* Sasc[2] */
  0.0980723, 1.067, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  1.10977, /* Smn_rec */
  10.2716, /* Stl_rec */
  0.0339335, /* ftl_rec */
  1.10977, /* Smn_asc */
  10.2716, /* Stl_asc */
  0.0339335, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.787241, 0.787241}, {0.742299, 0.742299}}, /* fd_np_sgl[2][2] */
  {0.774089, 0.769273}, /* fp_np_sgl[2] */
  {{1.31876, 0.0}, {0.934049, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{0.638336, 0.0}, {0.519405, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.46029, 0.46029}, {0.402308, 0.402308}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.773972, 0.761986}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.132456, 0.679806}, {0.122075, 0.572431}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.151427, 0.673289}, {0.107643, 0.629359}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.692797e-01, 1.336589e-02}, /* tau_k_1_p[2] */
  {6.620223e-02, 3.723578e-02}, /* tau_k_2_p[2] */
  {1.713830e-01, 1.307457e-02}, /* tau_k_1_n[2] */
  {6.675977e-02, 3.423165e-02}, /* tau_k_2_n[2] */
  {-1.247497e-02, 4.956063e-02}, /* sigma_k[2] */
  5.094241e-01, /* fk2 */
  /* outlier parameters */
  44.2129, /* sig_ol */
  0.0119662, /* fol_sgl */
  0.000314807, /* fol_mul */
};

dtres_param_t dtres_param_svd2_ichep04_mc = {
  /* R_det */
  /*  multiple track vertex */
  {0.627412, 0.106106}, /* Srec[2] */
  {0.100731, 0.0}, /* ftl_rec_mlt[2] */
  6.23692, /* Stl_rec_mlt */
  {0.862452, 0.0206427},/* Sasc[2] */
  0.0980723, 0.922342, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.994838, /* Smn_rec */
  4.55494, /* Stl_rec */
  0.068628, /* ftl_rec */
  0.994838, /* Smn_asc */
  4.55494, /* Stl_asc */
  0.068628, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.759525, 0.759525}, {0.750635, 0.750635}}, /* fd_np_sgl[2][2] */
  {0.774089, 0.769273}, /* fp_np_sgl[2] */
  {{1.31876, 0.0}, {0.934049, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{0.638336, 0.0}, {0.519405, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.340936, 0.340936}, {0.313528, 0.313528}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.773972, 0.761986}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.132456, 0.679806}, {0.122075, 0.572431}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.151427, 0.673289}, {0.107643, 0.629359}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.692797e-01, 1.336589e-02}, /* tau_k_1_p[2] */
  {6.620223e-02, 3.723578e-02}, /* tau_k_2_p[2] */
  {1.713830e-01, 1.307457e-02}, /* tau_k_1_n[2] */
  {6.675977e-02, 3.423165e-02}, /* tau_k_2_n[2] */
  {-1.247497e-02, 4.956063e-02}, /* sigma_k[2] */
  5.094241e-01, /* fk2 */
  /* outlier parameters */
  32.2588, /* sig_ol */
  0.0156605, /* fol_sgl */
  0.000270744, /* fol_mul */
};

dtres_param_t dtres_param_posi_error_svd2_ichep04_v0 = {
  /* R_det */
  /*  multiple track vertex */
    {2.941638E-01, 2.680165E-02}, /* Srec[2] */
  {1.934048E-01, 0.0}, /* ftl_rec_mlt[2] */
    5.805351E-01, /* Stl_rec_mlt */
  {4.926124E-01, 1.173910E-02},/* Sasc[2] */
  2.257050E-03, 9.250735E-02, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
    0.0000, /* Stl_asc_mlt */
    /*  single track vertex */
    3.651375E-01, /* Smn_rec */
    2.177915E+00, /* Stl_rec */
    2.581997E-02, /* ftl_rec */
    3.651375E-01, /* Smn_asc */
    2.177915E+00, /* Stl_asc */
    2.581997E-02, /* ftl_asc */
    /* Rnp parameters */
    /*  single track vertex */
    /* index ... 0 <--> B0_B0B */
    /*           1 <--> B+_B-  */
  {{2.675478E-02, 2.675478E-02}, {5.109274E-02, 5.109274E-02}}, /* fd_np_sgl[2][2] */
    {2.078770E-02, 5.195840E-02}, /* fp_np_sgl[2] */
    {{4.049970E-02, 0.0}, {7.348070E-02, 0.0}}, /* tau_np_p_sgl[2][2] */
    {{6.226030E-02, 0.0}, {1.414230E-01, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{6.035001E-02, 6.035001E-02}, {6.583662E-02, 6.583662E-02}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
    {4.424370E-03, 6.534310E-03}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
    {{4.567270E-03, 1.290090E-02}, {6.496420E-03, 1.597400E-02}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
    {{8.294310E-03, 1.911180E-02}, {1.183310E-02, 2.884200E-02}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {4.258054e-03, 1.489901e-03}, /* tau_k_1_p[2] */
  {2.882638e-03, 1.424035e-03}, /* tau_k_2_p[2] */
  {4.353473e-03, 1.597366e-03}, /* tau_k_1_n[2] */    
  {2.934835e-03, 1.574420e-03}, /* tau_k_2_n[2] */
  {1.922392e-03, 1.107931e-03}, /* sigma_k[2] */      
  3.112788e-02, /* fk2 */
    /* outlier parameters */
    6.763932E+00, /* sig_ol */
    3.837933E-03, /* fol_sgl */
    3.566260E-04, /* fol_mul */
};

dtres_param_t dtres_param_nega_error_svd2_ichep04_v0 = {
  /* R_det */
  /*  multiple track vertex */
    {-2.638597E-01, -5.283685E-02}, /* Srec[2] */
  {-4.831824E-02, 0.0}, /* ftl_rec_mlt[2] */
    -1.166409E+00, /* Stl_rec_mlt */
  {-6.293215E-02, -2.473299E-02},/* Sasc[2] */
  -2.209080E-03, -1.341872E-01, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
    0.0000, /* Stl_asc_mlt */
    /*  single track vertex */
    -4.165114E-02, /* Smn_rec */
    -2.847818E+00, /* Stl_rec */
    -4.952554E-03, /* ftl_rec */
    -4.165114E-02, /* Smn_asc */
    -2.847818E+00, /* Stl_asc */
    -4.952554E-03, /* ftl_asc */
    /* Rnp parameters */
    /*  single track vertex */
    /* index ... 0 <--> B0_B0B */
    /*           1 <--> B+_B-  */
  {{-2.815595E-02, -2.815595E-02}, {-5.495795E-02, -5.495795E-02}}, /* fd_np_sgl[2][2] */
    {-2.069210E-02, -5.420790E-02}, /* fp_np_sgl[2] */
    {{-3.898860E-02, 0.0}, {-6.915910E-02, 0.0}}, /* tau_np_p_sgl[2][2] */
    {{-5.636480E-02, 0.0}, {-1.103080E-01, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-6.768774E-02, -6.768774E-02}, {-9.760931E-02, -9.760931E-02}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
    {-4.399570E-03, -6.532750E-03}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
    {{-4.577700E-03, -1.261650E-02}, {-6.372320E-03, -1.575800E-02}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
    {{-8.113570E-03, -1.851060E-02}, {-1.167420E-02, -2.760620E-02}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {-3.942710e-03, -1.574283e-03}, /* tau_k_1_p[2] */
  {-2.904957e-03, -1.398102e-03}, /* tau_k_2_p[2] */
  {-4.036754e-03, -1.700615e-03}, /* tau_k_1_n[2] */
  {-2.948232e-03, -1.554997e-03}, /* tau_k_2_n[2] */
  {-2.097279e-03, -1.105820e-03}, /* sigma_k[2] */
  -3.048329e-02, /* fk2 */
  /* outlier parameters */
    -1.985724E+01, /* sig_ol */
    -1.727439E-03, /* fol_sgl */
    -1.130770E-04, /* fol_mul */
};

dtres_param_t dtres_param_svd2_ichep04 = {
  /* R_det */
  /*  multiple track vertex */
  {0.702034, 0.132134}, /* Srec[2] */
  {0.0948268, 0.0}, /* ftl_rec_mlt[2] */
  4.86672, /* Stl_rec_mlt */
  {0.69771, 0.0499168},/* Sasc[2] */
  0.0980723, 1.06446, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  1.10593, /* Smn_rec */
  10.0092, /* Stl_rec */
  0.0345081, /* ftl_rec */
  1.10593, /* Smn_asc */
  10.0092, /* Stl_asc */
  0.0345081, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.789743, 0.789743}, {0.741749, 0.741749}}, /* fd_np_sgl[2][2] */
  {0.774089, 0.769273}, /* fp_np_sgl[2] */
  {{1.31876, 0.0}, {0.934049, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{0.638336, 0.0}, {0.519405, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.457612, 0.457612}, {0.400228, 0.400228}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.773972, 0.761986}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.132456, 0.679806}, {0.122075, 0.572431}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.151427, 0.673289}, {0.107643, 0.629359}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.692797e-01, 1.336589e-02}, /* tau_k_1_p[2] */
  {6.620223e-02, 3.723578e-02}, /* tau_k_2_p[2] */
  {1.713830e-01, 1.307457e-02}, /* tau_k_1_n[2] */
  {6.675977e-02, 3.423165e-02}, /* tau_k_2_n[2] */
  {-1.247497e-02, 4.956063e-02}, /* sigma_k[2] */
  5.094241e-01, /* fk2 */
  /* outlier parameters */
  44.5257, /* sigma_ol */
  0.0120795, /* f_ol_sgl */
  0.000310578, /* f_ol_mul */
};


dtres_param_t dtres_param_posi_error_svd2_ichep04 = {
  /* R_det */
  /*  multiple track vertex */
  {0.29399208, 0.026666302}, /* Srec[2] */
  {0.190664696, 0.0}, /* ftl_rec_mlt[2] */
  0.582523132, /* Stl_rec_mlt */
  {0.492627468, 0.011737733},/* Sasc[2] */
  0.00225705, 0.092467982, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.365152525, /* Smn_rec */
  2.171407295, /* Stl_rec */
  0.025844094, /* ftl_rec */
  0.365152525, /* Smn_asc */
  2.171407295, /* Stl_asc */
  0.025844094, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.026743868, 0.026743868}, {0.050891012, 0.050891012}}, /* fd_np_sgl[2][2] */
  {0.0207877, 0.0519584}, /* fp_np_sgl[2] */
  {{0.0404997, 0.0}, {0.0734807, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{0.0622603, 0.0}, {0.141423, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.060445007, 0.060445007}, {0.065735575, 0.065735575}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {0.00442437, 0.00653431}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{0.00456727, 0.0129009}, {0.00649642, 0.015974}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{0.00829431, 0.0191118}, {0.0118331, 0.028842}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {4.258054e-03, 1.489901e-03}, /* tau_k_1_p[2] */
  {2.882638e-03, 1.424035e-03}, /* tau_k_2_p[2] */
  {4.353473e-03, 1.597366e-03}, /* tau_k_1_n[2] */    
  {2.934835e-03, 1.574420e-03}, /* tau_k_2_n[2] */
  {1.922392e-03, 1.107931e-03}, /* sigma_k[2] */      
  3.112788e-02, /* fk2 */
  /* outlier parameters */
  6.641605717, /* sigma_ol */
  0.00383841, /* f_ol_sgl */
  0.000355688, /* f_ol_mul */
};

dtres_param_t dtres_param_nega_error_svd2_ichep04 = {
  /* R_det */
  /*  multiple track vertex */
  {-0.261707956, -0.052112185}, /* Srec[2] */
  {-0.047594054, 0.0}, /* ftl_rec_mlt[2] */
  -1.1357572, /* Stl_rec_mlt */
  {-0.062999125, -0.024735712},/* Sasc[2] */
  -0.00220908, -0.134424635, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  -0.041645712, /* Smn_rec */
  -2.840824255, /* Stl_rec */
  -0.005034424, /* ftl_rec */
  -0.041645712, /* Smn_asc */
  -2.840824255, /* Stl_asc */
  -0.005034424, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{-0.028350522, -0.028350522}, {-0.055673858, -0.055673858}}, /* fd_np_sgl[2][2] */
  {-0.0206921, -0.0542079}, /* fp_np_sgl[2] */
  {{-0.0389886, 0.0}, {-0.0691591, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{-0.0563648, 0.0}, {-0.110308, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-0.068105943, -0.068105943}, {-0.098082123, -0.098082123}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-0.00439957, -0.00653275}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-0.0045777, -0.0126165}, {-0.00637232, -0.015758}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-0.00811357, -0.0185106}, {-0.0116742, -0.0276062}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {-3.942710e-03, -1.574283e-03}, /* tau_k_1_p[2] */
  {-2.904957e-03, -1.398102e-03}, /* tau_k_2_p[2] */
  {-4.036754e-03, -1.700615e-03}, /* tau_k_1_n[2] */
  {-2.948232e-03, -1.554997e-03}, /* tau_k_2_n[2] */
  {-2.097279e-03, -1.105820e-03}, /* sigma_k[2] */
  -3.048329e-02, /* fk2 */
  /* outlier parameters */
  -19.86010176, /* sigma_ol */
  -0.001728323, /* f_ol_sgl */
  -0.000112877, /* f_ol_mul */
};

dtres_param_t dtres_param_svd2_lp05 = {
  /* R_det */
  /*  multiple track vertex */
  {0.65244, 0.125157}, /* Srec[2] */
  {0.118444, 0.0}, /* ftl_rec_mlt[2] */
  4.78117, /* Stl_rec_mlt */
  {0.75968, 0.0453833},/* Sasc[2] */
  0.0980723, 0.950182, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  1.02207, /* Smn_rec */
  4.18482, /* Stl_rec */
  0.0870378, /* ftl_rec */
  1.02207, /* Smn_asc */
  4.18482, /* Stl_asc */
  0.0870378, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.776802, 0.776802}, {0.709836, 0.709836}}, /* fd_np_sgl[2][2] */
  {0.774089, 0.769273}, /* fp_np_sgl[2] */
  {{1.31876, 0.0}, {0.934049, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{0.638336, 0.0}, {0.519405, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.393355, 0.393355}, {0.326565, 0.326565}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {0.773972, 0.761986}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{0.132456, 0.679806}, {0.122075, 0.572431}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{0.151427, 0.673289}, {0.107643, 0.629359}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {1.692797e-01, 1.336589e-02}, /* tau_k_1_p[2] */
  {6.620223e-02, 3.723578e-02}, /* tau_k_2_p[2] */
  {1.713830e-01, 1.307457e-02}, /* tau_k_1_n[2] */
  {6.675977e-02, 3.423165e-02}, /* tau_k_2_n[2] */
  {-1.247497e-02, 4.956063e-02}, /* sigma_k[2] */
  5.094241e-01, /* fk2 */
  /* outlier parameters */
  35.4457, /* sigma_ol */
  0.0180001, /* f_ol_sgl */
  0.000361175, /* f_ol_mul */
};

dtres_param_t dtres_param_posi_error_svd2_lp05 = {
  /* R_det */
  /*  multiple track vertex */
  {0.230376977, 0.031248506}, /* Srec[2] */
  {0.047948743, 0.0}, /* ftl_rec_mlt[2] */
  0.459171004, /* Stl_rec_mlt */
  {0.499735846, 0.011027941},/* Sasc[2] */
  0.00225705, 0.099621476, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  0.448388223, /* Smn_rec */
  4.079661237, /* Stl_rec */
  0.028493377, /* ftl_rec */
  0.448388223, /* Smn_asc */
  4.079661237, /* Stl_asc */
  0.028493377, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{0.039552949, 0.039552949}, {0.094057747, 0.094057747}}, /* fd_np_sgl[2][2] */
  {0.0207877, 0.0519584}, /* fp_np_sgl[2] */
  {{0.0404997, 0.0}, {0.0734807, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{0.0622603, 0.0}, {0.141423, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{0.074804174, 0.074804174}, {0.120571726, 0.120571726}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {0.00442437, 0.00653431}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{0.00456727, 0.0129009}, {0.00649642, 0.015974}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{0.00829431, 0.0191118}, {0.0118331, 0.028842}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.0, 0.0}, /* tau_k_1_p[2] */
  {0.0, 0.0}, /* tau_k_2_p[2] */
  {0.0, 0.0}, /* tau_k_1_n[2] */
  {0.0, 0.0}, /* tau_k_2_n[2] */
  {0.0, 0.0}, /* sigma_k[2] */
  0.0, /* fk2 */
  /* outlier parameters */
  7.812762694, /* sigma_ol */
  0.002226117, /* f_ol_sgl */
  0.000203678, /* f_ol_mul */
};

dtres_param_t dtres_param_nega_error_svd2_lp05 = {
  /* R_det */
  /*  multiple track vertex */
  {-0.114260797, -0.017523041}, /* Srec[2] */
  {-0.028109606, 0.0}, /* ftl_rec_mlt[2] */
  -0.606904348, /* Stl_rec_mlt */
  {-0.105422757, -0.020903425},/* Sasc[2] */
  -0.00220908, -0.133797923, /* Snp, Snp_global */
  {0.0000, 0.0}, /* ftl_asc_mlt[2] */
  0.0000, /* Stl_asc_mlt */
  /*  single track vertex */
  -0.052976488, /* Smn_rec */
  -0.831459499, /* Stl_rec */
  -0.031392851, /* ftl_rec */
  -0.052976488, /* Smn_asc */
  -0.831459499, /* Stl_asc */
  -0.031392851, /* ftl_asc */
  /* Rnp parameters */
  /*  single track vertex */
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  {{-0.037217378, -0.037217378}, {-0.125742306, -0.125742306}}, /* fd_np_sgl[2][2] */
  {-0.0206921, -0.0542079}, /* fp_np_sgl[2] */
  {{-0.0389886, 0.0}, {-0.0691591, 0.0}}, /* tau_np_p_sgl[2][2] */
  {{-0.0563648, 0.0}, {-0.110308, 0.0}}, /* tau_np_n_sgl[2][2] */
  /*  multiple track vertex  */
  {{-0.087882144, -0.087882144}, {-0.158834765, -0.158834765}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-0.00439957, -0.00653275}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-0.0045777, -0.0126165}, {-0.00637232, -0.015758}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-0.00811357, -0.0185106}, {-0.0116742, -0.0276062}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (For partial reconstracted event) */
  {0.0, 0.0}, /* tau_k_1_p[2] */
  {0.0, 0.0}, /* tau_k_2_p[2] */
  {0.0, 0.0}, /* tau_k_1_n[2] */
  {0.0, 0.0}, /* tau_k_2_n[2] */
  {0.0, 0.0}, /* sigma_k[2] */
  0.0, /* fk2 */
  /* outlier parameters */
  -10.62940498, /* sigma_ol */
  -0.002825928, /* f_ol_sgl */
  -0.000114168, /* f_ol_mul */
};

dtres_param_t dtres_param_center_svd1_ichep06_mc = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+1.086620e+00, +1.359900e-02}, /* Srec[2] */
  {+0.000000e+00, +4.470700e-04}, /* ftl_rec_mlt[2] */
  +6.648120e+00, /* Stl_rec_mlt */
  {+8.967260e-01, +5.098590e-03}, /* Sasc[2] */
  +2.020690e-02, +9.910830e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +9.835630e-01, /* Smn_rec */
  +1.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +9.835630e-01, /* Smn_asc */
  +1.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+6.796420e-01, +1.000000e+00}, {+6.956790e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+7.339510e-01, +7.454320e-01}, /* fp_np_sgl[2] */
  {{+6.412400e-01, +1.044720e+00}, {+3.135950e-01, +1.057200e+00}}, /* tau_np_p_sgl[2][2] */
  {{+3.137720e-01, +1.199360e+00}, {+1.957760e-01, +1.307330e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+3.491480e-01, +5.123690e-01}, {+3.149430e-01, +4.385360e-01}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {+7.916810e-01, +7.748680e-01}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{+7.112550e-02, +6.072890e-01}, {+7.004320e-02, +5.297430e-01}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{+6.233280e-02, +6.267190e-01}, {+1.108210e-02, +6.169610e-01}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+1.675300e-01, +1.157618e-02}, /* tau_k_1_p[2] */
  {+6.642841e-02, +3.364998e-02}, /* tau_k_2_p[2] */
  {+1.674884e-01, +1.180294e-02}, /* tau_k_1_n[2] */
  {+6.625765e-02, +3.254376e-02}, /* tau_k_2_n[2] */
  {-1.955903e-02, +4.648172e-02}, /* sigma_k[2] */
  +5.696577e-01, /* fk2 */
  /* outlier parameters */
  +3.616460e+01, /* sigma_ol */
  +4.104210e-02, /* f_ol_sgl */
  +2.661410e-04, /* f_ol_mul */
};


dtres_param_t dtres_param_center_svd2_ichep06_mc = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+6.901070e-01, +2.948450e-02}, /* Srec[2] */
  {+9.984920e-02, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +5.575640e+00, /* Stl_rec_mlt */
  {+8.855180e-01, +4.990170e-03}, /* Sasc[2] */
  +3.356030e-02, +9.244060e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.025930e+00, /* Smn_rec */
  +4.304330e+00, /* Stl_rec */
  +6.702860e-02, /* ftl_rec */
  +1.025930e+00, /* Smn_asc */
  +4.304330e+00, /* Stl_asc */
  +6.702860e-02, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.172600e-01, +1.000000e+00}, {+7.201820e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.409870e-01, +8.340980e-01}, /* fp_np_sgl[2] */
  {{+1.450340e+00, +0.000000e+00}, {+8.861100e-01, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+1.001860e+00, +0.000000e+00}, {+5.190220e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+2.564670e-01, +4.161530e-01}, {+2.987660e-01, +4.261770e-01}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {+7.702640e-01, +7.619030e-01}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{+1.106250e-01, +6.400080e-01}, {+1.005250e-01, +5.708080e-01}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{+1.285160e-01, +5.893880e-01}, {+9.230920e-02, +5.605140e-01}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+1.658325e-01, +1.200647e-02}, /* tau_k_1_p[2] */
  {+6.344224e-02, +3.424377e-02}, /* tau_k_2_p[2] */
  {+1.675831e-01, +1.193076e-02}, /* tau_k_1_n[2] */
  {+6.454622e-02, +3.200147e-02}, /* tau_k_2_n[2] */
  {-1.711954e-02, +5.607785e-02}, /* sigma_k[2] */
  +5.123095e-01, /* fk2 */
  /* outlier parameters */
  +3.032360e+01, /* sigma_ol */
  +2.326310e-02, /* f_ol_sgl */
  +1.251280e-04, /* f_ol_mul */
};


dtres_param_t dtres_param_center_svd1_ichep06_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+9.103520e-01, +2.378820e-02}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+6.668790e-01, +9.948560e-03}, /* Sasc[2] */
  +2.020690e-02, +9.540690e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +9.062250e-01, /* Smn_rec */
  +1.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +9.062250e-01, /* Smn_asc */
  +1.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+5.975640e-01, +1.000000e+00}, {+6.349990e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+7.339510e-01, +7.454320e-01}, /* fp_np_sgl[2] */
  {{+6.412400e-01, +1.044720e+00}, {+3.135950e-01, +1.057200e+00}}, /* tau_np_p_sgl[2][2] */
  {{+3.137720e-01, +1.199360e+00}, {+1.957760e-01, +1.307330e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+2.505390e-01, +5.088590e-01}, {+1.657090e-01, +2.407730e-01}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {+7.916810e-01, +7.748680e-01}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{+7.112550e-02, +6.072890e-01}, {+7.004320e-02, +5.297430e-01}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{+6.233280e-02, +6.267190e-01}, {+1.108210e-02, +6.169610e-01}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+1.675300e-01, +1.157618e-02}, /* tau_k_1_p[2] */
  {+6.642841e-02, +3.364998e-02}, /* tau_k_2_p[2] */
  {+1.674884e-01, +1.180294e-02}, /* tau_k_1_n[2] */
  {+6.625765e-02, +3.254376e-02}, /* tau_k_2_n[2] */
  {-1.955903e-02, +4.648172e-02}, /* sigma_k[2] */
  +5.696577e-01, /* fk2 */
  /* outlier parameters */
  +3.742850e+01, /* sigma_ol */
  +4.216880e-02, /* f_ol_sgl */
  +2.515330e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_posi_error_svd1_ichep06_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+4.628298e-01, +4.328468e-03}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +0.000000e+00, /* Stl_rec_mlt */
  {+3.392652e-01, +1.483315e-03}, /* Sasc[2] */
  +3.846410e-04, +1.039467e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +2.780684e-01, /* Smn_rec */
  +0.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +2.780684e-01, /* Smn_asc */
  +0.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+4.145048e-02, +0.000000e+00}, {+7.314889e-02, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {+9.603540e-03, +1.671480e-02}, /* fp_np_sgl[2] */
  {{+4.720620e-02, +5.445670e-02}, {+6.584080e-02, +8.573460e-02}}, /* tau_np_p_sgl[2][2] */
  {{+6.505540e-02, +8.919230e-02}, {+9.677800e-02, +1.538420e-01}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+6.189260e-02, +8.952890e-02}, {+8.835669e-02, +1.417249e-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {+5.488520e-03, +6.839920e-03}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{+3.946010e-03, +9.076430e-03}, {+4.917510e-03, +1.035090e-02}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{+7.729440e-03, +1.610270e-02}, {+8.781900e-03, +2.003600e-02}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+4.104007e-03, +1.390391e-03}, /* tau_k_1_p[2] */
  {+2.176433e-03, +1.106957e-03}, /* tau_k_2_p[2] */
  {+4.061705e-03, +1.496087e-03}, /* tau_k_1_n[2] */
  {+2.249873e-03, +1.242811e-03}, /* tau_k_2_n[2] */
  {+1.823305e-03, +9.961778e-04}, /* sigma_k[2] */
  +2.676453e-02, /* fk2 */
  /* outlier parameters */
  +6.795858e+00, /* sigma_ol */
  +7.092991e-03, /* f_ol_sgl */
  +1.269280e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_nega_error_svd1_ichep06_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {-7.985202e-02, -3.893659e-03}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +0.000000e+00, /* Stl_rec_mlt */
  {-3.982427e-01, -4.158069e-03}, /* Sasc[2] */
  -3.717280e-04, -5.251425e-02, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  -6.925574e-02, /* Smn_rec */
  +0.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  -6.925574e-02, /* Smn_asc */
  +0.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{-1.412783e-01, +0.000000e+00}, {-8.364410e-02, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {-9.845610e-03, -1.758140e-02}, /* fp_np_sgl[2] */
  {{-4.655180e-02, -5.277370e-02}, {-6.301150e-02, -8.340570e-02}}, /* tau_np_p_sgl[2][2] */
  {{-6.044220e-02, -8.924750e-02}, {-8.530400e-02, -1.551850e-01}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{-6.689689e-02, -7.732681e-02}, {-8.485551e-02, -1.552965e-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-5.364410e-03, -6.722400e-03}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-3.898330e-03, -9.086010e-03}, {-4.883970e-03, -1.032830e-02}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-7.576050e-03, -1.541360e-02}, {-8.779740e-03, -1.921560e-02}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {-3.856776e-03, -1.453927e-03}, /* tau_k_1_p[2] */
  {-2.213212e-03, -1.091627e-03}, /* tau_k_2_p[2] */
  {-3.824060e-03, -1.574079e-03}, /* tau_k_1_n[2] */
  {-2.287365e-03, -1.223838e-03}, /* tau_k_2_n[2] */
  {-1.784327e-03, -1.008945e-03}, /* sigma_k[2] */
  -2.686850e-02, /* fk2 */
  /* outlier parameters */
  -9.011589e+00, /* sigma_ol */
  -7.611657e-03, /* f_ol_sgl */
  -8.150846e-05, /* f_ol_mul */
};


dtres_param_t dtres_param_center_svd2_ichep06_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+6.557680e-01, +3.504870e-02}, /* Srec[2] */
  {+1.014940e-01, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +4.752860e+00, /* Stl_rec_mlt */
  {+7.901270e-01, +1.382910e-02}, /* Sasc[2] */
  +3.356030e-02, +8.772090e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.022370e+00, /* Smn_rec */
  +3.595670e+00, /* Stl_rec */
  +1.088300e-01, /* ftl_rec */
  +1.022370e+00, /* Smn_asc */
  +3.595670e+00, /* Stl_asc */
  +1.088300e-01, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.055090e-01, +1.000000e+00}, {+7.077830e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.409870e-01, +8.340980e-01}, /* fp_np_sgl[2] */
  {{+1.450340e+00, +0.000000e+00}, {+8.861100e-01, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+1.001860e+00, +0.000000e+00}, {+5.190220e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+2.967970e-01, +3.846020e-01}, {+2.827400e-01, +2.543190e-01}}, /* fd_np_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_st[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fd_np_mlt_stxi[2] */
  {+7.702640e-01, +7.619030e-01}, /* fp_np_mlt[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* fn_np_mlt[2] */
  {{+1.106250e-01, +6.400080e-01}, {+1.005250e-01, +5.708080e-01}}, /* tau_np_p_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_p_mlt_stxi[2] */
  {{+1.285160e-01, +5.893880e-01}, {+9.230920e-02, +5.605140e-01}}, /* tau_np_n_mlt[2][2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_xi[2] */
  {DTRES_EXTERAM_VAL, DTRES_EXTERAM_VAL}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+1.658325e-01, +1.200647e-02}, /* tau_k_1_p[2] */
  {+6.344224e-02, +3.424377e-02}, /* tau_k_2_p[2] */
  {+1.675831e-01, +1.193076e-02}, /* tau_k_1_n[2] */
  {+6.454622e-02, +3.200147e-02}, /* tau_k_2_n[2] */
  {-1.711954e-02, +5.607785e-02}, /* sigma_k[2] */
  +5.123095e-01, /* fk2 */
  /* outlier parameters */
  +3.362320e+01, /* sigma_ol */
  +2.618360e-02, /* f_ol_sgl */
  +2.417640e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_posi_error_svd2_ichep06_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+4.525464e-01, +8.494203e-03}, /* Srec[2] */
  {+3.984694e-02, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +4.281897e-01, /* Stl_rec_mlt */
  {+4.943813e-01, +2.640080e-03}, /* Sasc[2] */
  +3.796450e-04, +2.007288e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +4.450155e-01, /* Smn_rec */
  +4.462291e+00, /* Stl_rec */
  +4.205170e-02, /* ftl_rec */
  +4.450155e-01, /* Smn_asc */
  +4.462291e+00, /* Stl_asc */
  +4.205170e-02, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+2.936571e-02, +0.000000e+00}, {+5.574668e-02, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {+1.150870e-02, +2.599210e-02}, /* fp_np_sgl[2] */
  {{+3.024160e-02, +0.000000e+00}, {+3.346480e-02, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+7.423130e-02, +0.000000e+00}, {+8.124480e-02, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+8.849158e-02, +8.238034e-02}, {+1.201289e-01, +1.431401e-01}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {+2.405770e-03, +2.947600e-03}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{+2.076280e-03, +6.172150e-03}, {+2.263810e-03, +6.357770e-03}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{+3.812610e-03, +8.710030e-03}, {+3.835400e-03, +9.824210e-03}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.694618e-03, +1.033360e-03}, /* tau_k_1_p[2] */
  {+1.799392e-03, +9.972194e-04}, /* tau_k_2_p[2] */
  {+2.773887e-03, +1.085271e-03}, /* tau_k_1_n[2] */
  {+1.768092e-03, +1.077810e-03}, /* tau_k_2_n[2] */
  {+1.375563e-03, +6.966386e-04}, /* sigma_k[2] */
  +1.958643e-02, /* fk2 */
  /* outlier parameters */
  +6.801315e+00, /* sigma_ol */
  +5.266311e-03, /* f_ol_sgl */
  +8.545992e-05, /* f_ol_mul */
};

dtres_param_t dtres_param_nega_error_svd2_ichep06_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {-9.792007e-02, -6.010608e-03}, /* Srec[2] */
  {-5.254031e-02, +0.000000e+00}, /* ftl_rec_mlt[2] */
  -1.028490e+00, /* Stl_rec_mlt */
  {-2.602271e-01, -6.141551e-03}, /* Sasc[2] */
  -3.708640e-04, -7.748610e-02, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  -7.903206e-02, /* Smn_rec */
  -8.034759e-01, /* Stl_rec */
  -4.464682e-02, /* ftl_rec */
  -7.903206e-02, /* Smn_asc */
  -8.034759e-01, /* Stl_asc */
  -4.464682e-02, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{-8.212332e-02, +0.000000e+00}, {-7.239110e-02, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {-1.191620e-02, -2.769920e-02}, /* fp_np_sgl[2] */
  {{-2.988920e-02, +0.000000e+00}, {-3.253030e-02, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{-7.026250e-02, +0.000000e+00}, {-7.253760e-02, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{-3.819023e-02, -5.236396e-02}, {-8.159194e-02, -7.999538e-02}}, /* fd_np_mlt[2][2] */
  {0., 0.}, /* fd_np_mlt_st[2] */
  {0., 0.}, /* fd_np_mlt_xi[2] */
  {0., 0.}, /* fd_np_mlt_stxi[2] */
  {-2.401830e-03, -2.935040e-03}, /* fp_np_mlt[2] */
  {0., 0.}, /* fn_np_mlt[2] */
  {{-2.063710e-03, -6.192170e-03}, {-2.250420e-03, -6.357270e-03}}, /* tau_np_p_mlt[2][2] */
  {0., 0.}, /* tau_np_p_mlt_xi[2] */
  {0., 0.}, /* tau_np_p_mlt_stxi[2] */
  {{-3.767710e-03, -8.628310e-03}, {-3.797410e-03, -9.692480e-03}}, /* tau_np_n_mlt[2][2] */
  {0., 0.}, /* tau_np_n_mlt_xi[2] */
  {0., 0.}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {-2.574345e-03, -1.069034e-03}, /* tau_k_1_p[2] */
  {-1.815618e-03, -9.913188e-04}, /* tau_k_2_p[2] */
  {-2.649232e-03, -1.128920e-03}, /* tau_k_1_n[2] */
  {-1.779127e-03, -1.073200e-03}, /* tau_k_2_n[2] */
  {-1.376546e-03, -7.009413e-04}, /* sigma_k[2] */
  -1.940598e-02, /* fk2 */
  /* outlier parameters */
  -8.970727e+00, /* sigma_ol */
  -4.412399e-03, /* f_ol_sgl */
  -6.665218e-05, /* f_ol_mul */
};

dtres_param_t dtres_param_center_svd1_2010mdlh_mc = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+9.626350e-01, +1.985560e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+7.290850e-01, +1.719270e-01}, /* Sasc[2] */
  +0.000000e+00, +1.000000e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.109750e+00, /* Smn_rec */
  +1.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +1.109750e+00, /* Smn_asc */
  +1.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.817180e-01, +1.000000e+00}, {+8.295620e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.186460e-01, +8.495510e-01}, /* fp_np_sgl[2] */
  {{+1.847670e+00, +0.000000e+00}, {+1.417050e+00, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+2.041140e+00, +0.000000e+00}, {+1.982780e+00, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+4.664410e-01, +6.371510e-01}, {+4.053710e-01, +6.110950e-01}}, /* fd_np_mlt[2][2] */
  {+2.706290e-01, +3.052200e-01}, /* fd_np_mlt_st[2] */
  {-2.204070e-01, -1.830920e-01}, /* fd_np_mlt_xi[2] */
  {+2.228050e-01, +1.983780e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.232800e-01, +1.371750e-01}, /* fn_np_mlt[2] */
  {{-5.202290e-03, +7.168080e-01}, {-6.925350e-03, +6.471710e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.966400e-02, -2.786130e-02}, /* tau_np_p_mlt_xi[2] */
  {+2.514870e-01, +2.420110e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+4.517970e-02, +5.151980e-01}, {-3.982540e-02, +6.046470e-01}}, /* tau_np_n_mlt[2][2] */
  {-7.831950e-02, -3.812870e-02}, /* tau_np_n_mlt_xi[2] */
  {+4.304680e-01, +3.372510e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.351098e-01, /* fk2 */
  /* outlier parameters */
  +3.319470e+01, /* sigma_ol */
  +3.784100e-02, /* f_ol_sgl */
  +2.146080e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_center_svd1_2010mdlh_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+7.046620e-01, +2.120840e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+4.834940e-01, +2.366200e-01}, /* Sasc[2] */
  +0.000000e+00, +1.057080e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +9.798010e-01, /* Smn_rec */
  +1.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +9.798010e-01, /* Smn_asc */
  +1.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.817180e-01, +1.000000e+00}, {+8.295620e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.186460e-01, +8.495510e-01}, /* fp_np_sgl[2] */
  {{+1.847670e+00, +0.000000e+00}, {+1.417050e+00, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+2.041140e+00, +0.000000e+00}, {+1.982780e+00, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+4.664410e-01, +6.371510e-01}, {+4.053710e-01, +6.110950e-01}}, /* fd_np_mlt[2][2] */
  {+2.706290e-01, +3.052200e-01}, /* fd_np_mlt_st[2] */
  {-2.204070e-01, -1.830920e-01}, /* fd_np_mlt_xi[2] */
  {+2.228050e-01, +1.983780e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.232800e-01, +1.371750e-01}, /* fn_np_mlt[2] */
  {{-5.202290e-03, +7.168080e-01}, {-6.925350e-03, +6.471710e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.966400e-02, -2.786130e-02}, /* tau_np_p_mlt_xi[2] */
  {+2.514870e-01, +2.420110e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+4.517970e-02, +5.151980e-01}, {-3.982540e-02, +6.046470e-01}}, /* tau_np_n_mlt[2][2] */
  {-7.831950e-02, -3.812870e-02}, /* tau_np_n_mlt_xi[2] */
  {+4.304680e-01, +3.372510e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.351098e-01, /* fk2 */
  /* outlier parameters */
  +4.369930e+01, /* sigma_ol */
  +3.700380e-02, /* f_ol_sgl */
  +1.141840e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_posi_error_svd1_2010mdlh_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+1.714002e-01, +4.347115e-02}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +0.000000e+00, /* Stl_rec_mlt */
  {+2.820590e-01, +3.743452e-02}, /* Sasc[2] */
  +0.000000e+00, +1.168773e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +2.777485e-01, /* Smn_rec */
  +0.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +2.777485e-01, /* Smn_asc */
  +0.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+1.033570e-02, +0.000000e+00}, {+9.976570e-03, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {+1.728010e-02, +1.775870e-02}, /* fp_np_sgl[2] */
  {{+8.514250e-02, +0.000000e+00}, {+8.078460e-02, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+2.773060e-01, +0.000000e+00}, {+2.876040e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+4.014220e-02, +3.985010e-02}, {+3.832010e-02, +3.764570e-02}}, /* fd_np_mlt[2][2] */
  {+7.042270e-02, +6.576100e-02}, /* fd_np_mlt_st[2] */
  {+1.581500e-02, +1.578260e-02}, /* fd_np_mlt_xi[2] */
  {+2.917370e-02, +2.775710e-02}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+7.256960e-03, +8.579120e-03}, /* fn_np_mlt[2] */
  {{+1.544320e-02, +3.880040e-02}, {+1.252400e-02, +3.221560e-02}}, /* tau_np_p_mlt[2][2] */
  {+5.370360e-03, +4.757190e-03}, /* tau_np_p_mlt_xi[2] */
  {+1.442400e-02, +1.300540e-02}, /* tau_np_p_mlt_stxi[2] */
  {{+2.696420e-02, +6.384480e-02}, {+1.982320e-02, +5.045020e-02}}, /* tau_np_n_mlt[2][2] */
  {+1.321680e-02, +1.033470e-02}, /* tau_np_n_mlt_xi[2] */
  {+3.442280e-02, +2.706580e-02}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+4.641307e-03, +0.000000e+00}, /* tau_k_1_p[2] */
  {+6.089655e-03, +0.000000e+00}, /* tau_k_2_p[2] */
  {+4.641307e-03, +0.000000e+00}, /* tau_k_1_n[2] */
  {+6.089655e-03, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +4.594755e-02, /* fk2 */
  /* outlier parameters */
  +1.496606e+01, /* sigma_ol */
  +6.532360e-03, /* f_ol_sgl */
  +5.683895e-05, /* f_ol_mul */
};

dtres_param_t dtres_param_nega_error_svd1_2010mdlh_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {-1.065471e-01, -4.704013e-02}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +0.000000e+00, /* Stl_rec_mlt */
  {-7.318439e-02, -5.605924e-02}, /* Sasc[2] */
  +0.000000e+00, -1.630550e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  -3.636229e-02, /* Smn_rec */
  +0.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  -3.636229e-02, /* Smn_asc */
  +0.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{-1.081850e-02, +0.000000e+00}, {-1.049090e-02, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {-1.847270e-02, -1.921890e-02}, /* fp_np_sgl[2] */
  {{-8.060650e-02, +0.000000e+00}, {-7.559220e-02, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{-2.431610e-01, +0.000000e+00}, {-2.570580e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{-4.109540e-02, -4.082290e-02}, {-3.895900e-02, -3.819170e-02}}, /* fd_np_mlt[2][2] */
  {-7.026490e-02, -6.626680e-02}, /* fd_np_mlt_st[2] */
  {-1.527740e-02, -1.453960e-02}, /* fd_np_mlt_xi[2] */
  {-3.001790e-02, -2.893270e-02}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {-7.267580e-03, -9.143810e-03}, /* fn_np_mlt[2] */
  {{-1.535860e-02, -3.861560e-02}, {-1.247460e-02, -3.175120e-02}}, /* tau_np_p_mlt[2][2] */
  {-5.304540e-03, -4.782510e-03}, /* tau_np_p_mlt_xi[2] */
  {-1.437050e-02, -1.281480e-02}, /* tau_np_p_mlt_stxi[2] */
  {{-2.683960e-02, -6.260460e-02}, {-1.994650e-02, -4.913800e-02}}, /* tau_np_n_mlt[2][2] */
  {-1.309970e-02, -1.037550e-02}, /* tau_np_n_mlt_xi[2] */
  {-3.318160e-02, -2.595050e-02}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {-4.091320e-03, +0.000000e+00}, /* tau_k_1_p[2] */
  {-6.138504e-03, +0.000000e+00}, /* tau_k_2_p[2] */
  {-4.091320e-03, +0.000000e+00}, /* tau_k_1_n[2] */
  {-6.138504e-03, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  -4.137452e-02, /* fk2 */
  /* outlier parameters */
  -4.845010e+00, /* sigma_ol */
  -1.079975e-02, /* f_ol_sgl */
  -6.836431e-05, /* f_ol_mul */
};


dtres_param_t dtres_param_center_svd2_2010mdlh_mc = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+9.271430e-01, +2.103700e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+8.210910e-01, +1.408000e-01}, /* Sasc[2] */
  +0.000000e+00, +1.000000e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.053040e+00, /* Smn_rec */
  +4.320600e+00, /* Stl_rec */
  +7.068970e-02, /* ftl_rec */
  +1.053040e+00, /* Smn_asc */
  +4.320600e+00, /* Stl_asc */
  +7.068970e-02, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.737830e-01, +1.000000e+00}, {+8.062880e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.013430e-01, +8.274010e-01}, /* fp_np_sgl[2] */
  {{+1.625960e+00, +0.000000e+00}, {+9.860720e-01, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+9.181350e-01, +0.000000e+00}, {+4.316230e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+5.600550e-01, +7.507230e-01}, {+5.345170e-01, +7.087410e-01}}, /* fd_np_mlt[2][2] */
  {+1.569090e-01, +1.870940e-01}, /* fd_np_mlt_st[2] */
  {-2.021350e-01, -1.879030e-01}, /* fd_np_mlt_xi[2] */
  {+2.324180e-01, +2.255560e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.224260e-01, +1.210540e-01}, /* fn_np_mlt[2] */
  {{+3.871670e-02, +7.653070e-01}, {+2.442560e-02, +7.411070e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.632670e-02, -1.976910e-02}, /* tau_np_p_mlt_xi[2] */
  {+3.214540e-01, +2.753010e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+8.292430e-02, +5.342620e-01}, {+6.553030e-02, +5.239120e-01}}, /* tau_np_n_mlt[2][2] */
  {-3.013640e-02, -1.906810e-02}, /* tau_np_n_mlt_xi[2] */
  {+3.899160e-01, +3.309140e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.107397e-01, /* fk2 */
  /* outlier parameters */
  +3.063190e+01, /* sigma_ol */
  +2.231740e-02, /* f_ol_sgl */
  +8.843970e-05, /* f_ol_mul */
};

dtres_param_t dtres_param_center_svd2_2010mdlh_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+8.075310e-01, +2.326080e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+6.429470e-01, +2.290660e-01}, /* Sasc[2] */
  +0.000000e+00, +1.014100e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.014710e+00, /* Smn_rec */
  +3.662940e+00, /* Stl_rec */
  +1.110420e-01, /* ftl_rec */
  +1.014710e+00, /* Smn_asc */
  +3.662940e+00, /* Stl_asc */
  +1.110420e-01, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.737830e-01, +1.000000e+00}, {+8.062880e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.013430e-01, +8.274010e-01}, /* fp_np_sgl[2] */
  {{+1.625960e+00, +0.000000e+00}, {+9.860720e-01, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+9.181350e-01, +0.000000e+00}, {+4.316230e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+5.600550e-01, +7.507230e-01}, {+5.345170e-01, +7.087410e-01}}, /* fd_np_mlt[2][2] */
  {+1.569090e-01, +1.870940e-01}, /* fd_np_mlt_st[2] */
  {-2.021350e-01, -1.879030e-01}, /* fd_np_mlt_xi[2] */
  {+2.324180e-01, +2.255560e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.224260e-01, +1.210540e-01}, /* fn_np_mlt[2] */
  {{+3.871670e-02, +7.653070e-01}, {+2.442560e-02, +7.411070e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.632670e-02, -1.976910e-02}, /* tau_np_p_mlt_xi[2] */
  {+3.214540e-01, +2.753010e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+8.292430e-02, +5.342620e-01}, {+6.553030e-02, +5.239120e-01}}, /* tau_np_n_mlt[2][2] */
  {-3.013640e-02, -1.906810e-02}, /* tau_np_n_mlt_xi[2] */
  {+3.899160e-01, +3.309140e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.107397e-01, /* fk2 */
  /* outlier parameters */
  +3.352530e+01, /* sigma_ol */
  +2.730960e-02, /* f_ol_sgl */
  +1.529520e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_posi_error_svd2_2010mdlh_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+2.793806e-01, +6.850826e-02}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +0.000000e+00, /* Stl_rec_mlt */
  {+3.862413e-01, +2.650982e-02}, /* Sasc[2] */
  +0.000000e+00, +7.507394e-02, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +4.407071e-01, /* Smn_rec */
  +3.616785e+00, /* Stl_rec */
  +1.547901e-02, /* ftl_rec */
  +4.407071e-01, /* Smn_asc */
  +3.616785e+00, /* Stl_asc */
  +1.547901e-02, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+8.226390e-03, +0.000000e+00}, {+1.278640e-02, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {+1.720800e-02, +3.618370e-02}, /* fp_np_sgl[2] */
  {{+4.998000e-02, +0.000000e+00}, {+4.034130e-02, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+8.457830e-02, +0.000000e+00}, {+8.966180e-02, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+1.285450e-02, +1.289750e-02}, {+1.200810e-02, +1.201540e-02}}, /* fd_np_mlt[2][2] */
  {+2.787120e-02, +2.669110e-02}, /* fd_np_mlt_st[2] */
  {+4.715760e-03, +4.766360e-03}, /* fd_np_mlt_xi[2] */
  {+1.194060e-02, +1.194550e-02}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+2.582800e-03, +2.907050e-03}, /* fn_np_mlt[2] */
  {{+5.570230e-03, +2.019060e-02}, {+4.425930e-03, +1.680450e-02}}, /* tau_np_p_mlt[2][2] */
  {+1.706220e-03, +1.529650e-03}, /* tau_np_p_mlt_xi[2] */
  {+7.185460e-03, +6.540220e-03}, /* tau_np_p_mlt_stxi[2] */
  {{+9.327340e-03, +3.109430e-02}, {+7.401770e-03, +2.599310e-02}}, /* tau_np_n_mlt[2][2] */
  {+4.305020e-03, +3.593850e-03}, /* tau_np_n_mlt_xi[2] */
  {+1.588790e-02, +1.355840e-02}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+1.548688e-03, +0.000000e+00}, /* tau_k_1_p[2] */
  {+2.349945e-03, +0.000000e+00}, /* tau_k_2_p[2] */
  {+1.548688e-03, +0.000000e+00}, /* tau_k_1_n[2] */
  {+2.349945e-03, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +1.653093e-02, /* fk2 */
  /* outlier parameters */
  +5.583821e+00, /* sigma_ol */
  +2.035564e-03, /* f_ol_sgl */
  +1.022429e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_nega_error_svd2_2010mdlh_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {-1.478871e-01, -5.859953e-02}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +0.000000e+00, /* Stl_rec_mlt */
  {-7.211574e-02, -5.239444e-02}, /* Sasc[2] */
  +0.000000e+00, -1.643656e-01, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +0.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  -3.897418e-02, /* Smn_rec */
  -3.911839e-01, /* Stl_rec */
  -3.946462e-02, /* ftl_rec */
  -3.897418e-02, /* Smn_asc */
  -3.911839e-01, /* Stl_asc */
  -3.946462e-02, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{-8.599940e-03, +0.000000e+00}, {-1.431760e-02, +0.000000e+00}}, /* fd_np_sgl[2][2] */
  {-1.769460e-02, -3.903760e-02}, /* fp_np_sgl[2] */
  {{-4.888840e-02, +0.000000e+00}, {-3.912560e-02, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{-7.819400e-02, +0.000000e+00}, {-7.529680e-02, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{-1.248950e-02, -1.255940e-02}, {-1.187880e-02, -1.189400e-02}}, /* fd_np_mlt[2][2] */
  {-2.910630e-02, -2.746740e-02}, /* fd_np_mlt_st[2] */
  {-4.845690e-03, -4.824440e-03}, /* fd_np_mlt_xi[2] */
  {-1.150930e-02, -1.166030e-02}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {-2.553030e-03, -2.867460e-03}, /* fn_np_mlt[2] */
  {{-5.430320e-03, -2.052420e-02}, {-4.336660e-03, -1.704030e-02}}, /* tau_np_p_mlt[2][2] */
  {-1.737370e-03, -1.545770e-03}, /* tau_np_p_mlt_xi[2] */
  {-7.070600e-03, -6.443980e-03}, /* tau_np_p_mlt_stxi[2] */
  {{-9.286360e-03, -3.092940e-02}, {-7.401010e-03, -2.571140e-02}}, /* tau_np_n_mlt[2][2] */
  {-4.304750e-03, -3.568940e-03}, /* tau_np_n_mlt_xi[2] */
  {-1.567680e-02, -1.343530e-02}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {-1.473479e-03, +0.000000e+00}, /* tau_k_1_p[2] */
  {-2.347952e-03, +0.000000e+00}, /* tau_k_2_p[2] */
  {-1.473479e-03, +0.000000e+00}, /* tau_k_1_n[2] */
  {-2.347952e-03, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  -1.577599e-02, /* fk2 */
  /* outlier parameters */
  -9.217143e+00, /* sigma_ol */
  -4.659891e-03, /* f_ol_sgl */
  -7.195453e-05, /* f_ol_mul */
};


dtres_param_t dtres_param_center_svd1_2010nn_mc = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+9.670470e-01, +2.006890e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+7.219140e-01, +1.733930e-01}, /* Sasc[2] */
  +0.000000e+00, +1.000000e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.109220e+00, /* Smn_rec */
  +1.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +1.109220e+00, /* Smn_asc */
  +1.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.821120e-01, +1.000000e+00}, {+8.300350e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.170990e-01, +8.500410e-01}, /* fp_np_sgl[2] */
  {{+1.846290e+00, +0.000000e+00}, {+1.417820e+00, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+2.064200e+00, +0.000000e+00}, {+1.988790e+00, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+4.671060e-01, +6.398390e-01}, {+4.028770e-01, +6.087960e-01}}, /* fd_np_mlt[2][2] */
  {+2.600710e-01, +3.061690e-01}, /* fd_np_mlt_st[2] */
  {-2.223320e-01, -1.825950e-01}, /* fd_np_mlt_xi[2] */
  {+2.300750e-01, +1.975340e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.236850e-01, +1.383780e-01}, /* fn_np_mlt[2] */
  {{-3.897030e-03, +7.127320e-01}, {-6.807620e-03, +6.468940e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.948990e-02, -2.773490e-02}, /* tau_np_p_mlt_xi[2] */
  {+2.514000e-01, +2.417630e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+4.662370e-02, +5.134400e-01}, {-3.966000e-02, +6.041600e-01}}, /* tau_np_n_mlt[2][2] */
  {-7.839280e-02, -3.795670e-02}, /* tau_np_n_mlt_xi[2] */
  {+4.332160e-01, +3.357360e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.351098e-01, /* fk2 */
  /* outlier parameters */
  +3.295330e+01, /* sigma_ol */
  +3.770840e-02, /* f_ol_sgl */
  +2.086000e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_center_svd1_2010nn_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+7.015380e-01, +2.113580e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+4.772600e-01, +2.362870e-01}, /* Sasc[2] */
  +0.000000e+00, +1.057680e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +9.769520e-01, /* Smn_rec */
  +1.000000e+00, /* Stl_rec */
  +0.000000e+00, /* ftl_rec */
  +9.769520e-01, /* Smn_asc */
  +1.000000e+00, /* Stl_asc */
  +0.000000e+00, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.821120e-01, +1.000000e+00}, {+8.300350e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.170990e-01, +8.500410e-01}, /* fp_np_sgl[2] */
  {{+1.846290e+00, +0.000000e+00}, {+1.417820e+00, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+2.064200e+00, +0.000000e+00}, {+1.988790e+00, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+4.671060e-01, +6.398390e-01}, {+4.028770e-01, +6.087960e-01}}, /* fd_np_mlt[2][2] */
  {+2.600710e-01, +3.061690e-01}, /* fd_np_mlt_st[2] */
  {-2.223320e-01, -1.825950e-01}, /* fd_np_mlt_xi[2] */
  {+2.300750e-01, +1.975340e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.236850e-01, +1.383780e-01}, /* fn_np_mlt[2] */
  {{-3.897030e-03, +7.127320e-01}, {-6.807620e-03, +6.468940e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.948990e-02, -2.773490e-02}, /* tau_np_p_mlt_xi[2] */
  {+2.514000e-01, +2.417630e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+4.662370e-02, +5.134400e-01}, {-3.966000e-02, +6.041600e-01}}, /* tau_np_n_mlt[2][2] */
  {-7.839280e-02, -3.795670e-02}, /* tau_np_n_mlt_xi[2] */
  {+4.332160e-01, +3.357360e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.086473e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.934304e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.351098e-01, /* fk2 */
  /* outlier parameters */
  +4.368020e+01, /* sigma_ol */
  +3.700480e-02, /* f_ol_sgl */
  +1.166780e-04, /* f_ol_mul */
};

dtres_param_t dtres_param_center_svd2_2010nn_mc = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+9.485760e-01, +2.074960e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+7.881960e-01, +1.170810e-01}, /* Sasc[2] */
  +0.000000e+00, +1.000000e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.047350e+00, /* Smn_rec */
  +4.394850e+00, /* Stl_rec */
  +6.804490e-02, /* ftl_rec */
  +1.047350e+00, /* Smn_asc */
  +4.394850e+00, /* Stl_asc */
  +6.804490e-02, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.714090e-01, +1.000000e+00}, {+8.055240e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.001330e-01, +8.281430e-01}, /* fp_np_sgl[2] */
  {{+1.619410e+00, +0.000000e+00}, {+9.816200e-01, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+9.166640e-01, +0.000000e+00}, {+4.314210e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+5.587240e-01, +7.479840e-01}, {+5.360980e-01, +7.103110e-01}}, /* fd_np_mlt[2][2] */
  {+1.631980e-01, +1.823800e-01}, /* fd_np_mlt_st[2] */
  {-2.016900e-01, -1.888730e-01}, /* fd_np_mlt_xi[2] */
  {+2.301290e-01, +2.294210e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.227150e-01, +1.208840e-01}, /* fn_np_mlt[2] */
  {{+3.949700e-02, +7.642590e-01}, {+2.628330e-02, +7.346910e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.630320e-02, -2.051890e-02}, /* tau_np_p_mlt_xi[2] */
  {+3.212710e-01, +2.781350e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+8.286050e-02, +5.311100e-01}, {+6.727790e-02, +5.206650e-01}}, /* tau_np_n_mlt[2][2] */
  {-2.979580e-02, -1.913500e-02}, /* tau_np_n_mlt_xi[2] */
  {+3.892550e-01, +3.308110e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.107397e-01, /* fk2 */
  /* outlier parameters */
  +3.076940e+01, /* sigma_ol */
  +2.225020e-02, /* f_ol_sgl */
  +7.313210e-05, /* f_ol_mul */
};

dtres_param_t dtres_param_center_svd2_2010nn_data = {
  /*** Rdet ***/
  /* multiple track vertex */
  {+7.990420e-01, +2.325800e-01}, /* Srec[2] */
  {+0.000000e+00, +0.000000e+00}, /* ftl_rec_mlt[2] */
  +1.000000e+00, /* Stl_rec_mlt */
  {+6.463950e-01, +2.287110e-01}, /* Sasc[2] */
  +0.000000e+00, +1.013030e+00, /* Snp, Snp_global */
  {+0.000000e+00, +0.000000e+00}, /* ftl_asc_mlt[2] */
  +1.000000e+00, /* Stl_asc_mlt */
  /* single track vertex */
  +1.011270e+00, /* Smn_rec */
  +3.633810e+00, /* Stl_rec */
  +1.125900e-01, /* ftl_rec */
  +1.011270e+00, /* Smn_asc */
  +3.633810e+00, /* Stl_asc */
  +1.125900e-01, /* ftl_asc */
  /*** Rnp ***/
  /* index ... 0 <--> B0_B0B */
  /*           1 <--> B+_B-  */
  /* single track vertex */
  {{+7.714090e-01, +1.000000e+00}, {+8.055240e-01, +1.000000e+00}}, /* fd_np_sgl[2][2] */
  {+8.001330e-01, +8.281430e-01}, /* fp_np_sgl[2] */
  {{+1.619410e+00, +0.000000e+00}, {+9.816200e-01, +0.000000e+00}}, /* tau_np_p_sgl[2][2] */
  {{+9.166640e-01, +0.000000e+00}, {+4.314210e-01, +0.000000e+00}}, /* tau_np_n_sgl[2][2] */
  /* multiple track vertex */
  {{+5.587240e-01, +7.479840e-01}, {+5.360980e-01, +7.103110e-01}}, /* fd_np_mlt[2][2] */
  {+1.631980e-01, +1.823800e-01}, /* fd_np_mlt_st[2] */
  {-2.016900e-01, -1.888730e-01}, /* fd_np_mlt_xi[2] */
  {+2.301290e-01, +2.294210e-01}, /* fd_np_mlt_stxi[2] */
  {+0.000000e+00, +0.000000e+00}, /* fp_np_mlt[2] */
  {+1.227150e-01, +1.208840e-01}, /* fn_np_mlt[2] */
  {{+3.949700e-02, +7.642590e-01}, {+2.628330e-02, +7.346910e-01}}, /* tau_np_p_mlt[2][2] */
  {-2.630320e-02, -2.051890e-02}, /* tau_np_p_mlt_xi[2] */
  {+3.212710e-01, +2.781350e-01}, /* tau_np_p_mlt_stxi[2] */
  {{+8.286050e-02, +5.311100e-01}, {+6.727790e-02, +5.206650e-01}}, /* tau_np_n_mlt[2][2] */
  {-2.979580e-02, -1.913500e-02}, /* tau_np_n_mlt_xi[2] */
  {+3.892550e-01, +3.308110e-01}, /* tau_np_n_mlt_stxi[2] */
  /* Rk parameters (for partial reconstracted event) */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_p[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_p[2] */
  {+2.035482e-01, +0.000000e+00}, /* tau_k_1_n[2] */
  {+9.808753e-02, +0.000000e+00}, /* tau_k_2_n[2] */
  {+0.000000e+00, +0.000000e+00}, /* sigma_k[2] */
  +3.107397e-01, /* fk2 */
  /* outlier parameters */
  +3.351010e+01, /* sigma_ol */
  +2.732950e-02, /* f_ol_sgl */
  +1.536610e-04, /* f_ol_mul */
};

/* Default Resolution Parameters */
/* (== parameter set for 2004 publication ) */

dtres_param_t dtres_param_default    = dtres_param_center_svd1_2010mdlh_data;
dtres_param_t dtres_param_default_mc = dtres_param_center_svd1_2010mdlh_mc;
dtres_param_t dtres_param_posi_error = dtres_param_posi_error_svd1_2010mdlh_data;
dtres_param_t dtres_param_nega_error = dtres_param_nega_error_svd1_2010mdlh_data;

dtres_param_t dtres_param_default_svd2    = dtres_param_center_svd2_2010mdlh_data;
dtres_param_t dtres_param_default_svd2_mc = dtres_param_center_svd2_2010mdlh_mc;
dtres_param_t dtres_param_posi_error_svd2 = dtres_param_posi_error_svd2_2010mdlh_data;
dtres_param_t dtres_param_nega_error_svd2 = dtres_param_nega_error_svd2_2010mdlh_data;

void dt_resol_global_set_bgc(const double bgc){
  dt_resol_global::inv_bgc = 1.0/bgc;
}

void dt_resol_global_set_bmass(const double x){
  dt_resol_global::inv_mb = 1.0/x;
}
void dt_resol_global_set_b0_mass(const double x){
  dt_resol_global::mbzero     = x;
}

void dt_resol_global_set_bplus_mass(const double x){
  dt_resol_global::mbplus     = x;
}


__CLOSE_LINKAGE__
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
