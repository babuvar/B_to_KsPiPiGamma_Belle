//               
//    conv_coef.h   
//       ---- Caliculation of coefficient at convolution calculation
//    $Id: conv_coef.h 9932 2006-11-12 14:26:53Z katayama $ 
//                
//    $Log$
//    Revision 1.1  2002/09/05 01:27:06  katayama
//    New package tatami from Nakadaira/Sumisawa san
//
//    Revision 1.4  2002/07/12 09:09:59  nakadair
//    Bug fix in add_MfEnp_coef(...).
//
//    Revision 1.3  2002/05/11 07:51:37  nakadair
//    Update for solaris and GCC3.
//
//    Revision 1.2  2002/05/09 06:22:00  nakadair
//    Bug fix in add_AnEn_coef(...).
//
//    Revision 1.1  2002/05/07 21:38:22  nakadair
//    Add convolution with Resolution function and fortran interface.
// 
//                

#ifndef __CONV_COEF_H__
#define __CONV_COEF_H__

#include "belle.h"
#include <cmath>
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

inline void add_EpEn_coef(double* fEp1, double* fEn2,
			  const double tau1, const double tau2,
			  const double weight = 1.0){
  const double inv_tausum = 1.0/(tau1+tau2);
  *fEp1 += weight*tau1*inv_tausum;
  *fEn2 += weight*tau2*inv_tausum;
  return;
}

inline void add_EnEp_coef(double* fEn1, double* fEp2,
			  const double tau1, const double tau2,
			  const double weight = 1.0){
  const double inv_tausum = 1.0/(tau1+tau2);
  *fEn1 += weight*tau1*inv_tausum;
  *fEp2 += weight*tau2*inv_tausum;
  return;
}

inline void add_EnEn_coef(double* fEn1, double* fEn2, double* fxEn1,
			  const double tau1, const double tau2,
			  const double weight = 1.0){
  if(tau1==tau2 /* ||!finite(1.0/(tau1-tau2)) */){
    *fxEn1 += weight; 
  }else{
    const double inv_tausub = 1.0/(tau1-tau2);
    *fEn1 +=  weight*tau1*inv_tausub;
    *fEn2 += -weight*tau2*inv_tausub;
  }
  return;
}

inline void add_EpEp_coef(double* fEp1, double* fEp2, double* fxEp1,
			  const double tau1, const double tau2,
			  const double weight = 1.0){
  if(tau1==tau2 /* ||!finite(1.0/(tau1-tau2)) */){
    *fxEp1 += weight; 
  }else{
    const double inv_tausub = 1.0/(tau1-tau2);
    *fEp1 +=  weight*tau1*inv_tausub;
    *fEp2 += -weight*tau2*inv_tausub;
  }
  return;
}

inline void add_EpEnp_coef(double* fEp1, double* fEn2, double* fEp2, double* fxEp1,
			   const double tau1, const double tau2_n, const double tau2_p,
			   const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_EpEn_coef(fEp1, fEn2,        tau1, tau2_n, weight2*tau2_n);
  add_EpEp_coef(fEp1, fEp2, fxEp1, tau1, tau2_p, weight2*tau2_p);
  return;
}

inline void add_EpEf_coef(double* fEp1, double* fEn2, double* fEp2, double* fxEp1,
			  const double tau1, const double tau2,
			  const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_EpEn_coef(fEp1, fEn2,        tau1, tau2, weight2);
  add_EpEp_coef(fEp1, fEp2, fxEp1, tau1, tau2, weight2);
  return;
}

inline void add_EnEnp_coef(double* fEn1, double* fEn2, double* fEp2, double* fxEn1,
			   const double tau1, const double tau2_n, const double tau2_p,
			   const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_EnEn_coef(fEn1, fEn2, fxEn1, tau1, tau2_n, weight2*tau2_n);
  add_EnEp_coef(fEn1, fEp2,        tau1, tau2_p, weight2*tau2_p);
  return;
}

inline void add_EnEf_coef(double* fEn1, double* fEn2, double* fEp2, double* fxEn1,
			  const double tau1, const double tau2,
			  const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_EnEn_coef(fEn1, fEn2, fxEn1, tau1, tau2, weight2);
  add_EnEp_coef(fEn1, fEp2, tau1, tau2, weight2);
  return;
}

inline void add_EnpEnp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fxEn1, double* fxEp1,
			    const double tau1_n, const double tau1_p, const double tau2_n, const double tau2_p,
			    const double weight = 1.0){
  const double norm1   = 1.0/(std::abs(tau1_n)+std::abs(tau1_p));
  const double norm2   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2   = norm1*norm2*weight;  
  const double weight2_n = weight2*tau1_n;  
  const double weight2_p = weight2*tau1_p;  
  add_EnEn_coef(fEn1, fEn2, fxEn1, tau1_n, tau2_n, weight2_n*tau2_n);
  add_EnEp_coef(fEn1, fEp2,        tau1_n, tau2_p, weight2_n*tau2_p);
  add_EpEn_coef(fEp1, fEn2,        tau1_p, tau2_n, weight2_p*tau2_n);
  add_EpEp_coef(fEp1, fEp2, fxEp1, tau1_p, tau2_p, weight2_p*tau2_p);
  return;
}

inline void add_EfEnp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fxEn1, double* fxEp1,
			   const double tau1, const double tau2_n, const double tau2_p,
			   const double weight = 1.0){
  const double norm2   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2   = 0.5*norm2*weight;  
  const double weight2_n = weight2*tau2_n;  
  const double weight2_p = weight2*tau2_p;  
  add_EnEn_coef(fEn1, fEn2, fxEn1, tau1, tau2_n, weight2_n);
  add_EnEp_coef(fEn1, fEp2,        tau1, tau2_p, weight2_p);
  add_EpEn_coef(fEp1, fEn2,        tau1, tau2_n, weight2_n);
  add_EpEp_coef(fEp1, fEp2, fxEp1, tau1, tau2_p, weight2_p);
  return;
}

inline void add_EnpEf_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fxEn1, double* fxEp1,
			   const double tau1_n, const double tau1_p, const double tau2,
			   const double weight = 1.0){
  const double norm1   = 1.0/(std::abs(tau1_n)+std::abs(tau1_p));
  const double weight2 = norm1*0.5*weight;  
  const double weight2_n = weight2*tau1_n;  
  const double weight2_p = weight2*tau1_p;  
  add_EnEn_coef(fEn1, fEn2, fxEn1, tau1_n, tau2, weight2_n);
  add_EnEp_coef(fEn1, fEp2,        tau1_n, tau2, weight2_n);
  add_EpEn_coef(fEp1, fEn2,        tau1_p, tau2, weight2_p);
  add_EpEp_coef(fEp1, fEp2, fxEp1, tau1_p, tau2, weight2_p);
  return;
}

inline void add_EfEf_coef(double* fEn1, double* fEp1,
			  double* fEn2, double* fEp2, double* fxEn1, double* fxEp1,
			  const double tau1, const double tau2, const double weight = 1.0){
  const double weight2 = 0.25*weight;  
  add_EnEn_coef(fEn1, fEn2, fxEn1, tau1, tau2, weight2);
  add_EnEp_coef(fEn1, fEp2, tau1, tau2, weight2);
  add_EpEn_coef(fEp1, fEn2, tau1, tau2, weight2);
  add_EpEp_coef(fEp1, fEp2, fxEp1, tau1, tau2, weight2);
  return;
}

inline void add_xEpEp_coef(double* fxEp1, double* fEp1, double* fEp2, double* fxxEp1,
			   const double tau1, const double tau2, const double weight = 1.0){
  if(tau1==tau2 /* ||!finite(1.0/(tau1-tau2)) */ ){
    *fxxEp1 += weight;
  }else{
    const double inv_tausub = 1.0/(tau1-tau2);
    const double inv_tausub2 = inv_tausub*inv_tausub;
    *fxEp1 +=  weight*tau1*inv_tausub;
    *fEp1  += -weight*tau1*tau2*inv_tausub2;
    *fEp2  +=  weight*tau2*tau2*inv_tausub2;
  }
  return;
}

inline void add_xEnEn_coef(double* fxEn1, double* fEn1, double* fEn2, double* fxxEn1,
			   const double tau1, const double tau2, const double weight = 1.0){
  if(tau1==tau2 /* ||!finite(1.0/(tau1-tau2)) */ ){
    *fxxEn1 += weight;
  }else{
    const double inv_tausub = 1.0/(tau1-tau2);
    const double inv_tausub2 = inv_tausub*inv_tausub;
    *fxEn1 +=  weight*tau1*inv_tausub;
    *fEn1  += -weight*tau1*tau2*inv_tausub2;
    *fEn2  +=  weight*tau2*tau2*inv_tausub2;
  }
  return;
}

inline void add_xEpEn_coef(double* fxEp1, double* fEp1, double* fEn2,
			   const double tau1, const double tau2, const double weight = 1.0){
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;
  *fxEp1 +=  weight*tau1*inv_tausum;
  *fEp1  +=  weight*tau1*tau2*inv_tausum2;
  *fEn2  +=  weight*tau2*tau2*inv_tausum2;
  return;
}

inline void add_xEnEp_coef(double* fxEn1, double* fEn1, double* fEp2,
			   const double tau1, const double tau2, const double weight = 1.0){
  const double inv_tausum = 1.0/(tau1+tau2);
  const double inv_tausum2 = inv_tausum*inv_tausum;
  *fxEn1 +=  weight*tau1*inv_tausum;
  *fEn1  +=  weight*tau1*tau2*inv_tausum2;
  *fEp2  +=  weight*tau2*tau2*inv_tausum2;
  return;
}

inline void add_xEpEnp_coef(double* fxEp1, double* fEp1, double* fEn2, double* fEp2, double* fxxEp1, 
			    const double tau1, const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1, tau2_n, weight2*tau2_n);
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1, tau2_p, weight2*tau2_p);
  return;
}

inline void add_xEnEnp_coef(double* fxEn1, double* fEn1, double* fEn2, double* fEp2, double* fxxEn1, 
			    const double tau1, const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1, tau2_n, weight2*tau2_n);
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1, tau2_p, weight2*tau2_p);
  return;
}

inline void add_xEpEf_coef(double* fxEp1, double* fEp1, double* fEn2, double* fEp2, double* fxxEp1, 
			   const double tau1, const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1, tau2, weight2);
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1, tau2, weight2);
  return;
}

inline void add_xEnEf_coef(double* fxEn1, double* fEn1, double* fEn2, double* fEp2, double* fxxEn1, 
			   const double tau1, const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1, tau2, weight2);
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1, tau2, weight2);
  return;
}

inline void add_xEnpEp_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1, double* fEp2, double* fxxEp1,
			    const double tau1_n, const double tau1_p, const double tau2, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau1_n)+std::abs(tau1_p));
  const double weight2 = norm*weight;  
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1_n, tau2, weight2*tau1_n);
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1_p, tau2, weight2*tau1_p);
  return;
}

inline void add_xEnpEn_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1, double* fEn2, double* fxxEn1,
			    const double tau1_n, const double tau1_p, const double tau2, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau1_n)+std::abs(tau1_p));
  const double weight2 = norm*weight;  
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1_n, tau2, weight2*tau1_n);
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1_p, tau2, weight2*tau1_p);
  return;
}

inline void add_xEfEp_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1, double* fEp2, double* fxxEp1,
			   const double tau1, const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1, tau2, weight2);
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1, tau2, weight2);
  return;
}

inline void add_xEfEn_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1, double* fEn2, double* fxxEn1,
			   const double tau1, const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1, tau2, weight2);
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1, tau2, weight2);
  return;
}



inline void add_xEnpEnp_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1,
			     double* fEn2, double* fEp2, double* fxxEn1, double* fxxEp1,
			     const double tau1_n, const double tau1_p,
			     const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm1   = 1.0/(std::abs(tau1_n)+std::abs(tau1_p));
  const double norm2   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm1*norm2*weight;  
  const double weight2_n = weight2*tau1_n;  
  const double weight2_p = weight2*tau1_p;  
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1_n, tau2_n, weight2_n*tau2_n);
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1_n, tau2_p, weight2_n*tau2_p);
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1_p, tau2_p, weight2_p*tau2_p);
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1_p, tau2_n, weight2_p*tau2_n);

  return;
}

inline void add_xEnpEf_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1,
			    double* fEn2, double* fEp2, double* fxxEn1, double* fxxEp1,
			    const double tau1_n, const double tau1_p,
			    const double tau2, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau1_n)+std::abs(tau1_p));
  const double weight2 = 0.5*norm*weight;  
  const double weight2_n = weight2*tau1_n;  
  const double weight2_p = weight2*tau1_p;  
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1_n, tau2, weight2_n);
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1_n, tau2, weight2_n);
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1_p, tau2, weight2_p);
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1_p, tau2, weight2_p);

  return;
}

inline void add_xEfEnp_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1,
			    double* fEn2, double* fEp2, double* fxxEn1, double* fxxEp1,
			    const double tau1, const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = 0.5*norm*weight;  
  const double weight2_n = weight2*tau2_n;  
  const double weight2_p = weight2*tau2_p;  
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1, tau2_n, weight2_n);
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1, tau2_p, weight2_p);
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1, tau2_n, weight2_n);
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1, tau2_p, weight2_p);

  return;
}

inline void add_xEfEf_coef(double* fxEp1, double* fEp1, double* fxEn1, double* fEn1,
			   double* fEn2, double* fEp2, double* fxxEn1, double* fxxEp1,
			   const double tau1, const double tau2, const double weight = 1.0){
  const double weight2 = 0.25*weight;  
  add_xEpEp_coef(fxEp1, fEp1, fEp2, fxxEp1, tau1, tau2, weight2);
  add_xEnEp_coef(fxEn1, fEn1, fEp2,         tau1, tau2, weight2);
  add_xEpEn_coef(fxEp1, fEp1, fEn2,         tau1, tau2, weight2);
  add_xEnEn_coef(fxEn1, fEn1, fEn2, fxxEn1, tau1, tau2, weight2);
  return;
}


inline void add_ApEp_coef(double* fAp1, double* fMp1, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fAp1 += -weight*A*inv_at2*inv_tau;
  *fMp1 += -weight*A*inv_at2*dm;
  *fEp2 +=  weight*A*dm;
  return;
}

inline void add_AnEn_coef(double* fAn1, double* fMn1, double* fEn2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fAn1 += -weight*A*inv_at2*inv_tau;
  *fMn1 += +weight*A*inv_at2*dm;
  *fEn2 += -weight*A*dm;
  return;
}

inline void add_ApEn_coef(double* fAp1, double* fMp1, double* fEn2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fAp1 += weight*A*inv_at2*inv_tau;
  *fMp1 += weight*A*inv_at2*dm;
  *fEn2 += weight*A*dm;
  return;
}

inline void add_AnEp_coef(double* fAn1, double* fMn1, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fAn1 +=  weight*A*inv_at2*inv_tau;
  *fMn1 += -weight*A*inv_at2*dm;
  *fEp2 += -weight*A*dm;
  return;
}

inline void add_MpEp_coef(double* fMp1, double* fAp1, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fMp1 += -weight*A*inv_at2*inv_tau;
  *fAp1 +=  weight*A*inv_at2*dm;
  *fEp2 +=  weight*A*inv_tau;
  return;
}

inline void add_MnEn_coef(double* fMn1, double* fAn1, double* fEn2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1-inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fMn1 += -weight*A*inv_at2*inv_tau;
  *fAn1 += -weight*A*inv_at2*dm;
  *fEn2 +=  weight*A*inv_tau;
  return;
}

inline void add_MpEn_coef(double* fMp1, double* fAp1, double* fEn2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fMp1 +=  weight*A*inv_at2*inv_tau;
  *fAp1 += -weight*A*inv_at2*dm;
  *fEn2 +=  weight*A*inv_tau;
  return;
}

inline void add_MnEp_coef(double* fMn1, double* fAn1, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double inv_at1 = 1.0/std::abs(tau1);
  const double inv_at2 = 1.0/std::abs(tau2);
  const double inv_tau = inv_at1+inv_at2;
  const double A = 1.0/(inv_tau*inv_tau+dm*dm);
  *fMn1 +=  weight*A*inv_at2*inv_tau;
  *fAn1 +=  weight*A*inv_at2*dm;
  *fEp2 +=  weight*A*inv_tau;
  return;
}

inline void add_ApEnp_coef(double* fAp1, double* fMp1,
			   double* fEn2, double* fEp2,
			   const double tau1, const double dm,
			   const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_ApEn_coef(fAp1, fMp1, fEn2, tau1, dm, tau2_n, weight2*tau2_n);
  add_ApEp_coef(fAp1, fMp1, fEp2, tau1, dm, tau2_p, weight2*tau2_p);
}  

inline void add_AnEnp_coef(double* fAn1, double* fMn1,
			   double* fEn2, double* fEp2,
			   const double tau1, const double dm,
			   const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_AnEn_coef(fAn1, fMn1, fEn2, tau1, dm, tau2_n, weight2*tau2_n);
  add_AnEp_coef(fAn1, fMn1, fEp2, tau1, dm, tau2_p, weight2*tau2_p);
}  

inline void add_ApEf_coef(double* fAp1, double* fMp1,
			  double* fEn2, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_ApEn_coef(fAp1, fMp1, fEn2, tau1, dm, tau2, weight2);
  add_ApEp_coef(fAp1, fMp1, fEp2, tau1, dm, tau2, weight2);
}  

inline void add_AnEf_coef(double* fAn1, double* fMn1,
			  double* fEn2, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_AnEn_coef(fAn1, fMn1, fEn2, tau1, dm, tau2, weight2);
  add_AnEp_coef(fAn1, fMn1, fEp2, tau1, dm, tau2, weight2);
}  

inline void add_AfEp_coef(double* fAn1, double* fAp1, 
			  double* fMn1, double* fMp1, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = /* 0.5* */ weight;  
  add_ApEp_coef(fAp1, fMp1, fEp2, tau1, dm, tau2, weight2);
  add_AnEp_coef(fAn1, fMn1, fEp2, tau1, dm, tau2, weight2);
  return;
}

inline void add_AfEn_coef(double* fAn1, double* fAp1, 
			  double* fMn1, double* fMp1, double* fEn2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = /* 0.5* */ weight;  
  add_ApEn_coef(fAp1, fMp1, fEn2, tau1, dm, tau2, weight2);
  add_AnEn_coef(fAn1, fMn1, fEn2, tau1, dm, tau2, weight2);
  return;
}

inline void add_AfEnp_coef(double* fAn1, double* fAp1, 
			   double* fMn1, double* fMp1,
			   double* fEn2, double* fEp2,
			   const double tau1, const double dm,
			   const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2   = /* 0.5* */ norm*weight;  
  const double weight2_n = weight2*tau2_n;  
  const double weight2_p = weight2*tau2_p;  
  add_ApEn_coef(fAp1, fMp1, fEn2, tau1, dm, tau2_n, weight2_n);
  add_AnEn_coef(fAn1, fMn1, fEn2, tau1, dm, tau2_n, weight2_n);
  add_ApEp_coef(fAp1, fMp1, fEp2, tau1, dm, tau2_p, weight2_p);
  add_AnEp_coef(fAn1, fMn1, fEp2, tau1, dm, tau2_p, weight2_p);
  return;
}

inline void add_AfEf_coef(double* fAn1, double* fAp1, 
			  double* fMn1, double* fMp1,
			  double* fEn2, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = /* 0.25* */ 0.5*weight;  
  add_ApEn_coef(fAp1, fMp1, fEn2, tau1, dm, tau2, weight2);
  add_AnEn_coef(fAn1, fMn1, fEn2, tau1, dm, tau2, weight2);
  add_ApEp_coef(fAp1, fMp1, fEp2, tau1, dm, tau2, weight2);
  add_AnEp_coef(fAn1, fMn1, fEp2, tau1, dm, tau2, weight2);
  return;
}

inline void add_MpEnp_coef(double* fMp1, double* fAp1,
			   double* fEn2, double* fEp2,
			   const double tau1, const double dm,
			   const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_MpEn_coef(fMp1, fAp1, fEn2, tau1, dm, tau2_n, weight2*tau2_n);
  add_MpEp_coef(fMp1, fAp1, fEp2, tau1, dm, tau2_p, weight2*tau2_p);
}  

inline void add_MnEnp_coef(double* fMn1, double* fAn1,
			   double* fEn2, double* fEp2,
			   const double tau1, const double dm,
			   const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2 = norm*weight;  
  add_MnEn_coef(fMn1, fAn1, fEn2, tau1, dm, tau2_n, weight2*tau2_n);
  add_MnEp_coef(fMn1, fAn1, fEp2, tau1, dm, tau2_p, weight2*tau2_p);
}  

inline void add_MpEf_coef(double* fMp1, double* fAp1,
			  double* fEn2, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_MpEn_coef(fMp1, fAp1, fEn2, tau1, dm, tau2, weight2);
  add_MpEp_coef(fMp1, fAp1, fEp2, tau1, dm, tau2, weight2);
}  

inline void add_MnEf_coef(double* fMn1, double* fAn1,
			  double* fEn2, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = 0.5*weight;  
  add_MnEn_coef(fMn1, fAn1, fEn2, tau1, dm, tau2, weight2);
  add_MnEp_coef(fMn1, fAn1, fEp2, tau1, dm, tau2, weight2);
}  

inline void add_MfEp_coef(double* fMn1, double* fMp1, 
			  double* fAn1, double* fAp1, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = /* 0.5* */ weight;  
  add_MpEp_coef(fMp1, fAp1, fEp2, tau1, dm, tau2, weight2);
  add_MnEp_coef(fMn1, fAn1, fEp2, tau1, dm, tau2, weight2);
  return;
}

inline void add_MfEn_coef(double* fMn1, double* fMp1, 
			  double* fAn1, double* fAp1, double* fEn2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = /* 0.5* */ weight;  
  add_MpEn_coef(fMp1, fAp1, fEn2, tau1, dm, tau2, weight2);
  add_MnEn_coef(fMn1, fAn1, fEn2, tau1, dm, tau2, weight2);
  return;
}

inline void add_MfEnp_coef(double* fMn1, double* fMp1, 
			   double* fAn1, double* fAp1,
			   double* fEn2, double* fEp2,
			   const double tau1, const double dm,
			   const double tau2_n, const double tau2_p, const double weight = 1.0){
  const double norm   = 1.0/(std::abs(tau2_n)+std::abs(tau2_p));
  const double weight2   = /* 0.5* */ norm*weight;  
  const double weight2_n = weight2*tau2_n;  
  const double weight2_p = weight2*tau2_p;  
  add_MpEn_coef(fMp1, fAp1, fEn2, tau1, dm, tau2_n, weight2_n);
  add_MnEn_coef(fMn1, fAn1, fEn2, tau1, dm, tau2_n, weight2_n);
  add_MpEp_coef(fMp1, fAp1, fEp2, tau1, dm, tau2_p, weight2_p);
  add_MnEp_coef(fMn1, fAn1, fEp2, tau1, dm, tau2_p, weight2_p);
  return;
}

inline void add_MfEf_coef(double* fMn1, double* fMp1, 
			  double* fAn1, double* fAp1,
			  double* fEn2, double* fEp2,
			  const double tau1, const double dm,
			  const double tau2, const double weight = 1.0){
  const double weight2 = /* 0.25* */ 0.5* weight;  
  add_MpEn_coef(fMp1, fAp1, fEn2, tau1, dm, tau2, weight2);
  add_MnEn_coef(fMn1, fAn1, fEn2, tau1, dm, tau2, weight2);
  add_MpEp_coef(fMp1, fAp1, fEp2, tau1, dm, tau2, weight2);
  add_MnEp_coef(fMn1, fAn1, fEp2, tau1, dm, tau2, weight2);
  return;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __CONV_COEF_H__ */

