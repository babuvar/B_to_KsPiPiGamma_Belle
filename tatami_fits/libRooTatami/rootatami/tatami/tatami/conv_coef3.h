//
//   conv_coef3.h
//     --- Description: Caliculation of coefficient at convolution calculation(3functions)
//   $Id: conv_coef3.h 9932 2006-11-12 14:26:53Z katayama $
//
//   $Log$
//   Revision 1.1  2002/09/05 01:27:06  katayama
//   New package tatami from Nakadaira/Sumisawa san
//
//   Revision 1.2  2002/05/10 13:27:58  nakadair
//   Omit definition of unused parameters.
//
//   Revision 1.1  2002/05/07 21:38:23  nakadair
//   Add convolution with Resolution function and fortran interface.
//
//

#ifndef __CONV_COEF3_H__
#define __CONV_COEF3_H__

#include "belle.h"
#include "tatami/conv_coef.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


inline void add_EpEnEn_coef(double* fEp1, double* fEn2,
			    double* fEn3, double* fxEn2,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEp1=0.0, _fEn2=0.0;
  add_EpEn_coef(&_fEp1, &_fEn2, tau1, tau2);
  add_EpEn_coef(fEp1, fEn3,        tau1, tau3, weight*_fEp1);
  add_EnEn_coef(fEn2, fEn3, fxEn2, tau2, tau3, weight*_fEn2);
  return;
}

inline void add_EnEpEn_coef(double* fEn1, double* fEp2,
			    double* fEn3, double* fxEn1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEp2=0.0;
  add_EnEp_coef(&_fEn1, &_fEp2, tau1, tau2);
  
  add_EnEn_coef(fEn1, fEn3, fxEn1, tau1, tau3, weight*_fEn1);
  add_EpEn_coef(fEp2, fEn3,        tau2, tau3, weight*_fEp2);
  return;
}

inline void add_EnEnEn_coef(double* fEn1, double* fEn2, double* fEn3,
			    double* fxEn1, double* fxEn2, double* fxxEn1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEn2=0.0, _fxEn1=0.0;
  add_EnEn_coef(&_fEn1, &_fEn2, &_fxEn1, tau1, tau2);
  
  add_EnEn_coef (fEn1,  fEn3, fxEn1,         tau1, tau3, weight*_fEn1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2,         tau2, tau3, weight*_fEn2);
  add_xEnEn_coef(fxEn1, fEn1, fEn3,  fxxEn1, tau1, tau3, weight*_fxEn1);
  return;
}

inline void add_EpEpEn_coef(double* fEp1, double* fEp2, double* fEn3, double* fxEp1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){ 
  double _fEp1=0.0, _fEp2=0.0, _fxEp1=0.0;
  add_EpEp_coef(&_fEp1, &_fEp2, &_fxEp1, tau1, tau2);

  add_EpEn_coef (fEp1,  fEn3,       tau1, tau3, weight*_fEp1);
  add_EpEn_coef (fEp2,  fEn3,       tau2, tau3, weight*_fEp2);
  add_xEpEn_coef(fxEp1, fEp1, fEn3, tau1, tau3, weight*_fxEp1);
  return;
}

inline void add_EpEnpEn_coef(double* fEp1, double* fEn2, double* fEp2, double* fEn3,
			     double* fxEp1, double* fxEn2,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){
  double _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEp1=0.0;
  add_EpEnp_coef(&_fEp1, &_fEn2, &_fEp2, &_fxEp1, tau1, tau2_n, tau2_p);

  add_EpEn_coef (fEp1,  fEn3,        tau1,   tau3, weight*_fEp1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2, tau2_n, tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,        tau2_p, tau3, weight*_fEp2);
  add_xEpEn_coef(fxEp1, fEp1, fEn3,  tau1,   tau3, weight*_fxEp1);
  return;
}

inline void add_EpEfEn_coef(double* fEp1, double* fEn2, double* fEp2, double* fEn3,
			    double* fxEp1, double* fxEn2,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEp1=0.0;
  add_EpEf_coef(&_fEp1, &_fEn2, &_fEp2, &_fxEp1, tau1, tau2);

  add_EpEn_coef (fEp1,  fEn3,        tau1, tau3, weight*_fEp1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2, tau2, tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,        tau2, tau3, weight*_fEp2);
  add_xEpEn_coef(fxEp1, fEp1, fEn3,  tau1, tau3, weight*_fxEp1);
  return;
}

inline void add_EnEnpEn_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn1, double* fxEn2, double* fxxEn1,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){ 
  double _fEn1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0;
  add_EnEnp_coef(&_fEn1, &_fEn2, &_fEp2, &_fxEn1, tau1, tau2_n, tau2_p);

  add_EnEn_coef (fEn1,  fEn3, fxEn1,         tau1,   tau3, weight*_fEn1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2,         tau2_n, tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,                tau2_p, tau3, weight*_fEp2);
  add_xEnEn_coef(fxEn1, fEn1, fEn3,  fxxEn1, tau1,   tau3, weight*_fxEn1);
  return;
}

inline void add_EnEfEn_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn1, double* fxEn2, double* fxxEn1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0;
  add_EnEf_coef(&_fEn1, &_fEn2, &_fEp2, &_fxEn1, tau1, tau2);

  add_EnEn_coef (fEn1,  fEn3, fxEn1,         tau1, tau3, weight*_fEn1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2,         tau2, tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,                tau2, tau3, weight*_fEp2);
  add_xEnEn_coef(fxEn1, fEn1, fEn3,  fxxEn1, tau1, tau3, weight*_fxEn1);
  return;
}

inline void add_EnpEnpEn_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3,
			      double* fxEn1, double* fxEp1, double* fxEn2, double* fxxEn1,
			      const double tau1_n, const double tau1_p, const double tau2_n, const double tau2_p,
			      const double tau3, const double weight = 1.0){

  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;

  add_EnpEnp_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1,
		  tau1_n, tau1_p, tau2_n, tau2_p);

  add_EnEn_coef (fEn1,  fEn3, fxEn1,         tau1_n, tau3, weight*_fEn1);
  add_EpEn_coef (fEp1,  fEn3,                tau1_p, tau3, weight*_fEp1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2,         tau2_n, tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,                tau2_p, tau3, weight*_fEp2);
  add_xEnEn_coef(fxEn1, fEn1, fEn3,  fxxEn1, tau1_n, tau3, weight*_fxEn1);
  add_xEpEn_coef(fxEp1, fEn1, fEn3,          tau1_p, tau3, weight*_fxEp1);
  return;
}

inline void add_EfEnpEn_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn1, double* fxEp1, double* fxEn2, double* fxxEn1,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){
  
  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;
  add_EfEnp_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1, tau1, tau2_n, tau2_p);

  add_EnEn_coef (fEn1,  fEn3, fxEn1,         tau1,   tau3, weight*_fEn1);
  add_EpEn_coef (fEp1,  fEn3,                tau1,   tau3, weight*_fEp1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2,         tau2_n, tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,                tau2_p, tau3, weight*_fEp2);
  add_xEnEn_coef(fxEn1, fEn1, fEn3,  fxxEn1, tau1,   tau3, weight*_fxEn1);
  add_xEpEn_coef(fxEp1, fEn1, fEn3,          tau1,   tau3, weight*_fxEp1);
  return;
  return;
}


inline void add_EnpEfEn_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn1, double* fxEp1, double* fxEn2, double* fxxEn1,
			     const double tau1_n, const double tau1_p, const double tau2,
			     const double tau3, const double weight = 1.0){

  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;
  add_EnpEf_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1,
		 tau1_n, tau1_p, tau2);

  add_EnEn_coef (fEn1,  fEn3, fxEn1,         tau1_n, tau3, weight*_fEn1);
  add_EpEn_coef (fEp1,  fEn3,                tau1_p, tau3, weight*_fEp1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2,         tau2,   tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,                tau2,   tau3, weight*_fEp2);
  add_xEnEn_coef(fxEn1, fEn1, fEn3,  fxxEn1, tau1_n, tau3, weight*_fxEn1);
  add_xEpEn_coef(fxEp1, fEn1, fEn3,          tau1_p, tau3, weight*_fxEp1);
  return;
}

inline void add_EfEfEn_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn1, double* fxEp1, double* fxEn2, double* fxxEn1,
			    const double tau1, const double tau2, const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;
  add_EfEf_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1, tau1, tau2);
  
  add_EnEn_coef (fEn1,  fEn3, fxEn1,         tau1, tau3, weight*_fEn1);
  add_EpEn_coef (fEp1,  fEn3,                tau1, tau3, weight*_fEp1);
  add_EnEn_coef (fEn2,  fEn3, fxEn2,         tau2, tau3, weight*_fEn2);
  add_EpEn_coef (fEp2,  fEn3,                tau2, tau3, weight*_fEp2);
  add_xEnEn_coef(fxEn1, fEn1, fEn3,  fxxEn1, tau1, tau3, weight*_fxEn1);
  add_xEpEn_coef(fxEp1, fEn1, fEn3,          tau1, tau3, weight*_fxEp1);
  return;
}

inline void add_ApEpEn_coef(double* fAp1, double* fMp1, double* fEp2, double* fEn3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAp1=0.0, _fMp1=0.0, _fEp2=0.0;
  add_ApEp_coef(&_fAp1, &_fMp1, &_fEp2, tau1, dm, tau2);

  add_ApEn_coef(fAp1, fMp1, fEn3, tau1, dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3, tau1, dm, tau3, weight*_fMp1);
  add_EpEn_coef(fEp2,       fEn3, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_AnEnEn_coef(double* fAn1, double* fMn1, double* fEn2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEn2=0.0;
  add_AnEn_coef(&_fAn1, &_fMn1, &_fEn2, tau1, dm, tau2);

  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_ApEnEn_coef(double* fAp1, double* fMp1, double* fEn2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAp1=0.0, _fMp1=0.0, _fEn2=0.0;
  add_ApEn_coef(&_fAp1, &_fMp1, &_fEn2, tau1, dm, tau2);

  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_AnEpEn_coef(double* fAn1, double* fMn1, double* fEp2, double* fEn3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEp2=0.0;
  add_AnEp_coef(&_fAn1, &_fMn1, &_fEp2, tau1, dm, tau2);

  add_AnEn_coef(fAn1, fMn1, fEn3, tau1, dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3, tau1, dm, tau3, weight*_fMn1);
  add_EpEn_coef(fEp2,       fEn3, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MpEpEn_coef(double* fMp1, double* fAp1, double* fEp2, double* fEn3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEp2=0.0;
  add_MpEp_coef(&_fMp1, &_fAp1, &_fEp2, tau1, dm, tau2);

  add_MpEn_coef(fMp1, fAp1, fEn3, tau1, dm, tau3, weight*_fMp1);
  add_ApEn_coef(fAp1, fMp1, fEn3, tau1, dm, tau3, weight*_fAp1);
  add_EpEn_coef(fEp2,       fEn3, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MnEnEn_coef(double* fMn1, double* fAn1, double* fEn2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEn2=0.0;
  add_MnEn_coef(&_fMn1, &_fAn1, &_fEn2, tau1, dm, tau2);

  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_MpEnEn_coef(double* fMp1, double* fAp1, double* fEn2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEn2=0.0;
  add_MpEn_coef(&_fMp1, &_fAp1, &_fEn2, tau1, dm, tau2);

  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_MnEpEn_coef(double* fMn1, double* fAn1, double* fEp2, double* fEn3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEp2=0.0;
  add_MnEp_coef(&_fMn1, &_fAn1, &_fEp2, tau1, dm, tau2);

  add_MnEn_coef(fMn1, fAn1, fEn3, tau1, dm, tau3, weight*_fMn1);
  add_AnEn_coef(fAn1, fMn1, fEn3, tau1, dm, tau3, weight*_fAn1);
  add_EpEn_coef(fEp2,       fEn3, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_ApEnpEn_coef(double* fAp1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){

  double _fAp1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_ApEnp_coef(&_fAp1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1,   dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1,   dm, tau3, weight*_fMp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2_n,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2_p,     tau3, weight*_fEp2);
  return;
}  

inline void add_AnEnpEn_coef(double* fAn1, double* fMn1,
			     double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AnEnp_coef(&_fAn1, &_fMn1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1,   dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1,   dm, tau3, weight*_fMn1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2_n,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2_p,     tau3, weight*_fEp2);
  return;
}  

inline void add_ApEfEn_coef(double* fAp1, double* fMp1,
			    double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAp1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_ApEf_coef(&_fAp1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2,     tau3, weight*_fEp2);
  return;
}  

inline void add_AnEfEn_coef(double* fAn1, double* fMn1,
			    double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AnEf_coef(&_fAn1, &_fMn1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2,     tau3, weight*_fEp2);
  return;
}  

inline void add_AfEpEn_coef(double* fAn1, double* fAp1, 
			    double* fMn1, double* fMp1, double* fEp2, double* fEn3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEp2=0.0;
  add_AfEp_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEp2, tau1, dm, tau2);

  add_AnEn_coef(fAn1, fMn1, fEn3, tau1, dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3, tau1, dm, tau3, weight*_fMn1);
  add_ApEn_coef(fAp1, fMp1, fEn3, tau1, dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3, tau1, dm, tau3, weight*_fMp1);
  add_EpEn_coef(fEp2,       fEn3, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_AfEnEn_coef(double* fAn1, double* fAp1, 
			    double* fMn1, double* fMp1, double* fEn2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEn2=0.0;
  add_AfEn_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEn2, tau1, dm, tau2);

  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_AfEnpEn_coef(double* fAn1, double* fAp1, 
			     double* fMn1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){

  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AfEnp_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1,   dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1,   dm, tau3, weight*_fMn1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1,   dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1,   dm, tau3, weight*_fMp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2_n,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2_p,     tau3, weight*_fEp2);
  return;
}

inline void add_AfEfEn_coef(double* fAn1, double* fAp1, 
			    double* fMn1, double* fMp1,
			    double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AfEf_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2);
  
  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MpEnpEn_coef(double* fMp1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MpEnp_coef(&_fMp1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1,   dm, tau3, weight*_fMp1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1,   dm, tau3, weight*_fAp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2_n,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2_p,     tau3, weight*_fEp2);

  return;
}  

inline void add_MnEnpEn_coef(double* fMn1, double* fAn1,
			     double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MnEnp_coef(&_fMn1, &_fAn1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1,   dm, tau3, weight*_fMn1);
  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1,   dm, tau3, weight*_fAn1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2_n,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2_p,     tau3, weight*_fEp2);
  return;
}  

inline void add_MpEfEn_coef(double* fMp1, double* fAp1,
			    double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MpEf_coef(&_fMp1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2,     tau3, weight*_fEp2);

  return;
}  

inline void add_MnEfEn_coef(double* fMn1, double* fAn1,
			    double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MnEf_coef(&_fMn1, &_fAn1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2,     tau3, weight*_fEp2);
  return;
}  

inline void add_MfEpEn_coef(double* fMn1, double* fMp1, 
			    double* fAn1, double* fAp1, double* fEp2, double* fEn3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEp2=0.0;
  add_MfEp_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEp2, tau1, dm, tau2);

  add_MnEn_coef(fMn1, fAn1, fEn3,  tau1, dm, tau3, weight*_fMn1);
  add_MpEn_coef(fMp1, fAp1, fEn3,  tau1, dm, tau3, weight*_fMp1);
  add_AnEn_coef(fAn1, fMn1, fEn3,  tau1, dm, tau3, weight*_fAn1);
  add_ApEn_coef(fAp1, fMp1, fEn3,  tau1, dm, tau3, weight*_fAp1);
  add_EpEn_coef(fEp2,       fEn3,  tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MfEnEn_coef(double* fMn1, double* fMp1, 
			    double* fAn1, double* fAp1, double* fEn2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEn2=0.0;
  add_MfEn_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEn2, tau1, dm, tau2);

  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_MfEnpEn_coef(double* fMn1, double* fMp1, 
			     double* fAn1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEn3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MfEnp_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1,   dm, tau3, weight*_fMn1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1,   dm, tau3, weight*_fMp1);
  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1,   dm, tau3, weight*_fAn1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1,   dm, tau3, weight*_fAp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2_n,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2_p,     tau3, weight*_fEp2);
  return;
}

inline void add_MfEfEn_coef(double* fMn1, double* fMp1, 
			    double* fAn1, double* fAp1,
			    double* fEn2, double* fEp2, double* fEn3,
			    double* fxEn2,
			    const double tau1, const double dm,	
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MfEf_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_MnEn_coef(fMn1, fAn1, fEn3,        tau1, dm, tau3, weight*_fMn1);
  add_MpEn_coef(fMp1, fAp1, fEn3,        tau1, dm, tau3, weight*_fMp1);
  add_AnEn_coef(fAn1, fMn1, fEn3,        tau1, dm, tau3, weight*_fAn1);
  add_ApEn_coef(fAp1, fMp1, fEn3,        tau1, dm, tau3, weight*_fAp1);
  add_EnEn_coef(fEn2,       fEn3, fxEn2, tau2,     tau3, weight*_fEn2);
  add_EpEn_coef(fEp2,       fEn3,        tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_EpEnEp_coef(double* fEp1, double* fEn2,
			    double* fEp3, double* fxEp1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEp1=0.0, _fEn2=0.0;
  add_EpEn_coef(&_fEp1, &_fEn2, tau1, tau2);
  add_EpEp_coef(fEp1, fEp3, fxEp1, tau1, tau3, weight*_fEp1);
  add_EnEp_coef(fEn2, fEp3,        tau2, tau3, weight*_fEn2);
  return;
}

inline void add_EnEpEp_coef(double* fEn1, double* fEp2,
			    double* fEp3, double* fxEp2,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEp2=0.0;
  add_EnEp_coef(&_fEn1, &_fEp2, tau1, tau2);
  
  add_EnEp_coef(fEn1, fEp3,        tau1, tau3, weight*_fEn1);
  add_EpEp_coef(fEp2, fEp3, fxEp2, tau2, tau3, weight*_fEp2);
  return;
}

inline void add_EnEnEp_coef(double* fEn1, double* fEn2, double* fEp3,
			    double* fxEn1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEn2=0.0, _fxEn1=0.0;
  add_EnEn_coef(&_fEn1, &_fEn2, &_fxEn1, tau1, tau2);
  
  add_EnEp_coef (fEn1,        fEp3, tau1, tau3, weight*_fEn1);
  add_EnEp_coef (fEn2,        fEp3, tau2, tau3, weight*_fEn2);
  add_xEnEp_coef(fxEn1, fEn1, fEp3, tau1, tau3, weight*_fxEn1);
  return;
}

inline void add_EpEpEp_coef(double* fEp1, double* fEp2, double* fEp3,
			    double* fxEp1, double* fxEp2, double* fxxEp1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){ 
  double _fEp1=0.0, _fEp2=0.0, _fxEp1=0.0;
  add_EpEp_coef(&_fEp1, &_fEp2, &_fxEp1, tau1, tau2);

  add_EpEp_coef (fEp1,  fEp3, fxEp1,         tau1, tau3, weight*_fEp1);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,         tau2, tau3, weight*_fEp2);
  add_xEpEp_coef(fxEp1, fEp1, fEp3,  fxxEp1, tau1, tau3, weight*_fxEp1);
  return;
}

inline void add_EpEnpEp_coef(double* fEp1, double* fEn2, double* fEp2, double* fEp3,
			     double* fxEp1, double* fxEp2, double* fxxEp1, 
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){
  double _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEp1=0.0;
  add_EpEnp_coef(&_fEp1, &_fEn2, &_fEp2, &_fxEp1, tau1, tau2_n, tau2_p);

  add_EpEp_coef (fEp1,  fEp3, fxEp1,        tau1,   tau3, weight*_fEp1);
  add_EnEp_coef (fEn2,  fEp3,               tau2_n, tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,        tau2_p, tau3, weight*_fEp2);
  add_xEpEp_coef(fxEp1, fEp1, fEp3, fxxEp1, tau1,   tau3, weight*_fxEp1);
  return;
}

inline void add_EpEfEp_coef(double* fEp1, double* fEn2, double* fEp2, double* fEp3,
			    double* fxEp1, double* fxEp2, double* fxxEp1, 
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEp1=0.0;
  add_EpEf_coef(&_fEp1, &_fEn2, &_fEp2, &_fxEp1, tau1, tau2);

  add_EpEp_coef (fEp1,  fEp3, fxEp1,        tau1, tau3, weight*_fEp1);
  add_EnEp_coef (fEn2,  fEp3,               tau2, tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,        tau2, tau3, weight*_fEp2);
  add_xEpEp_coef(fxEp1, fEp1, fEp3, fxxEp1, tau1, tau3, weight*_fxEp1);
  return;
}

inline void add_EnEnpEp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEp3,
			     double* fxEn1, double* fxEp2,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){ 
  double _fEn1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0;
  add_EnEnp_coef(&_fEn1, &_fEn2, &_fEp2, &_fxEn1, tau1, tau2_n, tau2_p);

  add_EnEp_coef (fEn1,  fEp3,         tau1,   tau3, weight*_fEn1);
  add_EnEp_coef (fEn2,  fEp3,         tau2_n, tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,  tau2_p, tau3, weight*_fEp2);
  add_xEnEp_coef(fxEn1, fEn1, fEp3,   tau1,   tau3, weight*_fxEn1);
  return;
}

inline void add_EnEfEp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEp3,
			    double* fxEn1, double* fxEp2,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0;
  add_EnEf_coef(&_fEn1, &_fEn2, &_fEp2, &_fxEn1, tau1, tau2);

  add_EnEp_coef (fEn1,  fEp3,         tau1, tau3, weight*_fEn1);
  add_EnEp_coef (fEn2,  fEp3,         tau2, tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,  tau2, tau3, weight*_fEp2);
  add_xEnEp_coef(fxEn1, fEn1, fEp3,   tau1, tau3, weight*_fxEn1);
  return;
}

inline void add_EnpEnpEp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEp3,
			      double* fxEn1, double* fxEp1, double* fxEp2, double* fxxEp1,
			      const double tau1_n, const double tau1_p, const double tau2_n, const double tau2_p,
			      const double tau3, const double weight = 1.0){

  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;
  add_EnpEnp_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1,
		  tau1_n, tau1_p, tau2_n, tau2_p);

  add_EnEp_coef (fEn1,  fEp3,                tau1_n, tau3, weight*_fEn1);
  add_EpEp_coef (fEp1,  fEp3, fxEp1,         tau1_p, tau3, weight*_fEp1);
  add_EnEp_coef (fEn2,  fEp3,                tau2_n, tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,         tau2_p, tau3, weight*_fEp2);
  add_xEnEp_coef(fxEn1, fEn1, fEp3,          tau1_n, tau3, weight*_fxEn1);
  add_xEpEp_coef(fxEp1, fEn1, fEp3,  fxxEp1, tau1_p, tau3, weight*_fxEp1);
  return;
}

inline void add_EfEnpEp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEp3,
			     double* fxEn1, double* fxEp1, double* fxEp2, double* fxxEp1,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){

  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;
  add_EfEnp_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1, tau1, tau2_n, tau2_p);

  add_EnEp_coef (fEn1,  fEp3,                tau1,   tau3, weight*_fEn1);
  add_EpEp_coef (fEp1,  fEp3, fxEp1,         tau1,   tau3, weight*_fEp1);
  add_EnEp_coef (fEn2,  fEp3,                tau2_n, tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,         tau2_p, tau3, weight*_fEp2);
  add_xEnEp_coef(fxEn1, fEn1, fEp3,          tau1,   tau3, weight*_fxEn1);
  add_xEpEp_coef(fxEp1, fEn1, fEp3,  fxxEp1, tau1,   tau3, weight*_fxEp1);
  return;
}

inline void add_EnpEfEp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEp3,
			     double* fxEn1, double* fxEp1, double* fxEp2, double* fxxEp1,
			     const double tau1_n, const double tau1_p, const double tau2,
			     const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;
  add_EnpEf_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1, tau1_n, tau1_p, tau2);

  add_EnEp_coef (fEn1,  fEp3,                tau1_n, tau3, weight*_fEn1);
  add_EpEp_coef (fEp1,  fEp3, fxEp1,         tau1_p, tau3, weight*_fEp1);
  add_EnEp_coef (fEn2,  fEp3,                tau2,   tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,         tau2,   tau3, weight*_fEp2);
  add_xEnEp_coef(fxEn1, fEn1, fEp3,          tau1_n, tau3, weight*_fxEn1);
  add_xEpEp_coef(fxEp1, fEn1, fEp3,  fxxEp1, tau1_p, tau3, weight*_fxEp1);
  return;
}

inline void add_EfEfEp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEp3,
			    double* fxEn1, double* fxEp1, double* fxEp2, double* fxxEp1,
			    const double tau1, const double tau2, const double tau3, const double weight = 1.0){
  double _fEn1=0.0, _fEp1=0.0, _fEn2=0.0, _fEp2=0.0, _fxEn1=0.0, _fxEp1=0.0;
  add_EfEf_coef(&_fEn1, &_fEp1, &_fEn2, &_fEp2, &_fxEn1, &_fxEp1, tau1, tau2);
  
  add_EnEp_coef (fEn1,  fEp3,                tau1, tau3, weight*_fEn1);
  add_EpEp_coef (fEp1,  fEp3, fxEp1,         tau1, tau3, weight*_fEp1);
  add_EnEp_coef (fEn2,  fEp3,                tau2, tau3, weight*_fEn2);
  add_EpEp_coef (fEp2,  fEp3, fxEp2,         tau2, tau3, weight*_fEp2);
  add_xEnEp_coef(fxEn1, fEn1, fEp3,          tau1, tau3, weight*_fxEn1);
  add_xEpEp_coef(fxEp1, fEn1, fEp3,  fxxEp1, tau1, tau3, weight*_fxEp1);
  return;
}

inline void add_ApEpEp_coef(double* fAp1, double* fMp1, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAp1=0.0, _fMp1=0.0, _fEp2=0.0;
  add_ApEp_coef(&_fAp1, &_fMp1, &_fEp2, tau1, dm, tau2);

  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_AnEnEp_coef(double* fAn1, double* fMn1, double* fEn2, double* fEp3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEn2=0.0;
  add_AnEn_coef(&_fAn1, &_fMn1, &_fEn2, tau1, dm, tau2);

  add_AnEp_coef(fAn1, fMn1, fEp3, tau1, dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3, tau1, dm, tau3, weight*_fMn1);
  add_EnEp_coef(fEn2,       fEp3, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_ApEnEp_coef(double* fAp1, double* fMp1, double* fEn2, double* fEp3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAp1=0.0, _fMp1=0.0, _fEn2=0.0;
  add_ApEn_coef(&_fAp1, &_fMp1, &_fEn2, tau1, dm, tau2);

  add_ApEp_coef(fAp1, fMp1, fEp3, tau1, dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3, tau1, dm, tau3, weight*_fMp1);
  add_EnEp_coef(fEn2,       fEp3, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_AnEpEp_coef(double* fAn1, double* fMn1, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEp2=0.0;
  add_AnEp_coef(&_fAn1, &_fMn1, &_fEp2, tau1, dm, tau2);

  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MpEpEp_coef(double* fMp1, double* fAp1, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEp2=0.0;
  add_MpEp_coef(&_fMp1, &_fAp1, &_fEp2, tau1, dm, tau2);

  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MnEnEp_coef(double* fMn1, double* fAn1, double* fEn2, double* fEp3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEn2=0.0;
  add_MnEn_coef(&_fMn1, &_fAn1, &_fEn2, tau1, dm, tau2);

  add_MnEp_coef(fMn1, fAn1, fEp3, tau1, dm, tau3, weight*_fMn1);
  add_AnEp_coef(fAn1, fMn1, fEp3, tau1, dm, tau3, weight*_fAn1);
  add_EnEp_coef(fEn2,       fEp3, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_MpEnEp_coef(double* fMp1, double* fAp1, double* fEn2, double* fEp3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEn2=0.0;
  add_MpEn_coef(&_fMp1, &_fAp1, &_fEn2, tau1, dm, tau2);

  add_MpEp_coef(fMp1, fAp1, fEp3, tau1, dm, tau3, weight*_fMp1);
  add_ApEp_coef(fAp1, fMp1, fEp3, tau1, dm, tau3, weight*_fAp1);
  add_EnEp_coef(fEn2,       fEp3, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_MnEpEp_coef(double* fMn1, double* fAn1, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEp2=0.0;
  add_MnEp_coef(&_fMn1, &_fAn1, &_fEp2, tau1, dm, tau2);

  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_ApEnpEp_coef(double* fAp1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){

  double _fAp1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_ApEnp_coef(&_fAp1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1,   dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1,   dm, tau3, weight*_fMp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2_n,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2_p,     tau3, weight*_fEp2);
  return;
}  

inline void add_AnEnpEp_coef(double* fAn1, double* fMn1,
			     double* fEn2, double* fEp2, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AnEnp_coef(&_fAn1, &_fMn1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1,   dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1,   dm, tau3, weight*_fMn1);
  add_EnEp_coef(fEn2,       fEp3,        tau2_n,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2_p,     tau3, weight*_fEp2);
  return;
}  

inline void add_ApEfEp_coef(double* fAp1, double* fMp1,
			    double* fEn2, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAp1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_ApEf_coef(&_fAp1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}  

inline void add_AnEfEp_coef(double* fAn1, double* fMn1,
			    double* fEn2, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fMn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AnEf_coef(&_fAn1, &_fMn1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_EnEp_coef(fEn2,       fEp3,        tau2,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}  

inline void add_AfEpEp_coef(double* fAn1, double* fAp1, 
			    double* fMn1, double* fMp1, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEp2=0.0;
  add_AfEp_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEp2, tau1, dm, tau2);

  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_AfEnEp_coef(double* fAn1, double* fAp1, 
			    double* fMn1, double* fMp1, double* fEn2, double* fEp3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEn2=0.0;
  add_AfEn_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEn2, tau1, dm, tau2);

  add_AnEp_coef(fAn1, fMn1, fEp3, tau1, dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3, tau1, dm, tau3, weight*_fMn1);
  add_ApEp_coef(fAp1, fMp1, fEp3, tau1, dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3, tau1, dm, tau3, weight*_fMp1);
  add_EnEp_coef(fEn2,       fEp3, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_AfEnpEp_coef(double* fAn1, double* fAp1, 
			     double* fMn1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){

  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AfEnp_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1,   dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1,   dm, tau3, weight*_fMn1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1,   dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1,   dm, tau3, weight*_fMp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2_n,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2_p,     tau3, weight*_fEp2);
  return;
}

inline void add_AfEfEp_coef(double* fAn1, double* fAp1, 
			    double* fMn1, double* fMp1,
			    double* fEn2, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fAn1=0.0, _fAp1=0.0, _fMn1=0.0, _fMp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_AfEf_coef(&_fAn1, &_fAp1, &_fMn1, &_fMp1, &_fEn2, &_fEp2, tau1, dm, tau2);
  
  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MpEnpEp_coef(double* fMp1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MpEnp_coef(&_fMp1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1,   dm, tau3, weight*_fMp1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1,   dm, tau3, weight*_fAp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2_n,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2_p,     tau3, weight*_fEp2);

  return;
}  

inline void add_MnEnpEp_coef(double* fMn1, double* fAn1,
			     double* fEn2, double* fEp2, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MnEnp_coef(&_fMn1, &_fAn1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1,   dm, tau3, weight*_fMn1);
  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1,   dm, tau3, weight*_fAn1);
  add_EnEp_coef(fEn2,       fEp3,        tau2_n,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2_p,     tau3, weight*_fEp2);
  return;
}  

inline void add_MpEfEp_coef(double* fMp1, double* fAp1,
			    double* fEn2, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMp1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MpEf_coef(&_fMp1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);

  return;
}  

inline void add_MnEfEp_coef(double* fMn1, double* fAn1,
			    double* fEn2, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fAn1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MnEf_coef(&_fMn1, &_fAn1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_EnEp_coef(fEn2,       fEp3,        tau2,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}  

inline void add_MfEpEp_coef(double* fMn1, double* fMp1, 
			    double* fAn1, double* fAp1, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEp2=0.0;
  add_MfEp_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEp2, tau1, dm, tau2);

  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_MfEnEp_coef(double* fMn1, double* fMp1, 
			    double* fAn1, double* fAp1, double* fEn2, double* fEp3,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEn2=0.0;
  add_MfEn_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEn2, tau1, dm, tau2);

  add_MnEp_coef(fMn1, fAn1, fEp3, tau1, dm, tau3, weight*_fMn1);
  add_MpEp_coef(fMp1, fAp1, fEp3, tau1, dm, tau3, weight*_fMp1);
  add_AnEp_coef(fAn1, fMn1, fEp3, tau1, dm, tau3, weight*_fAn1);
  add_ApEp_coef(fAp1, fMp1, fEp3, tau1, dm, tau3, weight*_fAp1);
  add_EnEp_coef(fEn2,       fEp3, tau2,     tau3, weight*_fEn2);
  return;
}

inline void add_MfEnpEp_coef(double* fMn1, double* fMp1, 
			     double* fAn1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MfEnp_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2_n, tau2_p);

  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1,   dm, tau3, weight*_fMn1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1,   dm, tau3, weight*_fMp1);
  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1,   dm, tau3, weight*_fAn1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1,   dm, tau3, weight*_fAp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2_n,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2_p,     tau3, weight*_fEp2);
  return;
}

inline void add_MfEfEp_coef(double* fMn1, double* fMp1, 
			    double* fAn1, double* fAp1,
			    double* fEn2, double* fEp2, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,	
			    const double tau2, const double tau3, const double weight = 1.0){
  double _fMn1=0.0, _fMp1=0.0, _fAn1=0.0, _fAp1=0.0, _fEn2=0.0, _fEp2=0.0;
  add_MfEf_coef(&_fMn1, &_fMp1, &_fAn1, &_fAp1, &_fEn2, &_fEp2, tau1, dm, tau2);

  add_MnEp_coef(fMn1, fAn1, fEp3,        tau1, dm, tau3, weight*_fMn1);
  add_MpEp_coef(fMp1, fAp1, fEp3,        tau1, dm, tau3, weight*_fMp1);
  add_AnEp_coef(fAn1, fMn1, fEp3,        tau1, dm, tau3, weight*_fAn1);
  add_ApEp_coef(fAp1, fMp1, fEp3,        tau1, dm, tau3, weight*_fAp1);
  add_EnEp_coef(fEn2,       fEp3,        tau2,     tau3, weight*_fEn2);
  add_EpEp_coef(fEp2,       fEp3, fxEp2, tau2,     tau3, weight*_fEp2);
  return;
}

inline void add_EpEnEnp_coef(double* fEp1, double* fEn2,
			     double* fEn3, double* fEp3,
			     double* fxEp1, double* fxEn2,
			     const double tau1, const double tau2,
			     const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EpEnEn_coef(fEp1, fEn2, fEn3, fxEn2, tau1, tau2, tau3_n, weight2*tau3_n);
  add_EpEnEp_coef(fEp1, fEn2, fEp3, fxEp1, tau1, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EnEpEnp_coef(double* fEn1, double* fEp2,
			     double* fEn3, double* fEp3,
			     double* fxEn1, double* fxEp2,
			     const double tau1, const double tau2,
			     const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EnEpEn_coef(fEn1, fEp2, fEn3, fxEn1, tau1, tau2, tau3_n, weight2*tau3_n);
  add_EnEpEp_coef(fEn1, fEp2, fEp3, fxEp2, tau1, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EnEnEnp_coef(double* fEn1, double* fEn2, double* fEn3, double* fEp3,
			     double* fxEn1, double* fxEn2, double* fxxEn1,
			     const double tau1, const double tau2,
			     const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EnEnEn_coef(fEn1, fEn2, fEn3, fxEn1, fxEn2, fxxEn1, tau1, tau2, tau3_n, weight2*tau3_n);
  add_EnEnEp_coef(fEn1, fEn2, fEp3, fxEn1,                tau1, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EpEpEnp_coef(double* fEp1, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp1, double* fxEp2, double* fxxEp1,
			     const double tau1, const double tau2,
			     const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EpEpEn_coef(fEp1, fEp2, fEn3, fxEp1,                tau1, tau2, tau3_n, weight2*tau3_n);
  add_EpEpEp_coef(fEp1, fEp2, fEp3, fxEp1, fxEp2, fxxEp1, tau1, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EpEnpEnp_coef(double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEp1,
			      const double tau1, const double tau2_n, const double tau2_p,
			      const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EpEnpEn_coef(fEp1, fEn2, fEp2, fEn3, fxEp1, fxEn2,         tau1, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_EpEnpEp_coef(fEp1, fEn2, fEp2, fEp3, fxEp1, fxEp2, fxxEp1, tau1, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EpEfEnp_coef(double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEp1,
			     const double tau1, const double tau2,
			     const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EpEfEn_coef(fEp1, fEn2, fEp2, fEn3, fxEp1, fxEn2,         tau1, tau2, tau3_n, weight2*tau3_n);
  add_EpEfEp_coef(fEp1, fEn2, fEp2, fEp3, fxEp1, fxEp2, fxxEp1, tau1, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EnEnpEnp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn1, double* fxEn2, double* fxEp2, double* fxxEn1,
			      const double tau1, const double tau2_n, const double tau2_p,
			      const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EnEnpEn_coef(fEn1, fEn2, fEp2, fEn3, fxEn1, fxEn2, fxxEn1, tau1, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_EnEnpEp_coef(fEn1, fEn2, fEp2, fEp3, fxEn1, fxEp2,         tau1, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EnEfEnp_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn1, double* fxEn2, double* fxEp2, double* fxxEn1,
			     const double tau1, const double tau2,
			     const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EnEfEn_coef(fEn1, fEn2, fEp2, fEn3, fxEn1, fxEn2, fxxEn1, tau1, tau2, tau3_n, weight2*tau3_n);
  add_EnEfEp_coef(fEn1, fEn2, fEp2, fEp3, fxEn1, fxEp2,         tau1, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EnpEnpEnp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			       double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			       const double tau1_n, const double tau1_p, const double tau2_n, const double tau2_p,
			       const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EnpEnpEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1_n, tau1_p, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_EnpEnpEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1_n, tau1_p, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EfEnpEnp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			      const double tau1, const double tau2_n, const double tau2_p,
			      const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EfEnpEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_EfEnpEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EnpEfEnp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			      const double tau1_n, const double tau1_p, const double tau2,
			      const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EnpEfEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1_n, tau1_p, tau2, tau3_n, weight2*tau3_n);
  add_EnpEfEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1_n, tau1_p, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_EfEfEnp_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			     const double tau1, const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_EfEfEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1, tau2, tau3_n, weight2*tau3_n);
  add_EfEfEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_ApEpEnp_coef(double* fAp1, double* fMp1, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_ApEpEn_coef(fAp1, fMp1, fEp2, fEn3,        tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_ApEpEp_coef(fAp1, fMp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AnEnEnp_coef(double* fAn1, double* fMn1, double* fEn2, double* fEn3, double* fEp3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AnEnEn_coef(fAn1, fMn1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_AnEnEp_coef(fAn1, fMn1, fEn2, fEp3,        tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_ApEnEnp_coef(double* fAp1, double* fMp1, double* fEn2, double* fEn3, double* fEp3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_ApEnEn_coef(fAp1, fMp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_ApEnEp_coef(fAp1, fMp1, fEn2, fEp3,        tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AnEpEnp_coef(double* fAn1, double* fMn1, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AnEpEn_coef(fAn1, fMn1, fEp2, fEn3,        tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_AnEpEp_coef(fAn1, fMn1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MpEpEnp_coef(double* fMp1, double* fAp1, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MpEpEn_coef(fMp1, fAp1, fEp2, fEn3,        tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MpEpEp_coef(fMp1, fAp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MnEnEnp_coef(double* fMn1, double* fAn1, double* fEn2, double* fEn3, double* fEp3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MnEnEn_coef(fMn1, fAn1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MnEnEp_coef(fMn1, fAn1, fEn2, fEp3,        tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MpEnEnp_coef(double* fMp1, double* fAp1, double* fEn2, double* fEn3, double* fEp3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MpEnEn_coef(fMp1, fAp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MpEnEp_coef(fMp1, fAp1, fEn2, fEp3,        tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MnEpEnp_coef(double* fMn1, double* fAn1, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MnEpEn_coef(fMn1, fAn1, fEp2,       fxEp2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MnEpEp_coef(fMn1, fAn1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_ApEnpEnp_coef(double* fAp1, double* fMp1,
			      double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn2, double* fxEp2,
			      const double tau1, const double dm,
			      const double tau2_n, const double tau2_p, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_ApEnpEn_coef(fAp1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_ApEnpEp_coef(fAp1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AnEnpEnp_coef(double* fAn1, double* fMn1,
			      double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn2, double* fxEp2,
			      const double tau1, const double dm,
			      const double tau2_n, const double tau2_p, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AnEnpEn_coef(fAn1, fMn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_AnEnpEp_coef(fAn1, fMn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_ApEfEnp_coef(double* fAp1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_ApEfEn_coef(fAp1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_ApEfEp_coef(fAp1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AnEfEnp_coef(double* fAn1, double* fMn1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AnEfEn_coef(fAn1, fMn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_AnEfEp_coef(fAn1, fMn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AfEpEnp_coef(double* fAn1, double* fAp1,
			     double* fMn1, double* fMp1, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AfEpEn_coef(fAn1, fAp1, fMn1, fMp1, fEp2, fEn3,        tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_AfEpEp_coef(fAn1, fAp1, fMn1, fMp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AfEnEnp_coef(double* fAn1, double* fAp1,
			     double* fMn1, double* fMp1, double* fEn2, double* fEn3, double* fEp3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AfEnEn_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_AfEnEp_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp3,        tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AfEnpEnp_coef(double* fAn1, double* fAp1,
			      double* fMn1, double* fMp1,
			      double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn2, double* fxEp2,
			      const double tau1, const double dm,
			      const double tau2_n, const double tau2_p,
			      const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AfEnpEn_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_AfEnpEp_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_AfEfEnp_coef(double* fAn1, double* fAp1,
			     double* fMn1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_AfEfEn_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_AfEfEp_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MpEnpEnp_coef(double* fMp1, double* fAp1,
			      double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn2, double* fxEp2,
			      const double tau1, const double dm,
			      const double tau2_n, const double tau2_p, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MpEnpEn_coef(fMp1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_MpEnpEp_coef(fMp1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MnEnpEnp_coef(double* fMn1, double* fAn1,
			      double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn2, double* fxEp2,
			      const double tau1, const double dm,
			      const double tau2_n, const double tau2_p, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MnEnpEn_coef(fMn1, fAn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_MnEnpEp_coef(fMn1, fAn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MpEfEnp_coef(double* fMp1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MpEfEn_coef(fMp1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MpEfEp_coef(fMp1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MnEfEnp_coef(double* fMn1, double* fAn1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MnEfEn_coef(fMn1, fAn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MnEfEp_coef(fMn1, fAn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MfEpEnp_coef(double* fMn1, double* fMp1,
			     double* fAn1, double* fAp1, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MfEpEn_coef(fMn1, fMp1, fAn1, fAp1, fEp2, fEn3,        tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MfEpEp_coef(fMn1, fMp1, fAn1, fAp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MfEnEnp_coef(double* fMn1, double* fMp1,
			     double* fAn1, double* fAp1, double* fEn2, double* fEn3, double* fEp3,
			     double* fxEn2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MfEnEn_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MfEnEp_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp3,        tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MfEnpEnp_coef(double* fMn1, double* fMp1,
			      double* fAn1, double* fAp1,
			      double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn2, double* fxEp2,
			      const double tau1, const double dm,
			      const double tau2_n, const double tau2_p, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MfEnpEn_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3_n, weight2*tau3_n);
  add_MfEnpEp_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3_p, weight2*tau3_p);
  return;
}

inline void add_MfEfEnp_coef(double* fMn1, double* fMp1,
			     double* fAn1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2, const double tau3_n, const double tau3_p, const double weight = 1.0){
  const double norm   = 1.0/(fabs(tau3_n)+fabs(tau3_p));
  const double weight2 = norm*weight;

  add_MfEfEn_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3_n, weight2*tau3_n);
  add_MfEfEp_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3_p, weight2*tau3_p);
  return;
}


inline void add_EpEnEf_coef(double* fEp1, double* fEn2,
			    double* fEn3, double* fEp3,
			    double* fxEp1, double* fxEn2,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EpEnEn_coef(fEp1, fEn2, fEn3, fxEn2, tau1, tau2, tau3, weight2);
  add_EpEnEp_coef(fEp1, fEn2, fEp3, fxEp1, tau1, tau2, tau3, weight2);
  return;
}

inline void add_EnEpEf_coef(double* fEn1, double* fEp2,
			    double* fEn3, double* fEp3,
			    double* fxEn1, double* fxEp2,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EnEpEn_coef(fEn1, fEp2, fEn3, fxEn1, tau1, tau2, tau3, weight2);
  add_EnEpEp_coef(fEn1, fEp2, fEp3, fxEp2, tau1, tau2, tau3, weight2);
  return;
}

inline void add_EnEnEf_coef(double* fEn1, double* fEn2, double* fEn3, double* fEp3,
			    double* fxEn1, double* fxEn2, double* fxxEn1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EnEnEn_coef(fEn1, fEn2, fEn3, fxEn1, fxEn2, fxxEn1, tau1, tau2, tau3, weight2);
  add_EnEnEp_coef(fEn1, fEn2, fEp3, fxEn1,                tau1, tau2, tau3, weight2);
  return;
}

inline void add_EpEpEf_coef(double* fEp1, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp1, double* fxEp2, double* fxxEp1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EpEpEn_coef(fEp1, fEp2, fEn3, fxEp1,                tau1, tau2, tau3, weight2);
  add_EpEpEp_coef(fEp1, fEp2, fEp3, fxEp1, fxEp2, fxxEp1, tau1, tau2, tau3, weight2);
  return;
}

inline void add_EpEnpEf_coef(double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEp1,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EpEnpEn_coef(fEp1, fEn2, fEp2, fEn3, fxEp1, fxEn2,         tau1, tau2_n, tau2_p, tau3, weight2);
  add_EpEnpEp_coef(fEp1, fEn2, fEp2, fEp3, fxEp1, fxEp2, fxxEp1, tau1, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_EpEfEf_coef(double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEp1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EpEfEn_coef(fEp1, fEn2, fEp2, fEn3, fxEp1, fxEn2,         tau1, tau2, tau3, weight2);
  add_EpEfEp_coef(fEp1, fEn2, fEp2, fEp3, fxEp1, fxEp2, fxxEp1, tau1, tau2, tau3, weight2);
  return;
}

inline void add_EnEnpEf_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn1, double* fxEn2, double* fxEp2, double* fxxEn1,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EnEnpEn_coef(fEn1, fEn2, fEp2, fEn3, fxEn1, fxEn2, fxxEn1, tau1, tau2_n, tau2_p, tau3, weight2);
  add_EnEnpEp_coef(fEn1, fEn2, fEp2, fEp3, fxEn1, fxEp2,         tau1, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_EnEfEf_coef(double* fEn1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn1, double* fxEn2, double* fxEp2, double* fxxEn1,
			    const double tau1, const double tau2,
			    const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EnEfEn_coef(fEn1, fEn2, fEp2, fEn3, fxEn1, fxEn2, fxxEn1, tau1, tau2, tau3, weight2);
  add_EnEfEp_coef(fEn1, fEn2, fEp2, fEp3, fxEn1, fxEp2,         tau1, tau2, tau3, weight2);
  return;
}

inline void add_EnpEnpEf_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			      double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			      const double tau1_n, const double tau1_p, const double tau2_n, const double tau2_p,
			      const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EnpEnpEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1_n, tau1_p, tau2_n, tau2_p, tau3, weight2);
  add_EnpEnpEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1_n, tau1_p, tau2_n, tau2_p, tau3, weight2);
  return;
}


inline void add_EfEnpEf_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			     const double tau1, const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EfEnpEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1, tau2_n, tau2_p, tau3, weight2);
  add_EfEnpEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_EnpEfEf_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			     const double tau1_n, const double tau1_p, const double tau2,
			     const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EnpEfEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1_n, tau1_p, tau2, tau3, weight2);
  add_EnpEfEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1_n, tau1_p, tau2, tau3, weight2);
  return;
}

inline void add_EfEfEf_coef(double* fEn1, double* fEp1, double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn1, double* fxEp1, double* fxEn2, double* fxEp2, double* fxxEn1, double* fxxEp1,
			    const double tau1, const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_EfEfEn_coef(fEn1, fEp1, fEn2, fEp2, fEn3, fxEn1, fxEp1, fxEn2, fxxEn1, tau1, tau2, tau3, weight2);
  add_EfEfEp_coef(fEn1, fEp1, fEn2, fEp2, fEp3, fxEn1, fxEp1, fxEp2, fxxEp1, tau1, tau2, tau3, weight2);
  return;
}

inline void add_ApEpEf_coef(double* fAp1, double* fMp1, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_ApEpEn_coef(fAp1, fMp1, fEp2, fEn3,        tau1, dm, tau2, tau3, weight2);
  add_ApEpEp_coef(fAp1, fMp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_AnEnEf_coef(double* fAn1, double* fMn1, double* fEn2, double* fEn3, double* fEp3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AnEnEn_coef(fAn1, fMn1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_AnEnEp_coef(fAn1, fMn1, fEn2, fEp3,        tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_ApEnEf_coef(double* fAp1, double* fMp1, double* fEn2, double* fEn3, double* fEp3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_ApEnEn_coef(fAp1, fMp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_ApEnEp_coef(fAp1, fMp1, fEn2, fEp3,        tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_AnEpEf_coef(double* fAn1, double* fMn1, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AnEpEn_coef(fAn1, fMn1, fEp2, fEn3,        tau1, dm, tau2, tau3, weight2);
  add_AnEpEp_coef(fAn1, fMn1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MpEpEf_coef(double* fMp1, double* fAp1, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MpEpEn_coef(fMp1, fAp1, fEp2, fEn3,        tau1, dm, tau2, tau3, weight2);
  add_MpEpEp_coef(fMp1, fAp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MnEnEf_coef(double* fMn1, double* fAn1, double* fEn2, double* fEn3, double* fEp3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MnEnEn_coef(fMn1, fAn1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_MnEnEp_coef(fMn1, fAn1, fEn2, fEp3,        tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MpEnEf_coef(double* fMp1, double* fAp1, double* fEn2, double* fEn3, double* fEp3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MpEnEn_coef(fMp1, fAp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_MpEnEp_coef(fMp1, fAp1, fEn2, fEp3,        tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MnEpEf_coef(double* fMn1, double* fAn1, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MnEpEn_coef(fMn1, fAn1, fEp2,       fxEp2, tau1, dm, tau2, tau3, weight2);
  add_MnEpEp_coef(fMn1, fAn1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_ApEnpEf_coef(double* fAp1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_ApEnpEn_coef(fAp1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  add_ApEnpEp_coef(fAp1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_AnEnpEf_coef(double* fAn1, double* fMn1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AnEnpEn_coef(fAn1, fMn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  add_AnEnpEp_coef(fAn1, fMn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_ApEfEf_coef(double* fAp1, double* fMp1,
			    double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn2, double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_ApEfEn_coef(fAp1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_ApEfEp_coef(fAp1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_AnEfEf_coef(double* fAn1, double* fMn1,
			    double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn2, double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AnEfEn_coef(fAn1, fMn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_AnEfEp_coef(fAn1, fMn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_AfEpEf_coef(double* fAn1, double* fAp1,
			    double* fMn1, double* fMp1, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AfEpEn_coef(fAn1, fAp1, fMn1, fMp1, fEp2, fEn3,        tau1, dm, tau2, tau3, weight2);
  add_AfEpEp_coef(fAn1, fAp1, fMn1, fMp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_AfEnEf_coef(double* fAn1, double* fAp1,
			    double* fMn1, double* fMp1, double* fEn2, double* fEn3, double* fEp3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AfEnEn_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_AfEnEp_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp3,        tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_AfEnpEf_coef(double* fAn1, double* fAp1,
			     double* fMn1, double* fMp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p,
			     const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AfEnpEn_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  add_AfEnpEp_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_AfEfEf_coef(double* fAn1, double* fAp1,
			    double* fMn1, double* fMp1,
			    double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn2, double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_AfEfEn_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_AfEfEp_coef(fAn1, fAp1, fMn1, fMp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MpEnpEf_coef(double* fMp1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MpEnpEn_coef(fMp1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  add_MpEnpEp_coef(fMp1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_MnEnpEf_coef(double* fMn1, double* fAn1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MnEnpEn_coef(fMn1, fAn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  add_MnEnpEp_coef(fMn1, fAn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_MpEfEf_coef(double* fMp1, double* fAp1,
			    double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn2, double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MpEfEn_coef(fMp1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_MpEfEp_coef(fMp1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MnEfEf_coef(double* fMn1, double* fAn1,
			    double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn2, double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MnEfEn_coef(fMn1, fAn1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_MnEfEp_coef(fMn1, fAn1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MfEpEf_coef(double* fMn1, double* fMp1,
			    double* fAn1, double* fAp1, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MfEpEn_coef(fMn1, fMp1, fAn1, fAp1, fEp2, fEn3,        tau1, dm, tau2, tau3, weight2);
  add_MfEpEp_coef(fMn1, fMp1, fAn1, fAp1, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MfEnEf_coef(double* fMn1, double* fMp1,
			    double* fAn1, double* fAp1, double* fEn2, double* fEn3, double* fEp3,
			    double* fxEn2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MfEnEn_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_MfEnEp_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp3,        tau1, dm, tau2, tau3, weight2);
  return;
}

inline void add_MfEnpEf_coef(double* fMn1, double* fMp1,
			     double* fAn1, double* fAp1,
			     double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			     double* fxEn2, double* fxEp2,
			     const double tau1, const double dm,
			     const double tau2_n, const double tau2_p, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MfEnpEn_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  add_MfEnpEp_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2_n, tau2_p, tau3, weight2);
  return;
}

inline void add_MfEfEf_coef(double* fMn1, double* fMp1,
			    double* fAn1, double* fAp1,
			    double* fEn2, double* fEp2, double* fEn3, double* fEp3,
			    double* fxEn2, double* fxEp2,
			    const double tau1, const double dm,
			    const double tau2, const double tau3, const double weight = 1.0){
  const double weight2 = 0.5*weight;

  add_MfEfEn_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEn3, fxEn2, tau1, dm, tau2, tau3, weight2);
  add_MfEfEp_coef(fMn1, fMp1, fAn1, fAp1, fEn2, fEp2, fEp3, fxEp2, tau1, dm, tau2, tau3, weight2);
  return;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __CONV_COEF3_H__ */

