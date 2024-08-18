//               
//    test_tatami.cc   
//       ---- Test program using tatami
//    $Id: test_tatami.cc 9948 2006-12-05 06:11:00Z katayama $ 
//                
//    $Log$
//    Revision 1.4  2004/10/19 07:09:54  kohji
//    1) SVD1 and SVD2 resolution parameters are updated for 2004 publication.
//    2) "expno" is added to add_outlier and dtres_systematics functions to treat
//    difference between SVD1 and SVD2 properly[cpfit_ml:0771].
//    3) New functions are added for cosh/sinh terms
//     (norm_)HAfRkRdetRnp_{fullrec,partial} for HAf=exp(-|t|Gamma)sinh(DGamma/2 t)
//     (norm_)HMfRkRdetRnp_{fullrec,partial} for HMf=exp(-|t|Gamma)cosh(DGamma/2 t)
//
//    Revision 1.3  2004/04/20 03:18:32  katayama
//    For 64bit machines
//
//    Revision 1.2  2002/10/02 00:40:53  katayama
//    use std::
//
//    Revision 1.1  2002/09/05 01:27:02  katayama
//    New package tatami from Nakadaira/Sumisawa san
//
//    Revision 1.9  2002/07/05 06:36:08  nakadair
//    Update for GCC3.
//
//    Revision 1.8  2002/07/05 06:17:35  nakadair
//    Add new checks.
//    Change to use CalcAkCk(...) instead of calc_akck(...).
//
//    Revision 1.7  2002/06/08 14:59:49  nakadair
//    Add code to check normalization and floating exception on ix86_linux.
//
//    Revision 1.6  2002/05/14 04:11:47  nakadair
//    Update for GCC3.
//
//    Revision 1.5  2002/05/11 07:51:42  nakadair
//    Update for solaris and GCC3.
//
//    Revision 1.4  2002/05/10 13:16:45  nakadair
//    Change output conversion specification of double variable.
//
//    Revision 1.3  2002/05/09 07:08:28  nakadair
//    Add add_outlier(...).
//
//    Revision 1.2  2002/05/09 06:26:11  nakadair
//    Add PDF examples.
//
//    Revision 1.1  2002/05/07 21:38:53  nakadair
//    Add test programs.
// 
//                
#include "belle.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "tatami/tatami.h"

#if (defined(__i386__)||defined(__x86_64__))&&defined(__linux__)
#include <ieee754.h>
#include <fpu_control.h>
#endif /* __ix86_linux__*/

#if defined(__sparc__)&&defined(__sun__)
#ifndef __GNUC__
#include <sunmath.h>
#endif
#endif /* __sparc_sun__ */

#include "belleutil/debugout.h"

using namespace Belle;

inline int
set_abort_by_sigfpe(){
//   signal(SIGFPE, crash_at_sigfpe);

#if (defined(__i386__)||defined(__x86_64__))&&defined(__linux__)
  unsigned int fpuflag = 0;
  _FPU_GETCW(fpuflag);
  dout(Debugout::INFO,"test_tatami") << "[set_abort_by_sigfpe] Chanege anction of SIGFPE" << std::endl;
  dout(Debugout::INFO,"test_tatami") << "FPU FLAG : " << std::hex << fpuflag ;
  fpuflag = 0x1372;
  _FPU_SETCW(fpuflag);
  _FPU_GETCW(fpuflag);
  dout(Debugout::INFO,"test_tatami") << " ==> " << std::hex << fpuflag << std::endl;
#endif /* __ix86_linux__*/
  return 0;
}

const char* test_tatami_id = "$Id: test_tatami.cc 9948 2006-12-05 06:11:00Z katayama $";

typedef double (*Xi_t)(const double t, const double tau, const double xd);

typedef double (*Ei_t)(const double t, const double tau);

extern "C" int setresp_(int* expmc);

inline double Xi_conv_Rk_full_by_int(Xi_t p_func,
				     const double t,
				     const double tau, const double xd,
				     const double ak, const double ck,
				     const double t_ini = -25.0,
				     const double t_end = +25.0,
				     const int ndiv = 5000){

  if(ck==0.0){
    return (*p_func)(t/ak, tau, xd);
  }
  double area=0.0;
  /* integral (-20sigma -- +20 sigma)*/
  const double t_del = (t_end-t_ini)/ndiv;
  for(int i=0;i<ndiv;++i){
    double t_run = t_ini+t_del*i;
    if(ck>0.0){
      area += (*p_func)(t_run, tau, xd)
	* Ep(t-t_run-((ak-1)*t_run+ck*std::abs(t_run)), ck*tau);
    }else{
      area += (*p_func)(t_run, tau, xd)
	* En(t-t_run-((ak-1)*t_run+ck*std::abs(t_run)), -ck*tau);
    }
  }
  if(ck>0.0){
    area -= 0.5*(*p_func)(t_ini, tau, xd)
      * Ep(t-t_ini-((ak-1)*t_ini+ck*std::abs(t_ini)), ck*tau);
    area -= 0.5*(*p_func)(t_end, tau, xd)
      * Ep(t-t_end-((ak-1)*t_end+ck*std::abs(t_end)), ck*tau);
  }else{
    area -= 0.5*(*p_func)(t_ini, tau, xd)
      * En(t-t_ini-((ak-1)*t_ini+ck*std::abs(t_ini)), -ck*tau);
    area -= 0.5*(*p_func)(t_end, tau, xd)
      * En(t-t_end-((ak-1)*t_end+ck*std::abs(t_end)), -ck*tau);
  }
  area *= t_del;
  return area;
}

inline double Xi_conv_Rk_full_by_int(Ei_t p_func,
				     const double t,
				     const double tau,
				     const double ak, const double ck,
				     const double t_ini = -25.0,
				     const double t_end = +25.0,
				     const int ndiv = 5000){
  if(ck==0.0){
    return (*p_func)(t/ak, tau);
  }
  double area=0.0;
  /* integral (-20sigma -- +20 sigma)*/
  const double t_del = (t_end-t_ini)/ndiv;
  for(int i=0;i<ndiv;++i){
    double t_run = t_ini+t_del*i;
    if(ck>0.0){
      area += (*p_func)(t_run, tau)
	* Ep(t-t_run-((ak-1)*t_run+ck*std::abs(t_run)), ck*tau);
    }else{
      area += (*p_func)(t_run, tau)
	* En(t-t_run-((ak-1)*t_run+ck*std::abs(t_run)), -ck*tau);
    }
  }
  if(ck>0.0){
    area -= 0.5*(*p_func)(t_ini, tau)
      * Ep(t-t_ini-((ak-1)*t_ini+ck*std::abs(t_ini)), ck*tau);
    area -= 0.5*(*p_func)(t_end, tau)
      * Ep(t-t_end-((ak-1)*t_end+ck*std::abs(t_end)), ck*tau);
  }else{
    area -= 0.5*(*p_func)(t_ini, tau)
      * En(t-t_ini-((ak-1)*t_ini+ck*std::abs(t_ini)), -ck*tau);
    area -= 0.5*(*p_func)(t_end, tau)
      * En(t-t_end-((ak-1)*t_end+ck*std::abs(t_end)), -ck*tau);
  }
  area *= t_del;
  return area;
}

extern "C" {
  extern dtres_param_t respar_;
}

#define M_SQRT2PI 2.50662827463100050240
#define M_SQRT1_PI 0.56418958354775628695
#define M_SQRTPI_2 1.25331413731550025120

const double g(const double x, const double mean, const double sigma)
{
  return std::exp(-0.5*(x-mean)*(x-mean)/(sigma*sigma))/(M_SQRT2PI*sigma);
}

const double gep(const double x, const double tau,
		 const double mean, const double sigma, const double offset)
{
  double x2 = x + offset;
  double a = M_SQRT1_2*(sigma/tau - (x2-mean)/sigma);
  double b = 0.5*sigma*sigma/(tau*tau) - (x2-mean)/tau;

  return (a >= 26 ? 0.5*M_SQRT1_PI*std::exp(b-a*a)*(1.0 - 0.5/(a*a))/(a*tau)
	  : 0.5*std::exp(b)*erfc(a)/tau);
}


const double gen(const double x, const double tau,
		 const double mean, const double sigma, const double offset)
{
  double x2 = x + offset;
  double a = M_SQRT1_2*(sigma/tau + (x2-mean)/sigma);
  double b = 0.5*sigma*sigma/(tau*tau) + (x2-mean)/tau;

  return (a >= 26 ? 0.5*M_SQRT1_PI*std::exp(b-a*a)*(1.0 - 0.5/(a*a))/(a*tau)
	  : 0.5*std::exp(b)*erfc(a)/tau);
}



const double gmp(const double x, const double tau, const double dm,
		 const double mean, const double sigma, const double offset)
{
  double x2, xx, yy;

  x2 = x + offset;
  xx = M_SQRT1_2*sigma*dm;
  yy = M_SQRT1_2*(sigma/tau-(x2-mean)/sigma);

  if(yy >= 0) return M_SQRTPI_2*sigma*g(x2, mean, sigma)*rewerf(xx, yy);
  else 
    return std::cos(2*xx*yy)*std::exp(yy*yy-xx*xx-0.5*(x2-mean)*(x2-mean)/(sigma*sigma))
      -M_SQRTPI_2*sigma*g(x2, mean, sigma)*rewerf(xx, -yy);
}


const double gmn(const double x, const double tau, const double dm,
		 const double mean, const double sigma, const double offset)
{
  double x2, xx, yy;

  x2 = x + offset;
  xx = M_SQRT1_2*sigma*dm;
  yy = M_SQRT1_2*(sigma/tau+(x2-mean)/sigma);

  if(yy >= 0) return M_SQRTPI_2*sigma*g(x2, mean, sigma)*rewerf(xx, yy);
  else
    return std::cos(2.0*xx*yy)
      *std::exp(yy*yy-xx*xx-0.5*(x2-mean)*(x2-mean)/(sigma*sigma))
      -M_SQRTPI_2*sigma*g(x2, mean, sigma)*rewerf(xx, -yy);
}

const double gap(const double x, const double tau, const double dm,
		 const double mean, const double sigma, const double offset)
{
  double x2, xx, yy;

  x2 = x + offset;
  xx = M_SQRT1_2*sigma*dm;
  yy = M_SQRT1_2*(sigma/tau-(x2-mean)/sigma);

  if(yy >= 0) return M_SQRTPI_2*sigma*g(x2, mean, sigma)*imwerf(xx, yy);
  else
    return M_SQRTPI_2*sigma*g(x2, mean, sigma)*imwerf(xx, -yy)
      -std::sin(2.0*xx*yy)*std::exp(yy*yy-xx*xx-0.5*(x2-mean)*(x2-mean)/(sigma*sigma));
}

const double gan(const double x, const double tau, const double dm,
		 const double mean, const double sigma, const double offset)
{
  double x2, xx, yy;

  x2 = x + offset;
  xx = M_SQRT1_2*sigma*dm;
  yy = M_SQRT1_2*(sigma/tau+(x2-mean)/sigma);

  if(yy >= 0) return -M_SQRTPI_2*sigma*g(x2, mean, sigma)*imwerf(xx, yy);
  else
    return -M_SQRTPI_2*sigma*g(x2, mean, sigma)*imwerf(xx, -yy)
      +std::sin(2.0*xx*yy)*std::exp(yy*yy-xx*xx-0.5*(x2-mean)*(x2-mean)/(sigma*sigma));
}

typedef double (*Ei_cg_t)(const double t, const double tau, const double m, const double s);

typedef double (*Xi_cg_t)(const double t, const double tau, const double xd, const double m, const double s);


inline double integral(Xi_cg_t func, const double tau, const double xd,
		       const double m, const double s,
		       const double b, const double u, const double dx){
  if(u==b)
    return 0.0;

  double area = 0.0;
  const int ndiv = (int)(double)(std::abs(u-b)/dx);
  const double t_ini = (u<b) ? u : b;
  for(int i=0;i<ndiv;++i){
    const double t_run = t_ini+dx*i;
    area += (*func)(t_run, tau, xd, m, s);
  }
  area -= (*func)(u, tau, xd, m, s)*0.5;
  area -= (*func)(b, tau, xd, m, s)*0.5;
  area *= dx;
if(u<b)
    area *= -1.0;
  return area;
}

inline double integral(Xi_cg_t func, const double tau, const double xd,
		       const double m, const double s,
		       const double b, const double u, const int n){
  const double dx = std::abs(u-b)/(double)n;
  return integral(func, tau, xd, m, s, b, u, dx);
};

inline double integral(Ei_cg_t func, const double tau,
		       const double m, const double s,
		       const double b, const double u, const double dx){
  if(u==b)
    return 0.0;

  double area = 0.0;
  const int ndiv = (int)(std::abs(u-b)/dx);
  const double t_ini = (u<b) ? u : b;
  for(int i=0;i<ndiv;++i){
    const double t_run = t_ini+dx*i;
    area += (*func)(t_run, tau, m, s);
  }
  area -= (*func)(u, tau, m, s)*0.5;
  area -= (*func)(b, tau, m, s)*0.5;
  area *= dx;
if(u<b)
    area *= -1.0;
  return area;
}

inline double integral(Ei_cg_t func, const double tau,
		       const double m, const double s,
		       const double b, const double u, const int n){
  const double dx = std::abs(u-b)/(double)n;
  return integral(func, tau, m, s, b, u, dx);
};

typedef double (*nEi_cg_t)(const double l, const double u, const double tau, const double m, const double s, const double o);

typedef double (*nXi_cg_t)(const double l, const double u, const double tau, const double xd, const double m, const double s, const double o);


void check_norm(Ei_cg_t func, nEi_cg_t func_norm){
  const double tau = 1.5;
  const double m   = 0.1;
  const double s   = 1.0;

  const double int_grad  = 0.001;
  
  const double l_lmt = -200.0;
  const double u_lmt =  200.0;
  const double grad  = 1.0;
  const int nl       = (int)((u_lmt-l_lmt)/grad);
  for(int ii=0;ii<nl;++ii){ 
    const double xl = l_lmt + grad * ii;
    const double lx = (xl>0) ? l_lmt : xl;
    const double ux = (xl>0) ? xl    : u_lmt;
    const double Int1 = integral(func, tau, m, s, lx, ux, int_grad);
    const double Int2 = func_norm(lx, ux, tau, m, s, 0.0);

    const double dInt = Int1-Int2;
    const double aInt = (Int1+Int2)/2.0;

    const double rdInt = (aInt==0.0) ? 0.0 : (dInt/aInt);
    
    if(debugout("INFO")) std::printf("% e % e % e % e % e\n",
				     xl, Int1, Int2, dInt, rdInt);
  }
  return;
}

void check_norm(Xi_cg_t func, nXi_cg_t func_norm){
  const double tau = 1.5;
  const double xd  = 0.425;
  const double m   = 0.1;
  const double s   = 1.0;

  const double int_grad  = 0.001;
  
  const double l_lmt = -200.0;
  const double u_lmt =  200.0;
  const double grad  = 1.0;
  const int nl       = (int)((u_lmt-l_lmt)/grad);
  for(int ii=0;ii<nl;++ii){ 
    const double xl = l_lmt + grad * ii;
    const double lx = (xl>0) ? l_lmt : xl;
    const double ux = (xl>0) ? xl    : u_lmt;
    const double Int1 = integral(func, tau, xd, m, s, lx, ux, int_grad);
    const double Int2 = func_norm(lx, ux, tau, xd, m, s, 0.0);

    const double dInt = Int1-Int2;
    const double aInt = (Int1+Int2)/2.0;

    const double rdInt = (aInt==0.0) ? 0.0 : (dInt/aInt);
    
    if(debugout("INFO")) std::printf("% e % e % e % e % e\n",
				     xl, Int1, Int2, dInt, rdInt);
  }
  return;
}

int main(int /* argc*/, char** /* argv */){
  set_abort_by_sigfpe();
//   check_norm(Ep_conv_gauss, norm_Ep_conv_gauss);
//   check_norm(En_conv_gauss, norm_En_conv_gauss);
//   check_norm(An_conv_gauss, norm_An_conv_gauss);
  check_norm(Ap_conv_gauss, norm_Ap_conv_gauss);
//   check_norm(Mn_conv_gauss, norm_Mn_conv_gauss);
//   check_norm(Mp_conv_gauss, norm_Mp_conv_gauss);
  exit(0);

#if 0
  /* Cross check with tomura-higuchi implementation */
  const double lll = -100.0;
  const double uuu =  100.0;
  const double dx  = 0.001;
  const int np  = (int)((uuu-lll)/dx);
  for(int ix=0;ix<np;++ix){
    const double x = lll+dx*ix;
//     double y1t = gep(x, 1.0, 0.0, 0.2, 0.0);
//     double y1m = nEp_conv_gauss(x, 0.0, 0.2);
//     double y2t = gen(x, 1.0, 0.0, 0.2, 0.0);
//     double y2m = nEn_conv_gauss(x, 0.0, 0.2);
#if 0
    double y1t = gmp(x, 1.0, 0.45, 0.0, 0.2, 0.0);
    double y1m = nMp_conv_gauss(x, 0.45, 0.0, 0.2);
    double y2t = gmn(x, 1.0, 0.45, 0.0, 0.2, 0.0);
    double y2m = nMn_conv_gauss(x, 0.45, 0.0, 0.2);

#else
    double y1t = gap(x, 1.0, 0.45, 0.0, 0.2, 0.0);
    double y1m = nAp_conv_gauss(x, 0.45, 0.0, 0.2);
    double y2t = gan(x, 1.0, 0.45, 0.0, 0.2, 0.0);
    double y2m = nAn_conv_gauss(x, 0.45, 0.0, 0.2);
#endif
    if(debugout("INFO")) std::printf("% e % e % e % e % e % e % e \n",
				     x, y1t, y1m, y2t, y2m,
				     (y1t+y1m==0.0) ? 0.0 : 2.0*(y1t-y1m)/(y1t+y1m),
				     (y2t+y2m==0.0) ? 0.0 : 2.0*(y2t-y2m)/(y2t+y2m));
  }
#endif

# if 0 
  /* Test of parameter exchange */
  /*  between C/C++ and fortran */
  int expmc = 1;
  int status = setresp_(&expmc);
  if(debugout("INFO")) std::printf("%d : % e  --- % e\n",
				   status, (double)respar_.Srec[0], dtres_param_default.Srec[0]);
  if(debugout("INFO")) std::printf("%d : % e  --- % e\n",
				   status, (double)respar_.tau_np_p_mlt[0][0],     dtres_param_default.tau_np_p_mlt[0][0]);
  if(debugout("INFO")) std::printf("%d : % e  --- % e\n",							 		   
				   status, (double)respar_.tau_np_p_mlt[0][1],     dtres_param_default.tau_np_p_mlt[0][1]);
  if(debugout("INFO")) std::printf("%d : % e  --- % e\n",							 		   
				   status, (double)respar_.tau_np_p_mlt[1][0],     dtres_param_default.tau_np_p_mlt[1][0]);
  if(debugout("INFO")) std::printf("%d : % e  --- % e\n",							 		   
				   status, (double)respar_.tau_np_p_mlt[1][1],     dtres_param_default.tau_np_p_mlt[1][1]);
  if(debugout("INFO")) std::printf("%d : % e  --- % e\n",
				   status, (double)respar_.fk2,     dtres_param_default.fk2);
  return 0;
#endif
  
  /* Test of Check XfRk_fullrec  */
  /* Comparison with numerical integral */
  const double llmt = -70.0;
  const double ulmt =  70.0;
  const double granularity = 0.1;
  const double nx = (ulmt-llmt)/granularity;


  const double ak = 0.98;
  const double ck = 0.5;

  const double dm  = 0.4;
  const double tau = 1.5; 

  const int bufsiz = 2048;
  char buf[bufsiz] = "";

  const char* fn1 = "mfrk.dat";
  std::ofstream fout1(fn1);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn1 << std::endl;
  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double y1 = Xi_conv_Rk_full_by_int(Mf, x, tau, dm, ak, ck);
    const double y2 = MfRk_fullrec(x, tau, dm, ak, ck);
   std::sprintf(buf, "% e % e % e % e", x, y1, y2, y2-y1);
    fout1 << buf << std::endl;
  }
  fout1.close();

  const char* fn2 = "afrk.dat";
  std::ofstream fout2(fn2);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn2 << std::endl;
  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double y1 = Xi_conv_Rk_full_by_int(Af, x, tau, dm, ak, ck);
    const double y2 = AfRk_fullrec(x, tau, dm, ak, ck);
   std::sprintf(buf, "% e % e % e % e", x, y1, y2, y2-y1);
    fout2 << buf << std::endl;
  }
  fout2.close();

  const char* fn3 = "efrk.dat";
  std::ofstream fout3(fn3);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn3 << std::endl;
  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double y1 = Xi_conv_Rk_full_by_int(Ef, x, tau, ak, ck);
    const double y2 = EfRk_fullrec(x, tau, ak, ck);
   std::sprintf(buf, "% e % e % e % e", x, y1, y2, y2-y1);
    fout3 << buf << std::endl;
  }
  fout3.close();

  /* Test of Check XfRk_partial  */
  /* Comparison with raw distribution */

  const char* fn4 = "efrkpar.dat";
  std::ofstream fout4(fn4);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn4 << std::endl;
  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double dz = 0.425*0.03*x;
    const double y1 = Ef(x, tau);
    const double y2 = EfRk_partial(x, tau, dz, &dtres_param_default);
   std::sprintf(buf, "% e % e % e % e", x, y1, y2, y2-y1);
    fout4 << buf << std::endl;
  }
  fout4.close();
  
  const char* fn5 = "afrkpar.dat";
  std::ofstream fout5(fn5);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn5 << std::endl;
  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double dz = 0.425*0.03*x;
    const double y1 = Af(x, tau, dm);
    const double y2 = AfRk_partial(x, tau, dm, dz, &dtres_param_default);
   std::sprintf(buf, "% e % e % e % e", x, y1, y2, y2-y1);
    fout5 << buf << std::endl;
  }
  fout5.close();
  
  const char* fn6 = "mfrkpar.dat";
  std::ofstream fout6(fn6);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn6 << std::endl;
  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double dz = 0.425*0.03*x;
#if 1
    const double y1 = Mf(x, tau, dm);
    const double y2 = MfRk_partial(x, tau, dm, dz, &dtres_param_default);
#else
    const double tau1 = 1.5;
    const double dm   = 0.45;
    const double tau2 = 0.3;
#if 0
    double fEn1=0.0, fEn2=0.0, fxEn1=0.0;
    add_EnEn_coef(&fEn1, &fEn2, &fxEn1, tau1, tau2);
    const double y1 = En_conv_En(x, 0.0, tau1, 0.0, tau2);
    const double y2 = fEn1*En(x, tau1)+fxEn1*xEn(x,tau1)+fEn2*En(x, tau2);
#endif
#if 0
    double fEp1=0.0, fEn2=0.0;
    add_EpEn_coef(&fEp1, &fEn2, tau1, tau2);
    const double y1 = Ep_conv_En(x, 0.0, tau1, 0.0, tau2);
    const double y2 = fEp1*Ep(x, tau1)+fEn2*En(x, tau2);
#endif

#if 0
    double fEn1=0.0, fEp2=0.0;
    add_EnEp_coef(&fEn1, &fEp2, tau1, tau2);
    const double y1 = En_conv_Ep(x, 0.0, tau1, 0.0, tau2);
    const double y2 = fEn1*En(x, tau1)+fEp2*Ep(x, tau2);
#endif
#if 1
    double fEp1=0.0, fEp2=0.0, fxEp1=0.0;
    add_EpEp_coef(&fEp1, &fEp2, &fxEp1, tau1, tau2);
    const double y1 = Ep_conv_Ep(x, 0.0, tau1, 0.0, tau2);
    const double y2 = fEp1*Ep(x, tau1)+fxEp1*xEp(x,tau1)+fEp2*Ep(x, tau2);
#endif

#if 1
    double fAp1=0.0, fMp1=0.0, fEn2=0.0;
    add_ApEn_coef(&fAp1, &fMp1,& fEn2, tau1, dm, tau2);
    const double y1 = Ap_conv_En(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fAp1*Ap(x, tau1, dm)+fMp1*Mp(x, tau1, dm)+fEn2*En(x, tau2);
#endif
#if 0
    double fAn1=0.0, fMn1=0.0, fEn2=0.0;
    add_AnEn_coef(&fAn1, &fMn1,& fEn2, tau1, dm, tau2);
    const double y1 = An_conv_En(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fAn1*An(x, tau1, dm)+fMn1*Mn(x, tau1, dm)+fEn2*En(x, tau2);
#endif
#if 0
    double fAp1=0.0, fMp1=0.0, fEp2=0.0;
    add_ApEp_coef(&fAp1, &fMp1,& fEp2, tau1, dm, tau2);
    const double y1 = Ap_conv_Ep(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fAp1*Ap(x, tau1, dm)+fMp1*Mp(x, tau1, dm)+fEp2*Ep(x, tau2);
#endif
#if 0
    double fAn1=0.0, fMn1=0.0, fEp2=0.0;
    add_AnEp_coef(&fAn1, &fMn1,& fEp2, tau1, dm, tau2);
    const double y1 = An_conv_Ep(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fAn1*An(x, tau1, dm)+fMn1*Mn(x, tau1, dm)+fEp2*Ep(x, tau2);
#endif
#if 0
    double fMp1=0.0, fAp1=0.0, fEn2=0.0;
    add_MpEn_coef(&fMp1, &fAp1,& fEn2, tau1, dm, tau2);
    const double y1 = Mp_conv_En(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fMp1*Mp(x, tau1, dm)+fAp1*Ap(x, tau1, dm)+fEn2*En(x, tau2);
#endif
#if 0
    double fMn1=0.0, fAn1=0.0, fEn2=0.0;
    add_MnEn_coef(&fMn1, &fAn1,& fEn2, tau1, dm, tau2);
    const double y1 = Mn_conv_En(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fMn1*Mn(x, tau1, dm)+fAn1*An(x, tau1, dm)+fEn2*En(x, tau2);
#endif
#if 0
    double fMp1=0.0, fAp1=0.0, fEp2=0.0;
    add_MpEp_coef(&fMp1, &fAp1,& fEp2, tau1, dm, tau2);
    const double y1 = Mp_conv_Ep(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fMp1*Mp(x, tau1, dm)+fAp1*Ap(x, tau1, dm)+fEp2*Ep(x, tau2);
#endif
#if 0
    double fMn1=0.0, fAn1=0.0, fEp2=0.0;
    add_MnEp_coef(&fMn1, &fAn1,& fEp2, tau1, dm, tau2);
    const double y1 = Mn_conv_Ep(x, 0.0, tau1, dm, 0.0, tau2);
    const double y2 = fMn1*Mn(x, tau1, dm)+fAn1*An(x, tau1, dm)+fEp2*Ep(x, tau2);
#endif
#endif
   std::sprintf(buf, "% e % e % e % e", x, y1, y2, y2-y1);
    fout6 << buf << std::endl;
  }
  fout6.close();
  
  /* CPfit PDG and so on */

  const double sin2phi1 = 0.8;
  const int bflavor  = 0; /* B0 := 0, B+ := 1; */
  const double w = 0.1;
  const double amix = 1.0 - 2.0 * w;
  const int ntrk_rec       = 2;
  const double sz_rec      = 0.008;
  const double chisq_z_rec = 1.200;
  const int ndf_z_rec      = 1;
  const int ntrk_asc       = 1;
  const double sz_asc      = 0.008;
  const double chisq_z_asc = 1.200;
  const int ndf_z_asc      = 2;
  const int keeptagl       = 1;
  const int expno          = 27;

  double a_k = 1.0, c_k = 0.0;
  const double cos_theta_b = 0.2;
  const double Eb_cms      = 5.29002;

  CalcAkCk(cos_theta_b, Eb_cms, &a_k, &c_k);

  const char* fn7 = "cpfitpdf.dat";
  std::ofstream fout7(fn7);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn7 << std::endl;

  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double L_exp = EfRkRdetRnp_fullrec(x, bflavor, tau, a_k, c_k,
					     ntrk_rec, sz_rec,
					     chisq_z_rec, ndf_z_rec,
					     ntrk_asc, sz_asc,
					     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default);

    const double L_sin = AfRkRdetRnp_fullrec(x, bflavor, tau, dm, a_k, c_k,
					     ntrk_rec, sz_rec,
					     chisq_z_rec, ndf_z_rec,
					     ntrk_asc, sz_asc,
					     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)/tau*0.5;
    const double y1 = add_outlier(expno, x, L_exp+amix*sin2phi1*L_sin, ntrk_rec, ntrk_asc, &dtres_param_default);
    const double y2 = add_outlier(expno, x, L_exp-amix*sin2phi1*L_sin, ntrk_rec, ntrk_asc, &dtres_param_default);
   std::sprintf(buf, "% e % e % e % e % e % e", x, y1, y2, y2-y1, (y2-y1)/(y2+y1),
	    L_sin);
    fout7 << buf << std::endl;
  }
  fout7.close();


  const char* fn8 = "dmfitpdf.dat";
  std::ofstream fout8(fn8);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn8 << std::endl;

  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double L_exp = EfRkRdetRnp_fullrec(x, bflavor, tau, a_k, c_k,
					     ntrk_rec, sz_rec,
					     chisq_z_rec, ndf_z_rec,
					     ntrk_asc, sz_asc,
					     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)*0.5;
    
    const double L_cos = MfRkRdetRnp_fullrec(x, bflavor, tau, dm, a_k, c_k,
					     ntrk_rec, sz_rec,
					     chisq_z_rec, ndf_z_rec,
					     ntrk_asc, sz_asc,
					     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)/tau*0.25;
    const double y1 = add_outlier(expno, x, L_exp+amix*L_cos, ntrk_rec, ntrk_asc, &dtres_param_default);
    const double y2 = add_outlier(expno, x, L_exp-amix*L_cos, ntrk_rec, ntrk_asc, &dtres_param_default);
   std::sprintf(buf, "% e % e % e % e % e", x, y1, y2, y2-y1, (y2-y1)/(y2+y1));
    fout8 << buf << std::endl;
  }
  fout8.close();


  const double Apipi = +0.7;
  const double Spipi = -0.7;


  const char* fn9 = "pipipdf.dat";
  std::ofstream fout9(fn9);
  dout(Debugout::ERR,"test_tatami") << "Output to file: " << fn9 << std::endl;

  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    memset(buf, '\0', bufsiz);
    const double L_exp = EfRkRdetRnp_fullrec(x, bflavor, tau, a_k, c_k,
					     ntrk_rec, sz_rec,
					     chisq_z_rec, ndf_z_rec,
					     ntrk_asc, sz_asc,
					     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)*0.5;

    const double L_sin = AfRkRdetRnp_fullrec(x, bflavor, tau, dm, a_k, c_k,
					     ntrk_rec, sz_rec,
					     chisq_z_rec, ndf_z_rec,
					     ntrk_asc, sz_asc,
					     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)/tau*0.25;

    const double L_cos = MfRkRdetRnp_fullrec(x, bflavor, tau, dm, a_k, c_k,
					     ntrk_rec, sz_rec,
					     chisq_z_rec, ndf_z_rec,
					     ntrk_asc, sz_asc,
					     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)/tau*0.25;
    const double y1 = add_outlier(expno, x, L_exp+amix*(Apipi*L_cos+Spipi*L_sin), ntrk_rec, ntrk_asc, &dtres_param_default);
    const double y2 = add_outlier(expno, x, L_exp-amix*(Apipi*L_cos+Spipi*L_sin), ntrk_rec, ntrk_asc, &dtres_param_default);
   std::sprintf(buf, "% e % e % e % e % e", x, y1, y2, y2-y1, (y2-y1)/(y2+y1));
    fout9 << buf << std::endl;
  }
  fout9.close();


  double area_EfRk        = 0.0;
  double area_EfRkRdet    = 0.0;
  double area_EfRkRdetRnp = 0.0;
  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
    area_EfRk += EfRk_fullrec(x, tau, ak, ck);

    area_EfRkRdet += EfRkRdet_fullrec(x, tau, a_k, c_k,
				      ntrk_rec, sz_rec,
				      chisq_z_rec, ndf_z_rec,
				      ntrk_asc, sz_asc,
				      chisq_z_asc, ndf_z_asc, &dtres_param_default);

    area_EfRkRdetRnp += EfRkRdetRnp_fullrec(x, bflavor, tau, a_k, c_k,
					    ntrk_rec, sz_rec,
					    chisq_z_rec, ndf_z_rec,
					    ntrk_asc, sz_asc,
					    chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default);
  }
  area_EfRk -= EfRk_fullrec(llmt, tau, ak, ck)*0.5;
  area_EfRk -= EfRk_fullrec(ulmt, tau, ak, ck)*0.5;
  area_EfRk *= granularity;
  

  area_EfRkRdet -= EfRkRdet_fullrec(llmt, tau, a_k, c_k,
				    ntrk_rec, sz_rec,
				    chisq_z_rec, ndf_z_rec,
				    ntrk_asc, sz_asc,
				    chisq_z_asc, ndf_z_asc, &dtres_param_default)*0.5;
  area_EfRkRdet -= EfRkRdet_fullrec(ulmt, tau, a_k, c_k,
				    ntrk_rec, sz_rec,
				    chisq_z_rec, ndf_z_rec,
				    ntrk_asc, sz_asc,
				    chisq_z_asc, ndf_z_asc, &dtres_param_default)*0.5;
  area_EfRkRdet *= granularity;

  area_EfRkRdetRnp -= EfRkRdetRnp_fullrec(llmt, bflavor, tau, a_k, c_k,
				     ntrk_rec, sz_rec,
				     chisq_z_rec, ndf_z_rec,
				     ntrk_asc, sz_asc,
				     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)*0.5;
  
  area_EfRkRdetRnp -= EfRkRdetRnp_fullrec(ulmt, bflavor, tau, a_k, c_k,
				     ntrk_rec, sz_rec,
				     chisq_z_rec, ndf_z_rec,
				     ntrk_asc, sz_asc,
				     chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)*0.5;
  area_EfRkRdetRnp *= granularity;
  dout(Debugout::INFO,"test_tatami") << "area_EfRk        : " << area_EfRk         << std::endl;
  dout(Debugout::INFO,"test_tatami") << "area_EfRkRdet    : " << area_EfRkRdet     << std::endl;
  dout(Debugout::INFO,"test_tatami") << "area_EfRkRdetRnp : " << area_EfRkRdetRnp  << std::endl;

  area_EfRk         = 0;
  area_EfRkRdet     = 0;
  area_EfRkRdetRnp  = 0;

  double area_Rk = 0.0;
  const double dz = 0.425*0.30*1.0;

  dtres_param_default.fk2 = 0.5;

  for(int i=0;i<nx;++i){
    double x = llmt + granularity*i;
//     const double dz = 0.425*0.03*x;
    area_Rk += Rk_partial(x, dz, &dtres_param_default);
    
    area_EfRk += EfRk_partial(x, tau, dz, &dtres_param_default);
    
    area_EfRkRdet += EfRkRdet_partial(x, tau, dz,
				      ntrk_rec, sz_rec,
				      chisq_z_rec, ndf_z_rec,
				      ntrk_asc, sz_asc,
				      chisq_z_asc, ndf_z_asc,
				      &dtres_param_default);

    area_EfRkRdetRnp += EfRkRdetRnp_partial(x, bflavor, tau, dz,
					    ntrk_rec, sz_rec,
					    chisq_z_rec, ndf_z_rec,
					    ntrk_asc, sz_asc,
					    chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default);
  }
//   const double ldz = 0.425*0.003*llmt;
//   const double udz = 0.425*0.003*ulmt;
  const double ldz = dz;
  const double udz = dz;

  area_Rk -= Rk_partial(llmt, ldz, &dtres_param_default)*0.5;
  area_Rk -= Rk_partial(ulmt, udz, &dtres_param_default)*0.5;
  area_Rk *= granularity;

  area_EfRk -= EfRk_partial(llmt, tau, ldz, &dtres_param_default)*0.5;
  area_EfRk -= EfRk_partial(ulmt, tau, udz, &dtres_param_default)*0.5;
  area_EfRk *= granularity;

  area_EfRkRdet -= EfRkRdet_partial(llmt, tau, ldz,
				    ntrk_rec, sz_rec,
				    chisq_z_rec, ndf_z_rec,
				    ntrk_asc, sz_asc,
				    chisq_z_asc, ndf_z_asc, &dtres_param_default)*0.5;
  area_EfRkRdet -= EfRkRdet_partial(ulmt, tau, udz,
				    ntrk_rec, sz_rec,
				    chisq_z_rec, ndf_z_rec,
				    ntrk_asc, sz_asc,
				    chisq_z_asc, ndf_z_asc, &dtres_param_default)*0.5;
  area_EfRkRdet *= granularity;

  area_EfRkRdetRnp -= EfRkRdetRnp_partial(llmt, bflavor, tau, ldz,
					  ntrk_rec, sz_rec,
					  chisq_z_rec, ndf_z_rec,
					  ntrk_asc, sz_asc,
					  chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)*0.5;
  
  area_EfRkRdetRnp -= EfRkRdetRnp_partial(ulmt, bflavor, tau, udz,
					  ntrk_rec, sz_rec,
					  chisq_z_rec, ndf_z_rec,
					  ntrk_asc, sz_asc,
					  chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default)*0.5;
  area_EfRkRdetRnp *= granularity;
  dout(Debugout::INFO,"test_tatami") << "area_Rk          : " << area_Rk           << std::endl;
  dout(Debugout::INFO,"test_tatami") << "area_EfRk        : " << area_EfRk         << std::endl;
  dout(Debugout::INFO,"test_tatami") << "area_EfRkRdet    : " << area_EfRkRdet     << std::endl;
  dout(Debugout::INFO,"test_tatami") << "area_EfRkRdetRnp : " << area_EfRkRdetRnp  << std::endl;

  const double nL_exp = norm_EfRkRdetRnp_fullrec(llmt, ulmt, bflavor, tau, a_k, c_k,
						 ntrk_rec, sz_rec,
						 chisq_z_rec, ndf_z_rec,
						 ntrk_asc, sz_asc,
						 chisq_z_asc, ndf_z_asc, keeptagl, &dtres_param_default);

  dout(Debugout::INFO,"test_tatami") << "nL_exp = " << nL_exp << std::endl;
  return 0;
}
