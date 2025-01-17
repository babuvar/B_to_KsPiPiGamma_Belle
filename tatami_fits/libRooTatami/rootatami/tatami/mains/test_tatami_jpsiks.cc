//
// $Id: test_tatami_jpsiks.cc 11208 2011-03-01 14:53:58Z hitoshi $
//
// test program using the package "tatami".
// extract sin2phi1 using J/psiKs mode (with ICHEP02 method).
//
// $Log$
// Revision 1.5  2004/10/19 07:09:54  kohji
// 1) SVD1 and SVD2 resolution parameters are updated for 2004 publication.
// 2) "expno" is added to add_outlier and dtres_systematics functions to treat
// difference between SVD1 and SVD2 properly[cpfit_ml:0771].
// 3) New functions are added for cosh/sinh terms
//  (norm_)HAfRkRdetRnp_{fullrec,partial} for HAf=exp(-|t|Gamma)sinh(DGamma/2 t)
//  (norm_)HMfRkRdetRnp_{fullrec,partial} for HMf=exp(-|t|Gamma)cosh(DGamma/2 t)
//
// Revision 1.4  2003/08/06 06:43:31  katayama
// New version from Nakadaira san
//
// Revision 1.3  2002/12/10 12:27:48  sumisawa
// Bug fix in (norm_)AfRkRdet(...), (norm_)MfRkRdet(...) and modify test_tatami_jpsiks.cc to speed up
//
// Revision 1.2  2002/09/05 05:24:35  sumisawa
// minor bug fix
//
// Revision 1.1  2002/09/05 01:27:03  katayama
// New package tatami from Nakadaira/Sumisawa san
//
//

#include "belle.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <cmath>
#include <sys/types.h>
#include <regex.h>

#include "tatami/tatami.h"
#include "belleutil/debugout.h"

#ifdef __sun
#include <ieeefp.h>
#endif

using namespace Belle;
//
// MINUIT interface
//

#define MINUIT_MINIMIZE (1)
#define MINUIT_MINOS    (2)

typedef void (*MNFUNC_PTR)(int*, double*, double*, double*, int*, void*);
typedef double (*mnfcn_t)(double*);

extern "C" {
  void mninit_(int*, int*, int*);
  void mnexcm_(MNFUNC_PTR, char*, double*, int*, int*, void*, int);
  void mnparm_(int*, char*, double*, double*, double*, double*, int*, int);
  void mnpout_(int*, char*, double*, double*, double*, double*, int*, int);
  void mnerrs_(int*, double*, double*, double*, double*);
}

typedef struct mnparam {
        /* inputs */
        char   mp_name[16];
        int    mp_fix;
        double mp_init;
        double mp_min, mp_max;
        double mp_acc;
        /* outputs */
        double mp_value;
        double mp_errneg, mp_errpos, mp_errpar;
} mnparam_t;

static mnparam_t MNpar[101];
static mnfcn_t MNfcn;
static int MNtouch = 1;

static void mnfcn(int *npar, double *grad, double *fval,
		  double *xval, int *iflag, void *futil)
{
  *fval = MNfcn ? (*MNfcn)(xval) : 0.0;
}

static void mnfunc(mnfcn_t fcn)
{
  MNtouch = 1;
  MNfcn = fcn;
}

static void mnreset(void)
{
  int i, one(1), ret;
  double opt;
  
  opt = 0.0;
  char command[] = "CLEAR";
  mnexcm_(mnfcn, command, &opt, &one, &ret, NULL, strlen(command));
  
  MNtouch = 1;
  MNfcn = NULL;
  for(i=0; i<100; i++) MNpar[i].mp_name[0] = '\0';
}

static void mninit(void)
{
  int ird(5), iwr(6), isav(7);
  double opt;
  int one = 1, ret;
  static int first = 1;
  
  MNtouch = 1;
  
  mnreset();
  
  if(first) mninit_(&ird, &iwr, &isav);
  first = 0;
  
  opt = -1.0;
  char command[] = "SET PRI";
  mnexcm_(mnfcn, command, &opt, &one, &ret, NULL, strlen(command));
}

static int mnpnum(void)
{
  int i(0);
  
  while( MNpar[i].mp_name[0] ) ++i;
  
  return i;
}

static void mnparam(char *name, double init, int fix,
		    double min, double max, double acc)
{
  int pnum = mnpnum();
  if( pnum>=100 ) return;
  
  mnparam_t &mnp(MNpar[pnum]);
  
  strncpy( mnp.mp_name, name, sizeof(mnp.mp_name) );
  mnp.mp_fix    = fix;
  mnp.mp_init   = init;
  mnp.mp_min    = min;
  mnp.mp_max    = max;
  mnp.mp_acc    = acc;
  
  mnp.mp_value  = init;
  mnp.mp_errneg = 0.0;
  mnp.mp_errpos = 0.0;
  mnp.mp_errpar = 0.0;

  MNtouch = 1;
}

static void mnfit(int which)
{
  int i, pnum=mnpnum(), num, ret;
  double opt[2];
  int zero(0), one(1), two (2);
  char command[100];
  
  if( MNtouch ){
    opt[0] = 1;
    strcpy(command, "SET PRI");
    mnexcm_(mnfcn, command, opt, &one, &ret, NULL, strlen(command));
    
    opt[0] = 0.0;
    strcpy(command, "SET WAR");
    mnexcm_(mnfcn, command, opt, &one, &ret, NULL, strlen(command));
    
    opt[0] = 1.0;
    strcpy(command, "SET STR");
    mnexcm_(mnfcn, command, opt, &one, &ret, NULL, strlen(command));
    
    for( i=0, num=1; i<pnum; i++, num++ ){
      int erflg;
      mnparm_(&num, MNpar[i].mp_name, &MNpar[i].mp_init, &MNpar[i].mp_acc,
	      &MNpar[i].mp_min, &MNpar[i].mp_max, &erflg,
	      strlen(MNpar[i].mp_name));
      if( MNpar[i].mp_fix ){
	opt[0] = (double)num;
	strcpy(command, "FIX");
	mnexcm_(mnfcn, command, opt, &one, &ret, NULL, strlen(command));
      }
    }
  }
  
  if( which==MINUIT_MINIMIZE ){
    strcpy(command, "MINI");
    mnexcm_(mnfcn, command, opt, &zero, &ret, NULL, strlen(command));
    mnexcm_(mnfcn, command, opt, &zero, &ret, NULL, strlen(command));
    strcpy(command, "HESS");
    mnexcm_(mnfcn, command, opt, &zero, &ret, NULL, strlen(command));
  }else{
    strcpy(command, "MINO");
    mnexcm_(mnfcn, command, opt, &zero, &ret, NULL, strlen(command));
  }

  for( i=0, num=1; i<pnum; i++, num++ ){
    char name[16];
    double min, max, globcc;
    int varbl;
    
    mnpout_(&num, name, &MNpar[i].mp_value, &MNpar[i].mp_errpar, &min, &max,
	    &varbl, strlen(name));
    mnerrs_(&num,
	    &MNpar[i].mp_errpos, &MNpar[i].mp_errneg, &MNpar[i].mp_errpar,
	    &globcc);
  }
  
  MNtouch = 0;
}


double mnvalue(int num, double *errneg, double *errpos, double *errpar)
{
  if( errneg ) *errneg = MNpar[num].mp_errneg;
  if( errpos ) *errpos = MNpar[num].mp_errpos;
  if( errpar ) *errpar = MNpar[num].mp_errpar;
  
  return MNpar[num].mp_value;
}

//------------------------------------------------------------
// calculate signal fraction

static const double argus (double x, double offset, 
			   double ebeam, double efact)
{
  if ( x > ebeam ) return 0.0;
  double dz1 = x + offset;
  double dz2 = 1.0 - ( dz1 / ebeam )* (dz1 / ebeam );
  double dz3 = efact*dz2;
  double argus = dz1 * sqrt(dz2) * exp(dz3 );
  return argus;
}

static const double int_argus(const double lower, const double upper,
			      const double offset,
			      const double e, const double a)
{
  double e2 = e*e;
  double ll2 = 1.0 - (lower+offset)*(lower+offset)/e2;
  double ul2 = 1.0 - (upper+offset)*(upper+offset)/e2;

  if(ll2 < 0.0) return 0.0;
  if(ul2 < 0.0) ul2 = 0.0;

  return 0.5*e2/a*((sqrt(ll2)*exp(a*ll2) - sqrt(ul2)*exp(a*ul2))
                   +0.5*sqrt(M_PI/fabs(a))
                   *(erfc(sqrt(fabs(a)*ll2)) - erfc(sqrt(fabs(a)*ul2))));
}


static const double sig_func(const int expno, const double mb, const double de,
			     const double ebeam, const int rbin)
{
  static const double de_min(-0.04);
  static const double de_max(0.04);
  static const double mb_min(5.2694);
  static const double mb_max(5.2894);


  double de_mean, de_sigma, de_fmain, de_mean2, de_sigma2,
    mb_mean, mb_sigma, argus_a, bg_de1;
  double sig_area[7];

  if(expno < 30){ //SVD1
    de_mean   = 6.93161e-04;
    de_sigma  = 6.44620e-03;
    de_fmain  = 6.57860e-01;
    de_mean2  = -9.97717e-04;
    de_sigma2 = 1.62821e-02;
    mb_mean   = 5.27951e+00;
    mb_sigma  = 2.48677e-03;
    argus_a   = -4.35127e+01;
    bg_de1    = -2.06461e+00;
    sig_area[0] = 0.955393;
    sig_area[1] = 0.973217;
    sig_area[2] = 0.964847;
    sig_area[3] = 0.978738;
    sig_area[4] = 0.979486;
    sig_area[5] = 0.986856;
    sig_area[6] = 0.988662;
  }else{ //SVD2
    de_mean   =  6.89758e-04;
    de_sigma  =  6.42571e-03;
    de_fmain  =  7.65624e-01;
    de_mean2  =  1.77623e-03;
    de_sigma2 =  2.03155e-02;
    mb_mean   =  5.27956e+00;
    mb_sigma  =  2.60387e-03;
    argus_a   = -4.40840e+01;
    bg_de1    = -2.48766e+00;
    sig_area[0] = 0.963777;
    sig_area[1] = 0.966783;
    sig_area[2] = 0.969343;
    sig_area[3] = 0.967370;
    sig_area[4] = 0.978108;
    sig_area[5] = 0.980252;
    sig_area[6] = 0.997117;
  }
  //signal
  const double int_sig =
    norm_gaussian(mb_min, mb_max, mb_mean, mb_sigma)
    * (de_fmain * norm_gaussian(de_min, de_max, de_mean, de_sigma)
       + (1-de_fmain) * norm_gaussian(de_min, de_max, de_mean2, de_sigma2));
  
  const double sig = sig_area[rbin] / int_sig
    * gaussian(mb, mb_mean, mb_sigma)
    * (de_fmain * gaussian(de, de_mean, de_sigma)
       + (1-de_fmain) * gaussian(de, de_mean2, de_sigma2));
  
  //background
  const double int_bg =
    int_argus(mb_min, mb_max, 0., ebeam, argus_a)
    * (1.0 + bg_de1*0.5*(de_min + de_max))*(de_max - de_min);

  const double bg = (1 - sig_area[rbin]) / int_bg
    * argus(mb, 0., ebeam, argus_a)
    * (1.0 + bg_de1*de);
  
  return sig/(sig + bg);
}

static const int set_rbin(const double r)
{
  return(0.   <=r && r<=0.1   ? 0 :
	 0.1  < r && r<=0.25  ? 1 :
	 0.25 < r && r<=0.5   ? 2 :
	 0.5  < r && r<=0.625 ? 3 :
	 0.625< r && r<=0.75  ? 4 :
	 0.75 < r && r<=0.875 ? 5 :
	 0.875< r && r<=1.0   ? 6 : 7);
}

static const double set_wtag(const int expno, const int rbin)
{
  static const double w_svd1[7]= {0.5, 0.418852, 0.329879, 0.233898, 0.170608,
				  0.099791,0.0228501};

  static const double w_svd2[7]= {0.5, 0.418826, 0.319303, 0.222948, 0.163191,
				  0.104085, 0.0251454};
  if(expno < 30){
    return w_svd1[rbin];
  }else{
    return w_svd2[rbin];
  }
}

static const double set_delta_wtag(const int expno, const int rbin)
{
  static const double dw_svd1[7] = {0., 0.0569661, 0.0126192,
				    -0.0147724, -0.000550289,
				    0.00887704, 0.00465683};

  static const double dw_svd2[7] = {0., -0.00877001, 0.0103515,
				    -0.0109253, -0.0186365,
				    0.00168037, -0.0036441};

  if(expno < 30){
    return dw_svd1[rbin];
  }else{
    return dw_svd2[rbin];
  }
}

//----------------------------------------------------

#define DECAY_PSI_KS            (1010)

#define Coeff_cm2ps   (dt_resol_global::inv_bgc)

struct vertex {
  double vt_mc;
  double vt_pos;
  double vt_time;
  double vt_err;
  int    vt_ntrk;
  double vt_chi2;
  int    vt_ndf;
  double vt_cl;
};

struct ICPV_data {
  // MC
  int    icpv_hepevt;
  double icpv_mc_delta_t;
  int    icpv_mc_flv;
  
  // event tags
  int    icpv_exp, icpv_run, icpv_evt, icpv_frm;
  int    icpv_decay;
  
  // vertex
  struct vertex icpv_rec;
  struct vertex icpv_asc;
  
  // kinematics
  double icpv_mbc;
  double icpv_deltae;
  double icpv_ecms;
  double icpv_ebeam;

  // flavor
  int    icpv_flv;
  int    icpv_rbin;
  int icpv_fbtg_mode;
  int icpv_keeptagl;
  double icpv_wtag, icpv_delta_wtag;
  double icpv_r;
  
  // CM system
  double icpv_mc_pbs;
  double icpv_pbs;
  double icpv_mc_costh;
  double icpv_costh;
  
  // for K*
  double icpv_costhtr;
  
  // signal/background fractions
  double icpv_sigfrac;

  //PDF
  mutable int need_calc;
  mutable double L_e, L_s, L_a, pdf_bg;
  mutable double int_L_e, int_L_s, int_L_a, int_pdf_bg;
};

static double FIT_S, FIT_A;
static std::vector<double *> MINUIT_par;

static std::vector<struct ICPV_data> ICPV_data;

static void read_param(void)
{
  const double min(-25.), max(25.);
  int fixed(0);

  FIT_S = 0.;
  double val(FIT_S);
  char name[] = "S";
  mnparam(name, val, fixed, min, max, 0.01);
  MINUIT_par.push_back(&FIT_S);

  FIT_A = 0.;
  val = FIT_A;
  strcpy(name, "A");
  mnparam(name, val, fixed, min, max, 0.01);
  MINUIT_par.push_back(&FIT_A);
}

static bool pattern_match(const char *buf)
{
  static regex_t preg;
  
  static bool first = true;
  if( first ){
    regcomp(&preg, "^[-+.0-9 \tEe]*\n$", 0);
    first = false;
  }
  if( regexec(&preg,buf,0,NULL,0)==REG_NOMATCH ){
    fprintf(stderr,"skipped bad event record (regexec): %s", buf);
    return false;
  }

  return true;
}

static bool valid_event(const struct ICPV_data &ev)
{
  
  if(abs(ev.icpv_hepevt)!=521 && abs(ev.icpv_hepevt)!=511)
    return false;

  if(ev.icpv_decay == DECAY_PSI_KS){
    if(fabs(ev.icpv_mbc - 5.2794) >= 0.01 ||
       fabs(ev.icpv_deltae) >= 0.04) return false;
  }else{
    return false;
  }

  if( ev.icpv_rbin==7 ) return false;
  //if( ev.icpv_flv==0 ) return false;
  
  const struct vertex &rec = ev.icpv_rec;
  const struct vertex &asc = ev.icpv_asc;
  
  if( rec.vt_err<=0 ) return false;
  if( asc.vt_err<=0 ) return false;
  
  if(rec.vt_time - asc.vt_time <= dt_resol_global::dt_llmt ) return false;
  if(rec.vt_time - asc.vt_time >= dt_resol_global::dt_ulmt ) return false;
 
  if(rec.vt_ndf > 0 && rec.vt_chi2/rec.vt_ndf >= 50) return false;
  if(asc.vt_ndf > 0 && asc.vt_chi2/asc.vt_ndf >= 50) return false;

  if(rec.vt_ntrk == 1 && rec.vt_err >= 0.05) return false;
  if(asc.vt_ntrk == 1 && asc.vt_err >= 0.05) return false;

  if(rec.vt_ntrk > 1 && rec.vt_err >= 0.02) return false;
  if(asc.vt_ntrk > 1 && asc.vt_err >= 0.02) return false;

  return true;
}

static void preset_event(struct ICPV_data &ev)
{
  struct vertex &rec = ev.icpv_rec;
  struct vertex &asc = ev.icpv_asc;
  
  rec.vt_time = rec.vt_pos * Coeff_cm2ps;
  asc.vt_time = asc.vt_pos * Coeff_cm2ps;

  const double wtag(ev.icpv_wtag);
  ev.icpv_r = 1 - 2*wtag;
  ev.icpv_rbin = set_rbin(ev.icpv_r);
  ev.icpv_wtag = set_wtag(ev.icpv_exp, ev.icpv_rbin);
  ev.icpv_delta_wtag = set_delta_wtag(ev.icpv_exp, ev.icpv_rbin);

  ev.icpv_sigfrac = sig_func(ev.icpv_exp, ev.icpv_mbc, ev.icpv_deltae,
			     ev.icpv_ebeam, ev.icpv_rbin);

  ev.icpv_ecms = sqrt(ev.icpv_pbs*ev.icpv_pbs +
		      dt_resol_global::mbzero*dt_resol_global::mbzero);

  ev.need_calc = 1;
}

static bool load_event(const char *buf, struct ICPV_data &ev)
{
  int dummy_i;
  double dummy_f;
  
  struct vertex &rec = ev.icpv_rec;
  struct vertex &asc = ev.icpv_asc;
  
  int ntoken;

  if( !pattern_match(buf) ) return false;
  
  ntoken = sscanf(
		  buf,
		  "%d %d %d %d "
		  "%d "
		  "%d "
		  "%le %le %le "
		  "%le %le "
		  "%d "
		  "%le %le %d %le %d %le "
		  "%le %le %d %le %d %le "
		  "%f "
		  "%d %d %le "
		  "%le %le %le %le %le %le\n",
		  &ev.icpv_exp, &ev.icpv_run, &ev.icpv_evt, &ev.icpv_frm,
		  &ev.icpv_hepevt,
		  &ev.icpv_mc_flv,
		  &rec.vt_mc, &asc.vt_mc, &ev.icpv_mc_delta_t,
		  &ev.icpv_mc_pbs, &ev.icpv_mc_costh,
		  &ev.icpv_decay,
		  &rec.vt_pos, &rec.vt_err, &rec.vt_ntrk, &rec.vt_chi2, &rec.vt_ndf, &rec.vt_cl,
		  &asc.vt_pos, &asc.vt_err, &asc.vt_ntrk, &asc.vt_chi2, &asc.vt_ndf, &asc.vt_cl,
		  &dummy_f,
		  //&ev.icpv_flv, &ev.icpv_fbtg_mode, &ev.icpv_wtag,
		  &ev.icpv_flv, &ev.icpv_keeptagl, &ev.icpv_wtag,
		  &ev.icpv_deltae, &ev.icpv_mbc, &ev.icpv_pbs, &ev.icpv_costh,
		  &ev.icpv_ebeam, &ev.icpv_costhtr);

  if( ntoken!=34 ) return false;
  preset_event(ev);    
  return valid_event(ev);
}

static void read_data(const char *file)
{
  int c=0;
  FILE *fp = fopen(file,"r");
  if( !fp ){
    perror("fopen");
    exit(1);
  }
  
  char buf[1024];
  while( fgets(buf,sizeof(buf),fp) ){
    struct ICPV_data ev;
    
    if( !load_event(buf,ev) ) continue;
    
    ICPV_data.push_back(ev);
  }
  
  fclose(fp);
  
  printf("number of fitted events: %d\n", ICPV_data.size());
}

static void load_fits(const double *par)
{
  for(unsigned int i=0; i<MINUIT_par.size(); ++i)
    *MINUIT_par[i] = par[i];
}

static double PDF_bg(const int expno,
		     const struct vertex &rec, const struct vertex &asc,
		     double &pdf, double &int_pdf)
{
  pdf = 0;
  int_pdf = 0;

  double s_main, s_tail, f_tail, tau_bg, f_delta, mu_l, mu_d;
  
  tau_bg = 1.4969;
  mu_l = -0.19137;
  mu_d = 0.72087E-02;
  
  if(rec.vt_ntrk > 1 && asc.vt_ntrk > 1){
   //multi track vertexing
    s_main  =  1.0920;
    s_tail  =  2.0939;
    f_tail  = 0.48098;
    f_delta = 0.61036;
  }else{
    s_main  = 1.1399;
    s_tail  = 23.106;
    f_tail  = 0.28176E-01;
    f_delta = 0.46589;
  }

  const double sigma_main(s_main
			  * sqrt(rec.vt_err*rec.vt_err + asc.vt_err*asc.vt_err)
			  * Coeff_cm2ps);
  const double sigma_tail(s_tail * sigma_main);

  const double dt_llmt(dt_resol_global::dt_llmt);
  const double dt_ulmt(dt_resol_global::dt_ulmt);

  const double delt(rec.vt_time - asc.vt_time);
  
  // lifetime compronent
  if(f_delta < 1){
    double pdf_l(0), int_pdf_l(0.);
    if(f_tail < 1){
      pdf_l += (1 - f_tail) *
	Enp_conv_gauss(delt, tau_bg, tau_bg, mu_l, sigma_main);
      int_pdf_l += (1 - f_tail) *
	norm_Enp_conv_gauss(dt_llmt, dt_ulmt,tau_bg,tau_bg,mu_l,sigma_main);
    }
    if(f_tail > 0){
      pdf_l += f_tail *
	Enp_conv_gauss(delt, tau_bg, tau_bg, mu_l, sigma_tail);
      int_pdf_l += f_tail *
	norm_Enp_conv_gauss(dt_llmt, dt_ulmt,tau_bg,tau_bg,mu_l,sigma_tail);
    }
    
    pdf += pdf_l * (1 - f_delta);
    int_pdf += int_pdf_l * (1 - f_delta);
  }

  // delta function
  if(f_delta > 0){
    double pdf_d(0), int_pdf_d(0.);
    if(f_tail < 1){
      pdf_d += (1 - f_tail) * gaussian(delt, mu_d, sigma_main);
      int_pdf_d += (1 - f_tail) *
	norm_gaussian(dt_llmt, dt_ulmt,mu_d,sigma_main);
    }
    if(f_tail > 0){
      pdf_d += f_tail * gaussian(delt, mu_d, sigma_tail);
      int_pdf_d += f_tail *
	norm_gaussian(dt_llmt, dt_ulmt,mu_d,sigma_tail);
    }
    
    pdf += pdf_d * f_delta;
    int_pdf += int_pdf_d * f_delta;
  }

  return (int_pdf > 0 ? pdf/int_pdf : 0);
}

static const double PDF_sin_cos(const struct ICPV_data &ev)
{
  static const double tau_b(1.525);
  static const double m_d(0.507);

  const struct vertex &rec = ev.icpv_rec;
  const struct vertex &asc = ev.icpv_asc;
  
  const double delt(rec.vt_time - asc.vt_time);

  const double dt_llmt(dt_resol_global::dt_llmt);
  const double dt_ulmt(dt_resol_global::dt_ulmt);

  // signal PDF

  double a_k, c_k;
  const double amix(ev.icpv_flv * (1 - 2*ev.icpv_wtag));
  const double cexp(1. - ev.icpv_flv * ev.icpv_delta_wtag);
  CalcAkCk(ev.icpv_costh, ev.icpv_ecms, &a_k, &c_k, 0);

  const int mc = 0 /* 1 for MC */;
  const dtres_param_t* dtres_param = get_dtres_param(ev.icpv_exp, mc);

  if(ev.need_calc){
    ev.L_e = 
      EfRkRdetRnp_fullrec(delt, 0, tau_b, a_k, c_k,
			  rec.vt_ntrk, rec.vt_err, rec.vt_chi2, rec.vt_ndf,
			  asc.vt_ntrk, asc.vt_err, asc.vt_chi2, asc.vt_ndf,
			  ev.icpv_keeptagl,
			  dtres_param);   
    ev.L_s =
      AfRkRdetRnp_fullrec(delt, 0, tau_b, m_d, a_k, c_k,
			  rec.vt_ntrk, rec.vt_err, rec.vt_chi2, rec.vt_ndf,
			  asc.vt_ntrk, asc.vt_err, asc.vt_chi2, asc.vt_ndf,
			  ev.icpv_keeptagl,
			  dtres_param) * 0.5 / tau_b;   

    ev.L_a =
      MfRkRdetRnp_fullrec(delt, 0, tau_b, m_d, a_k, c_k,
			  rec.vt_ntrk, rec.vt_err, rec.vt_chi2, rec.vt_ndf,
			  asc.vt_ntrk, asc.vt_err, asc.vt_chi2, asc.vt_ndf,
			  ev.icpv_keeptagl,
			  dtres_param) * 0.5 / tau_b;   

    ev.int_L_e =
      norm_EfRkRdetRnp_fullrec(dt_llmt, dt_ulmt, 0, tau_b, a_k, c_k,
			       rec.vt_ntrk, rec.vt_err, rec.vt_chi2, rec.vt_ndf,
			       asc.vt_ntrk, asc.vt_err, asc.vt_chi2, asc.vt_ndf,
			       ev.icpv_keeptagl,
			       dtres_param);   
//    ev.int_L_s =
//      norm_AfRkRdetRnp_fullrec(dt_llmt, dt_ulmt, 0, tau_b, m_d, a_k, c_k,
//			       rec.vt_ntrk, rec.vt_err, rec.vt_chi2, rec.vt_ndf,
//			       asc.vt_ntrk, asc.vt_err, asc.vt_chi2, asc.vt_ndf,
//			       dtres_param) * 0.5 / tau_b;   
  }
  const double pdf_sig(cexp * ev.L_e + amix * (FIT_S*ev.L_s + FIT_A*ev.L_a));
//  const double int_pdf_sig(cexp * ev.int_L_e +
//			   amix * (FIT_sin2phi1*ev.int_L_s));
  const double int_pdf_sig(2. * ev.int_L_e);

  //BG PDF
  if(ev.need_calc){
    PDF_bg(ev.icpv_exp, rec, asc, ev.pdf_bg, ev.int_pdf_bg);
    ev.int_pdf_bg *= 2.;
  }

  const double pdf(AddOutlierWithBkg(ev.icpv_exp,
				     delt, ev.icpv_sigfrac,
				     pdf_sig, ev.pdf_bg,
				     rec.vt_ntrk, asc.vt_ntrk,
				     dtres_param,
				     int_pdf_sig, ev.int_pdf_bg,
				     dt_llmt, dt_ulmt, 0.5));
//  const double pdf(pdf_sig/int_pdf_sig);

  ev.need_calc = 0;

  return pdf;
}

static const double PDF(const struct ICPV_data &ev)
{
  if(ev.icpv_decay == 1010)
    return PDF_sin_cos(ev);
  else
    return 0;
}

static double func_belle(double *par) // this is FCN
{
  load_fits(par);
  
  double sum_log_pdf(0.0);
  
  for(std::vector<struct ICPV_data>::const_iterator it=ICPV_data.begin();
      it!=ICPV_data.end(); ++it){
    double pdf = PDF(*it);
    sum_log_pdf += log(pdf);
    if( !finite(pdf)||pdf<=0 ){
      printf("%d %d %d %d: pdf ====> %le\n",
	     (*it).icpv_exp,
	     (*it).icpv_run,
	     (*it).icpv_evt,
	     (*it).icpv_frm, pdf);
      exit(7);
    }
  }
  
  if( !finite(sum_log_pdf) ){
    printf("sum_log_pdf = %lf\n",sum_log_pdf);
    exit(9);
  }
  
  return -2.0*sum_log_pdf; // to be minimized
}


void plot_data(const char *file)
{
  FILE *fp = fopen(file,"w");

  for(std::vector<struct ICPV_data>::const_iterator it=ICPV_data.begin();
      it!=ICPV_data.end(); ++it){
    const struct ICPV_data &ev(*it);

    const double delt = ev.icpv_rec.vt_time - ev.icpv_asc.vt_time;

    fprintf(fp,"%d %+6.4le\n", ev.icpv_flv, delt);
  }

  fclose(fp);
}

void plot_fit(const char *file)
{
  const double dt_lim(20);
  const unsigned nbin(200);

  std::vector<double> tot_p(nbin), tot_n(nbin);
  for(unsigned i=0; i<nbin; ++i){
    tot_p[i] = 0.;
    tot_n[i] = 0.;
  }

  for(std::vector<struct ICPV_data>::const_iterator it=ICPV_data.begin();
      it!=ICPV_data.end(); ++it){
    struct ICPV_data ev(*it);

    const double delt = ev.icpv_rec.vt_time - ev.icpv_asc.vt_time;

    for(unsigned j=0; j < nbin; ++j){
      ev.icpv_rec.vt_time = -dt_lim + (2*dt_lim / nbin)*(j+0.5);
      ev.icpv_asc.vt_time = 0;

      ev.icpv_flv = 1;
      ev.need_calc = 1;
      tot_p[j] += PDF_sin_cos(ev);

      ev.icpv_flv = -1;
      ev.need_calc = 1;
      tot_n[j] += PDF_sin_cos(ev);
    }
  }

  FILE *fp = fopen(file,"w");
  for(unsigned i=0; i<nbin; ++i){
    const double delt = -dt_lim + (2*dt_lim / nbin)*(i+0.5);
    fprintf(fp,"%+6.4le %+6.4le %+6.4le\n", delt, tot_p[i], tot_n[i]);
  }
  fclose(fp);
}

int main(int argc, char **argv)
{
  if( argc != 2){
    fprintf(stderr,"usage: test_tatami_jpsiks [datafile]\n");
    exit(3);
  }

  const char *datafile = argv[1];

  mninit();
  
  mnreset();

  read_param();
  
  read_data(datafile);

  if(ICPV_data.size() == 0){
    fprintf(stderr, "Error : No events in \"%s\" file!!\n",
            datafile);
    exit(12);
  }  

  mnfunc(func_belle);
  
  mnfit(MINUIT_MINIMIZE);
  mnfit(MINUIT_MINOS   );

  //load central value
  int ipar  =MINUIT_par.size();
  double errneg[ipar], errpos[ipar], errpar[ipar];
  for(unsigned int i=0; i<ipar; ++i){
    *MINUIT_par[i] = mnvalue(i, &errneg[i], &errpos[i], &errpar[i]);
  }
  plot_data("plot.d");
  plot_fit("plot.f");
  
  mnreset();
}
