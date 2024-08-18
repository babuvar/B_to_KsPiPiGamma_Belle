//
// $Id: test_tatami_toymc.cc 9948 2006-12-05 06:11:00Z katayama $
//
// test program using the package "tatami".
// example of toyMC generator for J/psiKs.
//
// $Log$
// Revision 1.3  2005/09/22 05:16:46  kohji
// 1. Fixed the out-of-range bug in test_tatami_toymc.cc
//        and norm_Mf() bug in libcnv.cc, pointed out by Sasa.
// 2. Added resolution parameters used for LP05,
//  dtres_param_svd2_lp05, dtres_param_posi_error_svd2_lp05,
//  	 dtres_param_nega_error_svd2_lp05.
//
// Revision 1.2  2003/08/06 06:43:31  katayama
// New version from Nakadaira san
//
// Revision 1.1  2002/10/24 12:35:58  katayama
// new from Sumisawa san
//
//

#include "belle.h"
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <cmath>
#include <sys/types.h>
#include <regex.h>

#include "tatami/tatami.h"

#ifdef __sun
#include <ieeefp.h>
#endif

#include "belleutil/debugout.h"

using namespace Belle;


extern "C" void ranecu_(float[], const int*, const int*);
extern "C" void ranecq_(int*, int*, int*, const char*, int);

/*** random number generator ***/

static int ranecu_seq(1);
static int counter(0);

static void rand_seed()
{
  static int first=1;

  if( !first ) return;

  first = 0;

  //srand48(getpid());
  int iseq(ranecu_seq);
  int seed1(getpid());
  int seed2((int)time(NULL));

  // dout(Debugout::DDEBUG,"test_tatami_toymc")<<"Seed1 = "<<seed1<<" : Seed2 = "<<seed2<<std::endl;

  ranecq_(&seed1, &seed2, &iseq, "S", 1);
  //srand48();
}

static double make_rand(void)
{
  float ran;
  int size(1);
  int iseq(ranecu_seq);

  ranecu_(&ran, &size, &iseq);

  if(++counter%100000 == 0){
    int s1, s2;
    ranecq_(&s1, &s2, &iseq, "R", 1);
    // dout(Debugout::DDEBUG,"test_tatami_toymc")<<"Counter = "<<counter
    //<<" Seed1 = "<<s1<<" : Seed2 = "<<s2<<std::endl;
  }
  
  return ran;
}

double rand_exp(double tau)
{
  rand_seed();

  double x = make_rand();
  double y = -log(x)*tau;
  
  return y;
}

double rand_gauss(const double mean, const double sigma)
{
  double x1, x2, y1, y2;

  rand_seed();

  x1 = make_rand();
  x2 = make_rand()*2.*M_PI;
  y1 = sqrt(-2.*log(x1))*cos(x2);
  y2 = sqrt(-2.*log(x1))*sin(x2);

  return y1*sigma+mean;
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


static const double sig_func(const double mb, const double de,
			     const double ebeam,
			     double* const psig=0x0, double* const pbg=0x0)
{
  static const double de_mean(0.54653E-03);
  static const double de_sigma(0.62154E-02);
  static const double de_fmain(0.63769);
  static const double de_sigma2(0.15816E-01);
  static const double mb_mean(5.2796);
  static const double mb_sigma(0.24917E-02);
  static const double sig_area(0.976245);
  static const double argus_a(-43.753);
  static const double bg_de1(-1.8516);

  static const double de_min(-0.04);
  static const double de_max(0.04);
  static const double mb_min(5.2694);
  static const double mb_max(5.2894);

  //signal
  const double int_sig =
    norm_gaussian(mb_min, mb_max, mb_mean, mb_sigma)
    * (de_fmain * norm_gaussian(de_min, de_max, de_mean, de_sigma)
       + (1-de_fmain) * norm_gaussian(de_min, de_max, de_mean, de_sigma2));
  
  const double sig = sig_area / int_sig
    * gaussian(mb, mb_mean, mb_sigma)
    * (de_fmain * gaussian(de, de_mean, de_sigma)
       + (1-de_fmain) * gaussian(de, de_mean, de_sigma2));
  
  //background
  const double int_bg =
    int_argus(mb_min, mb_max, 0., ebeam, argus_a)
    * (1.0 + bg_de1*0.5*(de_min + de_max))*(de_max - de_min);

  const double bg = (1 - sig_area) / int_bg
    * argus(mb, 0., ebeam, argus_a)
    * (1.0 + bg_de1*de);

  if(psig != 0x0) *psig = sig;
  if(pbg != 0x0) *pbg = bg;
  
  return sig/(sig + bg);
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
  int icpv_fbtg_mode;
  double icpv_wtag;
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
};

//--------------------------------
static double calc_pbstar(const double mbc, const double eb)
{
  return sqrt(eb*eb - mbc*mbc);
}

static void generate_kine(struct ICPV_data &ev, int &is_sig)
{
  rand_seed();

  const double ebeam(ev.icpv_ebeam);

  static int kine_flag(1);
  static double cmax(0);

  if(kine_flag == 1){
    kine_flag = 0;

    double s, b;
    sig_func(5.2794, 0., ebeam, &s, &b);
    cmax = 1.1 * (s + (b>0 ? b : 0.));
  }

  // signal box
  static const double sigbox_mb_min(5.2694);
  static const double sigbox_mb_max(5.2894);
  static const double sigbox_de_min(-0.04);
  static const double sigbox_de_max(0.04);

  double mbc, deltae;
  double c, x, sig, bg;
  do{
    mbc = sigbox_mb_min + (sigbox_mb_max-sigbox_mb_min)*make_rand();
    deltae = sigbox_de_min + (sigbox_de_max-sigbox_de_min)*make_rand();

    sig_func(mbc, deltae, ebeam, &sig, &bg);

    c = sig + (bg>0 ? bg : 0);
    x = cmax * make_rand();
  }while(x > c);

  if(x <= sig)
    is_sig = 1;
  else
    is_sig = 0;
  
  ev.icpv_mbc = mbc;
  ev.icpv_deltae = deltae;
  ev.icpv_pbs = calc_pbstar(mbc, ebeam);
}

static double generate_ebeam(void)
{
  return 5.290024912;
}

// # of tracks
static const int generate_ntrack_cp(void)
{
  rand_seed();

  return make_rand() < 143./1550. ? 1 : 2;
}


static const int generate_ntrack_tag(void)
{
  rand_seed();

  return make_rand() < 427./1550. ? 1 : 2;
}

//------------------------------------------
// cos(thetaB*)
static const double
generate_cos_thetaB_sig(void)
{
  rand_seed();

  double x = make_rand();
  double y = 2.0*cos((asin(sqrt(x))+2.0*M_PI)*2.0/3.0);
  
  return y;
}


static const double
generate_cos_thetaB_bg(void)
{
  rand_seed();

  return 1 - 2*make_rand();
}

//-------------------------------------------------------------------------
// vertex errors
static const double rand_pseudo_threshold(const double c1, const double c2,
                                          const double tau)
{
  rand_seed();

  const double max = (c2-c1)*.5+tau;
  const double x = make_rand()*max;

  double Fx;

  if( x<(c2-c1)*.5 ){
    Fx = sqrt(2*(c2-c1)*x)+c1;
  } else {
    Fx = c2-tau*log(1-(x+(c1-c2)*.5)/tau);
  }

  return Fx;
}

static const double generate_sigma_ntrk_cp(void)
{
  const double sigma_min  = 0.25821E-02;
  const double sigma_peak = 0.36795E-02;
  const double sigma_tau  = 0.19887E-02;

  return rand_pseudo_threshold(sigma_min, sigma_peak, sigma_tau);
}


static const double generate_sigma_ntrk_tag(void)
{
  const double sigma_min  = 0.27346E-02;
  const double sigma_peak = 0.51262E-02;
  const double sigma_tau  = 0.35819E-02;

  return rand_pseudo_threshold(sigma_min, sigma_peak, sigma_tau);
}

static const double generate_sigma_1trk(void)
{
  const double sigma_min  = 0.27974E-02;
  const double sigma_peak = 0.50779E-02;
  const double sigma_tau  = 0.13028E-01;

  return rand_pseudo_threshold(sigma_min, sigma_peak, sigma_tau);
}

//--------------------------------------------------------------------------
// xi

static const double rand_three_exp(const double f1, const double f2,
                                   const double t1, const double t2,
                                   const double t3)
{
  rand_seed();

  double tau;

  if( make_rand()<f1 )
    tau = t1;
  else
    if( make_rand()<f2 )
      tau = t2;
    else
      tau = t3;

  return rand_exp(tau);
}

static const double generate_xi_cp(void)
{
  const double f1 = 0.91058;
  const double f2 = 0.98164;
  const double t1 = 0.46210;
  const double t2 = 2.8456;
  const double t3 = 22.291;

  for(;;){
    const double xi = rand_three_exp(f1, f2, t1, t2, t3);
    if( xi<=100 ) return xi;
  }
}

static const double generate_xi_tag(void)
{
  const double f1 = 0.68444;
  const double f2 = 0.94250;
  const double t1 = 0.99754;
  const double t2 = 4.8307;
  const double t3 = 26.542;

  for(;;){
    const double xi = rand_three_exp(f1, f2, t1, t2, t3);
    if( xi<=100 ) return xi;
  }
}
//---------------------------------
// generate wrong tag fraction for flavor tagging.
//
inline double sum_in_range(const double* const hist,
                           const int i, const int j){
  const int ll = (i>j) ? j : i;
  const int ul = (i>j) ? i : j;
  double sum = 0.0;
  for(int n=ll;n<=ul;++n){
    sum += hist[n];
  }
  return sum;
}

static void generate_hamlet(struct ICPV_data &ev)
{
  rand_seed();

  ev.icpv_flv = ev.icpv_mc_flv;

  ev.icpv_fbtg_mode = 0;

  const double fhist[] = {
    0.399, 0.146, 0.104, 0.122, 0.094, 0.137
  };
  const int nhist = sizeof(fhist)/sizeof(double);

  const double qhist[] = {
    (0.000+0.250)*0.5,
    (0.250+0.500)*0.5,
    (0.500+0.625)*0.5,
    (0.625+0.750)*0.5,
    (0.750+0.875)*0.5,
    (0.875+1.000)*0.5
  };
  const double whist[] = {
    0.458, 0.336, 0.228, 0.160, 0.112, 0.020
  };

  double tagq = make_rand()*(sum_in_range(fhist, 0, nhist-1));

  ev.icpv_wtag = -1;
  for(int i=0;i<nhist;++i){
    if(tagq>sum_in_range(fhist, 0, i)){
      continue;
    }
    ev.icpv_wtag = (1.0-qhist[i])*0.5;
    // assumed wtag fraction in CP fitter;
    if( make_rand()<whist[i]) ev.icpv_flv *= -1;
    break;
  }
}

//---------------------------------
// generate Dt distribution (w/o convlutions of Rdet and Rnp)
//
static const double generate_dt_Phys_Rk(const double sin2phi1,
					const double tau, const double dm,
					struct ICPV_data &ev)
{
  rand_seed();

  double a_k, c_k;
  CalcAkCk(ev.icpv_costh, ev.icpv_ecms, &a_k, &c_k, 0);

  double dt, x, c, akck;
  dt = rand_exp(tau);

  x = 2*make_rand()-1;
  if(x > c_k/a_k){
    akck = a_k + c_k;
    dt *= akck;
  }else{
    akck = a_k - c_k;
    dt *= -akck;
  }

  const double a(c_k * tau * dm / a_k);
  c = 1/(1 + a*a)
    * sin2phi1 * (sin(dm*dt/akck) - a * cos(dm*dt/akck));
  x = 2*make_rand()-1;

  if(x<c){
    ev.icpv_mc_flv = 1;
  }else{
    ev.icpv_mc_flv = -1;
  }  

  return dt;
}

//
// RdetRnp component
//
static const double generate_sig_res(const struct vertex &rec,
				     const struct vertex &asc, 
				     const dtres_param_t* dtres_param)
{
  rand_seed();
  const int keeptagl(0);

  const double dt_llmt(dt_resol_global::dt_llmt);
  const double dt_ulmt(dt_resol_global::dt_ulmt);

  const double c_max(1.2 * 
		     RdetRnp(0., 0,
			     rec.vt_ntrk, rec.vt_err, rec.vt_chi2, rec.vt_ndf,
			     asc.vt_ntrk, asc.vt_err, asc.vt_chi2, asc.vt_ndf, keeptagl,
			     dtres_param));
  double res, c, x;
  do{
    res = dt_llmt + make_rand() * (dt_ulmt - dt_llmt);
    c = RdetRnp(res, 0,
		rec.vt_ntrk, rec.vt_err, rec.vt_chi2, rec.vt_ndf,
		asc.vt_ntrk, asc.vt_err, asc.vt_chi2, asc.vt_ndf, keeptagl,
		dtres_param);
    x = make_rand() * c_max;
  }while(x > c);

  return res;
}

//
// generate Dt distribution for signal (J/psi Ks)
//
static const double generate_dz_sig(const double sin2phi1,
				    struct ICPV_data &ev, 
				    const dtres_param_t* dtres_param)
{
  static const double tau_b(1.542);
  static const double m_d(0.489);

  struct vertex &rec = ev.icpv_rec;
  struct vertex &asc = ev.icpv_asc;

  const double dt_llmt(dt_resol_global::dt_llmt);
  const double dt_ulmt(dt_resol_global::dt_ulmt);

  double gen, res, dt;

  do{
    gen = generate_dt_Phys_Rk(sin2phi1, tau_b, m_d, ev);
    res = generate_sig_res(rec, asc, dtres_param);

    dt = gen + res;
  }while(dt <= dt_llmt || dt >= dt_ulmt);

  rec.vt_mc = gen / Coeff_cm2ps;
  
  return (dt / Coeff_cm2ps);
}

static const double PDF_bg(const double delt,
			   const struct vertex &rec, const struct vertex &asc)
{
  double pdf(0);

  double s_main, s_tail, f_tail, tau_bg, f_delta, mu_l, mu_d;
  
  mu_l = -0.12035;
  mu_d = mu_l;
  tau_bg = 1.7490;
  
  if(rec.vt_ntrk > 1 && asc.vt_ntrk > 1){
    s_main = 0.58894;
    s_tail = 3.2665;
    f_tail = 0.87145;
    f_delta = 0.68071;
  }else{
    s_main = 1.0683;
    s_tail = 3.0000;
    f_tail = 0.;
    f_delta = 0.55672;
  }

  const double sigma_main(s_main
			  * sqrt(rec.vt_err*rec.vt_err + asc.vt_err*asc.vt_err)
			  * Coeff_cm2ps);
  const double sigma_tail(s_tail * sigma_main);

  // lifetime compronent
  if(f_delta < 1){
    double pdf_l(0);
    if(f_tail < 1){
      pdf_l += (1 - f_tail) *
	Enp_conv_gauss(delt, tau_bg, tau_bg, mu_l, sigma_main);
    }
    if(f_tail > 0){
      pdf_l += f_tail *
	Enp_conv_gauss(delt, tau_bg, tau_bg, mu_l, sigma_tail);
    }
    
    pdf += pdf_l * (1 - f_delta);
  }

  // delta function
  if(f_delta > 0){
    double pdf_d(0);
    if(f_tail < 1){
      pdf_d += (1 - f_tail) * gaussian(delt, mu_d, sigma_main);
    }
    if(f_tail > 0){
      pdf_d += f_tail * gaussian(delt, mu_d, sigma_tail);
    }
    
    pdf += pdf_d * f_delta;
  }

  return pdf;
}

//
// generate Dt distribution for BG
//
static const double generate_dz_bg(const struct vertex &rec,
				   const struct vertex &asc)
{
  rand_seed();

  const double dt_llmt(dt_resol_global::dt_llmt);
  const double dt_ulmt(dt_resol_global::dt_ulmt);

  const double c_max(1.1 * PDF_bg(-0.12035, rec, asc));

  double dt, c, x;
  do{
    dt = dt_llmt + make_rand() * (dt_ulmt - dt_llmt);

    c = PDF_bg(dt, rec, asc);
    x = make_rand() * c_max;
  }while(x > c);

  return (dt / Coeff_cm2ps);
}

//
// Out lier component
//
static const double generate_dz_ol(const dtres_param_t* dtres_param)
{
  rand_seed();

  const double dt_llmt(dt_resol_global::dt_llmt);
  const double dt_ulmt(dt_resol_global::dt_ulmt);
  
  const double sigma(dtres_param->sig_ol);

  double dt;
  do{
    dt = rand_gauss(0., sigma);
  }while(dt <= dt_llmt || dt>= dt_ulmt);

  return (dt / Coeff_cm2ps);  
}

static void generate_event(const double sin2phi1, struct ICPV_data &ev)
{
  rand_seed();
  
  struct vertex &rec = ev.icpv_rec;
  struct vertex &asc = ev.icpv_asc;

  const int mc = 0;

  ev.icpv_exp = 13;
  ev.icpv_run = 0;
  ev.icpv_evt = 0;
  ev.icpv_frm = 0;
  ev.icpv_hepevt = 511;
  rec.vt_mc = 0;
  asc.vt_mc = 0;
  ev.icpv_mc_delta_t = 0;
  ev.icpv_mc_pbs = 0;
  ev.icpv_mc_costh = 0;
  ev.icpv_decay = DECAY_PSI_KS;
  asc.vt_pos = 0.;
  rec.vt_cl = 0.;
  asc.vt_cl = 0.;
  ev.icpv_costhtr = 0.;

  rec.vt_ntrk = generate_ntrack_cp();
  rec.vt_ndf = 2 * rec.vt_ntrk;
  asc.vt_ntrk = generate_ntrack_tag();
  asc.vt_ndf = 2 * asc.vt_ntrk;

  const dtres_param_t* dtres_param = get_dtres_param(ev.icpv_exp, mc, 
						     dt_resol_global::ichep2002_80);
  
  double f_ol(dtres_param->fol_mul);
  if(rec.vt_ntrk == 1){
    rec.vt_err = generate_sigma_1trk();
    rec.vt_chi2 = 1.;
    f_ol = dtres_param->fol_sgl;
  }else{
    rec.vt_err = generate_sigma_ntrk_cp();
    rec.vt_chi2 = generate_xi_cp() * rec.vt_ndf;
  }
  if(asc.vt_ntrk == 1){
    asc.vt_err = generate_sigma_1trk();
    asc.vt_chi2 = 1.;
    f_ol = dtres_param->fol_sgl;
  }else{
    asc.vt_err = generate_sigma_ntrk_tag();
    asc.vt_chi2 = generate_xi_tag() * asc.vt_ndf;
  }

  ev.icpv_ebeam = generate_ebeam();

  int is_sig;
  generate_kine(ev, is_sig);
  
  ev.icpv_ecms = sqrt(ev.icpv_pbs*ev.icpv_pbs +
		      dt_resol_global::mbzero*dt_resol_global::mbzero);

  if(is_sig){
    ev.icpv_costh = generate_cos_thetaB_sig();
  }else{
    ev.icpv_costh = generate_cos_thetaB_bg();
  }

  if(make_rand() < f_ol){
    ev.icpv_mc_flv = make_rand() < 0.5 ? +1 : -1;
    rec.vt_pos = generate_dz_ol(dtres_param);
  }else{
    if(is_sig){
      rec.vt_pos = generate_dz_sig(sin2phi1, ev, dtres_param);
    }else{
      ev.icpv_mc_flv = make_rand() < 0.5 ? +1 : -1;
      rec.vt_pos = generate_dz_bg(rec, asc);
    }
  }

  generate_hamlet(ev);
  ev.icpv_hepevt *= ev.icpv_mc_flv;
}

static void dump_event(const struct ICPV_data &ev)
{
  const struct vertex &rec = ev.icpv_rec;
  const struct vertex &asc = ev.icpv_asc;
  
  if(debugout("DUMP")) std::printf("%d %d %d %d "
				   "%d "
				   "%d "
				   "%le %le %le "
				   "%le %le "
				   "%d "
				   "%le %le %d %le %d %le "
				   "%le %le %d %le %d %le "
				   "0 "
				   "%d %d %le "
				   "%le %le %le %le %le %le\n",
				   ev.icpv_exp, ev.icpv_run, ev.icpv_evt, ev.icpv_frm,
				   ev.icpv_hepevt,
				   ev.icpv_mc_flv,
				   rec.vt_mc, asc.vt_mc, ev.icpv_mc_delta_t,
				   ev.icpv_mc_pbs, ev.icpv_mc_costh,
				   ev.icpv_decay,
				   rec.vt_pos, rec.vt_err, rec.vt_ntrk, rec.vt_chi2, rec.vt_ndf, rec.vt_cl,
				   asc.vt_pos, asc.vt_err, asc.vt_ntrk, asc.vt_chi2, asc.vt_ndf, asc.vt_cl,
				   ev.icpv_flv, ev.icpv_fbtg_mode, ev.icpv_wtag,
				   ev.icpv_deltae, ev.icpv_mbc, ev.icpv_pbs, ev.icpv_costh,
				   ev.icpv_ebeam, ev.icpv_costhtr);
}

int main(int argc, char **argv)
{
  if(argc != 3){
    dout(Debugout::ERR,"test_tatami_toymc") <<"usage: test_tatami_toymc [sin2phi1] [nevent]"<< std::endl;
    exit(3);
  }

  const double sin2phi1 = atof(argv[1]);
  const int nevent = atoi(argv[2]);

  for(int i=0; i<nevent; ++i){
    struct ICPV_data toy;
    generate_event(sin2phi1, toy);
    dump_event(toy);
  }
}
