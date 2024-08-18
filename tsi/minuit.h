/* minuit.h */

#ifndef MINUIT_H
#define MINUIT_H

#ifdef __cplusplus
#define LINKAGE "C"
#else
#define LINKAGE
#endif

typedef struct mnparam {
	/* inputs */
	char   mp_name[48];
	int    mp_fix;
	int    mp_seed;
	double mp_rnd;
	double mp_init;
	double mp_min, mp_max;
	double mp_acc;
	/* outputs */
	double mp_value;
	double mp_errneg, mp_errpos, mp_errpar;
} mnparam_t;

typedef double (*mnfcn_t)(double*);


extern LINKAGE void mninit(void);
extern LINKAGE void mnreset(void);
extern LINKAGE void mnend(void);
extern LINKAGE void mnparam(const char *name, const double init, const int fix, const double min, const double max, const double acc);
extern LINKAGE void mnblind(const char *name, const int blind);
extern LINKAGE double mnvalue(int num, double *errneg, double *errpos, double *errpar);
extern LINKAGE void mnpset(mnparam_t *mnp);
extern LINKAGE void mnpget(const int num, mnparam_t *mnp);
extern LINKAGE int mnpnum(void);
extern LINKAGE void mnfunc(mnfcn_t fcn);
extern LINKAGE void mnfit(const int which, ...);
extern LINKAGE int mnstat(void);
extern LINKAGE void mnemat(double *emat);
extern LINKAGE int mncont(const int p1, const int p2, const int np, double *x, double *y);
extern LINKAGE void mndump(const int warn, const int dump, const int step);
extern LINKAGE double mnminimum(void);
extern LINKAGE void mnresult(const char format);
extern LINKAGE void mnstrategy(const int opt);

#define MINUIT_FIXED  (1)
#define MINUIT_FLOAT  (0)

#define MINUIT_MINIMIZE	(1)
#define MINUIT_HESSE	(2)
#define MINUIT_MINOS	(3)
#define MINUIT_MINOS_WITH_NERVOUS_LIST	(4)

#define MINUIT_WARN_OFF		( 0)
#define MINUIT_WARN_ON		( 1)
#define MINUIT_DUMP_NONE	(-1)
#define MINUIT_DUMP_NORMAL	( 0)
#define MINUIT_DUMP_VERBOSE	( 1)
#define MINUIT_STEP_OFF		( 0)
#define MINUIT_STEP_ON		( 1)

#endif /* MINUIT_H */

/*
	$Id: minuit.h,v 1.4 2003/08/05 02:32:33 higuchit Exp $
	$Log: minuit.h,v $
	Revision 1.4  2003/08/05 02:32:33  higuchit
	2003 summer version.
	
	Revision 1.3  2001/09/21 08:54:15  higuchi
	typedef bug.
	
	Revision 1.2  2001/09/20 06:49:29  higuchi
	protoized
	
	Revision 1.1.1.1  2001/07/14 05:05:35  higuchi
	initial import
	
*/

