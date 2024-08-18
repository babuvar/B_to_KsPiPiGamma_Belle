/* minuit.c */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "./minuit.h"

extern void mninit_(), mnexcm_(), mncomd_(), mnparm_(), mnpout_(), mnerrs_(), mnemat_(), mncont_(), mnstat_();


typedef struct mnopt_t {
	double mo_warn;
	double mo_dump;
	int    mo_step;
	double mo_strategy;
} mnopt_t;


static mnopt_t MNopt = { 0.0, -1.0, 1, 0.0, };
static mnparam_t MNpar[101];
static mnfcn_t MNfcn;
static int MNtouch = 1;


static void mnfcn(int *npar, double *grad, double *fval, double *xval, int *iflag, void *futil)
{
	int i, pnum=mnpnum();
	double xval_copy[101];

	if( MNopt.mo_step ){
		for( i=0; i<pnum; i++ )
			printf("%-16s: %+le\n", MNpar[i].mp_name, xval[i]);
	}

	for( i=0; i<pnum; i++ )
		xval_copy[i] = !MNpar[i].mp_seed ? xval[i] : xval[i]-MNpar[i].mp_rnd;

	*fval = MNfcn ? (*MNfcn)(xval_copy) : 0.0;

	if( MNopt.mo_step ){
		printf("FCN             : %+le\n\n", *fval);
	}
}


int mnpfind(const char *name)
{
	int i, pnum=mnpnum();

	for( i=0; i<pnum; i++ )
		if( !strcmp(MNpar[i].mp_name, name) ) return i;

	return -1;
}


static mnparam_t *mnfind(const char *name)
{
	int i, pnum=mnpnum();

	for( i=0; i<pnum; i++ )
		if( !strcmp(MNpar[i].mp_name, name) ) return &MNpar[i];

	return NULL;
}


void mninit(void)
{
	int ird = 5, iwr = 6, isav = 7;
	double opt;
	int one = 1, ret;
	static int first = 1;

	MNtouch = 1;

	mnend();

	if( first ) mninit_(&ird, &iwr, &isav);
	first = 0;

	opt = -1.0;
	mnexcm_(mnfcn, "SET PRI", &opt, &one, &ret, NULL, 7);
}


void mnreset(void)
{
	int i, one = 1, ret;
	double opt;

	opt = 0.0;
	mnexcm_(mnfcn, "CLEAR", &opt, &one, &ret, NULL, 5);

	MNtouch = 1;
	MNfcn = NULL;
	for(i=0; i<100; i++) MNpar[i].mp_name[0] = '\0';
}


void mnparam(const char *name, const double init, const int fix, const double min, const double max, const double acc)
{
	mnparam_t mnp;

	strncpy( mnp.mp_name, name, sizeof(mnp.mp_name) );
	mnp.mp_fix    = fix;
	mnp.mp_seed   = 0;
	mnp.mp_rnd    = 0;
	mnp.mp_init   = init;
	mnp.mp_min    = min;
	mnp.mp_max    = max;
	mnp.mp_acc    = acc;

	mnp.mp_value  = init;
	mnp.mp_errneg = 0.0;
	mnp.mp_errpos = 0.0;
	mnp.mp_errpar = 0.0;

	mnpset(&mnp);
}


void mnblind(const char *name, const int seed)
{
	mnparam_t *mnp = mnfind(name);

	if( !mnp ) return;

	mnp->mp_seed = seed;
	if( seed ){
		double rnd = 0;
		const double min=mnp->mp_min, max=mnp->mp_max, del=max-min;

		if( del!=0 ){
			srand48((long)seed);
			rnd = 8.*del*drand48() + (min-4.*del);

			mnp->mp_rnd   = rnd;
			mnp->mp_init += rnd;
			mnp->mp_min  += rnd;
			mnp->mp_max  += rnd;
		} else {
			srand48((long)seed);
			rnd = 100.*drand48();

			mnp->mp_rnd   = rnd;
			mnp->mp_init += rnd;
		}
	}

	MNtouch = 1;
}


void mnpset(mnparam_t *mnp)
{
	int pnum = mnpnum();
	if( pnum>=100 ) return;

	memcpy(&MNpar[pnum],mnp,sizeof(mnparam_t));
	MNtouch = 1;
}


int mnpnum(void)
{
	int i = 0;

	while( MNpar[i].mp_name[0] ) i++;

	return i;
}


void mnfunc(mnfcn_t fcn)
{
	MNtouch = 1;

	MNfcn = fcn;
}


void mnfit(const int which, ...)
{
	int i, pnum=mnpnum(), num, ret;
	double opt[2];
	int zero = 0, one = 1, two = 2;

	if( MNtouch ){
		opt[0] = MNopt.mo_dump;
		mnexcm_(mnfcn, "SET PRI", opt, &one, &ret, NULL, 7);

		opt[0] = 0.0;
		if( MNopt.mo_warn ){
			mnexcm_(mnfcn, "SET WAR", opt, &one, &ret, NULL, 7);
		} else {
			mnexcm_(mnfcn, "SET NOW", opt, &one, &ret, NULL, 7);
		}

		opt[0] = MNopt.mo_strategy;
		mnexcm_(mnfcn, "SET STR", opt, &one, &ret, NULL, 7);

		for( i=0, num=1; i<pnum; i++, num++ ){
			int erflg;
			mnparm_(&num, MNpar[i].mp_name, &MNpar[i].mp_init, &MNpar[i].mp_acc,
				&MNpar[i].mp_min, &MNpar[i].mp_max, &erflg, strlen(MNpar[i].mp_name));
			if( MNpar[i].mp_fix ){
				opt[0] = (double)num;
				mnexcm_(mnfcn, "FIX", opt, &one, &ret, NULL, 3);
			}
		}
	}

	switch( which ){
		case MINUIT_MINIMIZE:
			mnexcm_(mnfcn, "MINI", opt, &zero, &ret, NULL, 4);
			break;

		case MINUIT_HESSE:
			mnexcm_(mnfcn, "HESS", opt, &zero, &ret, NULL, 4);
			break;

		case MINUIT_MINOS:
			mnexcm_(mnfcn, "MINO", opt, &zero, &ret, NULL, 4);
			break;

		case MINUIT_MINOS_WITH_NERVOUS_LIST:
		{
			va_list ap;
			int *nervous_list;
			int i, nopt, nnervous;
			double opt[1024];
			int pnum;

			va_start(ap, which);
			nervous_list = va_arg(ap, int*);
			nnervous     = va_arg(ap, int);
			va_end(ap);

			pnum = mnpnum();

			opt[0] = 10000;
			for( i=0; i<nnervous; i++ ){
				int idx;
				
				idx = nervous_list[i];
				if(idx<0 || idx>=pnum){
					fprintf(stderr, "mnfit(MINUIT_MINOS_WITH_NERVOUS_LIST): bad parno=%d\n", idx);
					exit(1);
				}
				idx++ /* ++ for C->FORTRAN conversion */;
				opt[i+1] = idx /* +1 for 10000 [maxcalls] */;
			}
			nopt = 1 + nnervous;

			mnexcm_(mnfcn, "MINO", opt, &nopt, &ret, NULL, 4);
			break;
		}

		default:
			1;
	}

	for( i=0, num=1; i<pnum; i++, num++ ){
		char name[16];
		double val, min, max, globcc;
		int varbl;

		mnpout_(&num,
			name, &MNpar[i].mp_value, &MNpar[i].mp_errpar, &min, &max, &varbl, 16);
		mnerrs_(&num,
			&MNpar[i].mp_errpos, &MNpar[i].mp_errneg, &MNpar[i].mp_errpar, &globcc);
	}

	MNtouch = 0;
}


int mnstat(void)
{
	double fmin, fedm, errdef;
	int npari, nparx, istat=0;

	mnstat_(&fmin, &fedm, &errdef, &npari, &nparx, &istat);

	return istat;
}


void mnpget(const int num, mnparam_t *mnp)
{
	memcpy(mnp, &MNpar[num], sizeof(mnparam_t));
}


double mnvalue(int num, double *errneg, double *errpos, double *errpar)
{
	if( errneg ) *errneg = MNpar[num].mp_errneg;
	if( errpos ) *errpos = MNpar[num].mp_errpos;
	if( errpar ) *errpar = MNpar[num].mp_errpar;

	return MNpar[num].mp_value;
}


void mnemat(double *emat)
{
	int pnum = mnpnum();

	mnemat_(emat, &pnum);
}


int mncont(const int p1, const int p2, const int np, double *x, double *y)
{
	int nfound;
	int zero=0;

	double *_x = x, *_y = y;

	if( !_x ) _x = (double*)malloc(sizeof(double)*np);
	if( !_y ) _y = (double*)malloc(sizeof(double)*np);

	mncont_(mnfcn, &p1, &p2, &np, _x, _y, &nfound, zero);

	if( !x ) free(_x);
	if( !y ) free(_y);

	return nfound;
}


double mnminimum(void)
{
	int i;
	int pnum = mnpnum();
	double *xval;
	double fcn;
	
	xval = malloc(sizeof(double)*pnum);

	for( i=0; i<pnum; i++ ) xval[i] = !MNpar[i].mp_seed ? MNpar[i].mp_value : MNpar[i].mp_value-MNpar[i].mp_rnd;

	fcn = MNfcn ? (*MNfcn)(xval) : 0.0;

	free(xval);

	return fcn;
}


void mnresult(const char format)
{
	int i;
	int pnum = mnpnum();

	printf("\n");
	printf(" FCN= %f\n\n", mnminimum());
	printf("  EXT PARAMETER                            PARABOLIC           MINOS ERRORS\n");
  printf("  NO.   NAME                  VALUE          ERROR        NEGATIVE       POSITIVE\n");

	for( i=0; i<pnum; i++ )
	{
		printf(">>  %2d  ", i+1);
		printf("%-18s  ", MNpar[i].mp_name);
		if( format == 'f' || format == 'F' ){
			printf("%+f  ",   MNpar[i].mp_value);
		} else {
			printf("%+e  ",   MNpar[i].mp_value);
		}
		if( MNpar[i].mp_fix ){
			printf("    fixed    ");
		} else if( MNpar[i].mp_acc==0 ){
			printf("   constant  ");
		} else {
			if( format == 'f' || format == 'F' ){
				MNpar[i].mp_errpar!=0 ? printf("%+f  ",   MNpar[i].mp_errpar) : printf("             ");
				MNpar[i].mp_errneg!=0 ? printf("%+f  ",   MNpar[i].mp_errneg) : printf("             ");
				MNpar[i].mp_errpos!=0 ? printf("%+f  ",   MNpar[i].mp_errpos) : printf("             ");
			} else {
				MNpar[i].mp_errpar!=0 ? printf("%+e  ",   MNpar[i].mp_errpar) : printf("             ");
				MNpar[i].mp_errneg!=0 ? printf("%+e  ",   MNpar[i].mp_errneg) : printf("             ");
				MNpar[i].mp_errpos!=0 ? printf("%+e  ",   MNpar[i].mp_errpos) : printf("             ");
			}
		}
		printf("\n");
	}

	printf("\n");
}


void mndump(const int warn, const int dump, const int step)
{
	MNopt.mo_warn = (double)warn;
	MNopt.mo_dump = (double)dump;
	MNopt.mo_step = step;
	MNtouch = 1;
}


void mnstrategy(const int opt)
{
	MNopt.mo_strategy = (double)opt;
	MNtouch = 1;
}


void mnend(void)
{
	mnreset();
}

