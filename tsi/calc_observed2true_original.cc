/* calc_observed2true.cc */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "minuit.h"

#include "calc_coeffs.h"


static double Phi1, Phi3, rTSI, deltaTSI;
static double observedScp, observedAcp;

static double N_cpod_SLtag_qpos  = 1.0, N_cpod_SLtag_qneg  = 1.0;
static double N_cpod_HADtag_qpos = 1.0, N_cpod_HADtag_qneg = 1.0;
static double N_cpev_SLtag_qpos  = 1.0, N_cpev_SLtag_qneg  = 1.0;
static double N_cpev_HADtag_qpos = 1.0, N_cpev_HADtag_qneg = 1.0;


/*------------------------------------------------------------------------------*/
/* read init file                                                               */
/*------------------------------------------------------------------------------*/
static void
read_init(const char *file)
{
	FILE *fp;
	int nparam = 0;
	char buf[256];


	fp = fopen(file,"r");
	assert(fp);

	while( fgets(buf,sizeof(buf),fp) )
	{
		char key[64];
		double val;

		const int ntoken = sscanf(buf, "%s %le", key, &val);

		if( ntoken<2 ) continue;
		if( key[0]=='#' || key[0]=='!' ) continue;

		if     ( !strcmp(key,"phi1")               ){ Phi1               = val; nparam++;}
		else if( !strcmp(key,"phi3")               ){ Phi3               = val; nparam++;}
		else if( !strcmp(key,"rTSI")               ){ rTSI               = val; nparam++;}
		else if( !strcmp(key,"deltaTSI")           ){ deltaTSI           = val; nparam++;}
		else if( !strcmp(key,"observedAcp")        ){ observedAcp        = val; nparam++;}
		else if( !strcmp(key,"observedScp")        ){ observedScp        = val; nparam++;}
		else if( !strcmp(key,"N_cpod_SLtag_qpos")  ){ N_cpod_SLtag_qpos  = val; nparam++;}
		else if( !strcmp(key,"N_cpod_SLtag_qneg")  ){ N_cpod_SLtag_qneg  = val; nparam++;}
		else if( !strcmp(key,"N_cpod_HADtag_qpos") ){ N_cpod_HADtag_qpos = val; nparam++;}
		else if( !strcmp(key,"N_cpod_HADtag_qneg") ){ N_cpod_HADtag_qneg = val; nparam++;}
		else if( !strcmp(key,"N_cpev_SLtag_qpos")  ){ N_cpev_SLtag_qpos  = val; nparam++;}
		else if( !strcmp(key,"N_cpev_SLtag_qneg")  ){ N_cpev_SLtag_qneg  = val; nparam++;}
		else if( !strcmp(key,"N_cpev_HADtag_qpos") ){ N_cpev_HADtag_qpos = val; nparam++;}
		else if( !strcmp(key,"N_cpev_HADtag_qneg") ){ N_cpev_HADtag_qneg = val; nparam++;}
		else continue;
	}

	fclose(fp);

	assert(nparam==14);

	if( observedAcp*observedAcp + observedScp*observedScp >= 1 ){
		fprintf(stderr, "Scp^2 + Acp^2 >=1\n");
		exit(1);
	}
}


/*------------------------------------------------------------------------------*/
/* FCN                                                                          */
/*------------------------------------------------------------------------------*/
void
true2tsi(
	const double trueAcp, const double trueScp,
	const double phi1, const double phi3, const double rtsi, const double deltatsi,
	double *tsiAcp, double *tsiScp
){
	class PDFcoeffs coeffs_cpod_true(+trueAcp, +trueScp);
	class PDFcoeffs coeffs_cpod_with_tsi = coeffs_cpod_true;
	coeffs_cpod_with_tsi.apply_TSI(phi1, phi3, rtsi, deltatsi);

	class PDFcoeffs coeffs_cpev_true(+trueAcp, -trueScp);
	class PDFcoeffs coeffs_cpev_with_tsi = coeffs_cpev_true;
	coeffs_cpev_with_tsi.apply_TSI(phi1, phi3, rtsi, deltatsi);


	const double Acp_cpod_SLtag_qpos  = +coeffs_cpod_true.cos(+1);
	const double Acp_cpod_SLtag_qneg  = -coeffs_cpod_true.cos(-1);
	const double Acp_cpod_HADtag_qpos = +coeffs_cpod_with_tsi.cos(+1);
	const double Acp_cpod_HADtag_qneg = -coeffs_cpod_with_tsi.cos(-1);

	const double Acp_cpev_SLtag_qpos  = +coeffs_cpev_true.cos(+1);
	const double Acp_cpev_SLtag_qneg  = -coeffs_cpev_true.cos(-1);
	const double Acp_cpev_HADtag_qpos = +coeffs_cpev_with_tsi.cos(+1);
	const double Acp_cpev_HADtag_qneg = -coeffs_cpev_with_tsi.cos(-1);


	const double Scp_cpod_SLtag_qpos  = +coeffs_cpod_true.sin(+1);
	const double Scp_cpod_SLtag_qneg  = -coeffs_cpod_true.sin(-1);
	const double Scp_cpod_HADtag_qpos = +coeffs_cpod_with_tsi.sin(+1);
	const double Scp_cpod_HADtag_qneg = -coeffs_cpod_with_tsi.sin(-1);

	const double Scp_cpev_SLtag_qpos  = -coeffs_cpev_true.sin(+1);
	const double Scp_cpev_SLtag_qneg  = +coeffs_cpev_true.sin(-1);
	const double Scp_cpev_HADtag_qpos = -coeffs_cpev_with_tsi.sin(+1);
	const double Scp_cpev_HADtag_qneg = +coeffs_cpev_with_tsi.sin(-1);


	const double N_tot =
		+ N_cpod_SLtag_qpos + N_cpod_SLtag_qneg + N_cpod_HADtag_qpos + N_cpod_HADtag_qneg
		+ N_cpev_SLtag_qpos + N_cpev_SLtag_qneg + N_cpev_HADtag_qpos + N_cpev_HADtag_qneg;

	const double AveAcp = 1. / N_tot * (
		+ Acp_cpod_SLtag_qpos  * N_cpod_SLtag_qpos
		+ Acp_cpod_SLtag_qneg  * N_cpod_SLtag_qneg
		+ Acp_cpod_HADtag_qpos * N_cpod_HADtag_qpos
		+ Acp_cpod_HADtag_qneg * N_cpod_HADtag_qneg
		+ Acp_cpev_SLtag_qpos  * N_cpev_SLtag_qpos
		+ Acp_cpev_SLtag_qneg  * N_cpev_SLtag_qneg
		+ Acp_cpev_HADtag_qpos * N_cpev_HADtag_qpos
		+ Acp_cpev_HADtag_qneg * N_cpev_HADtag_qneg
	);

	const double AveScp = 1. / N_tot * (
		+ Scp_cpod_SLtag_qpos  * N_cpod_SLtag_qpos
		+ Scp_cpod_SLtag_qneg  * N_cpod_SLtag_qneg
		+ Scp_cpod_HADtag_qpos * N_cpod_HADtag_qpos
		+ Scp_cpod_HADtag_qneg * N_cpod_HADtag_qneg
		+ Scp_cpev_SLtag_qpos  * N_cpev_SLtag_qpos
		+ Scp_cpev_SLtag_qneg  * N_cpev_SLtag_qneg
		+ Scp_cpev_HADtag_qpos * N_cpev_HADtag_qpos
		+ Scp_cpev_HADtag_qneg * N_cpev_HADtag_qneg
	);

	*tsiAcp = AveAcp;
	*tsiScp = AveScp;
}


static double
func_belle(double *par) /* this is FCN */
{
	const double trueAcp = par[0];
	const double trueScp = par[1];
	double tsiAcp, tsiScp;

	true2tsi(
		trueAcp, trueScp,
		Phi1, Phi3, rTSI, deltaTSI,
		&tsiAcp, &tsiScp
	);

	return (observedAcp-tsiAcp)*(observedAcp-tsiAcp) + (observedScp-tsiScp)*(observedScp-tsiScp);
}


/*------------------------------------------------------------------------------*/
/* MINUIT I/F                                                                   */
/*------------------------------------------------------------------------------*/
static const double
do_mnfit(const char *initfile)
{
	double fcn;

	mninit();
	mnreset();

	read_init(initfile);

	mndump(MINUIT_WARN_OFF, MINUIT_DUMP_NONE, MINUIT_STEP_OFF);
	mnstrategy(1);
	mnfunc(func_belle);

	mnparam("trueAcp", observedAcp, 0, observedAcp-0.10, observedAcp+0.10, 0.001);
	mnparam("trueScp", observedScp, 0, observedScp-0.10, observedScp+0.10, 0.001);

	mnfit(MINUIT_MINIMIZE);
	// mnfit(MINUIT_MINIMIZE);
	// mnfit(MINUIT_HESSE   );
	// mnfit(MINUIT_MINOS   );

	// mnresult('e');

	/* fill result */
	const double trueAcp = mnvalue(0,NULL,NULL,NULL);
	const double trueScp = mnvalue(1,NULL,NULL,NULL);

	mnend();


	printf("observedAcp = %f, observedScp = %f\n", observedAcp, observedScp);
	printf("trueAcp     = %f, trueScp     = %f\n", trueAcp,     trueScp    );
}


int
main(int argc, char *argv[])
{
	if( argc!=2 ){
		fprintf(stderr, "usage: %s <init-file>\n", argv[0]);
		exit(1);
	}


	do_mnfit(argv[1]);


	return 0;
}

