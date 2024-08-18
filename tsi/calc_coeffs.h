/* calc_coeffs.h */


#pragma once


#include <math.h>
#include <assert.h>
#include <complex>



using namespace std;
typedef complex<double> c_double;


inline c_double
calc_lambda(const double Acp, const double Scp)
// lambda = etaCP x exp(-2*i*phi1)
// Scp(q=+1) ~ +Im(lambda) = + ( -etaCP x sin(2phi1) )
// Scp(q=-1) ~ -Im(lambda) = - ( -etaCP x sin(2phi1) )
{
	const double mag_lambda_2 = (1+Acp)/(1-Acp);
	const double mag_lambda   = sqrt(mag_lambda_2);
	const double arg_lambda = asin(Scp*(1+mag_lambda_2)/(2.*mag_lambda));
	return c_double(polar(mag_lambda,arg_lambda));
}



/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------*/

class PDFcoeffs
{
public:
	PDFcoeffs(
		const double Acp, const double Scp
	){
		this->m_Acp = Acp;
		this->m_Scp = Scp;
		this->m_lambda = calc_lambda(Acp, Scp);
		this->m_rTSI = this->m_deltaTSI = 0;
		this->update_TSI();
	}

	PDFcoeffs(
		const c_double lambda
	){
		this->m_lambda = lambda;
		const double im_lambda  = this->m_lambda.imag();
		const double mag_lambda = abs(this->m_lambda);
		this->m_Acp = (mag_lambda*mag_lambda-1)/(mag_lambda*mag_lambda+1);
		this->m_Scp = 2.*im_lambda/(mag_lambda*mag_lambda+1);
		this->m_rTSI = this->m_deltaTSI = 0;
		this->update_TSI();
	}

public:
	void clear_TSI(void)
	{
		this->m_rTSI = this->m_deltaTSI = 0;
		this->update_TSI();
	}

	void apply_TSI(const double phi1, const double phi3, const double rTSI, const double deltaTSI)
	{
		this->m_phi1 = phi1;
		this->m_phi3 = phi3;
		this->m_rTSI = rTSI;
		this->m_deltaTSI = deltaTSI;

		this->update_TSI();
	}


private:
	void update_TSI(void)
	{
		const double re_lambda  = this->m_lambda.real();
		const double im_lambda  = this->m_lambda.imag();
		const double mag_lambda = abs(this->m_lambda);

		// q=+1
		{
			const double coeff_cosh =  .5*(mag_lambda*mag_lambda+1) - 2*m_rTSI*re_lambda*::cos(2*m_phi1+m_phi3-m_deltaTSI);
			const double coeff_sinh = 0;
			const double coeff_cos  = +.5*(mag_lambda*mag_lambda-1) + 2*m_rTSI*im_lambda*::sin(2*m_phi1+m_phi3-m_deltaTSI);
			const double coeff_sin  = +im_lambda + m_rTSI*(1-mag_lambda*mag_lambda)*::sin(2*m_phi1+m_phi3-m_deltaTSI);

			this->m_coeff_qpos_cosh = coeff_cosh / coeff_cosh;
			this->m_coeff_qpos_sinh = coeff_sinh / coeff_cosh;
			this->m_coeff_qpos_cos  = coeff_cos  / coeff_cosh;
			this->m_coeff_qpos_sin  = coeff_sin  / coeff_cosh;
		}
		
		// q=-1
		{
			const double coeff_cosh =  .5*(mag_lambda*mag_lambda+1) - 2*m_rTSI*re_lambda*::cos(2*m_phi1+m_phi3+m_deltaTSI);
			const double coeff_sinh = 0;
			const double coeff_cos  = -.5*(mag_lambda*mag_lambda-1) + 2*m_rTSI*im_lambda*::sin(2*m_phi1+m_phi3+m_deltaTSI);
			const double coeff_sin  = -im_lambda + m_rTSI*(1-mag_lambda*mag_lambda)*::sin(2*m_phi1+m_phi3+m_deltaTSI);

			this->m_coeff_qneg_cosh = coeff_cosh / coeff_cosh;
			this->m_coeff_qneg_sinh = coeff_sinh / coeff_cosh;
			this->m_coeff_qneg_cos  = coeff_cos  / coeff_cosh;
			this->m_coeff_qneg_sin  = coeff_sin  / coeff_cosh;
		}
	}

public:
	double cosh(const int flv_asc) const { assert(abs(flv_asc)==1); return flv_asc>0 ? m_coeff_qpos_cosh : m_coeff_qneg_cosh;}
	double sinh(const int flv_asc) const { assert(abs(flv_asc)==1); return flv_asc>0 ? m_coeff_qpos_sinh : m_coeff_qneg_sinh;}
	double cos (const int flv_asc) const { assert(abs(flv_asc)==1); return flv_asc>0 ? m_coeff_qpos_cos  : m_coeff_qneg_cos; }
	double sin (const int flv_asc) const { assert(abs(flv_asc)==1); return flv_asc>0 ? m_coeff_qpos_sin  : m_coeff_qneg_sin; }


private:
	c_double m_lambda;
	double m_phi1, m_phi3, m_rTSI, m_deltaTSI; // all for the TSI effect, effectless when rTSI=0
	double m_Acp, m_Scp;


private:
	double m_coeff_qpos_cosh, m_coeff_qneg_cosh;
	double m_coeff_qpos_sinh, m_coeff_qneg_sinh;
	double m_coeff_qpos_cos,  m_coeff_qneg_cos;
	double m_coeff_qpos_sin,  m_coeff_qneg_sin;
};

