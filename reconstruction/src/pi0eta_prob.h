//
// P.Koppenburg - 2003 08 29
//

/*
  2004.05.11 put under ushutil ( again after disk crash )
  const added in input Mdst_gammas
*/


#ifndef __PI0ETA_PROB__
#define __PI0ETA_PROB__
#include "belle.h"

#include "belle.h"
#include "panther/panther.h"
#include MDST_H
#include "belleCLHEP/Vector/LorentzVector.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
/*
    Probabilities given m     - the 2gamma mass, 
                        E    - the 2nd gamma E, 
	                theta - and the 2nd gamma polar angle
*/
double Pi0_Prob(double m, double E, double theta);             
double Eta_Prob(double m, double E, double theta);
/*
    Returns highest Probability for IN   - high energy photon 
                                    VETO - This Mdst_gamma is excluded from the search
                          returns:       - probability
			            OUT  - low energy photon
			            m    - mass for best combination
 */
double Closest_Pi0_Probability
( const Mdst_gamma& IN, Mdst_gamma& OUT, double& m );
double Closest_Pi0_Probability
( const Mdst_gamma& IN, Mdst_gamma& OUT, double& m, Mdst_gamma* VETO );
double Closest_Eta_Probability
( const Mdst_gamma& IN, Mdst_gamma& OUT, double& m ); 
double Closest_Eta_Probability
( const Mdst_gamma& IN, Mdst_gamma& OUT, double& m, Mdst_gamma* VETO ); 

// internal usage
double Pi0_Eta_Prob(int what, double m, double p2, double theta);
double Closest_Probability(int,const Mdst_gamma&, Mdst_gamma&, double&, Mdst_gamma* ); 

static const double ECL_F_THETA =  33 ; 
static const double ECL_B_THETA = 128 ;
static const double Degrees = 0.017453293 ;
#if defined(BELLE_NAMESPACE)
}// namecpace Belle
#endif

#endif
