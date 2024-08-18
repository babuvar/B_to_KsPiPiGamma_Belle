#include "belle.h"
//
// P.Koppenburg - 2003 08 29
//

#include <iostream>
#include <iomanip>
#include <math.h>
#include "pi0eta_prob.h"
//#include "ana_utils.h"
//#include "properties.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
/*==================================================================================================
* Pi0 and eta prob
*=================================================================================================*/
double Pi0_Prob(double m, double p2, double theta){   // mass of pair, p2 of second gamma. First gamma assumed p1>1.5
  return Pi0_Eta_Prob(1, m, p2, theta);
}
/*================================================================================================*/
double Eta_Prob(double m, double p2, double theta){   // mass of pair, p2 of second gamma. First gamma assumed p1>1.5
  return Pi0_Eta_Prob(2, m, p2, theta);
}
/*================================================================================================*/
double Pi0_Eta_Prob(int what, double m, double p2, double theta){   // mass of pair, p2 of second gamma. First gamma assumed p1>1.5
  double p_min, p_max, m_min, m_max;
  int p_bins, m_bins ; 
  double *prob ;
  if (what==1){
#include "pi0_prob.h"
  } else if (what==2){
#include "eta_prob.h"
  } else  return -1.;
// get row
  if (p2<=0.) return -2;
  double logp = log(1000*p2)/log(10.);
  if ((logp<p_min) || (logp>p_max)) return 0;
  double drow = p_bins*(logp-p_min)/(p_max-p_min);
  int row = int(drow);
// get column
  if ((m<m_min) || (m>m_max)) return 0;
  double dcol = m_bins*(m-m_min)/(m_max-m_min);
  int col = int(dcol);
  if ((row<0) || (row>=p_bins) || (col<0) || (col>=m_bins)){
    std::cout << "Pi0 Eta prob: " << what << " " << logp << " " << m << " get coords " << row << " " << col << std::endl ;
    return -99.;
  }
// get first order prob
  int pos = m_bins*row + col ;
  if (pos>m_bins*p_bins) {
    std::cout << "Pi0 Eta prob: " << what << " " << logp << " " << m << " get coords " << row << " " << col << " i.e. " << pos << std::endl ;
    return -99.;
  }
  double fop = prob[pos] ;
//  std::cout << "Pi0 Eta prob: " << what << " " << logp<< " "  << m 
//      << " get coords " << row << " " << col << " i.e. " << pos << " -> prob " << fop << std::endl ;
  return fop ;
}

/*==================================================================================================
* Closest probability Interface
*=================================================================================================*/
double Closest_Pi0_Probability
( const Mdst_gamma& IN,Mdst_gamma& OUT,double& mm )
{
  return Closest_Probability(1,IN,OUT,mm,NULL);
}
/*================================================================================================*/
double Closest_Pi0_Probability
( const Mdst_gamma& IN,Mdst_gamma& OUT,double& mm,Mdst_gamma* VETO )
{
  return Closest_Probability(1,IN,OUT,mm,VETO);
}
/*================================================================================================*/
double Closest_Eta_Probability
( const Mdst_gamma& IN,Mdst_gamma& OUT,double& mm ){
  return Closest_Probability(2,IN,OUT,mm,NULL);
} 
/*================================================================================================*/
double Closest_Eta_Probability
( const Mdst_gamma& IN,Mdst_gamma& OUT,double& mm,Mdst_gamma* VETO ){
  return Closest_Probability(2,IN,OUT,mm,VETO);
}
/*==================================================================================================
* Closest probability
*=================================================================================================*/
double Closest_Probability
(int what, const Mdst_gamma& IN, Mdst_gamma& OUT, double& mm, Mdst_gamma* VETO)
{
  
  HepLorentzVector Reco1( IN.px(), IN.py(), IN.pz()  ) ;
  Reco1= HepLorentzVector( Reco1.px(), Reco1.py(), Reco1.pz(), Reco1.rho()) ;
  double bp = -100.;
  mm = -1;
  
  Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();
  
  for (std::vector<Mdst_gamma>::iterator it1= GammaMgr.begin(); it1 !=GammaMgr.end(); ++it1){
    Mdst_gamma* D = &(*it1);
    if (VETO!=NULL){ if (D->get_ID()==VETO->get_ID()) continue ;}
    if (D->get_ID()==IN.get_ID()) continue ;
    
    HepLorentzVector PG( D->px(), D->py(), D->pz()  ) ;
    PG = HepLorentzVector( PG.px(), PG.py(), PG.pz(), PG.rho()) ;
    double m = (Reco1 + PG).m();
    double p = Pi0_Eta_Prob(what, m, PG.e(), PG.theta() );
    
    if (p>bp) {
      bp = p ;
      mm = m ;
      OUT = *D ;
    }
  }  
  return bp;  
}
#if defined(BELLE_NAMESPACE)
}// namecpace Belle
#endif

