//***************************************************************
// GammaUserInfo.cc     created 2019/02/19
//***************************************************************
#include "belle.h"
#include "GammaUserInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  //default constructor
  GammaUserInfo::GammaUserInfo(double pi0Veto, double etaVeto, double e9oe25) : m_pi0Veto(pi0Veto), m_etaVeto(etaVeto), m_e9oe25(e9oe25)
  { 
  }

  //copy constructor
  GammaUserInfo::GammaUserInfo( const GammaUserInfo &x ) : ParticleUserInfo(x) {
	m_pi0Veto = x.m_pi0Veto;
	m_etaVeto = x.m_etaVeto;
	m_e9oe25 = x.m_e9oe25;
  }

  // Operator "=".
  GammaUserInfo &GammaUserInfo::operator=( const GammaUserInfo &x ) {
	m_pi0Veto = x.m_pi0Veto;
	m_etaVeto = x.m_etaVeto;
	m_e9oe25 = x.m_e9oe25;
	return *this;
  }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
