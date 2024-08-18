//***************************************************************
// MYUserInfo.cc     created 2019/01/19 for B0->eta_c gamma K pi.
//***************************************************************
#include "belle.h"
#include "MYUserInfo.hh"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  //default constructor
  MYUserInfo::MYUserInfo() : m_ksmass(-1.0), m_etackpichi2(-1.0), 
			     m_pi0prob(-1.0), m_e9oe25(-1.0),
			     m_bkpichi2(-1.0){
  }

  //copy constructor
  MYUserInfo::MYUserInfo( const MYUserInfo &x ) : ParticleUserInfo(x) {
    m_ksmass = x.m_ksmass;
    m_etackpichi2 = x.m_etackpichi2;
    m_pi0prob = x.m_pi0prob;
    m_e9oe25 = x.m_e9oe25; 
    m_bkpichi2 = x.m_bkpichi2;
  }

  // Operator "=".
  MYUserInfo &MYUserInfo::operator=( const MYUserInfo &x ) {
    m_ksmass = x.m_ksmass;
    m_etackpichi2 = x.m_etackpichi2;
    m_pi0prob = x.m_pi0prob;
    m_e9oe25 = x.m_e9oe25; 
    m_bkpichi2 = x.m_bkpichi2;
    return *this;
  }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
