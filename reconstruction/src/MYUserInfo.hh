//***********************************************************
// MYUserInfo.hh            created 2013/01/16
//***********************************************************
#ifndef MYUSERINFO_HH
#define MYUSERINFO_HH

#include "belle.h"
#include "particle/ParticleUserInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  class MYUserInfo : public ParticleUserInfo {
  public:
    // double prob_kpi; // legacy of Nishida's example, not used.
  private:
    double m_ksmass;    // KS->pi+pi- mass before mass-vertex fit.
    double m_etackpichi2; // eta_c daughter Kpi vertex fit chi^2.
    double m_pi0prob; // gamma pi0 probability.
    double m_e9oe25; // gamma cluster E9/E25.
    double m_bkpichi2; // B daughter Kpi vertex fit chi^2.
  public:
    MYUserInfo();
    MYUserInfo( const MYUserInfo &x ); // copy constructor
    ~MYUserInfo() {}
    MYUserInfo *clone(void) const { return new MYUserInfo(*this); }
    MYUserInfo &operator = ( const MYUserInfo &x );
  public: // 
    void ksmass(double x) {m_ksmass = x;} // Set value.
    double ksmass() { return m_ksmass; }  // Return value.
    void etackpichi2(double x) { m_etackpichi2 = x; } 
    double etackpichi2() { return m_etackpichi2; }
    void pi0prob(double x){ m_pi0prob = x; }
    double pi0prob(){ return m_pi0prob; }
    void e9oe25(double x){ m_e9oe25 = x; }
    double e9oe25(){ return m_e9oe25; }
    void bkpichi2(double x){ m_bkpichi2 = x; }
    double bkpichi2(){ return m_bkpichi2; }
  };

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif
