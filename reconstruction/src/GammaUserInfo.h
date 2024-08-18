//***********************************************************
// GammaUserInfo.hh            created 2019/02/19
//***********************************************************
#ifndef GAMMAUSERINFO_HH
#define GAMMAUSERINFO_HH

#include "belle.h"
#include "particle/ParticleUserInfo.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "ksfwmoments.h" 

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  class GammaUserInfo : public ParticleUserInfo {
  public:
  private:
	double m_pi0Veto;    //
	double m_etaVeto;    //
	double m_e9oe25;    //

  public:
    GammaUserInfo(double pi0Veto = -1, double etaVeto = -1, double e9oe25 = -1);
    GammaUserInfo( const GammaUserInfo &x ); // copy constructor
    ~GammaUserInfo() {}
    GammaUserInfo *clone(void) const { return new GammaUserInfo(*this); }
    GammaUserInfo &operator = ( const GammaUserInfo &x );
  public: // 
	double getPi0Veto() {return m_pi0Veto;} // Set value.
	double getEtaVeto() {return m_etaVeto;} // Set value.
	double getE9E25() {return m_e9oe25;} // Set value.
  };

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif
