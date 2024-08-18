//***********************************************************
// BUserInfo.hh            created 2019/02/19
//***********************************************************
#ifndef MYUSERINFO_HH
#define MYUSERINFO_HH

#include "belle.h"
#include "particle/ParticleUserInfo.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "ksfwmoments.h" 
#include "hamlet/Hamlet.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  class BUserInfo : public ParticleUserInfo {
  public:
  private:
	double m_b0cl;    // B0 flavor
	double m_b0chi2;    // B0 flavor
	double m_b0ndf;    // B0 flavor

	double m_tagCL;    // B0 flavor
	double m_tagChi2;    // B0 flavor
	double m_tagNdf;    // B0 flavor
	double m_tagChi2wo;    // B0 flavor
	double m_tagNdfwo;    // B0 flavor
	int m_tagNtrk;    // B0 flavor
	int m_tagLepton;    // B0 flavor

	double m_cosT;

	double m_b0qr;    // B0 flavor
	double m_w;
	double m_dw;
	double m_wtag;
	int m_flavor;
	HepPoint3D m_b0Vtx;
	HepSymMatrix m_b0VtxErr;
	HepPoint3D m_tagVtx;
	HepSymMatrix m_tagVtxErr;
	ksfwmoments m_ksfw;
	//Hamlet * m_hamlet;
  public:
    BUserInfo(double b0cl = -1.0, double b0chi2 = -1.0, double b0ndf = -1.0,
	double tagCL = -1, double tagChi2 = -1, double tagNdf = -1, int tagNtrk = -1, int taglepton = 0,
	double cosT = -1, 
	double b0q = -2, double w = -1, double dw = -2, double wtag = -2, int flavor = -2);
    BUserInfo( const BUserInfo &x ); // copy constructor
    ~BUserInfo() 
	{
		//delete m_hamlet;
		//m_hamlet = NULL;
	}
    BUserInfo *clone(void) const { return new BUserInfo(*this); }
    BUserInfo &operator = ( const BUserInfo &x );
  public: // 
	double getB0qr() const { return m_b0qr; }  // Return value.
	double getW() const {return m_w;} // Set value.
	double getWtag() const {return m_wtag;} // Set value.
	int getFlavor() const {return m_flavor;} // Set value.
	double getdW() const {return m_dw;} // Set value.

	double getCosT() const {return m_cosT;} // Set value.

	double getB0CL() const { return m_b0cl; }  // Return value.
	double getB0Chi2() const { return m_b0chi2; }  // Return value.
	double getB0Ndf() const { return m_b0ndf; }  // Return value.
	double getTagCL() const { return m_tagCL; }  // Return value.
	double getTagChi2() const { return m_tagChi2; }  // Return value.
	double getTagNdf() const { return m_tagNdf; }  // Return value.
	int getTagNtrk() const { return m_tagNtrk; }  // Return value.
	int getTagLepton() const { return m_tagLepton; }  // Return value.
	//void setHamlet(Hamlet * hamlet) {m_hamlet = hamlet;} // Set value.
	void setB0tagwo(double ndf, double chi2) 
	{
		m_tagChi2wo = chi2; 
		m_tagNdfwo = ndf;
	} // Set value.
	double getTagChi2wo() const { return m_tagChi2wo; }  // Return value.
	double getTagNdfwo() const { return m_tagNdfwo; }  // Return value.
	void setB0qr(double x) {m_b0qr = x;} // Set value.
	void setB0Vtx(HepPoint3D b0Vtx, HepSymMatrix b0VtxErr) {m_b0Vtx = b0Vtx; m_b0VtxErr = b0VtxErr;};
	void setTagVtx(HepPoint3D tagVtx, HepSymMatrix tagVtxErr) {m_tagVtx = tagVtx; m_tagVtxErr = tagVtxErr;};
 	HepPoint3D getB0Vtx() {return m_b0Vtx;} 
 	HepPoint3D getTagVtx() {return m_tagVtx;} 
	HepSymMatrix getB0VtxErr() {return m_b0VtxErr;}
	HepSymMatrix getTagVtxErr() {return m_tagVtxErr;}
	ksfwmoments getKSFW(){return m_ksfw;}
	void setKSFW(ksfwmoments ksfw) {m_ksfw = ksfw;}
  };

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif
