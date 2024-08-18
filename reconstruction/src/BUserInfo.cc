//***************************************************************
// BUserInfo.cc     created 2019/02/19
//***************************************************************
#include "belle.h"
#include "BUserInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  //default constructor
  BUserInfo::BUserInfo(double b0cl, double b0chi2, double b0ndf, 
	double tagCL, double tagChi2, double tagNdf, int tagNtrk, int tagLepton,
	double cosT, double b0q, double w, double dw, double wtag, int flavor) : 
		m_b0qr(b0q), m_w(w), m_dw(dw), m_wtag(wtag), m_flavor(flavor), 
		m_b0cl(b0cl), m_b0chi2(b0chi2), m_b0ndf(b0ndf),
		m_tagCL(tagCL), m_tagChi2(tagChi2),  m_tagNdf(tagNdf), m_tagNtrk(tagNtrk), m_tagLepton(tagLepton),
		m_cosT(cosT){ 
	//m_hamlet = NULL;
  }

  //copy constructor
  BUserInfo::BUserInfo( const BUserInfo &x ) : ParticleUserInfo(x) {
	m_b0qr = x.m_b0qr;
	m_b0cl = x.m_b0cl;
	m_b0chi2 = x.m_b0chi2;
	m_b0ndf = x.m_b0ndf;
	m_tagCL = x.m_tagCL;
	m_tagChi2 = x.m_tagChi2;
	m_tagNdf = x.m_tagNdf;
	m_tagChi2wo = x.m_tagChi2wo;
	m_tagNdfwo = x.m_tagNdfwo;
	m_tagNtrk = x.m_tagNtrk;
	m_tagLepton = x.m_tagLepton;
	m_b0Vtx = x.m_b0Vtx;
	m_b0VtxErr = x.m_b0VtxErr;
	m_tagVtx = x.m_tagVtx;
	m_tagVtxErr = x.m_tagVtxErr;
	m_ksfw = x.m_ksfw;
	m_cosT = x.m_cosT;
	m_w = x.m_w;
	m_wtag = x.m_wtag;
	m_dw = x.m_dw;
	m_flavor = x.m_flavor;
  }

  // Operator "=".
  BUserInfo &BUserInfo::operator=( const BUserInfo &x ) {
	m_b0qr = x.m_b0qr;
	m_b0cl = x.m_b0cl;
	m_b0chi2 = x.m_b0chi2;
	m_b0ndf = x.m_b0ndf;
	m_tagCL = x.m_tagCL;
	m_tagChi2 = x.m_tagChi2;
	m_tagNdf = x.m_tagNdf;
	m_tagChi2wo = x.m_tagChi2wo;
	m_tagNdfwo = x.m_tagNdfwo;
	m_tagNtrk = x.m_tagNtrk;
	m_tagLepton = x.m_tagLepton;
	m_b0Vtx = x.m_b0Vtx;
	m_b0VtxErr = x.m_b0VtxErr;
	m_tagVtx = x.m_tagVtx;
	m_tagVtxErr = x.m_tagVtxErr;
	m_ksfw = x.m_ksfw;
	m_cosT = x.m_cosT;
	m_w = x.m_w;
	m_wtag = x.m_wtag;
	m_dw = x.m_dw;
	m_flavor = x.m_flavor;
	return *this;
  }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
