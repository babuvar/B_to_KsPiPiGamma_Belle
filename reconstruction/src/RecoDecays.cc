//********************************************************************
// RecoDecays.cc          .....created 2019/01/10 Kenkichi Miya. & S. Bilokin IPHC
// pi0/eta probability calculation added. 
// cos(theta_B) and KSFW variables calculation added. 2019 Jan. 9th.
//*******************************************************************
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
// belle lib. after 2007
#include "belle.h"

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "eid/eid.h"        // to use eid class.
#include "mdst/Muid_mdst.h" // to use Muid_mdst class.
#include "mdst/mdst.h"      // to use Benergy func.
#include "benergy/BeamEnergy.h" // to use another BeamEnergy class.
#include "kid/atc_pid.h"    // to use atc_pid class.
#include "helix/Helix.h"    // to use helix class.
#include "particle/Particle.h"
#include "particle/Relation.h"
#include "particle/Momentum.h"
#include "particle/utility.h"
#include "particle/combination.h"
#include "ip/IpProfile.h"  // to use IP Profile, run-dep. IP information.
#include "kfitter/kvertexfitter.h" // to use vertexfitter.
#include "kfitter/kmassfitter.h" // to use massfitter.
#include "tagv/TagV.h"
#include "icpv_skim/utility.h" // to use interface routine for kfitter call.
#include "nisKsFinder/nisKsFinder.h" // to use nisKsFinder
#include "pi0eta_prob.h" // to use pi0/eta probability caclulation routine.
#include "geninfo.h" // MC Mathcing
#include "ksfwmoments.h" 
#include "hamlet/Hamlet.h"
#include "hamlet/Fbtag_NN1.h"
#include "hamlet/Fbtag_MultDimLikelihood0.h"
#include HEPEVT_H
#include MDST_H
#include BELLETDF_H
#include EVTCLS_H

#include "BUserInfo.h" // to use my own UserInfo class.
#include "GammaUserInfo.h" // to use my own UserInfo class.

#if defined(BELLE_NAMESPACE)

using namespace std;

namespace Belle{
#endif

//=====
// Copied from Trabelsi Karim's example.
// Track's helix parameter with respect to the run-dep. IP.
//=====
HepVector param_at_ip(Particle &p){

  const Mdst_charged charged(p.mdstCharged());

  double thisMass = p.mass();

  int hyp = 4;
  if(thisMass < 0.005){ hyp = 0; // e = 0.000511  
  }else if(thisMass < 0.110){ hyp = 1; // mu = 0.1056
  }else if(thisMass < 0.200){ hyp = 2; // pi = 0.13956
  }else if(thisMass < 0.5){ hyp = 3; // K = 0.4936
  }
  const HepPoint3D pivot(charged.trk().mhyp(hyp).pivot_x(),
			 charged.trk().mhyp(hyp).pivot_y(),
			 charged.trk().mhyp(hyp).pivot_z());
  
  HepVector  a(5);
  a[0] = charged.trk().mhyp(hyp).helix(0);
  a[1] = charged.trk().mhyp(hyp).helix(1);
  a[2] = charged.trk().mhyp(hyp).helix(2);
  a[3] = charged.trk().mhyp(hyp).helix(3);
  a[4] = charged.trk().mhyp(hyp).helix(4);
  HepSymMatrix Ea(5,0);
  Ea[0][0] = charged.trk().mhyp(hyp).error(0);
  Ea[1][0] = charged.trk().mhyp(hyp).error(1);
  Ea[1][1] = charged.trk().mhyp(hyp).error(2);
  Ea[2][0] = charged.trk().mhyp(hyp).error(3);
  Ea[2][1] = charged.trk().mhyp(hyp).error(4);
  Ea[2][2] = charged.trk().mhyp(hyp).error(5);
  Ea[3][0] = charged.trk().mhyp(hyp).error(6);
  Ea[3][1] = charged.trk().mhyp(hyp).error(7);
  Ea[3][2] = charged.trk().mhyp(hyp).error(8);
  Ea[3][3] = charged.trk().mhyp(hyp).error(9);
  Ea[4][0] = charged.trk().mhyp(hyp).error(10);
  Ea[4][1] = charged.trk().mhyp(hyp).error(11);
  Ea[4][2] = charged.trk().mhyp(hyp).error(12);
  Ea[4][3] = charged.trk().mhyp(hyp).error(13);
  Ea[4][4] = charged.trk().mhyp(hyp).error(14);
  Helix helix(pivot, a, Ea);

  const Hep3Vector&   IP     = IpProfile::position();
  if (IP.mag())
    helix.pivot(IP);
  return helix.a();
}

//======
// Set proper error matrix for gamma. Note errCart length is 4.
// Copied from FindGamma in icpv_skim package.
//======
void errGam( HepSymMatrix& errCart, Mdst_gamma gamma )
{//          ^^^^^^^^^^^^^ output , ^^^^^^^^^^ input

  HepSymMatrix  errEcl( 3, 0 ); // 3x3 initialize to zero
  errEcl[ 0 ][ 0 ] = gamma.ecl().error( 0 ); // Energy
  errEcl[ 1 ][ 0 ] = gamma.ecl().error( 1 );
  errEcl[ 1 ][ 1 ] = gamma.ecl().error( 2 ); // Phi
  errEcl[ 2 ][ 0 ] = gamma.ecl().error( 3 );
  errEcl[ 2 ][ 1 ] = gamma.ecl().error( 4 );
  errEcl[ 2 ][ 2 ] = gamma.ecl().error( 5 ); // Theta

  HepMatrix  jacobian( 4, 3, 0 );
  double  cp = cos( gamma.ecl().phi() );
  double  sp = sin( gamma.ecl().phi() );
  double  ct = cos( gamma.ecl().theta() );
  double  st = sin( gamma.ecl().theta() );
  double   E = gamma.ecl().energy();

  jacobian[ 0 ][ 0 ] =       cp * st;
  jacobian[ 0 ][ 1 ] =  -E * sp * st;
  jacobian[ 0 ][ 2 ] =   E * cp * ct;
  jacobian[ 1 ][ 0 ] =       sp * st;
  jacobian[ 1 ][ 1 ] =   E * cp * st;
  jacobian[ 1 ][ 2 ] =   E * sp * ct;
  jacobian[ 2 ][ 0 ] =            ct;
  jacobian[ 2 ][ 1 ] =           0.0;
  jacobian[ 2 ][ 2 ] =  -E      * st;
  jacobian[ 3 ][ 0 ] =           1.0;
  jacobian[ 3 ][ 1 ] =           0.0;
  jacobian[ 3 ][ 2 ] =           0.0;
  errCart = errEcl.similarity( jacobian );
  return;
}

bool sortBFunction(const Particle & b1, const Particle & b2)
{
	const BUserInfo& binfo1 = dynamic_cast<const BUserInfo&>(b1.userInfo());	
	const BUserInfo& binfo2 = dynamic_cast<const BUserInfo&>(b2.userInfo());
	if (binfo1.getB0CL() > binfo2.getB0CL())
	{
		return true;
	}
	return false;
}

//===
// Select gamma exceeding proper threshold.
//===
int gamma_tight( Hep3Vector gamma_3v )
{
   // Energy threshold of gamma shower in GeV.
   const double eth_barrel = 0.05;
   const double eth_fwdec  = 0.10;
   const double eth_bwdec  = 0.10;

   // Gaps between barrel and endcaps in cos(theta).
   const double cos_fwd_gap =  0.8501;
   const double cos_bwd_gap = -0.6643;
   
   int ireturn = -1; // return value, -1:bad, 0:good.

   if( gamma_3v.cosTheta() > cos_fwd_gap )// In Fwd. Endcap.
   {
      if( gamma_3v.mag() > eth_fwdec )
      {
	 ireturn = 0;
      }
   }
   else if( gamma_3v.cosTheta() < cos_fwd_gap 
	    && gamma_3v.cosTheta() > cos_bwd_gap )// In Barrel.
   {
      if( gamma_3v.mag() > eth_barrel )
      {
	 ireturn = 0;
      }
   }
   else // In Bwd. Endcap.
   {
      if( gamma_3v.mag() > eth_bwdec )
      {
	 ireturn = 0;
      }
   }
   return ireturn;
}
// From Miyabayashi-san, corrected
// Creates GammaUserInfo objects and attaches them to corresponding photons
// Runs pi0 eta veto and computes e9e25
void fillPhotonProperties(std::vector<Particle> & gammas, double her, double ler, double cross)
{
	Mdst_ecl_aux_Manager& mdstauxmgr = Mdst_ecl_aux_Manager::get_manager();
	unsigned nPhotons = gammas.size();
	for (unsigned int i = 0; i < gammas.size(); i++)
	{
		double pi0_prob = -1.0;
		double eta_prob = -1.0;
		double temp_pi0_prob = -1.0;
		double temp_eta_prob = -1.0;
		double masspi0=0.0;
		double masseta=0.0;
		double e9oe25 = -1.0;
		if (gamma_tight( gammas[i].p().vect() ) != 0)
		{
			continue;
		}
		for (unsigned int j = 0; j < gammas.size(); j++)
		{
			if (i == j || gammas[j].p().e() < 0.02 ||  gammas[j].p().e() > 5.0)
			{
				continue;
			}
			HepLorentzVector gamgam_lv = gammas[i].p() + gammas[j].p();
			double m_gamgam = gamgam_lv.mag();
			if( 0.034976 < m_gamgam && m_gamgam < 0.234976 )
			{
				temp_pi0_prob = Pi0_Prob( m_gamgam, gammas[j].p().e(), gammas[j].p().theta());
				if( temp_pi0_prob > pi0_prob )
				{
					pi0_prob = temp_pi0_prob; 
					masspi0 = m_gamgam;
				}
			}
			if( 0.4473 < m_gamgam && m_gamgam < 0.6473 )
			{
				temp_eta_prob = Eta_Prob( m_gamgam, gammas[j].p().e(), gammas[j].p().theta());
				if( temp_eta_prob > eta_prob )
				{
					eta_prob = temp_eta_prob;
					masseta = m_gamgam;
				}
			}
		}
		int id_ecl =gammas[i].mdstGamma().ecl().get_ID();
		std::vector<Mdst_ecl_aux>::iterator itaux;
		for(itaux=mdstauxmgr.begin(); itaux!=mdstauxmgr.end(); itaux++)
		{
			Mdst_ecl_aux& aux = *itaux;
			if(id_ecl == aux.get_ID()) e9oe25 = aux.e9oe25();
		}

		GammaUserInfo GUI(pi0_prob, eta_prob, e9oe25);
		gammas[i].userInfo(GUI);
	}
}

// Shamelessly stolen from Luka sensei
double getThrustAngle(Particle& brec, double her, double ler, double cross){

  Mdst_charged_Manager &chgmgr = Mdst_charged_Manager::get_manager();
  Mdst_gamma_Manager &gammgr = Mdst_gamma_Manager::get_manager();
  std::vector<Particle> vec_asc_track;
  std::vector<Particle> vec_rec_track;

        // make used track list by rec-side
        unsigned int list_pi_mask = 0;
        unsigned int list_gam_mask = 0;  
        for( int i=0; i<brec.relation().nFinalStateParticles(); i++ ){
                const Particle &p(brec.relation().finalStateParticle(i));
               	if(p.mdstCharged()) list_pi_mask |= 1<<(int)p.mdstCharged().get_ID();
               	if(p.mdstGamma()) list_gam_mask |= 1<<(int)p.mdstGamma().get_ID();
                vec_rec_track.push_back(p);
        }

        // fill other tracks
        for( std::vector<Mdst_charged>::iterator chg=chgmgr.begin(); chg!=chgmgr.end(); chg++ ){
                // skip rec-side tracks...
               	if( 1<<(int)chg->get_ID() & list_pi_mask ) continue;
                const char *pname = chg->charge()>0 ? "PI+" : "PI-";
                Particle p(*chg,std::string(pname));
                vec_asc_track.push_back(p);
        }
        for( std::vector<Mdst_gamma>::iterator gam=gammgr.begin(); gam!=gammgr.end(); gam++ ){
               	// skip rec-side tracks...
                if( 1<<(int)gam->get_ID() & list_gam_mask ) continue;
                Particle p(*gam);
                vec_asc_track.push_back(p);
        }

        Hep3Vector bthrust = calcuThrust(vec_rec_track, her, ler, cross);
        Hep3Vector rodthrust = calcuThrust(vec_asc_track, her, ler, cross);      
        double tht = bthrust.angle(rodthrust);    
  
	return tht;
}
// Get information about event: run number, event number, experiment number
void get_event_id(int &no_exp, int &no_run, int &no_evt )
{
	no_exp=-1, no_run=-1, no_evt=-1;

	belle_event *belle_event;
	belle_event = (struct belle_event*)BsGetEnt(BELLE_EVENT,1,BBS_No_Index);
	if( belle_event )
	{
		no_exp = belle_event->m_ExpNo;
		no_run = belle_event->m_RunNo;
		no_evt = belle_event->m_EvtNo & 0x0fffffff;
	}
	return;
}
inline const double calc_mbc(const Particle &bcand)
{
	
	const double cross = BeamEnergy::Cross_angle();//0.022;  // Finite angle crossing, in rad. 
	const double her = BeamEnergy::E_HER();//7.998213; // HER beam energy.(New 2003)
	const double ler = BeamEnergy::E_LER();//3.499218; // LER beam energy.(New 2003)
	const double benergy = BeamEnergy::Ecm() / 2.;
	const double mass_e = 0.000510998902;

	Vector4 p4_her(0.0, 0.0, +sqrt(her*her-mass_e*mass_e), her); p4_her.rotateY(cross);
	Vector4 p4_ler(0.0, 0.0, -sqrt(ler*ler-mass_e*mass_e), ler);

	Vector4 lv = p4_her + p4_ler;
	Vector4 b_momentum(bcand.p());

	b_momentum.boost(-lv.boostVector());
	b_momentum.setT(benergy);

	return b_momentum.mag();
}


// Function which returns generated value of delta t and delta z
std::vector<double> getDtDz(Gen_hepevt_Manager & gen_mgr, int & genBtagFlavor, bool runOnData = true)
{
	//--------------------------
	//other true information
	//for decay time study
	//--------------------------
	std::vector<double> vect_dtdz;
	double zrec_gen = 0;
	double ztag_gen = 0;
	unsigned counter=0;
	double dz_gen=0;
	double dt_gen=-11;
	double trec_gen = 0;
	double ttag_gen = 0;
	double boost_trec_gen = 0;
	double boost_ttag_gen = 0;
	if (runOnData)
	{
		vect_dtdz.push_back(dt_gen);
		vect_dtdz.push_back(trec_gen);
		vect_dtdz.push_back(ttag_gen);
		vect_dtdz.push_back(dz_gen);
		vect_dtdz.push_back(zrec_gen);
		vect_dtdz.push_back(ttag_gen);
		return vect_dtdz;
	}
	int resonancePDG = 10313;
	int resonancePDG2 = 30343;
	int photonPDG = 22;
	int B_counter = 0;
	for (std::vector<Gen_hepevt>::iterator j = gen_mgr.begin(); j != gen_mgr.end(); ++j)
	{
		

		if (abs((*j).idhep()) == 511)  //511 for B0, 521 for B+, 531 for Bs
		{
			Gen_hepevt &child1 = *(gen_mgr.begin() - 1 + (*j).daFirst());
			Gen_hepevt &child2 = *(gen_mgr.begin() - 1 + (*j).daLast());
			// 111 for pi0 and 443 for jpsi
			if ( (((abs(child1.idhep()) == resonancePDG) && (child2.idhep() == photonPDG))||
				 ((child1.idhep() == photonPDG) && (abs(child2.idhep()) == resonancePDG))) || 
				 (((abs(child1.idhep()) == resonancePDG2) && (child2.idhep() == photonPDG))||
				 ((child1.idhep() == photonPDG) && (abs(child2.idhep()) == resonancePDG2))) )
			{
				counter = ++counter;
				if (counter == 1)
				{
					zrec_gen = child1.VZ();
					trec_gen = child1.T();
					double boost_gamma = (*j).E()/(*j).M();
					boost_trec_gen = child1.T()/boost_gamma;
				}
			}
			else if ( (!((abs(child1.idhep()) == resonancePDG) && (child2.idhep() == photonPDG))||
				 	!((child1.idhep() == photonPDG) && (abs(child2.idhep()) == resonancePDG))) || 
				 	(!((abs(child1.idhep()) == resonancePDG2) && (child2.idhep() == photonPDG))||
				 	!((child1.idhep() == photonPDG) && (abs(child2.idhep()) == resonancePDG2))) ||
					 counter == 2)
			{
				genBtagFlavor = (*j).idhep() / abs((*j).idhep());
				ztag_gen = child1.VZ();
				ttag_gen = child1.T();
				double boost_gamma_tag = (*j).E()/(*j).M();
				boost_ttag_gen = child1.T()/boost_gamma_tag;
			}
			double zrecdiff=0;
			double ztagdiff=0;
			dz_gen = (zrec_gen-ztag_gen)/10;   //in cm, z is in mm
			double light_speed=2.99792458;
			dt_gen = ((boost_trec_gen-boost_ttag_gen)/light_speed)*10;//in ps , boosted
		}
	}
	vect_dtdz.push_back(dt_gen);
	vect_dtdz.push_back(trec_gen/10);
	vect_dtdz.push_back(ttag_gen/10);
	vect_dtdz.push_back(dz_gen);
	vect_dtdz.push_back(zrec_gen/10);
	vect_dtdz.push_back(ztag_gen/10);
 	return vect_dtdz;
}

bool enough_svd_hit(const Particle &p)
{
	Mdst_trk_fit &trk_fit = p.mdstCharged().trk().mhyp(2);
	if( trk_fit.nhits(3)<1 ) return false;
	if( trk_fit.nhits(4)<2 ) return false;
	return true;
}


struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

// Check if the particle HeVMgr[mo_id] is a mo_pdg, and if it's daugther are the one specified in vect_true_da_pdg + any number of gamma
// Don't work if pdg in vect_true_da_pdg < 22 = gamma_pdg
bool good_decay(std::vector<Gen_hepevt> HevMgr, int mo_id, int mo_pdg, std::vector<int> vect_true_da_pdg)
{
	std::vector<int> vect_da_pdg;

	// Sort the list
	std::sort(vect_true_da_pdg.begin(), vect_true_da_pdg.end(), greater());

	if(abs(HevMgr[mo_id].idhep())==abs(mo_pdg))
	{
		for(int i=HevMgr[mo_id].daFirst(); i<=HevMgr[mo_id].daLast(); i++)
		{
			vect_da_pdg.push_back(abs(HevMgr[i-1].idhep()));
		}

		// Sort the list, gamma will be at the end as gamma_pdg=22
		std::sort(vect_da_pdg.begin(), vect_da_pdg.end(), greater());
		
		// std::cout<<vect_true_da_pdg[0]<<" "<<vect_true_da_pdg[0]<<std::endl;
		// for(int i=0;i<vect_da_pdg.size();i++)
		// {
		// 	std::cout<<vect_da_pdg[i]<<" ";
		// }
		// std::cout<<std::endl;
		// std::cout<<std::endl;

		// Check the pdg of the 2 first daughter
		for(int i=0; i<vect_true_da_pdg.size(); i++)
		{
			if(vect_true_da_pdg[i]!=vect_da_pdg[i])
			{
				return false;
			}
		}

		for(int i=vect_true_da_pdg.size(); i<vect_da_pdg.size(); i++)
		{
			// Check if the other daughters are gammas
			if(vect_da_pdg[i]!=22)
			{
				return false;
			}
		}

		return true;

	}
	else{return false;}
}
std::vector<int> loop_on_daugthers(std::vector<Gen_hepevt> HevMgr, int mo_id, std::vector<int> vect_true_da_pdg, std::vector<int> vect_da_pdg)
{
	std::cout<<"INFO Looping over daugthers of particle "<<HevMgr[mo_id].idhep()<<std::endl;
	// If particle id is zero put a zero directly on the temporal list
	if (HevMgr[mo_id].idhep()==0)
	{
	   vect_da_pdg.push_back(0);
	   return vect_da_pdg;
	}
	// define a list of final state particles to finish the loop with, 
	std::vector<int> vect_final_state_pdg = {
            11, 12, 13, 14, 15, 16, // leptons
            22, //  photons
            111, 211, // pions
            321, //  K+
            310, 130, // Kshort, Klong
            2112, 2212 }; // neutrons, protons
	
	// If the daugthers particles are 0, put a zero and stop the loop
	if (HevMgr[mo_id].daFirst()==0)
	{
	   vect_da_pdg.push_back(0);
	   return vect_da_pdg;
	}
	
	for(int i=HevMgr[mo_id].daFirst(); i<=HevMgr[mo_id].daLast(); i++)
		{
			// Set that a daugther is not a final state
			bool fstatepart = false;
			std::cout<<"INFO Working on "<<HevMgr[i-1].idhep()<<std::endl;
			
			// Loop over the expected final particles
			for(int j=0; j<vect_final_state_pdg.size(); j++)
			{
				// Check that the daugther matches one of the expected final particles
				if (vect_final_state_pdg[j] == abs(HevMgr[i-1].idhep()))
				{
					//  if its a photon and either its energy is very low or it comes from a particle that is not a B0/B0bar 
					if(HevMgr[i-1].idhep()==22 && (HevMgr[i-1].E() < 0.05 || abs(HevMgr[mo_id].idhep())!= 511) )
					{
						// Mark the particle as a final state particle but it is not added to the list
						// std::cout<<"INFO No signal photon "<<HevMgr[i-1].idhep() <<" with energy "<<HevMgr[i-1].E()<<" and the mother is "<< HevMgr[mo_id].idhep()<<std::endl;
						fstatepart = true;
						continue;
					}
					
					// Add the particle to the list
					vect_da_pdg.push_back(HevMgr[i-1].idhep());
					std::cout<<"INFO Added daugther to the list "<<HevMgr[i-1].idhep()<<std::endl;
					fstatepart = true;
				}
			}
			// Loop over the final particle list, if the daugther is a final state particle stop checking

			// If the particle is not a final state particle, check for its daugthers
			if(fstatepart != true)
			{
				vect_da_pdg = loop_on_daugthers(HevMgr, i-1, vect_true_da_pdg, vect_da_pdg);
			}
		}
	std::cout<<"INFO Ended loop over daugthers of particle "<<HevMgr[mo_id].idhep()<<std::endl;
	return vect_da_pdg;
}

bool finalstate(std::vector<Gen_hepevt> HevMgr, int mo_id, std::vector<int> vect_true_da_pdg)
{
	// Check that the variable is a B/B0, otherwise leave
	if(abs(HevMgr[mo_id].idhep())!=abs(511))
	{
		std::cout<<"INFO Particle is not a b0/b0bar "<<HevMgr[mo_id].idhep()<<std::endl;
		return false;
	}

	// Create vector to store temporal particles
	std::vector<int> vect_temporal_da_pdg;
	// Fill up the vector with function loop on daugthers
	vect_temporal_da_pdg= loop_on_daugthers(HevMgr, mo_id, vect_true_da_pdg, vect_temporal_da_pdg);

	// Order both vectors to match them later
	std::sort(vect_temporal_da_pdg.begin(), vect_temporal_da_pdg.end(), greater());
	std::sort(vect_true_da_pdg.begin(), vect_true_da_pdg.end(), greater());

	//  print the list
	for (int i = 0; i < vect_temporal_da_pdg.size(); ++i) {
        std::cout << vect_temporal_da_pdg[i] << " ";
	}

	// Check that the size of the temporal vector is the same as the size of the target vector
	if(vect_true_da_pdg.size() != vect_temporal_da_pdg.size())
	{
		return false;
	}

	// check that the i element of a the true vector matches with the i element of the vector obtained from the loop on daugthers
	for(int i=0; i<vect_temporal_da_pdg.size(); i++)
		{
			if(vect_true_da_pdg[i]!=vect_temporal_da_pdg[i])
			{
				return false;
			}
		}
	
	return true;

}

// The is_signal_finalstate checks matches final state particles to a common mother independently of the Kres in the middle
// uses finalstate

int is_signal_finalstate(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// If finalstate returns true, meaning that the final state particles of a B matches the list {}, then the variable is set to 1. 
				// This takes into account possible brem/photos photons
				if(finalstate(HevMgr, i-1, {310, 211, -211, 22}) == true)
				{
					f_signal+=1;
					std::cout<<"INFO This triggered a yes "<<f_signal<<std::endl;
				}
			}

		}
	
	}
	return f_signal;
}

int is_signal_mcgen(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> rho ks0
					if(good_decay(HevMgr, k1_id, 10313, {113, 310})==true || good_decay(HevMgr, k1_id, 30343, {113, 310})==true)
					{
						// Order is always the same in HevMgr
						int rho_id = HevMgr[k1_id].daFirst()-1;
						int ks_id = HevMgr[k1_id].daFirst();

						// Check rho -> pi pi
						if(good_decay(HevMgr, rho_id, 113, {211, 211})==true)
						{
							// Check Ks0 -> pi pi
							if(good_decay(HevMgr, ks_id, 310, {211, 211})==true)
							{
								f_signal+=1;
							}
						}
					}
				}
			}
		}
		// if(f_signal==0)
		// {
		// 	std::cout<<"##########################################"<<std::endl;
		// 	std::cout<<"get_ID() \t= "<<i_hep->get_ID()<<std::endl;
		// 	std::cout<<"abs(idhep) \t= "<<abs(i_hep->idhep())<<std::endl;
		// 	std::cout<<"daFirst() \t= "<<i_hep->daFirst()<<std::endl;
		// 	std::cout<<"daLast() \t= "<<i_hep->daLast()<<std::endl;
		// }
		// if(f_signal==0)
		// {
		// 	std::cout<<i_hep->get_ID()-1<<" "<<i_hep->idhep()<<" "<<i_hep->daFirst()-1<<" "<<i_hep->daLast()-1<<std::endl;
		// }
	}
	return f_signal;
}


int is_signal_k1ks(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> K*0 pi0
					if(good_decay(HevMgr, k1_id, 10313, {111, 313})==true || good_decay(HevMgr, k1_id, 30343, {111, 313})==true)
					{
						int K_id = -999;
						int id1 = HevMgr[k1_id].daFirst()-1;
						int id2 = HevMgr[k1_id].daFirst();

						if(abs(HevMgr[id1].idhep())==313){K_id = id1;}
						else{K_id = id2;}

						// Check K* -> K0 pi0
						if(good_decay(HevMgr, K_id, 313, {311, 111})==true)
						{
							f_signal+=1;
						}
					}

					// Check K1 -> K*+ pi+
					if(good_decay(HevMgr, k1_id, 10313, {211, 323})==true || good_decay(HevMgr, k1_id, 30343, {211, 323})==true)
					{
						int K_id = -999;
						int id1 = HevMgr[k1_id].daFirst()-1;
						int id2 = HevMgr[k1_id].daFirst();

						if(abs(HevMgr[id1].idhep())==323){K_id = id1;}
						else{K_id = id2;}

						// Check K*+ -> K0 pi+
						if(good_decay(HevMgr, K_id, 323, {311, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k1k0s(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> K*0 pi0
					if(good_decay(HevMgr, k1_id, 10313, {10311, 111})==true || good_decay(HevMgr, k1_id, 30343, {10311, 111})==true)
					{
						int K_id = -999;
						int id1 = HevMgr[k1_id].daFirst()-1;
						int id2 = HevMgr[k1_id].daFirst();

						if(abs(HevMgr[id1].idhep())==10311){K_id = id1;}
						else{K_id = id2;}

						// Check K*0 -> K0 pi0
						if(good_decay(HevMgr, K_id, 10311, {311, 111})==true)
						{
							f_signal+=1;
						}
					}

					// Check K1 -> K*0+ pi+
					if(good_decay(HevMgr, k1_id, 10313, {10321, 211})==true || good_decay(HevMgr, k1_id, 30343, {10321, 211})==true)
					{
						int K_id = -999;
						int id1 = HevMgr[k1_id].daFirst()-1;
						int id2 = HevMgr[k1_id].daFirst();

						if(abs(HevMgr[id1].idhep())==10321){K_id = id1;}
						else{K_id = id2;}

						// Check K*0+ -> K0 pi+
						if(good_decay(HevMgr, K_id, 10321, {311, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k1k2p(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k2_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> K20 pi
					if(good_decay(HevMgr, k2_id, 10313, {325, 211})==true || good_decay(HevMgr, k2_id, 30343, {325, 211})==true)
					{
						int K_id = -999;
						int id1 = HevMgr[k2_id].daFirst()-1;
						int id2 = HevMgr[k2_id].daFirst();

						if(abs(HevMgr[id1].idhep())==325){K_id = id1;}
						else{K_id = id2;}

						// Check K20 -> K0 pi
						if(good_decay(HevMgr, K_id, 325, {311, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k1rho(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> rho0 k0
					if(good_decay(HevMgr, k1_id, 10313, {113, 311})==true || good_decay(HevMgr, k1_id, 30343, {113, 311})==true)
					{
						// Order is always the same in HevMgr
						int rho_id = HevMgr[k1_id].daFirst()-1;
						int k0_id = HevMgr[k1_id].daFirst();

						// Check rho -> pi pi
						if(good_decay(HevMgr, rho_id, 113, {211, 211})==true)
						{
							f_signal+=1;
						}
					}

					// Check K1 -> rho+ k+
					if(good_decay(HevMgr, k1_id, 10313, {213, 321})==true || good_decay(HevMgr, k1_id, 30343, {213, 321})==true)
					{
						// Order is always the same in HevMgr
						int rho_id = HevMgr[k1_id].daFirst()-1;
						int k0_id = HevMgr[k1_id].daFirst();

						// Check rho -> pi pi
						if(good_decay(HevMgr, rho_id, 213, {211, 111})==true || good_decay(HevMgr, rho_id, 213, {111, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k2ks(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K2 gamma
				if(good_decay(HevMgr, i-1, 511, {315, 22})==true)
				{
					int k2_id = HevMgr[i-1].daFirst()-1;

					// Check K2 -> K* pi
					if(good_decay(HevMgr, k2_id, 315, {211, 323})==true)
					{
						int K_id = -999;
						int id1 = HevMgr[k2_id].daFirst()-1;
						int id2 = HevMgr[k2_id].daFirst();

						if(abs(HevMgr[id1].idhep())==323){K_id = id1;}
						else{K_id = id2;}

						// Check K* -> K0 pi
						if(good_decay(HevMgr, K_id, 323, {311, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k2rho(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K2 gamma
				if(good_decay(HevMgr, i-1, 511, {315, 22})==true)
				{
					int k2_id = HevMgr[i-1].daFirst()-1;

					// Check K2 -> rho k0
					if(good_decay(HevMgr, k2_id, 315, {113, 311})==true)
					{
						// Order is always the same in HevMgr
						int rho_id = HevMgr[k2_id].daFirst()-1;
						int k0_id = HevMgr[k2_id].daFirst();

						// Check rho -> pi pi
						if(good_decay(HevMgr, rho_id, 113, {211, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k1kpp(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> K0 pi pi
					if(good_decay(HevMgr, k1_id, 10313, {311, 211, 211})==true || good_decay(HevMgr, k1_id, 30343, {311, 211, 211})==true)
					{
						f_signal+=1;
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k1omg(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> omega k0
					if(good_decay(HevMgr, k1_id, 10313, {223, 311})==true || good_decay(HevMgr, k1_id, 30343, {223, 311})==true)
					{
						int omega_id = HevMgr[k1_id].daFirst()-1;
						int k0_id = HevMgr[k1_id].daFirst();

						// Check omega -> pi pi
						if(good_decay(HevMgr, omega_id, 223, {211, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}


int is_signal_k1f0(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0 -> K1 gamma
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1 -> f0 k0
					if(good_decay(HevMgr, k1_id, 10313, {10221, 311})==true || good_decay(HevMgr, k1_id, 30343, {10221, 311})==true)
					{
						int f0_id = HevMgr[k1_id].daFirst()-1;
						int k0_id = HevMgr[k1_id].daFirst();

						// Check f0 -> pi pi
						if(good_decay(HevMgr, f0_id, 10221, {211, 211})==true)
						{
							f_signal+=1;
						}
					}
				}
			}
		}
	}
	return f_signal;
}




void coutevtgen(Gen_hepevt_Manager & HevMgr)
{
	std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// std::cout<<"##########################################"<<std::endl;
		// std::cout<<"get_ID() \t= "<<i_hep->get_ID()<<std::endl;
		// std::cout<<"abs(idhep) \t= "<<abs(i_hep->idhep())<<std::endl;
		// std::cout<<"daFirst() \t= "<<i_hep->daFirst()<<std::endl;
		// std::cout<<"daLast() \t= "<<i_hep->daLast()<<std::endl;

		std::cout<<i_hep->get_ID()-1<<" "<<i_hep->idhep()<<" "<<i_hep->daFirst()-1<<" "<<i_hep->daLast()-1<<std::endl;
	}
}

int is_signal_no_gamma_res(Gen_hepevt_Manager & HevMgr)
{
	int f_signal = 0;
	int fsig = 0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Check if B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0
				if(abs(HevMgr[i-1].idhep())==511 && (HevMgr[i-1].daLast()-HevMgr[i-1].daFirst()+1)==2)
				{
					// Check if K1 & gamma
					int da0 = HevMgr[i-1].daFirst()-1;
					int da1 = HevMgr[i-1].daLast()-1;

					if(((abs(HevMgr[da0].idhep())==10313 && abs(HevMgr[da1].idhep())==22) ||
					   (abs(HevMgr[da1].idhep())==10313 && abs(HevMgr[da0].idhep())==22)) ||
					   ((abs(HevMgr[da0].idhep())==30343 && abs(HevMgr[da1].idhep())==22) ||
					   (abs(HevMgr[da1].idhep())==30343 && abs(HevMgr[da0].idhep())==22)))
					{
						for(int j=da0; j<=da1; j++)
						{
							// Check K1
							if((abs(HevMgr[j].idhep())==10313 && (HevMgr[j].daLast()-HevMgr[j].daFirst()+1)==2) ||
								(abs(HevMgr[j].idhep())==30343 && (HevMgr[j].daLast()-HevMgr[j].daFirst()+1)==2))
							{
								int da00 = HevMgr[j].daFirst()-1;
								int da11 = HevMgr[j].daLast()-1;

								// Check if rho & Ks
								if(abs(HevMgr[da00].idhep())==113 && abs(HevMgr[da11].idhep())==310 ||
					   			abs(HevMgr[da11].idhep())==310 && abs(HevMgr[da00].idhep())==113)
								{
									for(int k=da00; k<=da11; k++)
									{
										// Check rho
										if(abs(HevMgr[k].idhep())==113 && (HevMgr[k].daLast()-HevMgr[k].daFirst()+1)==2)
										{
											int da000 = HevMgr[k].daFirst()-1;
											int da111 = HevMgr[k].daLast()-1;

											// Check pi+ pi-
											if(abs(HevMgr[da000].idhep())==211 && abs(HevMgr[da111].idhep())==211)
											{
												fsig+=1;
											} 
										}
										// Check Ks0
										if(abs(HevMgr[k].idhep())==310 && (HevMgr[k].daLast()-HevMgr[k].daFirst()+1)==2)
										{
											int da000 = HevMgr[k].daFirst()-1;
											int da111 = HevMgr[k].daLast()-1;

											// Check pi+ pi-
											if(abs(HevMgr[da000].idhep())==211 && abs(HevMgr[da111].idhep())==211)
											{
												fsig+=1;
											} 
										}
									}
								}
							}
						}
					}
				}
			}
		}
		// {
		// 	std::cout<<"FSIGNAL = "<<f_signal<<std::endl;
		// 	std::cout<<"##########################################"<<std::endl;
		// 	std::cout<<"get_ID() \t= "<<i_hep->get_ID()<<std::endl;
		// 	std::cout<<"abs(idhep) \t= "<<abs(i_hep->idhep())<<std::endl;
		// 	std::cout<<"daFirst() \t= "<<i_hep->daFirst()<<std::endl;
		// 	std::cout<<"daLast() \t= "<<i_hep->daLast()<<std::endl;
		// }
	}
	if(fsig==2 || fsig==4)
	{
		f_signal+=1;
	}
	return f_signal;
}


std::vector<Particle> get_pions_FS(Gen_hepevt_Manager & HevMgr)
{
	std::vector<Particle> vect_pi(4);
	int filled=0;
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// Check B0
				if(good_decay(HevMgr, i-1, 511, {10313, 22})==true || good_decay(HevMgr, i-1, 511, {30343, 22})==true)
				{
					int k1_id = HevMgr[i-1].daFirst()-1;

					// Check K1
					if(good_decay(HevMgr, k1_id, 10313, {113, 310})==true || good_decay(HevMgr, k1_id, 30343, {113, 310})==true)
					{
						int rho_id = HevMgr[k1_id].daFirst()-1;
						int ks_id = HevMgr[k1_id].daFirst();

						// Check rho
						if(good_decay(HevMgr, rho_id, 113, {211, 211})==true)
						{
							// Check Ks0
							if(good_decay(HevMgr, ks_id, 310, {211, 211})==true)
							{
								int pi1_id = HevMgr[rho_id].daFirst()-1;
								int pi2_id = HevMgr[rho_id].daFirst();
								int pi3_id = HevMgr[ks_id].daFirst()-1;
								int pi4_id = HevMgr[ks_id].daFirst();

								vect_pi[0]=HevMgr[pi1_id];
								vect_pi[1]=HevMgr[pi2_id];
								vect_pi[2]=HevMgr[pi3_id];
								vect_pi[3]=HevMgr[pi4_id];

								filled=1;
							}
						}
					}
				}
			}
		}
		// if(filled==1)
		// {
		// 	std::cout<<"##########################################"<<std::endl;
		// 	std::cout<<"get_ID() \t= "<<i_hep->get_ID()<<std::endl;
		// 	std::cout<<"abs(idhep) \t= "<<abs(i_hep->idhep())<<std::endl;
		// 	std::cout<<"daFirst() \t= "<<i_hep->daFirst()<<std::endl;
		// 	std::cout<<"daLast() \t= "<<i_hep->daLast()<<std::endl;
		// }
	}
	if(filled==1){return vect_pi;}

	else{return {};}
}

std::vector<HepLorentzVector> quadri_momentum_pi(Gen_hepevt_Manager & HevMgr, Particle B0)
{
	std::vector<HepLorentzVector> quadri_momentum(8);
	std::vector<Particle> vect_pi = get_pions_FS(HevMgr);
	if (vect_pi.empty())
	{
		HepLorentzVector fake;
		fake.setX(-999);fake.setPx(-999);
		fake.setY(-999);fake.setPy(-999);
		fake.setZ(-999);fake.setPz(-999);
		fake.setT(-999);
		return {fake,fake,fake,fake,fake,fake,fake,fake};
	}

	Particle HevMgr_pi1_rho = vect_pi[0];
	Particle HevMgr_pi2_rho = vect_pi[1];
	Particle HevMgr_pi1_ks0 = vect_pi[2];
	Particle HevMgr_pi2_ks0 = vect_pi[3];
	Particle Cand_pi1_rho = B0.child(0).child(1);
	Particle Cand_pi2_rho = B0.child(0).child(2);
	Particle Cand_pi1_ks0 = B0.child(0).child(0).child(0);
	Particle Cand_pi2_ks0 = B0.child(0).child(0).child(1);

	if(HevMgr_pi1_rho.charge()==1)
	{
		quadri_momentum[0] = HevMgr_pi1_rho.p();
		quadri_momentum[1] = HevMgr_pi2_rho.p();
	}
	else
	{
		quadri_momentum[1] = HevMgr_pi1_rho.p();
		quadri_momentum[0] = HevMgr_pi2_rho.p();
	}

	if(HevMgr_pi1_ks0.charge()==1)
	{
		quadri_momentum[2] = HevMgr_pi1_ks0.p();
		quadri_momentum[3] = HevMgr_pi2_ks0.p();
	}
	else
	{
		quadri_momentum[3] = HevMgr_pi1_ks0.p();
		quadri_momentum[2] = HevMgr_pi2_ks0.p();
	}

	if(Cand_pi1_rho.charge()==1)
	{
		quadri_momentum[4] = Cand_pi1_rho.p();
		quadri_momentum[5] = Cand_pi2_rho.p();
	}
	else
	{
		quadri_momentum[5] = Cand_pi1_rho.p();
		quadri_momentum[4] = Cand_pi2_rho.p();
	}

	if(Cand_pi1_ks0.charge()==1)
	{
		quadri_momentum[6] = Cand_pi1_ks0.p();
		quadri_momentum[7] = Cand_pi2_ks0.p();
	}
	else
	{
		quadri_momentum[7] = Cand_pi1_ks0.p();
		quadri_momentum[6] = Cand_pi2_ks0.p();
	}
	
	return quadri_momentum;
}

std::vector<double> get_angles_pions_FS(Gen_hepevt_Manager & HevMgr, Particle B0)
{
	std::vector<double> angles_btw_rec_mc(4);
	std::vector<double> angles_btw_pi(4);
	std::vector<Particle> vect_pi = get_pions_FS(HevMgr);
	if (vect_pi.empty()){return {-999,-999,-999,-999};}

	Particle HevMgr_pi1_rho = vect_pi[0];
	Particle HevMgr_pi2_rho = vect_pi[1];
	Particle HevMgr_pi1_ks0 = vect_pi[2];
	Particle HevMgr_pi2_ks0 = vect_pi[3];
	Particle Cand_pi1_rho = B0.child(0).child(1);
	Particle Cand_pi2_rho = B0.child(0).child(2);
	Particle Cand_pi1_ks0 = B0.child(0).child(0).child(0);
	Particle Cand_pi2_ks0 = B0.child(0).child(0).child(1);

	if(HevMgr_pi1_rho.charge()==Cand_pi1_rho.charge())
	{
		if(HevMgr_pi1_rho.charge()==1)
		{
			angles_btw_rec_mc[0] = (HevMgr_pi1_rho.p3()).angle(Cand_pi1_rho.p3());
			angles_btw_rec_mc[1] = (HevMgr_pi2_rho.p3()).angle(Cand_pi2_rho.p3());
		}
		else
		{
			angles_btw_rec_mc[1] = (HevMgr_pi1_rho.p3()).angle(Cand_pi1_rho.p3());
			angles_btw_rec_mc[0] = (HevMgr_pi2_rho.p3()).angle(Cand_pi2_rho.p3());
		}
	}
	else
	{
		if(HevMgr_pi1_rho.charge()==1)
		{
			angles_btw_rec_mc[0] = (HevMgr_pi1_rho.p3()).angle(Cand_pi2_rho.p3());
			angles_btw_rec_mc[1] = (HevMgr_pi2_rho.p3()).angle(Cand_pi1_rho.p3());
		}
		else
		{
			angles_btw_rec_mc[1] = (HevMgr_pi1_rho.p3()).angle(Cand_pi2_rho.p3());
			angles_btw_rec_mc[0] = (HevMgr_pi2_rho.p3()).angle(Cand_pi1_rho.p3());
		}
	}

	if(HevMgr_pi1_ks0.charge()==Cand_pi1_ks0.charge())
	{
		if(HevMgr_pi1_ks0.charge()==1)
		{
			angles_btw_rec_mc[2] = (HevMgr_pi1_ks0.p3()).angle(Cand_pi1_ks0.p3());
			angles_btw_rec_mc[3] = (HevMgr_pi2_ks0.p3()).angle(Cand_pi2_ks0.p3());
		}
		else
		{
			angles_btw_rec_mc[3] = (HevMgr_pi1_ks0.p3()).angle(Cand_pi1_ks0.p3());
			angles_btw_rec_mc[2] = (HevMgr_pi2_ks0.p3()).angle(Cand_pi2_ks0.p3());
		}
	}
	else
	{
		if(HevMgr_pi1_ks0.charge()==1)
		{
			angles_btw_rec_mc[2] = (HevMgr_pi1_ks0.p3()).angle(Cand_pi2_ks0.p3());
			angles_btw_rec_mc[3] = (HevMgr_pi2_ks0.p3()).angle(Cand_pi1_ks0.p3());
		}
		else
		{
			angles_btw_rec_mc[3] = (HevMgr_pi1_ks0.p3()).angle(Cand_pi2_ks0.p3());
			angles_btw_rec_mc[2] = (HevMgr_pi2_ks0.p3()).angle(Cand_pi1_ks0.p3());
		}
	}
	return angles_btw_rec_mc;
}

int is_signal_angle(Gen_hepevt_Manager & HevMgr, Particle B0)
{
	std::vector<Particle> vect_pi = get_pions_FS(HevMgr);
	if (vect_pi.empty()){ return 0;}

	std::vector<double> angles = get_angles_pions_FS(HevMgr,B0);
	std::sort(angles.begin(), angles.end(), greater());

	if (angles[0]<0.5)
	{ 
		return 1;
	}
	else {return 0;}
}

int get_mcpdg(Particle B0)
{
	// Sometimes the relation is not set
	if(B0.relation().genHepevt()==1)
	{
		return B0.relation().genHepevt().idhep();
	}
	else
	{
		return -999;
	}
}

double get_mcP(Particle B0)
{
	Gen_hepevt B0gen = B0.relation().genHepevt();
	// Sometimes the relation is not set
	if(B0gen==1)
	{
		if(abs(B0gen.idhep())==511)
		{
			return sqrt(B0gen.P(0)*B0gen.P(0) + B0gen.P(1)*B0gen.P(1) + B0gen.P(2)*B0gen.P(2));
		}
	}
	else
	{
		return -999;
	}
}

std::vector<double> get_pgen(Gen_hepevt_Manager & HevMgr)
{
	std::vector<double> pgen = {-999,-999};
	for(std::vector<Gen_hepevt>::const_iterator i_hep = HevMgr.begin(); i_hep != HevMgr.end(); i_hep++) 
	{
		// Y(4S)
		if(abs(i_hep->idhep())==300553 && (HevMgr[i_hep->get_ID()-1].daLast()-HevMgr[i_hep->get_ID()-1].daFirst()+1)==2)
		{
			// Loop on the 2 B0
			for(int i=HevMgr[i_hep->get_ID()-1].daFirst(); i<=HevMgr[i_hep->get_ID()-1].daLast(); i++)
			{
				// std::cout<<HevMgr[i-1].idhep()<<std::endl;
				if((HevMgr[i-1].idhep())==511)
				{
					pgen[0] = sqrt(HevMgr[i-1].P(1)*HevMgr[i-1].P(1) + HevMgr[i-1].P(2)*HevMgr[i-1].P(2) + HevMgr[i-1].P(0)*HevMgr[i-1].P(0));
					// std::cout<<pgen[0]<<" "<<HevMgr[i-1].P(1)<<" "<<HevMgr[i-1].P(2)<<" "<<HevMgr[i-1].P(3)<<" "<<std::endl;
				}
				if((HevMgr[i-1].idhep())==-511)
				{
					pgen[1] = sqrt(HevMgr[i-1].P(1)*HevMgr[i-1].P(1) + HevMgr[i-1].P(2)*HevMgr[i-1].P(2) + HevMgr[i-1].P(0)*HevMgr[i-1].P(0));
					// std::cout<<pgen[1]<<" "<<HevMgr[i-1].P(1)<<" "<<HevMgr[i-1].P(2)<<" "<<HevMgr[i-1].P(3)<<" "<<std::endl;
				}
			}
		}
	}
	// std::cout<<pgen[0]<<" "<<pgen[1]<<std::endl;
	return pgen;
}

void varghese_troubleshoot(string prefix)
{
Gen_hepevt_Manager &gen_mgr          = Gen_hepevt_Manager::get_manager();

for (std::vector<Gen_hepevt>::iterator j = gen_mgr.begin(); j != gen_mgr.end(); ++j){//loop over mc particles

if (abs((*j).idhep()) == 511){//B or Bbar
Gen_hepevt &child1 =  *(gen_mgr.begin() - 1 + (*j).daFirst());
Gen_hepevt &child2 =  *(gen_mgr.begin() - 1 + (*j).daLast());
if(child1.idhep() == 10313 && child2.idhep() == 22)
{cout<<prefix<<" flavour1"<<endl;}
if(child1.idhep() == -10313 && child2.idhep() == 22)
{cout<<prefix<<" flavour2"<<endl;}
}//B or Bbar
}//loop over mc particles

}//varghese_troubleshoot


//recursive function
int ancestorBind(const Gen_hepevt& mcparticle){

if(!mcparticle){return -999;}//termination due to no more mother

if( abs(mcparticle.idhep()) == 511 || abs(mcparticle.idhep()) == 521 ){//expected termination at B
return mcparticle.get_ID();
}
else{
int ancBind = ancestorBind(mcparticle.mother());
return ancBind;
}

}//recursive function



//Getting MC information to find out what physics processes contribute to BB-backgrounds, 
//see if these have any CP asymmetry of their own, in order to assign a syst-uncertainty 
//on S, to account for this effect
void true_mcBancestor(Particle B0, 
int &b1pdg, int &b1ind, int (&b1daupdg)[6], 
int &b2pdg, int &b2ind, int (&b2daupdg)[6], 
int &ksbind, int &pi1bind, int &pi2bind, int &gambind)
{
Gen_hepevt_Manager &GenMgr = Gen_hepevt_Manager::get_manager();

//PART 1 :  b1pdg, b1ind, b1daupdg, b2pdg, b2ind, b2daupdg
int Bcount = 0;
for (std::vector<Gen_hepevt>::iterator j = GenMgr.begin(); j != GenMgr.end(); ++j){//loop over mc particles
if (abs((*j).idhep()) == 511 || abs((*j).idhep()) == 521){//B0, B0bar, B+ or B-
Bcount++;

if(Bcount == 1){//Bcount == 1
int daucount = -1;

b1pdg = (*j).idhep();
b1ind = (*j).get_ID();
//cout<<"B1-pdg = "<<b1pdg<<", B1-index = "<<b1ind<<endl;
for (int i = (*j).daFirst(); i <= (*j).daLast(); ++i) {//loop over B daughters
daucount++;
if(i==0 || daucount > 5)continue;
const Gen_hepevt& child = GenMgr(Panther_ID(i));
b1daupdg[daucount] = child.idhep();
//cout<<"B1 child"<<daucount<<"-pdg = "<<b1daupdg[daucount]<<"\t";
}//loop over B daughters
//cout<<endl;
}//Bcount == 1

else if(Bcount == 2){//Bcount == 2
int daucount = -1;
b2pdg = (*j).idhep();
b2ind = (*j).get_ID();
//cout<<"B2-pdg = "<<b2pdg<<", B2-index = "<<b2ind<<endl;
for (int i = (*j).daFirst(); i <= (*j).daLast(); ++i) {//loop over B daughters
daucount++;
if(i==0 || daucount > 5)continue;
const Gen_hepevt& child = GenMgr(Panther_ID(i));
b2daupdg[daucount] = child.idhep();
//cout<<"B2 child"<<daucount<<"-pdg = "<<b2daupdg[daucount]<<"\t";
}//loop over B daughters
//cout<<endl;
}//Bcount == 2

}//B0, B0bar, B+ or B-

}//loop over mc particles

//PART 2 : ksbind, pi1bind, pi2bind, gambind
ksbind = ancestorBind(B0.child(0).child(0).genHepevt());
//cout<<"ks-B-index = "<<ksbind<<endl;
pi1bind = ancestorBind(B0.child(0).child(1).genHepevt());
//cout<<"pi1-B-index = "<<pi1bind<<endl;
pi2bind = ancestorBind(B0.child(0).child(2).genHepevt());
//cout<<"pi2-B-index = "<<pi2bind<<endl;
gambind = ancestorBind(B0.child(1).genHepevt());
//cout<<"gamma-B-index = "<<gambind<<endl;
//cout<<endl;
//cout<<"------------------------------------------------------------------"<<endl;
//cout<<endl;
}//true_mcBancestor





//============================================
// BASF module class definition.
//=============================================
class RecoDecays : public Module {
   // public members whi'ch are fixed procedures.
   public:
      RecoDecays ( void );          // constructor
      ~RecoDecays ( void ){};       // destructor
      void init ( int* );        // to be compatible with 20000430
      void term( void ){};
      void disp_stat ( const char* ) {} ;// to be compatible with 20000430
      void hist_def ( void );
      void event ( BelleEvent*, int* );
      void begin_run ( BelleEvent*, int* );// to be compatible with 20000430
      void end_run ( BelleEvent*, int* ){};// to be compatible with 20000430
      void other ( int*, BelleEvent*, int* ){};
      void getWandDeltaW(int no_exp, double qr, double &w, double &dw);
   // histgrams and ntuples have to be
   // declared as private members.
   private:
  //BelleHistogram *H_masspipi, *H_masskskpi;
  BelleTuple *T_bcand;
  BelleTuple *T_kscand;
  BelleTuple *T_k10cand;
	// float m_limits[8];
	// float m_svd1ws[7];
	// float m_svd2ws[7];
	// float m_svd1dws[7];
	// float m_svd2dws[7];
	float m_limits[8] = { 0, 0.1, 0.25, 0.5, 0.625, 0.75, 0.875, 1 };
	float m_svd1ws[7] = { 0.418852, 0.329879, 0.233898, 0.170608, 0.099791, 0.0228501 };
	float m_svd2ws[7] = { 0.418826, 0.319303, 0.222948, 0.163191, 0.104085, 0.0251454 };
	float m_svd1dws[7] = { 0.0569661, 0.0126192, -0.0147724, -0.000550289, 0.00887704, 0.00465683 };
	float m_svd2dws[7] = { -0.00877001, 0.0103515, -0.0109253, -0.0186365, 0.00168037, -0.0036441 };

};
  
//========================
// BASF Interface Function: 
//======================== 
extern "C" Module_descr *mdcl_RecoDecays ()
{ // need both module decleration and module description.
   RecoDecays *module = new RecoDecays; 
   Module_descr *dscr = new Module_descr ( "RecoDecays", module );
   IpProfile::define_global(dscr);
   return dscr;
}

//====================
// Member Functions
//====================
//

void RecoDecays::getWandDeltaW(int no_exp, double qr, double &w, double &dw)
{
	const int nbins = 7;
	for (int startbin = 0; startbin < nbins; startbin++)
	{
		int endbin = startbin +1;
		if (std::fabs(qr) > m_limits[startbin] && std::fabs(qr) < m_limits[endbin])
		{
			if(no_exp < 31) //CASE A
			{
				w  = m_svd1ws[startbin];
				dw = m_svd1dws[startbin];
			}
			else
			{
				w  = m_svd2ws[startbin];
				dw = m_svd2dws[startbin];
			}
		}
	}	
}

//------------
// Constructor
//------------
RecoDecays::RecoDecays ( void ) 
{
}

//--------------
// init function
//--------------
void RecoDecays::init ( int* ) 
{
	Hamlet::init();
	// https://belle.kek.jp/group/indirectcp/cpfit/2010/wtag_nn_2010 
	// DATA
	// m_limits     = { 0, 0.1, 0.25, 0.5, 0.625, 0.75, 0.875, 1 };
	// m_svd1ws     = { 0.5, 0.429307, 0.337238, 0.246732, 0.189617, 0.108998, 0.0257188 };
	// m_svd2ws     = { 0.5, 0.433074, 0.334596, 0.23933, 0.185627, 0.111058, 0.0233918 };
	// m_svd1dws    = { 0.0, -0.0137777, 0.0160512, -0.0165945, 0.0148771, 0.0123041, -0.000992505 };
	// m_svd2dws    = { 0.0, -0.0164068, -0.00390874, -0.0278437, 0.00826071, -0.00357243, -0.00544135 };
	// Hardcoded official w and dw values taken from a summer school tutorials
	//m_svd1ws     = { 0.418852, 0.329879, 0.233898, 0.170608, 0.099791, 0.0228501 };
	//m_svd2ws     = { 0.418826, 0.319303, 0.222948, 0.163191, 0.104085, 0.0251454 };
	//m_svd1dws    = { 0.0569661, 0.0126192, -0.0147724, -0.000550289, 0.00887704, 0.00465683 };
	//m_svd2dws    = { -0.00877001, 0.0103515, -0.0109253, -0.0186365, 0.00168037, -0.0036441 };
}

//------------------
// hist_def function
//------------------
void RecoDecays::hist_def ( void )
{
   extern BelleTupleManager *BASF_Histogram; // TupleManager is a wrapper
   BelleTupleManager& tm = *BASF_Histogram;  // class to use hbook and ntuple.
   //T_kscand = tm.ntuple("kscand", "m p mcflags pcms costheta isSignal");
   //T_k10cand = tm.ntuple("k10cand", "m p mcflags pcms costheta isSignal m12 m23 m13");
   T_bcand = tm.ntuple("b", "nexp nrun eventid ecm ebeam "
   				"m mbc de p mcflags ecms pcms costheta coscms mcpdg b0pgen b0bpgen mcP "
				"isSignal issgnevt sigmc sigrec sigk1ks sigk1k0s sigk1k2p sigk1rho sigk2ks sigk2rho sigk1kpp sigk1omg sigk1f0 issgnang issgngam ncands ncand "
				"a_ro_pi1 a_ro_pi2 a_ks_pi1 a_ks_pi2 a_pi_max "
				"p1mc_px p1mc_py p1mc_pz p1mc_E "
				"p2mc_px p2mc_py p2mc_pz p2mc_E "
				"p3mc_px p3mc_py p3mc_pz p3mc_E "
				"p4mc_px p4mc_py p4mc_pz p4mc_E "
				"p1_px p1_py p1_pz p1_E "
				"p2_px p2_py p2_pz p2_E "
				"p3_px p3_py p3_pz p3_E "
				"p4_px p4_py p4_pz p4_E "
				"qr w dw wtag flavor rbin genbq "
				"zrec vtxcl vtxh vtxchi2 vtxndf vtxntrk vtxzerr vtxyerr vtxxerr "
				"ztag tagcl tagh tagchi2 tagchi2wo tagndf tagndfwo taglepton "
				"tagzerr tagyerr tagxerr tagntrk "
				"m12 m23 m13 "
				"dt dtgen dtres dtpull dterr trecgen ttaggen "
				"dz dzgen dzres dzpull dzerr zrecgen ztaggen "
   				"x_p x_pcms x_issignal x_m x_ks_p x_ks_m x_ks_pcms x_ks_cos g_p g_pcms g_issignal g_cos x_cos g_e9e25 " 
				"x_p0_p x_p1_p x_p0_cos x_p1_cos x_p0_kid x_p1_kid x_p0_eid x_p1_eid x_p0_mid x_p1_mid x_p0_cdc x_p1_cdc x_p0_svd x_p1_svd "
   				"k0mm2 k0et k0hso00 k0hso01 k0hso02 k0hso03 k0hso04 k0hso10 k0hso12 k0hso14 k0hso20 k0hso22 k0hso24 k0hoo0 k0hoo1 k0hoo2 k0hoo3 k0hoo4 costh "
				"pi0veto etaveto "
				"mcb1_pdg mcb1_ind mcb1_da1 mcb1_da2 mcb1_da3 mcb1_da4 mcb1_da5 mcb1_da6 mcb2_pdg mcb2_ind mcb2_da1 mcb2_da2 mcb2_da3 mcb2_da4 mcb2_da5 mcb2_da6 ksmcbind p1mcbind p2mcbind gmcbind");
				//"ntracks0 ntracks1 ntracks2");
}

//================
// Event function
//================
void RecoDecays::event ( BelleEvent* evptr, int* status )
{


	int eventID = std::rand();

	// Get Mdst tables.
	Mdst_charged_Manager& mdstchgmgr = Mdst_charged_Manager::get_manager();
	Mdst_trk_Manager& mdsttrkmgr = Mdst_trk_Manager::get_manager();
	Mdst_trk_fit_Manager& mdsttrkfitmgr = Mdst_trk_fit_Manager::get_manager();
	// For photons.
	Mdst_ecl_Manager& mdsteclmgr = Mdst_ecl_Manager::get_manager();
	Mdst_pi0_Manager& mdstpi0mgr = Mdst_pi0_Manager::get_manager(); 
	Mdst_gamma_Manager& mdstgammgr = Mdst_gamma_Manager::get_manager();
	Mdst_ecl_aux_Manager& mdstauxmgr = Mdst_ecl_aux_Manager::get_manager();
	Belle_runhead_Manager &runhead_m = Belle_runhead_Manager::get_manager();
	Gen_hepevt_Manager &gen_mgr          = Gen_hepevt_Manager::get_manager();
	Belle_runhead &runhead = runhead_m( ( Panther_ID) 1 );
	bool runOnData = (runhead.ExpMC())? 0:1;
	if (runOnData && !IpProfile::usable())
	{
		return;
	}
	int no_exp = 0;
	int no_run = 0;
	int no_evt = 0;
	get_event_id(no_exp, no_run, no_evt);


	// Vectors to store candidate tracks.
	std::vector<Particle> vect_pi_plus; 
	std::vector<Particle> vect_pi_minus ;
	std::vector<Particle> vect_K_short; 
	std::vector<Particle> vect_gamma; 

	// Particle type definition.
	Ptype ptype_pi_plus("PI+");
	Ptype ptype_pi_minus("PI-");
	Ptype ptype_gamma("GAMM");
	Ptype ptype_b0(511); // B^0
	Ptype ptype_b0b(-511); // \bar{B}^0
	// Set CM frame as HepLorentsVector(constant) to boost back.
	const double anglex = BeamEnergy::Cross_angle();//0.022;  // Finite angle crossing, in rad. 
	const double Eher = BeamEnergy::E_HER();//7.998213; // HER beam energy.(New 2003)
	const double Eler = BeamEnergy::E_LER();//3.499218; // LER beam energy.(New 2003)
	// LorentzVector of lab. respective to Upsiron's rest frame.

	double ecm = BeamEnergy::Ecm();
	double Ebeam = Benergy(); 

	//----
	// Set atc_pid to reject kaon.
	//----
	atc_pid selKpi(3,1,5,3,2);

	// Same for electrons and muons
	atc_pid selEpi(3,1,5,0,2);
	atc_pid selMupi(3,1,5,1,2);

	
	//======
	// Sort up charged tracks to K+, pi+, K- and pi-.
	//======
	for(std::vector<Mdst_charged>::iterator itc = mdstchgmgr.begin();
	itc!= mdstchgmgr.end(); itc++){
	Mdst_charged& track = *itc;

	// K or pi?
	if( track.charge()>0.0 ){ // pi+
	  Particle particle_pi_plus( track, ptype_pi_plus );
	  HepVector a_pi( param_at_ip(particle_pi_plus));
	  // dr and dz cut.
	  if( fabs( (double)a_pi[0] )<1.0 && fabs( (double)a_pi[3] )< 3.0 ){
	    vect_pi_plus.push_back( particle_pi_plus );
	  }
	}else{ // pi-
	  Particle particle_pi_minus( track, ptype_pi_minus );
	  HepVector a_pi( param_at_ip(particle_pi_minus));
	  // dr and dz cut.
	  if( fabs( (double)a_pi[0] )<1.0 && fabs( (double)a_pi[3] )< 3.0 ){
	    vect_pi_minus.push_back( particle_pi_minus );
	  }
	}
	}// Mdst_charged loop end.

	// Mdst_vee2 if kind=1, then register as KS -> pi+pi- candidates.
	Mdst_vee2_Manager& mdstvee2mgr = Mdst_vee2_Manager::get_manager();   
	Mdst_vee_daughters_Manager& mdstveedaumgr 
	= Mdst_vee_daughters_Manager::get_manager();   
	// Mass window to accept KS -> pi+ pi-.
	const double mpipi_ul = 0.4976 + 0.020; // PDG KS mass + 20 MeV.
	const double mpipi_ll = 0.4976 - 0.020; // PDG KS mass - 20 MeV.
	// nisKSFinder needs IP.
	const Hep3Vector& IP = IpProfile::position();
	//if (!IP.mag()){ IP.setX(0.0); IP.setY(0.0); IP.setZ(0.0);}

	for(std::vector<Mdst_vee2>::iterator it_vee2=mdstvee2mgr.begin(); it_vee2!=mdstvee2mgr.end(); it_vee2++){
		Mdst_vee2& vee2 = *it_vee2;
		if( vee2.kind()==1 ){
			// nisKsFinder selection.
			nisKsFinder ksnb;
			ksnb.candidates(vee2, IP);
			if( ksnb.standard()==1){
				Particle particle_ks( vee2 );
				// KS candidate mass cut within +- 20 MeV, as defined earlier.
				if( particle_ks.p().mag() > mpipi_ll && particle_ks.p().mag() < mpipi_ul ){
					vect_K_short.push_back( particle_ks );
				}// Final KS mass window
			}// nisKsFinder 
		}// if vee2.kind()==1 end.
	}// Mdst_vee2 loop end.
	//=====
	// Needed instructions for gamma.
	//=====
	// Set gamma position and its error.
	const HepPoint3D origin; // Default constructor sets (0,0,0).
	static double largeError = 1.0; // 1.0*1.0cm^2.
	HepSymMatrix dx( 3, 0);
	dx[0][0] = largeError;
	dx[1][1] = largeError;
	dx[2][2] = largeError;
	// Convert Mdst_gamma into Particle object.
	for(std::vector<Mdst_gamma>::iterator it_gamma=mdstgammgr.begin();
		it_gamma!=mdstgammgr.end(); ++it_gamma)
	{
		Mdst_gamma& gamma = *it_gamma;
		Particle particle_gam( gamma );
		// Calculate error matrix and set it.
		HepSymMatrix errCart(4,0);
		errGam(errCart, gamma);
		//// std::cout<<"errCart "<<errCart[0][0]<<std::endl;
		particle_gam.momentum().momentum( particle_gam.p(), errCart );

		// Set largeError of position for kmassfitter call later.
		particle_gam.momentum().position( origin, dx );
		vect_gamma.push_back( particle_gam );
	}
	
	std::vector<Particle> Xsd_particle_list;
	std::vector<Particle> B0_particle_list_raw; B0_particle_list_raw.clear(); 
	std::vector<Particle> B0_particle_list; B0_particle_list.clear();
	// run pi0 eta veto and e9e25 for every photon:
	fillPhotonProperties(vect_gamma, Eher, Eler, anglex);
	// Reconstruct Kaon resonance Xsd:
	combination(Xsd_particle_list, Ptype("K10"), vect_K_short, vect_pi_plus, vect_pi_minus);
	// Apply a loose cut on Xsd mass:
	withMassCut(Xsd_particle_list,0.9, 2.5);
	// Apply a hard cut on Xsd CMS energy
	withPSCut(Xsd_particle_list,1.,3.5);
	// Apply a hard cut on photon CMS energy
	withPSCut(vect_gamma,1.5,3.5);
	// Reconstruct B0:
	combination(B0_particle_list_raw, ptype_b0, Xsd_particle_list, vect_gamma);
	
	for(int i = 0; i < B0_particle_list_raw.size(); i++)
	{

		HepLorentzVector b04vector = B0_particle_list_raw[i].p();
		b04vector.boost(-BeamEnergy::CMBoost());
	
		float mbc = calc_mbc(B0_particle_list_raw[i]);
		float de = b04vector.t() - Ebeam;
		// Cuts on mbc and de:
		if (mbc > 5.2 && de > -0.2 && de < 0.2)
		{
			/************************************************************/
			/***********************FLAVOR TAG***************************/
			/************************************************************/
			Hamlet hamlet;
			hamlet.setBcp( B0_particle_list_raw[i] ); // set B_CP side tracks
			hamlet.setTagMethod(Hamlet::MULT_DIM_LH);
			//hamlet.setTagMethod(Hamlet::NN1);
			//Fbtag_NN1 tagger_NN = hamlet.fbtg_NN1();
  			Fbtag_MultDimLikelihood0 tagger_MDL = hamlet.fbtg_mult_dim_likelihood();

			//Fbtag_NN1 tagger_NN = hamlet_NN.fbtg_NN1();
			double fq_nn = hamlet.q(); // will return flavor *q
			double wtag = hamlet.wtag(); // 
			int flavor = hamlet.flavor();
			double w = -1;
			double dw = -2;
			getWandDeltaW(no_exp, fq_nn, w,dw);


			/************************************************************/
			/***********************VERTEX FIT***************************/
			/************************************************************/
			kvertexfitter kvf1;
			Particle pi_plus = B0_particle_list_raw[i].child(0).child(1);
			if(enough_svd_hit(pi_plus))
			{
				addTrack2fit(kvf1, pi_plus);
			}
			Particle pi_minus = B0_particle_list_raw[i].child(0).child(2);
			if(enough_svd_hit(pi_minus))
			{
				addTrack2fit(kvf1, pi_minus);
			}
			addTube2fit(kvf1);
			unsigned err1 = kvf1.fit();
			HepPoint3D signalVertex = kvf1.vertex();
			HepSymMatrix signalVertexError = kvf1.errVertex();
			double cl = kvf1.cl();
			double chi2 = kvf1.chisq();
			double ndf = kvf1.dgf();
			if( !IpProfile::usable() )
			{
				//cl = -1;
			};
			/************************************************************/
			/************************TAG V FIT***************************/
			/************************************************************/
			TagVK tagv;
			tagv.setdefault(B0_particle_list_raw[i], signalVertex);
			std::vector<Particle> taglist;
			std::vector<Particle> k_p, k_m, pi_p, pi_m;
			makeKPi(k_p, k_m, pi_p, pi_m, 0);
			taglist = pi_p;
			taglist.insert(taglist.end(), pi_m.begin(), pi_m.end());
			for(int j=0; j<B0_particle_list_raw[i].relation().nFinalStateParticles(); j++)
			{
				removeParticle(taglist,B0_particle_list_raw[i].relation().finalStateParticle(j));
			}
			for(std::vector<Particle>::iterator it=taglist.begin(); it!=taglist.end(); ++it)
			{
				Particle * particle = &*it;
				tagv.push_back(particle);
			}
			tagv.fit();
			HepPoint3D tagVertex = tagv.vtx();
			HepSymMatrix tagVertexError = tagv.errVtx();
			double tagCL = tagv.cl();
			double tagChi2 = tagv.chisq();
			double tagNdf = tagv.ndf();
			int tagNtrk = tagv.ntrk();
			int taglepton = tagv.isTagLeptonVertex() ? (int)tagv.VertexTagLepton().get_ID() : 0;
			/************************************************************/
			/********************CONTINUUM SUPRESSION********************/
			/************************************************************/
			ksfwmoments km(B0_particle_list_raw[i], Ebeam, - BeamEnergy::CMBoost());
			km.usefinal(0);
			double cosThrustAngle = std::fabs(std::cos(getThrustAngle(B0_particle_list_raw[i], Eher, Eler, anglex)));
			BUserInfo bUI(cl, chi2, ndf, tagCL, tagChi2, tagNdf, tagNtrk, taglepton, cosThrustAngle, fq_nn, w, dw, wtag, flavor);
			//bUI.setHamlet(NULL);
			bUI.setKSFW(km);
			bUI.setB0Vtx(signalVertex, signalVertexError);
			bUI.setTagVtx(tagVertex, tagVertexError);
			bUI.setB0tagwo(tagv.ndf_woip(), tagv.chisq_woip()); 
			B0_particle_list_raw[i].userInfo(bUI);
			BUserInfo & mybui =  dynamic_cast<BUserInfo&>(B0_particle_list_raw[i].userInfo());
			//mybui.setHamlet(hamlet_NN2);
			//
			B0_particle_list.push_back(B0_particle_list_raw[i]);
		}
	}



	/**************************************************/
	/**************OUTPUT IN ROOTFILE******************/
	/**************************************************/
	std::sort(B0_particle_list.begin(), B0_particle_list.end(), sortBFunction);
	if (!runOnData)
	{
		setMCtruth(vect_K_short);
		setMCtruth(Xsd_particle_list);
		setMCtruth(B0_particle_list);
	}
	// if(B0_particle_list.size() > 0) 
	for(int b0_i=0;b0_i<B0_particle_list.size();b0_i++) 
	{
                Particle b0 = B0_particle_list[b0_i];

	    //Looking at Mbc sideband only for data before unblinding
	    //if(calc_mbc(b0) < 5.26){
		T_bcand->column("nexp", no_exp );
		T_bcand->column("nrun", no_run );
		T_bcand->column("eventid", eventID );
		std::cout<<"INFO Nexp "<<no_exp<<" Nrun "<<no_run<<" EventID "<<eventID<<std::endl;

		T_bcand->column("ecm", ecm/2. );
		T_bcand->column("ebeam", Ebeam );
		T_bcand->column("p", b0.p().vect().mag() );
		T_bcand->column("m12", (b0.child(0).child(0).p()+b0.child(0).child(1).p()).mag() );
		T_bcand->column("m23", (b0.child(0).child(1).p()+b0.child(0).child(2).p()).mag() );
		T_bcand->column("m13", (b0.child(0).child(0).p()+b0.child(0).child(2).p()).mag() );
		T_bcand->column("g_p", b0.child(1).p().vect().mag() );
		T_bcand->column("x_p", b0.child(0).p().vect().mag() );
		T_bcand->column("x_m", b0.child(0).mass() );
		T_bcand->column("x_ks_m", b0.child(0).child(0).mass() );
		T_bcand->column("x_ks_p", b0.child(0).child(0).p().vect().mag() );
		T_bcand->column("x_ks_cos", b0.child(0).child(0).p().vect().cosTheta() );
		HepLorentzVector ks4vector = b0.child(0).child(0).p();
		ks4vector.boost(-BeamEnergy::CMBoost());
		T_bcand->column("x_ks_pcms", ks4vector.vect().mag() );
		T_bcand->column("x_issignal", (getMCtruthFlag(b0.child(0)) == 1)? 1:0 );
		T_bcand->column("g_issignal", (getMCtruthFlag(b0.child(1)) == 1)? 1:0 );
		T_bcand->column("g_cos", b0.child(1).p().vect().cosTheta());
		T_bcand->column("x_cos", b0.child(0).p().vect().cosTheta());
		T_bcand->column("x_p0_p", b0.child(0).child(1).p().vect().mag() );
		T_bcand->column("x_p1_p", b0.child(0).child(2).p().vect().mag() );
		T_bcand->column("x_p0_cos", b0.child(0).child(1).p().vect().cosTheta() );
		T_bcand->column("x_p1_cos", b0.child(0).child(2).p().vect().cosTheta() );
		T_bcand->column("x_p0_kid", selKpi.prob( b0.child(0).child(1).mdstCharged()) );
		T_bcand->column("x_p1_kid", selKpi.prob( b0.child(0).child(2).mdstCharged()) );
		T_bcand->column("x_p0_mid", selMupi.prob( b0.child(0).child(1).mdstCharged()) );
		T_bcand->column("x_p1_mid", selMupi.prob( b0.child(0).child(2).mdstCharged()) );
		T_bcand->column("x_p0_eid", selEpi.prob( b0.child(0).child(1).mdstCharged()) );
		T_bcand->column("x_p1_eid", selEpi.prob( b0.child(0).child(2).mdstCharged()) );
		T_bcand->column("x_p0_cdc", b0.child(0).child(1).mdstCharged().trk().mhyp(2).nhits(0)+b0.child(0).child(1).mdstCharged().trk().mhyp(2).nhits(1) );
		T_bcand->column("x_p0_svd", b0.child(0).child(1).mdstCharged().trk().mhyp(2).nhits(3)+b0.child(0).child(1).mdstCharged().trk().mhyp(2).nhits(4) );
		T_bcand->column("x_p1_cdc", b0.child(0).child(2).mdstCharged().trk().mhyp(2).nhits(0)+b0.child(0).child(2).mdstCharged().trk().mhyp(2).nhits(1) );
		T_bcand->column("x_p1_svd", b0.child(0).child(2).mdstCharged().trk().mhyp(2).nhits(3)+b0.child(0).child(2).mdstCharged().trk().mhyp(2).nhits(4) );
		//const atc_pid & xsd_pi0_PID =  b0.child(0).child(1).pId().kId();
		//const atc_pid & xsd_pi1_PID =  b0.child(0).child(2).pId().kId();
		HepLorentzVector xsd4vector = b0.child(0).p();
		xsd4vector.boost(-BeamEnergy::CMBoost());
		T_bcand->column("x_pcms", xsd4vector.vect().mag());
		HepLorentzVector gamma4vector = b0.child(1).p();
		gamma4vector.boost(-BeamEnergy::CMBoost());
		T_bcand->column("g_pcms", gamma4vector.vect().mag());
		HepLorentzVector b04vector = b0.p();
		b04vector.boost(-BeamEnergy::CMBoost());
		T_bcand->column("pcms", b04vector.vect().mag() );
		T_bcand->column("de", b04vector.t() - Ebeam );
		T_bcand->column("m", b0.mass() );
		float mbc = calc_mbc(b0);
		float m_B = 5.27963;
		T_bcand->column("mbc", mbc);
		T_bcand->column("ecms", std::sqrt(m_B*m_B + b04vector.rho()*b04vector.rho()));
		T_bcand->column("costheta", b0.p().vect().cosTheta() );
		T_bcand->column("coscms", b04vector.vect().cosTheta() );
		T_bcand->column("mcpdg", get_mcpdg(b0));
		T_bcand->column("mcflags", getMCtruthFlag(b0));
		T_bcand->column("isSignal", (getMCtruthFlag(b0) == 1)? 1:0);
		// sigk1ks sigk1k0s sigk1k2p sigk1rho sigk2ks sigk2rho sigk1kpp sigk1omg sigk1f0
		int sigmc = is_signal_mcgen(gen_mgr);
		int sigrec = is_signal_finalstate(gen_mgr);
		int sigk1ks = is_signal_k1ks(gen_mgr);
		int sigk1k0s = is_signal_k1k0s(gen_mgr);
		int sigk1k2p = is_signal_k1k2p(gen_mgr);
		int sigk1rho = is_signal_k1rho(gen_mgr);
		int sigk2ks = is_signal_k2ks(gen_mgr);
		int sigk2rho = is_signal_k2rho(gen_mgr);
		int sigk1kpp = is_signal_k1kpp(gen_mgr);
		int sigk1omg = is_signal_k1omg(gen_mgr);
		int sigk1f0 = is_signal_k1f0(gen_mgr);
		int issgnevt = 0;
		if (sigmc || sigk1ks || sigk1k0s || sigk1k2p || sigk1rho || sigk2ks || sigk2rho || sigk1kpp || sigk1omg || sigk1f0){issgnevt = 1;}
		//if(getMCtruthFlag(b0)==1 && issgnevt==0){coutevtgen(gen_mgr);}
		T_bcand->column("issgnevt", (issgnevt == 1)? 1:0);
		std::cout<<"INFO Issgnevt is "<<((issgnevt == 1)? 1:0)<<std::endl;		
		T_bcand->column("sigmc", (sigmc >= 1)? 1:0);
		T_bcand->column("sigrec", (sigrec >= 1)? 1:0);
		std::cout<<"INFO Total event is a yes "<<((sigrec >= 1)? 1:0)<<std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		T_bcand->column("sigk1ks", (sigk1ks >= 1)? 1:0);
		T_bcand->column("sigk1k0s", (sigk1k0s >= 1)? 1:0);
		T_bcand->column("sigk1k2p", (sigk1k2p >= 1)? 1:0);
		T_bcand->column("sigk1rho", (sigk1rho >= 1)? 1:0);
		T_bcand->column("sigk2ks", (sigk2ks >= 1)? 1:0);
		T_bcand->column("sigk2rho", (sigk2rho >= 1)? 1:0);
		T_bcand->column("sigk1kpp", (sigk1kpp >= 1)? 1:0);
		T_bcand->column("sigk1omg", (sigk1omg >= 1)? 1:0);
		T_bcand->column("sigk1f0", (sigk1f0 >= 1)? 1:0);
		T_bcand->column("issgnang", (is_signal_angle(gen_mgr,b0) == 1)? 1:0);
		T_bcand->column("issgngam", (is_signal_no_gamma_res(gen_mgr) == 1)? 1:0);
		std::vector<double> pgen = get_pgen(gen_mgr);
		T_bcand->column("b0pgen", pgen[0]);
		T_bcand->column("b0bpgen", pgen[1]);
		T_bcand->column("mcP", get_mcP(b0));
		if(b0_i==0){T_bcand->column("ncands", B0_particle_list.size());}
		T_bcand->column("ncand", b0_i);
		std::vector<double> angles = get_angles_pions_FS(gen_mgr,b0);
		T_bcand->column("a_ro_pi1", angles[0]);
		T_bcand->column("a_ro_pi2", angles[1]);
		T_bcand->column("a_ks_pi1", angles[2]);
		T_bcand->column("a_ks_pi2", angles[3]);
		std::sort(angles.begin(), angles.end(), greater());
		T_bcand->column("a_pi_max", angles[0]);
		std::vector<HepLorentzVector> quadri_momentum = quadri_momentum_pi(gen_mgr,b0);
		T_bcand->column("p1mc_px", quadri_momentum[0].px());
		T_bcand->column("p1mc_py", quadri_momentum[0].py());
		T_bcand->column("p1mc_pz", quadri_momentum[0].pz());
		T_bcand->column("p1mc_E", quadri_momentum[0].e());
		T_bcand->column("p2mc_px", quadri_momentum[1].px());
		T_bcand->column("p2mc_py", quadri_momentum[1].py());
		T_bcand->column("p2mc_pz", quadri_momentum[1].pz());
		T_bcand->column("p2mc_E", quadri_momentum[1].e());
		T_bcand->column("p3mc_px", quadri_momentum[2].px());
		T_bcand->column("p3mc_py", quadri_momentum[2].py());
		T_bcand->column("p3mc_pz", quadri_momentum[2].pz());
		T_bcand->column("p3mc_E", quadri_momentum[2].e());
		T_bcand->column("p4mc_px", quadri_momentum[3].px());
		T_bcand->column("p4mc_py", quadri_momentum[3].py());
		T_bcand->column("p4mc_pz", quadri_momentum[3].pz());
		T_bcand->column("p4mc_E", quadri_momentum[3].e());
		T_bcand->column("p1_px", quadri_momentum[4].px());
		T_bcand->column("p1_py", quadri_momentum[4].py());
		T_bcand->column("p1_pz", quadri_momentum[4].pz());
		T_bcand->column("p1_E", quadri_momentum[4].e());
		T_bcand->column("p2_px", quadri_momentum[5].px());
		T_bcand->column("p2_py", quadri_momentum[5].py());
		T_bcand->column("p2_pz", quadri_momentum[5].pz());
		T_bcand->column("p2_E", quadri_momentum[5].e());
		T_bcand->column("p3_px", quadri_momentum[6].px());
		T_bcand->column("p3_py", quadri_momentum[6].py());
		T_bcand->column("p3_pz", quadri_momentum[6].pz());
		T_bcand->column("p3_E", quadri_momentum[6].e());
		T_bcand->column("p4_px", quadri_momentum[7].px());
		T_bcand->column("p4_py", quadri_momentum[7].py());
		T_bcand->column("p4_pz", quadri_momentum[7].pz());
		T_bcand->column("p4_E", quadri_momentum[7].e());
		BUserInfo& binfo = dynamic_cast<BUserInfo&>(b0.userInfo());
		T_bcand->column("qr", binfo.getB0qr());
		T_bcand->column("w", binfo.getW());
		T_bcand->column("dw", binfo.getdW());
		T_bcand->column("wtag", binfo.getWtag());
		T_bcand->column("flavor", binfo.getFlavor());
		T_bcand->column("rbin", 1.-(2*binfo.getWtag()));

		T_bcand->column("vtxcl", binfo.getB0CL());
		T_bcand->column("vtxh", binfo.getB0Chi2()/binfo.getB0Ndf());
		T_bcand->column("vtxchi2", binfo.getB0Chi2());
		T_bcand->column("vtxndf", binfo.getB0Ndf());
		T_bcand->column("vtxntrk", binfo.getB0Ndf()+1);
		T_bcand->column("tagcl", binfo.getTagCL());
		T_bcand->column("tagh", binfo.getTagChi2()/binfo.getTagNdf());
		T_bcand->column("tagchi2", binfo.getTagChi2());
		T_bcand->column("tagndf", binfo.getTagNdf());
		T_bcand->column("tagchi2wo", binfo.getTagChi2wo());
		T_bcand->column("tagndfwo", binfo.getTagNdfwo());
		T_bcand->column("tagntrk", binfo.getTagNtrk());
		T_bcand->column("taglepton", binfo.getTagLepton());
		double dz = -11;
		double dzerr = -1;
		double vtxzerr = -1;
		double tagzerr = -1;
		double vtxxerr = -1;
		double tagxerr = -1;
		double vtxyerr = -1;
		double tagyerr = -1;
		double vtxz = -1;
		double tagz = -1;
		double dt = -11;
		double dterr = -1;
		int genBFlavor = 0;
		std::vector<double> dtdz = getDtDz(gen_mgr, genBFlavor, runOnData);
		double dtgen = dtdz[0];
		double trecgen = dtdz[1];
		double ttaggen = dtdz[2];
		double dzgen = dtdz[3];
		double zrecgen = dtdz[4];
		double ztaggen = dtdz[5];
		if (binfo.getB0CL() > 0.0 && binfo.getTagCL() > 0.0)
		{
			dz = binfo.getB0Vtx()(2) - binfo.getTagVtx()(2);
			tagz = binfo.getTagVtx()(2);
			vtxz = binfo.getB0Vtx()(2);
			dzerr = std::sqrt(binfo.getB0VtxErr()(2,2) + binfo.getTagVtxErr()(2,2));
			dt = dz*10.0/0.12741;
			dterr = dzerr*10.0/0.12741;
			tagzerr = std::sqrt(binfo.getTagVtxErr()(3,3));
			vtxzerr = std::sqrt(binfo.getB0VtxErr()(3,3));
			tagyerr = std::sqrt(binfo.getTagVtxErr()(2,2));
			vtxyerr = std::sqrt(binfo.getB0VtxErr()(2,2));
			tagxerr = std::sqrt(binfo.getTagVtxErr()(1,1));
			vtxxerr = std::sqrt(binfo.getB0VtxErr()(1,1));
		}
		T_bcand->column("genbq", genBFlavor);
		T_bcand->column("tagzerr", tagzerr);
		T_bcand->column("vtxzerr", vtxzerr);
		T_bcand->column("tagyerr", tagyerr);
		T_bcand->column("vtxyerr", vtxyerr);
		T_bcand->column("tagxerr", tagxerr);
		T_bcand->column("vtxxerr", vtxxerr);
		T_bcand->column("ztag", tagz);
		T_bcand->column("zrec", vtxz);
		T_bcand->column("dt", dt);
		T_bcand->column("dtgen", dtgen);
		T_bcand->column("dtres", dt - dtgen);
		T_bcand->column("dtpull", (dt - dtgen)/dterr);
		T_bcand->column("dterr", dterr);
		T_bcand->column("trecgen", trecgen);
		T_bcand->column("ttaggen", ttaggen);
		T_bcand->column("dz", dz);
		T_bcand->column("dzgen", dzgen);
		T_bcand->column("dzres", dz - dzgen);
		T_bcand->column("dzpull", (dz - dzgen)/dzerr);
		T_bcand->column("dzerr", dzerr);
		T_bcand->column("zrecgen", zrecgen);
		T_bcand->column("ztaggen", ztaggen);
		ksfwmoments km = binfo.getKSFW();
		T_bcand->column("k0mm2",  km.mm2());
		T_bcand->column("k0et",   km.et());
		T_bcand->column("k0hso00",km.Hso(0,0));
		T_bcand->column("k0hso01",km.Hso(0,1));
		T_bcand->column("k0hso02",km.Hso(0,2));
		T_bcand->column("k0hso03",km.Hso(0,3));
		T_bcand->column("k0hso04",km.Hso(0,4));
		T_bcand->column("k0hso10",km.Hso(1,0));
		T_bcand->column("k0hso12",km.Hso(1,2));
		T_bcand->column("k0hso14",km.Hso(1,4));
		T_bcand->column("k0hso20",km.Hso(2,0));
		T_bcand->column("k0hso22",km.Hso(2,2));
		T_bcand->column("k0hso24",km.Hso(2,4));
		T_bcand->column("k0hoo0", km.Hoo(0));
		T_bcand->column("k0hoo1", km.Hoo(1));
		T_bcand->column("k0hoo2", km.Hoo(2));
		T_bcand->column("k0hoo3", km.Hoo(3));
		T_bcand->column("k0hoo4", km.Hoo(4));
		T_bcand->column("costh", binfo.getCosT());
		GammaUserInfo& gammainfo = dynamic_cast<GammaUserInfo&>(b0.child(1).userInfo());
		T_bcand->column("pi0veto", gammainfo.getPi0Veto());
		T_bcand->column("etaveto", gammainfo.getEtaVeto());
		T_bcand->column("g_e9e25", gammainfo.getE9E25());


		int b1pdg = -999, b1ind = -999, b1daupdg[6] = {-999, -999, -999, -999, -999, -999},
		b2pdg = -999, b2ind = -999, b2daupdg[6] = {-999, -999, -999, -999, -999, -999},
		ksbind = -999, pi1bind = -999, pi2bind = -999, gambind = -999;
		true_mcBancestor(b0, b1pdg, b1ind, b1daupdg, b2pdg, b2ind, b2daupdg, ksbind, pi1bind, pi2bind, gambind);

		T_bcand->column("mcb1_pdg", b1pdg);
		T_bcand->column("mcb1_ind", b1ind);
		T_bcand->column("mcb1_da1", b1daupdg[0]);
		T_bcand->column("mcb1_da2", b1daupdg[1]);
                T_bcand->column("mcb1_da3", b1daupdg[2]);          
                T_bcand->column("mcb1_da4", b1daupdg[3]);
                T_bcand->column("mcb1_da5", b1daupdg[4]);          
                T_bcand->column("mcb1_da6", b1daupdg[5]);		

		T_bcand->column("mcb2_pdg", b2pdg);
                T_bcand->column("mcb2_ind", b2ind);
                T_bcand->column("mcb2_da1", b2daupdg[0]);          
                T_bcand->column("mcb2_da2", b2daupdg[1]);
                T_bcand->column("mcb2_da3", b2daupdg[2]);
                T_bcand->column("mcb2_da4", b2daupdg[3]);
                T_bcand->column("mcb2_da5", b2daupdg[4]);
                T_bcand->column("mcb2_da6", b2daupdg[5]);

		T_bcand->column("ksmcbind", ksbind);
                T_bcand->column("p1mcbind", pi1bind);
                T_bcand->column("p2mcbind", pi2bind);
                T_bcand->column("gmcbind", gambind);

		T_bcand->dumpData();
	    //}//data sideband cut : if(calc_mbc(b0) < 5.26)
	    //cout<<"****************** Signal Candidate END ******************"<<endl;
	}
	//cout<<"++++++++++++++++++++++++++++++ EVENT END ++++++++++++++++++++++++++++++"<<endl;
}// envent function end.

//================
// Begin_run function
//================
void RecoDecays::begin_run( BelleEvent*, int* )
{
  // IP Profile initialization is necessary to use.
	BeamEnergy::begin_run();
	IpProfile::begin_run();
	Hamlet::begin_run(Hamlet::MULT_DIM_LH ); // MDLH method
	Hamlet::begin_run(Hamlet::NN1); // Neural network
   // If you want to use eid, do its initialization here.
   // eid::init_data();
}

#if defined(BELLE_NAMESPACE)
}
#endif
