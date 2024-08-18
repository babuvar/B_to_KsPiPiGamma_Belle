#include "RooTatamiHelper.h"
#include "tatami/tatami.h"
#include <string>
#include <regex.h>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>

RooDataSet ascii2DataSet(const std::string& asciifilename,
		bool mc,
		RooRealVar& dt,
		RooRealVar& mbc,
		RooRealVar& deltaE,
		RooCategory& expno,
		RooRealVar& CosThetaB,
		RooRealVar& ntrk_sig,
		RooRealVar& ntrk_tag,
		RooRealVar& z_err_sig,
		RooRealVar& z_err_tag,
		RooRealVar& chisq_tracks_sig,
		RooRealVar& chisq_tracks_tag,
		RooRealVar& dgf_tracks_sig,
		RooRealVar& dgf_tracks_tag,
		RooRealVar& keeptagl_tag,
		RooRealVar& ecms,
		RooCategory& rbin,
		RooRealVar& wtag,
		RooRealVar& delta_wtag,
		RooCategory& flavor_tag,
		RooRealVar& ebeam) {

	FILE *fp = fopen(asciifilename.c_str(), "r");
	if (!fp) {
		std::cout << "Failed to open " << asciifilename << std::endl;
		perror("fopen");
		exit(1);
	} else {
		std::cout << "Creating new RooDataSet data set with " << asciifilename << " as input." << std::endl;
	}

	RooArgSet m_varset(dt, mbc, deltaE, expno, CosThetaB, ntrk_sig, ntrk_tag, z_err_sig, z_err_tag);
	m_varset.add(chisq_tracks_sig);
	m_varset.add(chisq_tracks_tag);
	m_varset.add(dgf_tracks_sig);
	m_varset.add(dgf_tracks_tag);
	m_varset.add(keeptagl_tag);
	m_varset.add(ecms);
	m_varset.add(rbin);
	m_varset.add(wtag);
	m_varset.add(delta_wtag);
	m_varset.add(flavor_tag);
	m_varset.add(ebeam);
	RooDataSet r_set("r_set", "r_set", m_varset);

	int m_icpv_exp, m_icpv_run, m_icpv_evt, m_icpv_frm, m_icpv_hepevt, m_icpv_mc_flv;
	double m_cos_theta_tr, m_cos_theta_1, m_icpv_mc_delta_t, m_icpv_mc_pbs, m_icpv_mc_costh;
	int m_icpv_decay, m_icpv_flv, m_icpv_rec_vt_ntrk, m_icpv_asc_vt_ntrk, m_icpv_rec_vt_ndf, m_icpv_asc_vt_ndf, m_icpv_keeptagl;
	double m_icpv_rec_vt_pos, m_icpv_rec_vt_err, m_icpv_rec_vt_chi2, m_icpv_rec_vt_cl;
	double m_icpv_asc_vt_pos, m_icpv_asc_vt_err, m_icpv_asc_vt_chi2, m_icpv_asc_vt_cl;
	double m_icpv_wtag_mc, m_icpv_deltae, m_icpv_mbc, m_icpv_pbs, m_icpv_costh, m_icpv_ebeam, m_icpv_costhtr;
	double m_icpv_dt, m_icpv_r, m_icpv_wtag, m_icpv_delta_wtag, m_icpv_ecms;
	int m_icpv_rbin;
	float m_dummy_f;

	char buffer[1024];

	while (fgets(buffer, sizeof (buffer), fp))
	{
		int ntoken = 0;

		if (!(pattern_match(buffer)) == true) continue;

		ntoken = sscanf(
				buffer,
				"%d %d %d %d "
				"%d "
				"%d "
				"%le %le %le "
				"%le %le "
				"%d "
				"%le %le %d %le %d %le "
				"%le %le %d %le %d %le "
				"%f "
				"%d  %d %le "
				"%le %le %le %le %le %le \n",

				&m_icpv_exp, &m_icpv_run, &m_icpv_evt, &m_icpv_frm,
				&m_icpv_hepevt,
				&m_icpv_mc_flv,
				&m_cos_theta_tr, &m_cos_theta_1,
				&m_icpv_mc_delta_t,
				&m_icpv_mc_pbs, &m_icpv_mc_costh,
				&m_icpv_decay,
				&m_icpv_rec_vt_pos, &m_icpv_rec_vt_err, &m_icpv_rec_vt_ntrk, &m_icpv_rec_vt_chi2, &m_icpv_rec_vt_ndf, &m_icpv_rec_vt_cl,
				&m_icpv_asc_vt_pos, &m_icpv_asc_vt_err, &m_icpv_asc_vt_ntrk, &m_icpv_asc_vt_chi2, &m_icpv_asc_vt_ndf, &m_icpv_asc_vt_cl,
				&m_dummy_f,
				&m_icpv_flv, &m_icpv_keeptagl, &m_icpv_wtag_mc,
				&m_icpv_deltae, &m_icpv_mbc, &m_icpv_pbs, &m_icpv_costh,
				&m_icpv_ebeam, &m_icpv_costhtr
		);

		if (ntoken != 34) continue; // have read ascii record for event, now calculate needed quantities (e.g. determine r bin, wtag, sigfrac, ...)

		//calculate variables
		m_icpv_dt = (m_icpv_rec_vt_pos - m_icpv_asc_vt_pos) * Belle::dt_resol_global::inv_bgc;

		m_icpv_r = 1 - 2 * m_icpv_wtag_mc;
		m_icpv_rbin = get_rbin(m_icpv_r);
		m_icpv_wtag = get_wtag(m_icpv_exp, m_icpv_rbin, mc);
		m_icpv_delta_wtag = get_delta_wtag(m_icpv_exp, m_icpv_rbin, mc);

		m_icpv_ecms = sqrt(m_icpv_pbs * m_icpv_pbs + Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero);

		//perform ICPV cuts
		if (m_icpv_rbin == 7) continue;
		if (m_icpv_rec_vt_err <= 0) continue;
		if (m_icpv_asc_vt_err <= 0) continue;

		if (m_icpv_rec_vt_ndf > 0 && m_icpv_rec_vt_chi2 / m_icpv_rec_vt_ndf >= 50) continue;
		if (m_icpv_asc_vt_ndf > 0 && m_icpv_asc_vt_chi2 / m_icpv_asc_vt_ndf >= 50) continue;

		if (m_icpv_rec_vt_ntrk == 1 && m_icpv_rec_vt_err >= 0.05) continue;
		if (m_icpv_asc_vt_ntrk == 1 && m_icpv_asc_vt_err >= 0.05) continue;

		if (m_icpv_rec_vt_ntrk > 1 && m_icpv_rec_vt_err >= 0.02) continue;
		if (m_icpv_asc_vt_ntrk > 1 && m_icpv_asc_vt_err >= 0.02) continue;

		typedef std::pair<double, RooRealVar*> valuevar;
		std::vector<valuevar> vars;
		vars.push_back(make_pair(m_icpv_dt, &dt));
		vars.push_back(make_pair(m_icpv_deltae, &deltaE));
		vars.push_back(make_pair(m_icpv_mbc, &mbc));
		//vars.push_back(make_pair(m_icpv_exp, &expno));
		vars.push_back(make_pair(m_icpv_costh, &CosThetaB));
		vars.push_back(make_pair(m_icpv_rec_vt_ntrk, &ntrk_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_ntrk, &ntrk_tag));
		vars.push_back(make_pair(m_icpv_rec_vt_err, &z_err_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_err, &z_err_tag));
		vars.push_back(make_pair(m_icpv_rec_vt_chi2, &chisq_tracks_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_chi2, &chisq_tracks_tag));
		vars.push_back(make_pair(m_icpv_rec_vt_ndf, &dgf_tracks_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_ndf, &dgf_tracks_tag));
		vars.push_back(make_pair(m_icpv_keeptagl, &keeptagl_tag));
		vars.push_back(make_pair(m_icpv_ecms, &ecms));
		vars.push_back(make_pair(m_icpv_wtag, &wtag));
		vars.push_back(make_pair(m_icpv_delta_wtag, &delta_wtag));
		vars.push_back(make_pair(m_icpv_ebeam, &ebeam));

		bool in_range = true;


		BOOST_FOREACH(valuevar &p, vars) {
			if (p.first < p.second->getMin() || p.first > p.second->getMax()) {
				in_range = false;
				break;
			} else
				p.second->setVal(p.first);
		}

		if (flavor_tag.setIndex(m_icpv_flv, false)) in_range = false;
		if (rbin.setIndex(m_icpv_rbin, false)) in_range = false;

		if (m_icpv_exp > 30)
			expno.setIndex(31, false);
		else
			expno.setIndex(7, false);

		//if (ntrk_sig.setIndex(m_icpv_rec_vt_ntrk, false)) in_range = false;
		//if (ntrk_tag.setIndex(m_icpv_asc_vt_ntrk, false)) in_range = false;
		//if (dgf_tracks_sig.setIndex(m_icpv_rec_vt_ndf, false)) in_range = false;
		//if (dgf_tracks_tag.setIndex(m_icpv_asc_vt_ndf, false)) in_range = false;
		//if (keeptagl_tag.setIndex(m_icpv_keeptagl, false)) in_range = false;

		if (in_range) r_set.add(m_varset);

	}
	cout << "number of events for fit: " << r_set.numEntries() << endl;
	fclose(fp);
	return r_set;
}

RooDataSet ascii2DataSet(const std::string& asciifilename,
		bool mc,
		RooRealVar& dt,
		RooRealVar& mbc,
		RooRealVar& deltaE,
		RooCategory& expno,
		RooRealVar& CosThetaB,
		RooRealVar& ntrk_sig,
		RooRealVar& ntrk_tag,
		RooRealVar& z_err_sig,
		RooRealVar& z_err_tag,
		RooRealVar& chisq_tracks_sig,
		RooRealVar& chisq_tracks_tag,
		RooRealVar& dgf_tracks_sig,
		RooRealVar& dgf_tracks_tag,
		RooRealVar& keeptagl_tag,
		RooRealVar& ecms,
		RooCategory& rbin,
		RooRealVar& wtag,
		RooRealVar& delta_wtag,
		RooCategory& flavor_tag,
		RooRealVar& ebeam,
		std::vector<read_in_variable> add_var) {

	FILE *fp = fopen(asciifilename.c_str(), "r");
	if (!fp) {
		std::cout << "Failed to open " << asciifilename << std::endl;
		perror("fopen");
		exit(1);
	} else {
		std::cout << "Creating new RooDataSet data set with " << asciifilename << " as input." << std::endl;
	}

	RooArgSet m_varset(dt, mbc, deltaE, expno, CosThetaB, ntrk_sig, ntrk_tag, z_err_sig, z_err_tag);
	m_varset.add(chisq_tracks_sig);
	m_varset.add(chisq_tracks_tag);
	m_varset.add(dgf_tracks_sig);
	m_varset.add(dgf_tracks_tag);
	m_varset.add(keeptagl_tag);
	m_varset.add(ecms);
	m_varset.add(rbin);
	m_varset.add(wtag);
	m_varset.add(delta_wtag);
	m_varset.add(flavor_tag);
	m_varset.add(ebeam);

	BOOST_FOREACH(read_in_variable& ri_var, add_var) {
		m_varset.add(*(ri_var.first));
	}

	RooDataSet r_set("r_set", "r_set", m_varset);

	double d_readin[34];
	int i_readin[34];
	//int m_icpv_exp, m_icpv_run, m_icpv_evt, m_icpv_frm, m_icpv_hepevt, m_icpv_mc_flv;
	//double m_cos_theta_tr, m_cos_theta_1, m_icpv_mc_delta_t, m_icpv_mc_pbs, m_icpv_mc_costh;
	//int m_icpv_decay, m_icpv_flv, m_icpv_rec_vt_ntrk, m_icpv_asc_vt_ntrk, m_icpv_rec_vt_ndf, m_icpv_asc_vt_ndf, m_icpv_keeptagl;
	//double m_icpv_rec_vt_pos, m_icpv_rec_vt_err, m_icpv_rec_vt_chi2, m_icpv_rec_vt_cl;
	//double m_icpv_asc_vt_pos, m_icpv_asc_vt_err, m_icpv_asc_vt_chi2, m_icpv_asc_vt_cl;
	//double m_icpv_wtag_mc, m_icpv_deltae, m_icpv_mbc, m_icpv_pbs, m_icpv_costh, m_icpv_ebeam, m_icpv_costhtr;
	double m_icpv_dt, m_icpv_r, m_icpv_wtag, m_icpv_delta_wtag, m_icpv_ecms;
	int m_icpv_rbin;
	//float m_dummy_f;

	char buffer[1024];

	while (fgets(buffer, sizeof (buffer), fp))
	{
		int ntoken = 0;

		if (!(pattern_match(buffer)) == true) continue;

		string scan_string = "%d %d %d %d %d %d %le %le %le %le %le %d %le %le %d %le %d %le %le %le %d %le %d %le %f %d %d %le %le %le %le %le %le %le \n";

		ntoken = sscanf(
				buffer,
				scan_string.c_str(),

				&i_readin[0], &d_readin[1], &d_readin[2], &d_readin[3], &d_readin[4], &d_readin[5], &d_readin[6], &d_readin[7], &d_readin[8], &d_readin[9],
				&d_readin[10], &d_readin[11], &d_readin[12], &d_readin[13], &i_readin[14], &d_readin[15], &i_readin[16], &d_readin[17], &d_readin[18], &d_readin[19],
				&i_readin[20], &d_readin[21], &i_readin[22], &d_readin[23], &d_readin[24], &i_readin[25], &i_readin[26], &d_readin[27], &d_readin[28], &d_readin[29],
				&d_readin[30], &d_readin[31], &d_readin[32], &d_readin[33]
				                                                       //&m_icpv_exp, &m_icpv_run, &m_icpv_evt, &m_icpv_frm,
				                                                       //&m_icpv_hepevt,
				                                                       //&m_icpv_mc_flv,
				                                                       //&m_cos_theta_tr, &m_cos_theta_1,
				                                                       //&m_icpv_mc_delta_t,
				                                                       //&m_icpv_mc_pbs, &m_icpv_mc_costh,
				                                                       //&m_icpv_decay,
				                                                       //&m_icpv_rec_vt_pos, &m_icpv_rec_vt_err, &m_icpv_rec_vt_ntrk, &m_icpv_rec_vt_chi2, &m_icpv_rec_vt_ndf, &m_icpv_rec_vt_cl,
				                                                       //&m_icpv_asc_vt_pos, &m_icpv_asc_vt_err, &m_icpv_asc_vt_ntrk, &m_icpv_asc_vt_chi2, &m_icpv_asc_vt_ndf, &m_icpv_asc_vt_cl,
				                                                       //&m_dummy_f,
				                                                       //&m_icpv_flv, &m_icpv_keeptagl, &m_icpv_wtag_mc,
				                                                       //&m_icpv_deltae, &m_icpv_mbc, &m_icpv_pbs, &m_icpv_costh,
				                                                       //&m_icpv_ebeam, &m_icpv_costhtr
		);

		if (ntoken != 34) continue; // have read ascii record for event, now calculate needed quantities (e.g. determine r bin, wtag, sigfrac, ...)

		//calculate variables
		double m_icpv_rec_vt_pos = d_readin[12];
		double m_icpv_asc_vt_pos = d_readin[18];
		double m_icpv_wtag_mc = d_readin[27];
		int m_icpv_exp = i_readin[0];
		double m_icpv_pbs = d_readin[30];
		double m_icpv_rec_vt_err = d_readin[13];
		double m_icpv_asc_vt_err = d_readin[19];
		int m_icpv_rec_vt_ndf = i_readin[16];
		int m_icpv_asc_vt_ndf = i_readin[22];
		double m_icpv_rec_vt_chi2 = d_readin[15];
		double m_icpv_asc_vt_chi2 = d_readin[21];
		int m_icpv_rec_vt_ntrk = i_readin[14];
		int m_icpv_asc_vt_ntrk = i_readin[20];
		double m_icpv_deltae = d_readin[28];
		double m_icpv_mbc = d_readin[29];
		double m_icpv_costh = d_readin[31];
		int m_icpv_keeptagl = i_readin[26];
		double m_icpv_ebeam = d_readin[32];
		int m_icpv_flv = i_readin[25];

		m_icpv_dt = (m_icpv_rec_vt_pos - m_icpv_asc_vt_pos) * Belle::dt_resol_global::inv_bgc;

		m_icpv_r = 1 - 2 * m_icpv_wtag_mc;
		m_icpv_rbin = get_rbin(m_icpv_r);
		m_icpv_wtag = get_wtag(m_icpv_exp, m_icpv_rbin, mc);
		m_icpv_delta_wtag = get_delta_wtag(m_icpv_exp, m_icpv_rbin, mc);

		m_icpv_ecms = sqrt(m_icpv_pbs * m_icpv_pbs + Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero);

		//perform ICPV cuts
		if (m_icpv_rbin == 7) continue;
		if (m_icpv_rec_vt_err <= 0) continue;
		if (m_icpv_asc_vt_err <= 0) continue;

		if (m_icpv_rec_vt_ndf > 0 && m_icpv_rec_vt_chi2 / m_icpv_rec_vt_ndf >= 50) continue;
		if (m_icpv_asc_vt_ndf > 0 && m_icpv_asc_vt_chi2 / m_icpv_asc_vt_ndf >= 50) continue;

		if (m_icpv_rec_vt_ntrk == 1 && m_icpv_rec_vt_err >= 0.05) continue;
		if (m_icpv_asc_vt_ntrk == 1 && m_icpv_asc_vt_err >= 0.05) continue;

		if (m_icpv_rec_vt_ntrk > 1 && m_icpv_rec_vt_err >= 0.02) continue;
		if (m_icpv_asc_vt_ntrk > 1 && m_icpv_asc_vt_err >= 0.02) continue;

		typedef std::pair<double, RooRealVar*> valuevar;
		std::vector<valuevar> vars;
		vars.push_back(make_pair(m_icpv_dt, &dt));
		vars.push_back(make_pair(m_icpv_deltae, &deltaE));
		vars.push_back(make_pair(m_icpv_mbc, &mbc));
		//vars.push_back(make_pair(m_icpv_exp, &expno));
		vars.push_back(make_pair(m_icpv_costh, &CosThetaB));
		vars.push_back(make_pair(m_icpv_rec_vt_ntrk, &ntrk_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_ntrk, &ntrk_tag));
		vars.push_back(make_pair(m_icpv_rec_vt_err, &z_err_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_err, &z_err_tag));
		vars.push_back(make_pair(m_icpv_rec_vt_chi2, &chisq_tracks_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_chi2, &chisq_tracks_tag));
		vars.push_back(make_pair(m_icpv_rec_vt_ndf, &dgf_tracks_sig));
		vars.push_back(make_pair(m_icpv_asc_vt_ndf, &dgf_tracks_tag));
		vars.push_back(make_pair(m_icpv_keeptagl, &keeptagl_tag));
		vars.push_back(make_pair(m_icpv_ecms, &ecms));
		vars.push_back(make_pair(m_icpv_wtag, &wtag));
		vars.push_back(make_pair(m_icpv_delta_wtag, &delta_wtag));
		vars.push_back(make_pair(m_icpv_ebeam, &ebeam));

		BOOST_FOREACH(read_in_variable& ri_var, add_var) {
			vars.push_back(make_pair(d_readin[ri_var.second - 1], ri_var.first));
		}

		bool in_range = true;

		BOOST_FOREACH(valuevar &p, vars) {
			if (p.first < p.second->getMin() || p.first > p.second->getMax()) {
				in_range = false;
				break;
			} else
				p.second->setVal(p.first);
		}

		if (flavor_tag.setIndex(m_icpv_flv, false)) in_range = false;
		if (rbin.setIndex(m_icpv_rbin, false)) in_range = false;

		if (m_icpv_exp > 30)
			expno.setIndex(31, false);
		else
			expno.setIndex(1, false);

		//if (ntrk_sig.setIndex(m_icpv_rec_vt_ntrk, false)) in_range = false;
		//if (ntrk_tag.setIndex(m_icpv_asc_vt_ntrk, false)) in_range = false;
		//if (dgf_tracks_sig.setIndex(m_icpv_rec_vt_ndf, false)) in_range = false;
		//if (dgf_tracks_tag.setIndex(m_icpv_asc_vt_ndf, false)) in_range = false;
		//if (keeptagl_tag.setIndex(m_icpv_keeptagl, false)) in_range = false;

		if (in_range) r_set.add(m_varset);

	}
	cout << "number of events for fit: " << r_set.numEntries() << endl;
	fclose(fp);
	return r_set;
}

bool pattern_match(const std::string& temp_buf) {
	static regex_t preg;
	const char *buf = temp_buf.c_str();

	static bool first = true;
	if (first) {
		regcomp(&preg, "^[-+.0-9 \tEe]*\n$", 0);
		first = false;
	}

	if (regexec(&preg, buf, 0, NULL, 0) == REG_NOMATCH) {
		fprintf(stderr, "skipped bad event record (regexec): %s", buf);
		return false;
	}

	return true;
}

inline int get_rbin(double r) {
	return (0. <= r && r <= 0.1 ? 0 :
	0.1 < r && r <= 0.25 ? 1 :
	0.25 < r && r <= 0.5 ? 2 :
	0.5 < r && r <= 0.625 ? 3 :
	0.625 < r && r <= 0.75 ? 4 :
	0.75 < r && r <= 0.875 ? 5 :
	0.875 < r && r <= 1.0 ? 6 : 7);
}

inline double get_wtag(int expno, int rbin, bool mc) {
	double w_svd1_data[7] = {0.5, 0.418852, 0.329879, 0.233898, 0.170608, 0.099791, 0.0228501};

	double w_svd2_data[7] = {0.5, 0.418826, 0.319303, 0.222948, 0.163191, 0.104085, 0.0251454};

	double w_svd1_mc[7] = {0.5, 0.420827, 0.300296, 0.219317, 0.154636, 0.0916131, 0.0228891};

	double w_svd2_mc[7] = {0.5, 0.412222, 0.307838, 0.212765, 0.149933, 0.0913264, 0.0218754};

	if (mc) {
		if (expno < 30) {
			return w_svd1_mc[rbin];
		} else 	{
			return w_svd2_mc[rbin];
		}
	} else {
		if (expno < 30) {
			return w_svd1_data[rbin];
		} else {
			return w_svd2_data[rbin];
		}
	}
}

inline double get_delta_wtag(int expno, int rbin, bool mc) {
	double dw_svd1_data[7] = {0., 0.0569661, 0.0126192, -0.0147724, -0.000550289, 0.00887704, 0.00465683};

	double dw_svd2_data[7] = {0., -0.00877001, 0.0103515, -0.0109253, -0.0186365, 0.00168037, -0.0036441};

	double dw_svd1_mc[7] = {0., 0.0583019, 0.00573998, -0.0392635, 0.00474508, -0.0118737, -0.00585326};

	double dw_svd2_mc[7] = {0., 0.00408778, 0.010326, -0.00479522, 0.00151989, 0.0143633, 0.00189979};

	if (mc) {
		if (expno < 30) {
			return dw_svd1_mc[rbin];
		} else {
			return dw_svd2_mc[rbin];
		}
	} else {
		if (expno < 30) {
			return dw_svd1_data[rbin];
		} else {
			return dw_svd2_data[rbin];
		}
	}
}

//Only for testing purposes
void print_dtres_params(int expno, int mc){

Belle::dtres_param_t* dtres_param = Belle::get_dtres_param(expno, mc);

cout<<"-------------------------------"<<endl;
cout<<"dtres_param->Srec[0] = "<<dtres_param->Srec[0]<<endl;
cout<<"dtres_param->Srec[1] = "<<dtres_param->Srec[1]<<endl;
cout<<"dtres_param->Sasc[0] = "<<dtres_param->Sasc[0]<<endl;
cout<<"dtres_param->Sasc[1] = "<<dtres_param->Sasc[1]<<endl;
cout<<"-------------------------------"<<endl;

}


void dtres_systematics_b0(int item, double sigma, int expno, int mc){

Belle::dtres_systematics(Belle::dtres_items_fullrecon_b0[item], sigma, expno, mc);


}

double get_dtres_parVal_b0(int item, int subitem, int expno, int mc, string& parname){

Belle::dtres_param_t* param = Belle::get_dtres_param(expno, mc);


switch(item){
  case 0: //par-0
  parname = "Rdet_Srec0";
    return (param->Srec)[0];
  case 1: //par-1
  parname = "Rdet_Srec1";
    return (param->Srec)[1];
  case 2: //par-2
  parname = "Rdet_Sasc0";
    return (param->Sasc)[0];
  case 3: //par-3
  parname = "Rdet_Sasc1";
    return (param->Sasc)[1];
  case 4: //par-4
  parname = "Rnp_Snp";
    return param->Snp;
  case 5: //par-5
  parname = "Rnp_Snp_global";
    return param->Snp_global;
  case 6: //par-6
    if(subitem == 0){
  parname = "Rdet_Smn_onetrk_rec";
    return param->Smn_rec;}
    if(subitem == 1) {
  parname = "Rdet_Smn_onetrk_asc";
    return param->Smn_asc;}
  case 7: //par-7
    if(subitem == 0) {
  parname = "Rdet_Stl_onetrk_rec";
    return param->Stl_rec;}
    if(subitem == 1) {
  parname = "Rdet_Stl_onetrk_asc";
    return param->Stl_asc;}
  case 8: //par-8
    if(subitem == 0) {
  parname = "Rdet_ftl_onetrk_rec";
    return param->ftl_rec;}
    if(subitem == 1) {
  parname = "Rdet_ftl_onetrk_asc";
    return param->ftl_asc;}
  case 9: //par-9
    if(subitem == 0) {parname = "Rdet_rec_mlt_ftail0";}
    if(subitem == 1) {parname = "Rdet_rec_mlt_ftail1";}
    return (param->ftl_rec_mlt)[subitem];
  case 10: //par-10
  parname = "Rdet_rec_mlt_stail";
    return param->Stl_rec_mlt;
  case 11: //par-11
    if(subitem == 0) {parname = "Rdet_asc_mlt_ftail0";}
    if(subitem == 1) {parname = "Rdet_asc_mlt_ftail1";}
    return (param->ftl_asc_mlt)[subitem];
  case 12: //par-12
  parname = "Rdet_asc_mlt_stail";
    return param->Stl_asc_mlt;
  case 13: //par-13
  parname = "Rnp_fd_np_sgl0_Bzero";
    return (param->fd_np_sgl)[0][0];
  case 14: //par-14
  parname = "Rnp_fd_np_sgl1_Bzero";
    return (param->fd_np_sgl)[0][1];
  case 15: //par-15
  parname = "Rnp_fp_np_sgl_Bzero";
    return (param->fp_np_sgl)[0];
  case 16: //par-16
  parname = "Rnp_tau_np_p_sgl0_Bzero";
    return (param->tau_np_p_sgl)[0][0];
  case 17: //par-17
  parname = "Rnp_tau_np_p_sgl1_Bzero";
    return (param->tau_np_p_sgl)[0][1];
  case 18: //par-18
  parname = "Rnp_tau_np_n_sgl0_Bzero";
    return (param->tau_np_n_sgl)[0][0];
  case 19: //par-19
  parname = "Rnp_tau_np_n_sgl1_Bzero";
    return (param->tau_np_n_sgl)[0][1];
  case 20: //par-20
  parname = "Rnp_fd_np_mlt0_Bzero";
    return (param->fd_np_mlt)[0][0];
  case 21: //par-21
  parname = "Rnp_fd_np_mlt1_Bzero";
    return (param->fd_np_mlt)[0][1];
  case 22: //par-22
  parname = "Rnp_fd_np_st_mlt_Bzero";
    return (param->fd_np_st_mlt)[0];
  case 23: //par-23
  parname = "Rnp_fd_np_xi_mlt_Bzero";
    return (param->fd_np_xi_mlt)[0];
  case 24: //par-24
  parname = "Rnp_fd_np_stxi_mlt_Bzero";
    return (param->fd_np_stxi_mlt)[0];
  case 25: //par-25
  parname = "Rnp_fp_np_mlt_Bzero";
    return (param->fp_np_mlt)[0];
  case 26: //par-26
  parname = "Rnp_fn_np_mlt_Bzero";
    return (param->fn_np_mlt)[0];
  case 27: //par-27
  parname = "Rnp_tau_np_p_mlt0_Bzero";
    return (param->tau_np_p_mlt)[0][0];
  case 28: //par-28
  parname = "Rnp_tau_np_p_mlt1_Bzero";
    return (param->tau_np_p_mlt)[0][1];
  case 29: //par-29
  parname = "Rnp_tau_np_p_xi_mlt_Bzero";
    return (param->tau_np_p_xi_mlt)[0];
  case 30: //par-30
  parname = "Rnp_tau_np_p_stxi_mlt_Bzero";
    return (param->tau_np_p_stxi_mlt)[0];
  case 31: //par-31
  parname = "Rnp_tau_np_n_mlt0_Bzero";
    return (param->tau_np_n_mlt)[0][0];
  case 32: //par-32
  parname = "Rnp_tau_np_n_mlt1_Bzero";
    return (param->tau_np_n_mlt)[0][1];
  case 33: //par-33
  parname = "Rnp_tau_np_n_xi_mlt_Bzero";
    return (param->tau_np_n_xi_mlt)[0];
  case 34: //par-34
  parname = "Rnp_tau_np_n_stxi_mlt_Bzero";
    return (param->tau_np_n_stxi_mlt)[0];
  case 35: //par-35
  parname = "Rol_sig_ol";
    return param->sig_ol;
  case 36: //par-36
  parname = "Rol_fol_sgl";
    return param->fol_sgl;
  case 37: //par-37
  parname = "Rol_fol_mul";
    return param->fol_mul;
  case 38: //par-38
    if(subitem == 0) {
  parname = "Dt_cutoff_llmt";
    return Belle::dt_resol_global::dt_llmt;}
    if(subitem == 1) {
  parname = "Dt_cutoff_ulmt";
    return Belle::dt_resol_global::dt_ulmt;}

  default:
    return -999;

}


}














