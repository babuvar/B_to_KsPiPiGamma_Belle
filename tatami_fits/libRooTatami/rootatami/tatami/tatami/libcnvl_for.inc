C
C       libcnvl_for.inc
C         ---- Fortran Interface of livcnvl package
C        $Id: libcnvl_for.inc 8536 2003-08-06 06:43:33Z katayama $
C
C        $Log$
C        Revision 1.2  2003/08/06 06:43:33  katayama
C        New version from Nakadaira san
C
C        Revision 1.1  2002/09/05 01:27:07  katayama
C        New package tatami from Nakadaira/Sumisawa san
C
C        Revision 1.4  2002/07/05 06:16:27  nakadair
C        Add interfece to class dt_resol_global.
C
C        Revision 1.3  2002/05/11 07:51:37  nakadair
C        Update for solaris and GCC3.
C
C        Revision 1.2  2002/05/09 07:08:23  nakadair
C        Add add_outlier(...).
C
C        Revision 1.1  2002/05/07 21:38:23  nakadair
C        Add convolution with Resolution function and fortran interface.
C
C
      DOUBLE PRECISION SREC
      DOUBLE PRECISION FTMREC
      DOUBLE PRECISION STMREC
      DOUBLE PRECISION SASC
      DOUBLE PRECISION FTMASC
      DOUBLE PRECISION STMASC
      DOUBLE PRECISION SMNREC
      DOUBLE PRECISION STLREC
      DOUBLE PRECISION FTLREC
      DOUBLE PRECISION SMNASC
      DOUBLE PRECISION STLASC
      DOUBLE PRECISION FTLASC
      DOUBLE PRECISION FDNPSG
      DOUBLE PRECISION FPNPSG
      DOUBLE PRECISION TNPPSG
      DOUBLE PRECISION TNPNSG
      DOUBLE PRECISION FDNPML
      DOUBLE PRECISION FPNPML
      DOUBLE PRECISION TNPPML
      DOUBLE PRECISION TNPNML
      DOUBLE PRECISION TAUK1P
      DOUBLE PRECISION TAUK2P
      DOUBLE PRECISION TAUK1N
      DOUBLE PRECISION TAUK2N
      DOUBLE PRECISION SIGMAK
      DOUBLE PRECISION FK2
      DOUBLE PRECISION SIGOL
      DOUBLE PRECISION FOLSG
      DOUBLE PRECISION FOLML
	


      COMMON /RESPAR/SREC(2), FTMREC, STMREC,
     \               SASC(4), FTMASC, STMASC,
     \               SMNREC, STLREC, FTLREC,
     \               SMNASC, STLASC, FTLASC,
     \               FDNPSG(2), FPNPSG(2), TNPPSG(2,2), TNPNSG(2,2),
     \               FDNPML(2), FPNPML(2), TNPPML(2,2), TNPNML(2,2),
     \               TAUK1P(2), TAUK2P(2), TAUK1N(2), TAUK2N(2),
     \               SIGMAK(2), FK2, SIGOL, FOLSG, FOLML
      
      SAVE /RESPAR/
      

c
c    Decralation of functions
c      
      DOUBLE PRECISION DiracDelta
      DOUBLE PRECISION gaussian
      DOUBLE PRECISION norm_gaussian
      DOUBLE PRECISION norm_gaussian_w_cutoff
      DOUBLE PRECISION Ep
      DOUBLE PRECISION En
      DOUBLE PRECISION Ef
      DOUBLE PRECISION Enp
      DOUBLE PRECISION nEp_conv_gauss
      DOUBLE PRECISION nEn_conv_gauss
      DOUBLE PRECISION Ep_conv_gauss
      DOUBLE PRECISION En_conv_gauss
      DOUBLE PRECISION Ef_conv_gauss
      DOUBLE PRECISION Enp_conv_gauss
      DOUBLE PRECISION nMp
      DOUBLE PRECISION Mp
      DOUBLE PRECISION nMn
      DOUBLE PRECISION Mn
      DOUBLE PRECISION nMf
      DOUBLE PRECISION Mf
      DOUBLE PRECISION nAp
      DOUBLE PRECISION Ap
      DOUBLE PRECISION nAn
      DOUBLE PRECISION An
      DOUBLE PRECISION nAf
      DOUBLE PRECISION Af
      DOUBLE PRECISION nMp_conv_gauss
      DOUBLE PRECISION nMn_conv_gauss
      DOUBLE PRECISION nAp_conv_gauss
      DOUBLE PRECISION nAn_conv_gauss
      DOUBLE PRECISION Mp_conv_gauss
      DOUBLE PRECISION Mn_conv_gauss
      DOUBLE PRECISION Mf_conv_gauss
      DOUBLE PRECISION Ap_conv_gauss
      DOUBLE PRECISION An_conv_gauss
      DOUBLE PRECISION Af_conv_gauss
      DOUBLE PRECISION CoefEEs
      DOUBLE PRECISION CoefEEo
      DOUBLE PRECISION CoefEEf
      DOUBLE PRECISION sum_sigma
      DOUBLE PRECISION sum_sigma3
      DOUBLE PRECISION norm_nEp
      DOUBLE PRECISION norm_nEn
      DOUBLE PRECISION norm_nEf
      DOUBLE PRECISION norm_Ep
      DOUBLE PRECISION norm_En
      DOUBLE PRECISION norm_Ef
      DOUBLE PRECISION norm_Ap
      DOUBLE PRECISION norm_Mp
      DOUBLE PRECISION norm_An
      DOUBLE PRECISION norm_Mn
      DOUBLE PRECISION norm_Af
      DOUBLE PRECISION norm_Mf
      DOUBLE PRECISION norm_nEp_conv_gauss
      DOUBLE PRECISION norm_nEn_conv_gauss
      DOUBLE PRECISION norm_Ep_conv_gauss
      DOUBLE PRECISION norm_En_conv_gauss
      DOUBLE PRECISION norm_nEf_conv_gauss
      DOUBLE PRECISION norm_Ef_conv_gauss
      DOUBLE PRECISION norm_Enp_conv_gauss
      DOUBLE PRECISION norm_nAn_conv_gauss
      DOUBLE PRECISION norm_nAp_conv_gauss
      DOUBLE PRECISION norm_nAf_conv_gauss
      DOUBLE PRECISION norm_nMn_conv_gauss
      DOUBLE PRECISION norm_nMp_conv_gauss
      DOUBLE PRECISION norm_nMf_conv_gauss
      DOUBLE PRECISION norm_An_conv_gauss
      DOUBLE PRECISION norm_Ap_conv_gauss
      DOUBLE PRECISION norm_Af_conv_gauss
      DOUBLE PRECISION norm_Mn_conv_gauss
      DOUBLE PRECISION norm_Mp_conv_gauss
      DOUBLE PRECISION norm_Mf_conv_gauss
      DOUBLE PRECISION Ap_conv_Ep
      DOUBLE PRECISION An_conv_En
      DOUBLE PRECISION Mp_conv_Ep
      DOUBLE PRECISION Mn_conv_En
      DOUBLE PRECISION Ap_conv_En
      DOUBLE PRECISION An_conv_Ep
      DOUBLE PRECISION Mp_conv_En
      DOUBLE PRECISION Mn_conv_Ep
      DOUBLE PRECISION Ap_conv_Ef
      DOUBLE PRECISION An_conv_Ef
      DOUBLE PRECISION Mp_conv_Ef
      DOUBLE PRECISION Mn_conv_Ef
      DOUBLE PRECISION Af_conv_Ep
      DOUBLE PRECISION Af_conv_En
      DOUBLE PRECISION Mf_conv_Ep
      DOUBLE PRECISION Mf_conv_En
      DOUBLE PRECISION Af_conv_Ef
      DOUBLE PRECISION Mf_conv_Ef
      DOUBLE PRECISION Ef_conv_EnpGauss
      INTEGER Ep_conv_En_coef
      INTEGER En_conv_Ep_coef
      INTEGER Ep_conv_Ep_coef
      INTEGER En_conv_En_coef
      DOUBLE PRECISION Ep_conv_En
      DOUBLE PRECISION En_conv_Ep
      DOUBLE PRECISION xEp
      DOUBLE PRECISION xEn
      DOUBLE PRECISION xEf
      DOUBLE PRECISION norm_xEp
      DOUBLE PRECISION norm_xEn
      DOUBLE PRECISION norm_xEf
      DOUBLE PRECISION xEp_conv_gauss
      DOUBLE PRECISION xEn_conv_gauss
      DOUBLE PRECISION xEf_conv_gauss
      DOUBLE PRECISION xEp_conv_gauss_w_meanshift
      DOUBLE PRECISION xEn_conv_gauss_w_meanshift
      DOUBLE PRECISION xEf_conv_gauss_w_meanshift
      DOUBLE PRECISION xxEp
      DOUBLE PRECISION xxEn
      DOUBLE PRECISION xxEf
      DOUBLE PRECISION norm_xxEp
      DOUBLE PRECISION norm_xxEn
      DOUBLE PRECISION norm_xxEf
      DOUBLE PRECISION xxEp_conv_gauss
      DOUBLE PRECISION xxEn_conv_gauss
      DOUBLE PRECISION xxEf_conv_gauss
      DOUBLE PRECISION xxEp_conv_gauss_w_meanshift
      DOUBLE PRECISION xxEn_conv_gauss_w_meanshift
      DOUBLE PRECISION xxEf_conv_gauss_w_meanshift
      DOUBLE PRECISION Ep_conv_Ep
      DOUBLE PRECISION En_conv_En
      DOUBLE PRECISION Ep_conv_Ef
      DOUBLE PRECISION En_conv_Ef
      DOUBLE PRECISION Ef_conv_Ep
      DOUBLE PRECISION Ef_conv_En
      DOUBLE PRECISION Ef_conv_Ef
      DOUBLE PRECISION xEp_conv_Ep
      DOUBLE PRECISION xEn_conv_En
      DOUBLE PRECISION xEp_conv_En
      DOUBLE PRECISION xEn_conv_Ep
      DOUBLE PRECISION xEf_conv_Ep
      DOUBLE PRECISION xEf_conv_En
      DOUBLE PRECISION xEf_conv_Ef
      DOUBLE PRECISION xEp_conv_Ef
      DOUBLE PRECISION xEn_conv_Ef
      DOUBLE PRECISION Ep_conv_Ep_Ep
      DOUBLE PRECISION En_conv_En_Ep
      DOUBLE PRECISION Ep_conv_En_Ep
      DOUBLE PRECISION En_conv_Ep_Ep
      DOUBLE PRECISION Ef_conv_Ep_Ep
      DOUBLE PRECISION Ef_conv_En_Ep
      DOUBLE PRECISION Ef_conv_Ef_Ep
      DOUBLE PRECISION Ep_conv_Ep_En
      DOUBLE PRECISION En_conv_En_En
      DOUBLE PRECISION Ep_conv_En_En
      DOUBLE PRECISION En_conv_Ep_En
      DOUBLE PRECISION Ef_conv_Ep_En
      DOUBLE PRECISION Ef_conv_En_En
      DOUBLE PRECISION Ef_conv_Ef_En
      DOUBLE PRECISION Ep_conv_Ep_Ef
      DOUBLE PRECISION En_conv_En_Ef
      DOUBLE PRECISION Ep_conv_En_Ef
      DOUBLE PRECISION En_conv_Ep_Ef
      DOUBLE PRECISION Ef_conv_Ep_Ef
      DOUBLE PRECISION Ef_conv_En_Ef
      DOUBLE PRECISION Ef_conv_Ef_Ef
      DOUBLE PRECISION xEp_conv_Ep_Ef
      DOUBLE PRECISION xEn_conv_En_Ef
      DOUBLE PRECISION xEp_conv_En_Ef
      DOUBLE PRECISION xEn_conv_Ep_Ef
      DOUBLE PRECISION xEf_conv_Ep_Ef
      DOUBLE PRECISION xEf_conv_En_Ef
      DOUBLE PRECISION xEf_conv_Ef_Ef
      DOUBLE PRECISION Ep_conv_Ep_gauss
      DOUBLE PRECISION En_conv_En_gauss
      DOUBLE PRECISION Ep_conv_En_gauss
      DOUBLE PRECISION En_conv_Ep_gauss
      DOUBLE PRECISION Ef_conv_Ep_gauss
      DOUBLE PRECISION Ef_conv_En_gauss
      DOUBLE PRECISION Ef_conv_Ef_gauss
      DOUBLE PRECISION xEp_conv_Ep_gauss
      DOUBLE PRECISION xEn_conv_En_gauss
      DOUBLE PRECISION xEp_conv_En_gauss
      DOUBLE PRECISION xEn_conv_Ep_gauss
      DOUBLE PRECISION xEf_conv_Ep_gauss
      DOUBLE PRECISION xEf_conv_En_gauss
      DOUBLE PRECISION xEf_conv_Ef_gauss
      DOUBLE PRECISION Ep_conv_Ep_Ep_gauss
      DOUBLE PRECISION En_conv_En_Ep_gauss
      DOUBLE PRECISION Ep_conv_En_Ep_gauss
      DOUBLE PRECISION En_conv_Ep_Ep_gauss
      DOUBLE PRECISION Ef_conv_Ep_Ep_gauss
      DOUBLE PRECISION Ef_conv_En_Ep_gauss
      DOUBLE PRECISION Ef_conv_Ef_Ep_gauss
      DOUBLE PRECISION Ep_conv_Ep_En_gauss
      DOUBLE PRECISION En_conv_En_En_gauss
      DOUBLE PRECISION Ep_conv_En_En_gauss
      DOUBLE PRECISION En_conv_Ep_En_gauss
      DOUBLE PRECISION Ef_conv_Ep_En_gauss
      DOUBLE PRECISION Ef_conv_En_En_gauss
      DOUBLE PRECISION Ef_conv_Ef_En_gauss
      DOUBLE PRECISION norm_xEp_conv_gauss
      DOUBLE PRECISION norm_xEn_conv_gauss
      DOUBLE PRECISION norm_xEf_conv_gauss
      DOUBLE PRECISION norm_xxEp_conv_gauss
      DOUBLE PRECISION norm_xxEn_conv_gauss
      DOUBLE PRECISION norm_xxEf_conv_gauss
      DOUBLE PRECISION norm_Ep_conv_Ep
      DOUBLE PRECISION norm_En_conv_En
      DOUBLE PRECISION norm_Ep_conv_En
      DOUBLE PRECISION norm_En_conv_Ep
      DOUBLE PRECISION norm_Ef_conv_Ep
      DOUBLE PRECISION norm_Ef_conv_En
      DOUBLE PRECISION norm_Ef_conv_Ef
      DOUBLE PRECISION norm_xEp_conv_Ep
      DOUBLE PRECISION norm_xEn_conv_En
      DOUBLE PRECISION norm_xEp_conv_En
      DOUBLE PRECISION norm_xEn_conv_Ep
      DOUBLE PRECISION norm_xEf_conv_Ep
      DOUBLE PRECISION norm_xEf_conv_En
      DOUBLE PRECISION norm_xEf_conv_Ef
      DOUBLE PRECISION norm_Ep_conv_Ep_Ep
      DOUBLE PRECISION norm_En_conv_En_Ep
      DOUBLE PRECISION norm_Ep_conv_En_Ep
      DOUBLE PRECISION norm_En_conv_Ep_Ep
      DOUBLE PRECISION norm_Ef_conv_Ep_Ep
      DOUBLE PRECISION norm_Ef_conv_En_Ep
      DOUBLE PRECISION norm_Ef_conv_Ef_Ep
      DOUBLE PRECISION norm_Ep_conv_Ep_En
      DOUBLE PRECISION norm_En_conv_En_En
      DOUBLE PRECISION norm_Ep_conv_En_En
      DOUBLE PRECISION norm_En_conv_Ep_En
      DOUBLE PRECISION norm_Ef_conv_Ep_En
      DOUBLE PRECISION norm_Ef_conv_En_En
      DOUBLE PRECISION norm_Ef_conv_Ef_En
      DOUBLE PRECISION norm_Ep_conv_Ep_gauss
      DOUBLE PRECISION norm_En_conv_En_gauss
      DOUBLE PRECISION norm_Ep_conv_En_gauss
      DOUBLE PRECISION norm_En_conv_Ep_gauss
      DOUBLE PRECISION norm_Ef_conv_Ep_gauss
      DOUBLE PRECISION norm_Ef_conv_En_gauss
      DOUBLE PRECISION norm_Ef_conv_Ef_gauss
      DOUBLE PRECISION norm_xEp_conv_Ep_gauss
      DOUBLE PRECISION norm_xEn_conv_En_gauss
      DOUBLE PRECISION norm_xEp_conv_En_gauss
      DOUBLE PRECISION norm_xEn_conv_Ep_gauss
      DOUBLE PRECISION norm_xEf_conv_Ep_gauss
      DOUBLE PRECISION norm_xEf_conv_En_gauss
      DOUBLE PRECISION norm_xEf_conv_Ef_gauss
      DOUBLE PRECISION norm_Ep_conv_Ep_Ep_gauss
      DOUBLE PRECISION norm_En_conv_En_Ep_gauss
      DOUBLE PRECISION norm_Ep_conv_En_Ep_gauss
      DOUBLE PRECISION norm_En_conv_Ep_Ep_gauss
      DOUBLE PRECISION norm_Ef_conv_Ep_Ep_gauss
      DOUBLE PRECISION norm_Ef_conv_En_Ep_gauss
      DOUBLE PRECISION norm_Ef_conv_Ef_Ep_gauss
      DOUBLE PRECISION norm_Ep_conv_Ep_En_gauss
      DOUBLE PRECISION norm_En_conv_En_En_gauss
      DOUBLE PRECISION norm_Ep_conv_En_En_gauss
      DOUBLE PRECISION norm_En_conv_Ep_En_gauss
      DOUBLE PRECISION norm_Ef_conv_Ep_En_gauss
      DOUBLE PRECISION norm_Ef_conv_En_En_gauss
      DOUBLE PRECISION norm_Ef_conv_Ef_En_gauss
      DOUBLE PRECISION Rk_fullrec
      DOUBLE PRECISION Rk_partial
      DOUBLE PRECISION Rnp
      DOUBLE PRECISION Rasc
      DOUBLE PRECISION RascRnp
      DOUBLE PRECISION Rdet
      DOUBLE PRECISION RdetRnp
      DOUBLE PRECISION EfRk_fullrec
      DOUBLE PRECISION AfRk_fullrec
      DOUBLE PRECISION MfRk_fullrec
      DOUBLE PRECISION EfRk_partial
      DOUBLE PRECISION AfRk_partial
      DOUBLE PRECISION MfRk_partial
      DOUBLE PRECISION EfRkRdet_fullrec
      DOUBLE PRECISION AfRkRdet_fullrec
      DOUBLE PRECISION MfRkRdet_fullrec
      DOUBLE PRECISION EfRkRdet_partial
      DOUBLE PRECISION AfRkRdet_partial
      DOUBLE PRECISION MfRkRdet_partial
      DOUBLE PRECISION EfRkRdetRnp_fullrec
      DOUBLE PRECISION AfRkRdetRnp_fullrec
      DOUBLE PRECISION MfRkRdetRnp_fullrec
      DOUBLE PRECISION EfRkRdetRnp_partial
      DOUBLE PRECISION AfRkRdetRnp_partial
      DOUBLE PRECISION MfRkRdetRnp_partial
      DOUBLE PRECISION norm_Rk_partial
      DOUBLE PRECISION norm_Rasc
      DOUBLE PRECISION norm_Rrec
      DOUBLE PRECISION norm_Rnp
      DOUBLE PRECISION norm_RascRnp
      DOUBLE PRECISION norm_Rdet
      DOUBLE PRECISION norm_RdetRnp
      DOUBLE PRECISION norm_EfRk_fullrec
      DOUBLE PRECISION norm_AfRk_fullrec
      DOUBLE PRECISION norm_MfRk_fullrec
      DOUBLE PRECISION norm_EfRk_partial
      DOUBLE PRECISION norm_AfRk_partial
      DOUBLE PRECISION norm_MfRk_partial
      DOUBLE PRECISION norm_EfRkRdet_fullrec
      DOUBLE PRECISION norm_AfRkRdet_fullrec
      DOUBLE PRECISION norm_MfRkRdet_fullrec
      DOUBLE PRECISION norm_EfRkRdet_partial
      DOUBLE PRECISION norm_AfRkRdet_partial
      DOUBLE PRECISION norm_MfRkRdet_partial
      DOUBLE PRECISION norm_EfRkRdetRnp_fullrec
      DOUBLE PRECISION norm_AfRkRdetRnp_fullrec
      DOUBLE PRECISION norm_MfRkRdetRnp_fullrec
      DOUBLE PRECISION norm_EfRkRdetRnp_partial
      DOUBLE PRECISION norm_AfRkRdetRnp_partial
      DOUBLE PRECISION norm_MfRkRdetRnp_partial
      INTEGER cpdefpar
      DOUBLE PRECISION add_outlier
      DOUBLE PRECISION add_norm_outlier
c
      EXTERNAL DiracDelta
      EXTERNAL gaussian
      EXTERNAL norm_gaussian
      EXTERNAL norm_gaussian_w_cutoff
      EXTERNAL Ep
      EXTERNAL En
      EXTERNAL Ef
      EXTERNAL Enp
      EXTERNAL nEp_conv_gauss
      EXTERNAL nEn_conv_gauss
      EXTERNAL Ep_conv_gauss
      EXTERNAL En_conv_gauss
      EXTERNAL Ef_conv_gauss
      EXTERNAL Enp_conv_gauss
      EXTERNAL nMp
      EXTERNAL Mp
      EXTERNAL nMn
      EXTERNAL Mn
      EXTERNAL nMf
      EXTERNAL Mf
      EXTERNAL nAp
      EXTERNAL Ap
      EXTERNAL nAn
      EXTERNAL An
      EXTERNAL nAf
      EXTERNAL Af
      EXTERNAL nMp_conv_gauss
      EXTERNAL nMn_conv_gauss
      EXTERNAL nAp_conv_gauss
      EXTERNAL nAn_conv_gauss
      EXTERNAL Mp_conv_gauss
      EXTERNAL Mn_conv_gauss
      EXTERNAL Mf_conv_gauss
      EXTERNAL Ap_conv_gauss
      EXTERNAL An_conv_gauss
      EXTERNAL Af_conv_gauss
      EXTERNAL CoefEEs
      EXTERNAL CoefEEo
      EXTERNAL CoefEEf
      EXTERNAL sum_sigma
      EXTERNAL sum_sigma3
      EXTERNAL norm_nEp
      EXTERNAL norm_nEn
      EXTERNAL norm_nEf
      EXTERNAL norm_Ep
      EXTERNAL norm_En
      EXTERNAL norm_Ef
      EXTERNAL norm_Ap
      EXTERNAL norm_Mp
      EXTERNAL norm_An
      EXTERNAL norm_Mn
      EXTERNAL norm_Af
      EXTERNAL norm_Mf
      EXTERNAL norm_nEp_conv_gauss
      EXTERNAL norm_nEn_conv_gauss
      EXTERNAL norm_Ep_conv_gauss
      EXTERNAL norm_En_conv_gauss
      EXTERNAL norm_nEf_conv_gauss
      EXTERNAL norm_Ef_conv_gauss
      EXTERNAL norm_Enp_conv_gauss
      EXTERNAL norm_nAn_conv_gauss
      EXTERNAL norm_nAp_conv_gauss
      EXTERNAL norm_nAf_conv_gauss
      EXTERNAL norm_nMn_conv_gauss
      EXTERNAL norm_nMp_conv_gauss
      EXTERNAL norm_nMf_conv_gauss
      EXTERNAL norm_An_conv_gauss
      EXTERNAL norm_Ap_conv_gauss
      EXTERNAL norm_Af_conv_gauss
      EXTERNAL norm_Mn_conv_gauss
      EXTERNAL norm_Mp_conv_gauss
      EXTERNAL norm_Mf_conv_gauss
      EXTERNAL Ap_conv_Ep
      EXTERNAL An_conv_En
      EXTERNAL Mp_conv_Ep
      EXTERNAL Mn_conv_En
      EXTERNAL Ap_conv_En
      EXTERNAL An_conv_Ep
      EXTERNAL Mp_conv_En
      EXTERNAL Mn_conv_Ep
      EXTERNAL Ap_conv_Ef
      EXTERNAL An_conv_Ef
      EXTERNAL Mp_conv_Ef
      EXTERNAL Mn_conv_Ef
      EXTERNAL Af_conv_Ep
      EXTERNAL Af_conv_En
      EXTERNAL Mf_conv_Ep
      EXTERNAL Mf_conv_En
      EXTERNAL Af_conv_Ef
      EXTERNAL Mf_conv_Ef
      EXTERNAL Ef_conv_EnpGauss
      EXTERNAL Ep_conv_En_coef
      EXTERNAL En_conv_Ep_coef
      EXTERNAL Ep_conv_Ep_coef
      EXTERNAL En_conv_En_coef
      EXTERNAL Ep_conv_En
      EXTERNAL En_conv_Ep
      EXTERNAL xEp
      EXTERNAL xEn
      EXTERNAL xEf
      EXTERNAL norm_xEp
      EXTERNAL norm_xEn
      EXTERNAL norm_xEf
      EXTERNAL xEp_conv_gauss
      EXTERNAL xEn_conv_gauss
      EXTERNAL xEf_conv_gauss
      EXTERNAL xEp_conv_gauss_w_meanshift
      EXTERNAL xEn_conv_gauss_w_meanshift
      EXTERNAL xEf_conv_gauss_w_meanshift
      EXTERNAL xxEp
      EXTERNAL xxEn
      EXTERNAL xxEf
      EXTERNAL norm_xxEp
      EXTERNAL norm_xxEn
      EXTERNAL norm_xxEf
      EXTERNAL xxEp_conv_gauss
      EXTERNAL xxEn_conv_gauss
      EXTERNAL xxEf_conv_gauss
      EXTERNAL xxEp_conv_gauss_w_meanshift
      EXTERNAL xxEn_conv_gauss_w_meanshift
      EXTERNAL xxEf_conv_gauss_w_meanshift
      EXTERNAL Ep_conv_Ep
      EXTERNAL En_conv_En
      EXTERNAL Ep_conv_Ef
      EXTERNAL En_conv_Ef
      EXTERNAL Ef_conv_Ep
      EXTERNAL Ef_conv_En
      EXTERNAL Ef_conv_Ef
      EXTERNAL xEp_conv_Ep
      EXTERNAL xEn_conv_En
      EXTERNAL xEp_conv_En
      EXTERNAL xEn_conv_Ep
      EXTERNAL xEf_conv_Ep
      EXTERNAL xEf_conv_En
      EXTERNAL xEf_conv_Ef
      EXTERNAL xEp_conv_Ef
      EXTERNAL xEn_conv_Ef
      EXTERNAL Ep_conv_Ep_Ep
      EXTERNAL En_conv_En_Ep
      EXTERNAL Ep_conv_En_Ep
      EXTERNAL En_conv_Ep_Ep
      EXTERNAL Ef_conv_Ep_Ep
      EXTERNAL Ef_conv_En_Ep
      EXTERNAL Ef_conv_Ef_Ep
      EXTERNAL Ep_conv_Ep_En
      EXTERNAL En_conv_En_En
      EXTERNAL Ep_conv_En_En
      EXTERNAL En_conv_Ep_En
      EXTERNAL Ef_conv_Ep_En
      EXTERNAL Ef_conv_En_En
      EXTERNAL Ef_conv_Ef_En
      EXTERNAL Ep_conv_Ep_Ef
      EXTERNAL En_conv_En_Ef
      EXTERNAL Ep_conv_En_Ef
      EXTERNAL En_conv_Ep_Ef
      EXTERNAL Ef_conv_Ep_Ef
      EXTERNAL Ef_conv_En_Ef
      EXTERNAL Ef_conv_Ef_Ef
      EXTERNAL xEp_conv_Ep_Ef
      EXTERNAL xEn_conv_En_Ef
      EXTERNAL xEp_conv_En_Ef
      EXTERNAL xEn_conv_Ep_Ef
      EXTERNAL xEf_conv_Ep_Ef
      EXTERNAL xEf_conv_En_Ef
      EXTERNAL xEf_conv_Ef_Ef
      EXTERNAL Ep_conv_Ep_gauss
      EXTERNAL En_conv_En_gauss
      EXTERNAL Ep_conv_En_gauss
      EXTERNAL En_conv_Ep_gauss
      EXTERNAL Ef_conv_Ep_gauss
      EXTERNAL Ef_conv_En_gauss
      EXTERNAL Ef_conv_Ef_gauss
      EXTERNAL xEp_conv_Ep_gauss
      EXTERNAL xEn_conv_En_gauss
      EXTERNAL xEp_conv_En_gauss
      EXTERNAL xEn_conv_Ep_gauss
      EXTERNAL xEf_conv_Ep_gauss
      EXTERNAL xEf_conv_En_gauss
      EXTERNAL xEf_conv_Ef_gauss
      EXTERNAL Ep_conv_Ep_Ep_gauss
      EXTERNAL En_conv_En_Ep_gauss
      EXTERNAL Ep_conv_En_Ep_gauss
      EXTERNAL En_conv_Ep_Ep_gauss
      EXTERNAL Ef_conv_Ep_Ep_gauss
      EXTERNAL Ef_conv_En_Ep_gauss
      EXTERNAL Ef_conv_Ef_Ep_gauss
      EXTERNAL Ep_conv_Ep_En_gauss
      EXTERNAL En_conv_En_En_gauss
      EXTERNAL Ep_conv_En_En_gauss
      EXTERNAL En_conv_Ep_En_gauss
      EXTERNAL Ef_conv_Ep_En_gauss
      EXTERNAL Ef_conv_En_En_gauss
      EXTERNAL Ef_conv_Ef_En_gauss
      EXTERNAL norm_xEp_conv_gauss
      EXTERNAL norm_xEn_conv_gauss
      EXTERNAL norm_xEf_conv_gauss
      EXTERNAL norm_xxEp_conv_gauss
      EXTERNAL norm_xxEn_conv_gauss
      EXTERNAL norm_xxEf_conv_gauss
      EXTERNAL norm_Ep_conv_Ep
      EXTERNAL norm_En_conv_En
      EXTERNAL norm_Ep_conv_En
      EXTERNAL norm_En_conv_Ep
      EXTERNAL norm_Ef_conv_Ep
      EXTERNAL norm_Ef_conv_En
      EXTERNAL norm_Ef_conv_Ef
      EXTERNAL norm_xEp_conv_Ep
      EXTERNAL norm_xEn_conv_En
      EXTERNAL norm_xEp_conv_En
      EXTERNAL norm_xEn_conv_Ep
      EXTERNAL norm_xEf_conv_Ep
      EXTERNAL norm_xEf_conv_En
      EXTERNAL norm_xEf_conv_Ef
      EXTERNAL norm_Ep_conv_Ep_Ep
      EXTERNAL norm_En_conv_En_Ep
      EXTERNAL norm_Ep_conv_En_Ep
      EXTERNAL norm_En_conv_Ep_Ep
      EXTERNAL norm_Ef_conv_Ep_Ep
      EXTERNAL norm_Ef_conv_En_Ep
      EXTERNAL norm_Ef_conv_Ef_Ep
      EXTERNAL norm_Ep_conv_Ep_En
      EXTERNAL norm_En_conv_En_En
      EXTERNAL norm_Ep_conv_En_En
      EXTERNAL norm_En_conv_Ep_En
      EXTERNAL norm_Ef_conv_Ep_En
      EXTERNAL norm_Ef_conv_En_En
      EXTERNAL norm_Ef_conv_Ef_En
      EXTERNAL norm_Ep_conv_Ep_gauss
      EXTERNAL norm_En_conv_En_gauss
      EXTERNAL norm_Ep_conv_En_gauss
      EXTERNAL norm_En_conv_Ep_gauss
      EXTERNAL norm_Ef_conv_Ep_gauss
      EXTERNAL norm_Ef_conv_En_gauss
      EXTERNAL norm_Ef_conv_Ef_gauss
      EXTERNAL norm_xEp_conv_Ep_gauss
      EXTERNAL norm_xEn_conv_En_gauss
      EXTERNAL norm_xEp_conv_En_gauss
      EXTERNAL norm_xEn_conv_Ep_gauss
      EXTERNAL norm_xEf_conv_Ep_gauss
      EXTERNAL norm_xEf_conv_En_gauss
      EXTERNAL norm_xEf_conv_Ef_gauss
      EXTERNAL norm_Ep_conv_Ep_Ep_gauss
      EXTERNAL norm_En_conv_En_Ep_gauss
      EXTERNAL norm_Ep_conv_En_Ep_gauss
      EXTERNAL norm_En_conv_Ep_Ep_gauss
      EXTERNAL norm_Ef_conv_Ep_Ep_gauss
      EXTERNAL norm_Ef_conv_En_Ep_gauss
      EXTERNAL norm_Ef_conv_Ef_Ep_gauss
      EXTERNAL norm_Ep_conv_Ep_En_gauss
      EXTERNAL norm_En_conv_En_En_gauss
      EXTERNAL norm_Ep_conv_En_En_gauss
      EXTERNAL norm_En_conv_Ep_En_gauss
      EXTERNAL norm_Ef_conv_Ep_En_gauss
      EXTERNAL norm_Ef_conv_En_En_gauss
      EXTERNAL norm_Ef_conv_Ef_En_gauss
      EXTERNAL Rk_fullrec
      EXTERNAL Rk_partial
      EXTERNAL Rnp
      EXTERNAL Rasc
      EXTERNAL RascRnp
      EXTERNAL Rdet
      EXTERNAL RdetRnp
      EXTERNAL EfRk_fullrec
      EXTERNAL AfRk_fullrec
      EXTERNAL MfRk_fullrec
      EXTERNAL EfRk_partial
      EXTERNAL AfRk_partial
      EXTERNAL MfRk_partial
      EXTERNAL EfRkRdet_fullrec
      EXTERNAL AfRkRdet_fullrec
      EXTERNAL MfRkRdet_fullrec
      EXTERNAL EfRkRdet_partial
      EXTERNAL AfRkRdet_partial
      EXTERNAL MfRkRdet_partial
      EXTERNAL EfRkRdetRnp_fullrec
      EXTERNAL AfRkRdetRnp_fullrec
      EXTERNAL MfRkRdetRnp_fullrec
      EXTERNAL EfRkRdetRnp_partial
      EXTERNAL AfRkRdetRnp_partial
      EXTERNAL MfRkRdetRnp_partial
      EXTERNAL norm_Rk_partial
      EXTERNAL norm_Rasc
      EXTERNAL norm_Rrec
      EXTERNAL norm_Rnp
      EXTERNAL norm_RascRnp
      EXTERNAL norm_Rdet
      EXTERNAL norm_RdetRnp
      EXTERNAL norm_EfRk_fullrec
      EXTERNAL norm_AfRk_fullrec
      EXTERNAL norm_MfRk_fullrec
      EXTERNAL norm_EfRk_partial
      EXTERNAL norm_AfRk_partial
      EXTERNAL norm_MfRk_partial
      EXTERNAL norm_EfRkRdet_fullrec
      EXTERNAL norm_AfRkRdet_fullrec
      EXTERNAL norm_MfRkRdet_fullrec
      EXTERNAL norm_EfRkRdet_partial
      EXTERNAL norm_AfRkRdet_partial
      EXTERNAL norm_MfRkRdet_partial
      EXTERNAL norm_EfRkRdetRnp_fullrec
      EXTERNAL norm_AfRkRdetRnp_fullrec
      EXTERNAL norm_MfRkRdetRnp_fullrec
      EXTERNAL norm_EfRkRdetRnp_partial
      EXTERNAL norm_AfRkRdetRnp_partial
      EXTERNAL norm_MfRkRdetRnp_partial
      EXTERNAL cpdefpar
      EXTERNAL add_outlier
      EXTERNAL add_norm_outlier
      EXTERNAL dt_resol_global_set_bgc
      EXTERNAL dt_resol_global_set_bmass
      EXTERNAL dt_resol_global_set_b0_mass
      EXTERNAL dt_resol_global_set_bplus_mass
