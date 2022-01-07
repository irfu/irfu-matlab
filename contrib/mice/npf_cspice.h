/*

-Header_File npf_cspice.h ( Mice routine name pointer mapping structure )

-Abstract

   List the Mice call name to CSPICE call mappings.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading

   MICE.REQ

-Keywords

-Brief_I/O

   None.

-Detailed_Input

   None.

-Detailed_Output

   None.

-Parameters
*/

#ifndef npf_cspice_h
#define npf_cspice_h

/*
-Exceptions

   None.

-Files

   None.

-Particulars

   Below, the "npf" structure and the assignments within. This
   structure defines mapping between the routine name field "name"
   to the pointer corresponding to that routine's interface.

   A "name" with format name_c maps to a cspice_name call, a "name"
   name_s maps to a mice_name and cspice_name call from MATLAB where
   "name" represents the call base name.

   Update this list to add a new interface call to Mice.

-Examples

   None.

-Restrictions

   None.

-Literature_References

   None.

-Author_and_Institution

   N.J. Bachman        (JPL)
   G. Chinn            (JPL)
   M. Costa Sitja      (JPL)
   J. Diaz del Rio     (ODC Space)
   E.D. Wright         (JPL)

-Version

   -Mice 1.5.0 05-NOV-2021 (JDR) (MCS)

      Added interfaces:

         cspice_azlcpo
         cspice_azlrec
         cspice_badkpv
         cspice_bltfrm
         cspice_chbder
         cspice_chbigr
         cspice_chbint
         cspice_chbval
         cspice_ckfrot
         cspice_ckfxfm
         cspice_ckgr02
         cspice_ckgr03
         cspice_cklpf
         cspice_ckmeta
         cspice_cknr02
         cspice_cknr03
         cspice_ckupf
         cspice_dafhsf
         cspice_dafps
         cspice_dafrs
         cspice_dasadc
         cspice_dasadd
         cspice_dasadi
         cspice_dashfs
         cspice_daslla
         cspice_dasllc
         cspice_dasonw
         cspice_dasops
         cspice_dasopw
         cspice_dasrdc
         cspice_dasrdd
         cspice_dasrdi
         cspice_dasudc
         cspice_dasudd
         cspice_dasudi
         cspice_daswbr
         cspice_dazldr
         cspice_dlabbs
         cspice_dlabns
         cspice_dlaens
         cspice_dlafps
         cspice_dlaopn
         cspice_dnearp
         cspice_dpmax
         cspice_dpmin
         cspice_drdazl
         cspice_ednmpt
         cspice_edpnt
         cspice_evsgp4
         cspice_expool
         cspice_getelm
         cspice_getfat
         cspice_getfvn
         cspice_hrmesp
         cspice_hrmint
         cspice_intmax
         cspice_intmin
         cspice_invstm
         cspice_kplfrm
         cspice_ldpool
         cspice_lgresp
         cspice_lgrind
         cspice_lgrint
         cspice_nextwd
         cspice_nthwd
         cspice_oscltx
         cspice_polyds
         cspice_prop2b
         cspice_qderiv
         cspice_qdq2av
         cspice_qxq
         cspice_recazl
         cspice_repmc
         cspice_repmct
         cspice_repmd
         cspice_repmf
         cspice_repmi
         cspice_repml
         cspice_repmot
         cspice_rotvec
         cspice_scfmt
         cspice_scpart
         cspice_spkapo (mice_spkapo)
         cspice_spkez  (mice_spkez)
         cspice_spkgeo (mice_spkgeo)
         cspice_spklef
         cspice_spkssb
         cspice_spkuef
         cspice_spkw09
         cspice_spkw10
         cspice_spkw13
         cspice_stelab
         cspice_stlabx
         cspice_surfpv
         cspice_szpool
         cspice_tangpt
         cspice_tipbod
         cspice_tisbod
         cspice_tkfram
         cspice_tparch
         cspice_tparse
         cspice_trgsep
         cspice_twovxf
         cspice_vprojg
         cspice_vupack

   -Mice 1.4.0 05-JAN-2017 (EDW) (NJB)

      Added interfaces:

         cspice_bodfnd
         cspice_ccifrm (mice_ccifrm)
         cspice_dascls
         cspice_dasec
         cspice_dasopr
         cspice_dlabfs
         cspice_dlafns
         cspice_dskb02
         cspice_dskcls
         cspice_dskd02
         cspice_dskgd
         cspice_dskgtl
         cspice_dski02
         cspice_dskmi2
         cspice_dskn02
         cspice_dskobj
         cspice_dskopn
         cspice_dskp02
         cspice_dskrb2
         cspice_dsksrf
         cspice_dskstl
         cspice_dskv02
         cspice_dskw02
         cspice_dskx02
         cspice_dskxsi
         cspice_dskz02
         cspice_illum_pl02
         cspice_illum_plid_pl02 (mice_illum_plid_pl02)
         cspice_illumf (mice_illumf)
         cspice_illumg (mice_illumg)
         cspice_latsrf
         cspice_limb_pl02
         cspice_limbpt
         cspice_llgrid_pl02
         cspice_lspcn
         cspice_pckcov
         cspice_pckfrm
         cspice_pltar
         cspice_pltexp
         cspice_pltnp
         cspice_pltnrm
         cspice_pltvol
         cspice_srfc2s (mice_srfc2s)
         cspice_srfcss (mice_srfcss)
         cspice_srfnrm
         cspice_srfrec
         cspice_srfs2c (mice_srfs2c)
         cspice_srfscc (mice_srfscc)
         cspice_subpt_pl02
         cspice_subsol_pl02
         cspice_term_pl02
         cspice_termpt

   -Mice 1.3.0 31-OCT-2012 (EDW)(SCK)

      Added interfaces:

         cspice_cidfrm (mice_cidfrm)
         cspice_cnmfrm (mice_cnmfrm)
         cspice_dafac
         cspice_dafbbs
         cspice_dafbfs
         cspice_dafcls
         cspice_dafcs
         cspice_dafdc
         cspice_dafec
         cspice_daffna
         cspice_daffpa
         cspice_dafgda
         cspice_dafgn
         cspice_dafgs
         cspice_dafopr
         cspice_dafopw
         cspice_dafus
         cspice_dcyldr
         cspice_dgeodr
         cspice_dlatdr
         cspice_dpgrdr
         cspice_drdcyl
         cspice_drdgeo
         cspice_drdlat
         cspice_drdpgr
         cspice_drdsph
         cspice_dsphdr
         cspice_dvpool
         cspice_edlimb (mice_edlimb)
         cspice_edterm
         cspice_fovray
         cspice_fovtrg
         cspice_frame
         cspice_frinfo (mice_frinfo)
         cspice_frmnam
         cspice_gfilum
         cspice_gfpa
         cspice_gfstol
         cspice_inedpl
         cspice_inelpl
         cspice_inrypl
         cspice_invort
         cspice_namfrm
         cspice_npedln (mice_npedln)
         cspice_npelpt (mice_npelpt)
         cspice_nplnpt (mice_nplnpt)
         cspice_occult
         cspice_phaseq
         cspice_pjelpl (mice_pjelpl)
         cspice_pl2nvc
         cspice_pl2nvp
         cspice_pl2psv
         cspice_psv2pl
         cspice_pxfrm2
         cspice_spkcls
         cspice_spkcpo (mice_spkcpo)
         cspice_spkcpt (mice_spkcpt)
         cspice_spkcvo (mice_spkcvo)
         cspice_spkcvt (mice_spkcvt)
         cspice_spkopn
         cspice_spkpvn
         cspice_spksfs
         cspice_spkw08
         cspice_surfpt (mice_surfpt)
         cspice_timdef_get
         cspice_timdef_set
         cspice_vprjp
         cspice_vprjpi
         cspice_vproj
         cspice_xfmsta

   -Mice 1.2.0 28-APR-2010 (EDW)

      Added interfaces:

         cspice_bodc2s (mice_bodc2s)
         cspice_ducrss
         cspice_dvcrss
         cspice_dvdot
         cspice_dvhat
         cspice_dvnorm
         cspice_dvsep
         cspice_ekfind
         cspice_ekgc
         cspice_ekgd
         cspice_ekgi
         cspice_eknelt
         cspice_gfrr
         cspice_vperp
         cspice_vnorm
         cspice_unitim

   -Mice 1.1.0 23-FEB-2009 (EDW)

      Added interfaces:

         cspice_el2cgv
         cspice_gfdist
         cspice_gfposc
         cspice_gfrfov
         cspice_gfsep
         cspice_gfsntc
         cspice_gfsubc
         cspice_gftfov
         cspice_lmpool
         cspice_nvc2pl
         cspice_nvp2pl
         cspice_saelgv
         cspice_vrotv
         cspice_wncard
         cspice_wnsumd (mice_wnsumd)
         cspice_wnvald

   -Mice 1.0.0 12-FEB-2008 (EDW)

-Index_Entries

*/

static struct npf
   {
   char  * name;
   void    (*pfunc)(int, mxArray*[], int, const mxArray*[]);
   }
   NPF[] =
     {
     { "axisar_c", &cspice_axisar },
     { "azlcpo_c", &cspice_azlcpo },
     { "azlrec_c", &cspice_azlrec },
     { "b1900_c",  &cspice_b1900  },
     { "b1950_c",  &cspice_b1950  },
     { "badkpv_c", &cspice_badkpv },
     { "bltfrm_c", &cspice_bltfrm },
     { "bodc2n_s", &mice_bodc2n   },
     { "bodc2s_s", &mice_bodc2s   },
     { "boddef_c", &cspice_boddef },
     { "bodfnd_c", &cspice_bodfnd },
     { "bodn2c_s", &mice_bodn2c   },
     { "bods2c_s", &mice_bods2c   },
     { "bodvcd_c", &cspice_bodvcd },
     { "bodvrd_c", &cspice_bodvrd },
     { "ccifrm_s", &mice_ccifrm   },
     { "cgv2el_c", &cspice_cgv2el },
     { "chbder_c", &cspice_chbder },
     { "chbigr_c", &cspice_chbigr },
     { "chbint_c", &cspice_chbint },
     { "chbval_c", &cspice_chbval },
     { "cidfrm_s", &mice_cidfrm   },
     { "ckcls_c",  &cspice_ckcls  },
     { "ckcov_c",  &cspice_ckcov  },
     { "ckfrot_c", &cspice_ckfrot },
     { "ckfxfm_c", &cspice_ckfxfm },
     { "ckgp_c",   &cspice_ckgp   },
     { "ckgpav_c", &cspice_ckgpav },
     { "ckgr02_c", &cspice_ckgr02 },
     { "ckgr03_c", &cspice_ckgr03 },
     { "cklpf_c",  &cspice_cklpf  },
     { "ckmeta_c", &cspice_ckmeta },
     { "cknr02_c", &cspice_cknr02 },
     { "cknr03_c", &cspice_cknr03 },
     { "ckobj_c",  &cspice_ckobj  },
     { "ckopn_c",  &cspice_ckopn  },
     { "ckupf_c",  &cspice_ckupf  },
     { "ckw01_c",  &cspice_ckw01  },
     { "ckw02_c",  &cspice_ckw02  },
     { "ckw03_c",  &cspice_ckw03  },
     { "clight_c", &cspice_clight },
     { "clpool_c", &cspice_clpool },
     { "cnmfrm_s", &mice_cnmfrm   },
     { "conics_c", &cspice_conics },
     { "convrt_c", &cspice_convrt },
     { "cspice_mice", &cspice_mice},
     { "cyllat_c", &cspice_cyllat },
     { "cylrec_c", &cspice_cylrec },
     { "cylsph_c", &cspice_cylsph },
     { "dafac_c",  &cspice_dafac  },
     { "dafbbs_c", &cspice_dafbbs },
     { "dafbfs_c", &cspice_dafbfs },
     { "dafcls_c", &cspice_dafcls },
     { "dafcs_c",  &cspice_dafcs  },
     { "dafdc_c",  &cspice_dafdc  },
     { "dafec_c",  &cspice_dafec  },
     { "daffna_c", &cspice_daffna },
     { "daffpa_c", &cspice_daffpa },
     { "dafgda_c", &cspice_dafgda },
     { "dafgn_c",  &cspice_dafgn  },
     { "dafgs_c",  &cspice_dafgs  },
     { "dafhsf_c", &cspice_dafhsf },
     { "dafopr_c", &cspice_dafopr },
     { "dafopw_c", &cspice_dafopw },
     { "dafps_c",  &cspice_dafps  },
     { "dafrs_c",  &cspice_dafrs  },
     { "dafus_c",  &cspice_dafus  },
     { "dasadc_c", &cspice_dasadc },
     { "dasadd_c", &cspice_dasadd },
     { "dasadi_c", &cspice_dasadi },
     { "dascls_c", &cspice_dascls },
     { "dasec_c",  &cspice_dasec  },
     { "dashfs_c", &cspice_dashfs },
     { "daslla_c", &cspice_daslla },
     { "dasllc_c", &cspice_dasllc },
     { "dasonw_c", &cspice_dasonw },
     { "dasopr_c", &cspice_dasopr },
     { "dasops_c", &cspice_dasops },
     { "dasopw_c", &cspice_dasopw },
     { "dasrdc_c", &cspice_dasrdc },
     { "dasrdd_c", &cspice_dasrdd },
     { "dasrdi_c", &cspice_dasrdi },
     { "dasudc_c", &cspice_dasudc },
     { "dasudd_c", &cspice_dasudd },
     { "dasudi_c", &cspice_dasudi },
     { "daswbr_c", &cspice_daswbr },
     { "dazldr_c", &cspice_dazldr },
     { "dcyldr_c", &cspice_dcyldr },
     { "deltet_c", &cspice_deltet },
     { "dgeodr_c", &cspice_dgeodr },
     { "dlabbs_c", &cspice_dlabbs },
     { "dlabfs_c", &cspice_dlabfs },
     { "dlabns_c", &cspice_dlabns },
     { "dlaens_c", &cspice_dlaens },
     { "dlafns_c", &cspice_dlafns },
     { "dlafps_c", &cspice_dlafps },
     { "dlaopn_c", &cspice_dlaopn },
     { "dlatdr_c", &cspice_dlatdr },
     { "dnearp_c", &cspice_dnearp },
     { "dpgrdr_c", &cspice_dpgrdr },
     { "dpmax_c",  &cspice_dpmax  },
     { "dpmin_c",  &cspice_dpmin  },
     { "dpr_c",    &cspice_dpr    },
     { "drdazl_c", &cspice_drdazl },
     { "drdcyl_c", &cspice_drdcyl },
     { "drdgeo_c", &cspice_drdgeo },
     { "drdlat_c", &cspice_drdlat },
     { "drdpgr_c", &cspice_drdpgr },
     { "drdsph_c", &cspice_drdsph },
     { "dskb02_c", &cspice_dskb02 },
     { "dskcls_c", &cspice_dskcls },
     { "dskd02_c", &cspice_dskd02 },
     { "dskgd_c",  &cspice_dskgd  },
     { "dskgtl_c", &cspice_dskgtl },
     { "dski02_c", &cspice_dski02 },
     { "dskmi2_c", &cspice_dskmi2 },
     { "dskn02_c", &cspice_dskn02 },
     { "dskobj_c", &cspice_dskobj },
     { "dskopn_c", &cspice_dskopn },
     { "dskp02_c", &cspice_dskp02 },
     { "dskrb2_c", &cspice_dskrb2 },
     { "dsksrf_c", &cspice_dsksrf },
     { "dskstl_c", &cspice_dskstl },
     { "dskv02_c", &cspice_dskv02 },
     { "dskw02_c", &cspice_dskw02 },
     { "dskx02_c", &cspice_dskx02 },
     { "dskxsi_c", &cspice_dskxsi },
     { "dskxv_c",  &cspice_dskxv  },
     { "dskz02_c", &cspice_dskz02 },
     { "dsphdr_c", &cspice_dsphdr },
     { "dtpool_s", &mice_dtpool   },
     { "ducrss_c", &cspice_ducrss },
     { "dvcrss_c", &cspice_dvcrss },
     { "dvdot_c",  &cspice_dvdot  },
     { "dvhat_c",  &cspice_dvhat  },
     { "dvnorm_c", &cspice_dvnorm },
     { "dvpool_c", &cspice_dvpool },
     { "dvsep_c",  &cspice_dvsep  },
     { "edlimb_s", &mice_edlimb   },
     { "ednmpt_c", &cspice_ednmpt },
     { "edpnt_c",  &cspice_edpnt  },
     { "edterm_c", &cspice_edterm },
     { "ekfind_c", &cspice_ekfind },
     { "ekgc_c",   &cspice_ekgc   },
     { "ekgd_c",   &cspice_ekgd   },
     { "ekgi_c",   &cspice_ekgi   },
     { "eknelt_c", &cspice_eknelt },
     { "el2cgv_c", &cspice_el2cgv },
     { "et2lst_c", &cspice_et2lst },
     { "et2utc_c", &cspice_et2utc },
     { "etcal_c",  &cspice_etcal  },
     { "eul2m_c",  &cspice_eul2m  },
     { "eul2xf_c", &cspice_eul2xf },
     { "evsgp4_c", &cspice_evsgp4 },
     { "expool_c", &cspice_expool },
     { "fovray_c", &cspice_fovray },
     { "fovtrg_c", &cspice_fovtrg },
     { "frame_c",  &cspice_frame  },
     { "frinfo_s", &mice_frinfo   },
     { "frmnam_c", &cspice_frmnam },
     { "furnsh_c", &cspice_furnsh },
     { "gcpool_c", &cspice_gcpool },
     { "gdpool_c", &cspice_gdpool },
     { "georec_c", &cspice_georec },
     { "getelm_c", &cspice_getelm },
     { "getfat_c", &cspice_getfat },
     { "getfov_c", &cspice_getfov },
     { "getfvn_c", &cspice_getfvn },
     { "gfdist_c", &cspice_gfdist },
     { "gfilum_c", &cspice_gfilum },
     { "gfoclt_c", &cspice_gfoclt },
     { "gfpa_c",   &cspice_gfpa   },
     { "gfposc_c", &cspice_gfposc },
     { "gfrfov_c", &cspice_gfrfov },
     { "gfrr_c",   &cspice_gfrr   },
     { "gfsep_c",  &cspice_gfsep  },
     { "gfsntc_c", &cspice_gfsntc },
     { "gfstol_c", &cspice_gfstol },
     { "gfsubc_c", &cspice_gfsubc },
     { "gftfov_c", &cspice_gftfov },
     { "gipool_c", &cspice_gipool },
     { "gnpool_c", &cspice_gnpool },
     { "halfpi_c", &cspice_halfpi },
     { "hrmesp_c", &cspice_hrmesp },
     { "hrmint_c", &cspice_hrmint },
     { "illum_c",  &cspice_illum  },
     { "illum_pl02",        &cspice_illum_pl02    },
     { "illum_plid_pl02_s", &mice_illum_plid_pl02 },
     { "illumf_s", &mice_illumf   },
     { "illumg_s", &mice_illumg   },
     { "ilumin_s", &mice_ilumin   },
     { "inedpl_c", &cspice_inedpl },
     { "inelpl_c", &cspice_inelpl },
     { "inrypl_c", &cspice_inrypl },
     { "intmax_c", &cspice_intmax },
     { "intmin_c", &cspice_intmin },
     { "invort_c", &cspice_invort },
     { "invstm_c", &cspice_invstm },
     { "j1900_c",  &cspice_j1900  },
     { "j1950_c",  &cspice_j1950  },
     { "j2000_c",  &cspice_j2000  },
     { "j2100_c",  &cspice_j2100  },
     { "jyear_c",  &cspice_jyear  },
     { "kclear_c", &cspice_kclear },
     { "kdata_c",  &cspice_kdata  },
     { "kinfo_c",  &cspice_kinfo  },
     { "kplfrm_c", &cspice_kplfrm },
     { "ktotal_c", &cspice_ktotal },
     { "latcyl_c", &cspice_latcyl },
     { "latrec_c", &cspice_latrec },
     { "latsph_c", &cspice_latsph },
     { "latsrf_c", &cspice_latsrf },
     { "ldpool_c", &cspice_ldpool },
     { "lgresp_c", &cspice_lgresp },
     { "lgrind_c", &cspice_lgrind },
     { "lgrint_c", &cspice_lgrint },
     { "limbpt_c", &cspice_limbpt },
     { "limb_pl02",   &cspice_limb_pl02   },
     { "llgrid_pl02", &cspice_llgrid_pl02 },
     { "lmpool_c", &cspice_lmpool },
     { "lspcn_c",  &cspice_lspcn  },
     { "ltime_c",  &cspice_ltime  },
     { "m2eul_c",  &cspice_m2eul  },
     { "m2q_c",    &cspice_m2q    },
     { "namfrm_c", &cspice_namfrm },
     { "nearpt_s", &mice_nearpt   },
     { "nextwd_c", &cspice_nextwd },
     { "npedln_s", &mice_npedln   },
     { "npelpt_s", &mice_npelpt   },
     { "nplnpt_s", &mice_nplnpt   },
     { "nthwd_c",  &cspice_nthwd  },
     { "nvc2pl_c", &cspice_nvc2pl },
     { "nvp2pl_c", &cspice_nvp2pl },
     { "occult_c", &cspice_occult },
     { "oscelt_c", &cspice_oscelt },
     { "oscltx_c", &cspice_oscltx },
     { "pckcov_c", &cspice_pckcov },
     { "pckfrm_c", &cspice_pckfrm },
     { "pcpool_c", &cspice_pcpool },
     { "pdpool_c", &cspice_pdpool },
     { "pgrrec_c", &cspice_pgrrec },
     { "phaseq_c", &cspice_phaseq },
     { "pi_c",     &cspice_pi     },
     { "pipool_c", &cspice_pipool },
     { "pjelpl_s", &mice_pjelpl   },
     { "pl2nvc_c", &cspice_pl2nvc },
     { "pl2nvp_c", &cspice_pl2nvp },
     { "pl2psv_c", &cspice_pl2psv },
     { "pltar_c",  &cspice_pltar  },
     { "pltexp_c", &cspice_pltexp },
     { "pltnp_c",  &cspice_pltnp  },
     { "pltnrm_c", &cspice_pltnrm },
     { "pltvol_c", &cspice_pltvol },
     { "polyds_c", &cspice_polyds },
     { "prop2b_c", &cspice_prop2b },
     { "psv2pl_c", &cspice_psv2pl },
     { "pxform_c", &cspice_pxform },
     { "pxfrm2_c", &cspice_pxfrm2 },
     { "q2m_c",    &cspice_q2m    },
     { "qderiv_c", &cspice_qderiv },
     { "qdq2av_c", &cspice_qdq2av },
     { "qxq_c",    &cspice_qxq    },
     { "radrec_c", &cspice_radrec },
     { "rav2xf_c", &cspice_rav2xf },
     { "raxisa_c", &cspice_raxisa },
     { "recazl_c", &cspice_recazl },
     { "reccyl_c", &cspice_reccyl },
     { "recgeo_c", &cspice_recgeo },
     { "reclat_c", &cspice_reclat },
     { "recpgr_c", &cspice_recpgr },
     { "recrad_c", &cspice_recrad },
     { "recsph_c", &cspice_recsph },
     { "repmc_c",  &cspice_repmc  },
     { "repmct_c", &cspice_repmct },
     { "repmd_c",  &cspice_repmd  },
     { "repmf_c",  &cspice_repmf  },
     { "repmi_c",  &cspice_repmi  },
     { "repml_c",  &cspice_repml  },
     { "repmot_c", &cspice_repmot },
     { "rotate_c", &cspice_rotate },
     { "rotmat_c", &cspice_rotmat },
     { "rotvec_c", &cspice_rotvec },
     { "rpd_c",    &cspice_rpd    },
     { "saelgv_c", &cspice_saelgv },
     { "scdecd_c", &cspice_scdecd },
     { "sce2c_c",  &cspice_sce2c  },
     { "sce2s_c",  &cspice_sce2s  },
     { "scencd_c", &cspice_scencd },
     { "scfmt_c",  &cspice_scfmt  },
     { "scpart_c", &cspice_scpart },
     { "scs2e_c",  &cspice_scs2e  },
     { "sct2e_c",  &cspice_sct2e  },
     { "sctiks_c", &cspice_sctiks },
     { "sincpt_s", &mice_sincpt   },
     { "spd_c",    &cspice_spd    },
     { "sphcyl_c", &cspice_sphcyl },
     { "sphlat_c", &cspice_sphlat },
     { "sphrec_c", &cspice_sphrec },
     { "spkapo_s", &mice_spkapo   },
     { "spkcls_c", &cspice_spkcls },
     { "spkcov_c", &cspice_spkcov },
     { "spkcpo_s", &mice_spkcpo   },
     { "spkcpt_s", &mice_spkcpt   },
     { "spkcvo_s", &mice_spkcvo   },
     { "spkcvt_s", &mice_spkcvt   },
     { "spkez_s",  &mice_spkez    },
     { "spkezr_s", &mice_spkezr   },
     { "spkgeo_s", &mice_spkgeo   },
     { "spklef_c", &cspice_spklef },
     { "spkobj_c", &cspice_spkobj },
     { "spkopn_c", &cspice_spkopn },
     { "spkpos_s", &mice_spkpos   },
     { "spkpvn_s", &mice_spkpvn   },
     { "spksfs_s", &mice_spksfs   },
     { "spkssb_c", &cspice_spkssb },
     { "spkuef_c", &cspice_spkuef },
     { "spkw08_c", &cspice_spkw08 },
     { "spkw09_c", &cspice_spkw09 },
     { "spkw10_c", &cspice_spkw10 },
     { "spkw13_c", &cspice_spkw13 },
     { "srfc2s_s", &mice_srfc2s   },
     { "srfcss_s", &mice_srfcss   },
     { "srfnrm_c", &cspice_srfnrm },
     { "srfrec_c", &cspice_srfrec },
     { "srfs2c_s", &mice_srfs2c   },
     { "srfscc_s", &mice_srfscc   },
     { "srfxpt_s", &mice_srfxpt   },
     { "stelab_c", &cspice_stelab },
     { "stlabx_c", &cspice_stlabx },
     { "stpool_c", &cspice_stpool },
     { "str2et_c", &cspice_str2et },
     { "subpnt_s", &mice_subpnt   },
     { "subpt_pl02", &cspice_subpt_pl02  },
     { "subpt_s",  &mice_subpt    },
     { "subslr_s", &mice_subslr   },
     { "subsol_c", &cspice_subsol },
     { "subsol_pl02", &cspice_subsol_pl02 },
     { "surfnm_c", &cspice_surfnm },
     { "surfpt_s", &mice_surfpt   },
     { "surfpv_c", &cspice_surfpv },
     { "sxform_c", &cspice_sxform },
     { "szpool_c", &cspice_szpool },
     { "tangpt_c", &cspice_tangpt },
     { "term_pl02", &cspice_term_pl02 },
     { "termpt_c", &cspice_termpt },
     { "timdef_get_c", &cspice_timdef_get  },
     { "timdef_set_c", &cspice_timdef_set  },
     { "timout_c", &cspice_timout },
     { "tipbod_c", &cspice_tipbod },
     { "tisbod_c", &cspice_tisbod },
     { "tkfram_c", &cspice_tkfram },
     { "tkvrsn_c", &cspice_tkvrsn },
     { "tparch_c", &cspice_tparch },
     { "tparse_c", &cspice_tparse },
     { "tpictr_c", &cspice_tpictr },
     { "trgsep_c", &cspice_trgsep },
     { "tsetyr_c", &cspice_tsetyr },
     { "twopi_c",  &cspice_twopi  },
     { "twovec_c", &cspice_twovec },
     { "twovxf_c", &cspice_twovxf },
     { "tyear_c",  &cspice_tyear  },
     { "unitim_c", &cspice_unitim },
     { "unload_c", &cspice_unload },
     { "unorm_c",  &cspice_unorm  },
     { "vdist_c",  &cspice_vdist  },
     { "vhat_c",   &cspice_vhat   },
     { "vnorm_c",  &cspice_vnorm  },
     { "vperp_c",  &cspice_vperp  },
     { "vprjp_c",  &cspice_vprjp  },
     { "vprjpi_c", &cspice_vprjpi },
     { "vproj_c",  &cspice_vproj  },
     { "vprojg_c", &cspice_vprojg },
     { "vrotv_c",  &cspice_vrotv  },
     { "vsep_c",   &cspice_vsep   },
     { "vupack_c", &cspice_vupack },
     { "wncard_c", &cspice_wncard },
     { "wncomd_c", &cspice_wncomd },
     { "wndifd_c", &cspice_wndifd },
     { "wnelmd_c", &cspice_wnelmd },
     { "wnexpd_c", &cspice_wnexpd },
     { "wnextd_c", &cspice_wnextd },
     { "wnfetd_c", &cspice_wnfetd },
     { "wnfild_c", &cspice_wnfild },
     { "wnfltd_c", &cspice_wnfltd },
     { "wnincd_c", &cspice_wnincd },
     { "wninsd_c", &cspice_wninsd },
     { "wnintd_c", &cspice_wnintd },
     { "wnreld_c", &cspice_wnreld },
     { "wnsumd_s", &mice_wnsumd   },
     { "wnunid_c", &cspice_wnunid },
     { "wnvald_c", &cspice_wnvald },
     { "xf2eul_c", &cspice_xf2eul },
     { "xf2rav_c", &cspice_xf2rav },
     { "xfmsta_c", &cspice_xfmsta },
     };
#endif


