
DSFN = DataSet FileNaming

Package for generating and interpreting official SolO dataset filenames. This
includes datasets from instruments other than SolO/RPW.

Supports unofficial filenaming used internally by SolO/RPW during early phases
of development of RPW.
Supports unofficial extension to official dataset filenaming. (This support will
probably be abolished eventually.)

NOTE: There is one esoteric official non-RPW case of dataset filenames which is
      not (yet) supported:
      solo.adm.dsfn.create_dataset_filename() will not (yet) work for dataset
      filenames with uppercase outside of the archiving level in the dataset ID.
      Ex: solo_L1_swa-eas2-NM3D_20201027T000007-20201027T030817_V01.cdf
NOTE: Code is used by e.g. BICAS, and BICAS's batch processing
      (bicas.tools.batch).
NOTE: This code is managed by Erik P G Johansson. There may be
      backwards-incompatible change made to it (though rarely)!


================================
SUPPORTED FILENAMING CONVENTIONS
================================

--------------------------------------------------------------------
"LES" filenaming convention (datasets from early RPW ground testing)
--------------------------------------------------------------------
solo_L1_rpw-tds-lfm-rswf_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
--
NOTE: Always time interval on form SECOND_TO_SECOND.
NOTE: Found in e.g. brain.irfu.se:/data/solo/necp_em1_run/.


--------------------------------------------------------------------
"CNE" filenaming convention (datasets from early RPW ground testing)
--------------------------------------------------------------------
ROC-SGSE_HK_RPW-BIA_19850de_CNE_V02.cdf
--
NOTE: Uppercase.
NOTE: No time interval string.
NOTE: "ROC-SGSE" instead of "SOLO".
NOTE: Found in e.g. brain.irfu.se:/data/solo/remote/roc_sgse/.


------------------------------
Inflight filenaming convention
------------------------------
solo_HK_rpw-bia_20200301_V01.cdf                   # NOTE: No -cdag.
solo_L2_rpw-lfr-surv-cwf-e-cdag_20200213_V01.cdf   # NOTE: -cdag.
solo_L1_rpw-bia-sweep-cdag_20200307T053018-20200307T053330_V01.cdf
solo_L1_rpw-bia-current-cdag_20200401T000000-20200421T000000_V01.cdf
    NOTE: Phased out for currents (not sweeps).
          /EJ+XB-mail 2020-05-27
solo_L1_rpw-bia-current-cdag_20200301-20200331_V01.cdf
--
NOTE: With and without "cdag". "cdag" is an RPW-internal extension to
      the official SolO filenaming convention.
NOTE: Multiple ways of specifying time interval.


============================================
UNSUPPORTED FILENAMING CONVENTIONS: EXAMPLES
============================================
Master CDFs:
  SOLO_HK_RPW-BIA_V01.cdf
  SOLO_L2_RPW-LFR-SURV-CWF-E_V04.cdf
Draft Master CDFs:
  SOLO_L1_RPW-BIA-CURRENT_V01.Draft.cdf
      (Probably not even official; probably ROC ad hoc filenames)
Summary plots (which BTW are NOT datasets!):
  solo_L3_rpw-lfr-surv-cwf-e_20200423_V01.png
  solo_L3_rpw-lfr-surv-swf-e_20200423_V01.png


/Erik P G Johansson 2024-08-13
