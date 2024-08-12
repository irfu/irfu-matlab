
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

/Erik P G Johansson 2024-08-13
