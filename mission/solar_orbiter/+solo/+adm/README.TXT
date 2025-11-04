
ADM = Automatic Dataset Management

Package for automatic "management" of SolO datasets which follow official file
naming conventions. E.g. for
- identifying datasets among non-dataset files (based on filenaming)
- iterating over datasets
- finding datasets which overlap with a selected time interval
- identifying datasets with only the latest version
- identifying datasets of specific types (dataset ID)
- grouping datasets
etc.

NOTE: The code is used by e.g. BICAS proper and BICAS's batch processing tool
      (bicas.tools.batch).
NOTE: This code is developed and managed by Erik P G Johansson.
      Backwards-incompatible changes may be made to it!

/Erik P G Johansson 2025-11-04


#########################
Miscellaneous conventions
#########################
DF
      Class solo.adm.dsfn.DatasetFilename
