
package for generating SolO quicklooks available at
https://www.space.irfu.se/soloql/.

QLI = QuickLooks IRFU ("IRFU quicklooks")

NOTE: Not to be confused with LESIA's summary plots or BIA_QL3.
NOTE: These quicklooks use both SolO/RPW datasets AND datasets from
other instruments.



=================================
 How to generate quicklook files
=================================
For an example of how to generate multiple quicklooks (files) locally, see
solo.qli.batch.generate_quicklooks_demo().



===========================================
 How to edit the content of the quicklooks
===========================================
generate_quicklooks_24h_6h_2h.m
generate_quicklook_7days.m
    Files which do the actual plotting.

generate_quicklooks_24h_6h_2h_using_DB_SPICE.m
generate_quicklook_7days_using_DB_SPICE.m
    Files which call generate_quicklooks_24h_6h_2h.m
    and generate_quicklook_7days.m
    respectively and supply them with bulk data. If you want the quicklooks to
    have access to more data e.g. more science data (zVariables), then you need
    to edit these.



==========================================================
 Various source code and documentation naming conventions
==========================================================
NOTE: The source code additionally tries to follow BICAS naming conventions. See
      irfu-matlab/mission/solar_orbiter/bicas/readme.txt
--
DASA
    Date Array Selection Algorithm. One of any of a set of algorithms for
    selecting which dates for which to generate quicklooks.
    Ex: DMRQ
DASAID
    DASA ID. String constant which represents a specific DASA.
DICT
    MATLAB dictionary object.
DMRQ
    "Datasets More Recent than QLI (file)" DASA.
    A day is included if (1) the most recent QLI FMD (for that day) is
    more recent than the most recent input dataset (for that day), or (2) there
    is no QLI for that day.
DSMD
    Class solo.adm.DatasetMetadata.
DT
    datetime (class).
FMD
    File Modification Date (timestamp).
    If represented as datetime, then always without TimeZone.
FSR
    Class solo.qli.batch.FileSystemReaderAbstract and subclasses thereof.
GQL
    Class solo.qli.batch.GenerateQuicklooksAbstract and subclasses thereof.
QLI
    QuickLooks IRFU ("IRFU quicklooks"). In the code it refers to the file
    representations of the quicklooks.
SolO DB
    Informal term for the irfu-matlab functionality behind solo.get_ts() and
    solo.db_list_files().
UFD
    Class solo.qli.batch.UmdFmdDictionary. Models a dictionary UMD-->FMD.
UMD
    UTC Midnight Date (timestamp). Timestamp in UTC at midnight. Used for
    uniquely representing an entire day as a whole.
