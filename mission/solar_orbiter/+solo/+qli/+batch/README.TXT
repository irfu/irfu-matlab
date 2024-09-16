
Package for *batch* generation of SolO quicklooks (in solo.qli).

The code is designed to generic in the sense that it can be reused for
any computer system. In practice, it is primarily intended for official
generation of quicklooks on a particular setup on a particular server at IRFU.
It is intended for both cron jobs and manually commanded generation of
quicklooks.

See irfu-matlab/mission/solar_orbiter/+solo/+qli/README.TXT for naming
conventions.


===================================================
How to do batch generation of quicklooks in general
===================================================
See solo.qli.batch.generate_quicklooks_demo().


=================================================
How to do OFFICIAL batch generation of quicklooks
=================================================
See solo.qli.batch.generate_quicklooks_shell().



solo.qli.batch.generate_quicklooks_shell()
    Intended to be called from the OS, e.g. bash scripts.
    Reads configuration file.
    Initializes "SolO DB".
    Calls solo.qli.batch.generate_quicklooks_syntax() to do most work.

solo.qli.batch.generate_quicklooks_syntax()
    (1) Obtains list of days to generate QLIs by calling
        solo.qli.batch.interface.get_days_from_DASA().
    (2) (Nominally) generates QLIs by calling
        solo.qli.batch.generate_quicklooks()

solo.qli.batch.generate_quicklooks()
    (1) (Nominally) iterates over submitted list dates to make 24h/6h/2h QLIs.
        Calls
        solo.qli.batch.GenerateQuicklooksImplementation.generate_quicklooks_24h_6h_2h_using_DB_SPICE().
    (2) (Nominally) constructs list of weeks and iterates over it to make weekly
        QLIs. Calls
        solo.qli.batch.GenerateQuicklooksImplementation.generate_quicklook_7days_using_DB_SPICE().

