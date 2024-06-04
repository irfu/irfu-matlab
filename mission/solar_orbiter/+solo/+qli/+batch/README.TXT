
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
