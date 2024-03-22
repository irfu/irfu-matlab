
NOTE: This code should NOT be of interest to the typical user.

offgen = Official Generation of quicklooks

This MATLAB package (directory) contains code which is associated
with the particular setup on IRFU servers which is used for generating official
IRFU quicklooks at IRFU (/data/solo/quicklooks/www/).

This code can be seen as a "user of" (or a particular "application of")
solo.qli and other code, but does not contribute any functionality/code that
is of interest outside of generating offical IRFU quicklooks on IRFU servers.
The actual code is intended to only contain thin wrappers which add very little
additional functionality to the outside-package code that is called.
It should specify particular hardcoded paths and particular server names when
needed.

Some of these functions are explicitly designed for being called from the OS,
e.g. bash scripts, and hence only accept string-valued arguments. Such
functions have "bash" in the name.
