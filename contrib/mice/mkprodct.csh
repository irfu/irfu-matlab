#!/bin/csh
#
#   Linux version.
#
#   This script builds the shared object library for the Matlab-CSPICE
#   interface (Mice). It assumes that it is executed from one of the
#   "product" directories in a tree that looks like the one displayed
#   below:
#
#                      package
#                         |
#                         |
#       +------+------+------+------+------+
#       |      |      |      |      |      |
#     data    doc    etc    exe    lib    src
#                                          |
#                                          |
#                         +----------+----------+------- ... ------+
#                         |          |          |                  |
#                     product_1  product_2  product_3    ...   product_n
#
#   Here's the basic strategy:
#
#     1)  Compile all of the .c files in the current directory to
#         LOCALLIB. The "mex" script include with the Matlab
#         distribution initiates the compile using the options defined
#         in the local mexopts.sh file. The "mex" is assumed to in
#         one of the directories in the current executable path.
#
#     2)  Move LOCALLIB to the 'lib' directory, set the library to the
#         correct name.
#
#   Change History:
#   ===============
#
#   Version 1.2.0,  Jan. 13, 2020, Boris Semenov
#
#      Added the -Wno-unused-but-set-variable option to suppress
#      warnings for set-but-unused "extra" returned by mice_checkargs
#      in non-vectorized wrappers.
#
#   Version 1.1.0,  Sep. 23, 2016, Ed Wright
#
#      Modified mex call to compensate for changes in MEX functionality.
#
#   Version 1.0.0,  Feb. 14, 2008, Ed Wright & Boris Semenov
#

#
# Check if "mex" is in the current executable path. If not, complain
# and stop.
#
which mex >& /dev/null

if ( $status != 0 ) then
   echo " "
   echo "      The 'mex' script include with the Matlab is not in the "
   echo "      current path:"
   echo " "
   echo $PATH
   echo " "
   echo "      As a result, I am unable to build this product."
   echo " "
   exit 1;
endif

#
# Set the names of the local and final libraries.
#
set LOCALLIB = "mice.mexa64"
set LIBRARY  = "../../lib/$LOCALLIB"

#
# Set compile options.
#
set MEXOPTIONS = '-Dunix -compatibleArrayDims'

#
# Compile all .c modules in the current directory into the local
# library.
#
mex $MEXOPTIONS *.c -output $LOCALLIB -I../../include       \
    CLIBS='$CLIBS ../../lib/cspice.a'                       \
    CFLAGS='$CFLAGS -ansi -Wall -Wundef -Wpointer-arith -Wcast-align -Wsign-compare -Wno-unused-but-set-variable'


#
# Move local library to its final location.
#
\mv $LOCALLIB $LIBRARY

#
# All done.
#
exit 0
