%***************************************************************************************************
% Copyright 2006, T.P. O'Brien
%
% This file is part of IRBEM-LIB.
%
%    IRBEM-LIB is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    IRBEM-LIB is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
%
INTRODUCTION

We provide a set of Matlab wrappers for the IRBEM FORTRAN library
through the use of a shared library file (.dll on windows, .so on
Unix). Matlab refers to this type of interface as a "mex file". In
addition to the shared library file, Matlab requires a C-style header
file that includes the C-style calling syntax for each FORTRAN
function. This .h file is provided with the IRBEM Matlab library.

INSTALLATION

The IRBEM compile scripts will create the shared library file for
use with IDL or Matlab. A copy of the shared library is created in the
"matlab" subdirectory, and is given a name "libirbem.dll" or
"libirbem.so" depending on the platform. Thus, in order to
access the library and the wrappers, the user needs only add the
"matlab" subdirectory (e.g. "c:\onera_desp\matlab") to the Matlab
path. In order for the library to easily locate the data files that it
needs, it is also advisable to add the "data" subdirectory to the
Matlab path. Matlab provides a GUI interface to add folders to the
search path from "Set Path" entry on the "File" menu. It is also
possible to add these two folders to the path using the Matlab
"addpath" command; however, this command does not make the change to
the path permanent, so that "addpath" will have to be issued once in
each new Matlab session. The "savepath" function can make the path
change permanent, as can the GUI path manager interface. Thus, one
could use the following commands, if the IRBEM library is located in
"c:\onera_desp":

addpath c:\onera_desp\matlab c:\onera_desp\data
savepath

USE OF THE MATLAB LIBRARY

The Matlab wrappers are built to provide Matlab-like function calls
into the library. These wrappers handle loading the library into
Matlab's function space, formatting inputs and outputs for proper
calls to the library, and providing "vectorized" functionality
wherever reasonable. Each wrapper function will determine whether the
FORTRAN library has been loaded, and, if not, attempt to load it from
the default location (anywhere in the Matlab path).

Each .m file provided includes a robust set of helps obtainable via
help <funtion_name> in the usual Matlab way. This help call will
provide details of how to call the function. 

Although the FORTRAN library often limits the size of arrays, the
Matlab wrappers typically handle arbitrarily large array inputs by
splicing together multiple calls to the library. Also, in many cases,
when a set of arrays is expected as input, the Matlab wrappers will
accept scalars for some, which will be repeated (via Matlab's repmat
function) to be the same size as the other array arguments.

The Matlab library also provides a single function for doing all
coordinate rotations. This function, onera_desp_lib_rotate will
perform rotations between any coordinate system supported by the
library. It will also compute the most efficient multi-step rotation
for those rotations not explicitly implemented in the FORTRAN library.

In many cases, the library requires integer inputs that represent
different options, e.g., kext=5 for the Olson-Pfitzer Quiet external
field model. In most cases, the Matlab wrapper supports string
(keyword) inputs in place of the integer values. This keyword approach
is implemented for kext, options, sysaxes, and whichm, among others.

Whenever date/time arguments are required by the FORTRAN library, the
Matlab library expects Matlab Date Numbers (construct argument
"matlabd" with Matlab's datenum function).

The FORTRAN library function fly_in_afrl_crres requires a set of text
files. If the path to these files is not specified, the wrapper will
attempt to guess it by locating one of these, 'crrespro_quiet.txt' in
the Matlab search path.

WRAPPER FUNCTIONS:

[Lm,Lstar,Blocal,Bmin,J,POSIT] = onera_desp_lib_drift_shell(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)

[Bmin,xGEO] = onera_desp_lib_find_magequator(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)

[Blocal,Bmirror,xGEO] = onera_desp_lib_find_mirror_point(kext,options,sysaxes,matlabd,x1,x2,x3,alpha,maginput)

Flux = onera_desp_lib_fly_in_afrl_crres(sysaxes,whichm,energy,x1,x2,x3,Ap15,crres_path)

Flux = onera_desp_lib_fly_in_nasa_aeap(sysaxes,whichm,energy,x1,x2,x3)

[Bgeo,B] = onera_desp_lib_get_field(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)

MLT = onera_desp_lib_get_mlt(matlabd,xGEO)

[Lm,Lstar,Blocal,Bmin,J,MLT] = onera_desp_lib_make_lstar(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)

[Lm,Lstar,Blocal,Bmin,J,MLT] = onera_desp_lib_make_lstar_shell_splitting(kext,options,sysaxes,matlabd,x1,x2,x3,alpha,maginput)

[Lm,Blocal,Bmin,J,POSIT] = onera_desp_lib_trace_field_line(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)

The following function handles all coordinate rotations. In some
cases, the FORTRAN library requests the coordinates as separate
arguments (e.g., alt, lat, lon). For these functions, the Matlab
wrapper expects a single Nx3 matrix, constructed by, e.g., X =
[alt(:), lat(:), lon(:)]. Also, some rotations are time invariant and
the matlabd argument can be omitted. The helps for
onera_desp_lib_rotate provide a list of available rotations. This
function will also perform multi-step rotations as needed to move
between any set of supported coordinate systems.

Y = onera_desp_lib_rotate(X,rotation,matlabd)
[Y,psi] = onera_desp_lib_rotate(X,rotation,matlabd)

HELPER FUNCTIONS

Helper functions that support the wrappers:

The following function creates the proper maginputs array for use with
the field models.

maginputs = onera_desp_lib_maginputs(Kp,Dst,Nsw,Vsw,Psw,ByGSM,BzGSM,G1,G2,G3,W1,W2,W3,W4,W5,W6,AL);

The following funnction will load the shared library file from the
non-default location.

onera_desp_lib_load(libfile,headerfile);

The following functions are used by the library to convert the Matlab
wrapper inputs into the inputs needed by the FORTRAN library. These
include as appropriate the look-up tables that convert keyword inputs
into the integer constants used by the FORTAN library.

kext = onera_desp_lib_kext(kext);
options  = onera_desp_lib_options(inoptions)
sysaxes = onera_desp_lib_sysaxes(sysaxes)
[iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd)


COORDINATE ROTATIONS
GEO2GSM
xGSM = onera_desp_lib_rotate(xGEO,'geo2gsm',matlabd);
[xGSM,psi] = onera_desp_lib_rotate(xGEO,'geo2gsm',matlabd);

GSM2GEO
xGEO = onera_desp_lib_rotate(xGSM,'gsm2geo',matlabd);
[xGEO,psi] = onera_desp_lib_rotate(xGSM,'gsm2geo',matlabd);

GEO2GSE
xGSE = onera_desp_lib_rotate(xGEO,'geo2gse',matlabd);

GSE2GEO
xGEO = onera_desp_lib_rotate(xGSE,'gse2geo',matlabd);

GDZ2GEO
******NOTE ORDER OF GDZ COORDINATES********
xGEO = onera_desp_lib_rotate([alti(:), lati(:), longi(:)],'gdz2geo');
xx = xGEO(:,1); yy = xGEO(:,2); zz = xGEO(:,3);

GEO2GDZ
******NOTE ORDER OF GDZ COORDINATES********
xGDZ = onera_desp_lib_rotate([xx(:) yy(:) zz(:)],'geo2gdz');
alti = xGDZ(:,1); lati = xGDZ(:,2); longi = xGDZ(:,3); 

GEO2GEI
xGEI = onera_desp_lib_rotate(xGEO,'geo2gei',matlabd);

GEI2GEO
xGEO = onera_desp_lib_rotate(xGEI,'gei2geo',matlabd);

GEO2SM
xSM = onera_desp_lib_rotate(xGEO,'geo2sm',matlabd);

SM2GEO
xGEO = onera_desp_lib_rotate(xSM,'sm2geo',matlabd);

GSM2SM
xSM = onera_desp_lib_rotate(xGSM,'gsm2sm',matlabd);

SM2GSM
xGSM = onera_desp_lib_rotate(xSM,'sm2gsm',matlabd);

GEO2MAG
xMAG = onera_desp_lib_rotate(xGEO,'geo2mag',matlabd);

MAG2GEO
xGEO = onera_desp_lib_rotate(xMAG,'mag2geo',matlabd);

SPH2CAR
xCAR = onera_desp_lib_rotate([r(:), lati(:), longi(:)],'sph2car');

CAR2SPH
xSPH = onera_desp_lib_rotate(xCAR,'car2sph');
r = xSPH(:,1); lati = xSPH(:,2); longi = xSPH(:,3);

RLL2GDZ
******NOTE ORDER OF GDZ COORDINATES********
xGDZ = onera_desp_lib_rotate([r(:), lati(:), longi(:)],'rll2gdz');
alti = xGDZ(:,1); lati = xGDZ(:,2); longi = xGDZ(:,3); 

GSE2HEE
xHEE = onera_desp_lib_rotate(xGSE,'gse2hee',matlabd);

HEE2GSE
xGSE = onera_desp_lib_rotate(xHEE,'hee2gse',matlabd);

HEE2HAE
xHAE = onera_desp_lib_rotate(xHEE,'hee2hae',matlabd);

HAE2HEE
xHEE = onera_desp_lib_rotate(xHAE,'hae2hee',matlabd);

HAE2HEEQ
xHEEQ = onera_desp_lib_rotate(xHAE,'hae2heeq',matlabd);

HEEQ2HAE
xHAE = onera_desp_lib_rotate(xHEEQ,'heeq2hae',matlabd);

This document is largely superseded by the project help pages
