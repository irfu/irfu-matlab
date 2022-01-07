%-Abstract
%
%   CSPICE_SCS2E converts a spacecraft clock string to ephemeris
%   seconds past J2000 (ET).
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      sc       a NAIF ID code for a spacecraft, one of whose clock values is
%               represented by `sclkch'.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%               The set of supported spacecraft clocks is listed in the
%               SCLK Required Reading.
%
%      sclkch   character string representation(s) of the spacecraft clock
%               value that corresponds to `et', for the spacecraft specified
%               by the input argument `sc'.
%
%               [n,c1] = size(sclkch); char = class(sclkch)
%
%                  or
%
%               [1,n] = size(sclkch); cell = class(sclkch)
%
%               `sclkch' is an absolute spacecraft clock time, so partition
%               information should be included in this string. The precise
%               format of `sclkch' is specified in the SCLK Required Reading.
%
%   the call:
%
%      [et] = cspice_scs2e( sc, sclkch )
%
%   returns:
%
%      et       the epoch(s), specified as ephemeris seconds past J2000, that
%               corresponds to `sclkch'.
%
%               [1,n] = size(et); double = class(et)
%
%               `et' returns with the same vectorization measure, N,
%               as `sclkch'.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Obtain from a series of CASSINI SCLK string representations
%      the corresponding ephemeris epoch and UTC time.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: scs2e_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     ------------
%            naif0012.tls                  Leapseconds
%            cas00071.tsc                  CASSINI SCLK
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                'cas00071.tsc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function scs2e_ex1()
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'scs2e_ex1.tm' )
%
%         %
%         % Assign values for the spacecraft ID (CASSINI).
%         %
%         SC     = -82;
%         sclkch = '1/1465644281.165';
%
%         %
%         % Convert 'sclkch' for spacecraft 'SC' to ephemeris time.
%         %
%         et = cspice_scs2e( SC, sclkch );
%
%         %
%         % Convert the ephemeris time to a UTC calendar string.
%         %
%         utc = cspice_et2utc( et, 'C', 3 );
%
%         disp( '  SCLK String           ET                  UTC Time')
%         disp( ['----------------  ----------------', ...
%                '  ------------------------']       )
%         disp( 'Scalar:' )
%         fprintf( '%s  %16.6f  %s\n', sclkch, et, utc  );
%
%         %
%         % Convert a vector of SCLK strings to ET and
%         % UTC.
%         %
%         sclkch =  strvcat( '1/1465644281.165' , ...
%                            '1/1465646281.165' , ...
%                            '1/1465647281'     , ...
%                            '1/1465647281.001' );
%
%         et  = cspice_scs2e( SC, sclkch );
%         utc = cspice_et2utc( et, 'C', 3 );
%
%         disp( 'Vector:' )
%         for i=1:4
%            fprintf( '%s  %16.6f  %s\n', sclkch(i,:), et(i), utc(i,:));
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%        SCLK String           ET                  UTC Time
%      ----------------  ----------------  ------------------------
%      Scalar:
%      1/1465644281.165  140223678.094534  2004 JUN 11 11:00:13.910
%      Vector:
%      1/1465644281.165  140223678.094534  2004 JUN 11 11:00:13.910
%      1/1465646281.165  140225678.080283  2004 JUN 11 11:33:33.896
%      1/1465647281      140226677.428631  2004 JUN 11 11:50:13.244
%      1/1465647281.001  140226677.432538  2004 JUN 11 11:50:13.248
%
%
%-Particulars
%
%   This routine is provided as a convenience; it is simply shorthand
%   for the code fragment
%
%      [sclkdp] = cspice_scencd( sc, sclkch );
%      [et]     = cspice_sct2e(  sc, sclkdp );
%
%-Exceptions
%
%   1)  If an SCLK kernel has not been loaded, does not contain all of
%       the required data, or contains invalid data, an error is
%       signaled by a routine in the call tree of this routine. The
%       output argument `et' will not be modified. This routine assumes
%       that that an SCLK kernel appropriate to the spacecraft clock
%       identified by the input argument `sc' has been loaded.
%
%   2)  If a leapseconds kernel is required for conversion between
%       SCLK and `et' but is not loaded, an error is signaled by a
%       routine in the call tree of this routine. The output argument
%       `et' will not be modified. When using SCLK kernels that map SCLK
%       to a time system other than `et' (also called barycentric
%       dynamical time---`TDB'), it is necessary to have a leapseconds
%       kernel loaded at the time this routine is called.
%
%       The time system that an SCLK kernel maps SCLK to is indicated
%       by the variable SCLK_TIME_SYSTEM_nn in the kernel, where nn
%       is the negative of the NAIF integer code for the spacecraft.
%       The time system used in a kernel is TDB if and only if the
%       variable is assigned the value 1.
%
%   3)  If the value of `sclkch' is invalid, an error is signaled by a
%       routine in the call tree of this routine. The output argument
%       `et' will not be modified.
%
%   4)  If any of the input arguments, `sc' or `sclkch', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `sc' or `sclkch', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  An SCLK kernel appropriate to the spacecraft clock identified
%       by `sc' must be loaded at the time this routine is called.
%
%   2)  If the SCLK kernel used with this routine does not map SCLK
%       directly to barycentric dynamical time, a leapseconds kernel
%       must be loaded at the time this routine is called.
%
%-Required_Reading
%
%   MICE.REQ
%   SCLK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and meta-kernel with CASSINI PDS
%       archived data. Reformatted code example output.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 06-JAN-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 18-APR-2006 (EDW)
%
%-Index_Entries
%
%   spacecraft_clock string to ephemeris time
%
%-&

function [et] = cspice_scs2e(sc, sclkch)

   switch nargin
      case 2

         sc     = zzmice_int(sc);
         sclkch = zzmice_str(sclkch);

      otherwise
         error ( 'Usage: [_et_] = cspice_scs2e(sc, _`sclkch`_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [et] = mice('scs2e_c',sc, sclkch);
   catch spiceerr
      rethrow(spiceerr)
   end





