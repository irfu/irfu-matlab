%-Abstract
%
%   CSPICE_SCT2E converts encoded spacecraft clock ("ticks")
%   to ephemeris seconds past J2000 (ET).
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
%      sc       a NAIF integer code for a spacecraft, one of whose encoded
%               clock values is represented by `sclkdp'.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%      sclkdp   the encoded spacecraft clock value(s).
%
%               [1,n] = size(sclkdp); double = class(sclkdp)
%
%               `sclkdp' represents time measured from spacecraft clock
%               start: partition information IS reflected in the encoded
%               value.
%
%   the call:
%
%      [et] = cspice_sct2e( sc, sclkdp )
%
%   returns:
%
%      et       the epoch(s), specified as ephemeris seconds past J2000, that
%               corresponds to `sclkdp'.
%
%               [1,n] = size(et); double = class(et)
%
%               `et' returns with the same vectorization measure, N,
%               as `sclkdp'.
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
%   1) Obtain from a series of CASSINI encoded spacecraft clock ("ticks")
%      values the corresponding ephemeris epochs and UTC times.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: sct2e_ex1.tm
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
%      function sct2e_ex1()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'sct2e_ex1.tm' )
%
%         %
%         % Assign values for the spacecraft ID (CASSINI).
%         %
%         SC     = -82;
%         sclkdp =  197483593295.0;
%
%         %
%         % Convert 'sclkdp' for spacecraft 'SC' to ephemeris time.
%         %
%         et = cspice_sct2e( SC, sclkdp );
%
%         %
%         % Convert the ephemeris time to a UTC calendar string.
%         %
%         utc = cspice_et2utc( et, 'C', 3 );
%
%         disp( '  Encoded SCLK          ET                  UTC Time')
%         disp( ['----------------  ----------------',  ...
%                '  ------------------------']       )
%         disp( 'Scalar:' )
%         fprintf( '%14.3f  %16.6f  %s\n', sclkdp, et, utc  );
%
%         %
%         % Convert a vector of SCLK values.
%         %
%         sclkdp = [ 197483587237.0,   ...
%                    197483587250.0,   ...
%                    197485901583.201, ...
%                    197486447183.0,   ...
%                    198136032015.400 ];
%
%         %
%         % Convert the `sclkdp' vector  for spacecraft `SC' to
%         % ephemeris time.
%         %
%         et = cspice_sct2e( SC, sclkdp );
%
%         %
%         % Convert the ephemeris time vector to a UTC calendar
%         % strings then output.
%         %
%         utc = cspice_et2utc( et, 'C', 3 );
%
%         disp( 'Vector:' )
%         for i=1:5
%            fprintf( '%14.3f  %16.6f  %s\n', sclkdp(i), et(i), utc(i,:));
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
%        Encoded SCLK          ET                  UTC Time
%      ----------------  ----------------  ------------------------
%      Scalar:
%      197483593295.000  140223701.758428  2004 JUN 11 11:00:37.574
%      Vector:
%      197483587237.000  140223678.094534  2004 JUN 11 11:00:13.910
%      197483587250.000  140223678.145315  2004 JUN 11 11:00:13.961
%      197485901583.201  140232718.444966  2004 JUN 11 13:30:54.260
%      197486447183.000  140234849.678995  2004 JUN 11 14:06:25.494
%      198136032015.400  142772272.350485  2004 JUL 10 22:56:48.167
%
%
%-Particulars
%
%   This routine operates on encoded SCLK values. These values
%   are obtained by calling the Mice routine cspice_scencd or other
%   SCLK conversion routines. The advantage of encoded SCLK, as
%   opposed to character string representations of SCLK is that
%   encoded SCLK values are easy to perform arithmetic operations on.
%   Additionally, working with encoded SCLK reduces the overhead of
%   repeated conversion of character strings to integers or double
%   precision numbers.
%
%   To convert a string representation of an SCLK value to ET, use
%   the Mice routine cspice_scs2e.
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
%
%   3)  If the clock type for the spacecraft clock identified by `sc' is
%       not supported by this routine, the error SPICE(NOTSUPPORTED)
%       is signaled by a routine in the call tree of this routine. The
%       output argument `et' will not be modified.
%
%   4)  If the input argument `sclkdp' is invalid, an error is signaled
%       by a routine in the call tree of this routine. The output
%       argument `et' will not be modified.
%
%   5)  If any of the input arguments, `sc' or `sclkdp', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   6)  If any of the input arguments, `sc' or `sclkdp', is not of the
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
%   TIME.REQ
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
%   -Mice Version 1.0.2, 06-JAN-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 04-SEP-2012 (EDW)
%
%       Edit to call example in -I/O to correct form.
%
%   -Mice Version 1.0.0, 18-APR-2006 (EDW)
%
%-Index_Entries
%
%   spacecraft_clock ticks to ephemeris time
%
%-&

function [et] = cspice_sct2e(sc,sclkdp)

   switch nargin
      case 2

         sc     = zzmice_int(sc);
         sclkdp = zzmice_dp(sclkdp);

      otherwise
         error( 'Usage: [_et_] = cspice_sct2e(sc, _sclkdp_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [et] = mice('sct2e_c',sc, sclkdp);
   catch spiceerr
      rethrow(spiceerr)
   end
