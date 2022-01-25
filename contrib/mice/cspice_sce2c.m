%-Abstract
%
%   CSPICE_SCE2C converts ephemeris seconds past J2000 (ET) to
%   continuous encoded spacecraft clock ("ticks"). Non-integral
%   tick values may be returned.
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
%      sc       the NAIF ID of the spacecraft clock whose
%               encoded SCLK value at the epoch `et' is desired.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%      et       the ephemeris time(s) expressed as ephemeris seconds
%               past J2000.
%
%               [1,n] = size(et); double = class(et)
%
%   the call:
%
%      sclkdp = cspice_sce2c( sc, et )
%
%   returns:
%
%      sclkdp   the encoded SCLK value(s) corresponding to `et'
%               for `sc'.
%
%               [1,n] = size(sclkdp); double = class(sclkdp)
%
%               `sclkdp' returns with the same vectorization measure, N,
%               as `et'.
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
%   1) Convert a series of UTC times to their corresponding CASSINI
%      continuous encoded spacecraft clock ("ticks") values.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: sce2c_ex1.tm
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
%            ---------                     --------
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
%      function sce2c_ex1()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'sce2c_ex1.tm' )
%
%         %
%         % Assign values for the spacecraft ID (CASSINI),
%         % and event time.
%         %
%         SC         = -82;
%         event_time = '2004 JUN 11 11:00:37.57200';
%
%         %
%         % Convert the time string to ephemeris time.
%         %
%         et = cspice_str2et( event_time );
%
%         %
%         % Convert the ephemeris time to the encoded SCLK
%         % format.
%         %
%         sclkdp = cspice_sce2c( SC, et );
%
%         disp( '         UTC Time               Encoded SCLK'    )
%         disp( '--------------------------  -------------------' )
%         disp( 'Scalar:' )
%         txt    = sprintf( '%s  %19.6f', event_time, sclkdp );
%         disp( txt )
%
%         %
%         % Vectorized use, a vector of UTC times.
%         %
%         event_time =  strvcat( '2004 JUN 11 12:00:37.57200', ...
%                                '2004 JUN 11 13:00:37.57200', ...
%                                '2004 JUN 11 14:00:37.57200' );
%
%         %
%         % Convert the time strings to ET.
%         %
%         et = cspice_str2et( event_time );
%
%         %
%         % Convert the 'et' array to the encoded
%         % spacecraft clock.
%         %
%         sclkdp = cspice_sce2c( SC, et );
%
%         disp( 'Vector:' )
%         for i=1:3
%            txt = sprintf( '%s  %19.6f', event_time(i,:), sclkdp(i) );
%            disp( txt )
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
%               UTC Time               Encoded SCLK
%      --------------------------  -------------------
%      Scalar:
%      2004 JUN 11 11:00:37.57200  197483593294.540863
%      Vector:
%      2004 JUN 11 12:00:37.57200  197484514901.107330
%      2004 JUN 11 13:00:37.57200  197485436507.673767
%      2004 JUN 11 14:00:37.57200  197486358114.240204
%
%
%-Particulars
%
%   This routine outputs continuous encoded SCLK values; unlike the
%   routine cspice_sce2t, the values output by this routine need not be
%   integral.
%
%   This routine supports use of non-integral encoded clock values in
%   C-kernels: non-integral clock values may be stored as pointing
%   time tags when a C-kernel is created, and they may be supplied
%   as request times to the C-kernel readers.
%
%   The advantage of encoded SCLK, as opposed to character string
%   representations of SCLK, is that encoded SCLK values are easy to
%   perform arithmetic operations on. Also, working with encoded SCLK
%   reduces the overhead of repeated conversion of  character strings
%   to integers or double precision numbers.
%
%   To convert ET to a string representation of an SCLK value, use
%   the Mice routine cspice_sce2s.
%
%   See the SCLK Required Reading for a list of the entire set of SCLK
%   conversion routines.
%
%-Exceptions
%
%   1)  If an SCLK kernel has not been loaded, does not contain all of
%       the required data, or contains invalid data, an error is
%       signaled by a routine in the call tree of this routine. The
%       output argument `sclkdp' will not be modified. This routine
%       assumes that that an SCLK kernel appropriate to the spacecraft
%       clock identified by the input argument `sc' has been loaded.
%
%   2)  If a leapseconds kernel is required for conversion between
%       SCLK and `et' but is not loaded, an error is signaled by a
%       routine in the call tree of this routine. The output argument
%       `sclkdp' will not be modified. When using SCLK kernels that map
%       SCLK to a time system other than `et' (also called barycentric
%       dynamical time---`TDB'), it is necessary to have a leapseconds
%       kernel loaded at the time this routine is called.
%
%       The time system that an SCLK kernel maps SCLK to is indicated
%       by the variable SCLK_TIME_SYSTEM_nn in the kernel, where nn
%       is the negative of the NAIF integer code for the spacecraft.
%       The time system used in a kernel is TDB if and only if the
%       variable is assigned the value 1.
%
%   3)  If the clock type for the spacecraft clock identified by `sc' is
%       not supported by this routine, the error SPICE(NOTSUPPORTED)
%       is signaled by a routine in the call tree of this routine. The
%       output argument `sclkdp' will not be modified.
%
%   4)  If the input `et' value is not representable as an encoded
%       spacecraft clock value for the spacecraft clock identified by
%       `sc', an error is signaled by a routine in the call tree of this
%       routine. The output argument `sclkdp' will not be modified.
%
%   5)  If any of the input arguments, `sc' or `et', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   6)  If any of the input arguments, `sc' or `et', is not of the
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
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
%   ephemeris time to continuous spacecraft_clock ticks
%
%-&

function [sclkdp] = cspice_sce2c(sc, et)

   switch nargin
      case 2

         sc = zzmice_int(sc);
         et = zzmice_dp(et);

      otherwise
         error ( 'Usage: [_sclkdp_] = cspice_sce2c(sc, _et_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [sclkdp] = mice('sce2c_c',sc, et);
   catch spiceerr
      rethrow(spiceerr)
   end



