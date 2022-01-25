%-Abstract
%
%   CSPICE_SCE2S converts an epoch specified as ephemeris seconds
%   past J2000 (ET) value describing a date to a character string
%   representation of a spacecraft clock value (SCLK).
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
%      sc       the NAIF ID of the spacecraft clock whose encoded SCLK value
%               at the epoch `et' is desired.
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
%      [sclkch] = cspice_sce2s( sc, et )
%
%   returns:
%
%      sclkch   the representation(s) of spacecraft `sc' clock count that
%               corresponds to `et'.
%
%               [n,c1] = size(sclkch); char = class(sclkch)
%
%               'sclkch' returns with the same vectorization measure (N)
%               as 'et'.
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
%   1) Convert a series of UTC times to the character string
%      representation of the CASSINI spacecraft clock values.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: sce2s_ex1.tm
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
%      function sce2s_ex1()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'sce2s_ex1.tm' )
%
%         %
%         % Assign values for the spacecraft ID (CASSINI).
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
%         % Convert the ephemeris time to the corresponding
%         % SCLK string appropriate for this spacecraft
%         %
%         sclkch = cspice_sce2s( SC, et );
%
%         disp( '         UTC Time               SCLK string'   )
%         disp( '--------------------------   ----------------' )
%         disp( 'Scalar:' )
%         txt    = sprintf( '%s   %s', event_time, sclkch );
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
%         % Convert the ephemeris time to the corresponding
%         % SCLK string appropriate for this spacecraft
%         %
%         sclkch = cspice_sce2s( SC, et );
%
%         disp( 'Vector:' )
%         for i=1:3
%            txt = sprintf( '%s   %s', event_time(i,:), sclkch(i,:));
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
%               UTC Time               SCLK string
%      --------------------------   ----------------
%      Scalar:
%      2004 JUN 11 11:00:37.57200   1/1465644305.079
%      Vector:
%      2004 JUN 11 12:00:37.57200   1/1465647905.085
%      2004 JUN 11 13:00:37.57200   1/1465651505.092
%      2004 JUN 11 14:00:37.57200   1/1465655105.098
%
%
%-Particulars
%
%   This routine is provided as a convenience; it is simply shorthand
%   for the code fragment
%
%      [sclkdp] = cspice_sce2t(  sc, et     );
%      [sclkch] = cspice_scdecd( sc, sclkdp );
%
%-Exceptions
%
%   1)  If an SCLK kernel has not been loaded, does not contain all of
%       the required data, or contains invalid data, an error is
%       signaled by a routine in the call tree of this routine. The
%       output argument `sclkch' will not be modified. This routine
%       assumes that an SCLK kernel appropriate to the spacecraft
%       clock identified by the input argument `sc' has been loaded.
%
%   2)  If a leapseconds kernel is required for conversion between
%       SCLK and `et' but is not loaded, an error is signaled by a
%       routine in the call tree of this routine. The output argument
%       `sclkch' will not be modified. When using an SCLK kernel that
%       maps SCLK to a time system other than `et' (also called
%       barycentric dynamical time---`TDB'), it is necessary to have a
%       leapseconds kernel loaded at the time this routine is called.
%
%       The time system to which an SCLK kernel maps SCLK epochs is
%       indicated by the variable SCLK_TIME_SYSTEM_nn in the kernel,
%       where nn is the negative of the NAIF integer code for the
%       spacecraft. The time system used in a kernel is TDB if and
%       only if the variable is assigned the value 1.
%
%   3)  If the input `et' value is not representable in the spacecraft
%       clock string format for the spacecraft clock identified by `sc',
%       an error is signaled by a routine in the call tree of this
%       routine. The output argument `sclkch' will not be modified.
%
%   4)  If any of the input arguments, `sc' or `et', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `sc' or `et', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   An SCLK kernel, appropriate to the spacecraft clock identified
%   by `sc', must be loaded at the time this routine is called.
%
%   If the SCLK kernel used with this routine does not map SCLK
%   directly to barycentric dynamical time, a leapseconds kernel
%   must be loaded at the time this routine is called.
%
%-Restrictions
%
%   None.
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
%   ephemeris time to spacecraft_clock string
%
%-&

function [sclkch] = cspice_sce2s(sc, et)

   switch nargin
      case 2

         sc = zzmice_int(sc);
         et = zzmice_dp(et);

      otherwise
         error ( 'Usage: [_`sclkch`_] = cspice_sce2s(sc, _et_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [sclkch] = mice('sce2s_c',sc, et);
   catch spiceerr
      rethrow(spiceerr)
   end



