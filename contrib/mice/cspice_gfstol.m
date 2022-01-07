%-Abstract
%
%   CSPICE_GFSTOL overrides the default GF convergence value used in the high
%   level GF routines.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      value    value to use as the GF subsystem convergence tolerance.
%
%               [1,1] = size(value); double = class(value)
%
%               This value will override the default tolerance,
%               SPICE_GF_CNVTOL, defined in MiceGF.m. Units are TDB seconds.
%
%   the call:
%
%      cspice_gfstol( value )
%
%   returns:
%
%      None.
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
%   1) In 14 A.D., the Roman princeps Tiberius sent his son Drusus to subdue
%      a revolt of a Roman Legion stationed in Pannonia. A Lunar eclipse
%      occurred during this mission.
%
%      Perform a search for occultation events of the sun by earth as
%      observed from the Moon center. Search during the interval from
%      14 A.D. SEP 1 to 14 A.D. SEP 30 (Julian).
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gfstol_ex1.tm
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
%            de408.bsp                     Planetary ephemeris covering
%                                          year 14 AD
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de408.bsp',
%                                'pck00009.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End meta-kernel
%
%
%      Example code begins here.
%
%
%      function gfstol_ex1()
%
%         TIMFMT  = 'YYYY ERA MON DD HR:MN:SC.#### ::JCAL';
%         MAXWIN  = 100;
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfstol_ex1.tm' )
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '14 A.D. SEP 1  00:00:00',                 ...
%                               '14 A.D. SEP 30 00:00:00'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % Select a 3-minute step. We'll ignore any occultations
%         % lasting less than 3 minutes.
%         %
%         step    = 180.;
%
%         occtyp  = 'any';
%         front   = 'earth';
%         fshape  = 'ellipsoid';
%         fframe  = 'iau_earth';
%         back    = 'sun';
%         bshape  = 'ellipsoid';
%         bframe  = 'iau_sun';
%         obsrvr  = 'moon';
%         abcorr  = 'lt';
%
%         %
%         % Perform the search. 'et(1)' and 'et(2)' have values ~-6*10^10,
%         % SPICE_GF_CNVTOL has value 10^-6, so double precision addition or
%         % subtraction of 'et(1)' and 'et(2)' with SPICE_GF_CNVTOL returns
%         % a result indistinguishable from 'et(1)' and 'et(2)'.
%         %
%         % Reduce the GF convergence tolerance by an order of magnitude
%         % to resolve this condition.
%         %
%
%         cspice_gfstol( 1e-5 )
%
%         result = cspice_gfoclt( occtyp, front,  fshape, fframe,          ...
%                                 back,   bshape, bframe, abcorr,          ...
%                                 obsrvr, step,   cnfine, MAXWIN);
%
%         %
%         % List the beginning and ending times in each interval
%         % if result contains data.
%         %
%         count = cspice_wncard(result);
%
%         for i=1:count
%
%            [left, right] = cspice_wnfetd( result, i );
%
%            output = cspice_timout( [left,right], TIMFMT );
%
%            if( isequal( left, right) )
%
%               disp( ['Event time: ' output(1,:)] )
%
%            else
%
%               disp( ['Start time :' output(1,:)] )
%               disp( ['Stop time  :' output(2,:)] )
%               disp( ' ')
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Start time :  14 A.D. SEP 27 05:02:02.8250
%      Stop time  :  14 A.D. SEP 27 09:33:31.6995
%
%
%-Particulars
%
%   The high level GF routines (see gf.req for a listing) use a default
%   value for the convergence tolerance, SPICE_GF_CNVTOL, defined in
%   MiceGF.m. It may occur that a GF search run needs a different
%   convergence tolerance. cspice_gfstol programmatically changes the
%   tolerance used by those routines.
%
%-Exceptions
%
%   1)  If `value' is not strictly greater-than-zero, the error
%       SPICE(INVALIDTOLERANCE) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the input argument `value' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `value' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   GF.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 21-JUL-2020 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added -Parameters,
%       -Exceptions, -Files, -Restrictions, -Literature_References and
%       -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 12-MAR-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   change default convergence tolerance for GF routines
%
%-&

function cspice_gfstol( value )

  switch nargin

      case 1

         value = zzmice_dp(value);

      otherwise

         error ( 'Usage: cspice_gfstol(value)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('gfstol_c', value);
   catch spiceerr
      rethrow(spiceerr)
   end


