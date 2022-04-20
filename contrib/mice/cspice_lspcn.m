%-Abstract
%
%   CSPICE_LSPCN computes L_s, the planetocentric longitude of the sun,
%   as seen from a specified body.
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
%      body     the name of the central body, typically a planet.
%
%               [1,c1] = size(body); char = class(body)
%
%                  or
%
%               [1,1] = size(body); cell = class(body)
%
%      et       the epoch(s) at which the longitude of the sun (L_s) is
%               to be computed. `et' is expressed as seconds past J2000
%               TDB (Barycentric Dynamical Time).
%
%               [1,n] = size(et); double = class(et)
%
%      abcorr   indicates the aberration corrections to be applied
%               when computing the longitude of the sun.
%
%               [1,c2] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               `abcorr' may be any of the following.
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       Correct the position of the sun,
%                             relative to the central body, for
%                             planetary (light time) aberration.
%
%                  'LT+S'     Correct the position of the sun,
%                             relative to the central body, for
%                             planetary and stellar aberrations.
%
%   the call:
%
%      [lspcn] = cspice_lspcn( body, et, abcorr )
%
%   returns:
%
%      lspcn    the planetocentric longitude(s) of the sun, often called
%               "L_s," for the specified body at the specified time.
%
%               [1,n] = size(et); double = class(et)
%
%               The longitude is defined in a right-handed frame whose
%               basis vectors are defined as follows:
%
%               - The positive Z direction is given by the instantaneous
%                 angular velocity vector of the orbit of the body about
%                 the sun.
%
%               - The positive X direction is that of the cross product
%                 of the instantaneous north spin axis of the body with
%                 the positive Z direction.
%
%               - The positive Y direction is Z x X.
%
%               Units are radians; the range is 0 to 2*pi. Longitudes are
%               positive to the east.
%
%               `lspcn' returns with the same vectorization measure (N) as
%               `et'.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Compute the planetocentric longitude of the Sun as seen from
%      the Earth on 21 March 2006. The result should be approximately
%      0, as the date corresponds to the Spring equinox.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: lspcn_ex1.tm
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
%            de421.bsp                     Planetary ephemeris
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            naif0011.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc',
%                                'naif0011.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function lspcn_ex1()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'lspcn_ex1.tm' )
%
%         et = cspice_str2et('21 march 2006');
%         lspcn = cspice_lspcn( 'earth', et, 'none' ) * cspice_dpr;
%         disp( lspcn )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%       0.23645
%
%
%   2) Compute the planetocentric longitude of the Sun as seen from
%      the Earth for a set of dates.
%
%      Use the meta-kernel from the first example load the required SPICE
%      kernels.
%
%      Example code begins here.
%
%
%      function lspcn_ex2()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'lspcn_ex1.tm' )
%
%         et = cspice_str2et('21 march 2005') + [0:1000000.:10000000.];
%         lspcn = cspice_lspcn( 'earth', et, 'none' ) * cspice_dpr;
%
%         utcstr = cspice_et2utc(et, 'C', 0);
%
%
%         disp( '      UTC time           lspcn' )
%         disp( '--------------------   --------')
%         for i = 1:11
%            fprintf( '%s   %8.4f\n', utcstr(i,:), lspcn(i) )
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%            UTC time           lspcn
%      --------------------   --------
%      2005 MAR 21 00:00:00     0.4815
%      2005 APR 01 13:46:40    11.9353
%      2005 APR 13 03:33:20    23.3193
%      2005 APR 24 17:20:00    34.6294
%      2005 MAY 06 07:06:40    45.8611
%      2005 MAY 17 20:53:20    57.0476
%      2005 MAY 29 10:40:00    68.1626
%      2005 JUN 10 00:26:40    79.2509
%      2005 JUN 21 14:13:20    90.3058
%      2005 JUL 03 04:00:00   101.3366
%      2005 JUL 14 17:46:40   112.3833
%
%
%-Particulars
%
%   The direction of the vernal equinox for the central body is
%   determined from the instantaneous equatorial and orbital planes
%   of the central body. This equinox definition is specified in
%   reference [1]. The "instantaneous orbital plane" is interpreted
%   in this routine as the plane normal to the cross product of the
%   position and velocity of the central body relative to the sun.
%   The geometric state of the central body relative to the sun is
%   used for this normal vector computation. The "instantaneous
%   equatorial plane" is normal to the central body's north pole
%   at the requested epoch. The pole direction is determined from
%   rotational elements loaded via a PCK file.
%
%   The result returned by this routine will depend on the
%   ephemeris data and rotational elements used. The result may
%   differ from that given in any particular version of the
%   Astronomical Almanac, due to differences in these input data,
%   and due to differences in precision of the computations.
%
%-Exceptions
%
%   1)  If the input body name cannot be translated to an ID code,
%       and if the name is not a string representation of an integer
%       (for example, '399'), the error SPICE(NOTRANSLATION) is
%       signaled by a routine in the call tree of this routine.
%
%   2)  If no SPK (ephemeris) file has been loaded prior to calling
%       this routine, or if the SPK data has insufficient coverage, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   3)  If a PCK file containing rotational elements for the central
%       body has not been loaded prior to calling this routine, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   4)  If the instantaneous angular velocity and spin axis of `body'
%       are parallel, an error is signaled by a
%       routine in the call tree of this routine.
%
%   5)  If any of the input arguments, `body', `et' or `abcorr', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `body', `et' or `abcorr', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   Appropriate SPICE kernels must be loaded by the calling program
%   before this routine is called.
%
%   The following data are required:
%
%   -  An SPK file (or files) containing ephemeris data sufficient to
%      compute the geometric state of the central body relative to
%      the sun at `et' must be loaded before this routine is called. If
%      light time correction is used, data must be available that
%      enable computation of the state the sun relative to the solar
%      system barycenter at the light-time corrected epoch. If
%      stellar aberration correction is used, data must be available
%      that enable computation of the state the central body relative
%      to the solar system barycenter at `et'.
%
%   -  A PCK file containing rotational elements for the central body
%      must be loaded before this routine is called.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   ABCORR.REQ
%
%-Literature_References
%
%   [1]  "The Astronomical Almanac for the Year 2005," page L9,
%        United States Naval Observatory, U.S. Government Printing
%        Office, Washington, D.C., 2004.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       meta-kernel to example #1. Updated code examples to produce
%       formatted output and added a call to cspice_kclear. Added the
%       problem statement to both examples.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       square brackets to output argument in function declaration.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 10-NOV-2015 (EDW)
%
%       Script rewritten to call CSPICE interface.
%
%   -Mice Version 0.9.0, 23-OCT-2015 (EDW)
%
%       Pure Matlab script.
%
%-Index_Entries
%
%   planetocentric longitude of sun
%   compute L_s
%   compute Ls
%   compute L_sub_s
%
%-&

function [lspcn] = cspice_lspcn( body, et, abcorr )

   switch nargin
      case 3

         body   = zzmice_str(body);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);

      otherwise

         error( 'Usage: [_lspcn_] = cspice_lspcn( `body`, _et_, `abcorr`)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [lspcn] = mice('lspcn_c', body, et, abcorr );
   catch spiceerr
      rethrow(spiceerr)
   end
