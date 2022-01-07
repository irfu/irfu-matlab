%-Abstract
%
%   CSPICE_NEARPT locates the point on the surface of an ellipsoid that is
%   nearest to a specified position. This routine also returns the altitude
%   of the position above the ellipsoid.
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
%      positn   the array(s) defining the Cartesian position of a point with
%               respect to the center of an ellipsoid.
%
%               [3,n] = size(positn); double = class(positn)
%
%               The vector is expressed in a body-fixed reference frame.
%               The semi-axes of the ellipsoid are aligned with the x, y, and
%               z-axes of the body-fixed frame.
%
%      a        the length of the semi-axis of the ellipsoid that is parallel
%               to the X-axis of the body-fixed reference frame.
%
%               [1,1] = size(a); double = class(a)
%
%      b        the length of the semi-axis of the ellipsoid that is parallel
%               to the Y-axis of the body-fixed reference frame.
%
%               [1,1] = size(b); double = class(b)
%
%      c        the length of the semi-axis of the ellipsoid that is parallel
%               to the Z-axis of the body-fixed reference frame.
%
%               [1,1] = size(c); double = class(c)
%
%   the call:
%
%      [npoint, alt] = cspice_nearpt( positn, a, b, c )
%
%   returns:
%
%      npoint   the array(s) defining the location on the ellipsoid closest
%               to `positn'.
%
%               [3,n] = size(npoint); double = class(npoint)
%
%               `npoint' is a 3-vector expressed in the body-fixed
%               reference frame.
%
%      alt      the value(s) of the altitudes of `positn' above the ellipsoid.
%
%               [1,1] = size(alt); double = class(alt)
%
%               If `positn' is inside the ellipsoid, `alt' will be negative
%               and have magnitude equal to the distance between `npoint' and
%               `positn'.
%
%               `npoint' and `alt' return with the same vectorization
%               measure, N, as `positn'.
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
%   1) Given a point outside an ellipsoid, compute the nearest point
%      on its surface.
%
%      Example code begins here.
%
%
%      function nearpt_ex1()
%
%         %
%         % Define the radii of an ellipsoid.
%         %
%         a  =  1.;
%         b  =  2.;
%         c  =  3.;
%
%         %
%         % Use point on the X axis, outside the ellipsoid.
%         %
%         point        = [ 3.5; 0.; 0. ];
%         [pnear, alt] = cspice_nearpt( point, a, b, c);
%
%         fprintf('Nearest point: %6.2f %6.2f %6.2f\n', pnear)
%         fprintf('Altitude:      %6.2f\n', alt              )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Nearest point:   1.00   0.00   0.00
%      Altitude:        2.50
%
%
%   2) Compute the point on the Earth nearest to the Moon at ephemeris
%      time 0.0 (Jan 1 2000, 12:00 TBD).
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: nearpt_ex2.tm
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
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function nearpt_ex2()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'nearpt_ex2.tm')
%
%         %
%         % Retrieve the position of the Moon wrt the Earth at
%         % ephemeris time 0.d (Jan 1 2000 12:00 TDB) in the Earth-fixed
%         % reference frame.
%         %
%         epoch       = 0.;
%         abcorr      = 'LT+S';
%         [pos,lt] = cspice_spkpos( 'moon', epoch, 'IAU_EARTH', ...
%                                       abcorr, 'earth');
%
%         %
%         % Retrieve the triaxial radii for Earth (body ID 399).
%         %
%         radii = cspice_bodvrd( 'EARTH', 'RADII', 3);
%
%         %
%         % Now calculate the point on the Earth nearest to the Moon
%         % given LT+S aberration correction at the epoch time.
%         %
%         [npoint,alt]= cspice_nearpt( pos, radii(1), radii(2), radii(3) );
%
%         fprintf('Epoch:         %15.8f\n', epoch               )
%         fprintf('Nearest point: %15.8f %15.8f %15.8f\n', npoint)
%         fprintf('Altitude:      %15.8f\n', alt                 )
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
%      Epoch:              0.00000000
%      Nearest point:   3347.08204098  -5294.53585186  -1198.28264121
%      Altitude:      396037.22419372
%
%
%-Particulars
%
%   A sister version of this routine exists named mice_nearpt that returns
%   the output arguments as fields in a single structure.
%
%   Many applications of this routine are more easily performed
%   using the higher-level Mice routine cspice_subpnt. This routine
%   is the mathematical workhorse on which cspice_subpnt relies.
%
%   This routine solves for the location, N, on the surface of an
%   ellipsoid nearest to an arbitrary location, P, relative to that
%   ellipsoid.
%
%-Exceptions
%
%   1)  If any of the axis lengths `a', `b' or `c' are non-positive, the
%       error SPICE(BADAXISLENGTH) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If the ratio of the longest to the shortest ellipsoid axis
%       is large enough so that arithmetic expressions involving its
%       squared value may overflow, the error SPICE(BADAXISLENGTH)
%       is signaled by a routine in the call tree of this routine.
%
%   3)  If any of the expressions
%
%          a * abs( positn(1) ) / m^2
%          b * abs( positn(2) ) / m^2
%          c * abs( positn(3) ) / m^2
%
%       where `m' is the minimum of { `a', `b', `c' }, is large enough so
%       that arithmetic expressions involving these sub-expressions
%       may overflow, the error SPICE(INPUTSTOOLARGE) is signaled by a
%       routine in the call tree of this routine.
%
%   4)  If the axes of the ellipsoid have radically different
%       magnitudes, for example if the ratios of the axis lengths vary
%       by 10 orders of magnitude, the results may have poor
%       precision. No error checks are done to identify this problem.
%
%   5)  If the axes of the ellipsoid and the input point `positn' have
%       radically different magnitudes, for example if the ratio of
%       the magnitude of `positn' to the length of the shortest axis is
%       1.E25, the results may have poor precision. No error checks
%       are done to identify this problem.
%
%   6)  If any of the input arguments, `positn', `a', `b' or `c', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `positn', `a', `b' or `c', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   See -Exceptions section.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 2.0.0, 01-NOV-2021 (NJB) (JDR) (EDW)
%
%       Edit to logic to reduce unneeded operations when
%       error or projection vectors equal zero. Addition
%       of details concerning the "ellipsoid near point"
%       problem and solution.
%
%       Fixed a typo in the -I/O section.
%
%       Edited the header to comply with NAIF standard. Added
%       examples' problem statement, and meta-kernel for code example #2.
%       Updated code examples to produce formatted output and added
%       cspice_kclear to code example #2.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions, -Particulars,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 03-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 21-DEC-2005 (EDW)
%
%-Index_Entries
%
%   distance from point to ellipsoid
%   nearest point on an ellipsoid
%
%-&

function [ npoint, alt ] = cspice_nearpt( positn, a, b, c )

   switch nargin
      case 4

         positn = zzmice_dp(positn);
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);

      otherwise

         error ( ['Usage: [_npoint(3)_, _alt_] = ' ...
                  'cspice_nearpt( _positn(3)_, a, b, c )'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [nearpt] = mice( 'nearpt_s', positn, a, b, c  );
      npoint   = reshape( [nearpt.pos], 3, [] );
      alt      = reshape( [nearpt.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end


