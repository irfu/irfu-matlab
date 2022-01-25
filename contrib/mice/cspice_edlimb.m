%-Abstract
%
%   CSPICE_EDLIMB calculates the limb of a triaxial ellipsoid
%   as viewed from a specified location.
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
%      a,
%      b,
%      c        the lengths of the semi-axes of a triaxial ellipsoid.
%
%               [1,1] = size(a); double = class(a)
%               [1,1] = size(b); double = class(b)
%               [1,1] = size(c); double = class(c)
%
%               The ellipsoid is centered at the origin and oriented so that
%               its axes lie on the x, y and z axes. `a', `b', and `c' are
%               the lengths of the semi-axes that respectively point in the
%               x, y, and z directions.
%
%      viewpt   a point from which the ellipsoid is viewed. `viewpt' must be
%               outside of the ellipsoid.
%
%               [3,1] = size(viewpt); double = class(viewpt)
%
%   the call:
%
%      [limb] = cspice_edlimb( a, b, c, viewpt )
%
%   returns:
%
%      limb   the SPICE ellipse that represents the limb of the ellipsoid
%             observed from `viewpt'.
%
%              [1,1] = size(limb); struct = class(limb)
%
%              The structure has the fields:
%
%                 center:    [3x1 double]
%                 semiMajor: [3x1 double]
%                 semiMinor: [3x1 double]
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
%   1) Given an ellipsoid and a viewpoint exterior to it, calculate
%      the limb ellipse as seen from that viewpoint.
%
%      Example code begins here.
%
%
%      function edlimb_ex1()
%
%         %
%         % Define an ellipsoid
%         %
%         a = sqrt(2.);
%         b = 2.*sqrt(2.);
%         c = sqrt(2.);
%
%         %
%         % Locate a viewpoint exterior to the ellipsoid.
%         %
%         viewpt = [ 2., 0.,  0. ]';
%
%         %
%         % Calculate the limb ellipse as seen by from the viewpoint.
%         %
%         limb = cspice_edlimb( a, b, c, viewpt );
%
%         %
%         % Output the structure components.
%         %
%         fprintf( 'Semiminor axis: %10.3f %10.3f %10.3f\n', ...
%                  limb.semiMinor                          );
%         fprintf( 'Semimajor axis: %10.3f %10.3f %10.3f\n', ...
%                  limb.semiMajor                          );
%         fprintf( 'Center:         %10.3f %10.3f %10.3f\n', ...
%                  limb.center                             );
%
%         %
%         % Check against expected values:
%         %
%         % Semiminor: 0., 0., -1.
%         % Semimajor: 0., 2.,  0.
%         % Center   : 1., 0.,  0.
%         %
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Semiminor axis:      0.000      0.000     -1.000
%      Semimajor axis:      0.000      2.000     -0.000
%      Center:              1.000      0.000      0.000
%
%
%   2) We'd like to find the apparent limb of Jupiter, corrected for
%      light time and stellar aberration, as seen from JUNO
%      spacecraft's position at a given UTC time.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: edlimb_ex2.tm
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
%            File name                           Contents
%            ---------                           --------
%            juno_rec_160522_160729_160909.bsp   JUNO s/c ephemeris
%            pck00010.tpc                        Planet orientation
%                                                and radii
%            naif0012.tls                        Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'juno_rec_160522_160729_160909.bsp',
%                                'pck00010.tpc',
%                                'naif0012.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function edlimb_ex2()
%
%         %
%         % Local parameters.
%         %
%         UTCSTR = '2016 Jul 14 19:45:00';
%
%         %
%         % Load the required kernels.
%         %
%         cspice_furnsh( 'edlimb_ex2.tm' );
%
%         %
%         % Find the viewing point in Jupiter-fixed coordinates. To
%         % do this, find the apparent position of Jupiter as seen
%         % from the spacecraft in Jupiter-fixed coordinates and
%         % negate this vector. In this case we'll use light time
%         % and stellar aberration corrections to arrive at the
%         % apparent limb. `jpos' is the Jupiter's position as seen
%         % from the spacecraft.  `scpos' is the spacecraft's position
%         % relative to Jupiter.
%         %
%         [et]       = cspice_str2et( UTCSTR );
%         [jpos, lt] = cspice_spkpos( 'JUPITER', et,    'J2000',           ...
%                                     'LT+S',    'JUNO'          );
%
%         scpos = -jpos;
%
%         %
%         % Get Jupiter's semi-axis lengths...
%         %
%         [rad] = cspice_bodvrd( 'JUPITER', 'RADII', 3 );
%
%         %
%         % ...and the transformation from J2000 to Jupiter
%         % equator and prime meridian coordinates. Note that we
%         % use the orientation of Jupiter at the time of
%         % emission of the light that arrived at the
%         % spacecraft at time `et'.
%         %
%         [tipm] = cspice_pxform( 'J2000', 'IAU_JUPITER', et-lt );
%
%         %
%         % Transform the spacecraft's position into Jupiter-
%         % fixed coordinates.
%         %
%         scpjfc = tipm * scpos;
%
%         %
%         % Find the apparent limb.  `limb' is a SPICE ellipse
%         % representing the limb.
%         %
%         [limb] = cspice_edlimb( rad(1), rad(2), rad(3), scpjfc );
%
%         %
%         % `center', `smajor', and `sminor' are the limb's center,
%         % semi-major axis of the limb, and a semi-minor axis
%         % of the limb.  We obtain these from `limb' using the
%         % Mice routine cspice_el2cgv ( Ellipse to center and
%         % generating vectors ).
%         %
%         [center, smajor, sminor] = cspice_el2cgv( limb );
%
%         %
%         % Output the structure components.
%         %
%         fprintf( 'Apparent limb of Jupiter as seen from JUNO:\n' )
%         fprintf( '   UTC time       : %s\n', UTCSTR )
%         fprintf( '   Semi-minor axis: %13.6f %13.6f %13.6f\n',           ...
%                               sminor(1), sminor(2), sminor(3) )
%         fprintf( '   Semi-major axis: %13.6f %13.6f %13.6f\n',           ...
%                               smajor(1), smajor(2), smajor(3) )
%         fprintf( '   Center         : %13.6f %13.6f %13.6f\n',           ...
%                               center(1), center(2), center(3) )
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
%      Apparent limb of Jupiter as seen from JUNO:
%         UTC time       : 2016 Jul 14 19:45:00
%         Semi-minor axis:  12425.547643  -5135.572410  65656.053303
%         Semi-major axis:  27305.667297  66066.222576     -0.000000
%         Center         :    791.732472   -327.228993   -153.408849
%
%
%-Particulars
%
%   The limb of a body, as seen from a viewing point, is the boundary
%   of the portion of the body's surface that is visible from that
%   viewing point. In this definition, we consider a surface point
%   to be `visible' if it can be connected to the viewing point by a
%   line segment that doesn't pass through the body. This is a purely
%   geometrical definition that ignores the matter of which portions
%   of the surface are illuminated, or whether the view is obscured by
%   any additional objects.
%
%   If a body is modeled as a triaxial ellipsoid, the limb is always
%   an ellipse. The limb is determined by its center, a semi-major
%   axis vector, and a semi-minor axis vector.
%
%   We note that the problem of finding the limb of a triaxial
%   ellipsoid is mathematically identical to that of finding its
%   terminator, if one makes the simplifying assumption that the
%   terminator is the limb of the body as seen from the vertex of the
%   umbra. So, this routine can be used to solve this simplified
%   version of the problem of finding the terminator.
%
%-Exceptions
%
%   1)  If the length of any semi-axis of the ellipsoid is
%       non-positive, the error SPICE(INVALIDAXISLENGTH) is signaled
%       by a routine in the call tree of this routine. `limb' is not
%       modified.
%
%   2)  If the length of any semi-axis of the ellipsoid is zero after
%       the semi-axis lengths are scaled by the reciprocal of the
%       magnitude of the longest semi-axis and then squared, the error
%       SPICE(DEGENERATECASE) is signaled by a routine in the call
%       tree of this routine. `limb' is not modified.
%
%   3)  If the viewing point `viewpt' is inside the ellipse, the error
%       SPICE(INVALIDPOINT) is signaled by a routine in the call tree
%       of this routine. `limb' is not modified.
%
%   4)  If the geometry defined by the input ellipsoid and viewing
%       point is so extreme that the limb cannot be found, the error
%       SPICE(DEGENERATECASE) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the shape of the ellipsoid and the viewing geometry are
%       such that the limb is an excessively flat ellipsoid, the
%       limb may be a degenerate ellipse. You must determine whether
%       this possibility poses a problem for your application.
%
%   6)  If any of the input arguments, `a', `b', `c' or `viewpt', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `a', `b', `c' or `viewpt', is
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
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   ELLIPSES.REQ
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
%   -Mice Version 1.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard.
%       Added example's problem statement and added second example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 09-NOV-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   ellipsoid limb
%
%-&

function [ limb ] = cspice_edlimb( a, b, c, viewpt )

   switch nargin
      case 4

         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
         viewpt = zzmice_dp(viewpt);

      otherwise

         error ( 'Usage: [ limb ] = cspice_edlimb( a, b, c, viewpt )' )

   end

   try
      [ limb ] = mice( 'edlimb_s', a, b, c, viewpt );
   catch spiceerr
      rethrow(spiceerr)
   end
