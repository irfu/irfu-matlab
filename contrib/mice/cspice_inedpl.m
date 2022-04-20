%-Abstract
%
%   CSPICE_INEDPL calculates the intercept of a triaxial ellipsoid
%   and a plane.
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
%      plane    a structure describing a SPICE plane.
%
%               [1,1] = size(plane); struct = class(plane)
%
%               The intersection of `plane' and the ellipsoid is sought.
%
%               The structure has the fields:
%
%                  normal:   [3,1] = size(normal); double = class(normal)
%                  constant: [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [ellips, found] = cspice_inedpl( a, b, c, plane )
%
%   returns:
%
%      ellips   a structure describing a SPICE ellipse that defines the
%               intersection of `plane' and the ellipsoid.
%
%               [1,1] = size(ellips); struct = class(ellips)
%
%               The structure has the fields:
%
%                  center:    [3,1] =  size(center); double = class(center)
%                  semiMajor: [3,1] =  size(semiMajor);
%                            double = class(semiMajor)
%
%      found    the boolean indicating whether `plane'
%               intersects the ellipsoid (true) or not (false).
%
%               [1,1] = size(found); logical = class(found)
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
%   1) Suppose we wish to find the limb of a body, as observed from
%      location `loc' in body-fixed coordinates. The Mice routine
%      cspice_edlimb solves this problem. Here's how cspice_inedpl is used in
%      that solution.
%
%      We assume `loc' is outside of the body. The body is modelled as
%      a triaxial ellipsoid with semi-axes of length `a', `b', and `c'.
%
%      The notation
%
%         < x, y >
%
%      indicates the inner product of the vectors `x' and `y'.
%
%      The limb lies on the plane defined by
%
%         < x,  n >  =  1,
%
%      where the vector `n' is defined as
%
%                     2              2              2
%         ( loc(1) / a ,   loc(2) / b ,   loc(3) / c  )
%
%      The assignments
%
%         n(1) = loc(1) / (a*a);
%         n(2) = loc(2) / (b*b);
%         n(3) = loc(3) / (c*c);
%
%      and the calls
%
%         [plane]                  = cspice_nvc2pl( n,  1.0 );
%
%         [limb, found]            = cspice_inedpl( a,  b,  c,  plane );
%
%         [center, smajor, sminor] = cspice_el2cgv( limb );
%
%      will return the center and semi-axes of the limb.
%
%
%      How do we know that  < x, n > = 1  for all `x' on the limb?
%      This is because all limb points `x' satisfy
%
%         < loc - x, surfnm(x) >  =  0,
%
%      where surfnm(x) is any surface normal at `x'. surfnm(x) is
%      parallel to the vector
%
%                        2            2            2
%         v = (  x(1) / a ,   x(2) / b ,   x(3) / c   )
%
%      so we have
%
%         < loc - x, v >  =  0,
%
%         < loc, v >      =  < x, v >  =  1  (from the original
%                                             ellipsoid
%                                             equation)
%      and finally
%
%         < x, n >  =  1
%
%      where `n' is as defined above.
%
%
%   2) We'd like to find the apparent limb of Jupiter, corrected for
%      light time and stellar aberration, as seen from JUNO
%      spacecraft's position at a given UTC time.
%
%      This example is equivalent to the one in cspice_edlimb, but it uses
%      cspice_inedpl to compute the limb.
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: inedpl_ex2.tm
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
%      function inedpl_ex2()
%
%         %
%         % Local parameters.
%         %
%         UTCSTR = '2016 Jul 14 19:45:00';
%
%         %
%         % Load the required kernels.
%         %
%         cspice_furnsh( 'inedpl_ex2.tm' );
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
%         scpos = tipm * scpos;
%
%         %
%         % Normalize the position to factors of the radii.
%         %
%         scpos = [ scpos(1)/rad(1)^2,                                     ...
%                   scpos(2)/rad(2)^2,                                     ...
%                   scpos(3)/rad(3)^2 ]';
%
%         %
%         % Find the apparent limb.  `limb' is a SPICE ellipse
%         % representing the limb.
%         %
%         [plane]       = cspice_nvc2pl( scpos, 1.0 );
%         [limb, found] = cspice_inedpl( rad(1), rad(2), rad(3), plane );
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
%   An ellipsoid and a plane can intersect in an ellipse, a single point, or
%   the empty set.
%
%-Exceptions
%
%   1)  If any of the lengths of the semi-axes of the input ellipsoid
%       are non-positive, the error SPICE(DEGENERATECASE) is signaled
%       by a routine in the call tree of this routine. `ellips' is not
%       modified. `found' is set to false.
%
%   2)  If the input plane in invalid, in other words, if the input
%       plane as the zero vector as its normal vector, the error
%       SPICE(INVALIDPLANE) is signaled by a routine in the call tree
%       of this routine. `ellips' is not modified. `found' is set to
%       false.
%
%   3)  If the input plane and ellipsoid are very nearly tangent,
%       roundoff error may cause this routine to give unreliable
%       results.
%
%   4)  If the input plane and ellipsoid are precisely tangent, the
%       intersection is a single point. In this case, the output
%       ellipse is degenerate, but `found' will still have the value
%       true. You must decide whether this output makes sense for
%       your application.
%
%   5)  If any of the input arguments, `a', `b', `c' or `plane', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `a', `b', `c' or `plane', is
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
%   PLANES.REQ
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
%   -Mice Version 1.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "ellipse" to "ellips".
%
%       Edited the header to comply with NAIF standard. Replaced
%       example with mathematical description of the algorithm in
%       cspice_edlimb, and added a second complete code example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 11-JUN-2013 (EDW)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 27-AUG-2012 (EDW)
%
%-Index_Entries
%
%   intersection of ellipsoid and plane
%
%-&

function [ellips, found] = cspice_inedpl( a, b, c, plane )

   switch nargin
      case 4

         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);
         plane = zzmice_pln(plane);

      otherwise

         error ( 'Usage: [ellips, found] = cspice_inedpl( a, b, c, plane )' )

   end

   try
      [ellips, found] = mice( 'inedpl_c', a, b, c, plane );

      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
