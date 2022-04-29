%-Abstract
%
%   CSPICE_NPEDLN calculates the nearest point on a triaxial
%   ellipsoid to a specified line, and the distance from the
%   ellipsoid point to the line.
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
%      linept,
%      linedr   respectively, a point and a direction vector that define a
%               line.
%
%               [3,1] = size(linept); double = class(linept)
%               [3,1] = size(linedr); double = class(linedr)
%
%               The line is the set of vectors
%
%                     linept   +   t * linedr
%
%               where `t' is any real number.
%
%   the call:
%
%      [pnear, dist] = cspice_npedln( a, b, c, linept, linedr )
%
%   returns:
%
%      pnear   the point on the ellipsoid closest to the line, if the line
%              doesn't intersect the ellipsoid.
%
%              [3,1] = size(pnear); double = class(pnear)
%
%              If the line intersects the ellipsoid, `pnear' will be a point
%              of intersection. If `linept' is outside of the ellipsoid,
%              `pnear' will be the closest point of intersection. If `linept'
%              is inside ellipsoid, `pnear' will not necessarily be the
%              closest point of intersection.
%
%      dist    the distance of the line from the ellipsoid.
%
%              [1,1] = size(dist); double = class(dist)
%
%              This is the minimum distance between any point on the line
%              and any point on the ellipsoid.
%
%              If the line intersects the ellipsoid, `dist' is zero.
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
%   1) Find the distance between an instrument optic axis ray and
%      the surface of a body modeled as a tri-axial ellipsoid, given the
%      instrument position and pointing unit vector, and the body
%      semi-axes lengths.
%
%      Example code begins here.
%
%
%      function npedln_ex1()
%
%         %
%         % We can find the distance between an instrument optic axis ray
%         % and the surface of a body modeled as a tri-axial ellipsoid
%         % using this routine.  If the instrument position and pointing
%         % unit vector in body-fixed coordinates are:
%         %
%         linept = [ 1.0e6,  2.0e6,  3.0e6 ]';
%         linedr = [ -4.472091234e-1, -8.944182469e-1, -4.472091234e-3 ]';
%
%         %
%         % The body semi-axes lengths:
%         %
%         a = 7.0e5;
%         b = 7.0e5;
%         c = 6.0e5;
%
%         %
%         % The call to cspice_npedln yields a value for `pnear', the nearest
%         % point on the body to the optic axis ray and a value for `dist',
%         % the distance to the ray.
%         %
%         [pnear, dist] = cspice_npedln( a, b, c, linept, linedr );
%
%         fprintf( 'Nearest point: %16.7f %16.7f %16.7f\n', pnear );
%         fprintf( 'Distance     : %16.7f\n',               dist  );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Nearest point:    -1633.3110792    -3266.6222158   599991.8335001
%      Distance     :  2389967.9338300
%
%
%-Particulars
%
%   For any ellipsoid and line, if the line does not intersect the
%   ellipsoid, there is a unique point on the ellipsoid that is
%   closest to the line. Therefore, the distance dist between
%   ellipsoid and line is well-defined. The unique line segment of
%   length dist that connects the line and ellipsoid is normal to
%   both of these objects at its endpoints.
%
%   If the line intersects the ellipsoid, the distance between the
%   line and ellipsoid is zero.
%
%-Exceptions
%
%   If this routine detects an error, the output arguments `pnear' and
%   `dist' are not modified.
%
%   1)  If the length of any semi-axis of the ellipsoid is
%       non-positive, the error SPICE(INVALIDAXISLENGTH) is signaled
%       by a routine in the call tree of this routine.
%
%   2)  If the line's direction vector is the zero vector, the error
%       SPICE(ZEROVECTOR) is signaled by a routine in the call tree of
%       this routine.
%
%   3)  If the length of any semi-axis of the ellipsoid is zero after
%       the semi-axis lengths are scaled by the reciprocal of the
%       magnitude of the longest semi-axis and then squared, the error
%       SPICE(DEGENERATECASE) is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If the input ellipsoid is extremely flat or needle-shaped
%       and has its shortest axis close to perpendicular to the input
%       line, numerical problems could cause this routine's algorithm
%       to fail, in which case, the error SPICE(DEGENERATECASE) is
%       signaled by a routine in the call tree of this routine.
%
%   5)  If any of the input arguments, `a', `b', `c', `linept' or
%       `linedr', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   6)  If any of the input arguments, `a', `b', `c', `linept' or
%       `linedr', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Reformatted example's output and added problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Edits to header documentation -I/O. The previous documentation
%       showed the routine would accept vectorized arguments. All
%       -I/O arguments are expected as non-vectorized.
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
%   distance between line and ellipsoid
%   distance between line of sight and body
%   nearest point on ellipsoid to line
%
%-&

function [pnear, dist] = cspice_npedln( a, b, c, linept, linedr )

   switch nargin
      case 5

         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
         linept = zzmice_dp(linept);
         linedr = zzmice_dp(linedr);

      otherwise

         error ( ['Usage: [ pnear(3), dist] = ' ...
                  'cspice_npedln( a, b, c, linept(3), linedr(3) )'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [npedln] = mice( 'npedln_s',a, b, c, linept, linedr );
      pnear    = reshape( [npedln.pos], 3, [] );
      dist     = reshape( [npedln.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end


