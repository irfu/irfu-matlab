%-Abstract
%
%   CSPICE_EDNMPT returns the unique point on an ellipsoid's surface where
%   the outward normal direction is a given vector.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      a        the length of the semi-axis of the ellipsoid that is parallel
%               to the X-axis of the body-fixed coordinate system.
%
%               [1,1] = size(a); double = class(a)
%
%      b        the length of the semi-axis of the ellipsoid that is parallel
%               to the Y-axis of the body-fixed coordinate system.
%
%               [1,1] = size(b); double = class(b)
%
%      c        the length of the semi-axis of the ellipsoid that is parallel
%               to the Z-axis of the body-fixed coordinate system.
%
%               [1,1] = size(c); double = class(c)
%
%      normal   a non-zero vector.
%
%               [3,1] = size(normal); double = class(normal)
%
%               The unique point on the ellipsoid at which `normal' is an
%               outward normal vector is sought.
%
%   the call:
%
%      [point] = cspice_ednmpt( a, b, c, normal )
%
%   returns:
%
%      point    the unique point on the ellipsoid at which `normal' is an
%               outward normal vector.
%
%               [3,1] = size(point); double = class(point)
%
%               `point' is a 3-vector giving the body-fixed coordinates
%               of a point on the ellipsoid. In body-fixed coordinates,
%               the semi-axes of the ellipsoid are aligned with the X,
%               Y, and Z-axes of the coordinate system.
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
%   1) Choose a triaxial ellipsoid with three unequal semi-axis
%      lengths. Pick several vectors; find the points on the
%      ellipsoid where the respective outward normals are parallel to
%      those vectors.
%
%      Check the results: at each point, a computed outward normal
%      vector should have very small angular separation from the
%      input vector. Also, the point should be on the surface of the
%      ellipsoid. The ellipsoid can be thought of as a level surface
%      of the function
%
%                           2        2         2
%         f(x, y, z) = (x/a)  + (y/b)  +  (z/c)
%
%      where `a', `b', `c' are the semi-axis lengths of the ellipsoid.
%      Specifically, the ellipsoid is the set
%
%         { (x, y, z) : f(x, y, z)  =  1 }
%
%      We can evaluate F at a point to determine whether that point
%      is close to the ellipsoid's surface.
%
%
%      Example code begins here.
%
%
%      function ednmpt_ex1()
%
%         %
%         % Initialize the ellipsoid semi-axes.
%         %
%         a = 10.0;
%         b =  5.0;
%         c =  2.0;
%
%         %
%         % Pick several vectors; find the points
%         % on the ellipsoid where the respective
%         % outward normals are parallel to those
%         % vectors; check the results.
%         %
%         xnorml   = [ 0.0, 0.0, 3.0 ]';
%         [point]  = cspice_ednmpt( a, b, c, xnorml );
%         [normal] = cspice_surfnm( a, b, c, point );
%
%         fprintf( ' \n' )
%         fprintf( 'Semi-axis lengths:    %13.8f %13.8f %13.8f\n', a, b, c )
%         fprintf( 'Input vector:         %13.8f %13.8f %13.8f\n',         ...
%                                  xnorml(1), xnorml(2), xnorml(3) )
%         fprintf( 'Point:                %13.8f %13.8f %13.8f\n',         ...
%                                     point(1), point(2), point(3) )
%         fprintf( 'Outward normal:       %13.8f %13.8f %13.8f\n',         ...
%                                  normal(1), normal(2), normal(3) )
%         fprintf( 'Angular error (rad):  %13.8f\n',                       ...
%                      cspice_vsep( normal, xnorml ) )
%         fprintf( 'Off-surface error:    %13.8f\n',                       ...
%                  (point(1)/a) ^ 2 + (point(2)/b) ^ 2 +                   ...
%                  (point(3)/c) ^ 2 - 1                 )
%         fprintf( ' \n' )
%
%         xnorml   = [ 15.0, -7.0, 3.0 ]';
%         [point]  = cspice_ednmpt( a, b, c, xnorml );
%         [normal] = cspice_surfnm( a, b, c, point );
%
%         fprintf( 'Semi-axis lengths:    %13.8f %13.8f %13.8f\n', a, b, c )
%         fprintf( 'Input vector:         %13.8f %13.8f %13.8f\n',         ...
%                                  xnorml(1), xnorml(2), xnorml(3) )
%         fprintf( 'Point:                %13.8f %13.8f %13.8f\n',         ...
%                                     point(1), point(2), point(3) )
%         fprintf( 'Outward normal:       %13.8f %13.8f %13.8f\n',         ...
%                                  normal(1), normal(2), normal(3) )
%         fprintf( 'Angular error (rad):  %13.8f\n',                       ...
%                      cspice_vsep( normal, xnorml ) )
%         fprintf( 'Off-surface error:    %13.8f\n',                       ...
%                  (point(1)/a) ^ 2 + (point(2)/b) ^ 2 +                   ...
%                  (point(3)/c) ^ 2 - 1                 )
%         fprintf( ' \n' )
%
%         xnorml   = [ 15.0, -7.0, 3.0 ]';
%         [point]  = cspice_ednmpt( a, b, c, xnorml );
%         [normal] = cspice_surfnm( a, b, c, point );
%
%         fprintf( 'Semi-axis lengths:    %13.8f %13.8f %13.8f\n', a, b, c )
%         fprintf( 'Input vector:         %13.8f %13.8f %13.8f\n',         ...
%                                  xnorml(1), xnorml(2), xnorml(3) )
%         fprintf( 'Point:                %13.8f %13.8f %13.8f\n',         ...
%                                     point(1), point(2), point(3) )
%         fprintf( 'Outward normal:       %13.8f %13.8f %13.8f\n',         ...
%                                  normal(1), normal(2), normal(3) )
%         fprintf( 'Angular error (rad):  %13.8f\n',                       ...
%                      cspice_vsep( normal, xnorml ) )
%         fprintf( 'Off-surface error:    %13.8f\n',                       ...
%                  (point(1)/a) ^ 2 + (point(2)/b) ^ 2 +                   ...
%                  (point(3)/c) ^ 2 - 1                 )
%         fprintf( ' \n' )
%
%         xnorml   = [ a/2, b/2, c/2 ]';
%         [point]  = cspice_ednmpt( a, b, c, xnorml );
%         [normal] = cspice_surfnm( a, b, c, point );
%
%         fprintf( 'Semi-axis lengths:    %13.8f %13.8f %13.8f\n', a, b, c )
%         fprintf( 'Input vector:         %13.8f %13.8f %13.8f\n',         ...
%                                  xnorml(1), xnorml(2), xnorml(3) )
%         fprintf( 'Point:                %13.8f %13.8f %13.8f\n',         ...
%                                     point(1), point(2), point(3) )
%         fprintf( 'Outward normal:       %13.8f %13.8f %13.8f\n',         ...
%                                  normal(1), normal(2), normal(3) )
%         fprintf( 'Angular error (rad):  %13.8f\n',                       ...
%                      cspice_vsep( normal, xnorml ) )
%         fprintf( 'Off-surface error:    %13.8f\n',                       ...
%                  (point(1)/a) ^ 2 + (point(2)/b) ^ 2 +                   ...
%                  (point(3)/c) ^ 2 - 1                 )
%         fprintf( ' \n' )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Semi-axis lengths:      10.00000000    5.00000000    2.00000000
%      Input vector:            0.00000000    0.00000000    3.00000000
%      Point:                   0.00000000    0.00000000    2.00000000
%      Outward normal:          0.00000000    0.00000000    1.00000000
%      Angular error (rad):     0.00000000
%      Off-surface error:       0.00000000
%
%      Semi-axis lengths:      10.00000000    5.00000000    2.00000000
%      Input vector:           15.00000000   -7.00000000    3.00000000
%      Point:                   9.73103203   -1.13528707    0.07784826
%      Outward normal:          0.89165745   -0.41610681    0.17833149
%      Angular error (rad):     0.00000000
%      Off-surface error:       0.00000000
%
%      Semi-axis lengths:      10.00000000    5.00000000    2.00000000
%      Input vector:           15.00000000   -7.00000000    3.00000000
%      Point:                   9.73103203   -1.13528707    0.07784826
%      Outward normal:          0.89165745   -0.41610681    0.17833149
%      Angular error (rad):     0.00000000
%      Off-surface error:       0.00000000
%
%      Semi-axis lengths:      10.00000000    5.00000000    2.00000000
%      Input vector:            5.00000000    2.50000000    1.00000000
%      Point:                   9.69412864    1.21176608    0.07755303
%      Outward normal:          0.88045091    0.44022545    0.17609018
%      Angular error (rad):     0.00000000
%      Off-surface error:       0.00000000
%
%
%-Particulars
%
%   This routine can be used to determine the distance between an
%   ellipsoid and a non-intersecting plane. This distance computation
%   supports computation of terminator points on an ellipsoid.
%
%-Exceptions
%
%   1)  If any of the semi-axis lengths is non-positive, the error
%       SPICE(BADAXISLENGTH) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If any of the semi-axis lengths underflows to zero when
%       divided by the largest semi-axis length, the error
%       SPICE(AXISUNDERFLOW) is signaled by a routine in the call tree
%       of this routine.
%
%   3)  If `normal' is the zero vector, the error SPICE(ZEROVECTOR)
%       is signaled by a routine in the call tree of this routine.
%
%   4)  If the input pass the above checks but lead to a
%       divide-by-zero error or to a computing an invalid argument of
%       a fractional exponential expression, the error
%       SPICE(DEGENERATECASE) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If any of the input arguments, `a', `b', `c' or `normal', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `a', `b', `c' or `normal', is
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
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   point on an ellipsoid having given surface normal
%
%-&
function [point] = cspice_ednmpt( a, b, c, normal )

   switch nargin
      case 4

         a = zzmice_dp(a);
         b = zzmice_dp(b);
         c = zzmice_dp(c);
         normal = zzmice_dp(normal);

      otherwise

         error ( 'Usage: [point(3)] = cspice_ednmpt( a, b, c, normal(3) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [point] = mice('ednmpt_c', a, b, c, normal);
   catch spiceerr
      rethrow(spiceerr)
   end
