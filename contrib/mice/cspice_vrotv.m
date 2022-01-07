%-Abstract
%
%   CSPICE_VROTV rotates a double precision 3-vector about a specified
%   axis vector by a specified angle (measured in radians) then
%   returns the rotated vector.
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
%      v        a vector to rotate.
%
%               [3,1] = size(v1); double = class(v1)
%
%      axis     a vector defining the axis about which to rotate `v'.
%
%               [3,1] = size(axis); double = class(axis)
%
%      theta    the value of the angle measured in radians through which
%               which rotate `v' about `axis'.
%
%               [1,1] = size(theta); double = class(theta)
%
%   the call:
%
%      r = cspice_vrotv( v, axis, theta )
%
%   returns:
%
%      r        the vector result of rotating `v' about `axis' through an
%               angle of `theta'.
%
%               [3,1] = size(r); double = class(r)
%
%               If axis is the zero vector, r = v.
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
%   1) Given an axis and an angle of rotation, perform rotations on
%      various vectors.
%
%      Example code begins here.
%
%
%      function vrotv_ex1()
%
%         %
%         % Given an axis of rotation and angle of rotation.
%         %
%         axis  = [ 0.; 0.; 1.];
%         theta = cspice_halfpi;
%
%         %
%         % Perform rotations on various vectors...
%         %
%         v1 = [ 1.; 2.; 3. ];
%         r1 = cspice_vrotv( v1, axis, theta );
%         fprintf( 'Input vector  : %12.6f %12.6f %12.6f\n',   v1 );
%         fprintf( 'Rotated vector: %12.6f %12.6f %12.6f\n\n', r1 );
%
%         v2 = [ 1.; 0.; 0. ];
%         r2 = cspice_vrotv( v2, axis, theta );
%         fprintf( 'Input vector  : %12.6f %12.6f %12.6f\n',   v2 );
%         fprintf( 'Rotated vector: %12.6f %12.6f %12.6f\n\n', r2 );
%
%         v3 = [ 0.; 1.; 0. ];
%         r3 = cspice_vrotv( v3, axis, theta );
%         fprintf( 'Input vector  : %12.6f %12.6f %12.6f\n', v3 );
%         fprintf( 'Rotated vector: %12.6f %12.6f %12.6f\n', r3 );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Input vector  :     1.000000     2.000000     3.000000
%      Rotated vector:    -2.000000     1.000000     3.000000
%
%      Input vector  :     1.000000     0.000000     0.000000
%      Rotated vector:     0.000000     1.000000     0.000000
%
%      Input vector  :     0.000000     1.000000     0.000000
%      Rotated vector:    -1.000000     0.000000     0.000000
%
%
%-Particulars
%
%   This routine computes the result of rotating (in a right handed
%   sense) the vector `v' about the axis represented by `axis' through
%   an angle of `theta' radians.
%
%   If `w' is a unit vector parallel to axis, then `r' is given by:
%
%       r = v + ( 1 - cos(theta) ) (w X(w X v)) + sin(theta) (w X v)
%
%   where "X" above denotes the vector cross product.
%
%-Exceptions
%
%   1)  If the input axis is the zero vector, `r' will be returned
%       as `v'.
%
%   2)  If any of the input arguments, `v', `axis' or `theta', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `v', `axis' or `theta', is not
%       of the expected type, or it does not have the expected
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
%   ROTATION.REQ
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
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and reformatted example's output.
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
%   -Mice Version 1.0.2, 18-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 14-JUL-2010 (EDW)
%
%       Corrected minor typo in header.
%
%   -Mice Version 1.0.0, 17-APR-2008 (EDW)
%
%-Index_Entries
%
%   vector rotation about an axis
%
%-&

function [r] = cspice_vrotv( v, axis, theta)

   switch nargin
      case 3

         v     = zzmice_dp(v);
         axis  = zzmice_dp(axis);
         theta = zzmice_dp(theta);

      otherwise

         error ( ['Usage: [r(3)] = cspice_vrotv( v(3), axis(3), theta)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [r] = mice('vrotv_c', v, axis, theta);
   catch spiceerr
      rethrow(spiceerr)
   end


