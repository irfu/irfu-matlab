%-Abstract
%
%   CSPICE_M2Q calculates a unit quaternion corresponding to a
%   specified rotation matrix.
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
%      r        a rotation matri(x|ces).
%
%               [3,3]   = size(r); double = class(r)
%
%                  or
%
%               [3,3,n] = size(r); double = class(r)
%
%   the call:
%
%      [q] = cspice_m2q( r )
%
%   returns:
%
%      q        an array of unit-length SPICE-style quaternion(s)
%               representing `r'.
%
%               If [3,3]   = size(r) then [4,1] = size(q)
%               If [3,3,n] = size(r) then [4,n] = size(q)
%                                        double = class(q)
%
%               See the discussion of quaternion styles in -Particulars
%               below.
%
%               `q' is a 4-dimensional vector. If `r' rotates vectors
%               in the counterclockwise sense by an angle of theta
%               radians about a unit vector `a', where
%
%                  0 < theta < pi
%                    -       -
%
%               then letting h = theta/2,
%
%                  q = ( cos(h), sin(h)a ,  sin(h)a ,  sin(h)a ).
%                                       1          2          3
%
%               The restriction that theta must be in the range
%               [0, pi] determines the output quaternion `q'
%               uniquely except when theta = pi; in this special
%               case, both of the quaternions
%
%                  q = ( 0,  a ,  a ,  a  )
%                             1    2    3
%               and
%
%                  q = ( 0, -a , -a , -a  )
%                             1    2    3
%
%               are possible outputs.
%
%               `q' returns with the same vectorization measure, N,
%               as `r' .
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
%   1) Create a 3-dimensional rotation matrix of 90 degrees about the
%      Z axis and convert it to a unit quaternion. Verify that the
%      norm of the quaternion is equal to 1.
%
%      Example code begins here.
%
%
%      function m2q_ex1()
%
%         %
%         % Create a rotation matrix of 90 degrees about the Z axis.
%         %
%         r = cspice_rotate( cspice_halfpi, 3);
%
%         fprintf('Rotation matrix:\n');
%         fprintf('%15.7f %15.7f %15.7f\n', r(1,:));
%         fprintf('%15.7f %15.7f %15.7f\n', r(2,:));
%         fprintf('%15.7f %15.7f %15.7f\n', r(3,:));
%
%         %
%         % Convert the matrix to a quaternion.
%         %
%         q = cspice_m2q( r );
%         fprintf('Unit quaternion:\n');
%         fprintf('%15.7f %15.7f %15.7f %15.7f\n', q);
%
%         %            _
%         % Confirm || q || = 1.
%         %
%         fprintf( '\n|| q || = %15.7f\n', q'  * q );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Rotation matrix:
%            0.0000000       1.0000000       0.0000000
%           -1.0000000       0.0000000       0.0000000
%            0.0000000       0.0000000       1.0000000
%      Unit quaternion:
%            0.7071068       0.0000000       0.0000000      -0.7071068
%
%      || q || =       1.0000000
%
%
%      Note, the call sequence:
%
%         q = cspice_m2q( r );
%         r = cspice_q2m( q );
%
%      preserves `r' except for round-off error. Yet, the call sequence:
%
%         r = cspice_q2m( q );
%         q = cspice_m2q( r );
%
%      may preserve `q' or return `-q'.
%
%-Particulars
%
%   About SPICE quaternions
%   =======================
%
%   There are (at least) two popular "styles" of quaternions; these
%   differ in the layout of the quaternion elements, the definition
%   of the multiplication operation, and the mapping between the set
%   of unit quaternions and corresponding rotation matrices.
%
%   SPICE-style quaternions have the scalar part in the first
%   component and the vector part in the subsequent components. The
%   SPICE convention, along with the multiplication rules for SPICE
%   quaternions, are those used by William Rowan Hamilton, the
%   inventor of quaternions.
%
%   Another common quaternion style places the scalar component
%   last. This style is often used in engineering applications.
%
%-Exceptions
%
%   1)  If `r' is not a rotation matrix, the error SPICE(NOTAROTATION)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If the input argument `r' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   3)  If the input argument `r' is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 09-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 10-JAN-2006 (EDW)
%
%-Index_Entries
%
%   matrix to quaternion
%
%-&

function [q] = cspice_m2q(r)

   switch nargin
      case 1

         r = zzmice_dp(r);

      otherwise

         error ( 'Usage: [_q(4)_] = cspice_m2q( _r(3,3)_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      [q] = mice('m2q_c', r);
   catch spiceerr
      rethrow(spiceerr)
   end



