%-Abstract
%
%   CSPICE_ROTATE calculates the 3x3 rotation matrix generated
%   by a rotation of a specified angle about a specified axis.
%   This rotation operates as a rotation of the coordinate
%   system.
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
%      angle    the angle(s), given in radians, through which the rotation is
%               performed.
%
%               [1,n] = size(angle); double = class(angle)
%
%      iaxis    the index of the axis of rotation.
%
%               [1,1] = size(iaxis); int32 = class(iaxis)
%
%               The X, Y, and Z axes have indices 1, 2 and 3 respectively.
%
%   the call:
%
%      [mout] = cspice_rotate( angle, iaxis )
%
%   returns:
%
%      mout     the rotation matri(x|ces) which describes the rotation of a
%               reference frame through `angle' radians about the axis whose
%               index is `iaxis'.
%
%               If [1,1] = size(angle) then [3,3]   = size(mout)
%               If [1,n] = size(angle) then [3,3,n] = size(mout)
%                                            double = class(mout)
%
%               `mout' returns with the same vectorization measure, N,
%               as `angle'.
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
%   1) Compute the 3x3 matrix that rotates vectors from one
%      frame to another frame rotated by pi/10 radians about
%      +Y with respect to the first frame, and use it to transform
%      an arbitrary vector from the first frame to the second frame.
%
%      Example code begins here.
%
%
%      function rotate_ex1()
%
%         %
%         % Let's pick an arbitrary vector.
%         %
%         vec1 = [ 0.2; 0.04; 1.0 ];
%         fprintf( 'Vector in base frame:\n' )
%         fprintf( '  %16.12f  %16.12f  %16.12f\n', vec1 );
%
%         %
%         % Compute Pi/10 frame rotation about the Y axis.
%         %
%         rotmat = cspice_rotate( 0.1*cspice_pi, 2 );
%
%         %
%         % Apply the coordinate rotation to the vector.
%         %
%         vec2 = rotmat * vec1;
%         fprintf( 'Vector in rotated frame:\n' )
%         fprintf( '  %16.12f  %16.12f  %16.12f\n', vec2 );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Vector in base frame:
%          0.200000000000    0.040000000000    1.000000000000
%      Vector in rotated frame:
%         -0.118805691116    0.040000000000    1.012859915170
%
%
%-Particulars
%
%   A rotation about the first, i.e. x-axis, is described by
%
%      |  1        0          0      |
%      |  0   cos(theta) sin(theta)  |
%      |  0  -sin(theta) cos(theta)  |
%
%   A rotation about the second, i.e. y-axis, is described by
%
%      |  cos(theta)  0  -sin(theta)  |
%      |      0       1        0      |
%      |  sin(theta)  0   cos(theta)  |
%
%   A rotation about the third, i.e. z-axis, is described by
%
%      |  cos(theta) sin(theta)   0   |
%      | -sin(theta) cos(theta)   0   |
%      |       0          0       1   |
%
%   cspice_rotate decides which form is appropriate according to the value
%   of `iaxis'.
%
%-Exceptions
%
%   1)  If the axis index is not in the range 1 to 3, it will be
%       treated the same as that integer 1, 2, or 3 that is congruent
%       to it mod 3.
%
%   2)  If any of the input arguments, `angle' or `iaxis', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `angle' or `iaxis', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%       example's problem statement and modified example code accordingly.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 10-MAR-2015 (EDW)
%
%      Edited -I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 10-JAN-2006 (EDW)
%
%-Index_Entries
%
%   generate a rotation matrix
%
%-&

function [mout] = cspice_rotate( angle, iaxis )

   switch nargin

      case 2

         angle = zzmice_dp(angle);
         iaxis = zzmice_int( iaxis );

      otherwise
         error ( 'Usage: [_mout(3,3)_] = cspice_rotate( _angle_, iaxis )' )

   end

   %
   % Call the MEX library.
   %
   try
      [mout] = mice('rotate_c', angle, iaxis );
   catch spiceerr
      rethrow(spiceerr)
   end


