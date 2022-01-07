%-Abstract
%
%   CSPICE_ROTVEC transforms a vector to a new coordinate system rotated by
%   `angle' radians about axis `iaxis'. This transformation rotates `v1' by
%   -angle radians about the specified axis.
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
%      v1       a vector (typically representing a vector fixed in inertial
%               space) which is to be expressed in another coordinate system.
%
%               [3,1] = size(v1); double = class(v1)
%
%               The vector remains fixed but the coordinate system changes.
%
%      angle    an angle given in radians, through which the rotation is
%               performed.
%
%               [1,1] = size(angle); double = class(angle)
%
%      iaxis    the index of the axis of rotation.
%
%               [1,1] = size(iaxis); int32 = class(iaxis)
%
%               The X, Y, and Z axes have indices 1, 2 and 3 respectively.
%
%   the call:
%
%      [vout] = cspice_rotvec( v1, angle, iaxis )
%
%   returns:
%
%      vout     the vector expressed in the new coordinate system specified
%               by the angle of rotation and axis.
%
%               [3,1] = size(vout); double = class(vout)
%
%               If
%
%                  m = [angle]
%                             iaxis
%
%               represents the rotation matrix described by the `angle'
%               and `iaxis', (refer to the routine cspice_rotate) then
%
%                  vout =  m * v1  = [angle]      * v1
%                                           iaxis
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
%   1) Apply a rotation of -45.0 degrees about the +Z axis to
%      a 3 dimensional vector.
%
%      Example code begins here.
%
%
%      function rotvec_ex1()
%
%         %
%         % Input values.
%         %
%         v1    = [1.414, 0.0, 0.0]';
%
%         angle = cspice_pi/4;
%         iaxis = 3;
%
%         %
%         % Rotate `v1' by `angle' radians about `iaxis'.
%         %
%         [vout] = cspice_rotvec( v1, angle, iaxis );
%
%         fprintf( 'Input vector  : %9.3f %9.3f %9.3f\n', ...
%                                     v1(1), v1(2), v1(3) )
%         fprintf( 'Rotated vector: %9.3f %9.3f %9.3f\n', ...
%                               vout(1), vout(2), vout(3) )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Input vector  :     1.414     0.000     0.000
%      Rotated vector:     1.000    -1.000     0.000
%
%
%-Particulars
%
%   A rotation about the first, i.e. X-axis, is described by
%
%      .-                              -.
%      |   1       0            0       |
%      |   0   cos(theta)   sin(theta)  |
%      |   0  -sin(theta)   cos(theta)  |
%      `-                              -'
%
%   A rotation about the second, i.e. Y-axis, is described by
%
%      .-                              -.
%      |   cos(theta)   0  -sin(theta)  |
%      |       0        1       0       |
%      |   sin(theta)   1   cos(theta)  |
%      `-                              -'
%
%   A rotation about the third, i.e. Z-axis, is described by
%
%      .-                              -.
%      |   cos(theta)   sin(theta)   0  |
%      |  -sin(theta)   cos(theta)   0  |
%      |       0            0        1  |
%      `-                              -'
%
%   cspice_rotvec decides which form is appropriate according to the value
%   of `iaxis' and applies the rotation to the input vector.
%
%-Exceptions
%
%   1)  If the `iaxis' index is not in the range 1 to 3, it will be
%       treated the same as that integer 1, 2, or 3 that is congruent
%       to it mod 3.
%
%   2)  If any of the input arguments, `v1', `angle' or `iaxis', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `v1', `angle' or `iaxis', is
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
%   rotate a vector
%
%-&
function [vout] = cspice_rotvec( v1, angle, iaxis )

   switch nargin
      case 3

         v1 = zzmice_dp(v1);
         angle = zzmice_dp(angle);
         iaxis = zzmice_int(iaxis);

      otherwise

         error ( 'Usage: [vout(3)] = cspice_rotvec( v1(3), angle, iaxis )' )

   end

   %
   % Call the MEX library.
   %
   try
      [vout] = mice('rotvec_c', v1, angle, iaxis);
   catch spiceerr
      rethrow(spiceerr)
   end
