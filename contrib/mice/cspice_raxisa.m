%-Abstract
%
%   CSPICE_RAXISA computes the axis of the rotation given by an input matrix
%   and the angle of the rotation about that axis.
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
%      matrix   rotation matrix that gives the transformation from
%               some frame "frame1" to another frame "frame2".
%
%               [3,3]   = size(matrix); double = class(matrix)
%
%   the call:
%
%      [axis, angle] = cspice_raxisa( matrix )
%
%   returns:
%
%      axis   the unit vector pointing along the axis of the rotation. In
%             other words, 'axis' is a unit eigenvector of the input matrix,
%             corresponding to the eigenvalue 1. If the input matrix is
%             the identity matrix, 'axis' will be the vector (0, 0, 1).
%             If the input rotation is a rotation by pi radians, both
%             'axis' and -'axis' may be regarded as the axis of the rotation.
%
%             [3,1] = size(axis); double = class(axis)
%
%      angle  the angle between 'v' and 'matrix'*'v' for any non-zero vector
%             'v' orthogonal to 'axis'. 'angle' is given in radians.
%             The angle returned will be in the range from 0 to pi radians.
%
%             [1,1] = size(angle); double = class(angle)
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
%   1) Given an axis and an angle of rotation about that axis,
%      determine the rotation matrix. Using this matrix as input,
%      compute the axis and angle of rotation, and verify that
%      the later are equivalent by subtracting the original matrix
%      and the one resulting from using the computed axis and angle
%      of rotation on the cspice_axisar call.
%
%      Example code begins here.
%
%
%      function raxisa_ex1()
%
%         %
%         % Define an axis and an angle for rotation.
%         %
%         axis = [ 1.; 2.; 3. ];
%         angle = .1 * cspice_twopi;
%
%         %
%         % Determine the rotation matrix.
%         %
%         rot_mat = cspice_axisar( axis, angle );
%
%         %
%         % Now calculate the rotation axis and angle based on the
%         % matrix as input.
%         %
%         [axout, angout] = cspice_raxisa( rot_mat );
%         fprintf( 'Axis : %12.8f %12.8f %12.8f\n', axout );
%         fprintf( 'Angle: %12.8f\n', angout );
%
%         %
%         % Now input the axout and angout to cspice_axisar to
%         % compare against the original rotation matrix rot_mat.
%         %
%         rot_out = cspice_axisar( axout, angout );
%         diff = rot_mat - rot_out;
%         fprintf( '\nDifference between input and output matrices:\n');
%         fprintf('%20.16f %20.16f %20.16f\n', diff(1,:));
%         fprintf('%20.16f %20.16f %20.16f\n', diff(2,:));
%         fprintf('%20.16f %20.16f %20.16f\n', diff(3,:));
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Axis :   0.26726124   0.53452248   0.80178373
%      Angle:   0.62831853
%
%      Difference between input and output matrices:
%        0.0000000000000001   0.0000000000000000   0.0000000000000000
%       -0.0000000000000001   0.0000000000000001   0.0000000000000000
%        0.0000000000000000  -0.0000000000000001   0.0000000000000000
%
%
%      Note, the zero matrix is accurate to round-off error. A numerical
%      demonstration of equality.
%
%
%   2) This routine can be used to numerically approximate the
%      instantaneous angular velocity vector of a rotating object.
%
%      Suppose that R(t) is the rotation matrix whose columns
%      represent the inertial pointing vectors of the body-fixed axes
%      of an object at time t.
%
%      Then the angular velocity vector points along the vector given
%      by:
%
%                              T
%          limit  axis( R(t+h)R )
%          h-->0
%
%      And the magnitude of the angular velocity at time t is given
%      by:
%
%                             T
%         d angle ( R(t+h)R(t) )
%         ----------------------   at   h = 0
%                   dh
%
%      This code example computes the instantaneous angular velocity
%      vector of the Earth at 2000 Jan 01 12:00:00 TDB.
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model and orientation data for the Earth.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function raxisa_ex2()
%
%         %
%         % Load a PCK file containing a triaxial
%         % ellipsoidal shape model and orientation
%         % data for the Earth.
%         %
%         cspice_furnsh( 'pck00010.tpc' );
%
%         %
%         % Load time into the double precision variable `t'
%         % and the delta time (1 ms) into the double precision
%         % variable TH
%         %
%         t = 0.0;
%         h = 1e-3;
%
%         %
%         % Get the rotation matrices from IAU_EARTH to J2000
%         % at `t' and TH.
%         %
%         [rt]  = cspice_pxform( 'IAU_EARTH', 'J2000', t );
%         [rth] = cspice_pxform( 'IAU_EARTH', 'J2000', t+h );
%
%         %
%         % Compute the infinitesimal rotation r(t+h)r(t)^T
%         %
%         infrot = rth * rt.';
%
%         %
%         % Compute the `axis' and `angle' of the infinitesimal rotation
%         %
%         [axis, angle] = cspice_raxisa( infrot );
%
%         %
%         % Scale `axis' to get the angular velocity vector
%         %
%         angvel = angle/h * axis;
%
%         %
%         % Output the results.
%         %
%         fprintf( 'Instantaneous angular velocity vector:\n' )
%         fprintf( '%15.10f %14.10f %14.10f\n',                            ...
%                  angvel(1), angvel(2), angvel(3) )
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
%      Instantaneous angular velocity vector:
%         0.0000000000   0.0000000000   0.0000729212
%
%
%-Particulars
%
%   Every rotation matrix has an axis `a' such any vector `v'
%   parallel to that axis satisfies the equation
%
%      v = matrix * v
%
%   This routine returns a unit vector `axis' parallel to the axis of
%   the input rotation matrix. Moreover for any vector `w' orthogonal
%   to the axis of the rotation, the two vectors
%
%       axis,
%       w x (matrix*w)
%
%      (where "x" denotes the cross product operation)
%
%   will be positive scalar multiples of one another (at least
%   to within the ability to make such computations with double
%   precision arithmetic, and under the assumption that `matrix'
%   does not represent a rotation by zero or pi radians).
%
%   The angle returned will be the angle between `w' and matrix*w
%   for any vector orthogonal to `axis'.
%
%   If the input matrix is a rotation by 0 or pi radians some
%   choice must be made for the axis returned. In the case of
%   a rotation by 0 radians, `axis' is along the positive z-axis.
%   In the case of a rotation by 180 degrees, two choices are
%   possible. The choice made this routine is unspecified.
%
%-Exceptions
%
%   1)  If the input matrix is not a rotation matrix (where a fairly
%       loose tolerance is used to check this), an error is signaled
%       by a routine in the call tree of this routine.
%
%   2)  If the input matrix is the identity matrix, this routine
%       returns an angle of 0.0, and an axis of ( 0.0, 0.0, 1.0 ).
%
%   3)  If the input argument `matrix' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `matrix' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  If the input matrix is not a rotation matrix but is close
%       enough to pass the tests this routine performs on it, no error
%       will be signaled, but the results may have poor accuracy.
%
%   2)  The input matrix is taken to be an object that acts on
%       (rotates) vectors---it is not regarded as a coordinate
%       transformation. To find the axis and angle of a coordinate
%       transformation, input the transpose of that matrix to this
%       routine.
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
%       Edited the -Examples section to comply with NAIF standard.
%       Reformatted example's output, added problem statement and
%       a second example.
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
%   -Mice Version 1.0.1, 09-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 29-NOV-2005 (EDW)
%
%-Index_Entries
%
%   rotation axis of a matrix
%
%-&

function [axis, angle] = cspice_raxisa(matrix)

   switch nargin
      case 1

         matrix = zzmice_dp(matrix);

      otherwise

         error ( [ 'Usage: [ axis(3), angle] = ' ...
                   'cspice_raxisa( matrix(3,3) )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [axis, angle] = mice('raxisa_c', matrix);
   catch spiceerr
      rethrow(spiceerr)
   end


