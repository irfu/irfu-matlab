%-Abstract
%
%   CSPICE_AXISAR constructs a rotation matrix that rotates vectors by a
%   specified angle about a specified axis.
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
%      axis,
%      angle    respectively, a rotation axis and a rotation angle.
%
%               [3,1] = size(axis); double = class(axis)
%               [1,1] = size(angle); double = class(angle)
%
%               `axis' and `angle' determine a coordinate transformation
%               whose effect on any vector V is to rotate V by `angle'
%               radians about the vector `axis'.
%
%   the call:
%
%      [r] = cspice_axisar( axis, angle )
%
%   returns:
%
%      r        a rotation matrix representing the coordinate transformation
%               determined by `axis' and `angle': for each vector `v', r*v is
%               the vector resulting from rotating `v' by `angle' radians
%               about `axis'.
%
%               [3,3] = size(r); double = class(r)
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
%   1) Compute a matrix that rotates vectors by pi/2 radians about
%      the Z-axis, and compute the rotation axis and angle based on
%      that matrix.
%
%
%      Example code begins here.
%
%
%      function axisar_ex1()
%
%         %
%         % Define an axis and an angle for rotation.
%         %
%         axis  = [ 0.0, 0.0, 1.0 ]';
%         angle = cspice_halfpi;
%
%         %
%         % Determine the rotation matrix.
%         %
%         [rotmat] = cspice_axisar( axis, angle );
%
%         %
%         % Now calculate the rotation axis and angle based on
%         % `rotmat' as input.
%         %
%         [axout, angout] = cspice_raxisa( rotmat );
%
%         %
%         % Display the results.
%         %
%         fprintf( 'Rotation matrix:\n' );
%         fprintf( '\n' );
%         fprintf( '%10.5f %9.5f %9.5f\n', rotmat' );
%         fprintf( '\n' );
%         fprintf( 'Rotation axis       : %9.5f %9.5f %9.5f\n', axout );
%         fprintf( 'Rotation angle (deg): %9.5f\n', angout * cspice_dpr );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Rotation matrix:
%
%         0.00000  -1.00000   0.00000
%         1.00000   0.00000   0.00000
%         0.00000   0.00000   1.00000
%
%      Rotation axis       :   0.00000   0.00000   1.00000
%      Rotation angle (deg):  90.00000
%
%
%   2) Linear interpolation between two rotation matrices.
%
%      Let r(t) be a time-varying rotation matrix; `r' could be
%      a C-matrix describing the orientation of a spacecraft
%      structure. Given two points in time t1 and t2 at which
%      r(t) is known, and given a third time t3, where
%
%         t1  <  t3  <  t2,
%
%      we can estimate r(t3) by linear interpolation. In other
%      words, we approximate the motion of `r' by pretending that
%      `r' rotates about a fixed axis at a uniform angular rate
%      during the time interval [t1, t2]. More specifically, we
%      assume that each column vector of `r' rotates in this
%      fashion. This procedure will not work if `r' rotates through
%      an angle of pi radians or more during the time interval
%      [t1, t2]; an aliasing effect would occur in that case.
%
%
%      Example code begins here.
%
%
%      function axisar_ex2()
%
%         %
%         % Lets assume that r(t) is the matrix that rotates
%         % vectors by pi/2 radians about the Z-axis every
%         % minute.
%         %
%         % Let
%         %
%         %    r1 = r(t0 - 1), for t1 =  0", and
%         %    r2 = r(t2), for t1 = 60".
%         %
%         % Define both matrices and times.
%         %
%         axis = [ 0.0, 0.0, 1.0 ]';
%
%         t1   =  0.0;
%         t2   = 60.0;
%         t3   = 30.0;
%
%         [r1] = eye(3);
%         [r2] = cspice_axisar( axis, cspice_halfpi );
%
%         q = r2 * transpose( r1 );
%         [axis, angle] = cspice_raxisa( q );
%
%         %
%         % Find the fraction of the total rotation angle that `r'
%         % rotates through in the time interval [t1, t3].
%         %
%         frac = ( t3 - t1 )  /  ( t2 - t1 );
%
%         %
%         % Finally, find the rotation `delta' that r(t) undergoes
%         % during the time interval [t1, t3], and apply that
%         % rotation to `r1', yielding r(t3), which we'll call `r3'.
%         %
%         [delta] = cspice_axisar( axis, frac * angle );
%         r3      = delta * r1;
%
%         %
%         % Display the results.
%         %
%         fprintf( 'Time (s)            : %9.5f\n', t3 );
%         fprintf( 'Rotation axis       : %9.5f %9.5f %9.5f\n', axis );
%         fprintf( 'Rotation angle (deg): %9.5f\n',                         ...
%                                           frac * angle * cspice_dpr  );
%         fprintf( 'Rotation matrix     :\n' );
%         fprintf( '\n' );
%         fprintf( '%10.5f %9.5f %9.5f\n', r3' );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Time (s)            :  30.00000
%      Rotation axis       :   0.00000   0.00000   1.00000
%      Rotation angle (deg):  45.00000
%      Rotation matrix     :
%
%         0.70711  -0.70711   0.00000
%         0.70711   0.70711   0.00000
%         0.00000   0.00000   1.00000
%
%
%-Particulars
%
%   cspice_axisar can be thought of as a partial inverse of cspice_raxisa.
%   cspice_axisar really is a `left inverse': the code fragment
%
%      [axis, angle] = cspice_raxisa( r );
%      [r]           = cspice_axisar( axis, angle );
%
%   preserves `r', except for round-off error, as long as `r' is a
%   rotation matrix.
%
%   On the other hand, the code fragment
%
%      [r]           = cspice_axisar( axis, angle );
%      [axis, angle] = cspice_raxisa( r );
%
%   preserves `axis' and `angle', except for round-off error, only if
%   `angle' is in the range (0, pi). So cspice_axisar is a right inverse
%   of cspice_raxisa only over a limited domain.
%
%-Exceptions
%
%   1)  If `axis' is the zero vector, the rotation generated is the
%       identity. This is consistent with the specification of cspice_vrotv.
%
%   2)  If any of the input arguments, `axis' or `angle', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `axis' or `angle', is not of
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Updated code
%       example #1 and added second example.
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
%   -Mice Version 1.0.1, 28-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 29-NOV-2005 (EDW)
%
%-Index_Entries
%
%   axis and angle to rotation
%
%-&

function [r] = cspice_axisar( axis, angle)

   switch nargin
      case 2

         axis  = zzmice_dp(axis);
         angle = zzmice_dp(angle);

      otherwise

         error ( 'Usage: [r(3,3)] = cspice_axisar( axis(3), angle)' )

   end

   %
   % Call the MEX library.
   %
   try
      [r] = mice( 'axisar_c', axis, angle );
   catch spiceerr
      rethrow(spiceerr)
   end


