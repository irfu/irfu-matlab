%-Abstract
%
%   CSPICE_Q2M calculates the rotation matrix corresponding to a
%   specified unit quaternion.
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
%      q        an array of unit-length SPICE-style quaternion(s).
%
%               [4,n] = size(q); double = class(q)
%
%               `q' has the property that
%
%                  || q ||  =  1
%
%               See the discussion of quaternion styles in
%               -Particulars below.
%
%   the call:
%
%      [r] = cspice_q2m( q )
%
%   returns:
%
%      r        the rotation matri(x|ces) representing the same rotation as
%               does `q'.
%
%               If [4,1] = size(q) then [3,3]   = size(r)
%               If [4,n] = size(q) then [3,3,n] = size(r)
%                                       double = class(r)
%
%               See the discussion titled "Associating SPICE Quaternions
%               with Rotation Matrices" in -Particulars below.
%
%               `r' returns with the same vectorization measure, N,
%               as `q'.
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
%   1) Define a unit quaternion, confirm that its norm is equal to 1.0
%      and convert it to a matrix form.
%
%      Example code begins here.
%
%
%      function q2m_ex1()
%
%         %
%         % Define a unit quaternion.
%         %
%         q = [ sqrt(2.0)/2.0, 0.0, 0.0, -sqrt(2.0)/2.0]';
%         fprintf( 'Quaternion : %12.8f %12.8f %12.8f %12.8f\n', q );
%
%         %
%         % Confirm q satisfies || q || = 1. Calculate q * q.
%         %
%         fprintf( 'Norm       : %12.8f\n', q' * q );
%
%         %
%         % Convert the quaternion to a matrix form.
%         %
%         m = cspice_q2m( q );
%         fprintf( 'Matrix form:\n')
%         fprintf('%15.7f %15.7f %15.7f\n', m(1,:));
%         fprintf('%15.7f %15.7f %15.7f\n', m(2,:));
%         fprintf('%15.7f %15.7f %15.7f\n', m(3,:));
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Quaternion :   0.70710678   0.00000000   0.00000000  -0.70710678
%      Norm       :   1.00000000
%      Matrix form:
%            0.0000000       1.0000000       0.0000000
%           -1.0000000       0.0000000      -0.0000000
%           -0.0000000       0.0000000       1.0000000
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
%   If a 4-dimensional vector `q' satisfies the equality
%
%      || q ||   =  1
%
%   or equivalently
%
%          2          2          2          2
%      q(1)   +   q(2)   +   q(3)   +   q(4)   =  1,
%
%   then we can always find a unit vector `a' and a scalar r such that
%
%      q = ( cos(r/2), sin(r/2)a(1), sin(r/2)a(2), sin(r/2)a(3) ).
%
%   We can interpret `a' and r as the axis and rotation angle of a
%   rotation in 3-space. If we restrict r to the range [0, pi],
%   then r and `a' are uniquely determined, except if r = pi. In this
%   special case, A and -an are both valid rotation axes.
%
%   Every rotation is represented by a unique orthogonal matrix; this
%   routine returns that unique rotation matrix corresponding to `q'.
%
%   The Mice routine cspice_m2q is an one-sided inverse of this routine:
%   given any rotation matrix `r', the calls
%
%      [q] = cspice_m2q( r );
%      [r] = cspice_q2m( q );
%
%   leave `r' unchanged, except for round-off error. However, the
%   calls
%
%      [r] = cspice_q2m( q );
%      [q] = cspice_m2q( r );
%
%   might preserve `q' or convert `q' to -q.
%
%
%   Quaternion Styles
%   -----------------
%
%   There are different "styles" of quaternions used in
%   science and engineering applications. Quaternion styles
%   are characterized by
%
%   -  The order of quaternion elements
%
%   -  The quaternion multiplication formula
%
%   -  The convention for associating quaternions
%      with rotation matrices
%
%   Two of the commonly used styles are
%
%      - "SPICE"
%
%         > Invented by Sir William Rowan Hamilton
%         > Frequently used in mathematics and physics textbooks
%
%      - "Engineering"
%
%         > Widely used in aerospace engineering applications
%
%
%   Mice routine interfaces ALWAYS use SPICE quaternions.
%   Quaternions of any other style must be converted to SPICE
%   quaternions before they are passed to SPICELIB routines.
%
%
%   Relationship between SPICE and Engineering Quaternions
%   ------------------------------------------------------
%
%   Let `m' be a rotation matrix such that for any vector `v',
%
%      m*v
%
%   is the result of rotating `v' by theta radians in the
%   counterclockwise direction about unit rotation axis vector `a'.
%   Then the SPICE quaternions representing `m' are
%
%      (+/-) (  cos(theta/2),
%               sin(theta/2) a(1),
%               sin(theta/2) a(2),
%               sin(theta/2) a(3)  )
%
%   while the engineering quaternions representing `m' are
%
%      (+/-) ( -sin(theta/2) a(1),
%              -sin(theta/2) a(2),
%              -sin(theta/2) a(3),
%               cos(theta/2)       )
%
%   For both styles of quaternions, if a quaternion q represents
%   a rotation matrix `m', then -q represents `m' as well.
%
%   Given an engineering quaternion
%
%      qeng   = ( q0,  q1,  q2,  q3 )
%
%   the equivalent SPICE quaternion is
%
%      qspice = ( q3, -q0, -q1, -q2 )
%
%
%   Associating SPICE Quaternions with Rotation Matrices
%   ----------------------------------------------------
%
%   Let `from' and `to' be two right-handed reference frames, for
%   example, an inertial frame and a spacecraft-fixed frame. Let the
%   symbols
%
%      v    ,   v
%       from     to
%
%   denote, respectively, an arbitrary vector expressed relative to
%   the `from' and `to' frames. Let `m' denote the transformation matrix
%   that transforms vectors from frame `from' to frame `to'; then
%
%      v   =  m * v
%       to         from
%
%   where the expression on the right hand side represents left
%   multiplication of the vector by the matrix.
%
%   Then if the unit-length SPICE quaternion q represents `m', where
%
%      q = (q0, q1, q2, q3)
%
%   the elements of `m' are derived from the elements of q as follows:
%
%        .-                                                         -.
%        |           2    2                                          |
%        | 1 - 2*( q2 + q3 )   2*(q1*q2 - q0*q3)   2*(q1*q3 + q0*q2) |
%        |                                                           |
%        |                                                           |
%        |                               2    2                      |
%    m = | 2*(q1*q2 + q0*q3)   1 - 2*( q1 + q3 )   2*(q2*q3 - q0*q1) |
%        |                                                           |
%        |                                                           |
%        |                                                   2    2  |
%        | 2*(q1*q3 - q0*q2)   2*(q2*q3 + q0*q1)   1 - 2*( q1 + q2 ) |
%        |                                                           |
%        `-                                                         -.
%
%   Note that substituting the elements of -q for those of q in the
%   right hand side leaves each element of `m' unchanged; this shows
%   that if a quaternion q represents a matrix `m', then so does the
%   quaternion -q.
%
%   To map the rotation matrix `m' to a unit quaternion, we start by
%   decomposing the rotation matrix as a sum of symmetric
%   and skew-symmetric parts:
%
%                                      2
%      m = [ i  +  (1-cos(theta)) omega  ] + [ sin(theta) omega ]
%
%                   symmetric                   skew-symmetric
%
%
%   `omega' is a skew-symmetric matrix of the form
%
%                 .-             -.
%                 |  0   -n3   n2 |
%                 |               |
%       omega  =  |  n3   0   -n1 |
%                 |               |
%                 | -n2   n1   0  |
%                 `-             -'
%
%   The vector N of matrix entries (n1, n2, n3) is the rotation axis
%   of `m' and theta is M's rotation angle. Note that N and theta
%   are not unique.
%
%   Let
%
%      C = cos(theta/2)
%      s = sin(theta/2)
%
%   Then the unit quaternions `q' corresponding to `m' are
%
%      q = +/- ( C, s*n1, s*n2, s*n3 )
%
%   The mappings between quaternions and the corresponding rotations
%   are carried out by the Mice routines
%
%      cspice_q2m {quaternion to matrix}
%      cspice_m2q {matrix to quaternion}
%
%   cspice_m2q always returns a quaternion with scalar part greater than
%   or equal to zero.
%
%
%   SPICE Quaternion Multiplication Formula
%   ---------------------------------------
%
%   Given a SPICE quaternion
%
%      q = ( q0, q1, q2, q3 )
%
%   corresponding to rotation axis `a' and angle theta as above, we can
%   represent `q' using "scalar + vector" notation as follows:
%
%      s =   q0           = cos(theta/2)
%
%      v = ( q1, q2, q3 ) = sin(theta/2) * a
%
%      q = s + v
%
%   Let `q1' and `q2' be SPICE quaternions with respective scalar
%   and vector parts s1, s2 and v1, v2:
%
%      q1 = s1 + v1
%      q2 = s2 + v2
%
%   We represent the dot product of v1 and v2 by
%
%      <v1, v2>
%
%   and the cross product of v1 and v2 by
%
%      v1 x v2
%
%   Then the SPICE quaternion product is
%
%      q1*q2 = s1*s2 - <v1,v2>  + s1*v2 + s2*v1 + (v1 x v2)
%
%   If `q1' and `q2' represent the rotation matrices `m1' and `m2'
%   respectively, then the quaternion product
%
%      q1*q2
%
%   represents the matrix product
%
%      m1*m2
%
%-Exceptions
%
%   1)  If `q' is not a unit quaternion, the output matrix `r' is
%       the rotation matrix that is the result of converting
%       normalized `q' to a rotation matrix.
%
%   2)  If `q' is the zero quaternion, the output matrix `r' is
%       the identity matrix.
%
%   3)  If the input argument `q' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   4)  If the input argument `q' is not of the expected type, or it
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
%       Edited the header to comply with NAIF standard. Adde complete code
%       example.
%
%       Extended -I/O and -Particulars sections.
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
%   quaternion to matrix
%
%-&

function [r] = cspice_q2m( q )

   switch nargin
      case 1

         q = zzmice_dp(q);

      otherwise

         error ( 'Usage: [_r(3,3)_] = cspice_q2m( _q(4)_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      [r] = mice('q2m_c', q );
   catch spiceerr
      rethrow(spiceerr)
   end


