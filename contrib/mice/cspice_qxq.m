%-Abstract
%
%   CSPICE_QXQ multiplies two quaternions.
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
%      q1       a 4-vector representing a SPICE-style quaternion.
%
%               [4,1] = size(q1); double = class(q1)
%
%               See the discussion of 'Quaternion Styles' in the
%               -Particulars section below.
%
%               Note that multiple styles of quaternions are in use.
%               This routine will not work properly if the input
%               quaternions do not conform to the SPICE convention.
%
%      q2       a second SPICE-style quaternion.
%
%               [4,1] = size(q2); double = class(q2)
%
%   the call:
%
%      [qout] = cspice_qxq( q1, q2 )
%
%   returns:
%
%      qout     4-vector representing the quaternion product
%
%                  q1 * q2
%
%               Representing q(i) as the sums of scalar (real)
%               part s(i) and vector (imaginary) part v(i)
%               respectively,
%
%                  q1 = s1 + v1
%                  q2 = s2 + v2
%
%               qout has scalar part s3 defined by
%
%                  s3 = s1 * s2 - <v1, v2>
%
%               and vector part v3 defined by
%
%                  v3 = s1 * v2  +  s2 * v1  +  v1 x v2
%
%               where the notation < , > denotes the inner
%               product operator and x indicates the cross
%               product operator.
%
%               [4,1] = size(qout); double = class(qout)
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
%   1) Given the 'basis' quaternions:
%
%         qid:  ( 1.0, 0.0, 0.0, 0.0 )
%         qi :  ( 0.0, 1.0, 0.0, 0.0 )
%         qj :  ( 0.0, 0.0, 1.0, 0.0 )
%         qk :  ( 0.0, 0.0, 0.0, 1.0 )
%
%      the following quaternion products give these results:
%
%          Product       Expected result
%         -----------   ----------------------
%          qi  * qj     ( 0.0, 0.0, 0.0, 1.0 )
%          qj  * qk     ( 0.0, 1.0, 0.0, 0.0 )
%          qk  * qi     ( 0.0, 0.0, 1.0, 0.0 )
%          qi  * qi     (-1.0, 0.0, 0.0, 0.0 )
%          qj  * qj     (-1.0, 0.0, 0.0, 0.0 )
%          qk  * qk     (-1.0, 0.0, 0.0, 0.0 )
%          qid * qi     ( 0.0, 1.0, 0.0, 0.0 )
%          qi  * qid    ( 0.0, 1.0, 0.0, 0.0 )
%          qid * qj     ( 0.0, 0.0, 1.0, 0.0 )
%
%      The following code example uses QXQ to produce these results.
%
%
%      Example code begins here.
%
%
%      function qxq_ex1()
%
%         %
%         % Let `qid', `qi', `qj', `qk' be the 'basis'
%         % quaternions.
%         %
%         qid = [1.0,  0.0,  0.0,  0.0]';
%         qi  = [0.0,  1.0,  0.0,  0.0]';
%         qj  = [0.0,  0.0,  1.0,  0.0]';
%         qk  = [0.0,  0.0,  0.0,  1.0]';
%
%         %
%         % Compute:
%         %
%         %    qi x qj = qk
%         %    qj x qk = qi
%         %    qk x qi = qj
%         %
%         [qout] = cspice_qxq( qi, qj );
%         fprintf( 'qi x qj  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '     qk  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                               qk(1), qk(2), qk(3), qk(4) )
%         fprintf( ' \n' )
%
%         [qout] = cspice_qxq( qj, qk );
%         fprintf( 'qj x qk  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '     qi  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                               qi(1), qi(2), qi(3), qi(4) )
%         fprintf( ' \n' )
%
%         [qout] = cspice_qxq( qk, qi );
%         fprintf( 'qk x qi  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '     qj  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                               qj(1), qj(2), qj(3), qj(4) )
%         fprintf( ' \n' )
%
%         %
%         % Compute:
%         %
%         %    qi x qi  ==  -qid
%         %    qj x qj  ==  -qid
%         %    qk x qk  ==  -qid
%         %
%         [qout] = cspice_qxq( qi, qi );
%         fprintf( 'qi x qi  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '     qid = %7.1f %7.1f %7.1f %7.1f\n', ...
%                           qid(1), qid(2), qid(3), qid(4) )
%         fprintf( ' \n' )
%
%         [qout] = cspice_qxq( qj, qj );
%         fprintf( 'qj x qj  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '     qid = %7.1f %7.1f %7.1f %7.1f\n', ...
%                           qid(1), qid(2), qid(3), qid(4) )
%         fprintf( ' \n' )
%
%         [qout] = cspice_qxq( qk, qk );
%         fprintf( 'qk x qk  = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '     qid = %7.1f %7.1f %7.1f %7.1f\n', ...
%                           qid(1), qid(2), qid(3), qid(4) )
%         fprintf( ' \n' )
%
%         %
%         % Compute:
%         %
%         %    qid x qi  = qi
%         %    qi  x qid = qi
%         %    qid x qj  = qj
%         %
%         [qout] = cspice_qxq( qid, qi );
%         fprintf( 'qid x qi = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '      qi = %7.1f %7.1f %7.1f %7.1f\n', ...
%                               qi(1), qi(2), qi(3), qi(4) )
%         fprintf( ' \n' )
%
%         [qout] = cspice_qxq( qi, qid );
%         fprintf( 'qi x qid = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '      qi = %7.1f %7.1f %7.1f %7.1f\n', ...
%                               qi(1), qi(2), qi(3), qi(4) )
%         fprintf( ' \n' )
%
%         [qout] = cspice_qxq( qid, qj );
%         fprintf( 'qid x qj = %7.1f %7.1f %7.1f %7.1f\n', ...
%                       qout(1), qout(2), qout(3), qout(4) )
%         fprintf( '      qj = %7.1f %7.1f %7.1f %7.1f\n', ...
%                               qj(1), qj(2), qj(3), qj(4) )
%         fprintf( ' \n' )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      qi x qj  =     0.0     0.0     0.0     1.0
%           qk  =     0.0     0.0     0.0     1.0
%
%      qj x qk  =     0.0     1.0     0.0     0.0
%           qi  =     0.0     1.0     0.0     0.0
%
%      qk x qi  =     0.0     0.0     1.0     0.0
%           qj  =     0.0     0.0     1.0     0.0
%
%      qi x qi  =    -1.0     0.0     0.0     0.0
%           qid =     1.0     0.0     0.0     0.0
%
%      qj x qj  =    -1.0     0.0     0.0     0.0
%           qid =     1.0     0.0     0.0     0.0
%
%      qk x qk  =    -1.0     0.0     0.0     0.0
%           qid =     1.0     0.0     0.0     0.0
%
%      qid x qi =     0.0     1.0     0.0     0.0
%            qi =     0.0     1.0     0.0     0.0
%
%      qi x qid =     0.0     1.0     0.0     0.0
%            qi =     0.0     1.0     0.0     0.0
%
%      qid x qj =     0.0     0.0     1.0     0.0
%            qj =     0.0     0.0     1.0     0.0
%
%
%   2) Compute the composition of two rotation matrices by
%      converting them to quaternions and computing their
%      product, and by directly multiplying the matrices.
%
%      Example code begins here.
%
%
%      function qxq_ex2()
%
%         %
%         % Local variables
%         %
%         cmat1 = [ [1.0,  0.0,  0.0]', ...
%                   [0.0, -1.0,  0.0]', ...
%                   [0.0,  0.0, -1.0]'  ]';
%
%         cmat2 = [ [0.0,  1.0,  0.0]', ...
%                   [1.0,  0.0,  0.0]', ...
%                   [0.0,  0.0, -1.0]'  ]';
%
%         %
%         % Convert the C-matrices to quaternions.
%         %
%         [q1] = cspice_m2q( cmat1 );
%         [q2] = cspice_m2q( cmat2 );
%
%         %
%         % Find the product.
%         %
%         [qout] = cspice_qxq( q1, q2 );
%
%         %
%         % Convert the result to a C-matrix.
%         %
%         [cmout] = cspice_q2m( qout );
%
%         fprintf( 'Using quaternion product:\n' )
%         fprintf( '%9.4f %9.4f %9.4f\n',             ...
%                  cmout(1,1), cmout(1,2), cmout(1,3) )
%         fprintf( '%9.4f %9.4f %9.4f\n',             ...
%                  cmout(2,1), cmout(2,2), cmout(2,3) )
%         fprintf( '%9.4f %9.4f %9.4f\n',             ...
%                  cmout(3,1), cmout(3,2), cmout(3,3) )
%
%         %
%         % Multiply `cmat1' and `cmat2' directly.
%         %
%         cmout = cmat1 * cmat2;
%
%         fprintf( 'Using matrix product:\n' )
%         fprintf( '%9.4f %9.4f %9.4f\n',             ...
%                  cmout(1,1), cmout(1,2), cmout(1,3) )
%         fprintf( '%9.4f %9.4f %9.4f\n',             ...
%                  cmout(2,1), cmout(2,2), cmout(2,3) )
%         fprintf( '%9.4f %9.4f %9.4f\n',             ...
%                  cmout(3,1), cmout(3,2), cmout(3,3) )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Using quaternion product:
%         0.0000    1.0000    0.0000
%        -1.0000    0.0000    0.0000
%         0.0000    0.0000    1.0000
%      Using matrix product:
%         0.0000    1.0000    0.0000
%        -1.0000    0.0000    0.0000
%         0.0000    0.0000    1.0000
%
%
%-Particulars
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
%   Mice function interfaces ALWAYS use SPICE quaternions.
%   Quaternions of any other style must be converted to SPICE
%   quaternions before they are passed to Mice functions.
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
%               sin(theta/2) * a(1),
%               sin(theta/2) * a(2),
%               sin(theta/2) * a(3)  )
%
%   while the engineering quaternions representing `m' are
%
%      (+/-) ( -sin(theta/2) * a(1),
%              -sin(theta/2) * a(2),
%              -sin(theta/2) * a(3),
%               cos(theta/2)         )
%
%   For both styles of quaternions, if a quaternion `q' represents
%   a rotation matrix `m', then -q represents `m' as well.
%
%   Given an engineering quaternion
%
%      qeng   = ( q1,  q2,  q3,  q4 )
%
%   the equivalent SPICE quaternion is
%
%      qspice = ( q4, -q1, -q2, -q3 )
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
%   Then if the unit-length SPICE quaternion `q' represents `m', where
%
%      q = (q1, q2, q3, q4)
%
%   the elements of `m' are derived from the elements of `q' as follows:
%
%        .-                                                           -.
%        |            2    2                                           |
%        |  1 - 2*( q3 + q4 )   2*(q2*q3 - q1*q4)   2*(q2*q4 + q1*q3)  |
%        |                                                             |
%        |                                                             |
%        |                                2    2                       |
%    m = |  2*(q2*q3 + q1*q4)   1 - 2*( q2 + q4 )   2*(q3*q4 - q1*q2)  |
%        |                                                             |
%        |                                                             |
%        |                                                    2    2   |
%        |  2*(q2*q4 - q1*q3)   2*(q3*q4 + q1*q2)   1 - 2*( q2 + q3 )  |
%        `-                                                           -'
%
%   Note that substituting the elements of -q for those of `q' in the
%   right hand side leaves each element of `m' unchanged; this shows
%   that if a quaternion `q' represents a matrix `m', then so does the
%   quaternion -q.
%
%   To map the rotation matrix `m' to a unit quaternion, we start by
%   decomposing the rotation matrix as a sum of symmetric
%   and skew-symmetric parts:
%
%                                        2
%      m = [ I  +  (1-cos(theta)) * omega  ] + [ sin(theta) * omega ]
%
%                       symmetric                 skew-symmetric
%
%
%   `omega' is a skew-symmetric matrix of the form
%
%                 .-               -.
%                 |   0   -n3   n2  |
%                 |                 |
%       omega  =  |   n3   0   -n1  |
%                 |                 |
%                 |  -n2   n1   0   |
%                 `-               -'
%
%   The vector `n' of matrix entries (n1, n2, n3) is the rotation axis
%   of `m' and `theta' is m's rotation angle. Note that `n' and `theta'
%   are not unique.
%
%   Let
%
%      cth = cos(theta/2)
%      sth = sin(theta/2)
%
%   Then the unit quaternions `q' corresponding to `m' are
%
%      q = +/- ( cth, sth*n1, sth*n2, sth*n3 )
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
%      q = ( q1, q2, q3, q4 )
%
%   corresponding to rotation axis `a' and angle `theta' as above, we can
%   represent `q' using 'scalar + vector' notation as follows:
%
%      s =   q1           = cos(theta/2)
%
%      v = ( q2, q3, q4 ) = sin(theta/2) * a
%
%      q = s + v
%
%   Let `quat1' and `quat2' be SPICE quaternions with respective scalar
%   and vector parts `s1', `s2' and `v1', `v2':
%
%      quat1 = s1 + v1
%      quat2 = s2 + v2
%
%   We represent the dot product of `v1' and `v2' by
%
%      <v1, v2>
%
%   and the cross product of `v1' and `v2' by
%
%      v1 x v2
%
%   Then the SPICE quaternion product is
%
%      quat1*quat2 = s1*s2 - <v1,v2>  + s1*v2 + s2*v1 + (v1 x v2)
%
%   If `quat1' and `quat2' represent the rotation matrices `m1' and `m2'
%   respectively, then the quaternion product
%
%      quat1*quat1
%
%   represents the matrix product
%
%      m1*m2
%
%-Exceptions
%
%   1)  If any of the input arguments, `q1' or `q2', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   2)  If any of the input arguments, `q1' or `q2', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
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
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   quaternion times quaternion
%   multiply quaternion by quaternion
%
%-&
function [qout] = cspice_qxq( q1, q2 )

   switch nargin
      case 2

         q1 = zzmice_dp(q1);
         q2 = zzmice_dp(q2);

      otherwise

         error ( 'Usage: [qout(4)] = cspice_qxq( q1(4), q2(4) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [qout] = mice('qxq_c', q1, q2);
   catch spiceerr
      rethrow(spiceerr)
   end
