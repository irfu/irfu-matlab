%-Abstract
%
%   CSPICE_EUL2M constructs a 3x3, double precision rotation matrix
%   from a set of Euler angles and the corresponding rotation axes.
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
%      angle3,
%      angle2,
%      angle1,
%      axis3,
%      axis2,
%      axis1    respectively, set(s) of three angles and three coordinate
%               axis numbers; each pair angleX and axisX specifies a
%               coordinate transformation consisting of a rotation by angleX
%               radians about the coordinate axis indexed by axisX.
%
%               [1,n] = size(angle3); double = class(angle3)
%               [1,n] = size(angle2); double = class(angle2)
%               [1,n] = size(angle1); double = class(angle1)
%               [1,1] = size(axis3); int32 = class(axis3)
%               [1,1] = size(axis2); int32 = class(axis2)
%               [1,1] = size(axis1); int32 = class(axis1)
%
%               These coordinate transformations are typically symbolized
%               by
%
%                  [ angleX ]     .
%                            axisX
%
%               See the -Particulars section below for details concerning
%               this notation.
%
%               Note that these coordinate transformations rotate vectors
%               by
%
%                  -angleX
%
%               radians about the axis indexed by axisX.
%
%               The values of axisX may be 1, 2, or 3, indicating the X,
%               Y, and Z axes respectively.
%
%   the call:
%
%      [r] = cspice_eul2m( angle3, angle2, angle1, axis3, axis2, axis1 )
%
%   returns:
%
%      r        a rotation matrix(es) representing the composition of the
%               rotations defined by the input angle-axis pairs.
%
%               If [1,1] = size(angle3) then [3,3]   = size(r)
%               If [1,n] = size(angle3) then [3,3,n] = size(r)
%                                             double = class(r)
%
%               Together, the three pairs specify a composite
%               transformation that is the result of performing the rotations
%               about the axes indexed by `axis1', `axis2', and `axis3', in
%               that order. So,
%
%                  r = [ angle3 ]      [ angle2 ]      [ angle1 ]
%                                axis3           axis2           axis1
%
%               See the -Particulars section below for details concerning
%               this notation.
%
%               The resulting matrix `r' may be thought of as a coordinate
%               transformation; applying it to a vector yields the
%               vector's coordinates in the rotated system.
%
%               Viewing `r' as a coordinate transformation matrix, the
%               basis that `r' transforms vectors to is created by rotating
%               the original coordinate axes first by `angle1' radians
%               about the coordinate axis indexed by `axis1', next by
%               `angle2' radians about the coordinate axis indexed by
%               `axis2', and finally by `angle3' radians about coordinate
%               axis indexed by `axis3'. At the second and third steps of
%               this process, the coordinate axes about which rotations
%               are performed belong to the bases resulting from the
%               previous rotations.
%
%               `r' returns with the same vectorization measure, N,
%               as `angle3', `angle2' and `angle1'.
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
%   1) Create a rotation matrix for a single coordinate rotation of
%      90 degrees about the Z axis, and compute the vector resulting
%      from applying that rotation to the +X unit vector.
%
%      Example code begins here.
%
%
%      function eul2m_ex1()
%         %
%         % Create the rotation matrix for a single coordinate
%         % rotation of 90 degrees about the Z axis. As the
%         % second and third angles are 0, the final two axes IDs,
%         % 1, 1, have no effect for in this example.
%         %
%         rot = cspice_eul2m( cspice_halfpi, 0, 0, 3, 1, 1 );
%
%         %
%         % Output the result of rotating the +x unit vector
%         % using the 'rot' matrix.
%         %
%         fprintf( '%9.3f %9.3f %9.3f\n', rot * [1; 0; 0 ]);
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%          0.000    -1.000     0.000
%
%
%-Particulars
%
%   A word about notation: the symbol
%
%      [ x ]
%           i
%
%   indicates a rotation of x radians about the ith coordinate axis.
%   To be specific, the symbol
%
%      [ x ]
%           1
%
%   indicates a coordinate system rotation of x radians about the
%   first, or x-, axis; the corresponding matrix is
%
%      .-                    -.
%      |  1      0       0    |
%      |                      |
%      |  0    cos(x)  sin(x) |
%      |                      |
%      |  0   -sin(x)  cos(x) |
%      `-                    -'
%
%   Remember, this is a COORDINATE SYSTEM rotation by x radians; this
%   matrix, when applied to a vector, rotates the vector by -x
%   radians, not x radians. Applying the matrix to a vector yields
%   the vector's representation relative to the rotated coordinate
%   system.
%
%   The analogous rotation about the second, or y-, axis is
%   represented by
%
%      [ x ]
%           2
%
%   which symbolizes the matrix
%
%      .-                    -.
%      | cos(x)   0   -sin(x) |
%      |                      |
%      |  0       1      0    |
%      |                      |
%      | sin(x)   0    cos(x) |
%      `-                    -'
%
%   and the analogous rotation about the third, or z-, axis is
%   represented by
%
%      [ x ]
%           3
%
%   which symbolizes the matrix
%
%      .-                    -.
%      |  cos(x)  sin(x)   0  |
%      |                      |
%      | -sin(x)  cos(x)   0  |
%      |                      |
%      |  0        0       1  |
%      `-                    -'
%
%   From time to time, (depending on one's line of work, perhaps) one
%   may happen upon a pair of coordinate systems related by a
%   sequence of rotations. For example, the coordinate system
%   defined by an instrument such as a camera is sometime specified
%   by RA, DEC, and twist; if alpha, delta, and phi are the rotation
%   angles, then the series of rotations
%
%      [ phi ]     [ pi/2  -  delta ]     [ alpha ]
%             3                      2             3
%
%   produces a transformation from inertial to camera coordinates.
%
%   This routine is related to the Mice routine cspice_m2eul, which
%   produces a sequence of Euler angles, given a rotation matrix.
%   This routine is a "left inverse" of cspice_m2eul: the sequence of
%   calls
%
%      [angle3, angle2, angle1] = cspice_m2eul( r, axis3, axis2, axis1 )
%      [r] = cspice_eul2m( angle3, angle2, angle1, axis3, axis2, axis1 )
%
%   preserves `r' to round-off error.
%
%   On the other hand, the sequence of calls
%
%      [r] = cspice_eul2m( angle3, angle2, angle1, axis3, axis2, axis1 )
%      [angle3, angle2, angle1] = cspice_m2eul( r, axis3, axis2, axis1 )
%
%   preserve `angle3', `angle2', and `angle1' only if these angles start
%   out in the ranges that cspice_m2eul's outputs are restricted to.
%
%-Exceptions
%
%   1)  If any of `axis3', `axis2', or `axis1' do not have values in
%
%          { 1, 2, 3 }
%
%       the error SPICE(BADAXISNUMBERS) is signaled by a routine in
%       the call tree of this routine.
%
%   2)  If any of the input arguments, `angle3', `angle2', `angle1',
%       `axis3', `axis2' or `axis1', is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `angle3', `angle2', `angle1',
%       `axis3', `axis2' or `axis1', is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%   4)  If the input vectorizable arguments `angle3', `angle2' and
%       `angle1' do not have the same measure of vectorization (N), an
%       error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  Beware: more than one definition of "RA, DEC and twist"
%       exists.
%
%-Required_Reading
%
%   MICE.REQ
%   ROTATION.REQ
%
%-Literature_References
%
%   [1]  W. Owen, "Galileo Attitude and Camera Models," JPL
%        Interoffice Memorandum 314-323, Nov. 11, 1983. NAIF document
%        number 204.0.
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
%       Edited the header to comply with NAIF standard. Extended arguments
%       description in -I/O section.
%
%       Added example's problem statement and reformatted example's
%       output.
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
%   -Mice Version 1.0.1, 06-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   euler angles to matrix
%
%-&

function [r] = cspice_eul2m(angle3, angle2, angle1, axis3, axis2, axis1)

   switch nargin
      case 6

         angle3 = zzmice_dp(angle3);
         angle2 = zzmice_dp(angle2);
         angle1 = zzmice_dp(angle1);
         axis3  = zzmice_int(axis3);
         axis2  = zzmice_int(axis2);
         axis1  = zzmice_int(axis1);

      otherwise

         error( ['Usage: [_r(3,3)_] = '                         ...
                 'cspice_eul2m(_angle3_, _angle2_, _angle1_, '  ...
                 'axis3, axis2, axis1)']  )

   end

   %
   % Call the MEX library.
   %
   try
      [r] = mice('eul2m_c',angle3,angle2,angle1,axis3,axis2,axis1);
   catch spiceerr
      rethrow(spiceerr)
   end
