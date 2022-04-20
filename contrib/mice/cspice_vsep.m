%-Abstract
%
%   CSPICE_VSEP finds the separation angle in radians between two double
%   precision, 3-dimensional vectors. This angle is defined as zero
%   if either vector is zero.
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
%      v1,
%      v2       two double precision 3-dimensional vectors.
%
%               [3,n] = size(v1); double = class(v1)
%               [3,n] = size(v2); double = class(v2)
%
%               Either `v1' or `v2', or both, may be the zero vector.
%
%               An implicit assumption exists that `v1' and `v2' are
%               specified in the same reference frame. If this is not
%               the case, the numerical result of this routine has no
%               meaning.
%
%   the call:
%
%      [vsep] = cspice_vsep( v1, v2 )
%
%   returns:
%
%      vsep     the value(s) of the angular separation between `v1' and `v2'
%               expressed in radians.
%
%               [1,n] = size(vsep); double = class(vsep)
%
%               cspice_vsep is strictly non-negative. If either `v1' or `v2'
%               is the zero vector, then cspice_vsep is defined to be 0
%               radians.
%
%               `vsep' returns with the same vectorization measure, N, as
%               `v1' and `v2'
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
%   1) Define two sets of 3-dimensional vectors and compute the
%      angular separation between each vector in first set and the
%      corresponding vector in the second set.
%
%
%      Example code begins here.
%
%
%      function vsep_ex1()
%
%         %
%         % Define a set of vectors, calculate angular
%         % separation as measured in radians.
%         %
%         v1 = [1; 0; 0];
%         v2 = [0; 1; 0];
%
%         sep = cspice_vsep( v1, v2 );
%         disp( 'Scalar:' )
%         fprintf( '   Vector 1:  %3.1f  %3.1f  %3.1f\n', v1  )
%         fprintf( '   Vector 2:  %3.1f  %3.1f  %3.1f\n', v2  )
%         fprintf( '   Angular separation: %10.6f\n\n', sep )
%
%         %
%         % Instead of two calls with 3-vectors,
%         % vectorize the input as two 3X2 array.
%         %
%         v1 = [ [1; 0; 0], [1; 0; 0] ];
%         v2 = [ [1; 0; 0], [0; 1; 0] ];
%
%         sep = cspice_vsep( v1, v2 );
%         disp( 'Vectorized:' )
%         for i=1:2
%            fprintf( '   Vector 1: %3.1f %3.1f %3.1f\n',  v1(:,i))
%            fprintf( '   Vector 2: %3.1f %3.1f %3.1f\n',  v2(:,i))
%            fprintf( '   Angular separation: %10.6f\n\n', sep(i) )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Scalar:
%         Vector 1:  1.0  0.0  0.0
%         Vector 2:  0.0  1.0  0.0
%         Angular separation:   1.570796
%
%      Vectorized:
%         Vector 1: 1.0 0.0 0.0
%         Vector 2: 1.0 0.0 0.0
%         Angular separation:   0.000000
%
%         Vector 1: 1.0 0.0 0.0
%         Vector 2: 0.0 1.0 0.0
%         Angular separation:   1.570796
%
%
%-Particulars
%
%   In the plane, it is a simple matter to calculate the angle
%   between two vectors once the two vectors have been made to be
%   unit length. Then, since the two vectors form the two equal
%   sides of an isosceles triangle, the length of the third side
%   is given by the expression
%
%      length = 2.0 * sin ( cspice_vsep/2.0 )
%
%   The length is given by the magnitude of the difference of the
%   two unit vectors
%
%      length = norm ( u1 - u2 )
%
%   Once the length is found, the value of cspice_vsep may be calculated
%   by inverting the first expression given above as
%
%      cspice_vsep = 2.0 * arcsin ( length/2.0 )
%
%   This expression becomes increasingly unstable when cspice_vsep gets
%   larger than pi/2 radians or 90 degrees. In this situation (which
%   is easily detected by determining the sign of the dot product of
%   `v1' and `v2') the supplementary angle is calculated first and
%   then cspice_vsep is given by
%
%         cspice_vsep = pi - supplementary_angle
%
%-Exceptions
%
%   1)  If any of the input arguments, `v1' or `v2', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   2)  If any of the input arguments, `v1' or `v2', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   3)  If the input vectorizable arguments `v1' and `v2' do not have
%       the same measure of vectorization (N), an error is signaled by
%       the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The user is required to insure that the input vectors will not
%       cause floating point overflow upon calculation of the vector
%       dot product since no error detection or correction code is
%       implemented. In practice, this is not a significant
%       restriction.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and modified code example to produce
%       formatted output.
%
%       Changed output argument name "sep" to "vsep" to comply with NAIF
%       standard.
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
%   -Mice Version 1.0.2, 17-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Corrections made to -Version section numbering. 10-APR-2010
%       notation now numbered as 1.0.1.
%
%   -Mice Version 1.0.1, 10-APR-2010 (EDW)
%
%       Edits to header -I/O section.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   angular separation of 3-dimensional vectors
%
%-&

function [vsep] = cspice_vsep(v1, v2)

   switch nargin
      case 2

         v1 = zzmice_dp(v1);
         v2 = zzmice_dp(v2);

      otherwise

         error ( 'Usage: [_vsep_] = cspice_vsep(_v1(3)_, _v2(3)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [vsep] = mice('vsep_c',v1, v2);
   catch spiceerr
      rethrow(spiceerr)
   end
