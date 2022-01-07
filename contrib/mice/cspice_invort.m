%-Abstract
%
%   CSPICE_INVORT constructs the inverse of a 3x3 matrix with orthogonal
%   columns and non-zero column norms using a numerically stable algorithm.
%   The rows of the output matrix are the columns of the input matrix divided
%   by the length squared of the corresponding columns.
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
%      m        a 3x3 matrix.
%
%               [3,3] = size(m); double = class(m)
%
%   the call:
%
%      [mit] = cspice_invort( m )
%
%   returns:
%
%      mit      the matrix obtained by transposing `m' and dividing the rows
%               by squares of their norms.
%
%               [3,3] = size(mit); double = class(mit)
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
%   1) Given a double precision 3x3 matrix with mutually orthogonal
%      rows of arbitrary length, compute its inverse. Check that the
%      original matrix times the computed inverse produces the
%      identity matrix.
%
%      Example code begins here.
%
%
%      function invort_ex1()
%
%         %
%         % Define a matrix to invert.
%         %
%         m = [ [0.0, -1.0, 0.0]', [0.5,  0.0, 0.0]', [0.0,  0.0, 1.0]' ]';
%
%         fprintf( 'Original Matrix:\n' )
%         for i=1:3
%
%            fprintf( '%16.7f %15.7f %15.7f\n', m(i,1), m(i,2), m(i,3) )
%
%         end
%
%         %
%         % Invert the matrix, then output.
%         %
%         [mout] = cspice_invort( m );
%
%         fprintf( ' \n' )
%         fprintf( 'Inverse Matrix:\n' )
%         for i=1:3
%
%            fprintf( '%16.7f %15.7f %15.7f\n',                            ...
%                     mout(i,1), mout(i,2), mout(i,3) )
%
%         end
%
%         %
%         % Check the `m' times `mout' produces the identity matrix.
%         %
%         imat = m * mout;
%
%         fprintf( ' \n' )
%         fprintf( 'Original times inverse:\n' )
%         for i=1:3
%
%            fprintf( '%16.7f %15.7f %15.7f\n',                            ...
%                     imat(i,1), imat(i,2), imat(i,3) )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Original Matrix:
%             0.0000000      -1.0000000       0.0000000
%             0.5000000       0.0000000       0.0000000
%             0.0000000       0.0000000       1.0000000
%
%      Inverse Matrix:
%             0.0000000       2.0000000       0.0000000
%            -1.0000000       0.0000000       0.0000000
%             0.0000000       0.0000000       1.0000000
%
%      Original times inverse:
%             1.0000000       0.0000000       0.0000000
%             0.0000000       1.0000000       0.0000000
%             0.0000000       0.0000000       1.0000000
%
%
%-Particulars
%
%   Suppose that m is the matrix
%
%          .-                      -.
%          |   A*u    B*v     C*w   |
%          |      1      1       1  |
%          |                        |
%          |   A*u    B*v     C*w   |
%          |      2      2       2  |
%          |                        |
%          |   A*u    B*v     C*w   |
%          |      3      3       3  |
%          `-                      -'
%
%   where the vectors (u , u , u ),  (v , v , v ),  and (w , w , w )
%                       1   2   3      1   2   3          1   2   3
%
%   are unit vectors. This routine produces the matrix:
%
%
%          .-                      -.
%          |   a*u    a*u     a*u   |
%          |      1      2       3  |
%          |                        |
%          |   b*v    b*v     b*v   |
%          |      1      2       3  |
%          |                        |
%          |   c*w    c*w     c*w   |
%          |      1      2       3  |
%          `-                      -'
%
%   where a = 1/A, b = 1/B, and c = 1/C.
%
%-Exceptions
%
%   1)  If any of the columns of `m' have zero length, the error
%       SPICE(ZEROLENGTHCOLUMN) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If any column is too short to allow computation of the
%       reciprocal of its length without causing a floating point
%       overflow, the error SPICE(COLUMNTOOSMALL) is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If the input argument `m' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   4)  If the input argument `m' is not of the expected type, or it
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
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Updated the header to comply with NAIF standard. Added
%       complete code example to -Examples section. Extended -Abstract
%       section.
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
%   -Mice Version 1.0.0, 14-NOV-2013 (EDW) (SCK)
%
%-Index_Entries
%
%   Transpose a matrix and invert the lengths of the rows
%   Invert a pseudo orthogonal matrix
%
%-&

function [mit] = cspice_invort( m )

   switch nargin
      case 1

         m = zzmice_dp(m);

      otherwise

         error( 'Usage: [mit(3,3)] = cspice_rotmat( m(3,3) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [mit] = mice('invort_c', m );
   catch spiceerr
      rethrow(spiceerr)
   end


