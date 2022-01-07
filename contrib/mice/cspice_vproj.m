%-Abstract
%
%   CSPICE_VPROJ computes the projection of one 3-dimensional vector onto
%   another 3-dimensional vector.
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
%      a        double precision, 3-dimensional vector(s).
%
%               [3,n] = size(a); double = class(a)
%
%               This vector is to be projected onto the vector `b'.
%
%      b        double precision, 3-dimensional vector(s).
%
%               [3,n] = size(b); double = class(b)
%
%               This vector is the vector which receives the projection.
%
%               An implicit assumption exists that `a' and `b' are specified
%               in the same reference frame. If this is not the case, the
%               numerical result has no meaning.
%
%   the call:
%
%      [p] = cspice_vproj( a, b )
%
%   returns:
%
%      p        the double precision, 3-dimensional vector(s) containing the
%               projection of `a' onto `b'.
%
%               [3,n] = size(p); double = class(p)
%
%               (`p' is necessarily parallel to `b'.) If `b' is the zero
%               vector then `p' will be returned as the zero vector.
%
%               `p' returns with the same vectorization measure, N, as
%               `a' and `b'.
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
%   1) Define two sets of vectors and compute the projection of
%      each vector of the first set on the corresponding vector of
%      the second set.
%
%      Example code begins here.
%
%
%      function vproj_ex1()
%
%         %
%         % Define two vector sets.
%         %
%         a = [ [ 6, 6, 6]', ...
%               [ 6, 6, 6]', ...
%               [ 6, 6, 0]', ...
%               [ 6, 0, 0]' ];
%
%         b = [ [ 2, 0, 0]', ...
%               [-3, 0, 0]', ...
%               [ 0, 7, 0]', ...
%               [ 0, 0, 9]' ];
%
%         %
%         % Calculate the projection.
%         %
%         p = cspice_vproj( a, b );
%
%         for i=1:4
%            fprintf( 'Vector A  : %5.1f %5.1f %5.1f\n', ...
%                             a(1,i), a(2,i), a(3,i) )
%            fprintf( 'Vector B  : %5.1f %5.1f %5.1f\n', ...
%                             b(1,i), b(2,i), b(3,i) )
%            fprintf( 'Projection: %5.1f %5.1f %5.1f\n\n', ...
%                             p(1,i), p(2,i), p(3,i) )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Vector A  :   6.0   6.0   6.0
%      Vector B  :   2.0   0.0   0.0
%      Projection:   6.0   0.0   0.0
%
%      Vector A  :   6.0   6.0   6.0
%      Vector B  :  -3.0   0.0   0.0
%      Projection:   6.0  -0.0  -0.0
%
%      Vector A  :   6.0   6.0   0.0
%      Vector B  :   0.0   7.0   0.0
%      Projection:   0.0   6.0   0.0
%
%      Vector A  :   6.0   0.0   0.0
%      Vector B  :   0.0   0.0   9.0
%      Projection:   0.0   0.0   0.0
%
%
%-Particulars
%
%   Given any vectors `a' and `b', there is a unique decomposition of `a' as
%   a sum v + p such that `v', the dot product of `v' and `b', is zero, and
%   the dot product of `p' with `b' is equal the product of the lengths of
%   `p' and `b'. `p' is called the projection of `a' onto `b'. It can be
%   expressed mathematically as
%
%      dot(a,b)
%      -------- * b
%      dot(b,b)
%
%   (This is not necessarily the prescription used to compute the
%   projection. It is intended only for descriptive purposes.)
%
%-Exceptions
%
%   1)  If any of the input arguments, `a' or `b', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   2)  If any of the input arguments, `a' or `b', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   3)  If the input vectorizable arguments `a' and `b' do not have
%       the same measure of vectorization (N), an error is signaled by
%       the Mice interface.
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
%   [1]  G. Thomas and R. Finney, "Calculus and Analytic Geometry,"
%        7th Edition, Addison Wesley, 1988.
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
%       Changed output argument name "vproj" to "p".
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and reformatted example's output.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 12-MAR-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   3-vector projection
%
%-&

function [p] = cspice_vproj( a, b)

   switch nargin
      case 2

         a = zzmice_dp(a);
         b = zzmice_dp(b);

      otherwise

         error ( 'Usage: [_p(3)_] = cspice_vproj(_a(3)_, _b(3)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [p] = mice('vproj_c', a, b);
   catch spiceerr
      rethrow(spiceerr)
   end

