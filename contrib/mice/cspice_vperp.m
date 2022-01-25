%-Abstract
%
%   CSPICE_VPERP calculates the component of a vector perpendicular to a
%   second vector.
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
%               It the vector whose component orthogonal to `b' is sought.
%               (There is a unique decomposition of `a' into a sum v + p,
%               where `v' is parallel to `b' and `p' is orthogonal to `b'. We
%               want the component `p'.)
%
%      b        double precision, 3-dimensional vector(s).
%
%               [3,n] = size(b); double = class(b)
%
%               This vector is the vector used as a reference for the
%               decomposition of `a'.
%
%   the call:
%
%      [p] = cspice_vperp( a, b )
%
%   returns:
%
%      p        the double precision, 3-dimensional vector(s) containing the
%               component of `a' that is orthogonal to `b'.
%
%               [3,n] = size(p); double = class(p)
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
%   1) Define two vector sets and compute the component of the vector
%      in the first set perpendicular to the vector in the second set.
%
%      Example code begins here.
%
%
%      function vperp_ex1()
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
%         % Calculate the decomposition.
%         %
%         p = cspice_vperp( a, b );
%
%         for i=1:4
%            fprintf( 'Vector A     : %5.1f %5.1f %5.1f\n', ...
%                                a(1,i), a(2,i), a(3,i) )
%            fprintf( 'Vector B     : %5.1f %5.1f %5.1f\n', ...
%                                b(1,i), b(2,i), b(3,i) )
%            fprintf( 'Perpendicular: %5.1f %5.1f %5.1f\n\n', ...
%                                p(1,i), p(2,i), p(3,i) )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Vector A     :   6.0   6.0   6.0
%      Vector B     :   2.0   0.0   0.0
%      Perpendicular:   0.0   6.0   6.0
%
%      Vector A     :   6.0   6.0   6.0
%      Vector B     :  -3.0   0.0   0.0
%      Perpendicular:   0.0   6.0   6.0
%
%      Vector A     :   6.0   6.0   0.0
%      Vector B     :   0.0   7.0   0.0
%      Perpendicular:   6.0   0.0   0.0
%
%      Vector A     :   6.0   0.0   0.0
%      Vector B     :   0.0   0.0   9.0
%      Perpendicular:   6.0   0.0   0.0
%
%
%-Particulars
%
%   Given and non-zero vector `b' and a vector `a', there is a unique
%   decomposition of `a' as a sum v + p such that `p' is orthogonal
%   to `b' and `v' is parallel to `b'. This routine finds the vector `p'.
%
%   If `b' is a zero vector, `p' will be identical to `a'.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "vperp" to "p".
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and reformatted code example's output.
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
%   -Mice Version 1.0.1, 09-NOV-2012 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-APR-2010 (EDW)
%
%-Index_Entries
%
%   perpendicular component of a 3-vector
%
%-&

function [p] = cspice_vperp( a, b)

   switch nargin
      case 2

         a = zzmice_dp(a);
         b = zzmice_dp(b);

      otherwise

         error ( 'Usage: [_p(3)_] = cspice_vperp(_a(3)_, _b(3)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [p] = mice('vperp_c', a, b);
   catch spiceerr
      rethrow(spiceerr)
   end



