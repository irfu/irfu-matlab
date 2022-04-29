%-Abstract
%
%   CSPICE_POLYDS computes the value of a polynomial and its first
%   `nderiv' derivatives at the value `t'.
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
%      coeffs   the coefficients of the polynomial that is to be evaluated.
%
%               [deg+1,1] = size(coeffs); double = class(coeffs)
%
%               The first element of this array should be the constant
%               term, the second element the linear coefficient, the third
%               term the quadratic coefficient, and so on. The number of
%               coefficients supplied should be one more than `deg'.
%
%                  f(x) =   coeffs(1) + coeffs(2)*x + coeffs(3)*x^2
%
%                         + coeffs(4)*x^4 + ... + coeffs(deg+1)*x^deg
%
%      deg      the degree of the polynomial to be evaluated.
%
%               [1,1] = size(deg); int32 = class(deg)
%
%               `deg' should be one less than the number of coefficients
%               supplied.
%
%      nderiv   the number of derivatives to compute.
%
%               [1,1] = size(nderiv); int32 = class(nderiv)
%
%               If `nderiv' is zero, only the polynomial will be evaluated.
%               If nderiv = 1, then the polynomial and its first derivative
%               will be evaluated, and so on. If the value of `nderiv' is
%               negative, the routine returns immediately.
%
%      t        the point at which the polynomial and its derivatives should
%               be evaluated.
%
%               [1,1] = size(t); double = class(t)
%
%   the call:
%
%      [p] = cspice_polyds( coeffs, deg, nderiv, t )
%
%   returns:
%
%      p        an array containing the value of the polynomial and its
%               derivatives evaluated at `t'.
%
%               [nderiv+1,1] = size(p); double = class(p)
%
%               The first element of the array contains the value of `p' at
%               `t'. The second element of the array contains the value of
%               the first derivative of `p' at `t' and so on. The `nderiv' +
%               1'st element of the array contains the nderiv'th derivative
%               of `p' evaluated at `t'.
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
%   1) For the polynomial
%
%         f(x) = 1 + 3*x + 0.5*x^2 + x^3 + 0.5*x^4 - x^5 + x^6
%
%      the coefficient set
%
%         Degree  coeffs
%         ------  ------
%         0       1
%         1       3
%         2       0.5
%         3       1
%         4       0.5
%         5      -1
%         6       1
%
%      Compute the value of the polynomial and it's first
%      3 derivatives at the value `t' = 1.0. We expect:
%
%         Derivative Number     t = 1
%         ------------------    -----
%         f(x)         0        6
%         f'(x)        1        10
%         f''(x)       2        23
%         f'''(x)      3        78
%
%
%      Example code begins here.
%
%
%      function polyds_ex1()
%
%         %
%         % Local constants.
%         %
%         NDERIV =   3;
%
%         %
%         % Local variables.
%         %
%
%         coeffs = [1.0,3.0,0.5,1.0,0.5,-1.0,1.0]';
%
%         t      = 1.0;
%         deg    = 6;
%
%         [p]    = cspice_polyds( coeffs, deg, NDERIV, t );
%
%         for i=1:NDERIV+1
%
%            fprintf( 'P = %f\n', p(i) )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      P = 6.000000
%      P = 10.000000
%      P = 23.000000
%      P = 78.000000
%
%
%-Particulars
%
%   This routine uses the user supplied coefficients (coeffs)
%   to evaluate a polynomial (having these coefficients) and its
%   derivatives at the point `t'. The zero'th derivative of the
%   polynomial is regarded as the polynomial itself.
%
%-Exceptions
%
%   1)  If `nderiv' is less than zero, an error is signaled by the Mice
%       interface.
%
%   2)  If the degree of the polynomial is less than 0, an error is
%       signaled by the Mice interface.
%
%   3)  If any of the input arguments, `coeffs', `deg', `nderiv' or
%       `t', is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   4)  If any of the input arguments, `coeffs', `deg', `nderiv' or
%       `t', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%   5)  If the number of elements in `coeffs' is less than deg+1, an error
%       is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  Depending on the coefficients the user should be careful when
%       taking high order derivatives. As the example shows, these
%       can get big in a hurry. In general the coefficients of the
%       derivatives of a polynomial grow at a rate greater
%       than N! (N factorial).
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
%   -Mice Version 1.0.0, 19-JUL-2021 (JDR)
%
%-Index_Entries
%
%   compute a polynomial and its derivatives
%
%-&
function [p] = cspice_polyds( coeffs, deg, nderiv, t )

   switch nargin
      case 4

         coeffs = zzmice_dp(coeffs);
         deg    = zzmice_int(deg);
         nderiv = zzmice_int(nderiv, [0, int32(inf)]);
         t      = zzmice_dp(t);

      otherwise

         error ( [ 'Usage: [p(nderiv+1)] = '                               ...
                   'cspice_polyds( coeffs(deg+1), deg, nderiv, t )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [p] = mice('polyds_c', coeffs, deg, nderiv, t);
   catch spiceerr
      rethrow(spiceerr)
   end
