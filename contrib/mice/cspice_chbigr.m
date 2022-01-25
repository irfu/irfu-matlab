%-Abstract
%
%   CSPICE_CHBIGR evaluates an indefinite integral of a Chebyshev expansion
%   at a specified point `x' and returns the value of the input expansion at
%   `x' as well. The constant of integration is selected to make the integral
%   zero when `x' equals the abscissa value x2s(1).
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
%      degp     the degree of the input Chebyshev expansion.
%
%               [1,1] = size(degp); int32 = class(degp)
%
%      cp       an array containing the coefficients of the input Chebyshev
%               expansion.
%
%               [degp+1,1] = size(cp); double = class(cp)
%
%               The coefficient of the i'th Chebyshev polynomial is
%               contained in element cp(i+1), for i = 0 : degp.
%
%      x2s      an array containing the "transformation parameters" of the
%               domain of the expansion.
%
%               [2,1] = size(x2s); double = class(x2s)
%
%               Element x2s(1) contains the midpoint of the interval on
%               which the input expansion is defined; x2s(2) is one-half of
%               the length of this interval; this value is called the
%               interval's "radius."
%
%               The input expansion defines a function f(x) on the
%               interval
%
%                  [ x2s(1)-x2s(2),  x2s(1)+x2s(2) ]
%
%               as follows:
%
%                  Define s = ( x - x2s(1) ) / x2s(2)
%
%
%                                    degp+1
%                                    __
%                                    \
%                     f(x) = g(s)  = /  cp(k)  T   (s)
%                                    --         k-1
%                                    k=1
%
%
%      x        the abscissa value at which the function defined by the input
%               expansion and its integral are to be evaluated.
%
%               [1,1] = size(x); double = class(x)
%
%               Normally `x' should lie in the closed interval
%
%                  [ x2s(1)-x2s(2),  x2s(1)+x2s(2) ]
%
%               See the -Restrictions section below.
%
%   the call:
%
%      [p, itgrlp] = cspice_chbigr( degp, cp, x2s, x )
%
%   returns:
%
%      p,
%      itgrlp   Define `s' and f(x) as above in the description of the
%               input argument `x2s'.
%
%               [1,1] = size(p); double = class(p)
%               [1,1] = size(itgrlp); double = class(itgrlp)
%
%               Then `p' is f(x), and `itgrlp' is an indefinite integral of
%               f(x), evaluated at `x'.
%
%               The indefinite integral satisfies
%
%                  d(itgrlp)/dx     = f(x)
%
%               and
%
%                  itgrlp( x2s(1) ) = 0
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
%   1) Let the domain of a polynomial to be evaluated be the
%      closed interval
%
%         [20, 30]
%
%      Let the input expansion represent the polynomial
%
%                           6
%         f(x)  = g(s) = 5*s
%
%      where
%
%         s     = (x - 20)/10
%
%      Let F(x) be an indefinite integral of f(x) such that
%
%         F(20) = 0
%
%      Evaluate
%
%         f(30) and F(30)
%
%
%      Example code begins here.
%
%
%      function chbigr_ex1()
%
%         %
%         % Let our domain be the interval [10, 30].
%         %
%         x2s = [ 20.0, 10.0 ]';
%
%         %
%         % Assign the expansion coefficients.
%         %
%         degp = 5;
%
%         cp   = [ 0.0, 3.75, 0.0, 1.875, 0.0, 0.375 ]';
%
%         %
%         % Evaluate the function and its integral at x = 30.
%         %
%         x = 30.0;
%
%         [p, itgrlp] = cspice_chbigr( degp, cp, x2s, x );
%
%         %
%         % We make the change of variables
%         %
%         %    s(x) = (1/10) * ( x - 20 )
%         %
%         % The expansion represents the polynomial
%         %
%         %                     5
%         %    f(x) = g(s) = 6*s
%         %
%         % An indefinite integral of the expansion is
%         %
%         %                                6
%         %    F(x) = G(s) * dx/ds = 10 * s
%         %
%         % where `G' is defined on the interval [-1, 1]. The result
%         % should be (due to the change of variables)
%         %
%         %      (G(1)  - G(0) ) * dx/ds
%         %
%         %    = (F(30) - F(20)) * 10
%         %
%         %    = 10
%         %
%         % The value of the expansion at `x' should be
%         %
%         %    f(30) = g(1) = 6
%         %
%         fprintf( 'ITGRLP = %f\n', itgrlp )
%         fprintf( 'P      = %f\n', p )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      ITGRLP = 10.000000
%      P      = 6.000000
%
%
%-Particulars
%
%   Let
%
%      T ,  n = 0, ...
%       n
%
%   represent the nth Chebyshev polynomial of the first kind:
%
%      T (x) = cos( n*arccos(x) )
%       n
%
%   The input coefficients represent the Chebyshev expansion
%
%                     degp+1
%                     __
%                     \
%      f(x) = g(s)  = /  cp(k)  T   (s)
%                     --         k-1
%                     k=1
%
%   where
%
%      s = ( x - x2s(1) ) / x2s(2)
%
%   This routine evaluates and returns the value at `x' of an
%   indefinite integral f(x), where
%
%      df(x)/dx    = f(x)  for all `x' in
%                          [x2s(1)-x2s(2), x2s(1)+x2s(2)]
%
%      f( x2s(1) ) = 0
%
%   The value at `x' of the input expansion
%
%      f(x)
%
%   is returned as well.
%
%   Note that numerical problems may result from applying this
%   routine to abscissa values outside of the interval defined
%   by the input parameters x2s(*). See the -Restrictions section.
%
%   To evaluate Chebyshev expansions and their derivatives, use the
%   Mice routines cspice_chbint or cspice_chbder.
%
%   This routine supports the SPICELIB SPK type 20 and PCK type 20
%   evaluators SPKE20 and PCKE20.
%
%-Exceptions
%
%   1)  If the input degree is negative, an error is signaled by the
%       Mice interface.
%
%   2)  If the input interval radius is non-positive, the error
%       SPICE(INVALIDRADIUS) is signaled by a routine in the call tree
%       of this routine.
%
%   3)  If any of the input arguments, `degp', `cp', `x2s' or `x', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `degp', `cp', `x2s' or `x', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   5)  If the number of elements in `cp' is less than degp+1, an error
%       is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The value (x-x2s(1)) / x2s(2) normally should lie within the
%       interval -1:1 inclusive, that is, the closed interval
%       [-1, 1]. Chebyshev polynomials increase rapidly in magnitude
%       as a function of distance of abscissa values from this
%       interval.
%
%       In typical SPICE applications, where the input expansion
%       represents position, velocity, or orientation, abscissa
%       values that map to points outside of [-1, 1] due to round-off
%       error will not cause numeric exceptions.
%
%   2)  No checks for floating point overflow are performed.
%
%   3)  Significant accumulated round-off error can occur for input
%       expansions of excessively high degree. This routine imposes
%       no limits on the degree of the input expansion; users must
%       verify that the requested computation provides appropriate
%       accuracy.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  W. Press, B. Flannery, S. Teukolsky and W. Vetterling,
%        "Numerical Recipes -- The Art of Scientific Computing,"
%        chapter 5.4, "Recurrence Relations and Clenshaw's Recurrence
%        Formula," p 161, Cambridge University Press, 1986.
%
%   [2]  "Chebyshev polynomials," Wikipedia, The Free Encyclopedia.
%        Retrieved 01:23, November 23, 2013, from
%        http://en.wikipedia.org/w/index.php?title=
%        Chebyshev_polynomials&oldid=574881046
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
%   integral of chebyshev_polynomial_expansion
%   integrate chebyshev_polynomial_expansion
%
%-&
function [p, itgrlp] = cspice_chbigr( degp, cp, x2s, x )

   switch nargin
      case 4

         degp = zzmice_int(degp);
         cp   = zzmice_dp(cp);
         x2s  = zzmice_dp(x2s);
         x    = zzmice_dp(x);

      otherwise

         error ( [ 'Usage: [p, itgrlp] = '                                 ...
                   'cspice_chbigr( degp, cp(degp+1), x2s(2), x )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [p, itgrlp] = mice('chbigr_c', degp, cp, x2s, x);
   catch spiceerr
      rethrow(spiceerr)
   end
