%-Abstract
%
%   CSPICE_LGRIND evaluates a Lagrange interpolating polynomial, for a
%   specified set of coordinate pairs, at a specified abscissa value.
%   This routine returns both the value of the polynomial and its
%   derivative.
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
%      n        the number of points defining the polynomial.
%
%               [1,1] = size(n); int32 = class(n)
%
%               The arrays `xvals' and `yvals' contain `n' elements.
%
%      xvals,
%      yvals    arrays of abscissa and ordinate values that together
%               define `n' ordered pairs.
%
%               [n,1] = size(xvals); double = class(xvals)
%               [n,1] = size(yvals); double = class(yvals)
%
%               The set of points
%
%                  ( xvals(i), yvals(i) )
%
%               define the Lagrange polynomial used for
%               interpolation. The elements of `xvals' must be
%               distinct and in increasing order.
%
%      x        the abscissa value at which the interpolating polynomial is
%               to be evaluated.
%
%               [1,1] = size(x); double = class(x)
%
%   the call:
%
%      [p, dp] = cspice_lgrind( n, xvals, yvals, x )
%
%   returns:
%
%      p        the value at `x' of the unique polynomial of degree n-1 that
%               fits the points in the plane defined by `xvals' and `yvals'.
%
%               [1,1] = size(p); double = class(p)
%
%      dp       the derivative at `x' of the interpolating polynomial
%               described above.
%
%               [1,1] = size(dp); double = class(dp)
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
%   1) Fit a cubic polynomial through the points
%
%          ( -1, -2 )
%          (  0, -7 )
%          (  1, -8 )
%          (  3, 26 )
%
%      and evaluate this polynomial at x = 2.
%
%      The returned value of `p' should be 1.0, since the
%      unique cubic polynomial that fits these points is
%
%                     3      2
%         f(x)   =   x  + 2*x  - 4*x - 7
%
%      The returned value of `dp' should be 16.0, since the
%      derivative of f(x) is
%
%          '           2
%         f (x)  =  3*x  + 4*x - 4
%
%
%      Example code begins here.
%
%
%      function lgrind_ex1()
%
%         n     =   4;
%
%         xvals = [ -1, 0, 1, 3 ]';
%
%         yvals = [ -2, -7, -8, 26 ]';
%
%         [p, dp] = cspice_lgrind( n, xvals, yvals, 2.0 );
%
%         fprintf( 'P, DP = %f %f\n', p, dp )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      P, DP = 1.000000 16.000000
%
%
%-Particulars
%
%   Given a set of `n' distinct abscissa values and corresponding
%   ordinate values, there is a unique polynomial of degree n-1, often
%   called the "Lagrange polynomial", that fits the graph defined by
%   these values. The Lagrange polynomial can be used to interpolate
%   the value of a function at a specified point, given a discrete
%   set of values of the function.
%
%   Users of this routine must choose the number of points to use
%   in their interpolation method. The authors of Reference [1] have
%   this to say on the topic:
%
%      Unless there is solid evidence that the interpolating function
%      is close in form to the true function `f', it is a good idea to
%      be cautious about high-order interpolation. We
%      enthusiastically endorse interpolations with 3 or 4 points, we
%      are perhaps tolerant of 5 or 6; but we rarely go higher than
%      that unless there is quite rigorous monitoring of estimated
%      errors.
%
%   The same authors offer this warning on the use of the
%   interpolating function for extrapolation:
%
%      ...the dangers of extrapolation cannot be overemphasized:
%      An interpolating function, which is perforce an extrapolating
%      function, will typically go berserk when the argument `x' is
%      outside the range of tabulated values by more than the typical
%      spacing of tabulated points.
%
%-Exceptions
%
%   1)  If any two elements of the array `xvals' are equal, the error
%       SPICE(DIVIDEBYZERO) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If `n' is less than 1, an error is signaled by the Mice interface.
%
%   3)  This routine does not attempt to ward off or diagnose
%       arithmetic overflows.
%
%   4)  If any of the input arguments, `n', `xvals', `yvals' or `x',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   5)  If any of the input arguments, `n', `xvals', `yvals' or `x',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   6)  If the number of elements in `xvals' or `yvals' is less than `n',
%       an error is signaled by the Mice interface.
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
%   [1]  W. Press, B. Flannery, S. Teukolsky and W. Vetterling,
%        "Numerical Recipes -- The Art of Scientific Computing,"
%        chapters 3.0 and 3.1, Cambridge University Press, 1986.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 03-NOV-2021 (JDR)
%
%-Index_Entries
%
%   interpolate function using Lagrange polynomial
%   Lagrange interpolation
%
%-&
function [p, dp] = cspice_lgrind( n, xvals, yvals, x )

   switch nargin
      case 4

         n     = zzmice_int(n, [1, int32(inf)]);
         xvals = zzmice_dp(xvals);
         yvals = zzmice_dp(yvals);
         x     = zzmice_dp(x);

      otherwise

         error ( [ 'Usage: [p, dp] = '                                     ...
                   'cspice_lgrind( n, xvals(n), yvals(n), x )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [p, dp] = mice('lgrind_c', n, xvals, yvals, x);
   catch spiceerr
      rethrow(spiceerr)
   end
