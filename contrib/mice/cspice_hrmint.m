%-Abstract
%
%   CSPICE_HRMINT evaluates a Hermite interpolating polynomial at a specified
%   abscissa value.
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
%               The arrays `xvals' and `yvals' contain `n' and 2*n elements
%               respectively.
%
%      xvals    an array of length `n' containing abscissa values.
%
%               [n,1] = size(xvals); double = class(xvals)
%
%      yvals    an array of length 2*n containing ordinate and derivative
%               values for each point in the domain defined by `xvals'.
%
%               [2*n,1] = size(yvals); double = class(yvals)
%
%               The elements
%
%                  yvals( 2*i - 1 )
%                  yvals( 2*i     )
%
%               give the value and first derivative of the output
%               polynomial at the abscissa value
%
%                  xvals(i)
%
%               where `i' ranges from 1 to `n'.
%
%      x        the abscissa value at which the interpolating polynomial and
%               its derivative are to be evaluated.
%
%               [1,1] = size(x); double = class(x)
%
%   the call:
%
%      [f, df] = cspice_hrmint( n, xvals, yvals, x )
%
%   returns:
%
%      f,
%      df       the value and derivative at `x' of the unique polynomial
%               of degree 2n-1 that fits the points and derivatives defined
%               by `xvals' and `yvals'.
%
%               [1,1] = size(f); double = class(f)
%               [1,1] = size(df); double = class(df)
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
%   1) Fit a 7th degree polynomial through the points ( x, y, y' )
%
%         ( -1,      6,       3 )
%         (  0,      5,       0 )
%         (  3,   2210,    5115 )
%         (  5,  78180,  109395 )
%
%      and evaluate this polynomial at x = 2.
%
%      The returned value should be 141.0, and the returned
%      derivative value should be 456.0, since the unique 7th degree
%      polynomial that fits these constraints is
%
%                   7       2
%         f(x)  =  x   +  2x  + 5
%
%
%      Example code begins here.
%
%
%      function hrmint_ex1()
%
%         n     =   4;
%
%         xvals = [ -1.0, 0.0, 3.0, 5.0 ]';
%
%         yvals = [ 6.0, 3.0, 5.0, 0.0, 2210.0, 5115.0, 78180.0, 109395.0 ]';
%
%         [answer, deriv] = cspice_hrmint( n, xvals, yvals, 2.0 );
%
%         fprintf( 'ANSWER = %f\n', answer )
%         fprintf( 'DERIV  = %f\n', deriv )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      ANSWER = 141.000000
%      DERIV  = 456.000000
%
%
%-Particulars
%
%   Users of this routine must choose the number of points to use
%   in their interpolation method. The authors of Reference [1] have
%   this to say on the topic:
%
%      Unless there is solid evidence that the interpolating function
%      is close in form to the true function f, it is a good idea to
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
%      function, will typically go berserk when the argument x is
%      outside the range of tabulated values by more than the typical
%      spacing of tabulated points.
%
%-Exceptions
%
%   1)  If two input abscissas are equal, the error SPICE(DIVIDEBYZERO)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If `n' is less than 1, the error SPICE(INVALIDSIZE) is
%       signaled by a routine in the call tree of this routine.
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
%   6)  If the number of elements in `xvals' is less than `n', an error
%       is signaled by the Mice interface.
%
%   7)  If the number of elements in `yvals' is less than 2*n, an error
%       is signaled by the Mice interface.
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
%   [2]  S. Conte and C. de Boor, "Elementary Numerical Analysis -- An
%        Algorithmic Approach," 3rd Edition, p 64, McGraw-Hill, 1980.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 01-JUL-2021 (JDR)
%
%-Index_Entries
%
%   interpolate function using Hermite polynomial
%   Hermite interpolation
%
%-&
function [f, df] = cspice_hrmint( n, xvals, yvals, x )

   switch nargin
      case 4

         n     = zzmice_int(n);
         xvals = zzmice_dp(xvals);
         yvals = zzmice_dp(yvals);
         x     = zzmice_dp(x);

      otherwise

         error ( [ 'Usage: [f, df] = '                                     ...
                   'cspice_hrmint( n, xvals(n), yvals(n+1), x )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [f, df] = mice('hrmint_c', n, xvals, yvals, x);
   catch spiceerr
      rethrow(spiceerr)
   end
