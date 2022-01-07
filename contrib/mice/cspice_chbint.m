%-Abstract
%
%   CSPICE_CHBINT returns the value of a polynomial and its derivative,
%   evaluated at the input `x', using the coefficients of the Chebyshev
%   expansion of the polynomial.
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
%      cp       an array of coefficients of a polynomial with respect to the
%               Chebyshev basis.
%
%               [degp+1,1] = size(cp); double = class(cp)
%
%               The polynomial to be evaluated is assumed to be of the
%               form:
%
%                  cp(degp+1)*T(degp,s) + cp(degp)*T(degp-1,s) + ...
%
%                                       + cp(2)*T(1,s) + cp(1)*T(0,s)
%
%               where T(i,s) is the i'th Chebyshev polynomial
%               evaluated at a number `s' whose double precision
%               value lies between -1 and 1. The value of `s' is
%               computed from the input variables x2s(1), x2s(2) and `x'.
%
%      degp     the degree of the Chebyshev polynomial to be evaluated.
%
%               [1,1] = size(degp); int32 = class(degp)
%
%      x2s      an array of two parameters.
%
%               [2,1] = size(x2s); double = class(x2s)
%
%               These parameters are used to transform the domain of the
%               input variable `x' into the standard domain of the Chebyshev
%               polynomial. x2s(1) should be a reference point in the domain
%               of `x'; x2s(2) should be the radius by which points are
%               allowed to deviate from the reference point and while
%               remaining within the domain of `x'. The value of `x' is
%               transformed into the value `s' given by
%
%                  s = ( x - x2s(1) ) / x2s(2)
%
%               Typically x2s(1) is the midpoint of the interval over
%               which `x' is allowed to vary and x2s(2) is the radius of
%               the interval.
%
%               The main reason for doing this is that a Chebyshev
%               expansion is usually fit to data over a span
%               from A to B where A and B are not -1 and 1
%               respectively. Thus to get the "best fit" the
%               data was transformed to the interval [-1,1] and
%               coefficients generated. These coefficients are
%               not rescaled to the interval of the data so that
%               the numerical "robustness" of the Chebyshev fit will
%               not be lost. Consequently, when the "best fitting"
%               polynomial needs to be evaluated at an intermediate
%               point, the point of evaluation must be transformed
%               in the same way that the generating points were
%               transformed.
%
%      x        the value for which the polynomial is to be evaluated.
%
%               [1,1] = size(x); double = class(x)
%
%   the call:
%
%      [p, dpdx] = cspice_chbint( cp, degp, x2s, x )
%
%   returns:
%
%      p        the value of the polynomial to be evaluated.
%
%               [1,1] = size(p); double = class(p)
%
%               It is given by
%
%                  cp(degp+1)*T(degp,s) + cp(degp)*T(degp-1,s) + ...
%
%                                       + cp(2)*T(1,s) + cp(1)*T(0,s)
%
%               where T(i,s) is the i'th Chebyshev polynomial
%               evaluated  at a number s = ( x - x2s(1) )/x2s(2)
%
%      dpdx     the value of the derivative of the polynomial at `x'.
%
%               [1,1] = size(dpdx); double = class(dpdx)
%
%               It is given by
%
%                  1/x2s(2) [    cp(degp+1)*T'(degp,s)
%
%                              + cp(degp)*T'(degp-1,s)
%
%                              + ...
%
%                              + cp(2)*T'(1,s)
%
%                              + cp(1)*T'(0,s) ]
%
%               where T(i,s) and T'(`i',S) are the i'th Chebyshev
%               polynomial and its derivative, respectively,
%               evaluated  at a number s = ( x - x2s(1) )/x2s(2)
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
%   1) Depending upon the user's needs, there are 3 routines
%      available for evaluating Chebyshev polynomials.
%
%         cspice_chbval   for evaluating a Chebyshev polynomial when no
%                         derivatives are desired.
%
%         cspice_chbint   for evaluating a Chebyshev polynomial and its
%                         first derivative.
%
%         cspice_chbder   for evaluating a Chebyshev polynomial and a user
%                         or application dependent number of derivatives.
%
%      Of these 3 the one most commonly employed by SPICE software
%      is cspice_chbint as it is used to interpolate ephemeris state
%      vectors; this requires the evaluation of a polynomial
%      and its derivative. When no derivatives are desired one
%      should use cspice_chbval, or when more than one or an unknown
%      number of derivatives are desired one should use cspice_chbder.
%
%      The code example below illustrates how this routine might
%      be used to obtain points for plotting a polynomial
%      and its derivative.
%
%
%      Example code begins here.
%
%
%      function chbint_ex1()
%
%         %
%         % Set the coefficients of the polynomial and its
%         % transformation parameters
%         %
%         cp   = [1.0,3.0,0.5,1.0,0.5,-1.0,1.0]';
%         x2s  = [0.5,3.0]';
%
%         degp = 6;
%         x    = 1.0;
%
%         [p, dpdx] = cspice_chbint( cp, degp, x2s, x );
%
%         fprintf( 'Value of the polynomial at x=1:  %9.6f\n', p )
%         fprintf( '   First derivative at x=1    :  %9.6f\n', dpdx )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Value of the polynomial at x=1:  -0.340878
%         First derivative at x=1    :   0.382716
%
%
%-Particulars
%
%   This routine computes the value of a Chebyshev polynomial
%   expansion and the derivative of the expansion with respect to `x'.
%   The polynomial is given by
%
%      cp(degp+1)*T(degp,s) + cp(degp)*T(degp-1,s) + ...
%
%                           + cp(2)*T(1,s) + cp(1)*T(0,s)
%
%   where
%
%      s  =  ( x - x2s(1) ) / x2s(2)
%
%   and
%
%      T(i,s) is the i'th Chebyshev polynomial of the first kind
%      evaluated at `s'.
%
%-Exceptions
%
%   1)  If any of the input arguments, `cp', `degp', `x2s' or `x', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `cp', `degp', `x2s' or `x', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   3)  If the input degree is negative, an error is signaled by the
%       Mice interface.
%
%   4)  If the number of elements in `cp' is less than degp+1, an error
%       is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  One needs to be careful that the value
%
%          (x-x2s(1)) / x2s(2)
%
%       lies between -1 and 1. Otherwise, the routine may fail
%       spectacularly (for example with a floating point overflow).
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
%   [2]  T. Rivlin, "The Chebyshev Polynomials," Wisley, 1974.
%
%   [3]  R. Weast and S. Selby, "CRC Handbook of Tables for
%        Mathematics," 4th Edition, CRC Press, 1976.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 01-NOV-2021 (JDR)
%
%-Index_Entries
%
%   interpolate a chebyshev expansion
%
%-&
function [p, dpdx] = cspice_chbint( cp, degp, x2s, x )

   switch nargin
      case 4

         cp   = zzmice_dp(cp);
         degp = zzmice_int(degp);
         x2s  = zzmice_dp(x2s);
         x    = zzmice_dp(x);

      otherwise

         error ( [ 'Usage: [p, dpdx] = '                                  ...
                   'cspice_chbint( cp(degp+1), degp, x2s(2), x )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [p, dpdx] = mice('chbint_c', cp, degp, x2s, x);
   catch spiceerr
      rethrow(spiceerr)
   end
