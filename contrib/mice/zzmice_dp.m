%-Abstract
%
%   ZZMICE_DP converts a numeric input to double precision format.
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
%      x        an input numeric to convert to double precision.
%
%               [n,m] = size(x); numeric = class(x)
%
%               `x' may have any size or shape, but must be of class
%               numeric.
%
%      nanok    an optional input logical flag indicating whether the input
%               argument `x' is allowed to be NaN.
%
%               [1,1] = size(nanok); logical = class(nanok)
%
%               If not provided, this routine defaults `nanok' to false.
%
%   the call:
%
%      [y] = zzmice_dp(x)
%
%         or
%
%      [y] = zzmice_dp(x, nanok)
%
%   returns:
%
%      y        the double precision representation of `x'.
%
%               [n,m] = size(y); numeric = class(y)
%
%               `y' returns with the same size and shape of `x' and
%               class double.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   None.
%
%-Particulars
%
%   This routine exists to support the NAIF MATLAB-CSPICE interface.
%
%-Exceptions
%
%   1)  If the number of arguments passed to the function is not 1 or 2,
%       the error MICE(USAGE) is signaled.
%
%   2)  If the argument `x' is not of numeric type, the error MICE(BADARG) is
%       signaled.
%
%   3)  If the argument `x' is not finite, i.e. it is "NaN", and the argument
%       `nanok' is not set to true, the error MICE(NOTFINITE) is signaled.
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
%   None.
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
%   -Mice Version 1.2.0, 30-JUN-2021 (JDR)
%
%       Added argument "nanok" in order to support non-finite numeric inputs.
%
%       Edited the header to comply with NAIF standard.
%
%   -Mice Version 1.1.1, 12-FEB-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.1.0, 27-JUL-2009 (EDW)
%
%       Added value check on 'nargin'. Incorrect input argument type/form
%       error tag changed from "MICE(BADVAL)" to "MICE(BADARG)."
%
%   -Mice Version 1.0.1, 30-DEC-2008 (EDW)
%
%       Function ensures all input values as finite.
%
%       Corrected misspellings.
%
%   -Mice Version 1.0.0, 30-JAN-2006 (EDW)
%
%-Index_Entries
%
%   input converstion to double precision
%   verification of numeric input arguments
%
%-&

function [y] = zzmice_dp( x, nanok )

   %
   % One argument: we do not have nanok flag. Set the default to
   % false.
   %
   if( isequal(nargin,1) )

      nanok = false;

   elseif( ~isequal(nargin,2) )

      error( 'MICE(USAGE): [_y_] = zzmice_dp( _x_, [nanok] )' )

   end

   %
   % Check if the input is numeric.
   %
   if( ~isnumeric(x) )

      error( [ 'MICE(BADARG): Improper type of input argument passed to '  ...
               'function. Value or values expected as double precision '   ...
               'or integer.' ] )

   end

   %
   % Check if we have NaN values in the input `x', if these are not allowed.
   %
   if( ~nanok && ~all( isfinite( x(:) ) ) )

      error( [ 'MICE(NOTFINITE): Improper type of input argument passed '  ...
               'to function. Value or values expected as finite double '   ...
               'precision or integer.' ] )

   end

   %
   % We know we have a numeric input. Convert it to double.
   %
   y = double(x);
