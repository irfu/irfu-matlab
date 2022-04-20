%-Abstract
%
%   CSPICE_WNRELD compares two double precision windows returning
%   a scalar boolean.
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
%      a        SPICE window containing zero or more intervals.
%
%               [2l,1] = size(a); double = class(a)
%
%      b        SPICE window containing zero or more intervals.
%
%               [2m,1] = size(b); double = class(b)
%
%      op       comparison operator, indicating the way to compare the input
%               windows.
%
%               [1,c1] = size(op); char = class(op)
%
%                  or
%
%               [1,1] = size(op); cell = class(op)
%
%               `op' may have any of the following values:
%
%                  Operator             Meaning
%                  --------  -----------------------------------------
%                    '='     a = b  is true if `a' and `b' are equal
%                                   (contain the same intervals).
%
%                    '<>'    a <> b is true if `a' and `b' are not
%                                   equal.
%
%                    '<='    a <= b is true if `a' is a subset of `b'.
%
%                    '<'     a < b  is true if `a' is a proper subset
%                                   of `b'.
%
%                    '>='    a >= b is true if `b' is a subset of `a'.
%
%                    '>'     a > b  is true if `b' is a proper subset
%                                   of `a'.
%
%   the call:
%
%      [wnreld] = cspice_wnreld( a, op, b )
%
%   returns:
%
%      wnreld   A scalar boolean with value of the comparison.
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
%   1) Given several double precision windows, compare them and output
%      the result of those comparisons.
%
%      Example code begins here.
%
%
%      function wnreld_ex1()
%
%         %
%         % Define the comparison operators.
%         %
%         ops = { '=', '<>', '<=', '<', '>=', '>' };
%
%         %
%         %  Let a contain the intervals
%         %
%         a = [ [ 1; 3 ];  [ 7; 11 ];  [ 23; 27 ] ];
%
%         %
%         %  Let b and c contain the intervals
%         %
%
%         b = [ [ 1; 2 ];  [  9; 9 ];  [ 24; 27 ] ];
%         c = b;
%
%         %
%         %  Let d contain the intervals
%         %
%         d = [ [ 5; 10 ];  [ 15; 25 ] ];
%
%         %
%         %  Finally, let e and f be empty windows (containing no intervals).
%         %
%         e = zeros(0,1);
%         f = e;
%
%         %
%         % Compare b and c, which contain the same intervals.
%         %
%         disp( 'b and c contain the same intervals:')
%         for i=1:numel(ops)
%            if ( cspice_wnreld( b, ops(i),  c ) )
%               fprintf( '  b %2s c  is True.\n', char(ops(i)))
%            else
%               fprintf( '  b %2s c  is False.\n', char(ops(i)))
%            end
%         end
%         disp( '' )
%
%         %
%         % Every point contained in b and c is also contained in a. Thus,
%         %
%         disp( 'Every point in b is also in contained in a:')
%         for i=1:numel(ops)
%            if ( cspice_wnreld( b, ops(i),  a ) )
%               fprintf( '  b %2s a  is True.\n', char(ops(i)))
%            else
%               fprintf( '  b %2s a  is False.\n', char(ops(i)))
%            end
%         end
%         disp( '' )
%
%         %
%         % Although a and d have points in common, neither contains the
%         % other. Thus
%         %
%         disp( 'a and d have points in common, neither contains the other:')
%         for i=1:numel(ops)
%            if ( cspice_wnreld( a, ops(i),  d ) )
%               fprintf( '  a %2s d  is True.\n', char(ops(i)))
%            else
%               fprintf( '  a %2s d  is False.\n', char(ops(i)))
%            end
%         end
%         disp( '' )
%
%         %
%         % In addition, any window is equal to itself, a subset of itself,
%         % and a superset of itself. Thus,
%         %
%         disp( 'A window compared to itself:')
%         for i=1:numel(ops)
%            if ( cspice_wnreld( a, ops(i),  a ) )
%               fprintf( '  a %2s a  is True.\n', char(ops(i)))
%            else
%               fprintf( '  a %2s a  is False.\n', char(ops(i)))
%            end
%         end
%         disp( '' )
%
%         %
%         % Finally, an empty window is a proper subset of any window
%         % except another empty window. Thus,
%         %
%         disp( 'A window compared to an empty window:')
%         for i=1:numel(ops)
%            if ( cspice_wnreld( a, ops(i),  e ) )
%               fprintf( '  a %2s e  is True.\n', char(ops(i)))
%            else
%               fprintf( '  a %2s e  is False.\n', char(ops(i)))
%            end
%         end
%         disp( '' )
%
%         disp( 'An empty window compared to another empty window:')
%         for i=1:numel(ops)
%            if ( cspice_wnreld( f, ops(i),  e ) )
%               fprintf( '  f %2s e  is True.\n', char(ops(i)))
%            else
%               fprintf( '  f %2s e  is False.\n', char(ops(i)))
%            end
%         end
%         disp( '' )
%
%         %
%         % is false.
%         %
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      b and c contain the same intervals:
%        b  = c  is True.
%        b <> c  is False.
%        b <= c  is True.
%        b  < c  is False.
%        b >= c  is True.
%        b  > c  is False.
%
%      Every point in b is also in contained in a:
%        b  = a  is False.
%        b <> a  is True.
%        b <= a  is True.
%        b  < a  is True.
%        b >= a  is False.
%        b  > a  is False.
%
%      a and d have points in common, neither contains the other:
%        a  = d  is False.
%        a <> d  is True.
%        a <= d  is False.
%        a  < d  is False.
%        a >= d  is False.
%        a  > d  is False.
%
%      A window compared to itself:
%        a  = a  is True.
%        a <> a  is False.
%        a <= a  is True.
%        a  < a  is False.
%        a >= a  is True.
%        a  > a  is False.
%
%      A window compared to an empty window:
%        a  = e  is False.
%        a <> e  is True.
%        a <= e  is False.
%        a  < e  is False.
%        a >= e  is True.
%        a  > e  is True.
%
%      An empty window compared to another empty window:
%        f  = e  is True.
%        f <> e  is False.
%        f <= e  is True.
%        f  < e  is False.
%        f >= e  is True.
%        f  > e  is False.
%
%
%-Particulars
%
%   This function returns True whenever the specified relationship
%   between the input windows `a' and `b' is satisfied. For example,
%   the expression
%
%      cspice_wnreld ( needed, '<=', avail )
%
%   is True whenever the window `needed' is a subset of the window
%   `avail'. One window is a subset of another window if each of
%   the intervals in the first window is included in one of the
%   intervals in the second window. In addition, the first window
%   is a proper subset of the second if the second window contains
%   at least one point not contained in the first window. (Thus,
%   '<' implies '<=', and '>' implies '>='.)
%
%   The following pairs of expressions are equivalent.
%
%      cspice_wnreld ( a, '>',  b )
%      cspice_wnreld ( b, '<',  a )
%
%      cspice_wnreld ( a, '>=', b )
%      cspice_wnreld ( b, '<=', a )
%
%-Exceptions
%
%   1)  If the relational operator is not recognized, the error
%       SPICE(INVALIDOPERATION) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  The cardinality of the input windows must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   3)  If any of the input arguments, `a', `op' or `b', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   4)  If any of the input arguments, `a', `op' or `b', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
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
%   WINDOWS.REQ
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
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and modified example code to generate
%       formatted output.
%
%       Added square brackets to the output argument in function
%       declaration, and renamed it to "wnreld".
%
%       Corrected error message format.
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
%   -Mice Version 1.0.2, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       "logical" call replaced with "zzmice_logical."
%
%       Corrected version ID in 23-JUL-2009 entry, "1.0.0" to "1.0.1."
%
%   -Mice Version 1.0.1, 23-JUL-2009 (EDW)
%
%       Replaced "boolean" calls with "logical" as "boolean" functionally
%       aliases "logical."
%
%   -Mice Version 1.0.0, 22-JUL-2007 (EDW)
%
%-Index_Entries
%
%   compare two d.p. windows
%
%-&

function [wnreld] = cspice_wnreld( a, op, b )

   switch nargin

      case 3

         a  = zzmice_win(a);
         op = zzmice_str(op);
         b  = zzmice_win(b);

      otherwise

         error( 'Usage: [wnreld] = cspice_wnreld( a, `op`, b )' )

      end

   try
      [wnreld] = mice( 'wnreld_c', [zeros(6,1); a], op, [zeros(6,1); b] );
      [wnreld] = zzmice_logical(wnreld);
   catch spiceerr
      rethrow(spiceerr)
   end



