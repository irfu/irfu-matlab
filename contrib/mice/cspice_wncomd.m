%-Abstract
%
%   CSPICE_WNCOMD returns the complement of a double precision
%   window with respect to a specified interval.
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
%      left,
%      right    values defining the left and right endpoints of the complement
%               interval.
%
%               [1,1] = size(right); double = class(right)
%
%      window   SPICE window containing zero or more intervals
%
%               [2m,1] = size(window); double = class(window)
%
%   the call:
%
%      [result] = cspice_wncomd( left, right, window )
%
%   returns:
%
%      result   SPICE window containing the complement of `window' with
%               respect to the interval `left' to `right'
%
%               [2n,1] = size(result); double = class(result)
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
%   1) Given a double precision window, containing the following three
%      intervals:
%
%         [ 1.0, 3.0 ], [ 7.0, 11.0 ], [ 23.0, 27.0 ]
%
%      compute its complement with respect to the intervals [2.0, 20.0]
%      and [0.0, 100.0]
%
%      Example code begins here.
%
%
%      function wncomd_ex1()
%
%         %
%         % Let `window' contain the intervals
%         %
%         window = [ [ 1; 3 ];  [ 7; 11 ];  [ 23; 27 ];  ];
%
%         %
%         % The floating point complement of window with respect
%         % to [2,20]
%         %
%         [win1] = cspice_wncomd( 2, 20, window );
%
%         fprintf( 'Complement window with respect to [2.0, 20.0]\n' );
%         for i=1:cspice_wncard(win1)
%
%            [left, right] = cspice_wnfetd( win1, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         %
%         % The complement with respect to [ 0, 100 ]
%         %
%         [win2] = cspice_wncomd( 0, 100, window );
%
%         fprintf( '\nComplement window with respect to [0.0, 100.0]\n' );
%         for i=1:cspice_wncard(win2)
%
%            [left, right] = cspice_wnfetd( win2, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Complement window with respect to [2.0, 20.0]
%              3.000000         7.000000
%             11.000000        20.000000
%
%      Complement window with respect to [0.0, 100.0]
%              0.000000         1.000000
%              3.000000         7.000000
%             11.000000        23.000000
%             27.000000       100.000000
%
%
%-Particulars
%
%   Mathematically, the complement of a window contains those
%   points that are not contained in the window. That is, the
%   complement of the set of closed intervals
%
%      [ a(0), b(0) ], [ a(1), b(1) ], ..., [ a(n), b(n) ]
%
%   is the set of open intervals
%
%      ( -inf, a(0) ), ( b(0), a(1) ), ..., ( b(n), +inf )
%
%   Because Matlab offers no satisfactory representation of
%   infinity, we must take the complement with respect to a
%   finite interval.
%
%   In addition, Matlab offers no satisfactory floating point
%   representation of open intervals. Therefore, the complement
%   of a floating point window is closure of the set theoretical
%   complement. In short, the floating point complement of the
%   window
%
%      [ a(0), b(0) ], [ a(1), b(1) ], ..., [ a(n), b(n) ]
%
%   with respect to the interval from left to right is the
%   intersection of the windows
%
%      ( -inf, a(0) ), ( b(0), a(1) ), ..., ( b(n), +inf )
%
%   and
%
%      [ left, right ]
%
%   Note that floating point intervals of measure zero (singleton
%   intervals) in the original window are replaced by gaps of
%   measure zero, which are filled. Thus, complementing a floating
%   point window twice does not necessarily yield the original window.
%
%-Exceptions
%
%   1)  If `left' is greater than `right', the error SPICE(BADENDPOINTS)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  The cardinality of the input `window' must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   3)  If any of the input arguments, `left', `right' or `window', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `left', `right' or `window', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and reformatted example's output.
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
%   -Mice Version 1.0.1, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 24-JUL-2007 (EDW)
%
%-Index_Entries
%
%   complement a d.p. window
%
%-&

function [result] = cspice_wncomd( left, right, window)

   switch nargin

      case 3

         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [result] = cspice_wncomd( left, right, window)' )

   end

   %
   % Call the windows routine, add to 'a' and 'b' the space needed for
   % the control segments.
   %
   try
      [result] = mice('wncomd_c', left, right, [zeros(6,1); window] );
   catch spiceerr
      rethrow(spiceerr)
   end


