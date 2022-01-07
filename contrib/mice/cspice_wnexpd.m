%-Abstract
%
%   CSPICE_WNEXPD expands each of the intervals of a double
%   precision window.
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
%      left     amount to subtract from each left interval endpoint in the
%               input 'window'
%
%               [1,1] = size(left); double = class(left)
%
%      right    amount to add to each right interval endpoint in the input
%               'window'
%
%               [1,1] = size(right); double = class(right)
%
%      window   SPICE window containing zero or more intervals.
%
%               [2m,1] = size(window); double = class(window)
%
%   the call:
%
%      window_f = cspice_wnexpd( left, right, window)
%
%   returns:
%
%      window_f   SPICE window containing zero or more intervals, representing
%                 the original 'window' with each of its intervals expanded by
%                 'left' units on the left and 'right' units on the right
%
%                 [2n,1] = size(window_f); double = class(window_f)
%
%                 'window_f' can overwrite 'window'.
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
%   1) Given a double precision window, containing the following four
%      intervals:
%
%         [ 1.0, 3.0 ], [ 7.0, 11.0 ], [ 23.0, 27.0 ], [ 29.0, 29.0 ]
%
%      expand each interval by 2.0 units on the left and 1.0 on the
%      right endpoints, then by -2.0 units on the left and 2.0 units on
%      the right, and finally -2.0 units on the left and -1.0 units on
%      the right.
%
%      Example code begins here.
%
%
%      function wnexpd_ex1()
%
%         %
%         % Let 'window' contain the intervals
%         %
%         window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; [ 29; 29 ]; ];
%
%         %
%         % Apply the following series of calls
%         %
%         window = cspice_wnexpd(  2,  1, window );
%         fprintf( '1: Expanded window by  2 (left) and  1 (right)\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wnexpd( -2,  2, window );
%         fprintf( '2: Expanded window by -2 (left) and  2 (right)\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wnexpd( -2, -1, window );
%         fprintf( '3: Expanded window by -2 (left) and -1 (right)\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      1: Expanded window by  2 (left) and  1 (right)
%             -1.000000         4.000000
%              5.000000        12.000000
%             21.000000        30.000000
%      2: Expanded window by -2 (left) and  2 (right)
%              1.000000         6.000000
%              7.000000        14.000000
%             23.000000        32.000000
%      3: Expanded window by -2 (left) and -1 (right)
%              3.000000         5.000000
%              9.000000        13.000000
%             25.000000        31.000000
%
%
%      Note that intervals may be "expanded" by negative amounts.
%      In the example above, the second call shifts each interval to
%      the right, while the third call undoes the effect of the first
%      call (without restoring the merged intervals).
%
%      Note also that the third call is exactly equivalent to the
%      call:
%
%         cspice_wncond( 2, 1, window )
%
%-Particulars
%
%   This function expands (lengthens) each of the intervals in
%   the input window. The adjustments are not necessarily symmetric.
%   That is, left units are subtracted from the left endpoint of
%   each interval, and right units are added to the right endpoint
%   of each interval, where left and right may be different.
%
%   Intervals are merged when expansion causes them to overlap.
%
%-Exceptions
%
%   1)  The cardinality of the input `window' must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   2)  If any of the input arguments, `left', `right' or `window', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `left', `right' or `window', is
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
%   expand the intervals of a d.p. window
%
%-&

function [window_f] = cspice_wnexpd( left, right, window)

   switch nargin

      case 3

         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [window_f] = cspice_wnexpd( left, right, window)' )

   end

   %
   % Please note, this call does not require addition of space for the window
   % control segment as needed by other windows interfaces. The interface
   % copies the data in 'window' to a work variable rather than directly
   % pass 'window' into a CSPICE call.
   %
   try
      [window_f] = mice('wnexpd_c', left, right, window );
   catch spiceerr
      rethrow(spiceerr)
   end


