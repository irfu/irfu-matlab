%-Abstract
%
%   CSPICE_WNINSD inserts an interval into a double precision window.
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
%      right      left and right endpoints of the interval to insert.
%
%                 [1,1] = size(right); double = class(right)
%
%      window_i   optional input SPICE window containing zero or more
%                 intervals. Inclusion of this window argument results in
%                 an output window consisting of a union of the data in
%                 `window_i' and the window defined as [left,right].
%
%                 [2m,1] = size(window_i); double = class(window_i)
%
%   the call:
%
%      [window_f] = cspice_wninsd( left, right )
%
%         or
%
%      [window_f] = cspice_wninsd( left, right, window_i )
%
%   returns:
%
%      window_f   SPICE Window consisting of either the interval [left,right]
%                 or the window union of `window_i' and [left,right]
%
%                 [2n,1] = size(window); double = class(window)
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
%   1) Create an empty window structure based on a cell containing a
%      double precision 6-vector and insert the following three
%      intervals into that window:
%
%         [ 1, 3 ]  [ 7, 11 ]  [ 23, 27 ]
%
%      Then, insert the following intervals into the window, printing
%      the resulting window after each insertion:
%
%         [ 5, 5 ],  [ 4, 8 ],  [ 0, 30 ]  and [ 31, 32 ]
%
%      Example code begins here.
%
%
%      function wninsd_ex1()
%
%         %
%         % Create an empty `window' to hold the data.
%         %
%         window = zeros(0,1);
%
%         %  Let `window' contain the intervals
%         %
%         %  [ 1, 3 ]  [ 7, 11 ]  [ 23, 27 ]
%         %
%         data = [ 1, 3, 7, 11, 23, 27];
%
%         %
%         % Add the data to `window'.
%         %
%         for i=1:numel(data)/2
%            window = cspice_wninsd( data(2*i - 1), data(2*i), window );
%         end
%
%         %
%         % Note, the direct assignment:
%         %
%         %   window = [ [1; 3]; [7; 11]; [23; 27] ];
%         %
%         % will also perform the assignment of `data' to `window' but
%         % NAIF recommends using cspice_wninsd when possible.
%         %
%         % Print out the original window.
%         %
%         fprintf( 'Original window:\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         %
%         % Perform a series of cspice_wninsd calls, adding intervals
%         % to `window':
%         %
%         window = cspice_wninsd( 5,5, window );
%         fprintf( '1: Window after inserting the interval [ 5, 5 ]:\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wninsd( 4,8, window );
%         fprintf( '2: Window after inserting the interval [ 4, 8 ]:\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wninsd( 0,30, window );
%         fprintf( '3: Window after inserting the interval [ 0, 30 ]:\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wninsd( 31,32, window );
%         fprintf( '4: Window after inserting the interval [ 31, 32 ]:\n' );
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
%      Original window:
%              1.000000         3.000000
%              7.000000        11.000000
%             23.000000        27.000000
%      1: Window after inserting the interval [ 5, 5 ]:
%              1.000000         3.000000
%              5.000000         5.000000
%              7.000000        11.000000
%             23.000000        27.000000
%      2: Window after inserting the interval [ 4, 8 ]:
%              1.000000         3.000000
%              4.000000        11.000000
%             23.000000        27.000000
%      3: Window after inserting the interval [ 0, 30 ]:
%              0.000000        30.000000
%      4: Window after inserting the interval [ 31, 32 ]:
%              0.000000        30.000000
%             31.000000        32.000000
%
%
%-Particulars
%
%   This routine inserts the interval from `left' to `right' into the
%   input window. If the new interval overlaps any of the intervals
%   in the window, the intervals are merged. Thus, the cardinality
%   of the input window can actually decrease as the result of an
%   insertion. However, because inserting an interval that is
%   disjoint from the other intervals in the window can increase the
%   cardinality of the window, the routine signals an error.
%
%-Exceptions
%
%   1)  If `left' is greater than `right', the error SPICE(BADENDPOINTS)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  The cardinality of the input `window_i' must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   3)  If any of the input arguments, `left', `right' or `window_i',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   4)  If any of the input arguments, `left', `right' or `window_i',
%       is not of the expected type, or it does not have the expected
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
%       Changed output argument name "window" to "window_f" for consistency
%       with other routines.
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
%       Removed irrelevant information related to other unary window
%       routines from -Particulars section.
%
%   -Mice Version 1.0.2, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 21-OCT-2008 (EDW)
%
%       Edited example to demonstrate creation of, and loading of,
%       an empty window.
%
%   -Mice Version 1.0.0, 10-JUL-2007 (EDW)
%
%-Index_Entries
%
%   insert an interval into a d.p. window
%
%-&

function [window_f] = cspice_wninsd( left, right, window_i )

   switch nargin
      case 2

         left   = zzmice_dp(left);
         right  = zzmice_dp(right);

      case 3

         left     = zzmice_dp(left);
         right    = zzmice_dp(right);
         window_i = zzmice_win(window_i);

      otherwise

         error ( ['Usage: [window_f] = cspice_wninsd( left, right, ', ...
                                                     '[window_i] )'] )

   end

   %
   % The call passed either two or three arguments. Branch accordingly.
   %
   if ( nargin == 2 )

      try
         [window_f] = mice( 'wninsd_c', left, right );
      catch
         rethrow(lasterror)
      end

   else

      try
         [window_f] = mice( 'wninsd_c', left, right );
      catch
         rethrow(lasterror)
      end

      %
      % Perform a window union with the window defined by [left,right]
      % and the input window `window_i'.
      %
      window_f = cspice_wnunid( window_f, window_i );

   end
