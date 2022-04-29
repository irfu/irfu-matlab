%-Abstract
%
%   CSPICE_WNFLTD filters (removes) small intervals from a
%   double precision window.
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
%      smlint   scalar limiting measure of the small intervals to be
%               filtered. Intervals of measure less than or equal to
%                `smlint' are removed from the window.
%
%               [1,1] = size(smlint); double = class(smlint)
%
%      window_i SPICE window containing zero or more intervals.
%
%               [2m,1] = size(window_i); double = class(window_i)
%
%   the call:
%
%      [window_f] = cspice_wnfltd( smlint, window_i )
%
%   returns:
%
%      window_f SPICE window containing zero or more intervals representing
%               the original `window_i', after `smlint' intervals have been
%               removed.
%
%               [2n,1] = size(window_f); double = class(window_f)
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
%      remove any intervals smaller than or equal to 0.0. From the
%      resulting window, remove any intervals smaller than or equal to
%      2.0; then those of measures less than or equal to 3; and,
%      finally, those smaller than or equal to 4.0.
%
%      Example code begins here.
%
%
%      function wnfltd_ex1()
%
%         %
%         % Let `window' contain the intervals
%         %
%         window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; [ 29; 29 ]; ];
%
%         %
%         % Apply the following series of calls
%         %
%         window = cspice_wnfltd(  0, window );
%         fprintf( '1: Removing intervals smaller than or equal to 0.0\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wnfltd(  2, window );
%         fprintf( '2: Removing intervals smaller than or equal to 2.0\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wnfltd(  3, window );
%         fprintf( '3: Removing intervals smaller than or equal to 3.0\n' );
%         for i=1:cspice_wncard(window)
%
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%         end
%
%         window = cspice_wnfltd(  4, window );
%         fprintf( '4: Removing intervals smaller than or equal to 4.0\n' );
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
%      1: Removing intervals smaller than or equal to 0.0
%              1.000000         3.000000
%              7.000000        11.000000
%             23.000000        27.000000
%      2: Removing intervals smaller than or equal to 2.0
%              7.000000        11.000000
%             23.000000        27.000000
%      3: Removing intervals smaller than or equal to 3.0
%              7.000000        11.000000
%             23.000000        27.000000
%      4: Removing intervals smaller than or equal to 4.0
%
%
%      Note that the last result (4) represents an empty set.
%
%-Particulars
%
%   This routine removes from the input window every interval with
%   measure less than or equal to the limiting measure `smlint'.
%
%-Exceptions
%
%   1)  The cardinality of the input `window_i' must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   2)  If `smlint' is less than or equal to zero, this routine has
%       no effect on the window.
%
%   3)  If any of the input arguments, `smlint' or `window_i', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `smlint' or `window_i', is not
%       of the expected type, or it does not have the expected
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Changed argument names "window" and "sml" to "window_i" and
%       "smlint" for consistency with other routines.
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
%   -Mice Version 1.0.1, 11-JUN-2013 (EDW)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 24-JUL-2007 (EDW)
%
%-Index_Entries
%
%   filter small intervals from a d.p. window
%
%-&

function [window_f] = cspice_wnfltd( smlint, window_i )

   switch nargin

      case 2

         smlint   = zzmice_dp(smlint);
         window_i = zzmice_win(window_i);

      otherwise

         error ( 'Usage: [window_f] = cspice_wnfltd( smlint, window_i )' )

   end

%
% Please note, this call does not require addition of space for the window
% control segment as needed by other windows interfaces. The interface
% copies the data in `window_i' to a work variable rather than directly
% pass `window_i' into a CSPICE call.
%
   try
      [window_f] = mice( 'wnfltd_c', smlint, window_i );
   catch spiceerr
      rethrow(spiceerr)
   end



