%-Abstract
%
%   CSPICE_WNEXTD extracts the left or right endpoints from
%   a double precision window.
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
%      side     indicating whether to extract the left or right endpoints of
%               the intervals in 'window'
%
%               [1,1] = size(side); char = class(side)
%
%                  or
%
%               [1,1] = size(side); cell = class(side)
%
%               'L', 'l'       Left endpoints.
%               'R', 'r'       Right endpoints.
%
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%   the call:
%
%      window_f = cspice_wnextd( side, window)
%
%   returns:
%
%      window   SPICE window containing zero or more intervals, representing
%               the collection of singleton intervals containing either the
%               left or the right endpoints of the intervals in the original
%               'window'
%
%               [2n,1] = size(window); double = class(window)
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Given a double precision window, containing the following four
%      intervals:
%
%         [ 1.0, 3.0 ], [ 7.0, 11.0 ], [ 23.0, 27.0 ],  [29.0, 29.0]
%
%      extract all its left and all its right endpoints and output the
%      resulting singleton intervals.
%
%      Example code begins here.
%
%
%      function wnextd_ex1()
%
%         %
%         % Let 'window' contain the intervals
%         %
%         window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; [ 29; 29 ]; ];
%
%         %
%         %
%         %
%         left  = cspice_wnextd( 'L', window );
%         fprintf( '1: Singletons from extracting left endpoints\n' );
%         for i=1:cspice_wncard(left)
%
%            [point1, point2] = cspice_wnfetd( left, i );
%            fprintf( '%16.6f %16.6f\n', point1, point2  );
%
%         end
%
%         right = cspice_wnextd( 'R', window );
%         fprintf( '2: Singletons from extracting right endpoints\n' );
%         for i=1:cspice_wncard(right)
%
%            [point1, point2] = cspice_wnfetd( left, i );
%            fprintf( '%16.6f %16.6f\n', point1, point2  );
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      1: Singletons from extracting left endpoints
%              1.000000         1.000000
%              7.000000         7.000000
%             23.000000        23.000000
%             29.000000        29.000000
%      2: Singletons from extracting right endpoints
%              1.000000         1.000000
%              7.000000         7.000000
%             23.000000        23.000000
%             29.000000        29.000000
%
%
%   2) Repeat the example 1, using MATLAB native functionality.
%
%      Example code begins here.
%
%
%      function wnextd_ex2()
%
%         %
%         % Let 'window' contain the intervals
%         %
%         window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; [ 29; 29 ]; ];
%
%         %
%         % A similar operation using MATLAB native functionality,
%         % though the returned arrays are not SPICE windows identical
%         % to the cspice_wnextd result.
%         %
%         index_left = 1: 2 : numel(window);
%         index_right= index_left + 1;
%
%         left  = window( index_left );
%         fprintf( 'Left  endpoints: %8.2f %8.2f %8.2f %8.2f\n', left  );
%         right = window( index_right);
%         fprintf( 'Right endpoints: %8.2f %8.2f %8.2f %8.2f\n', right );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Left  endpoints:     1.00     7.00    23.00    29.00
%      Right endpoints:     3.00    11.00    27.00    29.00
%
%
%      Note that the returned arrays are not SPICE windows identical
%      to the cspice_wnextd result.
%
%-Particulars
%
%   This function returns a window composed of singleton intervals
%   containing one of the endpoints of the intervals in `window'.
%
%-Exceptions
%
%   1)  If the endpoint specification, `side', is not recognized, the
%       error SPICE(INVALIDENDPNTSPEC) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  The cardinality of the input `window' must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   3)  If any of the input arguments, `side' or `window', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `side' or `window', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%       example's problem statements and reformatted example's output.
%       Created a second example using the existing MATLAB native
%       functionality.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 12-MAR-2012 (EDW) (SCK)
%
%      Edited -I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 24-JUL-2007 (EDW)
%
%-Index_Entries
%
%   extract the endpoints from a d.p. window
%
%-&

function [window_f] = cspice_wnextd( side, window)

   switch nargin

      case 2

         side   = zzmice_str(side);
         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [window_f] = cspice_wnextd( side, window )' )

   end

%
% Please note, this call does not require addition of space for the window
% control segment as needed by other windows interfaces. The interface
% copies the data in 'window' to a work variable rather than directly
% pass 'window' into a CSPICE call.
%
   try
      [window_f] = mice('wnextd_c', side, window );
   catch spiceerr
      rethrow(spiceerr)
   end
