%-Abstract
%
%   CSPICE_WNVALD forms a valid double precision window from the contents
%   of an 2Mx1 array.
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
%      window_i   2M endpoints of (possibly unordered and non-disjoint)
%                 intervals.
%
%                 [2m,1] = size(window_i); double = class(window_i)
%
%   the call:
%
%      [window_f] = cspice_wnvald( window_i )
%
%   returns:
%
%      window_f   SPICE window containing the ordered union of the intervals
%                 in the input array `window_i'.
%
%                 [2n,1] = size(window_f); double = class(window_f)
%
%                 The `window_f' may overwrite the input `window_i'.
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
%   1) Define an array containing a set of unordered and possibly
%      overlapping intervals, and validate the array as a SPICE
%      window.
%
%
%      Example code begins here.
%
%
%      function wnvald_ex1()
%
%         %
%         % Local variables
%         %
%         windat = [ [0;   0]; ...
%                    [10; 12]; ...
%                    [2;   7]; ...
%                    [13; 15]; ...
%                    [1;   5]; ...
%                    [23; 29] ];
%
%         %
%         % Validate the input `windat' array as a SPICE window.
%         %
%         window = cspice_wnvald(windat);
%
%         fprintf( 'Current intervals: %d\n', cspice_wncard(window) )
%         fprintf( 'Maximum intervals: %d\n', numel(window)/2       )
%         fprintf( '\nIntervals\n\n' );
%
%         for i=1:cspice_wncard(window)
%            [left, right] = cspice_wnfetd( window, i );
%            fprintf( '%10.6f   %10.6f\n', left, right )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Current intervals: 5
%      Maximum intervals: 5
%
%      Intervals
%
%        0.000000     0.000000
%        1.000000     7.000000
%       10.000000    12.000000
%       13.000000    15.000000
%       23.000000    29.000000
%
%
%      Note that SPICE windows lack a constant size as the windows
%      interfaces dynamically adjust window size before return, therefore
%      the current number of intervals equals the maximum number.
%
%-Particulars
%
%   On input, 'window' is a 2Mx1 array. During validation, the intervals
%   are ordered, and overlapping intervals are merged. The size of the output
%   window equals the number of endpoints remaining, and window is ready for
%   use with any of the window routines.
%
%-Exceptions
%
%   1)  If the original number of endpoints `n' is odd, the error
%       SPICE(UNMATCHENDPTS) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the left endpoint is greater than the right endpoint, the
%       error SPICE(BADENDPOINTS) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the input argument `window_i' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `window_i' is not of the expected type,
%       or it does not have the expected dimensions and size, an error
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
%   -Mice Version 1.2.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and modified code example to produce
%       formatted output.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.1.1, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.1.0, 27-JUL-2009 (EDW)
%
%       'zzmice_cell' modified to include the cell type identifier 'double'.
%
%   -Mice Version 1.0.0, 17-DEC-2008 (EDW)
%
%-Index_Entries
%
%   validate a d.p. window
%
%-&

function [window_f] = cspice_wnvald( window_i )

   switch nargin

      case 1

         %
         % In this case, the interface requires only a numeric
         % `window_i' with 2Nx1 dimension.
         %

         window_i = zzmice_cell( window_i, 'double');

      otherwise

         error ( 'Usage: [window_f] = cspice_wnvald( window_i )' )

   end

   %
   % Please note, this call does not require addition of space for the window
   % control segment as needed by other windows interfaces. The interface
   % copies the data in `window_i' to a work variable rather than directly
   % pass `window_i' into a CSPICE call.
   %
   try
      [window_f] = mice('wnvald_c', window_i );
   catch spiceerr
      rethrow(spiceerr)
   end



