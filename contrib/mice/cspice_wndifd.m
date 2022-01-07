%-Abstract
%
%   CSPICE_WNDIFD returns the difference of two double precision
%   windows.
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
%      a   SPICE window
%
%          [2l,1] = size(a); double = class(a)
%
%      b   SPICE window
%
%          [2m,1] = size(b); double = class(b)
%
%          Two SPICE windows containing zero or more intervals.
%
%   the call:
%
%      c = cspice_wndifd( a, b )
%
%   returns:
%
%      c   the window difference (in the SPICE sense) of 'a' and 'b', every
%          point contained in 'a', but not contained in 'b'.
%
%          [2n,1] = size(c); double = class(c)
%
%          'c' can overwrite 'a' or 'b'.
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
%   1) Let `a' contain the intervals
%
%         [ 1, 3 ]  [ 7, 11 ]  [ 23, 27 ]
%
%      and `b' contain the intervals
%
%         [ 2, 4 ]  [ 8, 10 ]  [ 16, 18 ]
%
%      Then the difference of `a' and `b' contains the intervals
%
%         [ 1, 2 ]  [ 7, 8 ]  [ 10, 11 ]  [ 23, 27 ]
%
%      The following code example demonstrates how to compute this
%      difference of `a' and `b' using SPICE.
%
%
%      Example code begins here.
%
%
%      function wndifd_ex1()
%
%         %
%         % Let 'a' contain the intervals
%         %
%         a = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; ];
%
%         %
%         % and b contain the intervals
%         %
%         b = [ [ 2; 4 ]; [ 8; 10 ]; [ 16; 18 ]; ];
%
%         %
%         % Then the difference of`'a' and `b', `c':
%         %
%         c = cspice_wndifd(a, b);
%
%         %
%         % Output the difference `c'
%         %
%         for i=1:cspice_wncard(c)
%
%            [left, right] = cspice_wnfetd( c, i );
%            fprintf( '%16.6f %16.6f\n', left, right  );
%
%      end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%              1.000000         2.000000
%              7.000000         8.000000
%             10.000000        11.000000
%             23.000000        27.000000
%
%
%-Particulars
%
%   Mathematically, the difference of two windows contains every
%   point contained in the first window but not contained in the
%   second window.
%
%   Matlab offers no satisfactory floating point representation
%   of open intervals. Thus, for floating point windows we must
%   return the closure of the set theoretical difference: that is,
%   the difference plus the endpoints of the first window that are
%   contained in the second window.
%
%-Exceptions
%
%   1)  The cardinality of the input windows must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   2)  If any of the input arguments, `a' or `b', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `a' or `b', is not of the
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
%       Improved -Particulars section.
%
%   -Mice Version 1.0.1, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 23-JUL-2007 (EDW)
%
%-Index_Entries
%
%   difference two d.p. windows
%
%-&

function [c] = cspice_wndifd( a, b )

   switch nargin

      case 2

         a    = zzmice_win(a);
         b    = zzmice_win(b);

      otherwise

         error ( 'Usage: [c] = cspice_wndifd( a, b )' )

   end

%
% Call the windows routine, add to 'a' and 'b' the space needed for
% the control segments.
%
   try
      [c] = mice('wndifd_c', [zeros(6,1); a], [zeros(6,1); b] );
   catch spiceerr
      rethrow(spiceerr)
   end



