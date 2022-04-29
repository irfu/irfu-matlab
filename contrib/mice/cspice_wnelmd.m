%-Abstract
%
%   CSPICE_WNELMD determines whether a point is an element of a double
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
%      point    value which may or may not exist in one of the intervals in
%               window.
%
%               [1,1] = size(point); double = class(point)
%
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%   the call:
%
%      [wnelmd] = cspice_wnelmd( point, window )
%
%   returns:
%
%      A boolean with value true if `point' exists as an element of
%      `window'.
%
%         a(i)  <  point  <  b(i)
%               -         -
%
%      for some interval [ a(i), b(i) ] in `window', false
%      otherwise.
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
%   1) Given a set of double precision numbers, determine whether they
%      are elements of a double precision window.
%
%      Example code begins here.
%
%
%      function wnelmd_ex1()
%
%         %
%         % Let `window' contain the intervals
%         %
%         window = [ [ 1; 3 ];  [ 7; 11 ];  [ 23; 27 ]; ];
%
%         points = [ 0.0, 1.0, 9.0, 13.0, 29.0 ];
%
%         %
%         % Loop over the points.
%         %
%         for i=1:numel(points)
%            if ( cspice_wnelmd( points(i), window ) )
%               fprintf( 'Point %8.3f - an element of the window\n', ...
%                         points(i) )
%            else
%               fprintf('Point %8.3f - not an element of the window\n', ...
%                         points(i) )
%            end
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Point    0.000 - not an element of the window
%      Point    1.000 - an element of the window
%      Point    9.000 - an element of the window
%      Point   13.000 - not an element of the window
%      Point   29.000 - not an element of the window
%
%
%-Particulars
%
%   None.
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
%   2)  If any of the input arguments, `point' or `window', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `point' or `window', is not of
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
%   -Mice Version 1.1.0, 27-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and modified example code to produce
%       formatted output.
%
%       Added square brackets to output argument in function declaration,
%       and renamed it to "wnelmd".
%
%       Corrected error message format.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
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
%   -Mice Version 1.0.0, 17-JUL-2007 (EDW)
%
%-Index_Entries
%
%   element of a d.p. window
%
%-&

function [wnelmd] = cspice_wnelmd( point, window )

   switch nargin

      case 2

         point  = zzmice_dp(point);
         window = zzmice_win(window);


      otherwise

         error( 'Usage: [wnelmd] = cspice_wnelmd( point, window )' )

      end

   try
      [wnelmd] = mice( 'wnelmd_c', point, [zeros(6,1); window] );
      [wnelmd] = zzmice_logical(wnelmd);
   catch spiceerr
      rethrow(spiceerr)
   end

