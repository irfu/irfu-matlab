%-Abstract
%
%   CSPICE_WNCARD returns the cardinality (number of intervals) of a
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
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%   the call:
%
%      [wncard] = cspice_wncard( window )
%
%   returns:
%
%      wncard   the cardinality (number of intervals) of the input `window'.
%
%               [1,1] = size(wncard); int32 = class(wncard)
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
%   1) Given a double precision window, determine its cardinality
%      (number of intervals).
%
%      Example code begins here.
%
%
%      function wncard_ex1()
%
%         %
%         % Define a window with three intervals (six values).
%         %
%         window = [ [ 1.; 3.]; [ 7.; 11.]; [23.; 27.] ];
%
%         %
%         % What's the window cardinality of `window'?
%         %
%         card = cspice_wncard(window);
%
%         %
%         % Print the window cardinality (this ought to show "3" ).
%         %
%         fprintf('Number of intervals in the window: %d\n', card )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Number of intervals in the window: 3
%
%
%-Particulars
%
%   This function returns the value numel(window)/2.
%
%-Exceptions
%
%   1)  If the number of elements in `window' is not even, the error
%       SPICE(INVALIDSIZE) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the input argument `window' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `window' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
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
%       example's problem statement.
%
%       Changed output argument name "card" to "wncard" to comply with
%       NAIF standard.
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
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 30-DEC-2008 (EDW)
%
%-Index_Entries
%
%   cardinality of a d.p. window
%
%-&

function [wncard] = cspice_wncard(window)

   switch nargin

      case 1

         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [wncard] = cspice_wncard( window )' )

   end

   %
   % Call the windows routine, add to `window' the space needed for
   % the control segments.
   %
   try
      [wncard] = mice('wncard_c', [zeros(6,1); window] );
   catch spiceerr
      rethrow(spiceerr)
   end


