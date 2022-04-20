%-Abstract
%
%   CSPICE_NPLNPT calculates the location on a defined line
%   nearest to a specified point, then determines the distance
%   between the two points.
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
%      linpt,
%      lindir   are, respectively, a point and a direction vector that define
%               a line.
%
%               [3,1] = size(linpt); double = class(linpt)
%               [3,1] = size(lindir); double = class(lindir)
%
%               The line is the set of vectors
%
%                  linept   +   t * linedr
%
%               where `t' is any real number.
%
%      point    a point in 3-dimensional space.
%
%               [3,n] = size(point); double = class(point)
%
%   the call:
%
%      [pnear, dist] = cspice_nplnpt( linpt, lindir, point )
%
%   returns:
%
%      pnear    the nearest point on the input line to the input `point'.
%
%               [3,1] = size(pnear); double = class(pnear)
%
%      dist     distance between the input line and input point.
%
%               [1,n] = size(dist); double = class(dist)
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
%   1) Define a line, given a point and a direction for the line, and an
%      arbitrary point in space, and calculate the location on the line
%      nearest to the arbitray point, and the distance between these two
%      points.
%
%      Example code begins here.
%
%
%      function nplnpt_ex1()
%
%         %
%         % Define a point on a line, a direction for the line, and
%         % an arbitrary point in space.
%         %
%         linept = [  1, 2,  3 ]';
%         linedr = [  0, 1,  1 ]';
%         point  = [ -6, 9, 10 ]';
%
%         %
%         % Calculate the location on the line nearest the point
%         % and the distance between the location and the defined
%         % point.
%         %
%         [pnear, dist] = cspice_nplnpt( linept, linedr, point  );
%         fprintf('Nearest point: %15.8f %15.8f %15.8f\n', pnear)
%         fprintf('Distance     : %15.8f\n', dist               )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Nearest point:      1.00000000      9.00000000     10.00000000
%      Distance     :      7.00000000
%
%
%-Particulars
%
%   For every line L and point P, there is a unique closest point
%   on L to P. Call this closest point C. It is always true that
%   P - C  is perpendicular to L, and the length of P - C is called
%   the "distance" between P and L.
%
%-Exceptions
%
%   1)  If the line direction vector `lindir' is the zero vector, the
%       error SPICE(ZEROVECTOR) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If any of the input arguments, `linpt', `lindir' or `point',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   3)  If any of the input arguments, `linpt', `lindir' or `point',
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
%   -Mice Version 1.1.0, 07-AUG-2020 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement, and updated code example to produce
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
%   -Mice Version 1.0.0, 14-NOV-2013 (EDW) (SCK)
%
%-Index_Entries
%
%   distance between point and line
%   nearest point on line to point
%
%-&

function [pnear, dist] = cspice_nplnpt( linpt, lindir, point )

   switch nargin
      case 3

         linpt  = zzmice_dp(linpt);
         lindir = zzmice_dp(lindir);
         point  = zzmice_dp(point);

      otherwise

         error ( ['Usage: [pnear(3), dist] = ' ...
                  'cspice_nplnpt( linpt(3), lindir(3), point(3) )'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [nplnpt] = mice( 'nplnpt_s', linpt, lindir, point );
      pnear    = reshape( [nplnpt.pos], 3, [] );
      dist     = reshape( [nplnpt.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end


