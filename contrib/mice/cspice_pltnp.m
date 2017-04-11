%-Abstract
%
%   CSPICE_PLTNP finds the nearest point on a triangular plate to a
%   given point.
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
%      point      is an arbitrary point in 3-dimensional space.
%
%                 [3,1] = size(point); double = class(point)
%
%      v1,
%      v2,
%      v3         are 3-vectors constituting the vertices of
%                 a triangular plate.
%
%                 The plate is allowed to be degenerate: it may
%                 consist of a line segment or of a single point.
%
%                 [3,1] = size(v1); double = class(v1)
%                 [3,1] = size(v2); double = class(v2)
%                 [3,1] = size(v3); double = class(v3)
%
%   the call:
%
%       [pnear, dist] = cspice_pltnp(point, v1, v2, v3)
%
%   returns:
%
%      pnear      is the closest point on the plate to `point'.
%                 `pnear' is unique, since the plate is convex.
%
%                 [3,1] = size(pnear); double = class(pnear)
%
%      dist       is the distance between `point' and `pnear'.
%
%                 [1,1] = size(dist); double = class(dist)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Find the nearest point to the point (2,2,2) on a plate having
%      vertices at the unit basis vectors that lie along the positive
%      X, Y, and Z coordinate axes.
%
%      function pltnp_t
%
%         point = [2.0, 2.0, 2.0]';
%         v1    = [1.0, 0.0, 0.0]';
%         v2    = [0.0, 1.0, 0.0]';
%         v3    = [0.0, 0.0, 1.0]';
%
%         [pnear, dist] = cspice_pltnp(point, v1, v2, v3);
%
%
%         fprintf ( [ '\n' ...
%                     'Plate vertex 1 = %14.7e %14.7e %14.7e\n' ...
%                     'Plate vertex 2 = %14.7e %14.7e %14.7e\n' ...
%                     'Plate vertex 3 = %14.7e %14.7e %14.7e\n' ...
%                     'Input point    = %14.7e %14.7e %14.7e\n' ...
%                     '\n'                                      ...
%                     'Near point     = %14.7e %14.7e %14.7e\n' ...
%                     'Distance       = %14.7e\n'               ...
%                     '\n'],                                    ...
%                     v1(1),    v1(2),    v1(3),                ...
%                     v2(1),    v2(2),    v2(3),                ...
%                     v3(1),    v3(2),    v3(3),                ...
%                     point(1), point(2), point(3),             ...
%                     pnear(1), pnear(2), pnear(3),             ...
%                     dist                                    )
%
%   Matlab outputs:
%
%      Plate vertex 1 =  1.0000000e+00  0.0000000e+00  0.0000000e+00
%      Plate vertex 2 =  0.0000000e+00  1.0000000e+00  0.0000000e+00
%      Plate vertex 3 =  0.0000000e+00  0.0000000e+00  1.0000000e+00
%      Input point    =  2.0000000e+00  2.0000000e+00  2.0000000e+00
%
%      Near point     =  3.3333333e-01  3.3333333e-01  3.3333333e-01
%      Distance       =  2.8867513e+00
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pltnp_c.
%
%   MICE.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 16-MAR-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   nearest point on triangular plate
%
%-&

function [pnear, dist] = cspice_pltnp(point, v1, v2, v3)

   switch nargin
      case 4

         point = zzmice_dp(point);
         v1    = zzmice_dp(v1);
         v2    = zzmice_dp(v2);
         v3    = zzmice_dp(v3);

      otherwise

         error ( 'Usage: [pnear, dist] = cspice_pltnp(point, v1, v2, v3)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [pnear, dist] = mice( 'pltnp_c', point, v1, v2, v3);
   catch
      rethrow(lasterror)
   end



