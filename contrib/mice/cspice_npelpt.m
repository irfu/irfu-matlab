%-Abstract
%
%   CSPICE_NPELPT finds the nearest point on an ellipse to a specified point,
%   both in three-dimensional space, and finds the distance between the
%   ellipse and the point.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      point    a point in 3-dimensional space.
%
%               [3,1] = size(point); double = class(point)
%
%      ellips   a SPICE ellipse that represents an ellipse in
%               three-dimensional space.
%
%               [1,1] = size(ellips); struct = class(ellips)
%
%               The structure has the fields:
%
%                 center:    [3,1] = size(center); double = class(center)
%                 semiMinor: [3,1] = size(semiMinor); double =
%                 class(semiMinor) semiMajor: [3,1] = size(semiMajor); double
%                 = class(semiMajor)
%
%   the call:
%
%      [pnear, dist] = cspice_npelpt( point, ellips )
%
%   returns:
%
%      pnear    the nearest point on `ellips' to `point'.
%
%               [3,1] = size(pnear); double = class(pnear)
%
%      dist     the distance between `point' and `pnear'.
%
%               [1,1] = size(dist); double = class(dist)
%
%               This is the distance between `point' and the ellipse.
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
%   1) Create a SPICE ellipse given its center and its semi-major and
%      semi-minor axes, and calculate the location on that ellipse
%      which is closest to a defined point, and the distance between
%      those points.
%
%      Example code begins here.
%
%
%      function npelpt_ex1()
%
%         %
%         % Define a center, and semimajor and semiminor axes for
%         % an ellipse.
%         %
%         % Also define an arbitrary point in space.
%         %
%         center = [  1.;  2.; -3. ];
%         smajor = [  3.;  0.;  0. ];
%         sminor = [  0.;  2.;  0. ];
%         point  = [ -4.;  2.;  1. ];
%
%         %
%         % Create an ellipse structure using `center', `smajor',
%         % and `sminor'.
%         %
%         ellips = cspice_cgv2el( center, smajor, sminor );
%
%         fprintf( 'Input SPICE ellipse:\n' );
%         fprintf( '  Semi-minor axis: %10.6f %10.6f %10.6f\n',   ...
%                                              ellips.semiMinor);
%         fprintf( '  Semi-major axis: %10.6f %10.6f %10.6f\n',   ...
%                                              ellips.semiMajor);
%         fprintf( '  Center         : %10.6f %10.6f %10.6f\n\n',   ...
%                                              ellips.center   );
%
%         %
%         % Calculate the location on the ellipse closest to
%         % the defined point.
%         %
%         [pnear, dist] = cspice_npelpt( point, ellips );
%         fprintf( 'Nearest point    : %10.6f %10.6f %10.6f\n', pnear )
%         fprintf( 'Distance         : %10.6f\n', dist                )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Input SPICE ellipse:
%        Semi-minor axis:   0.000000   2.000000   0.000000
%        Semi-major axis:   3.000000   0.000000   0.000000
%        Center         :   1.000000   2.000000  -3.000000
%
%      Nearest point    :  -2.000000   2.000000  -3.000000
%      Distance         :   4.472136
%
%
%-Particulars
%
%   Given an ellipse and a point in 3-dimensional space, if the
%   orthogonal projection of the point onto the plane of the ellipse
%   is on or outside of the ellipse, then there is a unique point on
%   the ellipse closest to the original point. This routine finds
%   that nearest point on the ellipse. If the projection falls inside
%   the ellipse, there may be multiple points on the ellipse that are
%   at the minimum distance from the original point. In this case,
%   one such closest point will be returned.
%
%   This routine returns a distance, rather than an altitude, in
%   contrast to the Mice routine cspice_nearpt. Because our ellipse is
%   situated in 3-space and not 2-space, the input point is not
%   "inside" or "outside" the ellipse, so the notion of altitude does
%   not apply to the problem solved by this routine. In the case of
%   cspice_nearpt, the input point is on, inside, or outside the ellipsoid,
%   so it makes sense to speak of its altitude.
%
%-Exceptions
%
%   1)  If the input ellipse `ellips' has one or both semi-axis lengths
%       equal to zero, the error SPICE(DEGENERATECASE) is signaled by
%       a routine in the call tree of this routine.
%
%   2)  If the geometric ellipse represented by `ellips' does not
%       have a unique point nearest to the input point, any point
%       at which the minimum distance is attained may be returned
%       in `pnear'.
%
%   3)  If a ratio of non-zero ellipse radii violates the constraints
%       imposed by cspice_nearpt, an error is signaled by a routine in the
%       call tree of this routine.
%
%   4)  The routine does not check for overflow when scaling or
%       translating the input point.
%
%   5)  If any of the input arguments, `point' or `ellips', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `point' or `ellips', is not of
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
%   ELLIPSES.REQ
%   MICE.REQ
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
%   -Mice Version 1.1.0, 27-AUG-2021 (EDW) (JDR)
%
%       Changed input argument name "ellipse" to "ellips".
%
%       Edited the header to comply with NAIF standard. Added
%       code example to -Examples section.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 14-NOV-2014 (EDW)
%
%-Index_Entries
%
%   nearest point on ellipse to point
%
%-&

function [ pnear, dist ] = cspice_npelpt( point, ellips )

   switch nargin
      case 2

         point  = zzmice_dp(point);
         ellips = zzmice_ell(ellips);

      otherwise

         error ( 'Usage: [pnear, dist] = cspice_npelpt( point, ellips )' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [npelpt] = mice( 'npelpt_s', point, ellips );
      pnear    = reshape( [npelpt.pos], 3, [] );
      dist     = reshape( [npelpt.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end
