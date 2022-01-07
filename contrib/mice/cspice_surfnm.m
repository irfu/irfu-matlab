%-Abstract
%
%   CSPICE_SURFNM computes the outward-pointing, unit normal vector at a
%   point on the surface of an ellipsoid.
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
%      a        the length of the semi-axis of the ellipsoid that is parallel
%               to the X-axis of the body-fixed reference frame.
%
%               [1,1] = size(a); double = class(a)
%
%      b        the length of the semi-axis of the ellipsoid that is parallel
%               to the Y-axis of the body-fixed reference frame.
%
%               [1,1] = size(b); double = class(b)
%
%      c        the length of the semi-axis of the ellipsoid that is parallel
%               to the Z-axis of the body-fixed reference frame.
%
%               [1,1] = size(c); double = class(c)
%
%      point    3-vector(s) giving the bodyfixed coordinates of a point on the
%               ellipsoid.
%
%               [3,n] = size(point); double = class(point)
%
%               In bodyfixed coordinates, the semi-axes of the ellipsoid
%               are aligned with the X, Y, and Z-axes of the reference frame.
%
%   the call:
%
%      [normal] = cspice_surfnm( a, b, c, point )
%
%   returns:
%
%      normal   the unit vector(s) pointing away from the ellipsoid and normal
%               to the ellipsoid at `point'.
%
%               If [3,1] = size(point) then [3,3]   = size(normal)
%               If [3,n] = size(point) then [3,3,n] = size(normal)
%                                            double = class(normal)
%
%               `normal' returns with the same vectorization measure, N,
%               as `point'.
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
%   1) Compute the surface normal to an ellipsoid defined by its three
%      radii at a set of locations.
%
%      Example code begins here.
%
%
%      function surfnm_ex1()
%
%         %
%         % Define the radii of an ellipsoid.
%         %
%         a  =  1.;
%         b  =  2.;
%         c  =  3.;
%
%         %
%         % Select a set of locations, three 3-vectors.
%         %
%         point = [ [ 0.; 0.; 3.], [ 0.; 2.; 0.], [-1; 0; 0] ];
%
%         %
%         % Calculate the surface normal to the ellipsoid at 'point'.
%         %
%         out_norm = cspice_surfnm( a, b, c, point);
%
%         n_elements = size(out_norm,2);
%         for i=1:n_elements
%            fprintf( ['The normal at (%4.1f,%4.1f,%4.1f)'                 ...
%                      ' equals (%4.1f,%4.1f,%4.1f)\n'],                   ...
%                     point(:,i), out_norm(:,i) );
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      The normal at ( 0.0, 0.0, 3.0) equals ( 0.0, 0.0, 1.0)
%      The normal at ( 0.0, 2.0, 0.0) equals ( 0.0, 1.0, 0.0)
%      The normal at (-1.0, 0.0, 0.0) equals (-1.0, 0.0, 0.0)
%
%
%-Particulars
%
%   This routine computes the outward pointing unit normal vector to
%   the ellipsoid having semi-axes of length `a', `b', and `c' from the
%   point `point'.
%
%-Exceptions
%
%   1)  If any of the axes are non-positive, the error
%       SPICE(BADAXISLENGTH) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If any of the input arguments, `a', `b', `c' or `point', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `a', `b', `c' or `point', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  It is assumed that the input `point' is indeed on the ellipsoid.
%       No checking for this is done.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Reformatted example's
%       output and added problem statement.
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
%   -Mice Version 1.0.1, 17-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 15-JUN-2006 (EDW)
%
%-Index_Entries
%
%   surface normal vector on an ellipsoid
%
%-&

function [normal] = cspice_surfnm(a, b, c, point)

   switch nargin
      case 4

         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);
         point = zzmice_dp(point);

      otherwise

         error ( ['Usage: [_normal(3)_] = '             ...
                  'cspice_surfnm( a, b, c, _point(3)_ )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [normal] = mice( 'surfnm_c',  a, b, c, point);
   catch spiceerr
      rethrow(spiceerr)
   end


