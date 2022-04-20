%-Abstract
%
%   CSPICE_INRYPL finds the intersection of a ray and a plane.
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
%      vertex   the position of a ray vertex.
%
%               [3,1] = size(normal); double = class(normal)
%
%      dir      the direction of a ray from `vertex'.
%
%               [3,1] = size(normal); double = class(normal)
%
%      plane    a structure describing a SPICE plane.
%
%               [1,1] = size(plane); struct = class(plane)
%
%               The structure has the fields:
%
%                  normal:   [3,1] = size(normal); double = class(normal)
%                  constant: [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [nxpts, xpt] = cspice_inrypl( vertex, dir, plane )
%
%   returns:
%
%      nxpts   the number of points of intersection of the
%              input ray and plane. Values and meanings of nxpts are:
%
%                 0     No intersection.
%
%                 1     One point of intersection. Note that
%                       this case may occur when the ray's
%                       vertex is in the plane.
%
%                 -1    An infinite number of points of
%                       intersection; the ray lies in the plane.
%
%              [1,1] = size(nxpts); int32 = class(nxpts)
%
%      xpt     the point of intersection of the input ray and plane, when
%              there is exactly one point of intersection.
%
%              If the ray lies in the plane, `xpt' is set equal to
%              vertex.
%
%              If there is no intersection, `xpt' is the zero vector.
%
%              [3,1] = size(xpt); double = class(xpt)
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
%   1) Determine the intersection between the Saturn ring plane and
%      a look direction as seen from a position in the Saturn
%      body-fixed frame. For this extremely simplistic example,
%      we take the equatorial plane as the ring plane.
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model and orientation data for Saturn.
%
%         pck00008.tpc
%
%
%      Example code begins here.
%
%
%      function inrypl_ex1()
%
%         %
%         % Load the PCK.
%         %
%         cspice_furnsh( 'pck00008.tpc' )
%
%         %
%         % Retrieve the triaxial radii of Saturn (699)
%         %
%         radii = cspice_bodvrd( 'SATURN', 'RADII', 3 );
%
%         %
%         % Define a position in the IAU_SATURN frame at three equatorial
%         % radius out along the x axis, a half radius above the
%         % equatorial plane. For this example, we'll assume 'vertex'
%         % represents the light-time corrected position of a vehicle
%         % to the Saturn ring plane.
%         %
%         vertex = [ 3.0 * radii(1), 0.0, radii(3) *.50 ]';
%
%         %
%         % Define a look vector in the y-z plane from 'vertex'.
%         %
%         %   'vertex'
%         %      *______ y
%         %     /|\
%         %    / | \  30 degrees
%         %   /  |  \
%         %  x  -z  'dir'
%         %
%         dir = [ 0.,
%                 cos( 30. *cspice_rpd() ),
%                -sin( 30. *cspice_rpd() )
%               ];
%
%         %
%         % Define the equatorial plane as a SPICE plane. The Z
%         % axis is normal to the plane, the origin lies in the
%         % plane.
%         %
%         normal = [ 0., 0., 1.]';
%         point  = [ 0., 0., 0.]';
%         plane  = cspice_nvp2pl( normal, point );
%
%         %
%         % Determine the intersection point of 'dir' and 'plane', if
%         % such an intersection exists.
%         %
%         [ nxpts, xpt ] =cspice_inrypl( vertex, dir, plane );
%
%
%         %
%         % Do we have an intersection?
%         %
%         if ( nxpts == 1 )
%            fprintf(['Vector intersects plane at: ',    ...
%                     '%12.9g %12.9g %12.9g\n'], xpt(:) )
%         end
%
%         %
%         % No intersection
%         %
%         if ( nxpts == 0 )
%            fprintf( 'No intersection between vector and plane.\n' )
%         end
%
%         %
%         % No intersection
%         %
%         if ( nxpts == -1 )
%            fprintf( 'Vector lies in plane, degenerate case.\n' )
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Vector intersects plane at:       180804   47080.6051            0
%
%
%-Particulars
%
%   The intersection of a ray and plane in three-dimensional space
%   can be a the empty set, a single point, or the ray itself.
%
%-Exceptions
%
%   1)  If the ray's direction vector is the zero vector, the error
%       SPICE(ZEROVECTOR) is signaled by a routine in the call tree of
%       this routine. `nxpts' and `xpt' are not modified.
%
%   2)  If the ray's vertex is further than cspice_dpmax() / 3 from the
%       origin, the error SPICE(VECTORTOOBIG) is signaled by a routine
%       in the call tree of this routine. `nxpts' and `xpt' are not
%       modified.
%
%   3)  If the input plane is further than cspice_dpmax() / 3 from the
%       origin, the error SPICE(VECTORTOOBIG) is signaled by a routine
%       in the call tree of this routine. `nxpts' and `xpt' are not
%       modified.
%
%   4)  The input plane should be created by one of the Mice
%       routines
%
%          cspice_nvc2pl
%          cspice_nvp2pl
%          cspice_psv2pl
%
%       Invalid input planes will cause unpredictable results.
%
%   5)  In the interest of good numerical behavior, in the case
%       where the ray's vertex is not in the plane, this routine
%       considers that an intersection of the ray and plane occurs
%       only if the distance between the ray's vertex and the
%       intersection point is less than cspice_dpmax() / 3.
%
%       If `vertex' is not in the plane and this condition is not
%       met, then `nxpts' is set to 0 and `xpt' is set to the zero
%       vector.
%
%   6)  If any of the input arguments, `vertex', `dir' or `plane', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `vertex', `dir' or `plane', is
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
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   PLANES.REQ
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
%       Edited -Examples section to comply with NAIF standard. Added
%       -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 27-AUG-2012 (EDW)
%
%-Index_Entries
%
%   intersection of ray and plane
%
%-&

function [nxpts, xpt] = cspice_inrypl( vertex, dir, plane )

   switch nargin
      case 3

         vertex  = zzmice_dp(vertex);
         dir     = zzmice_dp(dir);
         plane   = zzmice_pln(plane);

      otherwise

         error ( 'Usage: [nxpts, xpt] = cspice_inrypl( vertex, dir, plane )' )

   end

   try
      [nxpts, xpt] = mice( 'inrypl_c', vertex, dir, plane );
   catch spiceerr
      rethrow(spiceerr)
   end


