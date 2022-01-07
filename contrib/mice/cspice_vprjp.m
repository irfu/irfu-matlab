%-Abstract
%
%   CSPICE_VPRJP projects orthogonally a vector onto a specified plane.
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
%      vin      the vector to orthogonally project onto a specified plane.
%
%               [3,1] = size(vin); double = class(vin)
%
%      plane    a structure describing a SPICE plane onto which to
%               project `vin'.
%
%               [1,1] = size(plane); struct = class(plane)
%
%               The structure has the fields:
%
%                  normal:   [3,1] = size(normal);   double = class(normal)
%                  constant: [1,1] = size(constant); double = class(constant)
%
%               The normal vector component of a SPICE plane has unit length.
%
%   the call:
%
%      [vout] = cspice_vprjp( vin, plane )
%
%   returns:
%
%      vout     the vector resulting from the orthogonal projection of `vin'
%               onto `plane'.
%
%               [3,1] = size(vout); double = class(vout)
%
%               `vout' is the closest point in the specified plane to `vin'.
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
%   1) Find the closest point in the ring plane of a planet to a
%      spacecraft located at a point (in body-fixed coordinates).
%
%
%      Example code begins here.
%
%
%      function vprjp_ex1()
%
%         %
%         % Set the spacecraft location and define the normal
%         % vector as the normal to the equatorial plane, and
%         % the origin at the body/ring center.
%         %
%         scpos = [-29703.16955, 879765.72163, -137280.21757]';
%
%         norm  = [0.0, 0.0, 1.0]';
%
%         orig  = [0.0, 0.0, 0.0]';
%
%         %
%         % Create the plane structure.
%         %
%         [ringpl] = cspice_nvp2pl( norm, orig );
%
%         %
%         % Project the position vector onto the ring plane.
%         %
%         [proj] = cspice_vprjp( scpos, ringpl );
%
%         fprintf( 'Projection of S/C position onto ring plane:\n'     )
%         fprintf( '%17.5f %16.5f %16.5f\n', proj(1), proj(2), proj(3) )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Projection of S/C position onto ring plane:
%           -29703.16955     879765.72163          0.00000
%
%
%-Particulars
%
%   Projecting a vector `vin' orthogonally onto a plane can be thought of
%   as finding the closest vector in the plane to `vin'. This "closest
%   vector" always exists; it may be coincident with the original
%   vector.
%
%   Two related routines are cspice_vprjpi, which inverts an orthogonal
%   projection of a vector onto a plane, and cspice_vproj, which projects
%   a vector orthogonally onto another vector.
%
%-Exceptions
%
%   1)  If the normal vector of the input plane does not have unit
%       length (allowing for round-off error), the error
%       SPICE(NONUNITNORMAL) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If any of the input arguments, `vin' or `plane', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `vin' or `plane', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  It is recommended that the input plane be created by one of
%       the Mice routines
%
%          cspice_nvc2pl ( Normal vector and constant to plane )
%          cspice_nvp2pl ( Normal vector and point to plane    )
%          cspice_psv2pl ( Point and spanning vectors to plane )
%
%       In any case the input plane must have a unit length normal
%       vector and a plane constant consistent with the normal
%       vector.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  G. Thomas and R. Finney, "Calculus and Analytic Geometry,"
%        7th Edition, Addison Wesley, 1988.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 27-AUG-2021 (EDW) (JDR) (NJB)
%
%       Added error check for non-unit plane normal vector.
%
%       Edited the header to comply with NAIF standard.
%       Modified code example to produce formatted output and use actual
%       data.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 18-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 27-AUG-2012 (EDW)
%
%-Index_Entries
%
%   vector projection onto plane
%
%-&

function [vout] = cspice_vprjp( vin, plane )

   switch nargin

      case 2

         vin   = zzmice_dp( vin );
         plane = zzmice_pln( plane );

      otherwise

         error ( 'Usage: [vout(3)] = cspice_vprjp( vin(3), plane )' )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [vout] = mice('vprjp_c', vin, plane );
   catch spiceerr
      rethrow(spiceerr)
   end
