%-Abstract
%
%   CSPICE_RECLAT converts rectangular (Cartesian) coordinates to
%   latitudinal coordinates. All coordinates are expressed as
%   double precision values.
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
%      rectan   the array(s) containing the rectangular coordinates of the
%               position or set of positions.
%
%               [3,n] = size(rectan); double = class(rectan)
%
%   the call:
%
%      [radius, lon, lat] = cspice_reclat(rectan)
%
%   returns:
%
%      radius   the value(s) describing the distance of the position
%               from the origin.
%
%               [1,n] = size(radius); double = class(radius)
%
%      lon      the value(s) describing the angle of the position from
%               the XZ plane measured in radians.
%
%               [1,n] = size(lon); double = class(lon)
%
%      lat      the value(s) describing the angle of the position from the
%               XY plane measured in radians.
%
%               [1,n] = size(lat); double = class(lat)
%
%               The argument `radius' returns in the same units associated
%               with `rectan'.
%
%               `radius', `lon', and `lat' return with
%               the same vectorization measure, N, as `rectan'.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Compute the latitudinal coordinates of the position of the Moon
%      as seen from the Earth, and convert them to rectangular
%      coordinates.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: reclat_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            naif0012.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'naif0012.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function reclat_ex1()
%
%         %
%         % Load an SPK and leapseconds kernels.
%         %
%         cspice_furnsh( 'reclat_ex1.tm' )
%
%         %
%         % Convert the time to ET.
%         %
%         et = cspice_str2et( '2017 Mar 20' );
%
%         %
%         % Retrieve the position of the moon seen from earth at `et'
%         % in the J2000 frame without aberration correction.
%         %
%         [pos, et] = cspice_spkpos( 'MOON', et, 'J2000', 'NONE', 'EARTH' );
%
%         fprintf( 'Original rectangular coordinates:\n' )
%         fprintf( '   X          (km): %20.8f\n', pos(1) )
%         fprintf( '   Y          (km): %20.8f\n', pos(2) )
%         fprintf( '   Z          (km): %20.8f\n', pos(3) )
%
%         %
%         % Convert the position vector `pos' to latitudinal
%         % coordinates.
%         %
%         [radius, lon, lat] = cspice_reclat(pos);
%         fprintf( '\n' )
%         fprintf( 'Latitudinal coordinates:\n' )
%         fprintf( '   Radius     (km): %20.8f\n', radius )
%         fprintf( '   Longitude (deg): %20.8f\n', lon * cspice_dpr )
%         fprintf( '   Latitude  (deg): %20.8f\n', lat * cspice_dpr )
%
%         %
%         % Convert the latitudinal to rectangular.
%         %
%         [rectan] = cspice_latrec( radius, lon, lat);
%         fprintf( '\n' )
%         fprintf( 'Rectangular coordinates from cspice_latrec:\n' )
%         fprintf( '   X          (km): %20.8f\n', rectan(1) )
%         fprintf( '   Y          (km): %20.8f\n', rectan(2) )
%         fprintf( '   Z          (km): %20.8f\n', rectan(3) )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Original rectangular coordinates:
%         X          (km):      -55658.44323296
%         Y          (km):     -379226.32931475
%         Z          (km):     -126505.93063865
%
%      Latitudinal coordinates:
%         Radius     (km):      403626.33912495
%         Longitude (deg):         -98.34959789
%         Latitude  (deg):         -18.26566077
%
%      Rectangular coordinates from cspice_latrec:
%         X          (km):      -55658.44323296
%         Y          (km):     -379226.32931475
%         Z          (km):     -126505.93063865
%
%
%   2) Create a table showing a variety of rectangular coordinates
%      and the corresponding latitudinal coordinates.
%
%      Corresponding rectangular and latitudinal coordinates are
%      listed to three decimal places. Output angles are in degrees.
%
%
%      Example code begins here.
%
%
%      function reclat_ex2()
%
%         %
%         % Define eleven sets of rectangular coordinates.
%         %
%         rec = [ [ 0., 1., 0., 0., -1., 0., 0., 1., 1., 0., 1. ];         ...
%                 [ 0., 0., 1., 0., 0., -1., 0., 1., 0., 1., 1. ];         ...
%                 [ 0., 0., 0., 1., 0., 0., -1., 0., 1., 1., 1. ]    ];
%
%         %
%         % ...convert the rectangular coordinates to latitudinal coordinates
%         %
%         [radius, lon, lat] = cspice_reclat(rec);
%
%         %
%         % Convert `lon' and `lat' to degrees.
%         %
%         lon = lon * cspice_dpr;
%         lat  = lat  * cspice_dpr;
%
%         %
%         % Create an array of values for output.
%         %
%         output = [ radius; lon; lat; rec ];
%
%         %
%         % Output banner.
%         %
%         disp('    r       lon      lat    rect(1)  rect(2)  rect(3)')
%         disp(' -------  -------  -------  -------  -------  -------')
%
%         txt = sprintf( '%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',  output );
%         disp( txt )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%          r       lon      lat    rect(1)  rect(2)  rect(3)
%       -------  -------  -------  -------  -------  -------
%         0.000    0.000    0.000    0.000    0.000    0.000
%         1.000    0.000    0.000    1.000    0.000    0.000
%         1.000   90.000    0.000    0.000    1.000    0.000
%         1.000    0.000   90.000    0.000    0.000    1.000
%         1.000  180.000    0.000   -1.000    0.000    0.000
%         1.000  -90.000    0.000    0.000   -1.000    0.000
%         1.000    0.000  -90.000    0.000    0.000   -1.000
%         1.414   45.000    0.000    1.000    1.000    0.000
%         1.414    0.000   45.000    1.000    0.000    1.000
%         1.414   90.000   45.000    0.000    1.000    1.000
%         1.732   45.000   35.264    1.000    1.000    1.000
%
%
%-Particulars
%
%   This routine returns the latitudinal coordinates of a point
%   whose position is input in rectangular coordinates.
%
%   Latitudinal coordinates are defined by a distance from a central
%   reference point, an angle from a reference meridian, and an angle
%   above the equator of a sphere centered at the central reference
%   point.
%
%-Exceptions
%
%   1)  If the X and Y components of `rectan' are both zero, the
%       longitude is set to zero.
%
%   2)  If `rectan' is the zero vector, longitude and latitude are
%       both set to zero.
%
%   3)  If the input argument `rectan' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `rectan' is not of the expected type, or
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
%       Changed output arguments "longitude" and "latitude" to
%       "lon" and "lat" for consistency with other routines.
%
%       Edited the header to comply with NAIF standard. Added
%       meta-kernel to example #1. Updated code example #1 to produce
%       formatted output and added a call to cspice_kclear. Added the
%       problem statement to both examples.
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
%   -Mice Version 1.0.1, 01-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   rectangular to latitudinal coordinates
%
%-&

function [radius, lon, lat] = cspice_reclat(rectan)

   switch nargin
      case 1

         rectan = zzmice_dp(rectan);

      otherwise
         error ( ['Usage: [_radius_, _lon_, _lat_] = ' ...
                  'cspice_reclat(_rectan(3)_)'] )
   end

   %
   % Call the MEX library.
   %
   try
      [radius, lon, lat] = mice('reclat_c',rectan);
   catch spiceerr
      rethrow(spiceerr)
   end

