%-Abstract
%
%   CSPICE_SRFREC converts planetocentric latitude and longitude
%   of a surface point on a specified body to rectangular
%   coordinates.
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
%      body       the NAIF integer code of an extended body
%                 on which a surface point of interest is located.
%                 The body is modeled as a triaxial ellipsoid.
%
%                 [1,1] = size(body); int32 = class(body)
%
%      longitude  Longitude of the input point.  This is the angle between
%                 the prime meridian and the meridian containing `rectan'.
%                 The direction of increasing longitude is from the +X axis
%                 towards the +Y axis.
%
%                 Longitude is measured in radians.  On input, the range
%                 of longitude is unrestricted.
%
%                 [1,n] = size(lon); double = class(lon)
%
%      latitude   Latitude of the input point. This is the angle from
%                 the XY plane of the ray from the origin through the
%                 point.
%
%                 [1,n] = size(lat); double = class(lat)
%
%   the call:
%
%      rectan = cspice_srfrec( radius, lon, lat)
%
%   returns:
%
%      rectan   the array(s) containing the rectangular coordinates of the
%               position or set of positions
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               'rectan' returns with the same units associated with 'lon'.
%
%               'rectan' returns with the vectorization measure, N, as
%               'lon', and 'lat'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example (1):
%
%      %
%      % NAIF ID for our body of interest.
%      %
%      EARTH =  399;
%
%      %
%      % Load the kernel pool with a PCK file that contains
%      % values for the radii of the Earth.
%      %
%      cspice_furnsh( '/kernels/standard.tm' )
%
%      %
%      % Find 'x', the rectangular coordinates of the surface point
%      % defined by `lat' and `long'.  The NAIF integer code for
%      % the Earth is 399. (See the NAIF_IDS required reading file
%      % for the complete set of codes.)
%      %
%      lon =  100.;
%      lat =   35.;
%
%      fprintf( 'Original latitudinal coordinates: \n' )
%      fprintf( '                 Longitude (deg): %f\n', lon )
%      fprintf( '                 Latitude  (deg): %f\n\n', lat )
%
%
%      %
%      % Convert angles to radians forr input to cspice_srfrec.
%      %
%      x = cspice_srfrec( EARTH, lon*cspice_rpd(), lat*cspice_rpd() );
%
%      fprintf( 'Rectangular coordinates: \n')
%      fprintf( '                 X (km): %f\n', x(1) )
%      fprintf( '                 Y (km): %f\n', x(2) )
%      fprintf( '                 Z (km): %f\n\n', x(3) )
%
%
%      %
%      %
%      % Now try to recover the original latitudinal coordinates
%      % from the rectangular coordinates found by cspice_srfrec.
%      %
%      [radius, lon1, lat1] = cspice_reclat( x);
%
%      %
%      % Convert angles back to degree for display.
%      %
%      fprintf( 'Latitudinal coordinates recovered from \n' )
%      fprintf( 'rectangular coordinates: \n' )
%      fprintf( '                 Longitude (deg): %f\n', lon1*cspice_dpr() )
%      fprintf( '                 Latitude  (deg): %f\n', lat1*cspice_dpr() )
%
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear()
%
%
%   MATLAB outputs:
%
%      Original latitudinal coordinates:
%                       Longitude (deg): 100.000000
%                       Latitude  (deg): 35.000000
%
%      Rectangular coordinates:
%                       X (km): -906.249195
%                       Y (km): 5139.594582
%                       Z (km): 3654.299896
%
%      Latitudinal coordinates recovered from
%      rectangular coordinates:
%                       Longitude (deg): 100.000000
%                       Latitude  (deg): 35.000000
%
%
%   Example (2):
%
%      %
%      % Define ten sets of latitudinal coordinates.
%      %
%      longitudes = [ 0., 90., 0. 180., -90., ...
%                                     0., 45., 0., 90., 45. ];
%      latitudes  = [ 0., 0., 90., 0., 0.,    ...
%                                     -90., 0., 45., 45., 35.2643 ];
%
%      %
%      % Convert angles to radians forr input to cspice_srfrec.
%      %
%      rectan = cspice_srfrec( EARTH, longitudes*cspice_rpd(), ...
%                                     latitudes*cspice_rpd() );
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ longitudes; latitudes; rectan ];
%
%      %
%      % Output banner.
%      %
%      disp('  longitude  latitude       x         y           z   ')
%      disp('  --------   --------   --------   --------   --------')
%
%      txt = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f\n', output );
%      disp( txt )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear()
%
%
%   MATLAB outputs:
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear()
%
%        longitude  latitude       x         y           z   
%        --------   --------   --------   --------   --------
%          0.0000     0.0000  6378.1366     0.0000     0.0000
%         90.0000     0.0000     0.0000  6378.1366     0.0000
%          0.0000    90.0000     0.0000     0.0000  6356.7519
%        180.0000     0.0000 -6378.1366     0.0000     0.0000
%        -90.0000     0.0000     0.0000 -6378.1366     0.0000
%          0.0000   -90.0000     0.0000     0.0000 -6356.7519
%         45.0000     0.0000  4510.0236  4510.0236     0.0000
%          0.0000    45.0000  4502.4440     0.0000  4502.4440
%         90.0000    45.0000     0.0000  4502.4440  4502.4440
%         45.0000    35.2643  3678.2937  3678.2937  3678.2814
%         
%-Particulars
%
%   This routine returns the rectangular coordinates of a surface
%   point on an extended body with known radii, where the location
%   of the surface point is specified in planetocentric latitudinal
%   coordinates.
%
%   Latitudinal coordinates are defined by a distance from a central
%   reference point, an angle from a reference meridian, and an angle
%   above the equator of a sphere centered at the central reference
%   point.  In this case, the distance from the central reference
%   point is not required as an input because the fact that the
%   point is on the body's surface allows one to deduce this quantity.
%
%   Below are two tables that demonstrate by example the relationship
%   between rectangular and latitudinal coordinates.
%
%   Listed in the first table (under r, longitude and latitude ) are
%   latitudinal coordinate triples that approximately represent
%   points whose rectangular coordinates are taken from the set
%   {-1, 0, 1}.  (Angular quantities are given in degrees.)
%
%
%    r       longitude  latitude      rectan(1)  rectan(2) rectan(3)
%   ----------------------------      -------------------------------
%    0.0000    0.0000    0.0000         0.0000     0.0000   0.0000
%    1.0000    0.0000    0.0000         1.0000     0.0000   0.0000
%    1.0000   90.0000    0.0000         0.0000     1.0000   0.0000
%    1.0000    0.0000   90.0000         0.0000     0.0000   1.0000
%    1.0000  180.0000    0.0000        -1.0000     0.0000   0.0000
%    1.0000  -90.0000    0.0000         0.0000    -1.0000   0.0000
%    1.0000    0.0000  -90.0000         0.0000     0.0000  -1.0000
%    1.4142   45.0000    0.0000         1.0000     1.0000   0.0000
%    1.4142    0.0000   45.0000         1.0000     0.0000   1.0000
%    1.4142   90.0000   45.0000         0.0000     1.0000   1.0000
%    1.7320   45.0000   35.2643         1.0000     1.0000   1.0000
%
%
%   This routine is related to the CSPICE routine cspice_latrec, which
%   accepts a radius, longitude, and latitude as inputs and produces
%   equivalent rectangular coordinates as outputs.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine srfrec_c.
%
%   MICE.REQ
%   KERNEL.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 01-DEC-2016, EDW (JPL)
%
%-Index_Entries
%
%   convert bodyfixed latitudinal coordinates to rectangular
%   convert surface latitudinal coordinates to rectangular
%   surface point latitudinal coordinates to rectangular
%
%-&

function [rectan] = cspice_srfrec(body, longitude, latitude)

   switch nargin
      case 3

         body      = zzmice_int(body);
         longitude = zzmice_dp(longitude);
         latitude  = zzmice_dp(latitude);

      otherwise

         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_srfrec(body, _longitude_, _latitude_)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice('srfrec_c', body,longitude,latitude);
   catch
      rethrow(lasterror)
   end

