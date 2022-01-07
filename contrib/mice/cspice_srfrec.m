%-Abstract
%
%   CSPICE_SRFREC converts planetocentric latitude and longitude of a surface
%   point on a specified body to rectangular coordinates.
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
%      body     the NAIF integer code of an extended body on which a surface
%               point of interest is located.
%
%               [1,1] = size(body); int32 = class(body)
%
%               The body is modeled as a triaxial ellipsoid.
%
%      lon      the longitude of the input point(s).
%
%               [1,n] = size(lon); double = class(lon)
%
%               This is the angle between the prime meridian and the
%               meridian containing the point. The direction of increasing
%               longitude is from the +X axis towards the +Y axis.
%
%               Longitude is measured in radians. On input, the
%               range of longitude is unrestricted.
%
%      lat      the latitude of the input point(s).
%
%               [1,n] = size(lat); double = class(lat)
%
%               This is the angle from the XY plane of the ray from the
%               origin through the point.
%
%               Latitude is measured in radians. On input, the range
%               of latitude is unrestricted.
%
%   the call:
%
%      [rectan] = cspice_srfrec( body, lon, lat )
%
%   returns:
%
%      rectan   the rectangular coordinates of the input surface point(s).
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               Units are the same as those used to define the radii of
%               `body'. Normally, these units are km.
%
%               `rectan' returns with the vectorization measure, N, as
%               `lon', and `lat'.
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
%   1) Find the rectangular coordinates of the point
%
%         100 degrees planetocentric longitude
%         -35 degrees planetocentric latitude
%
%      on the Earth; then convert these coordinates back to
%      latitudinal coordinates. We should be able to recover
%      our original longitude and latitude values.
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model and orientation data for the Earth.
%
%         pck00008.tpc
%
%
%      Example code begins here.
%
%
%      function srfrec_ex1()
%
%
%         %
%         % NAIF ID for our body of interest.
%         %
%         EARTH =  399;
%
%         %
%         % Load the kernel pool with a PCK file that contains
%         % values for the radii of the Earth.
%         %
%         cspice_furnsh( 'pck00008.tpc' )
%
%         %
%         % Find `x', the rectangular coordinates of the surface point
%         % defined by `lat' and `lon'.  The NAIF integer code for
%         % the Earth is 399. (See the NAIF_IDS required reading file
%         % for the complete set of codes.)
%         %
%         lon =  100.;
%         lat =   35.;
%
%         fprintf( 'Original latitudinal coordinates: \n' )
%         fprintf( '                 Longitude (deg): %f\n', lon )
%         fprintf( '                 Latitude  (deg): %f\n\n', lat )
%
%
%         %
%         % Convert angles to radians forr input to cspice_srfrec.
%         %
%         x = cspice_srfrec( EARTH, lon*cspice_rpd(), lat*cspice_rpd() );
%
%         fprintf( 'Rectangular coordinates: \n')
%         fprintf( '                 X (km): %f\n', x(1) )
%         fprintf( '                 Y (km): %f\n', x(2) )
%         fprintf( '                 Z (km): %f\n\n', x(3) )
%
%
%         %
%         %
%         % Now try to recover the original latitudinal coordinates
%         % from the rectangular coordinates found by cspice_srfrec.
%         %
%         [radius, lon1, lat1] = cspice_reclat( x);
%
%         %
%         % Convert angles back to degree for display.
%         %
%         fprintf( 'Latitudinal coordinates recovered from \n' )
%         fprintf( 'rectangular coordinates: \n' )
%         fprintf( '                 Longitude (deg): %f\n',  ...
%                                          lon1*cspice_dpr() )
%         fprintf( '                 Latitude  (deg): %f\n',  ...
%                                          lat1*cspice_dpr() )
%
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear()
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Original latitudinal coordinates:
%                       Longitude (deg): 100.000000
%                       Latitude  (deg): 35.000000
%
%      Rectangular coordinates:
%                       X (km): -906.249429
%                       Y (km): 5139.595909
%                       Z (km): 3654.300840
%
%      Latitudinal coordinates recovered from
%      rectangular coordinates:
%                       Longitude (deg): 100.000000
%                       Latitude  (deg): 35.000000
%
%
%   2) Create a table showing a variety of Earth latitudinal coordinates
%      and the corresponding rectangular coordinates.
%
%      Corresponding latitudinal and rectangular coordinates are
%      listed to four decimal places.
%
%      Use the PCK file from example 1 above.
%
%
%      Example code begins here.
%
%
%      function srfrec_ex2()
%
%         %
%         % NAIF ID for our body of interest.
%         %
%         EARTH =  399;
%
%         %
%         % Load the kernel pool with a PCK file that contains
%         % values for the radii of the Earth.
%         %
%         cspice_furnsh( 'pck00008.tpc' )
%
%         %
%         % Define ten sets of latitudinal coordinates.
%         %
%         longitudes = [ 0., 90., 0. 180., -90., ...
%                                        0., 45., 0., 90., 45. ];
%         latitudes  = [ 0., 0., 90., 0., 0.,    ...
%                                        -90., 0., 45., 45., 35.2643 ];
%
%         %
%         % Convert angles to radians for input to cspice_srfrec.
%         %
%         rectan = cspice_srfrec( EARTH, longitudes*cspice_rpd(), ...
%                                        latitudes*cspice_rpd() );
%
%         %
%         % Create an array of values for output.
%         %
%         output = [ longitudes; latitudes; rectan ];
%
%         %
%         % Output banner.
%         %
%         disp('  longitude  latitude       x         y           z   ')
%         disp('  --------   --------   --------   --------   --------')
%
%         txt = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f\n', output );
%         disp( txt )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear()
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%        longitude  latitude       x         y           z
%        --------   --------   --------   --------   --------
%          0.0000     0.0000  6378.1400     0.0000     0.0000
%         90.0000     0.0000     0.0000  6378.1400     0.0000
%          0.0000    90.0000     0.0000     0.0000  6356.7500
%        180.0000     0.0000 -6378.1400     0.0000     0.0000
%        -90.0000     0.0000     0.0000 -6378.1400     0.0000
%          0.0000   -90.0000     0.0000     0.0000 -6356.7500
%         45.0000     0.0000  4510.0260  4510.0260     0.0000
%          0.0000    45.0000  4502.4445     0.0000  4502.4445
%         90.0000    45.0000     0.0000  4502.4445  4502.4445
%         45.0000    35.2643  3678.2946  3678.2946  3678.2824
%
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
%   point. In this case, the distance from the central reference
%   point is not required as an input because the fact that the
%   point is on the body's surface allows one to deduce this quantity.
%
%   Below are two tables that demonstrate by example the relationship
%   between rectangular and latitudinal coordinates.
%
%   Listed in the first table (under R, `lon' and `lat') are
%   latitudinal coordinate triples that approximately represent
%   points whose rectangular coordinates are taken from the set
%   {-1, 0, 1}. (Angular quantities are given in degrees.)
%
%
%        R           lon       lat    rectan(1)   rectan(2)  rectan(3)
%       --------------------------    --------------------------------
%       0.0000    0.0000    0.0000      0.0000      0.0000     0.0000
%       1.0000    0.0000    0.0000      1.0000      0.0000     0.0000
%       1.0000   90.0000    0.0000      0.0000      1.0000     0.0000
%       1.0000    0.0000   90.0000      0.0000      0.0000     1.0000
%       1.0000  180.0000    0.0000     -1.0000      0.0000     0.0000
%       1.0000  -90.0000    0.0000      0.0000     -1.0000     0.0000
%       1.0000    0.0000  -90.0000      0.0000      0.0000    -1.0000
%       1.4142   45.0000    0.0000      1.0000      1.0000     0.0000
%       1.4142    0.0000   45.0000      1.0000      0.0000     1.0000
%       1.4142   90.0000   45.0000      0.0000      1.0000     1.0000
%       1.7320   45.0000   35.2643      1.0000      1.0000     1.0000
%
%
%   This routine is related to the Mice routine cspice_latrec, which
%   accepts a radius, longitude, and latitude as inputs and produces
%   equivalent rectangular coordinates as outputs.
%
%-Exceptions
%
%   1)  If radii for `body' are not found in the kernel pool, an error
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If the size of the `body' body radii kernel variable is not
%       three, an error is signaled by a routine in the call tree of
%       this routine.
%
%   3)  If any of the three `body' body radii is less-than or equal to
%       zero, an error is signaled by a routine in the call tree of
%       this routine.
%
%   4)  If any of the input arguments, `body', `lon' or `lat', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `body', `lon' or `lat', is not
%       of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   6)  If the input vectorizable arguments `lon' and `lat' do not
%       have the same measure of vectorization (N), an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  A PCK text kernel containing the body radius definitions
%       required by this routine must be loaded into the kernel
%       pool prior to any calls to this routine.
%
%-Required_Reading
%
%   KERNEL.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
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
%       Changed input argument names "longitude" and "latitude" to "lon" and
%       "lat". Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections. Edited
%       the header to comply with NAIF standard.
%
%       Added example's task description.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 01-DEC-2016 (EDW)
%
%-Index_Entries
%
%   convert bodyfixed latitudinal coordinates to rectangular
%   convert surface latitudinal coordinates to rectangular
%   surface point latitudinal coordinates to rectangular
%
%-&

function [rectan] = cspice_srfrec(body, lon, lat)

   switch nargin
      case 3

         body = zzmice_int(body);
         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);

      otherwise

         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_srfrec(body, _lon_, _lat_)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice('srfrec_c', body, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end

