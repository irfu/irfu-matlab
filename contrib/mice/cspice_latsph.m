%-Abstract
%
%   CSPICE_LATSPH converts latitudinal coordinates to spherical
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
%   the call:
%
%      [rho, colat, slon] = cspice_latsph( radius, lon, lat)
%
%   returns:
%
%      rho      the value(s) describing the distance of the position
%               from the origin.
%
%               [1,n] = size(rho); double = class(rho)
%
%      colat    the value(s) describing the angle between the point and the
%               positive z-axis, measured in radians (also referred to
%               as the polar angle).
%
%               [1,n] = size(colat); double = class(colat)
%
%      slon     the value(s) describing the angle of the projection of the
%               point to the XY plane from the positive X-axis, measured
%               in radians, with range:
%
%                   -pi < slon <= pi
%
%               The positive Y-axis is at longitude PI/2 radians.
%
%               [1,n] = size(slon); double = class(slon)
%
%
%               The argument `rho' returns in the same units associated
%               with `radius'.
%
%               `rho', `colat', and `slon' return with the same vectorization
%               measure, N, as `radius', `lon', and `lat'.
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
%   1) Co-latitude is obtained by subtracting latitude from cspice_halfpi
%      Radius and longitude mean the same thing in both latitudinal
%      and spherical coordinates. The table below lists `lat' and
%      corresponding `colat' in terms of degrees.
%
%           lat     colat
%          -----    -----
%             0        90
%            20        70
%            45        45
%           -30       120
%            90         0
%           -45       135
%
%
%   2) Compute the latitudinal coordinates of the position of the Moon
%      as seen from the Earth, and convert them to spherical and
%      rectangular coordinates.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: latsph_ex2.tm
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
%      function latsph_ex2()
%
%         %
%         % Load an SPK, leapseconds, and PCK kernel set.
%         %
%         cspice_furnsh( 'latsph_ex2.tm' )
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
%         fprintf( '   X            (km): %20.8f\n', pos(1) )
%         fprintf( '   Y            (km): %20.8f\n', pos(2) )
%         fprintf( '   Z            (km): %20.8f\n', pos(3) )
%
%         %
%         % Convert the position vector `pos' to latitudinal
%         % coordinates.
%         %
%         [radius, lon, lat] = cspice_reclat(pos);
%         fprintf( '\n' )
%         fprintf( 'Latitudinal coordinates:\n' )
%         fprintf( '   Radius       (km): %20.8f\n', radius )
%         fprintf( '   Longitude   (deg): %20.8f\n', lon * cspice_dpr )
%         fprintf( '   Latitude    (deg): %20.8f\n', lat * cspice_dpr )
%
%         %
%         % Convert the latitudinal coords to spherical.
%         %
%         [rho, colat, slon] = cspice_latsph( radius, lon, lat);
%         fprintf( '\n' )
%         fprintf( 'Spherical coordinates:\n' )
%         fprintf( '   Radius       (km): %20.8f\n', rho )
%         fprintf( '   Polar Angle (deg): %20.8f\n', colat * cspice_dpr )
%         fprintf( '   Longitude   (deg): %20.8f\n', slon  * cspice_dpr )
%
%         %
%         % Convert the spherical to rectangular.
%         %
%         [rectan] = cspice_sphrec(rho, colat, slon);
%         fprintf( '\n' )
%         fprintf( 'Rectangular coordinates from cspice_sphrec:\n' )
%         fprintf( '   X            (km): %20.8f\n', rectan(1) )
%         fprintf( '   Y            (km): %20.8f\n', rectan(2) )
%         fprintf( '   Z            (km): %20.8f\n', rectan(3) )
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
%         X            (km):      -55658.44323296
%         Y            (km):     -379226.32931475
%         Z            (km):     -126505.93063865
%
%      Latitudinal coordinates:
%         Radius       (km):      403626.33912495
%         Longitude   (deg):         -98.34959789
%         Latitude    (deg):         -18.26566077
%
%      Spherical coordinates:
%         Radius       (km):      403626.33912495
%         Polar Angle (deg):         108.26566077
%         Longitude   (deg):         -98.34959789
%
%      Rectangular coordinates from cspice_sphrec:
%         X            (km):      -55658.44323296
%         Y            (km):     -379226.32931475
%         Z            (km):     -126505.93063865
%
%
%   3) Create a table showing a variety of latitudinal coordinates
%      and the corresponding spherical coordinates.
%
%      Corresponding latitudinal and spherical coordinates are
%      listed to four decimal places. Input and output angles are
%      in degrees.
%
%
%      Example code begins here.
%
%
%      function latsph_ex3()
%
%         %
%         % Define six sets of latitudinal coordinates, `lon' and `lat'
%         % expressed in degrees - converted to radians by use
%         % of cspice_rpd.
%         %
%         rad = [ 1.,  1., sqrt(2.), sqrt(2.),   1.,  0. ];
%         lon = [ 0., 90.,     180.,     180., 180., 33. ] * cspice_rpd;
%         lat = [ 0.,  0.,      45.,      -45., 90.,  0. ] * cspice_rpd;
%
%         %
%         % ...convert the latitudinal coordinates to spherical coordinates
%         %
%         [rho, colat, slon] = cspice_latsph(rad, lon, lat);
%
%         %
%         % ...convert angular measure to degrees.
%         %
%         colat = colat * cspice_dpr;
%         lon   = lon   * cspice_dpr;
%         slon  = slon  * cspice_dpr;
%         lat   = lat   * cspice_dpr;
%
%         %
%         % Output banner.
%         %
%         disp('   radius     lon       lat        r       colat      slon ')
%         disp('  -------  --------  --------   -------   -------   -------')
%
%         %
%         % Create an array of values for output.
%         %
%         output = [ rad; lon; lat; rho; colat; slon];
%         txt    = sprintf( '%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', output );
%         disp( txt )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%         radius     lon       lat        r       colat      slon
%        -------  --------  --------   -------   -------   -------
%          1.000     0.000     0.000     1.000    90.000     0.000
%          1.000    90.000     0.000     1.000    90.000    90.000
%          1.414   180.000    45.000     1.414    45.000   180.000
%          1.414   180.000   -45.000     1.414   135.000   180.000
%          1.000   180.000    90.000     1.000     0.000   180.000
%          0.000    33.000     0.000     0.000    90.000    33.000
%
%
%-Particulars
%
%   This routine returns the spherical coordinates of a point
%   whose position is input in latitudinal coordinates.
%
%   Latitudinal coordinates are defined by a distance from a central
%   reference point, an angle from a reference meridian, and an angle
%   above the equator of a sphere centered at the central reference
%   point.
%
%   Spherical coordinates are defined by a distance from a central
%   reference point, an angle from a reference meridian, and an angle
%   from the z-axis.
%
%-Exceptions
%
%   1)  If any of the input arguments, `radius', `lon' or `lat', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `radius', `lon' or `lat', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   3)  If the input vectorizable arguments `radius', `lon' and `lat'
%       do not have the same measure of vectorization (N), an error is
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
%       Edited the header to comply with NAIF standard. Added
%       meta-kernel to example #2. Updated code example #2 to produce
%       formatted output and added a call to cspice_kclear. Added the
%       problem statement to existing examples and added third example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Changed the input arguments "longitude", "latitude" to "lon" and
%       "lat", and the output arguments "radius" and "lon" to "rho" and
%       "slon" for consistency with other functions.
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
%   latitudinal to spherical coordinates
%
%-&

function [rho, colat, slon] = cspice_latsph( radius, lon, lat)

   switch nargin
      case 3

         radius = zzmice_dp(radius);
         lon    = zzmice_dp(lon);
         lat    = zzmice_dp(lat);

      otherwise

         error ( ['Usage: [_rho_, _colat_, _slon_] = '...
                  'cspice_latsph( _radius_, _lon_, _lat_)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [rho, colat, slon] = mice('latsph_c', radius, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end





