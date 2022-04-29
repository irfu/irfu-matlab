%-Abstract
%
%   CSPICE_RECGEO converts rectangular coordinates to geodetic
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
%      rectan   the array(s) containing the rectangular coordinates of the
%               position or set of positions.
%
%               [3,n] = size(rectan); double = class(rectan)
%
%      re       the value describing equatorial radius of the body
%               of interest.
%
%               [1,1] = size(re); double = class(re)
%
%      f        the value describing flattening coefficient of the body,
%               a dimensionless value defined as:
%
%                    equatorial_radius - polar_radius
%                    --------------------------------
%                           equatorial_radius
%
%               [1,1] = size(f); double = class(f)
%
%   the call:
%
%      [lon, lat, alt] = cspice_recgeo( rectan, re, f )
%
%   returns:
%
%      lon      the value(s) describing the geodetic longitude
%               measured in radians.
%
%               [1,n] = size(lon); double = class(lon)
%
%      lat      the value(s) describing the geodetic latitude
%               measured in radians.
%
%               [1,n] = size(lat); double = class(lat)
%
%      alt      the value(s) describing the altitude above the
%               reference spheroid.
%
%               [1,n] = size(alt); double = class(alt)
%
%               `lon', `lat', and `alt' return with the same vectorization
%               measure, N, as `rectan'.
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
%   1) Find the geodetic coordinates of the point having Earth
%      rectangular coordinates:
%
%         X (km) =  -2541.748162
%         Y (km) =   4780.333036
%         Z (km) =   3360.428190
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model and orientation data for the Earth.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function recgeo_ex1()
%
%         %
%         % Load a PCK file containing a triaxial
%         % ellipsoidal shape model and orientation
%         % data for the Earth.
%         %
%         cspice_furnsh( 'pck00010.tpc' )
%
%         %
%         % Retrieve the triaxial radii of the earth
%         %
%         radii = cspice_bodvrd( 'EARTH', 'RADII', 3 );
%
%         %
%         % Calculate the flatness coefficient. Set a bodyfixed
%         % position.
%         %
%         flat = (radii(1) - radii(3))/radii(1);
%         x    = [ -2541.748162; 4780.333036; 3360.428190];
%
%         [ lon, lat, alt] = cspice_recgeo( x, radii(1), flat );
%
%         %
%         % Output, convert the angular values to degrees.
%         %
%         lon = lon * cspice_dpr;
%         lat = lat * cspice_dpr;
%
%         %
%         % Print out the results.
%         %
%         fprintf( 'Rectangular coordinates in km (x, y, z)\n' )
%         fprintf( '   %13.6f %13.6f %13.6f\n', x );
%         fprintf( 'Geodetic coordinates in deg and km (lon, lat, alt)\n' )
%         fprintf( '   %13.6f %13.6f %13.6f\n', lon, lat, alt );
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
%      Rectangular coordinates in km (x, y, z)
%          -2541.748162   4780.333036   3360.428190
%      Geodetic coordinates in deg and km (lon, lat, alt)
%            118.000000     31.999957      0.001916
%
%
%   2) Create a table showing a variety of rectangular coordinates
%      and the corresponding Earth geodetic coordinates. The
%      values are computed using the equatorial radius of the Clark
%      66 spheroid and the Clark 66 flattening factor:
%
%         radius: 6378.2064
%         flattening factor: 1./294.9787
%
%      Note: the values shown above may not be current or suitable
%            for your application.
%
%
%      Corresponding rectangular and geodetic coordinates are
%      listed to three decimal places. Output angles are in degrees.
%
%
%      Example code begins here.
%
%
%      function recgeo_ex2()
%
%         %
%         % Using the equatorial radius of the Clark66 spheroid
%         % (CLARKR = 6378.2064 km) and the Clark 66 flattening
%         % factor (CLARKF = 1.0 / 294.9787 ) convert from
%         % body fixed rectangular coordinates.
%         %
%         CLARKR = 6378.2064;
%         CLARKF = 1./294.9787;
%
%         x = [ [    0.0,         0.0,         0.0     ]',                 ...
%               [ 6378.20640,     0.0,         0.0     ]',                 ...
%               [    0.0,      6378.20640,     0.0     ]',                 ...
%               [    0.0,         0.0,      6378.20640 ]',                 ...
%               [-6378.20640,     0.0,         0.0     ]',                 ...
%               [    0.0,     -6378.20640,     0.0     ]',                 ...
%               [    0.0,         0.0,     -6378.20640 ]',                 ...
%               [ 6378.20640,  6378.20640,     0.0     ]',                 ...
%               [ 6378.20640,     0.0,      6378.20640 ]',                 ...
%               [    0.0,      6378.20640,  6378.20640 ]',                 ...
%               [ 6378.20640,  6378.20640,  6378.20640 ]' ];
%
%         [lon, lat, alt] = cspice_recgeo(  x, CLARKR, CLARKF );
%
%         %
%         % Output, convert the angular values to degrees.
%         %
%         lon = lon * cspice_dpr;
%         lat = lat * cspice_dpr;
%
%         %
%         % Output banner.
%         %
%         disp(' rectan[1]  rectan[2]  rectan[3]    lon      lat       alt')
%         fprintf([' ---------  ---------  --------- ',                    ...
%                  ' -------  -------  ---------\n' ])
%
%         for i = 1:10
%            fprintf( '%10.3f %10.3f %10.3f', x(:,i) )
%            fprintf( ' %8.3f %8.3f %10.3f\n', lon(i), lat(i), alt(i));
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%       rectan[1]  rectan[2]  rectan[3]    lon      lat       alt
%       ---------  ---------  ---------  -------  -------  ---------
%           0.000      0.000      0.000    0.000   90.000  -6356.584
%        6378.206      0.000      0.000    0.000    0.000      0.000
%           0.000   6378.206      0.000   90.000    0.000      0.000
%           0.000      0.000   6378.206    0.000   90.000     21.623
%       -6378.206      0.000      0.000  180.000    0.000      0.000
%           0.000  -6378.206      0.000  -90.000    0.000      0.000
%           0.000      0.000  -6378.206    0.000  -90.000     21.623
%        6378.206   6378.206      0.000   45.000    0.000   2641.940
%        6378.206      0.000   6378.206    0.000   45.137   2652.768
%           0.000   6378.206   6378.206   90.000   45.137   2652.768
%
%
%-Particulars
%
%   Given the body-fixed rectangular coordinates of a point, and the
%   constants describing the reference spheroid,  this routine
%   returns the geodetic coordinates of the point. The body-fixed
%   rectangular frame is that having the x-axis pass through the
%   0 degree latitude 0 degree longitude point. The y-axis passes
%   through the 0 degree latitude 90 degree longitude. The z-axis
%   passes through the 90 degree latitude point. For some bodies
%   this coordinate system may not be a right-handed coordinate
%   system.
%
%-Exceptions
%
%   1)  If the equatorial radius is non-positive, the error
%       SPICE(VALUEOUTOFRANGE) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the flattening coefficient is greater than or equal to one,
%       the error SPICE(VALUEOUTOFRANGE) is signaled by a routine in
%       the call tree of this routine.
%
%   3)  For points inside the reference ellipsoid, the nearest
%       point on the ellipsoid to `rectan' may not be unique, so
%       latitude may not be well-defined.
%
%   4)  If any of the input arguments, `rectan', `re' or `f', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `rectan', `re' or `f', is not
%       of the expected type, or it does not have the expected
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
%   [1]  R. Bate, D. Mueller, and J. White, "Fundamentals of
%        Astrodynamics," Dover Publications Inc., 1971.
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
%       Edited the header to comply with NAIF standard. Split the existing
%       code example into two separate examples.
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
%   rectangular to geodetic
%
%-&

function [lon, lat, alt] = cspice_recgeo(rectan, re, f)

   switch nargin
      case 3

         rectan = zzmice_dp(rectan);
         re     = zzmice_dp(re);
         f      = zzmice_dp(f);

      otherwise
         error ( ['Usage: [_lon_, _lat_, _alt_] = '...
                  'cspice_recgeo(_rectan(3)_, re, f)' ] )
   end

   %
   % Call the MEX library.
   %
   try
      [lon, lat, alt] = mice( 'recgeo_c', rectan, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end






