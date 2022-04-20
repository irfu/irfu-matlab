%-Abstract
%
%   CSPICE_GEOREC converts geodetic coordinates to rectangular
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
%      lon      the value(s) describing the geodetic longitude measured
%               in radians.
%
%               [1,n] = size(lon); double = class(lon)
%
%      lat      the value(s) describing the geodetic latitude measured
%               in radians.
%
%               [1,n] = size(lat); double = class(lat)
%
%      alt      the value(s) describing the altitude above the reference
%               spheroid.
%
%               [1,n] = size(alt); double = class(alt)
%
%      re       the equatorial radius of the body of interest.
%
%               [1,1] = size(re); double = class(re)
%
%      f        the flattening coefficient of the body, a dimensionless
%               value defined as:
%
%                  equatorial_radius - polar_radius
%                  --------------------------------
%                         equatorial_radius
%
%               [1,1] = size(f); double = class(f)
%
%   the call:
%
%      [rectan] = cspice_georec( lon, lat, alt, re, f )
%
%   returns:
%
%      rectan   the array(s) containing the rectangular coordinates of the
%               position or set of positions.
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               `rectan' returns with the same units associated with
%               `alt' and `re'.
%
%               `rectan' returns with the same vectorization measure,
%               N, as `lon', `lat', and `alt'
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
%   1) Find the rectangular coordinates of the point having Earth
%      geodetic coordinates:
%
%         lon (deg) =  118.0
%         lat (deg) =   32.0
%         alt (km)  =    0.0
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
%      function georec_ex1()
%
%         %
%         % Load a PCK file containing a triaxial
%         % ellipsoidal shape model and orientation
%         % data for the Earth.
%         %
%         cspice_furnsh( 'pck00010.tpc' );
%
%         %
%         % Retrieve the triaxial radii of the Earth
%         %
%         [radii] = cspice_bodvrd( 'EARTH', 'RADII', 3 );
%
%         %
%         % Compute flattening coefficient.
%         %
%         re =  radii(1);
%         rp =  radii(3);
%         f  =  ( re - rp ) / re;
%
%         %
%         % Set a geodetic position.
%         %
%         lon = 118.0 * cspice_rpd;
%         lat =  30.0 * cspice_rpd;
%         alt =   0.0;
%
%         %
%         % Do the conversion.
%         %
%         [rectan] = cspice_georec( lon, lat, alt, radii(1), f );
%
%         fprintf( 'Geodetic coordinates in deg and km (lon, lat, alt)\n' )
%         fprintf( '%14.6f %13.6f %13.6f\n',                               ...
%                  lon * cspice_dpr, lat * cspice_dpr, alt )
%         fprintf( 'Rectangular coordinates in km (x, y, z)\n' )
%         fprintf( '%14.6f %13.6f %13.6f\n',                               ...
%                  rectan(1), rectan(2), rectan(3) )
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
%      Geodetic coordinates in deg and km (lon, lat, alt)
%          118.000000     30.000000      0.000000
%      Rectangular coordinates in km (x, y, z)
%        -2595.359123   4881.160589   3170.373523
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
%      listed to three decimal places. Input angles are in degrees.
%
%
%      Example code begins here.
%
%
%      function georec_ex2()
%
%         %
%         % Local parameters.
%         %
%         NREC = 11;
%
%         %
%         % Define the input geodetic coordinates. Angles in
%         % degrees.
%         %
%         lon = [ 0.0,   0.0, 90.0,  0.0, 180.0, -90.0,                    ...
%                 0.0,  45.0,  0.0, 90.0,  45.0 ];
%
%         lat = [  90.0, 0.0,  0.0,   90.0,    0.0,   0.0,                 ...
%                 -90.0, 0.0, 88.707, 88.707, 88.1713 ];
%
%         alt = [ -6356.5838, 0.0,     0.0,        0.0,        0.0,   0.0, ...
%                     0.0,    0.0, -6355.5725, -6355.5725, -6355.5612 ];
%
%         %
%         % Using the equatorial radius of the Clark66 spheroid
%         % (clarkr = 6378.2064 km) and the Clark 66 flattening
%         % factor (clarkf = 1.0 / 294.9787 ) convert from
%         % body fixed rectangular coordinates.
%         %
%         clarkr = 6378.2064;
%         clarkf = 1.0 / 294.9787;
%
%         %
%         % Print the banner.
%         %
%         fprintf( [ '   lon      lat       alt     rectan(1)  rectan(2)', ...
%                    '  rectan(3)\n' ]                                     )
%         fprintf( [ ' -------  -------  ---------  ---------  ---------', ...
%                    '  ---------\n' ]                                     )
%
%         %
%         % Do the conversion.
%         %
%         rlon     = lon * cspice_rpd;
%         rlat     = lat * cspice_rpd;
%
%         [rectan] = cspice_georec( rlon, rlat, alt, clarkr, clarkf );
%
%         for i=1:NREC
%
%            fprintf( '%8.3f %8.3f %10.3f', lon(i), lat(i), alt(i) )
%            fprintf( '%11.3f %10.3f %10.3f\n', rectan(:,i) )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%         lon      lat       alt     rectan(1)  rectan(2)  rectan(3)
%       -------  -------  ---------  ---------  ---------  ---------
%         0.000   90.000  -6356.584      0.000      0.000      0.000
%         0.000    0.000      0.000   6378.206      0.000      0.000
%        90.000    0.000      0.000      0.000   6378.206      0.000
%         0.000   90.000      0.000      0.000      0.000   6356.584
%       180.000    0.000      0.000  -6378.206      0.000      0.000
%       -90.000    0.000      0.000      0.000  -6378.206      0.000
%         0.000  -90.000      0.000      0.000      0.000  -6356.584
%        45.000    0.000      0.000   4510.073   4510.073      0.000
%         0.000   88.707  -6355.573      1.000      0.000      1.000
%        90.000   88.707  -6355.573      0.000      1.000      1.000
%        45.000   88.171  -6355.561      1.000      1.000      1.000
%
%
%-Particulars
%
%   Given the geodetic coordinates of a point, and the constants
%   describing the reference spheroid,  this routine returns the
%   bodyfixed rectangular coordinates of the point. The bodyfixed
%   rectangular frame is that having the X-axis pass through the
%   0 degree latitude 0 degree longitude point. The Y-axis passes
%   through the 0 degree latitude 90 degree longitude. The Z-axis
%   passes through the 90 degree latitude point. For some bodies
%   this coordinate system may not be a right-handed coordinate
%   system.
%
%-Exceptions
%
%   1)  If the flattening coefficient is greater than or equal to one,
%       the error SPICE(VALUEOUTOFRANGE) is signaled by a routine in
%       the call tree of this routine.
%
%   2)  If the equatorial radius is less than or equal to zero, the
%       error SPICE(VALUEOUTOFRANGE) is signaled by a routine in the
%       call tree of this routine.
%
%   3)  If any of the input arguments, `lon', `lat', `alt', `re' or
%       `f', is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   4)  If any of the input arguments, `lon', `lat', `alt', `re' or
%       `f', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%   5)  If the input vectorizable arguments `lon', `lat' and `alt' do
%       not have the same measure of vectorization (N), an error is
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
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
%   -Mice Version 1.0.1, 06-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   geodetic to rectangular coordinates
%
%-&

function [rectan] = cspice_georec( lon, lat, alt, re, f)

   switch nargin
      case 5

         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
         alt  = zzmice_dp(alt);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);

      otherwise

         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_georec( _lon_, _lat_, _alt_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice( 'georec_c', lon, lat, alt, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end
