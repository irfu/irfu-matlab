%-Abstract
%
%   CSPICE_LATREC converts latitudinal coordinates to rectangular
%   (Cartesian) coordinates.
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
%      [rectan] = cspice_latrec( radius, lon, lat )
%
%   returns:
%
%      rectan   the array(s) containing the rectangular coordinates of the
%               position or set of positions
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               `rectan' returns with the same units associated with `radius'.
%
%               `rectan' returns with the vectorization measure, N, as
%               `radius', `lon', and `lat'.
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
%         File name: latrec_ex1.tm
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
%      function latrec_ex1()
%
%         %
%         % Load an SPK and leapseconds kernels.
%         %
%         cspice_furnsh( 'latrec_ex1.tm' )
%
%         %
%         % Convert the time to ET.
%         %
%         et = cspice_str2et( '2017 Mar 20' );
%
%         %
%         % Retrieve the position of the moon seen from earth at 'et'
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
%         % Convert the position vector 'pos' to latitudinal
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
%   2) Create a table showing a variety of latitudinal coordinates
%      and the corresponding rectangular coordinates.
%
%      Corresponding latitudinal and rectangular coordinates are
%      listed to four decimal places. Input angles are in degrees.
%
%
%      Example code begins here.
%
%
%      function latrec_ex2()
%
%         %
%         % Define eleven sets of latitudinal coordinates.
%         %
%         r         = [ 0., 1., 1., 1., 1., 1., 1., ...
%                       sqrt(2), sqrt(2), sqrt(2), sqrt(3) ];
%         longitude = [ 0., 0., 90., 0. 180., -90., ...
%                       0., 45., 0., 90., 45. ];
%         latitude  = [ 0., 0., 0., 90., 0., 0.,    ...
%                       -90., 0., 45., 45., 35.2643 ];
%
%         %
%         % ...convert the latitudinal coordinates to rectangular coordinates
%         %
%         longitude = longitude * cspice_rpd;
%         latitude  = latitude  * cspice_rpd;
%
%         rectan = cspice_latrec(r, longitude, latitude);
%
%         %
%         % Loop over each set of coordinates for output, convert `longitude'
%         % and `latitude' to degrees...
%         %
%         longitude = longitude * cspice_dpr;
%         latitude  = latitude  * cspice_dpr;
%
%         %
%         % Create an array of values for output.
%         %
%         output = [ r; longitude; latitude; rectan ];
%
%         %
%         % Output banner.
%         %
%         disp('     r        lon       lat     rect(1)   rect(2)   rect(3)')
%         disp('  -------  --------  --------   -------   -------   -------')
%
%         txt = sprintf( '%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', output );
%         disp( txt )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%           r        lon       lat     rect(1)   rect(2)   rect(3)
%        -------  --------  --------   -------   -------   -------
%          0.000     0.000     0.000     0.000     0.000     0.000
%          1.000     0.000     0.000     1.000     0.000     0.000
%          1.000    90.000     0.000     0.000     1.000     0.000
%          1.000     0.000    90.000     0.000     0.000     1.000
%          1.000   180.000     0.000    -1.000     0.000     0.000
%          1.000   -90.000     0.000     0.000    -1.000     0.000
%          1.000     0.000   -90.000     0.000     0.000    -1.000
%          1.414    45.000     0.000     1.000     1.000     0.000
%          1.414     0.000    45.000     1.000     0.000     1.000
%          1.414    90.000    45.000     0.000     1.000     1.000
%          1.732    45.000    35.264     1.000     1.000     1.000
%
%
%-Particulars
%
%   This routine returns the rectangular coordinates of a point
%   whose position is input in latitudinal coordinates.
%
%   Latitudinal coordinates are defined by a distance from a central
%   reference point, an angle from a reference meridian, and an angle
%   above the equator of a sphere centered at the central reference
%   point.
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Changed input argument names "longitude" and "latitude" to "lon" and
%       "lat".
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
%   latitudinal to rectangular coordinates
%
%-&

function [rectan] = cspice_latrec(radius, lon, lat)

   switch nargin
      case 3

         radius = zzmice_dp(radius);
         lon    = zzmice_dp(lon);
         lat    = zzmice_dp(lat);

      otherwise

         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_latrec(_radius_, _lon_, _lat_)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice('latrec_c', radius, lon, lat);
   catch spiceerr
      rethrow(spiceerr)
   end

