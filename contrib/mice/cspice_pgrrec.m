%-Abstract
%
%   CSPICE_PGRREC converts planetographic coordinates to
%   rectangular coordinates.
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
%      body     the name of the body with which the planetographic coordinate
%               system is associated.
%
%               [1,c1] = size(body); char = class(body)
%
%                  or
%
%               [1,1] = size(body); cell = class(body)
%
%               `body' is used by this routine to look up from the
%               kernel pool the prime meridian rate coefficient giving
%               the body's spin sense. See the -Files and -Particulars
%               header sections below for details.
%
%      lon      the planetographic longitude(s) of the input point(s).
%
%               [1,n] = size(lon); double = class(lon)
%
%               This is the angle between the prime meridian and the
%               meridian containing the input point. For bodies having
%               prograde (aka direct) rotation, the direction of increasing
%               longitude is positive west: from the +X axis of the
%               rectangular coordinate system toward the -Y axis. For bodies
%               having retrograde rotation, the direction of increasing
%               longitude is positive east: from the +X axis toward the +Y
%               axis.
%
%               The earth, moon, and sun are exceptions:
%               planetographic longitude is measured positive east for
%               these bodies.
%
%               The default interpretation of longitude by this
%               and the other planetographic coordinate conversion
%               routines can be overridden; see the discussion in
%               -Particulars below for details.
%
%               `lon' is measured in radians. On input, the range
%               of longitude is unrestricted.
%
%      lat      the planetographic latitude(s) of the input point(s).
%
%               [1,n] = size(lat); double = class(lat)
%
%               For a point P on the reference spheroid, this is the angle
%               between the XY plane and the outward normal vector at P. For
%               a point P not on the reference spheroid, the planetographic
%               latitude is that of the closest point to P on the spheroid.
%
%               `lat' is measured in radians. On input, the
%               range of latitude is unrestricted.
%
%      alt      the altitude(s) of point(s) above the reference spheroid.
%
%               [1,n] = size(alt); double = class(alt)
%
%               Units of `alt' must match those of `re'.
%
%      re       the equatorial radius of a reference spheroid.
%
%               [1,1] = size(re); double = class(re)
%
%               This spheroid is a volume of revolution: its horizontal
%               cross sections are circular. The shape of the spheroid is
%               defined by an equatorial radius `re' and a polar radius RP.
%               Units of `re' must match those of `alt'.
%
%      f        the flattening coefficient of the body, a
%               dimensionless value defined as:
%
%                  (re - rp) / re
%
%               where `rp' is the polar radius of the spheroid, and the
%               units of `rp' match those of `re'.
%
%               [1,1] = size(f); double = class(f)
%
%   the call:
%
%      [rectan] = cspice_pgrrec( body, lon, lat, alt, re, f )
%
%   returns:
%
%      rectan   the rectangular coordinates of the input point(s).
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               See the discussion below in the -Particulars header section
%               for details.
%
%               The units associated with `rectan' are those associated
%               with the inputs `alt' and `re'.
%
%               `rectan' returns with the same vectorization measure, N,
%               as `lon', `lat' and `alt'.
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
%   1) Find the rectangular coordinates of the point having Mars
%      planetographic coordinates:
%
%         longitude = 90 degrees west
%         latitude  = 45 degrees north
%         altitude  = 300 km
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model and orientation data for Mars.
%
%         pck00008.tpc
%
%
%      Example code begins here.
%
%
%      function pgrrec_ex1()
%
%         %
%         % Load a PCK file containing a triaxial
%         % ellipsoidal shape model and orientation
%         % data for Mars.
%         %
%         cspice_furnsh( 'pck00008.tpc' )
%
%         %
%         % Example 1: convert a single set of planetographic
%         %            coordinates to rectangular bodyfixed
%         %            coordinates.
%         %
%         % Look up the radii for Mars.  Although we
%         % omit it here, we could check the kernel pool
%         % to make sure the variable BODY499_RADII
%         % has three elements and numeric data type.
%         % If the variable is not present in the kernel
%         % pool, cspice_bodvrd will signal an error.
%         %
%         body = 'MARS';
%         radii = cspice_bodvrd( body, 'RADII', 3 );
%
%         %
%         %
%         % Calculate the flatness coefficient. Set a bodyfixed
%         % position vector, `x'.
%         %
%         re   = radii(1);
%         rp   = radii(3);
%         flat = ( re - rp ) / re;
%
%         % Set a longitude, latitude, altitude position.
%         % Note that we must provide longitude and
%         % latitude in radians.
%         %
%         lon  = 90. * cspice_rpd;
%         lat  = 45.  * cspice_rpd;
%         alt  = 3.d2;
%
%         %
%         % Do the conversion.
%         %
%         x = cspice_pgrrec( body, lon, lat, alt, re, flat );
%
%         %
%         % Output.
%         %
%         disp( 'Rectangular coordinates in km (x, y, z)' )
%         fprintf( '%9.3f   %9.3f   %9.3f\n', x' )
%
%         disp( 'Planetographic coordinates in degs and km (lon, lat, alt)' )
%         fprintf( '%9.3f   %9.3f   %9.3f\n', lon *cspice_dpr() ...
%                                           , lat *cspice_dpr() ...
%                                           , alt               )
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
%      Rectangular coordinates in km (x, y, z)
%          0.000   -2620.679    2592.409
%      Planetographic coordinates in degs and km (lon, lat, alt)
%         90.000      45.000     300.000
%
%
%   2) Create a table showing a variety of rectangular coordinates
%      and the corresponding Mars planetographic coordinates. The
%      values are computed using the reference spheroid having radii
%
%         Equatorial radius:    3396.190
%         Polar radius:         3376.200
%
%      Note: the values shown above may not be current or suitable
%             for your application.
%
%
%      Corresponding rectangular and planetographic coordinates are
%      listed to three decimal places.
%
%      Use the PCK file from example 1 above.
%
%
%      Example code begins here.
%
%
%      %
%      % Example 2: convert a vectorized set of planetographic coordinates
%      %            to rectangular bodyfixed coordinates.
%      %
%      function pgrrec_ex2()
%
%         %
%         % Load a PCK file containing a triaxial
%         % ellipsoidal shape model and orientation
%         % data for Mars.
%         %
%         cspice_furnsh( 'pck00008.tpc' )
%
%         % Look up the radii for Mars.  Although we
%         % omit it here, we could check the kernel pool
%         % to make sure the variable BODY499_RADII
%         % has three elements and numeric data type.
%         % If the variable is not present in the kernel
%         % pool, cspice_bodvrd will signal an error.
%         %
%         body = 'MARS';
%         radii = cspice_bodvrd( body, 'RADII', 3 );
%
%         %
%         %
%         % Calculate the flatness coefficient. Set a bodyfixed
%         % position vector, `x'.
%         %
%         re   = radii(1);
%         rp   = radii(3);
%         flat = ( re - rp ) / re;
%
%         %
%         % Define 1xN arrays of longitudes, latitudes, and altitudes.
%         %
%         lon = [ 0.,   ...
%                 180., ...
%                 180., ...
%                 180., ...
%                 90.,  ...
%                 270., ...
%                 0.,   ...
%                 0.,   ...
%                 0. ];
%
%         lat = [ 0.,  ...
%                 0.,  ...
%                 0.,  ...
%                 0.,  ...
%                 0.,  ...
%                 0.,  ...
%                 90., ...
%                -90., ...
%                 90. ];
%
%         alt = [ 0., ...
%                 0., ...
%                 10., ...
%                 10., ...
%                 0., ...
%                 0., ...
%                 0., ...
%                 0., ...
%                -3376.200 ];
%
%         %
%         % Convert angular measures to radians.
%         %
%         lon = lon*cspice_rpd;
%         lat = lat*cspice_rpd;
%
%         %
%         % Using the same Mars parameters, convert the `lon', `lat', `alt'
%         % vectors to bodyfixed rectangular coordinates.
%         %
%         x = cspice_pgrrec( body, lon, lat, alt, re, flat);
%
%         disp( ['rectan(1)   rectan(2)   rectan(3)' ...
%                '         lon         lat         alt'] )
%         disp( ['---------------------------------' ...
%                '------------------------------------'] )
%
%         %
%         % Create an array of values for output.
%         %
%         output = [  x(1,:);         x(2,:);         x(3,:); ...
%                     lon*cspice_dpr; lat*cspice_dpr; alt ];
%
%         txt = sprintf( '%9.3f   %9.3f   %9.3f   %9.3f   %9.3f   %9.3f\n', ...
%                                                                      output);
%         disp( txt )
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
%      rectan(1)   rectan(2)   rectan(3)         lon         lat         alt
%      ---------------------------------------------------------------------
%       3396.190      -0.000       0.000       0.000       0.000       0.000
%      -3396.190      -0.000       0.000     180.000       0.000       0.000
%      -3406.190      -0.000       0.000     180.000       0.000      10.000
%      -3406.190      -0.000       0.000     180.000       0.000      10.000
%          0.000   -3396.190       0.000      90.000       0.000       0.000
%         -0.000    3396.190       0.000     270.000       0.000       0.000
%          0.000      -0.000    3376.200       0.000      90.000       0.000
%          0.000      -0.000   -3376.200       0.000     -90.000       0.000
%          0.000       0.000       0.000       0.000      90.000   -3376.200
%
%
%   3) Below we show the analogous relationships for the earth,
%      using the reference ellipsoid radii
%
%         Equatorial radius:    6378.140
%         Polar radius:         6356.750
%
%      Note the change in longitudes for points on the +/- Y axis
%      for the earth vs the Mars values.
%
%      rectan(1)   rectan(2)   rectan(3)         lon         lat         alt
%      ---------------------------------------------------------------------
%       6378.140       0.000       0.000       0.000       0.000       0.000
%      -6378.140       0.000       0.000     180.000       0.000       0.000
%      -6388.140       0.000       0.000     180.000       0.000      10.000
%      -6368.140       0.000       0.000     180.000       0.000     -10.000
%          0.000   -6378.140       0.000     270.000       0.000       0.000
%          0.000    6378.140       0.000      90.000       0.000       0.000
%          0.000       0.000    6356.750       0.000      90.000       0.000
%          0.000       0.000   -6356.750       0.000     -90.000       0.000
%          0.000       0.000       0.000       0.000      90.000   -6356.750
%
%-Particulars
%
%   Given the planetographic coordinates of a point, this routine
%   returns the body-fixed rectangular coordinates of the point. The
%   body-fixed rectangular frame is that having the X-axis pass
%   through the 0 degree latitude 0 degree longitude direction, the
%   Z-axis pass through the 90 degree latitude direction, and the
%   Y-axis equal to the cross product of the unit Z-axis and X-axis
%   vectors.
%
%   The planetographic definition of latitude is identical to the
%   planetodetic (also called "geodetic" in SPICE documentation)
%   definition. In the planetographic coordinate system, latitude is
%   defined using a reference spheroid. The spheroid is
%   characterized by an equatorial radius and a polar radius. For a
%   point P on the spheroid, latitude is defined as the angle between
%   the X-Y plane and the outward surface normal at P. For a point P
%   off the spheroid, latitude is defined as the latitude of the
%   nearest point to P on the spheroid. Note if P is an interior
%   point, for example, if P is at the center of the spheroid, there
%   may not be a unique nearest point to P.
%
%   In the planetographic coordinate system, longitude is defined
%   using the spin sense of the body. Longitude is positive to the
%   west if the spin is prograde and positive to the east if the spin
%   is retrograde. The spin sense is given by the sign of the first
%   degree term of the time-dependent polynomial for the body's prime
%   meridian Euler angle "W": the spin is retrograde if this term is
%   negative and prograde otherwise. For the sun, planets, most
%   natural satellites, and selected asteroids, the polynomial
%   expression for W may be found in a SPICE PCK kernel.
%
%   The earth, moon, and sun are exceptions: planetographic longitude
%   is measured positive east for these bodies.
%
%   If you wish to override the default sense of positive longitude
%   for a particular body, you can do so by defining the kernel
%   variable
%
%      BODY<body ID>_PGR_POSITIVE_LON
%
%   where <body ID> represents the NAIF ID code of the body. This
%   variable may be assigned either of the values
%
%      'WEST'
%      'EAST'
%
%   For example, you can have this routine treat the longitude
%   of the earth as increasing to the west using the kernel
%   variable assignment
%
%      BODY399_PGR_POSITIVE_LON = 'WEST'
%
%   Normally such assignments are made by placing them in a text
%   kernel and loading that kernel via cspice_furnsh.
%
%   The definition of this kernel variable controls the behavior of
%   the Mice planetographic routines
%
%      cspice_pgrrec
%      cspice_recpgr
%      cspice_dpgrdr
%      cspice_drdpgr
%
%   It does not affect the other Mice coordinate conversion
%   routines.
%
%-Exceptions
%
%   1)  If the body name `body' cannot be mapped to a NAIF ID code, and
%       if `body' is not a string representation of an integer, the
%       error SPICE(IDCODENOTFOUND) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If the kernel variable
%
%          BODY<ID code>_PGR_POSITIVE_LON
%
%       is present in the kernel pool but has a value other
%       than one of
%
%           'EAST'
%           'WEST'
%
%       the error SPICE(INVALIDOPTION) is signaled by a routine in the
%       call tree of this routine. Case and blanks are ignored when
%       these values are interpreted.
%
%   3)  If polynomial coefficients for the prime meridian of `body' are
%       not available in the kernel pool, and if the kernel variable
%       BODY<ID code>_PGR_POSITIVE_LON is not present in the kernel
%       pool, the error SPICE(MISSINGDATA) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If the equatorial radius is non-positive, the error
%       SPICE(VALUEOUTOFRANGE) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the flattening coefficient is greater than or equal to one,
%       the error SPICE(VALUEOUTOFRANGE) is signaled by a routine in
%       the call tree of this routine.
%
%   6)  If any of the input arguments, `body', `lon', `lat', `alt',
%       `re' or `f', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   7)  If any of the input arguments, `body', `lon', `lat', `alt',
%       `re' or `f', is not of the expected type, or it does not have
%       the expected dimensions and size, an error is signaled by the
%       Mice interface.
%
%   8)  If the input vectorizable arguments `lon', `lat' and `alt' do
%       not have the same measure of vectorization (N), an error is
%       signaled by the Mice interface.
%
%-Files
%
%   This routine expects a kernel variable giving BODY's prime
%   meridian angle as a function of time to be available in the
%   kernel pool. Normally this item is provided by loading a PCK
%   file. The required kernel variable is named
%
%      BODY<body ID>_PM
%
%   where <body ID> represents a string containing the NAIF integer
%   ID code for `body'. For example, if `body' is 'JUPITER', then
%   the name of the kernel variable containing the prime meridian
%   angle coefficients is
%
%      BODY599_PM
%
%   See the PCK Required Reading for details concerning the prime
%   meridian kernel variable.
%
%   The optional kernel variable
%
%      BODY<body ID>_PGR_POSITIVE_LON
%
%   also is normally defined via loading a text kernel. When this
%   variable is present in the kernel pool, the prime meridian
%   coefficients for `body' are not required by this routine. See the
%   -Particulars section below for details.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   KERNEL.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections. Fixed
%       typos in header. Completed list of Mice planetographic routines in
%       -Particulars section. Extended Required Reading document list.
%
%       Edited the header to comply with NAIF standard.
%       Split the existing code example into two separate examples and
%       added example 3.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 12-MAR-2012 (EDW) (SCK)
%
%       Corrected misspellings.
%
%   -Mice Version 1.0.0, 22-JAN-2008 (EDW)
%
%-Index_Entries
%
%   convert planetographic to rectangular coordinates
%
%-&

function [rectan] = cspice_pgrrec( body, lon, lat, alt, re, f)

   switch nargin
      case 6

         body = zzmice_str(body);
         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
         alt  = zzmice_dp(alt);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);

      otherwise

         error ( ['Usage: [_rectan(3)_] = ' ...
                  'cspice_pgrrec( body, _lon_, _lat_, _alt_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice( 'pgrrec_c', body, lon, lat, alt, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end
