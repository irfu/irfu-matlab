%-Abstract
%
%   CSPICE_DRDPGR computes the Jacobian matrix of the transformation
%   from planetographic to rectangular coordinates.
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
%               `body' is used by this routine to look up from the kernel
%               pool the prime meridian rate coefficient giving the body's
%               spin sense.
%
%      lon      the planetographic longitude(s) of the input point(s).
%
%               [1,n] = size(lon); double = class(lon)
%
%               This is the angle between the prime meridian and the
%               meridian containing the input point.For bodies having
%               prograde (aka direct) rotation, the direction of increasing
%               longitude is positive west:from the +X axis of the
%               rectangular coordinate system toward the -Y axis.For bodies
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
%               Longitude is measured in radians. On input, the range
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
%               Latitude is measured in radians. On input, the
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
%               This spheroid is a volume of revolution:its horizontal
%               cross sections are circular.The shape of the spheroid is
%               defined by an equatorial radius `re' and a polar radius
%               `rp'. Units of `re' must match those of `alt'.
%
%      f        the flattening coefficient
%
%                  f = (re-rp) / re
%
%               where `rp' is the polar radius of the spheroid.
%
%               [1,1] = size(f); double = class(f)
%
%               The units of `rp' match those of `re'. (More importantly
%               rp = re*(1-f) )
%
%   the call:
%
%      [jacobi] = cspice_drdpgr( body, lon, lat, alt, re, f )
%
%   returns:
%
%      jacobi   the matrix(es) of partial derivatives of the conversion from
%               planetographic to rectangular coordinates.
%
%               If [1,1] = size(lon) then [3,3]   = size(jacobi)
%               If [1,n] = size(lon) then [3,3,n] = size(jacobi)
%                                          double = class(jacobi)
%
%               It has the form
%
%                  .-                         -.
%                  | dx/dlon  dx/dlat  dx/dalt |
%                  | dy/dlon  dy/dlat  dy/dalt |
%                  | dz/dlon  dz/dlat  dz/dalt |
%                  `-                         -'
%
%               evaluated at the input values of `lon', `lat' and `alt'.
%
%               `jacobi' returns with the same vectorization measure (N)
%               as `lon', `lat' and `alt'.
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
%   1) Find the planetographic state of the Earth as seen from
%      Mars in the J2000 reference frame at January 1, 2005 TDB.
%      Map this state back to rectangular coordinates as a check.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: drdpgr_ex1.tm
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
%            de405.bsp                     Planetary ephemeris
%            pck00008.tpc                  Planet orientation and
%                                          radii
%            naif0007.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de405.bsp',
%                                'pck00008.tpc',
%                                'naif0007.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function drdpgr_ex1()
%
%         %
%         % Load SPK, PCK, and LSK kernels, use a meta kernel for
%         % convenience.
%         %
%         cspice_furnsh( 'drdpgr_ex1.tm' );
%
%         %
%         % Look up the radii for Mars.  Although we
%         % omit it here, we could first call badkpv_c
%         % to make sure the variable BODY499_RADII
%         % has three elements and numeric data type.
%         % If the variable is not present in the kernel
%         % pool, bodvrd_c will signal an error.
%         %
%         [radii] = cspice_bodvrd( 'MARS', 'RADII', 3 );
%
%         %
%         % Compute flattening coefficient.
%         %
%         re  =  radii(1);
%         rp  =  radii(3);
%         f   =  ( re - rp ) / re;
%
%         %
%         % Look up the geometric state of earth as seen from Mars at
%         % January 1, 2005 TDB, relative to the J2000 reference
%         % frame.
%         %
%         [et] = cspice_str2et( 'January 1, 2005 TDB' );
%
%         [state, lt] = cspice_spkezr( 'Earth', et, 'J2000', ...
%                                        'LT+S', 'Mars'        );
%
%         %
%         % Convert position to planetographic coordinates.
%         %
%         [lon, lat, alt] = cspice_recpgr( 'mars', state(1:3), re, f );
%
%         %
%         % Convert velocity to planetographic coordinates.
%         %
%
%         [jacobi] = cspice_dpgrdr( 'MARS', state(1), state(2), ...
%                                           state(3), re,       f );
%
%         pgrvel = jacobi * state(4:6);
%
%         %
%         % As a check, convert the planetographic state back to
%         % rectangular coordinates.
%         %
%         [rectan] = cspice_pgrrec( 'mars', lon, lat, alt, re, f );
%         [jacobi] = cspice_drdpgr( 'mars', lon, lat, alt, re, f );
%
%         drectn = jacobi * pgrvel;
%
%         fprintf( '\n' )
%         fprintf( 'Rectangular coordinates:\n' )
%         fprintf( '\n' )
%         fprintf( '  X (km)                 =  %17.8e\n', state (1) )
%         fprintf( '  Y (km)                 =  %17.8e\n', state (2) )
%         fprintf( '  Z (km)                 =  %17.8e\n', state (3) )
%         fprintf( '\n' )
%         fprintf( 'Rectangular velocity:\n' )
%         fprintf( '\n' )
%         fprintf( '  dX/dt (km/s)           =  %17.8e\n', state (4) )
%         fprintf( '  dY/dt (km/s)           =  %17.8e\n', state (5) )
%         fprintf( '  dZ/dt (km/s)           =  %17.8e\n', state (6) )
%         fprintf( '\n' )
%         fprintf( 'Ellipsoid shape parameters:\n' )
%         fprintf( '\n' )
%         fprintf( '  Equatorial radius (km) =  %17.8e\n', re )
%         fprintf( '  Polar radius      (km) =  %17.8e\n', rp )
%         fprintf( '  Flattening coefficient =  %17.8e\n', f )
%         fprintf( '\n' )
%         fprintf( 'Planetographic coordinates:\n' )
%         fprintf( '\n' )
%         fprintf( '  Longitude (deg)        =  %17.8e\n', ...
%                                                   lon / cspice_rpd )
%         fprintf( '  Latitude  (deg)        =  %17.8e\n', ...
%                                                   lat / cspice_rpd )
%         fprintf( '  Altitude  (km)         =  %17.8e\n', alt )
%         fprintf( '\n' )
%         fprintf( 'Planetographic velocity:\n' )
%         fprintf( '\n' )
%         fprintf( '  d Longitude/dt (deg/s) =  %17.8e\n', ...
%                                               pgrvel(1)/cspice_rpd )
%         fprintf( '  d Latitude/dt  (deg/s) =  %17.8e\n', ...
%                                               pgrvel(2)/cspice_rpd )
%         fprintf( '  d Altitude/dt  (km/s)  =  %17.8e\n', pgrvel(3) )
%         fprintf( '\n' )
%         fprintf( 'Rectangular coordinates from inverse mapping:\n' )
%         fprintf( '\n' )
%         fprintf( '  X (km)                 =  %17.8e\n', rectan (1) )
%         fprintf( '  Y (km)                 =  %17.8e\n', rectan (2) )
%         fprintf( '  Z (km)                 =  %17.8e\n', rectan (3) )
%         fprintf( '\n' )
%         fprintf( 'Rectangular velocity from inverse mapping:\n' )
%         fprintf( '\n' )
%         fprintf( '  dX/dt (km/s)           =  %17.8e\n', drectn (1) )
%         fprintf( '  dY/dt (km/s)           =  %17.8e\n', drectn (2) )
%         fprintf( '  dZ/dt (km/s)           =  %17.8e\n', drectn (3) )
%         fprintf( '\n' )
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
%      Rectangular coordinates:
%
%        X (km)                 =     1.46039732e+08
%        Y (km)                 =     2.78546607e+08
%        Z (km)                 =     1.19750315e+08
%
%      Rectangular velocity:
%
%        dX/dt (km/s)           =    -4.70432720e+01
%        dY/dt (km/s)           =     9.07326134e+00
%        dZ/dt (km/s)           =     4.75791694e+00
%
%      Ellipsoid shape parameters:
%
%        Equatorial radius (km) =     3.39619000e+03
%        Polar radius      (km) =     3.37620000e+03
%        Flattening coefficient =     5.88600756e-03
%
%      Planetographic coordinates:
%
%        Longitude (deg)        =     2.97667659e+02
%        Latitude  (deg)        =     2.08445040e+01
%        Altitude  (km)         =     3.36531825e+08
%
%      Planetographic velocity:
%
%        d Longitude/dt (deg/s) =    -8.35770664e-06
%        d Latitude/dt  (deg/s) =     1.59355667e-06
%        d Altitude/dt  (km/s)  =    -1.12116008e+01
%
%      Rectangular coordinates from inverse mapping:
%
%        X (km)                 =     1.46039732e+08
%        Y (km)                 =     2.78546607e+08
%        Z (km)                 =     1.19750315e+08
%
%      Rectangular velocity from inverse mapping:
%
%        dX/dt (km/s)           =    -4.70432720e+01
%        dY/dt (km/s)           =     9.07326134e+00
%        dZ/dt (km/s)           =     4.75791694e+00
%
%
%-Particulars
%
%   It is often convenient to describe the motion of an object in the
%   planetographic coordinate system. However, when performing
%   vector computations it's hard to beat rectangular coordinates.
%
%   To transform states given with respect to planetographic
%   coordinates to states with respect to rectangular coordinates,
%   one makes use of the Jacobian of the transformation between the
%   two systems.
%
%   Given a state in planetographic coordinates
%
%      ( lon, lat, alt, dlon, dlat, dalt )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation:
%
%                  t          |                                  t
%      (dx, dy, dz)   = jacobi|              * (dlon, dlat, dalt)
%                             |(lon,lat,alt)
%
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(lon,lat,alt)
%
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
%   the CSPICE planetographic routines
%
%      cspice_pgrrec
%      cspice_recpgr
%      cspice_dpgrdr
%      cspice_drdpgr
%
%   It does not affect the other SPICE coordinate conversion
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
%   -Particulars section for details.
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
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Added complete example to the -Examples section.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 09-NOV-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   Jacobian of rectangular w.r.t. planetographic coordinates
%
%-&

function [jacobi] = cspice_drdpgr( body, lon, lat, alt, re, f)

   switch nargin
      case 6

         body = zzmice_str(body);
         lon  = zzmice_dp(lon);
         lat  = zzmice_dp(lat);
         alt  = zzmice_dp(alt);
         re   = zzmice_dp(re);
         f    = zzmice_dp(f);

      otherwise

         error( ['Usage: [_jacobi(3,3)_] = '...
                 'cspice_drdpgr( `body`, _lon_, _lat_, _alt_, re, f)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdpgr_c', body, lon, lat, alt, re, f);
   catch spiceerr
      rethrow(spiceerr)
   end
