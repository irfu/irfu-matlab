%-Abstract
%
%   CSPICE_ET2LST computes the local solar time for a given ephemeris epoch
%   `et' for an object on the surface of a body at a specified longitude.
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
%      et       the ephemeris time(s) expressed as ephemeris seconds past
%               J2000 at which a local time is desired.
%
%               [1,n] = size(et); double = class(et)
%
%      body     the SPICE ID-code of the body on which to measure local time.
%
%               [1,1] = size(body); int32 = class(body)
%
%      lon      longitude (either planetocentric or planetographic)
%               in radians of the site on the surface
%               of body for which local time should be computed.
%
%               [1,1] = size(lon); double = class(lon)
%
%      type     the name for the form of longitude supplied by the
%               variable `lon'.
%
%               [1,c1] = size(type); char = class(type)
%
%                  or
%
%               [1,1] = size(type); cell = class(type)
%
%               Allowed values are 'PLANETOCENTRIC' and 'PLANETOGRAPHIC'. Note
%               the case of the letters in type is insignificant. Both
%               'PLANETOCENTRIC' and 'planetocentric' are recognized. Leading
%               and trailing blanks in type are not significant.
%
%   the call:
%
%      [hr, mn, sc, time, ampm] = cspice_et2lst( et, body, lon, type )
%
%   returns:
%
%      hr       the value(s) describing the integral number of the local
%               "hour" of the site specified at epoch `et'.
%
%               [1,n] = size(et); double = class(et)
%
%               Note that an "hour" of local time does not have the same
%               duration as an hour measured by conventional clocks. It is
%               simply a representation of an angle.
%
%      mn       the value(s) describing the integral number of "minutes" past
%               the hour of the local time of the site at the epoch `et'.
%
%               [1,n] = size(et); double = class(et)
%
%               Again note that a "local minute" is not the same as a minute
%               you would measure with conventional clocks.
%
%      sc       the value(s) describing the integral number of "seconds" past
%               the minute of the local time of the site at the epoch `et'.
%
%               [1,n] = size(et); double = class(et)
%
%               Again note that a "local second" is not the same as a second
%               you would measure with conventional clocks.
%
%      time     the array of local time(s) on a "24 hour" local clock.
%
%               [n,c2] = size(time); char = class(time)
%
%      ampm     array of local time(s) on a "12 hour" local clock together
%               with the traditional AM/PM label to indicate whether the Sun
%               has crossed the local zenith meridian.
%
%               [n,c3] = size(ampm); char = class(ampm)
%
%               All output arguments return with the same measure of
%               vectorization, N, as `et'.
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
%   1) Compute the local time at particular location on Mars.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: et2lst_ex1.tm
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
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            naif0012.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc',
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
%      function et2lst_ex1()
%
%         %
%         % Load a leapseconds kernel.
%         %
%         cspice_furnsh( 'et2lst_ex1.tm' )
%
%         %
%         % Define two UTC time strings to `utc'
%         %
%         utc                        = strvcat( '2002 SEP 02 00:00:00',    ...
%                                               '2002 SEP 30 00:00:00' );
%
%         %
%         % Convert `utc' the ephemeris time, 'et'
%         %
%         et                          = cspice_str2et(utc);
%
%         %
%         % Define a planetographic longitude in degrees, convert the
%         % value to radians
%         %
%         dlon                       =  326.17;
%         rlon                       =  dlon * cspice_rpd;
%
%         %
%         % Convert inputs to Local Solar Time.
%         %
%         [hr, mn, sc, time, ampm] = cspice_et2lst( et,                    ...
%                                                   499,                   ...
%                                                   rlon,                  ...
%                                                   'PLANETOGRAPHIC');
%
%         fprintf( ['The local time at Mars %6.2f degrees '                ...
%                  'planetographic longitude:\n'],                         ...
%                  dlon )
%         fprintf( '   at UTC %s, LST = %s\n', utc(1,:), ampm(1,:) )
%         fprintf( '   at UTC %s, LST = %s\n', utc(2,:), ampm(2,:) )
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
%      The local time at Mars 326.17 degrees planetographic longitude:
%         at UTC 2002 SEP 02 00:00:00, LST = 03:25:35 A.M.
%         at UTC 2002 SEP 30 00:00:00, LST = 09:33:00 A.M.
%
%
%-Particulars
%
%   This routine returns the local solar time at a user
%   specified location on a user specified body.
%
%   Let `sunlng' be the planetocentric longitude (in degrees) of
%   the sun as viewed from the center of the body of interest.
%
%   Let `sitlng' be the planetocentric longitude (in degrees) of
%   the site for which local time is desired.
%
%   We define local time to be 12 + (sitlng - sunlng)/15
%
%   (where appropriate care is taken to map ( sitlng - sunlng )
%   into the range from -180 to 180).
%
%   Using this definition, we see that from the point of view
%   of this routine, local solar time is simply a measure of angles
%   between meridians on the surface of a body. Consequently,
%   this routine is not appropriate for computing "local times"
%   in the sense of Pacific Standard Time. For computing times
%   relative to standard time zones on earth, see the routines
%   cspice_timout and cspice_str2et.
%
%
%   Regarding planetographic longitude
%   ----------------------------------
%
%   In the planetographic coordinate system, longitude is defined
%   using the spin sense of the body. Longitude is positive to the
%   west if the spin is prograde and positive to the east if the spin
%   is retrograde. The spin sense is given by the sign of the first
%   degree term of the time-dependent polynomial for the body's prime
%   meridian Euler angle "W":  the spin is retrograde if this term is
%   negative and prograde otherwise. For the sun, planets, most
%   natural satellites, and selected asteroids, the polynomial
%   expression for W may be found in a SPICE PCK kernel.
%
%   The earth, moon, and sun are exceptions: planetographic longitude
%   is measured positive east for these bodies.
%
%   If you wish to override the default sense of positive
%   planetographic longitude for a particular body, you can do so by
%   defining the kernel variable
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
%-Exceptions
%
%   1)  This routine defines local solar time for any point on the
%       surface of the Sun to be 12:00:00 noon.
%
%   2)  If the `type' of the coordinates is not recognized, the error
%       SPICE(UNKNOWNSYSTEM) is signaled by a routine in the call tree
%       of this routine.
%
%   3)  If the body-fixed frame to associate with `body' cannot be
%       determined, the error SPICE(CANTFINDFRAME) is signaled by a
%       routine in the call tree of this routine.
%
%   4)  If insufficient data is available to compute the location of
%       the sun in body-fixed coordinates, an error is signaled by a
%       routine in the call tree of this routine.
%
%   5)  If the BODY#_PM keyword required to determine the body
%       rotation sense is not found in the POOL or if it is found but
%       is not a numeric keyword with at least two elements, the error
%       SPICE(CANTGETROTATIONTYPE) is signaled by a routine in the
%       call tree of this routine.
%
%   6)  If any of the input arguments, `et', `body', `lon' or `type',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   7)  If any of the input arguments, `et', `body', `lon' or `type',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   Suitable SPK and PCK files must be loaded prior to calling this
%   routine so that the body-fixed position of the sun relative to
%   `body' can be computed. The PCK files must contain the standard
%   BODY#_PM keyword need by this routine to determine the body
%   rotation sense.
%
%   When the input longitude is planetographic, the default
%   interpretation of this value can be overridden using the optional
%   kernel variable
%
%      BODY<body ID>_PGR_POSITIVE_LON
%
%   which is normally defined via loading a text kernel.
%
%-Restrictions
%
%   1)  This routine relies on being able to determine the name
%       of the body-fixed frame associated with `body' through the
%       frames subsystem. If the `body' specified is NOT one of the
%       nine planets or their satellites, you will need to load
%       an appropriate frame definition kernel that contains
%       the relationship between the body id and the body-fixed frame
%       name. See frames.req required reading for more details
%       on specifying this relationship.
%
%   2)  The routine determines the body rotation sense using the PCK
%       keyword BODY#_PM. Therefore, you will need to a text PCK file
%       defining the complete set of the standard PCK body rotation
%       keywords for the body of interest. The text PCK file must be
%       loaded independently of whether a binary PCK file providing
%       rotation data for the same body is loaded or not.
%
%   3)  Although it is not currently the case for any of the Solar
%       System bodies, it is possible that the retrograde rotation
%       rate of a body would be slower than the orbital rate of the
%       body rotation around the Sun. The routine does not account for
%       such cases; for them it will compute incorrect the local time
%       progressing backwards.
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
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Changed output argument names "min" and "sec" to "mn" and "sc".
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and meta-kernel.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.2, 05-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 07-MAR-2007 (EDW)
%
%-Index_Entries
%
%   Compute the local time for a point on a body.
%
%-&

function [hr, mn, sc, time, ampm] = cspice_et2lst( et, body, lon, type)

   switch nargin
      case 4

         et    = zzmice_dp(et);
         body  = zzmice_int(body);
         lon   = zzmice_dp(lon);
         type  = zzmice_str(type);

      otherwise

         error ( ['Usage: [ _hr_, _mn_, _sc_, _`time`_, _`ampm`_] = ' ...
                 'cspice_et2lst( _et_, body, lon, `type`)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [hr, mn, sc, time, ampm] = mice('et2lst_c', et, body, lon, type);

      %
      % Convert the integers returned from the interface to double precision
      % in case a user includes the return arguments in a calculation
      % with other doubles.
      %
      hr = zzmice_dp(hr);
      mn = zzmice_dp(mn);
      sc = zzmice_dp(sc);

   catch spiceerr
      rethrow(spiceerr)
   end





