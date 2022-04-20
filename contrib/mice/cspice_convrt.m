%-Abstract
%
%   CSPICE_CONVRT performs a conversion from a measurement in
%   one unit set to the corresponding measure in another unit
%   set.
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
%      x        value(s) representing a measurement in the units specified
%               by `in'.
%
%               [1,n] = size(x); double = class(x)
%
%      in       the string specifying the units associated with measurement
%               `x'.
%
%               [1,c1] = size(in); char = class(in)
%
%                  or
%
%               [1,1] = size(in); cell = class(in)
%
%      out      the string specifying the units desired for the measurement
%               `x'.
%
%               [1,c2] = size(out); char = class(out)
%
%                  or
%
%               [1,1] = size(out); cell = class(out)
%
%               Acceptable units for `in' and `out':
%
%                  Angles:                 'RADIANS'
%                                          'DEGREES'
%                                          'ARCMINUTES'
%                                          'ARCSECONDS'
%                                          'HOURANGLE'
%                                          'MINUTEANGLE'
%                                          'SECONDANGLE'
%
%                  Metric Distances:       'METERS'
%                                          'M'
%                                          'KILOMETERS'
%                                          'KM'
%                                          'CENTIMETERS'
%                                          'CM'
%                                          'MILLIMETERS'
%                                          'MM'
%
%                  English Distances:      'FEET'
%                                          'INCHES'
%                                          'YARDS'
%                                          'STATUTE_MILES'
%                                          'NAUTICAL_MILES'
%
%                  Astrometric Distances:  'AU'
%                                          'PARSECS'
%                                          'LIGHTSECS'
%                                          'LIGHTYEARS' julian lightyears
%
%                  Time:                   'SECONDS'
%                                          'MINUTES'
%                                          'HOURS'
%                                          'DAYS'
%                                          'JULIAN_YEARS'
%                                          'TROPICAL_YEARS'
%                                          'YEARS' (same as julian years)
%
%               Neither `in' nor `out' are case sensitive.
%
%   the call:
%
%      [y] = cspice_convrt( x, in, out )
%
%   returns:
%
%      y        the value(s) representing the input `x' measurement converted
%               to the units defined by `out'.
%
%               [1,n] = size(y); double = class(y)
%
%               `y' returns with the same vectorization measure, N, as `x'.
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
%   1) Convert 300 miles (statute miles) to kilometers and determine
%      the number of lightyears in a parsec.
%
%      Example code begins here.
%
%
%      function convrt_ex1()
%
%         %
%         % Convert 300 miles (statute miles) to kilometers.
%         %
%         dist_sm = 300;
%         dist_km = cspice_convrt( dist_sm, 'statute_miles', 'km' );
%
%         fprintf( '300 miles in km  : %15.6f\n', dist_km )
%
%         %
%         % Determine the number of lightyears in a vector of parsec values.
%         %
%         parsec     = 1;
%         lightyears = cspice_convrt( parsec, 'parsecs', 'lightyears' );
%
%         fprintf( 'Lightyears/parsec: %15.6f\n', lightyears )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      300 miles in km  :      482.803200
%      Lightyears/parsec:        3.261564
%
%
%   2) Determine the number of lightyears for different parsec values
%      using the vectorized capability of the cspice_convert call.
%
%      Example code begins here.
%
%
%      function convrt_ex2()
%
%         %
%         % Determine the number of lightyears in a vector of parsec
%         % values.
%         %
%         parsec     = [1, 3, 5];
%         lightyears = cspice_convrt( parsec, 'parsecs', 'lightyears' );
%
%         fprintf( '# parsec  lightyears\n' );
%         fprintf( '--------  ----------\n' );
%
%         for i=1:3
%
%            fprintf( '%8.3f  %10.6f\n', parsec(i), lightyears(i) )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      # parsec  lightyears
%      --------  ----------
%         1.000    3.261564
%         3.000    9.784691
%         5.000   16.307819
%
%
%-Particulars
%
%   This routine converts a measurement x given in units specified by
%   in to the equivalent value y in units specified by out.
%
%   If a unit is not recognized, an error message is produced that
%   indicates which one was not recognized.
%
%   If input and output units are incompatible (for example angle
%   and distance units) and error message will be produced stating
%   the requested units and associated types.
%
%-Exceptions
%
%   1)  If the input units, output units, or both input and output
%       units are not recognized, the error SPICE(UNITSNOTREC) is
%       signaled by a routine in the call tree of this routine.
%
%   2)  If the units being converted between are incompatible, the
%       error SPICE(INCOMPATIBLEUNITS) is signaled by a routine in the
%       call tree of this routine.
%
%   3)  If any of the input arguments, `x', `in' or `out', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `x', `in' or `out', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  This routine does not do any checking for overflow. The caller
%       is required to make sure that the units used for the
%       measurement are such that no floating point overflow will
%       occur when the conversion is performed.
%
%   2)  Some of the units are not "defined" quantities. In such a case
%       a best estimate is provided as of the date of the current
%       version of this routine. Those estimated quantities are:
%
%          AU               The astronomical unit. The value was taken
%                           from the JPL ephemeris DE125. This value
%                           is an approximation and should not be used
%                           for high-accuracy work. It agrees with the
%                           value used in the JPL planetary ephemeris
%                           DE430 (149597870.700 km) at the 100m
%                           level.
%
%          TROPICAL_YEARS   The tropical year is the time from equinox
%                           to equinox. This varies slightly with
%                           time.
%
%          PARSECS          The parsec is the distance to an object
%                           whose parallax angle is one arcsecond. Its
%                           value is dependent upon the value of the
%                           astronomical unit.
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
%   -Mice Version 1.1.0, 23-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Added examples' problem statements and reformatted examples'
%       output.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 05-APR-2017 (EDW)
%
%       Header update to correspond to current SPICELIB/CSPICE version.
%
%   -Mice Version 1.0.2, 30-OCT-2014 (EDW)
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
%   convert units
%
%-&

function [y] = cspice_convrt( x, in, out)

   switch nargin
      case 3

         x   = zzmice_dp(x);
         in  = zzmice_str(in);
         out = zzmice_str(out);

      otherwise

         error( 'Usage: [_y_] = cspice_convrt( _x_, `in`, `out`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [y] = mice('convrt_c', x, in, out);
   catch spiceerr
      rethrow(spiceerr)
   end
