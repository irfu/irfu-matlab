%-Abstract
%
%   CSPICE_RECRAD converts rectangular (Cartesian) coordinates to
%   right ascension, declination coordinates.
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
%               position(s).
%
%               [3,n] = size(rectan); double = class(rectan)
%
%   the call:
%
%      [range, ra, dec] = cspice_recrad(rectan)
%
%   returns:
%
%      range    the value(s) describing the distance of the position
%               from the origin.
%
%               [1,n] = size(range); double = class(range)
%
%      ra       the value(s) describing the right ascension of the position
%               as measured in radians.
%
%               [1,n] = size(ra); double = class(ra)
%
%      dec      the value(s) describing the declination of the position as
%               measured in radians.
%
%               [1,n] = size(dec); double = class(dec)
%
%               'range', 'ra', and 'dec' return with the same
%               vectorization measure, N, as 'rectan'.
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
%   1) Output the right ascension and declination of the earth's pole
%      in the J2000 frame approximately every six months for the time
%      interval January 1, 2000 to January 1, 2005 (UTC).
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: recrad_ex1.tm
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
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            naif0012.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'pck00010.tpc',
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
%      function recrad_ex1()
%
%         %
%         % Load a standard kernel set.
%         %
%         cspice_furnsh( 'recrad_ex1.tm' )
%
%         %
%         % Define the time bounds for the time interval,
%         % 5 years,  convert to ephemeris time J2000.
%         %
%         utc_bounds = [ '1 Jan 2000'; '1 Jan 2005' ];
%         et_bounds = cspice_str2et( utc_bounds);
%
%         %
%         % Step in units of 6 months. 5 years ~ 10 steps.
%         %
%         step = (et_bounds(2) - et_bounds(1)) / 10.;
%
%         %
%         % Create an array of 10 ephemeris times starting at
%         % et_bounds(1) in intervals of 'step'.
%         %
%         et = [0:9]*step + et_bounds(1);
%
%         %
%         % Set the conversion constant "radians to degrees."
%         %
%         r2d = cspice_dpr;
%
%         %
%         % Convert the 10-vector of 'et' to an array of corresponding
%         % transformation matrices (dimensions (3,3,10) ).
%         %
%         mat = cspice_pxform( 'IAU_EARTH', 'J2000', et);
%
%         %
%         % Extract the pole vector from the transformation matrix,
%         % convert to RA and DEC expressed in degrees.
%         %
%         % The last column in each matrix is the pole vector (z = (0,0,1))
%         % of the earth in IAU expressed in J2000. We need to copy the
%         % set of pole vectors to a 3xN array. Use reshape to do this.
%         %
%         pole = reshape( mat(:,3,:), 3,[] );
%
%         [radius, ra, dec] = cspice_recrad(pole);
%
%         ra  = ra * r2d;
%         dec = dec * r2d;
%
%         %
%         % Convert ephemeris times to UTC strings.
%         %
%         utcstr = cspice_et2utc( et, 'C', 0 );
%
%         disp( '      UTC time        Right Ascension    Declination')
%         disp( '--------------------  ---------------  ---------------')
%         for i=1:10
%            fprintf( '%s  %15.9f  %15.9f\n' , utcstr(i,:), ra(i), dec(i))
%         end
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
%            UTC time        Right Ascension    Declination
%      --------------------  ---------------  ---------------
%      2000 JAN 01 00:00:00    180.000008762     89.999992386
%      2000 JUL 01 16:48:00    359.996802446     89.997221470
%      2000 DEC 31 09:36:00    359.993596129     89.994435326
%      2001 JUL 02 02:24:00    359.990389813     89.991649182
%      2001 DEC 31 19:12:00    359.987183497     89.988863039
%      2002 JUL 02 12:00:00    359.983977181     89.986076895
%      2003 JAN 01 04:48:00    359.980770864     89.983290751
%      2003 JUL 02 21:36:00    359.977564548     89.980504607
%      2004 JAN 01 14:24:00    359.974358232     89.977718464
%      2004 JUL 02 07:12:00    359.971151916     89.974932320
%
%
%-Particulars
%
%   This routine returns the range, right ascension, and declination
%   of a point specified in rectangular coordinates.
%
%   The output is defined by a distance from a central reference
%   point, an angle from a reference meridian, and an angle above
%   the equator of a sphere centered at the central reference
%   point.
%
%-Exceptions
%
%   1)  If the X and Y components of `rectan' are both zero, the
%       right ascension is set to zero.
%
%   2)  If `rectan' is the zero vector, right ascension and declination
%       are both set to zero.
%
%   3)  If the input argument `rectan' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `rectan' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
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
%   -Mice Version 1.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and meta-kernel. Reformatted code
%       example output.
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
%   rectangular coordinates to ra and dec
%   rectangular to right_ascension and declination
%
%-&

function [range, ra, dec] = cspice_recrad(rectan)

   switch nargin
      case 1

         rectan = zzmice_dp(rectan);

      otherwise
         error ( 'Usage: [_range_, _ra_, _dec_] = cspice_recrad(_rectan(3)_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [range, ra, dec] = mice('recrad_c',rectan);
   catch spiceerr
      rethrow(spiceerr)
   end
