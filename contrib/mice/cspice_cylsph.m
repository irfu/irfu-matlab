%-Abstract
%
%   CSPICE_CYLSPH converts cylindrical coordinates to spherical
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
%      r        the value(s) describing the distance of the point of
%               interest from z axis.
%
%               [1,n] = size(r); double = class(r)
%
%      clon     the value(s) describing the cylindrical angle of the point of
%               interest from the XZ plane measured in radians.
%
%               [1,n] = size(clon); double = class(clon)
%
%      z        the value(s) describing the height of the point above
%               the XY plane.
%
%               [1,n] = size(z); double = class(z)
%
%   the call:
%
%      [radius, colat, slon] = cspice_cylsph( r, clon, z)
%
%   returns:
%
%      radius   the value(s) describing the distance of the position
%               from the origin.
%
%               [1,n] = size(radius); double = class(radius)
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
%               The argument `radius' returns in the same units associated
%               with `r' and `z'.
%
%               `radius', `colat', and `slon' return with the vectorization
%               measure (N) as `r', `clon', and `z'.
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
%   1) Compute the cylindrical coordinates of the position of the Moon
%      as seen from the Earth, and convert them to spherical and
%      rectangular coordinates.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: cylsph_ex1.tm
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
%      function cylsph_ex1()
%
%         %
%         % Load an SPK and leapseconds kernels.
%         %
%         cspice_furnsh( 'cylsph_ex1.tm' )
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
%         % Convert the position vector 'pos' to cylindrical
%         % coordinates.
%         %
%         [r, clon, z]           = cspice_reccyl(pos);
%         fprintf( '\n' )
%         fprintf( 'Cylindrical coordinates:\n' )
%         fprintf( '   Radius       (km): %20.8f\n', r )
%         fprintf( '   Longitude   (deg): %20.8f\n', clon * cspice_dpr )
%         fprintf( '   Z            (km): %20.8f\n', z )
%
%         %
%         % Convert the cylindrical coords to spherical.
%         %
%         [radius, colat, slon] = cspice_cylsph(r, clon, z);
%         fprintf( '\n' )
%         fprintf( 'Spherical coordinates:\n' )
%         fprintf( '   Radius       (km): %20.8f\n', radius )
%         fprintf( '   Polar Angle (deg): %20.8f\n', colat * cspice_dpr )
%         fprintf( '   Longitude   (deg): %20.8f\n', slon  * cspice_dpr )
%
%         %
%         % Convert the spherical to rectangular.
%         %
%         [rectan]              = cspice_sphrec(radius, colat, slon);
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
%      Cylindrical coordinates:
%         Radius       (km):      383289.01777726
%         Longitude   (deg):         261.65040211
%         Z            (km):     -126505.93063865
%
%      Spherical coordinates:
%         Radius       (km):      403626.33912495
%         Polar Angle (deg):         108.26566077
%         Longitude   (deg):         261.65040211
%
%      Rectangular coordinates from cspice_sphrec:
%         X            (km):      -55658.44323296
%         Y            (km):     -379226.32931475
%         Z            (km):     -126505.93063865
%
%
%   2) Create a table showing a variety of cylindrical coordinates
%      and the corresponding spherical coordinates.
%
%      Corresponding cylindrical and spherical coordinates are
%      listed to three decimal places. All input and output angles
%      are in degrees.
%
%
%      Example code begins here.
%
%
%      function cylsph_ex2()
%
%         %
%         % Define six sets of cylindrical coordinates, `clon' expressed
%         % in degrees - converted to radians by use of cspice_rpd.
%         %
%         r     = [ 1.,  1.,   1.,   1.,   0.,  0. ];
%         clon  = [ 0., 90., 180., 180., 180., 33. ] * cspice_rpd;
%         z     = [ 0.,  0.,   1.,  -1.,   1.,  0. ];
%
%         %
%         % ...convert the cylindrical coordinates to spherical coordinates
%         %
%         [radius, colat, slon] = cspice_cylsph(r, clon, z);
%
%         %
%         % ...convert angular measure to degrees.
%         %
%         colat = colat * cspice_dpr;
%         slon  = slon    * cspice_dpr;
%         clon = clon   * cspice_dpr;
%
%         %
%         % Output banner.
%         %
%         disp('     r       clon        z       radius    slon      colat ')
%         disp('  -------  --------  --------   -------  --------  --------')
%
%         %
%         % Create an array of values for output.
%         %
%         output = [ r; clon; z; radius; slon; colat ];
%
%         txt = sprintf( '%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', output );
%         disp( txt )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%           r       clon        z       radius    slon      colat
%        -------  --------  --------   -------  --------  --------
%          1.000     0.000     0.000     1.000     0.000    90.000
%          1.000    90.000     0.000     1.000    90.000    90.000
%          1.000   180.000     1.000     1.414   180.000    45.000
%          1.000   180.000    -1.000     1.414   180.000   135.000
%          0.000   180.000     1.000     1.000   180.000     0.000
%          0.000    33.000     0.000     0.000    33.000     0.000
%
%
%-Particulars
%
%   This returns the spherical coordinates of a point whose position
%   is input through cylindrical coordinates.
%
%-Exceptions
%
%   1)  If any of the input arguments, `r', `clon' or `z', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `r', `clon' or `z', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%   3)  If the input vectorizable arguments `r', `clon' and `z' do not
%       have the same measure of vectorization (N), an error is
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
%       Changed the argument names "lonc" and "lon" to "clon" and "slon"
%       for consistency with other routines.
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       meta-kernel to example #1. Updated code example #1 to produce
%       formatted output and added a call to cspice_kclear. Added the
%       problem statement to both examples.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 30-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 01-DEC-2005 (EDW)
%
%-Index_Entries
%
%   cylindrical to spherical
%
%-&

function [radius, colat, slon] = cspice_cylsph( r, clon, z )

   switch nargin
      case 3

         r    = zzmice_dp(r);
         clon = zzmice_dp(clon);
         z    = zzmice_dp(z);

      otherwise

         error ( ['Usage: [_radius_, _colat_, _slon_] = '...
                                     'cspice_cylsph(_r_, _clon_, _z_)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [radius, colat, slon] = mice('cylsph_c', r, clon, z);
   catch spiceerr
      rethrow(spiceerr)
   end


