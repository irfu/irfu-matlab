%-Abstract
%
%   CSPICE_SPHCYL converts spherical coordinates to cylindrical
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
%   the call:
%
%      [r, clon, z] = cspice_sphcyl( radius, colat, slon )
%
%   returns:
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
%               The arguments `r' and `z' return in the same units associated
%               with `radius'.
%
%               `r', `clon', and `z' return with the same vectorization
%               measure, N, as `radius', `colat', and `slon'.
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
%   1) Compute the spherical coordinates of the position of the Moon
%      as seen from the Earth, and convert them to cylindrical and
%      rectangular coordinates.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: sphcyl_ex1.tm
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
%      function sphcyl_ex1()
%
%         %
%         % Load an SPK and leapseconds kernels.
%         %
%         cspice_furnsh( 'sphcyl_ex1.tm' )
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
%         % Convert the position vector `pos' to spherical
%         % coordinates.
%         %
%         [radius, colat, slon]  = cspice_recsph(pos);
%         fprintf( '\n' )
%         fprintf( 'Spherical coordinates:\n' )
%         fprintf( '   Radius       (km): %20.8f\n', radius )
%         fprintf( '   Polar Angle (deg): %20.8f\n', colat * cspice_dpr )
%         fprintf( '   Longitude   (deg): %20.8f\n', slon  * cspice_dpr )
%
%         %
%         % Convert the spherical coords to cylindrical.
%         %
%         [r, lon, z]           = cspice_sphcyl(radius, colat, slon);
%         fprintf( '\n' )
%         fprintf( 'Cylindrical coordinates:\n' )
%         fprintf( '   Radius       (km): %20.8f\n', r )
%         fprintf( '   Longitude   (deg): %20.8f\n', lon * cspice_dpr )
%         fprintf( '   Z            (km): %20.8f\n', z )
%
%         %
%         % Convert the cylindrical to rectangular.
%         %
%         [rectan]              = cspice_cylrec(r, lon, z);
%         fprintf( '\n' )
%         fprintf( 'Rectangular coordinates from cspice_cylrec:\n' )
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
%      Spherical coordinates:
%         Radius       (km):      403626.33912495
%         Polar Angle (deg):         108.26566077
%         Longitude   (deg):         -98.34959789
%
%      Cylindrical coordinates:
%         Radius       (km):      383289.01777726
%         Longitude   (deg):         -98.34959789
%         Z            (km):     -126505.93063865
%
%      Rectangular coordinates from cspice_cylrec:
%         X            (km):      -55658.44323296
%         Y            (km):     -379226.32931475
%         Z            (km):     -126505.93063865
%
%
%   2) Create a table showing a variety of spherical coordinates
%      and the corresponding cylindrical coordinates.
%
%      Corresponding spherical and cylindrical coordinates are
%      listed to three decimal places. Input and output angles are
%      in degrees.
%
%
%      Example code begins here.
%
%
%      function sphcyl_ex2()
%
%         %
%         % Define six sets of spherical coordinates, `slon' and `colat'
%         % expressed in degrees - converted to radians by use of cspice_rpd.
%         %
%         radius = [  1.,  1., 1.4142, 1.4142, 1.  , 0. ];
%         colat  = [ 90., 90., 45.   , 135.  , 0.  , 0. ] * cspice_rpd;
%         slon   = [  0., 90., 180.  , 180.  , 180., 33.] * cspice_rpd;
%
%         %
%         % ...convert the spherical coordinates to cylindrical coordinates
%         %
%         [r, clon, z] = cspice_sphcyl(radius, colat, slon);
%
%        %
%         % ...convert angular measure to degrees.
%         %
%         colat = colat * cspice_dpr;
%         clon = clon   * cspice_dpr;
%         slon = slon   * cspice_dpr;
%
%         %
%         % Output banner.
%         %
%         disp('    r       clon      z     radius     slon    colat' )
%         disp(' -------  -------  -------  -------  -------  -------')
%
%
%         %
%         % Create an array of values for output.
%         %
%         output = [ r; clon; z; radius; slon; colat ];
%         txt   = sprintf( '%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', output);
%         disp( txt )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%          r       clon      z     radius     slon    colat
%       -------  -------  -------  -------  -------  -------
%         1.000    0.000    0.000    1.000    0.000   90.000
%         1.000   90.000    0.000    1.000   90.000   90.000
%         1.000  180.000    1.000    1.414  180.000   45.000
%         1.000  180.000   -1.000    1.414  180.000  135.000
%         0.000  180.000    1.000    1.000  180.000    0.000
%         0.000   33.000    0.000    0.000   33.000    0.000
%
%
%   3) Other than the obvious conversion between coordinate systems
%      this routine could be used to obtain the axial projection
%      from a sphere to a cylinder about the z-axis that contains
%      the equator of the sphere.
%
%      Such a projection is valuable because it preserves the
%      areas between regions on the sphere and their projections to
%      the cylinder.
%
%
%      Example code begins here.
%
%
%      function sphcyl_ex3()
%
%         %
%         % Define the point whose projection is to be
%         % computed.
%         %
%         radius =   100.0;
%         slon   =    45.0  * cspice_rpd;
%         colat  =   102.5 * cspice_rpd;
%
%         %
%         % Convert the spherical coordinates to cylindrical.
%         %
%         [r, clon, z] = cspice_sphcyl( radius, colat, slon );
%
%         fprintf( 'Coordinates of the projected point on cylinder:\n' )
%         fprintf( ' \n' )
%         fprintf( ' Radius     (km):  %22.11f\n', r )
%         fprintf( ' Longitude (deg):  %22.11f\n', clon*cspice_dpr )
%         fprintf( ' Z          (km):  %22.11f\n', z )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Coordinates of the projected point on cylinder:
%
%       Radius     (km):          97.62960071199
%       Longitude (deg):          45.00000000000
%       Z          (km):         -21.64396139381
%
%
%-Particulars
%
%   This returns the cylindrical coordinates of a point whose
%   position is input through spherical coordinates.
%
%-Exceptions
%
%   1)  If any of the input arguments, `radius', `colat' or `slon', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `radius', `colat' or `slon', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   3)  If the input vectorizable arguments `radius', `colat' and
%       `slon' do not have the same measure of vectorization (N), an
%       error is signaled by the Mice interface.
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
%       Changed the output argument name "lonc" to "clon".
%
%       Edited the header to comply with NAIF standard. Added
%       meta-kernel to example #1. Updated code example #1 to produce
%       formatted output and added a call to cspice_kclear. Added the
%       problem statement to both examples and a third example.
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
%   -Mice Version 1.0.0, 12-DEC-2005 (EDW)
%
%-Index_Entries
%
%   spherical to cylindrical coordinates
%
%-&

function [ r, clon, z] = cspice_sphcyl(radius, colat, slon)

   switch nargin
      case 3

         radius = zzmice_dp(radius);
         colat  = zzmice_dp(colat);
         slon   = zzmice_dp(slon);

      otherwise

         error ( ['Usage: [ _r_, _clon_, _z_] = '...
                  'cspice_sphcyl(_radius_, _colat_, _slon_)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [ r, clon, z] = mice('sphcyl_c', radius, colat, slon );
   catch spiceerr
      rethrow(spiceerr)
   end


