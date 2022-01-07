%-Abstract
%
%   CSPICE_DRDCYL computes the Jacobian matrix of the transformation from
%   cylindrical to rectangular coordinates.
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
%      r        distance of the point of interest from z axis.
%
%               [1,n] = size(r); double = class(r)
%
%      clon     cylindrical angle (in radians) of the point of interest from
%               the XZ plane.
%
%               [1,n] = size(clon); double = class(clon)
%
%               The angle increases in the counterclockwise sense about the
%               +Z axis.
%
%      z        height of the point above XY plane.
%
%               [1,n] = size(z); double = class(z)
%
%   the call:
%
%      [jacobi] = cspice_drdcyl( r, clon, z )
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion between
%               cylindrical and rectangular coordinates.
%
%               If [1,1] = size(r) then [3,3]   = size(jacobi)
%               If [1,n] = size(r) then [3,3,n] = size(jacobi)
%                                        double = class(jacobi)
%
%               It has the form
%
%                  .-                                -.
%                  |  dx/dr     dx/dclon     dx/dz    |
%                  |                                  |
%                  |  dy/dr     dy/dclon     dy/dz    |
%                  |                                  |
%                  |  dz/dr     dz/dclon     dz/dz    |
%                  `-                                -'
%
%               evaluated at the input values of `r', `clon' and `z'.
%               Here `x', `y', and `z' are given by the familiar formulae
%
%                  x = r*cos(clon)
%                  y = r*sin(clon)
%                  z = z
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
%   1) Find the cylindrical state of the Earth as seen from
%      Mars in the IAU_MARS reference frame at January 1, 2005 TDB.
%      Map this state back to rectangular coordinates as a check.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: drdcyl_ex1.tm
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
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function drdcyl_ex1()
%
%         %
%         % Load SPK, PCK and LSK kernels, use a meta kernel for
%         % convenience.
%         %
%         cspice_furnsh( 'drdcyl_ex1.tm' );
%
%         %
%         % Look up the apparent state of earth as seen from Mars
%         % at January 1, 2005 TDB, relative to the IAU_MARS reference
%         % frame.
%         %
%         [et] = cspice_str2et( 'January 1, 2005 TDB' );
%
%         [state, lt] = cspice_spkezr( 'Earth', et,    'IAU_MARS',         ...
%                                      'LT+S',  'Mars'          );
%
%         %
%         % Convert position to cylindrical coordinates.
%         %
%         [r, clon, z] = cspice_reccyl( state(1:3) );
%
%         %
%         % Convert velocity to cylindrical coordinates.
%         %
%         [jacobi] = cspice_dcyldr( state(1), state(2), state(3) );
%
%         cylvel   = jacobi * state(4:6);
%
%         %
%         % As a check, convert the cylindrical state back to
%         % rectangular coordinates.
%         %
%         [rectan] = cspice_cylrec( r, clon, z );
%
%         [jacobi] = cspice_drdcyl( r, clon, z );
%
%         drectn   = jacobi * cylvel;
%
%         fprintf( ' \n' )
%         fprintf( 'Rectangular coordinates:\n' )
%         fprintf( ' \n' )
%         fprintf( ' X (km)                 =  %17.8e\n', state(1) )
%         fprintf( ' Y (km)                 =  %17.8e\n', state(2) )
%         fprintf( ' Z (km)                 =  %17.8e\n', state(3) )
%         fprintf( ' \n' )
%         fprintf( 'Rectangular velocity:\n' )
%         fprintf( ' \n' )
%         fprintf( ' dX/dt (km/s)           =  %17.8e\n', state(4) )
%         fprintf( ' dY/dt (km/s)           =  %17.8e\n', state(5) )
%         fprintf( ' dZ/dt (km/s)           =  %17.8e\n', state(6) )
%         fprintf( ' \n' )
%         fprintf( 'Cylindrical coordinates:\n' )
%         fprintf( ' \n' )
%         fprintf( ' Radius    (km)         =  %17.8e\n', r )
%         fprintf( ' Longitude (deg)        =  %17.8e\n', clon/cspice_rpd )
%         fprintf( ' Z         (km)         =  %17.8e\n', z )
%         fprintf( ' \n' )
%         fprintf( 'Cylindrical velocity:\n' )
%         fprintf( ' \n' )
%         fprintf( ' d Radius/dt    (km/s)  =  %17.8e\n', cylvel(1) )
%         fprintf( ' d Longitude/dt (deg/s) =  %17.8e\n',                  ...
%                                                     cylvel(2)/cspice_rpd )
%         fprintf( ' d Z/dt         (km/s)  =  %17.8e\n', cylvel(3) )
%         fprintf( ' \n' )
%         fprintf( 'Rectangular coordinates from inverse mapping:\n' )
%         fprintf( ' \n' )
%         fprintf( ' X (km)                 =  %17.8e\n', rectan(1) )
%         fprintf( ' Y (km)                 =  %17.8e\n', rectan(2) )
%         fprintf( ' Z (km)                 =  %17.8e\n', rectan(3) )
%         fprintf( ' \n' )
%         fprintf( 'Rectangular velocity from inverse mapping:\n' )
%         fprintf( ' \n' )
%         fprintf( ' dX/dt (km/s)           =  %17.8e\n', drectn(1) )
%         fprintf( ' dY/dt (km/s)           =  %17.8e\n', drectn(2) )
%         fprintf( ' dZ/dt (km/s)           =  %17.8e\n', drectn(3) )
%         fprintf( ' \n' )
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
%       X (km)                 =    -7.60961826e+07
%       Y (km)                 =     3.24363805e+08
%       Z (km)                 =     4.74704840e+07
%
%      Rectangular velocity:
%
%       dX/dt (km/s)           =     2.29520749e+04
%       dY/dt (km/s)           =     5.37601112e+03
%       dZ/dt (km/s)           =    -2.08811490e+01
%
%      Cylindrical coordinates:
%
%       Radius    (km)         =     3.33170387e+08
%       Longitude (deg)        =     1.03202903e+02
%       Z         (km)         =     4.74704840e+07
%
%      Cylindrical velocity:
%
%       d Radius/dt    (km/s)  =    -8.34966283e+00
%       d Longitude/dt (deg/s) =    -4.05392876e-03
%       d Z/dt         (km/s)  =    -2.08811490e+01
%
%      Rectangular coordinates from inverse mapping:
%
%       X (km)                 =    -7.60961826e+07
%       Y (km)                 =     3.24363805e+08
%       Z (km)                 =     4.74704840e+07
%
%      Rectangular velocity from inverse mapping:
%
%       dX/dt (km/s)           =     2.29520749e+04
%       dY/dt (km/s)           =     5.37601112e+03
%       dZ/dt (km/s)           =    -2.08811490e+01
%
%
%-Particulars
%
%   It is often convenient to describe the motion of an object in
%   the cylindrical coordinate system. However, when performing
%   vector computations its hard to beat rectangular coordinates.
%
%   To transform states given with respect to cylindrical coordinates
%   to states with respect to rectangular coordinates, one uses
%   the Jacobian of the transformation between the two systems.
%
%   Given a state in cylindrical coordinates
%
%      ( r, clon, z, dr, dclon, dz )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation:
%                  t          |                           t
%      (dx, dy, dz)   = jacobi|          * (dr, dclon, dz)
%                             |(r,clon,z)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(r,clon,z)
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
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Added complete code example.
%
%       Changed the input argument name "lon" to "clon" for consistency
%       with other routines.
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
%   Jacobian of rectangular w.r.t. cylindrical coordinates
%
%-&

function [jacobi] = cspice_drdcyl( r, clon, z )

   switch nargin
      case 3

         r    = zzmice_dp(r);
         clon = zzmice_dp(clon);
         z    = zzmice_dp(z);

      otherwise

         error( 'Usage: [_jacobi(3,3)_] = cspice_drdcyl( _r_, _clon_, _z_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdcyl_c', r, clon, z);
   catch spiceerr
      rethrow(spiceerr)
   end




