%-Abstract
%
%   CSPICE_DRDSPH computes the Jacobian matrix of the transformation from
%   spherical to rectangular coordinates.
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
%      r        the distance(s) of a point(s) from the origin.
%
%               [1,n] = size(r); double = class(r)
%
%      colat    the angle(s) between the point(s) and the positive z-axis, in
%               radians.
%
%               [1,n] = size(colat); double = class(colat)
%
%      slon     the angle(s) of the point(s) measured from the xz plane in
%               radians.
%
%               [1,n] = size(slon); double = class(slon)
%
%               The angle increases in the counterclockwise sense about
%               the +z axis.
%
%   the call:
%
%      [jacobi] = cspice_drdsph( r, colat, slon )
%
%   returns:
%
%      jacobi   the matrix(es) of partial derivatives of the conversion between
%               spherical and rectangular coordinates, evaluated at the input
%               coordinates.
%
%               If [1,1] = size(r) then [3,3]   = size(jacobi)
%               If [1,n] = size(r) then [3,3,n] = size(jacobi)
%                                        double = class(jacobi)
%
%               This matrix has the form:
%
%                  .-                                  -.
%                  |  dx/dr     dx/dcolat     dx/dslon  |
%                  |                                    |
%                  |  dy/dr     dy/dcolat     dy/dslon  |
%                  |                                    |
%                  |  dz/dr     dz/dcolat     dz/dslon  |
%                  `-                                  -'
%
%               evaluated at the input values of `r', `slon' and `lat'.
%               Here `x', `y', and `z' are given by the familiar formulae
%
%                  x = r*cos(slon)*sin(colat)
%                  y = r*sin(slon)*sin(colat)
%                  z = r*cos(colat)
%
%               `jacobi' returns with the same vectorization measure (N)
%               as `r', `colat' and `slon'.
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
%   1) Find the spherical state of the Earth as seen from
%      Mars in the IAU_MARS reference frame at January 1, 2005 TDB.
%      Map this state back to rectangular coordinates as a check.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: drdsph_ex1.tm
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
%      function drdsph_ex1()
%
%         %
%         % Load SPK, PCK and LSK kernels, use a meta kernel for
%         % convenience.
%         %
%         cspice_furnsh( 'drdsph_ex1.tm' );
%
%         %
%         % Look up the apparent state of earth as seen from Mars
%         % at January 1, 2005 TDB, relative to the IAU_MARS reference
%         % frame.
%         %
%         [et] = cspice_str2et( 'January 1, 2005 TDB' );
%
%         [state, lt] = cspice_spkezr( 'Earth', et,    'IAU_MARS',         ...
%                                      'LT+S',  'Mars'           );
%
%         %
%         % Convert position to spherical coordinates.
%         %
%         [r, colat, slon] = cspice_recsph( state(1:3) );
%
%         %
%         % Convert velocity to spherical coordinates.
%         %
%         [jacobi] = cspice_dsphdr( state(1), state(2), state(3) );
%
%         sphvel   = jacobi * state(4:6);
%
%         %
%         % As a check, convert the spherical state back to
%         % rectangular coordinates.
%         %
%         [rectan] = cspice_sphrec( r, colat, slon );
%
%         [jacobi] = cspice_drdsph( r, colat, slon );
%
%         drectn   = jacobi * sphvel;
%
%         fprintf( ' \n' )
%         fprintf( 'Rectangular coordinates:\n' )
%         fprintf( ' \n' )
%         fprintf( ' X (km)                  =  %17.8e\n', state(1) )
%         fprintf( ' Y (km)                  =  %17.8e\n', state(2) )
%         fprintf( ' Z (km)                  =  %17.8e\n', state(3) )
%         fprintf( ' \n' )
%         fprintf( 'Rectangular velocity:\n' )
%         fprintf( ' \n' )
%         fprintf( ' dX/dt (km/s)            =  %17.8e\n', state(4) )
%         fprintf( ' dY/dt (km/s)            =  %17.8e\n', state(5) )
%         fprintf( ' dZ/dt (km/s)            =  %17.8e\n', state(6) )
%         fprintf( ' \n' )
%         fprintf( 'Spherical coordinates:\n' )
%         fprintf( ' \n' )
%         fprintf( ' Radius     (km)         =  %17.8e\n', r )
%         fprintf( ' Colatitude (deg)        =  %17.8e\n', colat/cspice_rpd )
%         fprintf( ' Longitude  (deg)        =  %17.8e\n', slon/cspice_rpd )
%         fprintf( ' \n' )
%         fprintf( 'Spherical velocity:\n' )
%         fprintf( ' \n' )
%         fprintf( ' d Radius/dt     (km/s)  =  %17.8e\n', sphvel(1) )
%         fprintf( ' d Colatitude/dt (deg/s) =  %17.8e\n',                 ...
%                                     sphvel(2)/cspice_rpd )
%         fprintf( ' d Longitude/dt  (deg/s) =  %17.8e\n',                 ...
%                                     sphvel(3)/cspice_rpd )
%         fprintf( ' \n' )
%         fprintf( 'Rectangular coordinates from inverse mapping:\n' )
%         fprintf( ' \n' )
%         fprintf( ' X (km)                  =  %17.8e\n', rectan(1) )
%         fprintf( ' Y (km)                  =  %17.8e\n', rectan(2) )
%         fprintf( ' Z (km)                  =  %17.8e\n', rectan(3) )
%         fprintf( ' \n' )
%         fprintf( 'Rectangular velocity from inverse mapping:\n' )
%         fprintf( ' \n' )
%         fprintf( ' dX/dt (km/s)            =  %17.8e\n', drectn(1) )
%         fprintf( ' dY/dt (km/s)            =  %17.8e\n', drectn(2) )
%         fprintf( ' dZ/dt (km/s)            =  %17.8e\n', drectn(3) )
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
%       X (km)                  =    -7.60961826e+07
%       Y (km)                  =     3.24363805e+08
%       Z (km)                  =     4.74704840e+07
%
%      Rectangular velocity:
%
%       dX/dt (km/s)            =     2.29520749e+04
%       dY/dt (km/s)            =     5.37601112e+03
%       dZ/dt (km/s)            =    -2.08811490e+01
%
%      Spherical coordinates:
%
%       Radius     (km)         =     3.36535219e+08
%       Colatitude (deg)        =     8.18910134e+01
%       Longitude  (deg)        =     1.03202903e+02
%
%      Spherical velocity:
%
%       d Radius/dt     (km/s)  =    -1.12116011e+01
%       d Colatitude/dt (deg/s) =     3.31899303e-06
%       d Longitude/dt  (deg/s) =    -4.05392876e-03
%
%      Rectangular coordinates from inverse mapping:
%
%       X (km)                  =    -7.60961826e+07
%       Y (km)                  =     3.24363805e+08
%       Z (km)                  =     4.74704840e+07
%
%      Rectangular velocity from inverse mapping:
%
%       dX/dt (km/s)            =     2.29520749e+04
%       dY/dt (km/s)            =     5.37601112e+03
%       dZ/dt (km/s)            =    -2.08811490e+01
%
%
%-Particulars
%
%   It is often convenient to describe the motion of an object in
%   the spherical coordinate system. However, when performing
%   vector computations its hard to beat rectangular coordinates.
%
%   To transform states given with respect to spherical coordinates
%   to states with respect to rectangular coordinates, one uses
%   the Jacobian of the transformation between the two systems.
%
%   Given a state in spherical coordinates
%
%      ( r, colat, slon, dr, dcolat, dslon )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation:
%                  t          |                                    t
%      (dx, dy, dz)   = jacobi|              * (dr, dcolat, dslon )
%                             |(r,colat,slon)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(r,colat,slon)
%
%-Exceptions
%
%   1)  If any of the input arguments, `r', `colat' or `slon', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `r', `colat' or `slon', is not
%       of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   3)  If the input vectorizable arguments `r', `colat' and `slon' do
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
%       Changed the input argument name "lon" to "slon" for consistency
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
%   Jacobian of rectangular w.r.t. spherical coordinates
%
%-&

function [jacobi] = cspice_drdsph( r, colat, slon )

   switch nargin
      case 3

         r     = zzmice_dp(r);
         colat = zzmice_dp(colat);
         slon  = zzmice_dp(slon);

      otherwise

         error( ['Usage: [_jacobi(3,3)_] = ' ...
                         'cspice_drdsph( _r_, _colat_, _slon_ )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdsph_c', r, colat, slon);
   catch spiceerr
      rethrow(spiceerr)
   end
