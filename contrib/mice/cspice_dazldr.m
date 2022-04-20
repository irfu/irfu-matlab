%-Abstract
%
%   CSPICE_DAZLDR computes the Jacobian matrix of the transformation from
%   rectangular to azimuth/elevation coordinates.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      x,
%      y,
%      z        the rectangular coordinates of the point at which the
%               Jacobian matrix of the map from rectangular to
%               azimuth/elevation coordinates is desired.
%
%               [1,1] = size(x); double = class(x)
%               [1,1] = size(y); double = class(y)
%               [1,1] = size(z); double = class(z)
%
%      azccw    a flag indicating how the azimuth is measured.
%
%               [1,1] = size(azccw); logical = class(azccw)
%
%               If `azccw' is true, the azimuth increases in the
%               counterclockwise direction; otherwise it increases
%               in the clockwise direction.
%
%      elplsz   a flag indicating how the elevation is measured.
%
%               [1,1] = size(elplsz); logical = class(elplsz)
%
%               If `elplsz' is true, the elevation increases from the
%               XY plane toward +Z; otherwise toward -Z.
%
%   the call:
%
%      [jacobi] = cspice_dazldr( x, y, z, azccw, elplsz )
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the transformation from
%               rectangular to azimuth/elevation coordinates.
%
%               [3,3] = size(jacobi); double = class(jacobi)
%
%               It has the form
%
%                  .-                            -.
%                  |  dr/dx     dr/dy     dr/dz   |
%                  |  daz/dx    daz/dy    daz/dz  |
%                  |  del/dx    del/dy    del/dz  |
%                  `-                            -'
%
%                evaluated at the input values of `x', `y', and `z'.
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
%   1) Find the azimuth/elevation state of Venus as seen from the
%      DSS-14 station at a given epoch. Map this state back to
%      rectangular coordinates as a check.
%
%      Task description
%      ================
%
%      In this example, we will obtain the apparent state of Venus as
%      seen from the DSS-14 station in the DSS-14 topocentric
%      reference frame. We will use a station frames kernel and
%      transform the resulting rectangular coordinates to azimuth,
%      elevation and range and its derivatives using cspice_recazl and
%      cspice_dazldr.
%
%      We will map this state back to rectangular coordinates using
%      cspice_azlrec and cspice_drdazl.
%
%      In order to introduce the usage of the logical flags `azccw'
%      and `elplsz', we will request the azimuth to be measured
%      clockwise and the elevation positive towards +Z
%      axis of the DSS-14_TOPO reference frame.
%
%      Kernels
%      =======
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: dazldr_ex1.tm
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
%            File name                        Contents
%            ---------                        --------
%            de430.bsp                        Planetary ephemeris
%            naif0011.tls                     Leapseconds
%            earth_720101_070426.bpc          Earth historical
%                                                binary PCK
%            earthstns_itrf93_050714.bsp      DSN station SPK
%            earth_topo_050714.tf             DSN station FK
%
%         \begindata
%
%         KERNELS_TO_LOAD = ( 'de430.bsp',
%                             'naif0011.tls',
%                             'earth_720101_070426.bpc',
%                             'earthstns_itrf93_050714.bsp',
%                             'earth_topo_050714.tf'         )
%
%         \begintext
%
%         End of meta-kernel.
%
%
%      Example code begins here.
%
%
%      function dazldr_ex1()
%
%         %
%         % Local parameters
%         %
%         META =   'dazldr_ex1.tm';
%
%         %
%         % Load SPICE kernels.
%         %
%         cspice_furnsh( META );
%
%         %
%         % Convert the observation time to seconds past J2000 TDB.
%         %
%         obstim = '2003 OCT 13 06:00:00.000000 UTC';
%
%         [et]   = cspice_str2et( obstim );
%
%         %
%         % Set the target, observer, observer frame, and
%         % aberration corrections.
%         %
%         target = 'VENUS';
%         obs    = 'DSS-14';
%         ref    = 'DSS-14_TOPO';
%         abcorr = 'CN+S';
%
%         %
%         % Compute the observer-target state.
%         %
%         [state, lt] = cspice_spkezr( target, et, ref, abcorr, obs );
%
%         %
%         % Convert position to azimuth/elevation coordinates,
%         % with azimuth increasing clockwise and elevation
%         % positive towards +Z axis of the DSS-14_TOPO
%         % reference frame
%         %
%         azccw  = false;
%         elplsz = true;
%
%         [r, az, el] = cspice_recazl( state(1:3), azccw, elplsz );
%
%         %
%         % Convert velocity to azimuth/elevation coordinates.
%         %
%         [jacobi] = cspice_dazldr( state(1), state(2), state(3),          ...
%                                   azccw,    elplsz              );
%
%         azlvel   = jacobi * state(4:6);
%
%         %
%         % As a check, convert the azimuth/elevation state back to
%         % rectangular coordinates.
%         %
%         [rectan] = cspice_azlrec( r, az, el, azccw, elplsz );
%
%         [jacobi] = cspice_drdazl( r, az, el, azccw, elplsz );
%
%         drectn   = jacobi * azlvel;
%
%         fprintf( '\n' )
%         fprintf( 'AZ/EL coordinates:\n' )
%         fprintf( '\n' )
%         fprintf( '   Range      (km)        =  %19.8f\n', r )
%         fprintf( '   Azimuth    (deg)       =  %19.8f\n', az * cspice_dpr )
%         fprintf( '   Elevation  (deg)       =  %19.8f\n', el * cspice_dpr )
%         fprintf( '\n' )
%         fprintf( 'AZ/EL velocity:\n' )
%         fprintf( '\n' )
%         fprintf( '   d Range/dt     (km/s)  =  %19.8f\n', azlvel(1) )
%         fprintf( '   d Azimuth/dt   (deg/s) =  %19.8f\n',                ...
%                                    azlvel(2) * cspice_dpr )
%         fprintf( '   d Elevation/dt (deg/s) =  %19.8f\n',                ...
%                                    azlvel(3) * cspice_dpr )
%         fprintf( '\n' )
%         fprintf( 'Rectangular coordinates:\n' )
%         fprintf( '\n' )
%         fprintf( '   X (km)                 =  %19.8f\n', state(1) )
%         fprintf( '   Y (km)                 =  %19.8f\n', state(2) )
%         fprintf( '   Z (km)                 =  %19.8f\n', state(3) )
%         fprintf( '\n' )
%         fprintf( 'Rectangular velocity:\n' )
%         fprintf( '\n' )
%         fprintf( '   dX/dt (km/s)           =  %19.8f\n', state(4) )
%         fprintf( '   dY/dt (km/s)           =  %19.8f\n', state(5) )
%         fprintf( '   dZ/dt (km/s)           =  %19.8f\n', state(6) )
%         fprintf( '\n' )
%         fprintf( 'Rectangular coordinates from inverse mapping:\n' )
%         fprintf( '\n' )
%         fprintf( '   X (km)                 =  %19.8f\n', rectan(1) )
%         fprintf( '   Y (km)                 =  %19.8f\n', rectan(2) )
%         fprintf( '   Z (km)                 =  %19.8f\n', rectan(3) )
%         fprintf( '\n' )
%         fprintf( 'Rectangular velocity from inverse mapping:\n' )
%         fprintf( '\n' )
%         fprintf( '   dX/dt (km/s)           =  %19.8f\n', drectn(1) )
%         fprintf( '   dY/dt (km/s)           =  %19.8f\n', drectn(2) )
%         fprintf( '   dZ/dt (km/s)           =  %19.8f\n', drectn(3) )
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
%      AZ/EL coordinates:
%
%         Range      (km)        =   245721478.99272084
%         Azimuth    (deg)       =         294.48543372
%         Elevation  (deg)       =         -48.94609726
%
%      AZ/EL velocity:
%
%         d Range/dt     (km/s)  =          -4.68189834
%         d Azimuth/dt   (deg/s) =           0.00402256
%         d Elevation/dt (deg/s) =          -0.00309156
%
%      Rectangular coordinates:
%
%         X (km)                 =    66886767.37916667
%         Y (km)                 =   146868551.77222887
%         Z (km)                 =  -185296611.10841590
%
%      Rectangular velocity:
%
%         dX/dt (km/s)           =        6166.04150307
%         dY/dt (km/s)           =      -13797.77164550
%         dZ/dt (km/s)           =       -8704.32385654
%
%      Rectangular coordinates from inverse mapping:
%
%         X (km)                 =    66886767.37916658
%         Y (km)                 =   146868551.77222890
%         Z (km)                 =  -185296611.10841590
%
%      Rectangular velocity from inverse mapping:
%
%         dX/dt (km/s)           =        6166.04150307
%         dY/dt (km/s)           =      -13797.77164550
%         dZ/dt (km/s)           =       -8704.32385654
%
%
%-Particulars
%
%   When performing vector calculations with velocities it is
%   usually most convenient to work in rectangular coordinates.
%   However, once the vector manipulations have been performed
%   it is often desirable to convert the rectangular representations
%   into azimuth/elevation coordinates to gain insights about
%   phenomena in this coordinate system.
%
%   To transform rectangular velocities to derivatives of coordinates
%   in an azimuth/elevation coordinate system, one uses the Jacobian
%   matrix of the transformation between the two systems.
%
%   Given a state in rectangular coordinates
%
%      ( x, y, z, dx, dy, dz )
%
%   the corresponding azimuth/elevation coordinate derivatives are
%   given by the matrix equation:
%
%                    t          |                      t
%      (dr, daz, del)   = jacobi|        * (dx, dy, dz)
%                               |(x,y,z)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(x, y, z)
%
%   In the azimuth/elevation coordinate system, several conventions
%   exist on how azimuth and elevation are measured. Using the `azccw'
%   and `elplsz' flags, users indicate which conventions shall be used.
%   See the descriptions of these input arguments for details.
%
%-Exceptions
%
%   1)  If the input point is on the Z-axis ( x = 0 and y = 0 ), the
%       Jacobian matrix is undefined and therefore, the error
%       SPICE(POINTONZAXIS) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If any of the input arguments, `x', `y', `z', `azccw' or
%       `elplsz', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   3)  If any of the input arguments, `x', `y', `z', `azccw' or
%       `elplsz', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
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
%
%-Version
%
%   -Mice Version 1.0.0, 08-FEB-2021 (JDR)
%
%-Index_Entries
%
%   Jacobian matrix of AZ/EL w.r.t. rectangular coordinates
%   Rectangular to range, azimuth and elevation derivative
%   Rectangular to range, AZ and EL velocity conversion
%
%-&
function [jacobi] = cspice_dazldr( x, y, z, azccw, elplsz )

   switch nargin
      case 5

         x      = zzmice_dp(x);
         y      = zzmice_dp(y);
         z      = zzmice_dp(z);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);

      otherwise

         error ( [ 'Usage: [jacobi(3,3)] = '                                ...
                   'cspice_dazldr( x, y, z, azccw, elplsz )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('dazldr_c', x, y, z, azccw, elplsz);
   catch spiceerr
      rethrow(spiceerr)
   end
