%-Abstract
%
%   CSPICE_SXFORM returns the state transformation matrix from one
%   frame to another at a specified epoch.
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
%      from     the name of a reference frame in which a state is known.
%
%               [1,c1] = size(from); char = class(from)
%
%                  or
%
%               [1,1] = size(from); cell = class(from)
%
%      to       the name of a reference frame in which it is desired to
%               represent the state.
%
%               [1,c2] = size(to); char = class(to)
%
%                  or
%
%               [1,1] = size(to); cell = class(to)
%
%      et       the epoch(s) in ephemeris seconds past the epoch of J2000
%               (TDB) at which the state transformation matrix should be
%               evaluated.
%
%               [1,n] = size(et); double = class(et)
%
%   the call:
%
%      [xform] = cspice_sxform( from, to, et )
%
%   returns:
%
%      xform    the state transformation matri(x|ces) that transforms states
%               from the reference frame `from' to the frame `to' at epoch
%               `et'.
%
%               If [1,1] = size(et) then [6,6]   = size(xform)
%               If [1,n] = size(et) then [6,6,n] = size(xform)
%                                        double = class(xform)
%
%               If (x, y, z, dx, dy, dz) is a state relative to the frame
%               `from' then the vector ( x', y', z', dx', dy', dz' ) is the
%               same state relative to the frame `to' at epoch `et'. Here the
%               vector ( x', y', z', dx', dy', dz' ) is defined by the
%               equation:
%
%                  .-   -.     .-          -.   .-  -.
%                  | x'  |     |            |   | x  |
%                  | y'  |     |            |   | y  |
%                  | z'  |  =  |   xform    |   | z  |
%                  | dx' |     |            |   | dx |
%                  | dy' |     |            |   | dy |
%                  | dz' |     |            |   | dz |
%                  `-   -'     `-          -'   `-  -'
%
%               `xform' returns with the same vectorization measure, N,
%               as `et'.
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
%   1) Suppose you have geodetic coordinates of a station on the
%      surface of Earth and that you need the inertial (J2000)
%      state of this station. The following code example
%      illustrates how to transform the geodetic state of the
%      station to a J2000 state.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: sxform_ex1.tm
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
%            pck00008.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'pck00008.tpc',
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
%      function sxform_ex1()
%
%         %
%         % Load the PCK and LSK kernels.
%         %
%         cspice_furnsh( 'sxform_ex1.tm' )
%
%         %
%         % Define a geodetic longitude, latitude, altitude
%         % coordinate set. These coordinates are defined in the
%         % non-inertial, earth fixed frame "IAU_EARTH".
%         %
%         lon = 118.25 * cspice_rpd;
%         lat = 34.05  * cspice_rpd;
%         alt = 0.;
%
%         %
%         % Define a UTC time of interest. Convert the 'utc' string
%         % to ephemeris time J2000.
%         %
%         utc = 'January 1, 1990';
%         et = cspice_str2et( utc );
%
%         %
%         % Retrieve the equatorial and polar axis of the earth (body 399).
%         %
%         abc = cspice_bodvrd( 'EARTH', 'RADII', 3 );
%         equatr =  abc(1);
%         polar  =  abc(3);
%
%         %
%         % Calculate the flattening factor for earth.
%         %
%         f =  ( equatr - polar  ) / equatr;
%
%         %
%         % Calculate the Cartesian coordinates on earth for the
%         % location at 'lon', 'lat', 'alt'.
%         %
%         estate = cspice_georec( lon, lat, alt, equatr, f);
%
%         %
%         % cspice_georec returned the position vector of the geodetic
%         % coordinates, but we want the state vector. Since it is a fixed
%         % location referenced in the "IAU_EARTH" frame, the location has
%         % no velocity. We need to extend estate to a 6-vector, the final
%         % three elements with value 0.d.
%         %
%         estate = [ estate; [0.; 0.; 0.] ];
%
%         %
%         % Retrieve the transformation matrix from "IAU_EARTH"
%         % to "J2000" at epoch 'et'.
%         %
%         xform = cspice_sxform( 'IAU_EARTH', 'J2000', et );
%
%         jstate = xform * estate;
%
%         utcstr = cspice_et2utc( et, 'C', 3 );
%         fprintf( 'Epoch                         : %s\n', utcstr)
%         fprintf(['Position in J2000 frame   (km):',      ...
%                  ' %10.4f  %10.4f  %10.4f\n'], jstate(1:3) );
%         fprintf(['Velocity in J2000 frame (km/s):',      ...
%                  ' %10.4f  %10.4f  %10.4f\n'], jstate(4:6) );
%
%         %
%         % Return the state transformation matrices from "IAU_EARTH"
%         % to "J2000" approximately every three months for the time
%         % interval February 1, 1990 to February 1, 1991 (UTC).
%         %
%         %
%         % Define the time bounds for the time interval,
%         % 1 year,  convert to ephemeris time J2000.
%         %
%         utc_bounds = strvcat( '1 Feb 1990', '1 Feb 1991' );
%         et_bounds = cspice_str2et( utc_bounds );
%
%         %
%         % Step in units of 3 months. 1 year -> 4 steps.
%         %
%         step = (et_bounds(2) - et_bounds(1) ) / 4.;
%
%         %
%         % Create an array of 4 ephemeris times starting at
%         % et_bound(1) in intervals of 'step'.
%         %
%         et = [0:3]*step + et_bounds(1);
%
%         %
%         % Convert the 4-vector of 'et' to an array of corresponding
%         % transformation matrices (dimensions (6,6,4) ).
%         %
%         xform = cspice_sxform( 'IAU_EARTH', 'J2000', et );
%
%         %
%         % Transform the Cartesian state vector from "IAU_EARTH"
%         % to "J2000".
%         %
%         utcstr = cspice_et2utc( et, 'C', 3 );
%         for i=1:4
%            jstate = xform(:,:,i) * estate;
%            disp (' ' )
%            fprintf( 'Epoch                         : %s\n', utcstr(i,:))
%            fprintf(['Position in J2000 frame   (km):',      ...
%                     ' %10.4f  %10.4f  %10.4f\n'], jstate(1:3) );
%            fprintf(['Velocity in J2000 frame (km/s):',      ...
%                     ' %10.4f  %10.4f  %10.4f\n'], jstate(4:6) );
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
%      Epoch                         : 1990 JAN 01 00:00:00.000
%      Position in J2000 frame   (km): -4131.4630  -3308.3707   3547.0215
%      Velocity in J2000 frame (km/s):     0.2412     -0.3010      0.0002
%
%      Epoch                         : 1990 FEB 01 00:00:00.000
%      Position in J2000 frame   (km): -1876.4713  -4947.4745   3549.2275
%      Velocity in J2000 frame (km/s):     0.3608     -0.1366      0.0003
%
%      Epoch                         : 1990 MAY 03 06:00:00.249
%      Position in J2000 frame   (km):  1875.1001   4945.4239   3552.8083
%      Velocity in J2000 frame (km/s):    -0.3606      0.1370     -0.0003
%
%      Epoch                         : 1990 AUG 02 12:00:00.502
%      Position in J2000 frame   (km): -1887.0731  -4943.3818   3549.3093
%      Velocity in J2000 frame (km/s):     0.3605     -0.1374      0.0003
%
%      Epoch                         : 1990 NOV 01 18:00:00.752
%      Position in J2000 frame   (km):  1886.0424   4941.3201   3552.7263
%      Velocity in J2000 frame (km/s):    -0.3603      0.1378     -0.0003
%
%
%-Particulars
%
%   This routine provides the user level interface for computing
%   state transformations from one reference frame to another.
%
%   Note that the reference frames may be inertial or non-inertial.
%   However, the user must take care that sufficient SPICE kernel
%   information is loaded to provide a complete state transformation
%   path from the `from' frame to the `to' frame.
%
%-Exceptions
%
%   1)  If sufficient information has not been supplied via loaded
%       SPICE kernels to compute the transformation between the two
%       frames, an error is signaled by a routine in the call tree of
%       this routine.
%
%   2)  If either frame `from' or `to' is not recognized, the error
%       SPICE(UNKNOWNFRAME) is signaled by a routine in the call tree
%       of this routine.
%
%   3)  If any of the input arguments, `from', `to' or `et', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `from', `to' or `et', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%   ROTATION.REQ
%   FRAMES.REQ
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
%       Edited the header to comply with NAIF standard. Extended -I/O section.
%
%       Added example's meta-kernel. Reduced the size of the array of times
%       used to generate the example's output and reformatted the output.
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
%   -Mice Version 1.0.2, 05-FEB-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   Find a state transformation matrix
%
%-&

function [xform] = cspice_sxform(from, to, et)

   switch nargin
      case 3

         from = zzmice_str(from);
         to   = zzmice_str(to);
         et   = zzmice_dp(et);

      otherwise

         error( ['Usage: [_xform(6,6)_] = cspice_sxform( ',           ...
                                              '`from`, `to`, _et_ )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [xform] = mice('sxform_c', from, to, et);
   catch spiceerr
      rethrow(spiceerr)
   end


