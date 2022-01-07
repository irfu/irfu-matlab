%-Abstract
%
%   CSPICE_AZLREC converts from range, azimuth and elevation of a point to
%   rectangular coordinates.
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
%      range    the distance of the point from the origin.
%
%               [1,1] = size(range); double = class(range)
%
%               The input should be in terms of the same units in which the
%               output is desired.
%
%               Although negative values for `range' are allowed, its
%               use may lead to undesired results. See the -Exceptions
%               section for a discussion on this topic.
%
%      az       the azimuth of the point.
%
%               [1,1] = size(az); double = class(az)
%
%               This is the angle between the projection onto the XY plane
%               of the vector from the origin to the point and the +X axis of
%               the reference frame. `az' is zero at the +X axis.
%
%               The way azimuth is measured depends on the value of
%               the logical flag `azccw'. See the descriptions of the
%               argument `azccw' for details.
%
%               The range (i.e., the set of allowed values) of `az' is
%               unrestricted. See the -Exceptions section for a
%               discussion on the `az' range.
%
%               Units are radians.
%
%      el       the elevation of the point.
%
%               [1,1] = size(el); double = class(el)
%
%               This is the angle between the vector from the origin to the
%               point and the XY plane. `el' is zero at the XY plane.
%
%               The way elevation is measured depends on the value of
%               the logical flag `elplsz'. See the descriptions of the
%               argument `elplsz' for details.
%
%               The range (i.e., the set of allowed values) of `el' is
%               [-pi/2, pi/2], but no error checking is done to ensure
%               that `el' is within this range. See the -Exceptions
%               section for a discussion on the `el' range.
%
%               Units are radians.
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
%               If `elplsz' is true, the elevation increases from
%               the XY plane toward +Z; otherwise toward -Z.
%
%   the call:
%
%      [rectan] = cspice_azlrec( range, az, el, azccw, elplsz )
%
%   returns:
%
%      rectan   an array containing the rectangular coordinates of the point.
%
%               [3,1] = size(rectan); double = class(rectan)
%
%               The units associated with the point are those
%               associated with the input `range'.
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
%   1) Create four tables showing a variety of azimuth/elevation
%      coordinates and the corresponding rectangular coordinates,
%      resulting from the different choices of the `azccw' and `elplsz'
%      flags.
%
%      Corresponding azimuth/elevation and rectangular coordinates
%      are listed to three decimal places. Input angles are in
%      degrees.
%
%
%      Example code begins here.
%
%
%      function azlrec_ex1()
%
%         %
%         % Local parameters.
%         %
%         NREC =   11;
%
%         %
%         % Define the input azimuth/elevation coordinates and the
%         % different choices of the `azccw' and `elplsz' flags.
%         %
%         range  = [ 0.0, 1.0,   1.0,   1.0,   1.0,   1.0,                 ...
%                    1.0, 1.414, 1.414, 1.414, 1.732 ]';
%
%         az     = [ 0.0,   0.0, 270.0,   0.0, 180.0, 90.0,                ...
%                    0.0, 315.0,   0.0, 270.0, 315.0 ]';
%
%         el     = [  0.0, 0.0,   0.0, -90.0,   0.0,   0.0,                ...
%                    90.0, 0.0, -45.0, -45.0, -35.264 ]';
%
%         azccw  = [false,  true]';
%         elplsz = [false,  true]';
%
%         %
%         % Create a table for each combination of `azccw' and `elplsz'.
%         %
%         for i=1:2
%
%            for j=1:2
%
%               %
%               % Display the flag settings.
%               %
%               msg   = 'AZCCW = #; ELPLSZ = #';
%               [msg] = cspice_repml( msg, '#', azccw(i), 'C' );
%               [msg] = cspice_repml( msg, '#', elplsz(j), 'C' );
%
%               fprintf( '\n' )
%               fprintf( '%s\n', msg )
%
%               %
%               % Print the banner.
%               %
%               fprintf( '\n' )
%               fprintf( [ '   RANGE      AZ       EL    rect(1) ',        ...
%                          ' rect(2)  rect(3)\n' ]                  )
%               fprintf( [ '  -------  -------  -------  ------- ',        ...
%                          ' -------  -------\n' ]                  )
%
%               %
%               % Do the conversion. Input angles in degrees.
%               %
%               for n=1:NREC
%
%                  raz      = az(n) * cspice_rpd;
%                  rel      = el(n) * cspice_rpd;
%
%                  [rectan] = cspice_azlrec( range(n), raz,       rel,     ...
%                                            azccw(i), elplsz(j)       );
%
%                  fprintf( '%9.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',       ...
%                           range(n),  az(n),     el(n),                   ...
%                           rectan(1), rectan(2), rectan(3)          )
%
%               end
%
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      AZCCW = False; ELPLSZ = False
%
%         RANGE      AZ       EL    rect(1)  rect(2)  rect(3)
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          1.000  270.000    0.000   -0.000    1.000    0.000
%          1.000    0.000  -90.000    0.000    0.000    1.000
%          1.000  180.000    0.000   -1.000   -0.000    0.000
%          1.000   90.000    0.000    0.000   -1.000    0.000
%          1.000    0.000   90.000    0.000    0.000   -1.000
%          1.414  315.000    0.000    1.000    1.000    0.000
%          1.414    0.000  -45.000    1.000    0.000    1.000
%          1.414  270.000  -45.000   -0.000    1.000    1.000
%          1.732  315.000  -35.264    1.000    1.000    1.000
%
%      AZCCW = False; ELPLSZ = True
%
%         RANGE      AZ       EL    rect(1)  rect(2)  rect(3)
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          1.000  270.000    0.000   -0.000    1.000    0.000
%          1.000    0.000  -90.000    0.000    0.000   -1.000
%          1.000  180.000    0.000   -1.000   -0.000    0.000
%          1.000   90.000    0.000    0.000   -1.000    0.000
%          1.000    0.000   90.000    0.000    0.000    1.000
%          1.414  315.000    0.000    1.000    1.000    0.000
%          1.414    0.000  -45.000    1.000    0.000   -1.000
%          1.414  270.000  -45.000   -0.000    1.000   -1.000
%          1.732  315.000  -35.264    1.000    1.000   -1.000
%
%      AZCCW = True; ELPLSZ = False
%
%         RANGE      AZ       EL    rect(1)  rect(2)  rect(3)
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          1.000  270.000    0.000   -0.000   -1.000    0.000
%          1.000    0.000  -90.000    0.000    0.000    1.000
%          1.000  180.000    0.000   -1.000    0.000    0.000
%          1.000   90.000    0.000    0.000    1.000    0.000
%          1.000    0.000   90.000    0.000    0.000   -1.000
%          1.414  315.000    0.000    1.000   -1.000    0.000
%          1.414    0.000  -45.000    1.000    0.000    1.000
%          1.414  270.000  -45.000   -0.000   -1.000    1.000
%          1.732  315.000  -35.264    1.000   -1.000    1.000
%
%      AZCCW = True; ELPLSZ = True
%
%         RANGE      AZ       EL    rect(1)  rect(2)  rect(3)
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          1.000  270.000    0.000   -0.000   -1.000    0.000
%          1.000    0.000  -90.000    0.000    0.000   -1.000
%          1.000  180.000    0.000   -1.000    0.000    0.000
%          1.000   90.000    0.000    0.000    1.000    0.000
%          1.000    0.000   90.000    0.000    0.000    1.000
%          1.414  315.000    0.000    1.000   -1.000    0.000
%          1.414    0.000  -45.000    1.000    0.000   -1.000
%          1.414  270.000  -45.000   -0.000   -1.000   -1.000
%          1.732  315.000  -35.264    1.000   -1.000   -1.000
%
%
%   2) Compute the right ascension and declination of the pointing
%      direction of DSS-14 station at a given epoch.
%
%      Task Description
%      ================
%
%      In this example, we will obtain the right ascension and
%      declination of the pointing direction of the DSS-14 station at
%      a given epoch, by converting the station's pointing direction
%      given in azimuth and elevation to rectangular coordinates
%      in the station topocentric reference frame and applying a
%      frame transformation from DSS-14_TOPO to J2000, in order to
%      finally obtain the corresponding right ascension and
%      declination of the pointing vector.
%
%      In order to introduce the usage of the logical flags `azccw'
%      and `elplsz', we will assume that the azimuth is measured
%      counterclockwise and the elevation negative towards +Z
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
%         File name: azlrec_ex2.tm
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
%           File name                        Contents
%           ---------                        --------
%           naif0011.tls                     Leapseconds
%           earth_720101_070426.bpc          Earth historical
%                                            binary PCK
%           earth_topo_050714.tf             DSN station FK
%
%         \begindata
%
%         KERNELS_TO_LOAD = ( 'naif0011.tls',
%                             'earth_720101_070426.bpc',
%                             'earth_topo_050714.tf'     )
%
%         \begintext
%
%         End of meta-kernel.
%
%
%      Example code begins here.
%
%
%      function azlrec_ex2()
%
%         %
%         % Local parameters
%         %
%         META =   'azlrec_ex2.tm';
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
%         % Set the local topocentric frame
%         %
%         ref = 'DSS-14_TOPO';
%
%         %
%         % Set the station's pointing direction in azimuth and
%         % elevation. Set arbitrarily the range to 1.0. Azimuth
%         % and elevation shall be given in radians. Azimuth
%         % increases counterclockwise and elevation is negative
%         % towards +Z (above the local horizon)
%         %
%         az     =   75.00;
%         el     =  -27.25;
%         azr    =   az * cspice_rpd;
%         elr    =   el * cspice_rpd;
%         r      =    1.00;
%         azccw  = true;
%         elplsz = false;
%
%         %
%         % Obtain the rectangular coordinates of the station's
%         % pointing direction.
%         %
%         [ptarg] = cspice_azlrec( r, azr, elr, azccw, elplsz );
%
%         %
%         % Transform the station's pointing vector from the
%         % local topocentric frame to J2000.
%         %
%         [rotate] = cspice_pxform( ref, 'J2000', et );
%         jpos     = rotate * ptarg;
%
%         %
%         % Compute the right ascension and declination.
%         % Express both angles in degrees.
%         %
%         [range, ra, dec] = cspice_recrad( jpos );
%         ra  =   ra * cspice_dpr;
%         dec =   dec * cspice_dpr;
%
%         %
%         % Display the computed pointing vector, the input
%         % data and resulting the angles.
%         %
%         fprintf( '\n' )
%         fprintf( 'Pointing azimuth    (deg):  %14.8f\n', az )
%         fprintf( 'Pointing elevation  (deg):  %14.8f\n', el )
%
%         [msg] = cspice_repml( 'Azimuth counterclockwise?: #', '#',       ...
%                               azccw, 'C'                           );
%         fprintf( '%s\n', msg )
%
%         [msg] = cspice_repml( 'Elevation positive +Z?   : #', '#',       ...
%                               elplsz, 'C'                          );
%         fprintf( '%s\n', msg )
%
%         fprintf( 'Observation epoch        : %s\n', obstim )
%         fprintf( '\n' )
%         fprintf( 'Pointing direction (normalized):  \n' )
%         fprintf( '   %14.8f %14.8f %14.8f\n',                            ...
%                  ptarg(1), ptarg(2), ptarg(3) )
%         fprintf( '\n' )
%         fprintf( 'Pointing right ascension (deg):  %14.8f\n', ra )
%         fprintf( 'Pointing declination (deg):      %14.8f\n', dec )
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
%      Pointing azimuth    (deg):     75.00000000
%      Pointing elevation  (deg):    -27.25000000
%      Azimuth counterclockwise?: True
%      Elevation positive +Z?   : False
%      Observation epoch        : 2003 OCT 13 06:00:00.000000 UTC
%
%      Pointing direction (normalized):
%             0.23009457     0.85872462     0.45787392
%
%      Pointing right ascension (deg):    280.06179939
%      Pointing declination (deg):         26.92826084
%
%
%-Particulars
%
%   This routine converts the azimuth, elevation, and range
%   of a point into the associated rectangular coordinates.
%
%   The input is defined by the distance from the center of
%   the reference frame (range), the angle from a reference
%   vector (azimuth), and the angle above the XY plane of the
%   reference frame (elevation).
%
%   The way azimuth and elevation are measured depends on the
%   values given by the user to the `azccw' and `elplsz' logical
%   flags. See the descriptions of these input arguments
%   for details.
%
%-Exceptions
%
%   1)  If the value of the input argument `range' is negative
%       the output rectangular coordinates will be negated, i.e.
%       the resulting array will be of the same length
%       but opposite direction to the one that would be obtained
%       with a positive input argument `range' of value ||RANGE||.
%
%   2)  If the value of the input argument `el' is outside the
%       range [-pi/2, pi/2], the results may not be as
%       expected.
%
%   3)  If the value of the input argument `az' is outside the
%       range [0, 2*pi], the value will be mapped to a value
%       inside the range that differs from the input value by an
%       integer multiple of 2*pi.
%
%   4)  If any of the input arguments, `range', `az', `el', `azccw' or
%       `elplsz', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   5)  If any of the input arguments, `range', `az', `el', `azccw' or
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
%   range, az and el to rectangular coordinates
%   range, azimuth and elevation to rectangular
%   convert range, az and el to rectangular coordinates
%   convert range, azimuth and elevation to rectangular
%
%-&
function [rectan] = cspice_azlrec( range, az, el, azccw, elplsz )

   switch nargin
      case 5

         range  = zzmice_dp(range);
         az     = zzmice_dp(az);
         el     = zzmice_dp(el);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);

      otherwise

         error ( [ 'Usage: [rectan(3)] = '                                  ...
                   'cspice_azlrec( range, az, el, azccw, elplsz )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice('azlrec_c', range, az, el, azccw, elplsz);
   catch spiceerr
      rethrow(spiceerr)
   end
