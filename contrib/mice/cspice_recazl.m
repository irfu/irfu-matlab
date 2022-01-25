%-Abstract
%
%   CSPICE_RECAZL converts rectangular coordinates of a point to range,
%   azimuth and elevation.
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
%      rectan   the rectangular coordinates of a point.
%
%               [3,1] = size(rectan); double = class(rectan)
%
%      azccw    a flag indicating how azimuth is measured.
%
%               [1,1] = size(azccw); logical = class(azccw)
%
%               If `azccw' is true, azimuth increases in the
%               counterclockwise direction; otherwise it increases in
%               the clockwise direction.
%
%      elplsz   a flag indicating how elevation is measured.
%
%               [1,1] = size(elplsz); logical = class(elplsz)
%
%               If `elplsz' is true, elevation increases from
%               the XY plane toward +Z; otherwise toward -Z.
%
%   the call:
%
%      [range, az, el] = cspice_recazl( rectan, azccw, elplsz )
%
%   returns:
%
%      range    the distance of the point from the origin.
%
%               [1,1] = size(range); double = class(range)
%
%               The units associated with `range' are those associated
%               with the input point.
%
%      az       the azimuth of the point.
%
%               [1,1] = size(az); double = class(az)
%
%               This is the angle between the projection onto the XY plane
%               of the vector from the origin to the point and the +X axis of
%               the reference frame. `az' is zero at the +X axis.
%
%               The way azimuth is measured depends on the value of the
%               logical flag `azccw'. See the description of the argument
%               `azccw' for details.
%
%               `az' is output in radians. The range of `az' is [0, 2*pi].
%
%      el       the elevation of the point.
%
%               [1,1] = size(el); double = class(el)
%
%               This is the angle between the vector from the origin to the
%               point and the XY plane. `el' is zero at the XY plane.
%
%               The way elevation is measured depends on the value of
%               the logical flag `elplsz'. See the description of the
%               argument `elplsz' for details.
%
%               `el' is output in radians. The range of `el' is [-pi/2,
%               pi/2].
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
%   1) Create four tables showing a variety of rectangular
%      coordinates and the corresponding range, azimuth and
%      elevation, resulting from the different choices of the `azccw'
%      and `elplsz' flags.
%
%      Corresponding rectangular coordinates and azimuth, elevation
%      and range are listed to three decimal places. Output angles
%      are in degrees.
%
%
%      Example code begins here.
%
%
%      function recazl_ex1()
%
%         %
%         % Local parameters.
%         %
%         NREC =   11;
%
%         %
%         % Define the input rectangular coordinates and the
%         % different choices of the `azccw' and `elplsz' flags.
%         %
%         rectan = [ [0.0,0.0,0.0]',  [1.0,0.0,0.0]',  [0.0,1.0,0.0]',     ...
%                    [0.0,0.0,1.0]',  [-1.0,0.0,0.0]', [0.0,-1.0,0.0]',    ...
%                    [0.0,0.0,-1.0]', [1.0,1.0,0.0]',  [1.0,0.0,1.0]',     ...
%                    [0.0,1.0,1.0]',  [1.0,1.0,1.0]',                      ...
%                     ];
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
%               fprintf( [ '  rect(1)  rect(2)  rect(3)   RANGE      AZ ', ...
%                          '      EL\n' ]                                  )
%               fprintf( [ '  -------  -------  -------  ------- ',        ...
%                          ' -------  -------\n' ]                  )
%
%               %
%               % Do the conversion. Output angles in degrees.
%               %
%               for n=1:NREC
%
%                  [range, az, el] = cspice_recazl( rectan(:, n),          ...
%                                                   azccw(i),              ...
%                                                   elplsz(j) );
%
%                  fprintf( '%9.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',       ...
%                           rectan(1,n), rectan(2,n),    rectan(3,n),      ...
%                           range,        az * cspice_dpr, el * cspice_dpr )
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
%        rect(1)  rect(2)  rect(3)   RANGE      AZ       EL
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          0.000    1.000    0.000    1.000  270.000    0.000
%          0.000    0.000    1.000    1.000    0.000  -90.000
%         -1.000    0.000    0.000    1.000  180.000    0.000
%          0.000   -1.000    0.000    1.000   90.000    0.000
%          0.000    0.000   -1.000    1.000    0.000   90.000
%          1.000    1.000    0.000    1.414  315.000    0.000
%          1.000    0.000    1.000    1.414    0.000  -45.000
%          0.000    1.000    1.000    1.414  270.000  -45.000
%          1.000    1.000    1.000    1.732  315.000  -35.264
%
%      AZCCW = False; ELPLSZ = True
%
%        rect(1)  rect(2)  rect(3)   RANGE      AZ       EL
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          0.000    1.000    0.000    1.000  270.000    0.000
%          0.000    0.000    1.000    1.000    0.000   90.000
%         -1.000    0.000    0.000    1.000  180.000    0.000
%          0.000   -1.000    0.000    1.000   90.000    0.000
%          0.000    0.000   -1.000    1.000    0.000  -90.000
%          1.000    1.000    0.000    1.414  315.000    0.000
%          1.000    0.000    1.000    1.414    0.000   45.000
%          0.000    1.000    1.000    1.414  270.000   45.000
%          1.000    1.000    1.000    1.732  315.000   35.264
%
%      AZCCW = True; ELPLSZ = False
%
%        rect(1)  rect(2)  rect(3)   RANGE      AZ       EL
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          0.000    1.000    0.000    1.000   90.000    0.000
%          0.000    0.000    1.000    1.000    0.000  -90.000
%         -1.000    0.000    0.000    1.000  180.000    0.000
%          0.000   -1.000    0.000    1.000  270.000    0.000
%          0.000    0.000   -1.000    1.000    0.000   90.000
%          1.000    1.000    0.000    1.414   45.000    0.000
%          1.000    0.000    1.000    1.414    0.000  -45.000
%          0.000    1.000    1.000    1.414   90.000  -45.000
%          1.000    1.000    1.000    1.732   45.000  -35.264
%
%      AZCCW = True; ELPLSZ = True
%
%        rect(1)  rect(2)  rect(3)   RANGE      AZ       EL
%        -------  -------  -------  -------  -------  -------
%          0.000    0.000    0.000    0.000    0.000    0.000
%          1.000    0.000    0.000    1.000    0.000    0.000
%          0.000    1.000    0.000    1.000   90.000    0.000
%          0.000    0.000    1.000    1.000    0.000   90.000
%         -1.000    0.000    0.000    1.000  180.000    0.000
%          0.000   -1.000    0.000    1.000  270.000    0.000
%          0.000    0.000   -1.000    1.000    0.000  -90.000
%          1.000    1.000    0.000    1.414   45.000    0.000
%          1.000    0.000    1.000    1.414    0.000   45.000
%          0.000    1.000    1.000    1.414   90.000   45.000
%          1.000    1.000    1.000    1.732   45.000   35.264
%
%
%   2) Compute the apparent azimuth and elevation of Venus as seen
%      from the DSS-14 station.
%
%      Task Description
%      ================
%
%      In this example, we will obtain the apparent position of
%      Venus as seen from the DSS-14 station in the DSS-14 topocentric
%      reference frame. We will use a station frames kernel and
%      transform the resulting rectangular coordinates to azimuth,
%      elevation and range using cspice_azlrec.
%
%      In order to introduce the usage of the logical flags `azccw'
%      and `elplsz', we will request the azimuth to be measured
%      clockwise and the elevation positive towards the +Z
%      axis of the DSS-14_TOPO reference frame.
%
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
%         File name: recazl_ex2.tm
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
%                                             binary PCK
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
%      function recazl_ex2()
%
%         %
%         % Local parameters
%         %
%         META =   'recazl_ex2.tm';
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
%         % Compute the observer-target position.
%         %
%         [ptarg, lt] = cspice_spkpos( target, et, ref, abcorr, obs );
%
%         %
%         % Compute azimuth, elevation and range of Venus
%         % as seen from DSS-14, with azimuth increasing
%         % clockwise and elevation positive towards +Z
%         % axis of the DSS-14_TOPO reference frame
%         %
%         azccw  = false;
%         elplsz = true;
%
%         [r, az, el] = cspice_recazl( ptarg, azccw, elplsz );
%
%         %
%         % Express both angles in degrees.
%         %
%         el =   el * cspice_dpr;
%         az =   az * cspice_dpr;
%
%         %
%         % Display the computed position, the range and
%         % the angles.
%         %
%         fprintf( '\n' )
%         fprintf( 'Target:                %s\n', target )
%         fprintf( 'Observation time:      %s\n', obstim )
%         fprintf( 'Observer center:       %s\n', obs )
%         fprintf( 'Observer frame:        %s\n', ref )
%         fprintf( 'Aberration correction: %s\n', abcorr )
%         fprintf( '\n' )
%         fprintf( 'Observer-target position (km):\n' )
%         fprintf( '%21.8f %20.8f %20.8f\n', ptarg(1), ptarg(2), ptarg(3) )
%         fprintf( 'Light time (s):        %19.8f\n', lt )
%         fprintf( '\n' )
%         fprintf( 'Target azimuth          (deg):  %19.8f\n', az )
%         fprintf( 'Target elevation        (deg):  %19.8f\n', el )
%         fprintf( 'Observer-target distance (km):  %19.8f\n', r )
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
%      Target:                VENUS
%      Observation time:      2003 OCT 13 06:00:00.000000 UTC
%      Observer center:       DSS-14
%      Observer frame:        DSS-14_TOPO
%      Aberration correction: CN+S
%
%      Observer-target position (km):
%          66886767.37916669   146868551.77222887  -185296611.10841593
%      Light time (s):               819.63862811
%
%      Target azimuth          (deg):         294.48543372
%      Target elevation        (deg):         -48.94609726
%      Observer-target distance (km):   245721478.99272084
%
%
%-Particulars
%
%   This routine returns the range, azimuth, and elevation of a point
%   specified in rectangular coordinates.
%
%   The output is defined by the distance from the center of the
%   reference frame (range), the angle from a reference vector
%   (azimuth), and the angle above the XY plane of the reference
%   frame (elevation).
%
%   The way azimuth and elevation are measured depends on the values
%   given by the user to the `azccw' and `elplsz' logical flags. See the
%   descriptions of these input arguments for details.
%
%-Exceptions
%
%   1)  If the X and Y components of `rectan' are both zero, the
%       azimuth is set to zero.
%
%   2)  If `rectan' is the zero vector, azimuth and elevation
%       are both set to zero.
%
%   3)  If any of the input arguments, `rectan', `azccw' or `elplsz',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   4)  If any of the input arguments, `rectan', `azccw' or `elplsz',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
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
%   -Mice Version 1.0.0, 01-NOV-2021 (JDR)
%
%-Index_Entries
%
%   rectangular coordinates to range, az and el
%   rectangular to range, azimuth and elevation
%   convert rectangular coordinates to range, az and el
%   convert rectangular to range, azimuth and elevation
%
%-&
function [range, az, el] = cspice_recazl( rectan, azccw, elplsz )

   switch nargin
      case 3

         rectan = zzmice_dp(rectan);
         azccw  = zzmice_int(azccw);
         elplsz = zzmice_int(elplsz);

      otherwise

         error ( [ 'Usage: [range, az, el] = '                              ...
                   'cspice_recazl( rectan(3), azccw, elplsz )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [range, az, el] = mice('recazl_c', rectan, azccw, elplsz);
   catch spiceerr
      rethrow(spiceerr)
   end
