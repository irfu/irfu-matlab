%-Abstract
%
%   CSPICE_GETELM parses the "lines" of a two-line element set, returning the
%   elements in units suitable for use in SPICE software.
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
%      frstyr   the first year possible for two-line elements.
%
%               [1,1] = size(frstyr); int32 = class(frstyr)
%
%               Since two-line elements allow only two digits for the year,
%               some conventions must be followed concerning which century
%               the two digits refer to. `frstyr' is the year of the earliest
%               representable elements. The two-digit year is mapped to the
%               year in the interval from `frstyr' to frstyr + 99 that has
%               the same last two digits as the two digit year in the element
%               set. For example if `frstyr' is set to 1960 then the two
%               digit years are mapped as shown in the table below:
%
%                  Two-line         Maps to
%                  element year
%                     00            2000
%                     01            2001
%                     02            2002
%                      .              .
%                      .              .
%                      .              .
%                     58            2058
%                     59            2059
%                    -------------------
%                     60            1960
%                     61            1961
%                     62            1962
%                      .              .
%                      .              .
%                      .              .
%                     99            1999
%
%               Note that if Space Command should decide to represent
%               years in 21st century as 100 + the last two digits
%               of the year (for example: 2015 is represented as 115)
%               instead of simply dropping the first two digits of
%               the year, this routine will correctly map the year
%               as long as you set `frstyr' to some value between 1900
%               and 1999.
%
%      lines    a pair of lines of text that comprise a Space command
%               "two-line element" set.
%
%               [2,c1] = size(lines); char = class(lines)
%
%                  or
%
%               [2,1] = size(lines); cell = class(lines)
%
%               These text lines should be the same as they are presented
%               in the two-line element files available from Space Command
%               (formerly NORAD). See -Particulars for a detailed description
%               of the format.
%
%   the call:
%
%      [epoch, elems] = cspice_getelm( frstyr, lines )
%
%   returns:
%
%      epoch    the epoch of the two-line elements supplied via the input
%               array `lines'.
%
%               [1,1] = size(epoch); double = class(epoch)
%
%               `epoch' is returned in TDB seconds past J2000.
%
%      elems    an array containing the elements from the two-line set
%               supplied via the array `lines'.
%
%               [10,1] = size(elems); double = class(elems)
%
%               The elements are in units suitable for use by the SPICE
%               routines EV2LIN and cspice_spkw10.
%
%               Also note that the elements XNDD6O and BSTAR
%               incorporate the exponential factor present in the
%               input two-line elements in `lines'. (See -Particulars
%               below).
%
%                  elems(  1 ) = NDT20 in radians/minute**2
%                  elems(  2 ) = NDD60 in radians/minute**3
%                  elems(  3 ) = BSTAR
%                  elems(  4 ) = INCL  in radians
%                  elems(  5 ) = NODE0 in radians
%                  elems(  6 ) = ECC
%                  elems(  7 ) = OMEGA in radians
%                  elems(  8 ) = M0    in radians
%                  elems(  9 ) = N0    in radians/minute
%                  elems( 10 ) = `epoch' of the elements in seconds
%                                past ephemeris epoch J2000.
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
%   1) Suppose that you have collected the two-line element data
%      for a spacecraft with NORAD ID 18123. The following example
%      code demonstrates how you could go about creating a type 10
%      SPK segment.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: getelm_ex1.tm
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
%            File name           Contents
%            ---------           ------------------------------------
%            naif0012.tls        Leapseconds
%            geophysical.ker     geophysical constants for evaluation
%                                of two-line element sets.
%
%         The geophysical.ker is a PCK file that is provided with the
%         Mice toolkit under the "/data" directory.
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                'geophysical.ker'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function getelm_ex1()
%
%         %
%         % Local parameters.
%         %
%         SPK10 =   'getelm_ex1.bsp';
%
%         %
%         % The SPK type 10 segment will contain 18 two-line
%         % elements sets for the `norad' spacecraft 18123 with
%         % respect to the Earth (ID 399) in the J2000 reference
%         % frame.
%         %
%         % As stated in the naif_ids required reading, for Earth
%         % orbiting spacecraft lacking a DSN identification code,
%         % the NAIF ID is derived from the tracking ID assigned to
%         % it by `norad' via:
%         %
%         %    NAIF ID = -100000 - norad ID code
%         %
%         TLESSZ =   9;
%         BODY   =   -118123;
%         CENTER =   399;
%         FRMNAM =   'J2000';
%
%         %
%         % Local variables.
%         %
%         consts = zeros( 8,1         );
%         elems  = zeros( 10*TLESSZ,1 );
%         epochs = zeros( TLESSZ,1    );
%
%         %
%         % These are the variables that will hold the constants
%         % required by SPK type 10. These constants are available
%         % from the loaded PCK file, which provides the actual
%         % values and units as used by `norad' propagation model.
%         %
%         %    Constant   Meaning
%         %    --------   ------------------------------------------
%         %    J2         J2 gravitational harmonic for Earth.
%         %    J3         J3 gravitational harmonic for Earth.
%         %    J4         J4 gravitational harmonic for Earth.
%         %    KE         Square root of the GM for Earth.
%         %    QO         High altitude bound for atmospheric model.
%         %    SO         Low altitude bound for atmospheric model.
%         %    ER         Equatorial radius of the Earth.
%         %    AE         Distance units/earth radius.
%         %
%         noadpn = {'J2','J3','J4','KE','QO','SO','ER','AE'};
%
%         %
%         % Define the Two-Line Element sets.
%         %
%         tle    = [ '1 18123U 87 53  A 87324.61041692 -.00000023  '  ...
%                                         '00000-0 -75103-5 0 00675',
%                      '2 18123  98.8296 152.0074 0014950 168.7820 '  ...
%                                       '191.3688 14.12912554 21686',
%                    '1 18123U 87 53  A 87326.73487726  .00000045  '  ...
%                                         '00000-0  28709-4 0 00684',
%                      '2 18123  98.8335 154.1103 0015643 163.5445 '  ...
%                                       '196.6235 14.12912902 21988',
%                    '1 18123U 87 53  A 87331.40868801  .00000104  '  ...
%                                         '00000-0  60183-4 0 00690',
%                      '2 18123  98.8311 158.7160 0015481 149.9848 '  ...
%                                       '210.2220 14.12914624 22644',
%                    '1 18123U 87 53  A 87334.24129978  .00000086  '  ...
%                                         '00000-0  51111-4 0 00702',
%                      '2 18123  98.8296 161.5054 0015372 142.4159 '  ...
%                                       '217.8089 14.12914879 23045',
%                    '1 18123U 87 53  A 87336.93227900 -.00000107  '  ...
%                                         '00000-0 -52860-4 0 00713',
%                      '2 18123  98.8317 164.1627 0014570 135.9191 '  ...
%                                       '224.2321 14.12910572 23425',
%                    '1 18123U 87 53  A 87337.28635487  .00000173  '  ...
%                                         '00000-0  10226-3 0 00726',
%                      '2 18123  98.8284 164.5113 0015289 133.5979 '  ...
%                                       '226.6438 14.12916140 23475',
%                    '1 18123U 87 53  A 87339.05673569  .00000079  '  ...
%                                         '00000-0  47069-4 0 00738',
%                      '2 18123  98.8288 166.2585 0015281 127.9985 '  ...
%                                       '232.2567 14.12916010 24908',
%                    '1 18123U 87 53  A 87345.43010859  .00000022  '   ...
%                                         '00000-0  16481-4 0 00758',
%                      '2 18123  98.8241 172.5226 0015362 109.1515 '  ...
%                                       '251.1323 14.12915487 24626',
%                    '1 18123U 87 53  A 87349.04167543  .00000042  '  ...
%                                         '00000-0  27370-4 0 00764',
%                      '2 18123  98.8301 176.1010 0015565 100.0881 '  ...
%                                       '260.2047 14.12916361 25138' ];
%
%         %
%         % Load the PCK file that provides the geophysical
%         % constants required for the evaluation of the two-line
%         % elements sets. Load also an LSK, as it is required by
%         % cspice_getelm to perform time conversions. Use a metakernel for
%         % convenience.
%         %
%         cspice_furnsh( 'getelm_ex1.tm' );
%
%         %
%         % Retrieve the data from the kernel, and place it on
%         % the `consts' array.
%         %
%         for i=1:8
%
%            [consts(i)] = cspice_bodvcd( CENTER, noadpn(i), 1 );
%
%         end
%
%         %
%         % Convert the Two Line Elements lines to the
%         % element sets.
%         %
%         j = 0;
%         for i=1:TLESSZ
%
%            [tmpEpochs, tmpElems] = cspice_getelm( 1950,             ...
%                                                   tle(1+2*(i-1):i*2,:) );
%
%            epochs(i)             = tmpEpochs;
%            elems(1+j*10:10+j*10) = tmpElems;
%            j                     = j + 1;
%
%         end
%
%         %
%         % Define the beginning and end of the segment to be
%         % -/+ 12 hours from the first and last epochs,
%         % respectively.
%         %
%         first = epochs(1) - 0.5 * cspice_spd;
%         last  = epochs(TLESSZ-1) + 0.5 * cspice_spd;
%
%         %
%         % `ncomch' is the number of characters to reserve for the
%         % kernel's comment area. This example doesn't write
%         % comments, so set to zero.
%         %
%         ncomch = 0;
%
%         %
%         % Internal file name and segment ID.
%         %
%         ifname = 'Test for type 10 SPK internal file name';
%         segid  = 'SPK type 10 test segment';
%
%         %
%         % Open a new SPK file.
%         %
%         [handle] = cspice_spkopn( SPK10, ifname, ncomch );
%
%         %
%         % Now add the segment.
%         %
%         cspice_spkw10( handle, BODY,   CENTER, FRMNAM, first,  last,     ...
%                        segid,  consts, TLESSZ, elems,  epochs        );
%
%         %
%         % Close the SPK file.
%         %
%         cspice_spkcls( handle );
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program is executed, no output is presented on
%      screen. After run completion, a new SPK type 10 exists in
%      the output directory.
%
%
%   2) Suppose you have a set of two-line elements for the LUME 1
%      cubesat. This example shows how you can use this routine
%      together with the routine cspice_evsgp4 to propagate a state to an
%      epoch of interest.
%
%      Use the meta-kernel from the previous example to load the
%      required SPICE kernels.
%
%
%      Example code begins here.
%
%
%      function getelm_ex2()
%
%         %
%         % Local parameters.
%         %
%         TIMSTR =   '2020-05-26 02:25:00';
%
%         %
%         % The lume-1 cubesat is an Earth orbiting object; set
%         % the center ID to the Earth ID.
%         %
%         CENTER =   399;
%
%         %
%         % Local variables.
%         %
%
%         geophs = zeros(8,1);
%
%         %
%         % These are the variables that will hold the constants
%         % required by cspice_evsgp4. These constants are available from
%         % the loaded PCK file, which provides the actual values
%         % and units as used by NORAD propagation model.
%         %
%         %    Constant   Meaning
%         %    --------   ------------------------------------------
%         %    J2         J2 gravitational harmonic for Earth.
%         %    J3         J3 gravitational harmonic for Earth.
%         %    J4         J4 gravitational harmonic for Earth.
%         %    KE         Square root of the GM for Earth.
%         %    QO         High altitude bound for atmospheric model.
%         %    SO         Low altitude bound for atmospheric model.
%         %    ER         Equatorial radius of the Earth.
%         %    AE         Distance units/earth radius.
%         %
%         noadpn = {'J2','J3','J4','KE','QO','SO','ER','AE'};
%
%         %
%         % Define the Two-Line Element set for LUME-1.
%         %
%         tle = [ '1 43908U 18111AJ  20146.60805006  .00000806'            ...
%                                  '  00000-0  34965-4 0  9999',
%                 '2 43908  97.2676  47.2136 0020001 220.6050 '            ...
%                                  '139.3698 15.24999521 78544' ];
%
%         %
%         % Load the PCK file that provides the geophysical
%         % constants required for the evaluation of the two-line
%         % elements sets. Load also an LSK, as it is required by
%         % cspice_getelm to perform time conversions. Use a metakernel for
%         % convenience.
%         %
%         cspice_furnsh( 'getelm_ex1.tm' );
%
%         %
%         % Retrieve the data from the kernel, and place it on
%         % the `geophs' array.
%         %
%         for i=1:8
%
%            [geophs(i)] = cspice_bodvcd( CENTER, noadpn(i), 1 );
%
%         end
%
%         %
%         % Convert the Two Line Elements lines to the element sets.
%         % Set the lower bound for the years to be the beginning
%         % of the space age.
%         %
%         [epoch, elems] = cspice_getelm( 1957, tle );
%
%         %
%         % Now propagate the state using cspice_evsgp4 to the epoch of
%         % interest.
%         %
%         [et]    = cspice_str2et( TIMSTR );
%         [state] = cspice_evsgp4( et, geophs, elems );
%
%         %
%         % Display the results.
%         %
%         fprintf( 'Epoch   : %s\n', TIMSTR )
%         fprintf( 'Position: %15.8f %15.8f %15.8f\n',                     ...
%                         state(1), state(2), state(3) )
%         fprintf( 'Velocity: %15.8f %15.8f %15.8f\n',                     ...
%                         state(4), state(5), state(6) )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a PC/Linux/Matlab9.x/32-bit
%      platform, the output was:
%
%
%      Epoch   : 2020-05-26 02:25:00
%      Position:  -4644.60403398  -5038.95025539   -337.27141116
%      Velocity:     -0.45719025      0.92884817     -7.55917355
%
%
%-Particulars
%
%         This routine parses a Space Command Two-line element set and
%         returns the orbital elements properly scaled and in units
%         suitable for use by other SPICE software. Input elements
%         are provided in two-lines in accordance with the format
%         required by the two-line element sets available from Space
%         Command (formerly NORAD). See [1] and [2] for details.
%
%         Each of these lines is 69 characters long. The following table
%         defines each of the individual fields for lines 1 and 2.
%
%         Line  Column  Type  Description
%         ----  ------  ----  ------------------------------------------
%           1      01     N   Line number of Element Data (always 1)
%           1    03-07    N   Satellite number (NORAD catalog number)
%           1      08     A   Classification (U:Unclassified; S:Secret)
%           1    10-11    N   International designator (last two digits
%                             of launch year).
%           1    12-14    N   International designator (launch number of
%                             the year).
%           1    15-17    A   International designator (piece of the
%                             launch)
%           1    19-20    N   Epoch year (last two digits of year).
%           1    21-32    N   Epoch (day of the year and portion of the
%                             day)
%           1    34-43    N   NDT20: first time derivative of Mean
%                                    Motion
%           1    45-52    N   NDD60: Second time derivative of Mean
%                                    Motion (decimal point assumed)
%           1    54-61    N   BSTAR drag term (decimal point assumed)
%           1      63     N   Ephemeris type
%           1    65-68    N   Element number
%           1      69     N   Checksum.
%
%           2      01     N   Line number of Element Data (always 2)
%           2    03-07    N   Satellite number (must be the same as in
%                             line 1)
%           2    09-16    N   INCL: Inclination, in degrees
%           2    18-25    N   NODE0: Right Ascension of the Ascending
%                                    Node, in degrees
%           2    27-33    N   ECC: Eccentricity (decimal point assumed)
%           2    35-42    N   OMEGA: Argument of Perigee, in degrees
%           2    44-51    N   M0: Mean Anomaly, in degrees
%           2    53-63    N   N0: Mean Motion (revolutions per day)
%           2    64-68    N   Revolution number at epoch
%           2      69     N   Checksum
%
%         The column type A indicates "characters A-Z", the type N means
%         "numeric."
%
%         Column refers to the substring within the line, e.g.
%
%
%   1 22076U 92052A   97173.53461370 -.00000038  00000-0  10000-3 0   594
%   2 22076  66.0378 163.4372 0008359 278.7732  81.2337 12.80930736227550
%   ^
%   123456789012345678901234567890123456789012345678901234567890123456789
%            1         2         3         4         5         6
%
%
%         In this example, the satellite number (column 03-07) is 22076.
%
%         The "raw" elements used by this routine in the first lines are
%         described in detail below as in several instances exponents and
%         decimal points are implied. Note that the input units are
%         degrees, degrees/day**n and revolutions/day.
%
%         The epoch (column 19-32; line 1) has a format NNNNN.NNNNNNNN,
%         where:
%
%                    Fraction
%                DOY  of day
%                --- --------
%              NNNNN.NNNNNNNN
%              --
%            Year
%
%         An epoch of 00001.00000000 corresponds to 00:00:00 UTC on
%         2000 January 01.
%
%         The first derivative of Mean Motion (column 34-43, line 1), has
%         a format +.NNNNNNNN, where the first character could be either
%         a plus sign, a minus sign or a space.
%
%         The second derivative of Mean Motion (column 45-52, line 1) and
%         the BSTAR drag term (see [1] for details -- column 54-61, line
%         1) have a format +NNNNN-N, where the first character could be
%         either a plus sign, a minus sign or a space, the decimal point
%         is assumed, and the exponent is marked by the sign (+/-) at
%         character 6 (column 51 and 60 for the second derivative and
%         BSTAR drag term respectively).
%
%         The "raw" elements in the second line consists primarily of
%         mean elements calculated using the sgp4/sdp4 orbital model (See
%         [1]). The Inclination, the Right Ascension of the Ascending
%         Node, the Argument of Perigee and the Mean Anomaly have units
%         of degrees and can range from 0 up to 360 degrees, except for
%         the Inclination that ranges from 0 to 180 degrees. The
%         Eccentricity value is provided with an assumed leading decimal
%         point. For example, a value of 9790714 corresponds to an
%         eccentricity of 0.9790714. The Mean motion is measured in
%         revolutions per day and its format is NN.NNNNNNN.
%
%         This routine extracts these values, "inserts" the implied
%         decimal points and exponents and then converts the inputs
%         to units of radians, radians/minute, radians/minute**2, and
%         radians/minute**3
%
%-Exceptions
%
%   1)  If an error occurs while trying to parse the two-line element
%       set, the error SPICE(BADTLE) is signaled by a routine in the
%       call tree of this routine and a description of the detected
%       issue in the "two-line element" set is reported on the long
%       error message.
%
%   2)  If any of the input arguments, `frstyr' or `lines', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `frstyr' or `lines', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   You must have loaded a SPICE leapseconds kernel into the
%   kernel pool prior to calling this routine.
%
%-Restrictions
%
%   1)  The format of the two-line elements suffer from a "millennium"
%       problem --- only two digits are used for the year of the
%       elements. It is not clear how Space Command will deal with
%       this problem. NAIF hopes that by adjusting the input `frstyr'
%       you should be able to use this routine well into the 21st
%       century.
%
%       The approach taken to mapping the two-digit year to the
%       full year is given by the code below. Here, YR is the
%       integer obtained by parsing the two-digit year from the first
%       line of the elements.
%
%          begyr = (frstyr/100)*100;
%          year  = begyr + yr;
%
%          if ( year < frstyr )
%
%             year = year + 100;
%
%          end
%
%       This mapping will be changed if future two-line element
%       representations make this method of computing the full year
%       inaccurate.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  F. Hoots and R. Roehrich, "Spacetrack Report #3: Models for
%        Propagation of the NORAD Element Sets," U.S. Air Force
%        Aerospace Defense Command, Colorado Springs, CO, 1980.
%
%   [2]  "SDC/SCC Two Card Element Set - Transmission Format,"
%        ADCOM/DO Form 12.
%
%   [3]  F. Hoots, "Spacetrack Report #6: Models for Propagation of
%        Space Command Element Sets,"  U.S. Air Force Aerospace
%        Defense Command, Colorado Springs, CO, 1986.
%
%   [4]  F. Hoots, P. Schumacher and R. Glover, "History of Analytical
%        Orbit Modeling in the U. S. Space Surveillance System,"
%        Journal of Guidance, Control, and Dynamics. 27(2):174-185,
%        2004.
%
%   [5]  D. Vallado, P. Crawford, R. Hujsak and T. Kelso, "Revisiting
%        Spacetrack Report #3," paper AIAA 2006-6753 presented at the
%        AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO.,
%        August 21-24, 2006.
%
%-Author_and_Institution
%
%   M. Costa Sitja      (JPL)
%
%-Version
%
%   -Mice Version 1.0.0, 06-NOV-2021 (MCS)
%
%-Index_Entries
%
%   Parse two-line elements
%
%-&
function [epoch, elems] = cspice_getelm( frstyr, lines )

   switch nargin
      case 2

         frstyr = zzmice_int(frstyr);
         lines  = zzmice_str(lines);

      otherwise

         error ( ['Usage: [epoch, elems(10) ] = ' ...
                 'cspice_getelm( frstyr, `lines(2)` )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [epoch, elems] = mice( 'getelm_c', frstyr, lines );
   catch spiceerr
      rethrow(spiceerr)
   end
