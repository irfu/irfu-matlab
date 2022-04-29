%-Abstract
%
%   CSPICE_SPKW10 writes an SPK type 10 segment to the file specified by
%   the input `handle'.
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
%      handle   the file handle of an SPK file that has been opened for
%               writing.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      body     the NAIF ID for the body whose states are to be recorded in
%               an SPK file.
%
%               [1,1] = size(body); int32 = class(body)
%
%      center   the NAIF ID for the center of motion associated with `body'.
%
%               [1,1] = size(center); int32 = class(center)
%
%      frame    the reference frame that states are referenced to, for
%               example 'J2000'.
%
%               [1,c1] = size(frame); char = class(frame)
%
%                  or
%
%               [1,1] = size(frame); cell = class(frame)
%
%      first,
%      last     the bounds on the ephemeris times, expressed as
%               seconds past J2000, for which the states can be used to
%               interpolate a state for `body'.
%
%               [1,1] = size(first); double = class(first)
%               [1,1] = size(last); double = class(last)
%
%      segid    the segment identifier.
%
%               [1,c2] = size(segid); char = class(segid)
%
%                  or
%
%               [1,1] = size(segid); cell = class(segid)
%
%               An SPK segment identifier may contain up to 40 characters.
%
%      consts   the geophysical constants needed for evaluation of the two
%               line elements sets.
%
%               [8,1] = size(consts); double = class(consts)
%
%               The order of these constants must be:
%
%                  consts(1) = J2 gravitational harmonic for Earth.
%                  consts(2) = J3 gravitational harmonic for Earth.
%                  consts(3) = J4 gravitational harmonic for Earth.
%
%               These first three constants are dimensionless.
%
%                  consts(4) = KE: Square root of the GM for Earth where
%                              GM is expressed in Earth radii cubed
%                              per minutes squared.
%
%                  consts(5) = QO: High altitude bound for atmospheric
%                              model in km.
%
%                  consts(6) = SO: Low altitude bound for atmospheric
%                              model in km.
%
%                  consts(7) = RE: Equatorial radius of the earth in km.
%
%                  consts(8) = AE: Distance units/earth radius
%                             (normally 1).
%
%               Below are currently recommended values for these
%               items:
%
%                  J2 =    1.082616e-3
%                  J3 =   -2.53881e-6
%                  J4 =   -1.65597e-6
%
%               The next item is the square root of GM for the Earth
%               given in units of earth-radii**1.5/Minute
%
%                  KE =    7.43669161e-2
%
%               The next two items define the top and bottom of the
%               atmospheric drag model used by the type 10 ephemeris
%               type. Don't adjust these unless you understand the full
%               implications of such changes.
%
%                  QO =  120.0e0
%                  SO =   78.0e0
%
%               The ER value is the equatorial radius in km of the Earth
%               as used by NORAD.
%
%                  ER = 6378.135e0
%
%               The value of AE is the number of distance units per
%               Earth radii used by the NORAD state propagation
%               software. The value should be 1 unless you've got a very
%               good understanding of the NORAD routine SGP4 and the
%               affect of changing this value.
%
%                  AE =    1.0e0
%
%      n        the number of "two-line" element sets and epochs to be stored
%               in the segment.
%
%               [1,1] = size(n); int32 = class(n)
%
%      elems    a time-ordered array of two-line elements as supplied in
%               NORAD two-line element files.
%
%               [10*n,1] = size(elems); double = class(elems)
%
%               The i'th set of elements should be stored as shown here:
%
%                  `base' = (i-1)*10
%
%                  elems( base + 1  )  = NDT20 in radians/minute**2
%                  elems( base + 2  )  = NDD60 in radians/minute**3
%                  elems( base + 3  )  = BSTAR
%                  elems( base + 4  )  = INCL  in radians
%                  elems( base + 5  )  = NODE0 in radians
%                  elems( base + 6  )  = ECC
%                  elems( base + 7  )  = OMEGA in radians
%                  elems( base + 8  )  = M0    in radians
%                  elems( base + 9  )  = N0    in radians/minute
%                  elems( base + 10 )  = `epoch' of the elements in seconds
%                                        past ephemeris epoch J2000.
%
%               The meaning of these variables is defined by the
%               format of the two-line element files available from
%               NORAD.
%
%      epochs   an n-dimensional array that contains the epochs (ephemeris
%               seconds past J2000) corresponding to the elements in `elems'.
%
%               [n,1] = size(epochs); double = class(epochs)
%
%               The i'th `epoch' must equal the `epoch' of the i'th element
%               set. Epochs must form a strictly increasing sequence.
%
%   the call:
%
%      cspice_spkw10( handle, body,   center, frame,  first,  last, ...
%                     segid,  consts, n,      elems,  epochs  )
%
%   returns:
%
%      None.
%
%      The routine writes an SPK type 10 segment to the file attached to
%      `handle'.
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
%         File name: spkw10_ex1.tm
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
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                'geophysical.ker'  )
%
%         \begintext
%
%         The geophysical.ker is a PCK file that is provided with the
%         SPICE toolkit under the "/data" directory.
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function spkw10_ex1()
%
%         %
%         % Local parameters.
%         %
%         SPK10 =   'spkw10_ex1.bsp';
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
%                   '1 18123U 87 53  A 87345.43010859  .00000022  '   ...
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
%         cspice_furnsh( 'spkw10_ex1.tm' );
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
%         cspice_spkw10( handle, BODY,  CENTER, FRMNAM, first,             ...
%                        last,   segid, consts, TLESSZ, elems,             ...
%                        epochs );
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
%-Particulars
%
%   This routine writes a type 10 SPK segment to the SPK file open
%   for writing that is attached to `handle'.
%
%   The routine cspice_getelm reads two-line element sets, as those
%   distributed by NORAD, and converts them to the elements in units
%   suitable for use in this routine.
%
%-Exceptions
%
%   1)  If the structure or content of the inputs are invalid, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   2)  If any file access error occurs, the error is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If any of the input arguments, `handle', `body', `center',
%       `frame', `first', `last', `segid', `consts', `n', `elems' or
%       `epochs', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   4)  If any of the input arguments, `handle', `body', `center',
%       `frame', `first', `last', `segid', `consts', `n', `elems' or
%       `epochs', is not of the expected type, or it does not have the
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
%   NAIF_IDS.REQ
%   SPK.REQ
%
%-Literature_References
%
%   [1]  F. Hoots and R. Roehrich, "Spacetrack Report #3: Models for
%        Propagation of the NORAD Element Sets," U.S. Air Force
%        Aerospace Defense Command, Colorado Springs, CO, 1980.
%
%   [2]  F. Hoots, "Spacetrack Report #6: Models for Propagation of
%        Space Command Element Sets,"  U.S. Air Force Aerospace
%        Defense Command, Colorado Springs, CO, 1986.
%
%   [3]  F. Hoots, P. Schumacher and R. Glover, "History of Analytical
%        Orbit Modeling in the U. S. Space Surveillance System,"
%        Journal of Guidance, Control, and Dynamics. 27(2):174-185,
%        2004.
%
%   [4]  D. Vallado, P. Crawford, R. Hujsak and T. Kelso, "Revisiting
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
%   -Mice Version 1.0.0, 05-NOV-2021 (MCS)
%
%-Index_Entries
%
%   Write a type 10 SPK segment
%
%-&
function cspice_spkw10( handle, body,   center, frame, first,  ...
                        last,   segid,  consts, n,     elems,  ...
                        epochs )

   switch nargin
      case 11

         handle = zzmice_int(handle);
         body   = zzmice_int(body);
         center = zzmice_int(center);
         frame  = zzmice_str(frame);
         first  = zzmice_dp(first);
         last   = zzmice_dp(last);
         segid  = zzmice_str(segid);
         consts = zzmice_dp(consts);
         n      = zzmice_int(n);
         elems  = zzmice_dp(elems);
         epochs = zzmice_dp(epochs);

      otherwise

         error ( [ 'Usage: cspice_spkw10( handle, body, center, `frame`, '  ...
                   'first, last, `segid`, consts(8), n, elems(10*n), '      ...
                   'epochs(n) )=' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'spkw10_c', handle, body,   center, frame,  first,              ...
            last,       segid,  consts, n,      elems,  epochs );
   catch spiceerr
      rethrow(spiceerr)
   end
