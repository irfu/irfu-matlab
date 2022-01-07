%-Abstract
%
%   CSPICE_SPKW09 writes a type 9 segment to an SPK file.
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
%      body     the NAIF integer code for an ephemeris object whose state
%               relative to another body is described by the segment to be
%               created.
%
%               [1,1] = size(body); int32 = class(body)
%
%      center   the NAIF integer code for the center of motion of the object
%               identified by body.
%
%               [1,1] = size(center); int32 = class(center)
%
%      frame    the NAIF name for a reference frame relative to which the
%               state information for body is specified.
%
%               [1,c1] = size(frame); char = class(frame)
%
%                  or
%
%               [1,1] = size(frame); cell = class(frame)
%
%      first,
%      last     are, respectively, the start and stop times of the time
%               interval over which the segment defines the state of body.
%
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
%      degree   the degree of the Lagrange polynomials used to interpolate
%               the states.
%
%               [1,1] = size(degree); int32 = class(degree)
%
%               All components of the state vectors are interpolated by
%               polynomials of fixed degree.
%
%      states   contains a time-ordered array of geometric states (x, y, z,
%               dx/dt, dy/dt, dz/dt, in kilometers and kilometers per second)
%               of body relative to center, specified relative to frame.
%
%               [6,n] = size(states); double = class(states)
%
%      epochs   an array of epochs corresponding to the members of the state
%               array.
%
%               [1,n] = size(epochs); double = class(epochs)
%
%               The epochs are specified as seconds past J2000, TDB.
%
%   the call:
%
%      cspice_spkw09( handle, body,   center, frame,  first,  last, ...
%                     segid,  degree, states, epochs  )
%
%   returns:
%
%      None.
%
%      See -Particulars for a description of the effect of this routine.
%
%-Parameters
%
%   MAXDEG      is the maximum allowed degree of the interpolating
%               polynomial.
%
%               The current value of MAXDEG is 15.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Suppose that you have a time-ordered array of geometric states
%      of a new object that follows Phobos, with a delay of 1 hour,
%      in its orbit around Mars and are prepared to produce a segment
%      of type 09 in an SPK file. Create a new SPK file with this
%      segment. Use an existing SPK to create the input data for the
%      SPK segment.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: spkw09_ex1.tm
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
%            mar097.bsp                       Mars satellite ephemeris
%            naif0012.tls                     Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'mar097.bsp',
%                                'naif0012.tls' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function spkw09_ex1()
%
%         %
%         % Local parameters.
%         %
%         SPKNAM = 'spkw09_ex1.bsp';
%         DEGREE = 3;
%         MARS   = 499;
%         NEPOCS = 800;
%         NOBJ   = 403;
%
%         %
%         % Local variables.
%         %
%         epochs = zeros(NEPOCS,1);
%         states = zeros(6,NEPOCS);
%
%         %
%         % Load the input SPK file.
%         %
%         cspice_furnsh( 'spkw09_ex1.tm' );
%
%         %
%         % Convert the input UTC to ephemeris time
%         %
%         [et] = cspice_str2et( '2018 Apr 03 08:35' );
%
%         %
%         % Create the time-ordered array of geometric states,
%         % at unequal time steps.
%         %
%         time  = et;
%         step  = 60.0;
%         delta = 10.0;
%
%         for i=1:NEPOCS
%
%            [states(:,i), lt] = cspice_spkezr( 'PHOBOS', time,  'J2000', ...
%                                               'NONE',   'MARS'          );
%
%            epochs(i) = time + 3600.0;
%            time      =   time + step                             ...
%                        + sin( cspice_halfpi * i / 2.0 ) * delta;
%
%         end
%
%         %
%         % Open a new SPK file, with 5000 characters reserved
%         % for comments.
%         %
%         ifname   = 'Test SPK type 9 internal filename.';
%         [handle] = cspice_spkopn( SPKNAM, ifname, 5000 );
%
%         %
%         % Create a segment identifier.
%         %
%         segid = 'MY_SAMPLE_SPK_TYPE_9_SEGMENT';
%
%         %
%         % Write the segment.
%         %
%         cspice_spkw09( handle,    NOBJ,           MARS,  'J2000', ...
%                        epochs(1), epochs(NEPOCS), segid, DEGREE,  ...
%                        states,    epochs                          );
%
%         %
%         % Close the new SPK file.
%         %
%         cspice_spkcls( handle );
%
%         %
%         % Compute the state of Phobos as seen from Mars,
%         % 12 hours after the input UTC time.
%         %
%         et = et + 43200.0;
%         [state, lt] = cspice_spkezr( 'PHOBOS', et,    'J2000', ...
%                                      'NONE',   'MARS'          );
%
%         fprintf( 'Phobos as seen from Mars at t0:\n' )
%         fprintf( '   Epoch       (s): %19.5f\n', et )
%         fprintf( '   Position   (km): %13.5f %13.5f %13.5f\n', ...
%                                   state(1), state(2), state(3) )
%         fprintf( '   Velocity (km/s): %13.5f %13.5f %13.5f\n', ...
%                                   state(4), state(5), state(6) )
%         fprintf( '\n' )
%
%         %
%         % Load the newly created kernel, and compute the state
%         % of the new object as seen from Mars, 13 hours after
%         % the input UTC time.
%         %
%         cspice_furnsh( SPKNAM );
%         et = et + 3600.0;
%
%         [state, lt] = cspice_spkezr( '403', et, 'J2000', 'NONE', 'MARS' );
%
%         fprintf( 'Object 403 as seen from Mars at t0 + 1h:\n' )
%         fprintf( '   Epoch       (s): %19.5f\n', et )
%         fprintf( '   Position   (km): %13.5f %13.5f %13.5f\n', ...
%                                   state(1), state(2), state(3) )
%         fprintf( '   Velocity (km/s): %13.5f %13.5f %13.5f\n', ...
%                                   state(4), state(5), state(6) )
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
%      Phobos as seen from Mars at t0:
%         Epoch       (s):     576059769.18566
%         Position   (km):   -7327.26277    2414.32655    5207.10638
%         Velocity (km/s):      -0.94289      -1.89473      -0.39671
%
%      Object 403 as seen from Mars at t0 + 1h:
%         Epoch       (s):     576063369.18566
%         Position   (km):   -7327.26277    2414.32655    5207.10638
%         Velocity (km/s):      -0.94289      -1.89473      -0.39671
%
%
%      Note that after run completion, a new SPK file exists in
%      the output directory.
%
%-Particulars
%
%   This routine writes an SPK type 09 data segment to the open SPK
%   file according to the format described in the type 09 section of
%   the SPK Required Reading. The SPK file must have been opened with
%   write access.
%
%-Exceptions
%
%   If any of the following exceptions occur, this routine will return
%   without creating a new segment.
%
%   1)  If `frame' is not a recognized name, the error
%       SPICE(INVALIDREFFRAME) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the last non-blank character of `segid' occurs past index 40,
%       the error SPICE(SEGIDTOOLONG) is signaled by a routine in the
%       call tree of this routine.
%
%   3)  If `segid' contains any nonprintable characters, the error
%       SPICE(NONPRINTABLECHARS) is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If `degree' is not at least 1 or is greater than MAXDEG, the
%       error SPICE(INVALIDDEGREE) is signaled by a routine in the
%       call tree of this routine.
%
%   5)  If the number of states `n' is not at least degree+1, the
%       error SPICE(TOOFEWSTATES) is signaled by a routine in the call
%       tree of this routine.
%
%   6)  If `first' is greater than or equal to `last', the error
%       SPICE(BADDESCRTIMES) is signaled by a routine in the call tree
%       of this routine.
%
%   7)  If the elements of the array `epochs' are not in strictly
%       increasing order, the error SPICE(TIMESOUTOFORDER) is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If the first epoch, epochs(1), is greater than `first', the
%       error SPICE(BADDESCRTIMES) is signaled by a routine in the
%       call tree of this routine.
%
%   9)  If the last epoch, epochs(n), where `n' is the number of
%       states, is less than `last', the error SPICE(BADDESCRTIMES) is
%       signaled by a routine in the call tree of this routine.
%
%   10) If any of the input arguments, `handle', `body', `center',
%       `frame', `first', `last', `segid', `degree', `states' or
%       `epochs', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   11) If any of the input arguments, `handle', `body', `center',
%       `frame', `first', `last', `segid', `degree', `states' or
%       `epochs', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%   12) If the input vector arguments `states' and `epochs' do not
%       have the same dimension (n), an error is signaled by the Mice
%       interface.
%
%-Files
%
%   A new type 9 SPK segment is written to the SPK file attached
%   to `handle'.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   NAIF_IDS.REQ
%   SPC.REQ
%   SPK.REQ
%   TIME.REQ
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
%   -Mice Version 1.0.0, 26-NOV-2021 (JDR)
%
%-Index_Entries
%
%   write SPK type_9 ephemeris data segment
%
%-&
function cspice_spkw09( handle, body,  center, frame,  first, ...
                        last,   segid, degree, states, epochs )

   switch nargin
      case 10

         handle = zzmice_int(handle);
         body = zzmice_int(body);
         center = zzmice_int(center);
         frame = zzmice_str(frame);
         first = zzmice_dp(first);
         last = zzmice_dp(last);
         segid = zzmice_str(segid);
         degree = zzmice_int(degree);
         states = zzmice_dp(states);
         epochs = zzmice_dp(epochs);

      otherwise

         error ( ['Usage: cspice_spkw09( handle, body, center, `frame`,' ...
                  ' first, last, `segid`, degree, states(6,n), epochs(n) )'] )

   end

   %
   % Call the MEX library.
   %
   try
      mice('spkw09_c', handle, body,  center, frame,  first,  ...
                       last,   segid, degree, states, epochs);
   catch spiceerr
      rethrow(spiceerr)
   end
