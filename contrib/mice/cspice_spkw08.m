%-Abstract
%
%   CSPICE_SPKW08 writes a type 8 segment to an SPK file.
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
%      handle   the file handle of an SPK file open for writing.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      body     the SPICE ID code for an ephemeris object
%               whose state relative to another body is described
%               by the segment to be created.
%
%               [1,1] = size(body); int32 = class(body)
%
%      center   the SPICE ID code for the center of motion
%               of the object identified by body.
%
%               [1,1] = size(center); int32 = class(center)
%
%      frame    the name for a reference frame relative to which the state
%               information for body is specified.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%      first,
%      last     are, respectively, the start and stop times of
%               the time interval over which the segment defines
%               the state of body.
%
%               [1,1] = size(first); double = class(first)
%               [1,1] = size(last);  double = class(last)
%
%      segid    is the segment identifier. An SPK segment
%               identifier may contain up to 40 characters.
%
%      degree   the degree of the Lagrange polynomials used to
%               interpolate the states. All components of the
%               state vectors are interpolated by polynomials of
%               fixed degree.
%
%               [1,1] = size(degree); int32 = class(degree)
%
%      states   contains a time-ordered array of geometric states
%               ( x, y, z, dx/dt, dy/dt, dz/dt, in kilometers and
%               kilometers per second ) of body relative to center,
%               specified relative to frame.
%
%               [6,m] = size(states); double = class(states)
%
%      begtim   the epoch corresponding to the first state in
%               the state array. Because extra states are needed
%               at the beginning and end of the segment in order
%               for the interpolation method to work, begtim will
%               normally precede first.
%
%               [1,1] = size(begtim); double = class(begtim)
%
%      step     the time step separating the epochs of adjacent
%               states in the input state array. step is specified
%               in TDB seconds.
%
%               [1,1] = size(step); double = class(step)
%
%   the call:
%
%      cspice_spkw08( handle, body,  center, frame,  first,  ...
%                     last,   segid, degree, states, begtim, ...
%                     step )
%
%   returns:
%
%   The routine writes to the SPK file referred to by `handle' a type 8 SPK
%   segment containing the data listed in `states'.
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
%   1) This example demonstrates how to create an SPK type 8 kernel
%      containing only one segment, given a time-ordered set of
%      discrete states and epochs.
%
%
%      Example code begins here.
%
%
%      function spkw08_ex1()
%
%         %
%         % Define the segment identifier parameters.
%         %
%         BODY       = 3;
%         CENTER     = 10;
%         REF        = 'J2000';
%         POLY_DEG   = 3;
%         SPK8       = 'spkw08_ex1.bsp';
%         N_DISCRETE = 9;
%
%         %
%         % A set of epochs.
%         %
%         DISCRETEEPOCHS = (1:9)*100;
%
%         %
%         % An array of discrete states to write to the SPK segment.
%         %
%         base = [ (1:6)*100 ]';
%
%         %
%         % Create the 6xN array of states.
%         %
%         DISCRETESTATES = [(base+1), (base+2), (base+3), ...
%                           (base+4), (base+5), (base+6), ...
%                           (base+7), (base+8), (base+9) ];
%
%         %
%         % Create a segment identifier.
%         %
%         segid = 'SPK type 8 test segment';
%
%         %
%         % Open a new SPK file.
%         %
%         handle = cspice_spkopn( SPK8, segid, 4 );
%
%         step   = DISCRETEEPOCHS(2) - DISCRETEEPOCHS(1);
%
%         %
%         % Create a type 8 segment.
%         %
%         cspice_spkw08( handle,                       ...
%                        BODY,                         ...
%                        CENTER,                       ...
%                        REF,                          ...
%                        DISCRETEEPOCHS(1),            ...
%                        DISCRETEEPOCHS(N_DISCRETE),   ...
%                        segid,                        ...
%                        POLY_DEG,                     ...
%                        DISCRETESTATES,               ...
%                        DISCRETEEPOCHS(1),            ...
%                        step )
%
%         %
%         % Close the SPK file.
%         %
%         cspice_spkcls( handle )
%
%
%      When this program is executed, no output is presented on
%      screen. After run completion, a new SPK type 8 exists in
%      the output directory.
%
%-Particulars
%
%   This routine writes an SPK type 08 data segment to the open SPK
%   file according to the format described in the type 08 section of
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
%   6)  If `first' is greater than `last', the error SPICE(BADDESCRTIMES)
%       is signaled by a routine in the call tree of this routine.
%
%   7)  If `step' is non-positive, the error SPICE(INVALIDSTEPSIZE) is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If the start time of the first record exceeds the descriptor
%       begin time by more than a computed tolerance, or if the end
%       time of the last record precedes the descriptor end time by
%       more than a computed tolerance, the error SPICE(COVERAGEGAP)
%       is signaled by a routine in the call tree of this routine. See
%       the -Parameters section above for a description of the
%       tolerance.
%
%   9)  If any of the input arguments, `handle', `body', `center',
%       `frame', `first', `last', `segid', `degree', `states',
%       `begtim' or `step', is undefined, an error is signaled by the
%       Matlab error handling system.
%
%   10) If any of the input arguments, `handle', `body', `center',
%       `frame', `first', `last', `segid', `degree', `states',
%       `begtim' or `step', is not of the expected type, or it does
%       not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   A new type 8 SPK segment is written to the SPK file attached
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Changed the input argument name "epoch1" to "begtim" for
%       consistency with other routines.
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section. Updated the Required Reading section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 23-MAY-2012 (EDW)
%
%-Index_Entries
%
%   write SPK type_8 ephemeris data segment
%
%-&

function cspice_spkw08( handle, ...
                        body,   ...
                        center, ...
                        frame,  ...
                        first,  ...
                        last,   ...
                        segid,  ...
                        degree, ...
                        states, ...
                        begtim, ...
                        step )

   switch nargin
      case 11

         handle  = zzmice_int(handle);
         body    = zzmice_int(body);
         center  = zzmice_int( center);
         frame   = zzmice_str(frame);
         first   = zzmice_dp(first);
         last    = zzmice_dp(last);
         segid   = zzmice_str(segid);
         degree  = zzmice_int(degree);
         states  = zzmice_dp(states);
         begtim  = zzmice_dp(begtim);
         step    = zzmice_dp(step);

      otherwise

         error ( [ 'Usage: '                ...
                   'cspice_spkw08( handle, '...
                                  'body, '  ...
                                  'center, '...
                                  'frame, ' ...
                                  'first, ' ...
                                  'last, '  ...
                                  'segid, ' ...
                                  'degree, '...
                                  'states, '...
                                  'begtim, '...
                                  'step )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'spkw08_c', handle, ...
                        body,   ...
                        center, ...
                        frame,  ...
                        first,  ...
                        last,   ...
                        segid,  ...
                        degree, ...
                        states, ...
                        begtim, ...
                        step )

   catch spiceerr
      rethrow(spiceerr)
   end
