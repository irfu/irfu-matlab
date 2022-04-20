%-Abstract
%
%   CSPICE_CKCOV finds the coverage window for a specified object in a
%   specified CK file.
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
%      ckfnm    the name(s) of SPICE CK files.
%
%               [n,c1] = size(ckfnm); char = class(ck)
%
%                 or
%
%               [1,n] = size(ckfnm); cell = class(ck)
%
%      idcode   the CK ID code of an object, normally a spacecraft
%               structure or instrument, for which pointing data are expected
%               to exist in the specified CK file.
%
%               [1,1] = size(idcode); int32 = class(idcode)
%
%      needav   a flag indicating whether only segments having
%               angular velocity are to be considered when determining
%               coverage.
%
%               [1,1] = size(needav); logical = class(needav)
%
%               When `needav' is true, segments without angular velocity
%               don't contribute to the coverage window; when `needav' is
%               false, all segments for `idcode' may contribute to the
%               coverage window.
%
%      level    the string defining the level (granularity) at which
%               the coverage is examined.
%
%               [1,c2] = size(level); char = class(level)
%
%                  or
%
%               [1,1] = size(level); cell = class(level)
%
%               Allowed values and corresponding meanings are:
%
%                  'SEGMENT'    The output coverage window contains
%                               intervals defined by the start and stop
%                               times of segments for the object
%                               designated by `idcode'.
%
%                  'INTERVAL'   The output coverage window contains
%                               interpolation intervals of segments for
%                               the object designated by `idcode'. For type
%                               1 segments, which don't have
%                               interpolation intervals, each epoch
%                               associated with a pointing instance is
%                               treated as a singleton interval; these
%                               intervals are added to the coverage
%                               window.
%
%                               All interpolation intervals are
%                               considered to lie within the segment
%                               bounds for the purpose of this summary:
%                               if an interpolation interval extends
%                               beyond the segment coverage interval,
%                               only its intersection with the segment
%                               coverage interval is considered to
%                               contribute to the total coverage.
%
%      tol      a tolerance value expressed in ticks of the spacecraft clock
%               associated with `idcode'.
%
%               [1,1] = size(tol); double = class(tol)
%
%               Before each interval is inserted into the coverage window,
%               the interval is intersected with the segment coverage
%               interval, then if the intersection is non-empty, it is
%               expanded by `tol': the left endpoint of the intersection
%               interval is reduced by `tol' and the right endpoint is
%               increased by `tol'. Adjusted interval endpoints, when
%               expressed as encoded SCLK, never are less than zero ticks.
%               Any intervals that overlap as a result of the expansion are
%               merged.
%
%               The coverage window returned when tol > 0 indicates the
%               coverage provided by the file to the CK readers cspice_ckgpav
%               and cspice_ckgp when that value of `tol' is passed to them as
%               an input.
%
%      timsys   the name of the time system to use in the output coverage
%               window.
%
%               [1,c3] = size(timsys); char = class(timsys)
%
%                  or
%
%               [1,1] = size(timsys); cell = class(timsys)
%
%               `timsys' may have the values:
%
%                  'SCLK'    Elements of `cover' are expressed in encoded
%                            SCLK ("ticks"), where the clock is associated
%                            with the object designated by `idcode'.
%
%                  'TDB'     Elements of `cover' are expressed as seconds
%                            past J2000 TDB.
%
%      room     a parameter specifying the maximum number of intervals that
%               can be accommodated by the dynamically allocated workspace
%               window used internally by this routine.
%
%               [1,1] = size(room); int32 = class(room)
%
%               It's not necessary to compute an accurate estimate of how
%               many intervals will be returned in `cover'; rather, the
%               user can pick a size considerably larger than what's
%               really required.
%
%      cover_i  an optional input describing a either an empty window or a
%               window array created from a previous cspice_ckcov call.
%
%               [2m,1] = size(cover_i), double = class(cover_i)
%
%                  or
%
%               [0,0] = size(cover_i), double = class(cover_i)
%
%               Inclusion of this window argument results in an output window
%               consisting of a union of the data retrieved from the `ck'
%               kernels and the data in `cover_i'.
%
%   the call:
%
%      [cover] = cspice_ckcov( ckfmn, idcode, needav, level,               ...
%                              tol,   timsys, room,   cover_i )
%
%         or
%
%      [cover] = cspice_ckcov( ckfnm, idcode, needav, level,               ...
%                              tol,   timsys, room           )
%
%   returns:
%
%      cover    the window containing the coverage for `idcode' available
%               from `ckfnm'.
%
%               [2p,1] = size(cover), double = class(cover)
%
%                  or
%
%               [0,1] = size(cover), double = class(cover)
%
%               When the coverage level is 'INTERVAL', this is the set
%               of time intervals for which data for `idcode' are present
%               in the file `ckfnm', merged with the set of time intervals
%               present in `cover_i' if provided on input. The merged
%               coverage is represented as the union of one or more
%               disjoint time intervals. The array `cover' contains the
%               pairs of endpoints of these intervals.
%
%               Each window defined as a pair of endpoints such that:
%
%                  window 1 = cover(1:2)
%                  window 2 = cover(3:4)
%                  window 3 = cover(5:6)
%                          ...
%                  window p = cover(2p-1,2p)
%
%               When the coverage level is 'SEGMENT', `cover' is computed
%               in a manner similar to that described above, but the coverage
%               intervals used in the computation are those of segments
%               rather than interpolation intervals within segments.
%
%               When `tol' is > 0, the intervals comprising the coverage
%               window for `idcode' are expanded by `tol' and any intervals
%               overlapping as a result are merged. The resulting window is
%               returned in `cover'. The expanded window in no case
%               extends beyond the segment bounds in either direction by
%               more than `tol'.
%
%               The interval endpoints contained in `cover' are encoded
%               spacecraft clock times if `timsys' is 'SCLK'; otherwise the
%               times are converted from encoded spacecraft clock to seconds
%               past J2000 TDB.
%
%               `cover' returns an empty set if `ckfnm' lacks coverage for
%               `idcode' and `cover_i' is not provided on input. `cover'
%               can overwrite `cover_i'.
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
%   1) Display the interval-level coverage for each object in a
%      specified CK file. Use tolerance of zero ticks. Do not
%      request angular velocity. Express the results in the TDB time
%      system.
%
%      Find the set of objects in the file. Loop over the contents
%      of the ID code set: find the coverage for each item in the
%      set and display the coverage.
%
%
%      Example code begins here.
%
%
%      function ckcov_ex1( )
%
%         %
%         % Local constants.
%         %
%         MAXIV  = 100000;
%         WINSIZ = 2 * MAXIV;
%         MAXOBJ = 1000;
%         LEVEL  = { 'SEGMENT', 'INTERVAL' };
%
%         %
%         % Load a leapseconds kernel and the SCLK corresponding to the
%         % input CK. Prompt the user for the file names.
%         %
%         % Note, neither cspice_ckcov or cspice_ckobj require these
%         % kernels to function. We need these data for output time
%         % conversion.
%         %
%         LSK  = input( 'Name of leapseconds kernel > ', 's' );
%         SCLK = input( 'Name of SCLK kernel        > ', 's' );
%         cspice_furnsh( LSK  )
%         cspice_furnsh( SCLK )
%
%         %
%         % Get the name of the CK file.
%         %
%         CK   = input( 'Name of CK file            > ', 's' );
%
%         %
%         % Find the set of objects in the CK file.
%         %
%         ids = cspice_ckobj( CK, MAXOBJ );
%
%         %
%         % Find the coverage for 'INTERVAL' and 'SEGMENT' levels.
%         %
%         for l=1:numel(LEVEL)
%
%            %
%            % We want to display the coverage for each object. Loop over
%            % the contents of the ID code set, find the coverage for
%            % each item in the set, and display the coverage.
%            %
%            for i=1:numel(ids)
%
%               %
%               % Extract the coverage data for object `ids(i)'.
%               %
%               cover    = cspice_ckcov( CK,  ids(i), 0,     LEVEL(l), ...
%                                        0.0, 'TDB',  WINSIZ         );
%               [row,col]= size(cover);
%
%               %
%               % Display a simple banner.
%               %
%               fprintf( '========================================\n')
%               fprintf( '%s level coverage for object %d\n', ...
%                                             char(LEVEL(l)), ids(i) )
%
%               %
%               %  `cover' has dimension 2Nx1, where 'row' has the value
%               %  2N with each window defined as a pair of endpoints
%               %  such that:
%               %
%               %  window 1 = cover(1:2)
%               %  window 2 = cover(3:4)
%               %  window 3 = cover(5:6)
%               %        ...
%               %  window N = cover(2N-1,2N)
%               %
%               % Loop from 1 to 1row' with stepsize 2.
%               %
%               for j=1:2:row
%
%                  %
%                  % Convert the endpoints to TDB calendar format time
%                  % strings and display them. Pass the endpoints in an
%                  % array, so cspice_timout returns an array of time
%                  % strings.
%                  %
%                  % Recall a vectorized input has dimension 1xM so
%                  % transpose the `cover' slice.
%                  %
%                  timstr = cspice_timout( cover(j:j+1)', ...
%                               'YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );
%                  fprintf('Interval: %d\n'  , (j+1)/2 )
%                  fprintf('   Start: %s\n'  , timstr(1,:) )
%                  fprintf('    Stop: %s\n\n', timstr(2,:) )
%
%               end
%
%            end
%
%         end
%
%         %
%         % Empty the kernel pool.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the LSK file named naif0010.tls, the SCLK file
%      named cas00145.tsc and the CK file named 08052_08057ra.bc, the
%      output was:
%
%
%      Name of leapseconds kernel > naif0010.tls
%      Name of SCLK kernel        > cas00145.tsc
%      Name of CK file            > 08052_08057ra.bc
%      ========================================
%      SEGMENT level coverage for object -82000
%      Interval: 1
%         Start: 2008 FEB 21 00:01:07.771 (TDB)
%          Stop: 2008 FEB 26 00:01:04.752 (TDB)
%
%      ========================================
%      INTERVAL level coverage for object -82000
%      Interval: 1
%         Start: 2008 FEB 21 00:01:07.771 (TDB)
%          Stop: 2008 FEB 23 22:53:30.001 (TDB)
%
%      Interval: 2
%         Start: 2008 FEB 23 22:58:13.999 (TDB)
%          Stop: 2008 FEB 24 02:22:25.913 (TDB)
%
%      Interval: 3
%         Start: 2008 FEB 24 02:27:49.910 (TDB)
%          Stop: 2008 FEB 24 19:46:33.470 (TDB)
%
%      Interval: 4
%         Start: 2008 FEB 24 19:49:33.469 (TDB)
%          Stop: 2008 FEB 25 04:25:21.250 (TDB)
%
%      Interval: 5
%         Start: 2008 FEB 25 04:29:33.248 (TDB)
%          Stop: 2008 FEB 25 15:23:44.971 (TDB)
%
%      Interval: 6
%         Start: 2008 FEB 25 15:24:12.971 (TDB)
%          Stop: 2008 FEB 25 20:25:04.843 (TDB)
%
%      Interval: 7
%         Start: 2008 FEB 25 20:25:48.843 (TDB)
%          Stop: 2008 FEB 26 00:01:04.752 (TDB)
%
%
%-Particulars
%
%   This routine provides an API via which applications can determine
%   the coverage a specified CK file provides for a specified
%   object.
%
%-Exceptions
%
%   1)  If the input file has transfer format, the error
%       SPICE(INVALIDFORMAT) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the input file is not a transfer file but has architecture
%       other than DAF, the error SPICE(INVALIDARCHTYPE) is signaled
%       by a routine in the call tree of this routine.
%
%   3)  If the input file is a binary DAF file of type other than CK,
%       the error SPICE(INVALIDFILETYPE) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If the CK file cannot be opened or read, an error is signaled
%       by a routine in the call tree of this routine. The output
%       window will not be modified.
%
%   5)  If the size of the output window argument `cover' is
%       insufficient to contain the actual number of intervals in the
%       coverage window for `idcode', an error is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If `tol' is negative, the error SPICE(VALUEOUTOFRANGE) is
%       signaled by a routine in the call tree of this routine.
%
%   7)  If `level' is not recognized, the error SPICE(INVALIDOPTION)
%       is signaled by a routine in the call tree of this routine.
%
%   8)  If `timsys' is not recognized, the error SPICE(NOTSUPPORTED)
%       is signaled by a routine in the call tree of this routine.
%
%   9)  If a time conversion error occurs, the error is signaled by a
%       routine in the call tree of this routine.
%
%   10) If the output time system is TDB, the CK subsystem must be
%       able to map `idcode' to the ID code of the associated spacecraft
%       clock. If this mapping cannot be performed, an error is
%       signaled by a routine in the call tree of this routine.
%
%   11) If the input CK type is not one of the supported CK types, the
%       error SPICE(NOTSUPPORTED) is signaled by a routine in the call
%       tree of this routine. This problem may indicate the version of
%       the SPICE Toolkit being used is outdated and a new version is
%       required.
%
%   12) If any of the input arguments, `ckfnm', `idcode', `needav',
%       `level', `tol', `timsys', `room' or `cover_i', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   13) If any of the input arguments, `ckfnm', `idcode', `needav',
%       `level', `tol', `timsys', `room' or `cover_i', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   This routine reads a C-kernel.
%
%   If the output time system is 'TDB', then a leapseconds kernel
%   and an SCLK kernel for the spacecraft clock associated with
%   `idcode' must be loaded before this routine is called.
%
%   If the ID code of the clock associated with `idcode' is not
%   equal to
%
%      idcode / 1000
%
%   then the kernel variable
%
%      CK_<idcode>_SCLK
%
%   must be present in the kernel pool to identify the clock
%   associated with `idcode'. This variable must contain the ID code
%   to be used for conversion between SCLK and TDB. Normally this
%   variable is provided in a text kernel loaded via cspice_furnsh.
%
%-Restrictions
%
%   1)  When this routine is used to accumulate coverage for `idcode'
%       provided by multiple CK files, the inputs `needav', `level', `tol',
%       and `timsys'  must have the same values for all files in order
%       for the result to be meaningful.
%
%-Required_Reading
%
%   CK.REQ
%   DAF.REQ
%   MICE.REQ
%   TIME.REQ
%   WINDOWS.REQ
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
%   -Mice Version 1.3.0, 26-NOV-2021 (EDW) (JDR)
%
%       Changed argument names "ck" and "cov" to "ckfnm" and "cover" for
%       consistency with other routines.
%
%       Edited the header to comply with NAIF standard. Modified existing
%       example to prompt the user for the required input data and hardcoded
%       the coverage levels within the code. Extended -Index_Entries.
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
%       Updated description of argument "room".
%
%   -Mice Version 1.2.1, 10-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.2.0, 03-APR-2012 (EDW)
%
%       Edits to Example code and comments. No change to Example code
%       functionality.
%
%       Renamed the argument "size" to "room". "size" is a Matlab function
%       name and it's seriously dumb to use a function name word as an
%       argument name.
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Explicitly described ID variables as "CK IDs."
%
%   -Mice Version 1.1.0, 29-DEC-2008 (EDW)
%
%       Edited description of "size"; "size" now defines the maximum
%       number of intervals for the internal workspace window.
%
%       The "cover_i" argument may now have the empty array value, [],
%       on input.
%
%       Added range restriction on size.
%
%       Corrected misspellings.
%
%   -Mice Version 1.0.0, 18-JUN-2007 (EDW)
%
%-Index_Entries
%
%   get coverage window for ck_object
%   get coverage start and stop time for ck_object
%   get coverage start and stop time for CK frame
%   get coverage start and stop time for CK instrument
%
%-&

function [cover] = cspice_ckcov( ckfnm, idcode, needav, level,             ...
                                 tol,   timsys, room,   cover_i )

   switch nargin
      case 7

         ckfnm  = zzmice_str(ckfnm);
         idcode = zzmice_int(idcode);
         needav = zzmice_int(needav);
         level  = zzmice_str(level);
         tol    = zzmice_dp(tol);
         timsys = zzmice_str(timsys);
         room   = zzmice_int(room, [1, int32(inf)/2] );

      case 8

         ckfnm  = zzmice_str(ckfnm);
         idcode = zzmice_int(idcode);
         needav = zzmice_int(needav);
         level  = zzmice_str(level);
         tol    = zzmice_dp(tol);
         timsys = zzmice_str(timsys);
         room   = zzmice_int(room, [1, int32(inf)/2] );
         cover_i= zzmice_win(cover_i);

      otherwise

         error ( ['Usage: [cover] = cspice_ckcov( _`ckfnm`_, idcode, '     ...
                          'needav, level, tol, timsys, room, [cover_i] )'] )

   end

%
% The call passed either seven or eight arguments. Branch accordingly.
%
if ( nargin == 7 )

   %
   % Call the MEX library.
   %
   try
      cover = mice('ckcov_c', ckfnm, idcode, needav,                       ...
                   level,     tol,   timsys, room    );
   catch spiceerr
      rethrow(spiceerr)
   end


else

   %
   % Call the MEX library.
   %
   try
      cover = mice('ckcov_c', ckfnm, idcode, needav,                       ...
                   level,     tol,   timsys, room );
   catch spiceerr
      rethrow(spiceerr)
   end

   %
   % Union the coverage window `cover' returned from the interface to the
   % input coverage window, `cover_i'.
   %
   cover = cspice_wnunid( cover, cover_i );

end
