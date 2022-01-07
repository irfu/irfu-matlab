%-Abstract
%
%   CSPICE_SPKCOV finds the coverage window for a specified ephemeris object
%   in a specified SPK file.
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
%      spkfnm   the name(s) of SPICE SPK file(s).
%
%               [n,c1] = size(spkfnm); char = class(spkfnm)
%
%                 or
%
%               [1,n] = size(spkfnm); cell = class(spkfnm)
%
%      idcode   the integer ID code of an object for which ephemeris data are
%               expected to exist in the specified SPK file.
%
%               [1,1] = size(idcode); int32 = class(idcode)
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
%               window array created from a previous cspice_spkcov call.
%
%               [2m,1] = size(cover_i), double = class(cover_i)
%
%                  or
%
%               [0,0] = size(cover_i), double = class(cover_i)
%
%               Inclusion of this window argument results in an output window
%               consisting of a union of the data retrieved from the `spkfnm'
%               kernels and the data in `cover_i'.
%
%   the call:
%
%      [cover] = cspice_spkcov( spkfnm, idcode, room, cover_i )
%
%         or
%
%      [cover] = cspice_spkcov( spkfnm, idcode, room )
%
%   returns:
%
%      cover    the window containing the coverage for `idcode', i.e.
%               the set of time intervals for which `idcode'
%               data exist in the file `spkfnm'.
%
%               [2p,1] = size(cover), double = class(cover)
%
%                  or
%
%               [0,1] = size(cover), double = class(cover)
%
%               The array `cover' contains the pairs of endpoints of these
%               intervals.
%
%               Each window defined as a pair of endpoints such that:
%
%                  window 1 = cover(1:2)
%                  window 2 = cover(3:4)
%                  window 3 = cover(5:6)
%                           ...
%                  window p = cover(2p-1,2p)
%
%               The interval endpoints contained in `cover' are ephemeris
%               times, expressed as seconds past J2000 TDB.
%
%               `cover' returns an empty set if `spkfnm' lacks coverage for
%               `idcode'. If `cover_i' exists in the argument list, `cover'
%               returns as a union of the coverage data found in 'spkfnm' and
%               the data in `cover_i'. `cover' can overwrite `cover_i'.
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
%   1) Use a simple function to display the SPK IDs found in an SPK or set of
%      SPKs, and the time coverage of the data corresponding to those IDs.
%
%      This example calls both cspice_spkobj and cspice_spkcov. In practice,
%      algorithms using cspice_spkobj will also use cspice_spkcov and
%      vice-versa.
%
%      Use the LSK kernel below to load the leap seconds and time
%      constants required for the time conversions.
%
%         naif0012.tls
%
%
%      Example code begins here.
%
%
%      function spkcov_ex1( SPK )
%
%         MAXIV  = 1000;
%         WINSIZ = 2 * MAXIV;
%         LSK    = 'naif0012.tls';
%
%         %
%         % Note, neither cspice_spkcov or cspice_spkobj requires this
%         % kernel to function. We need the data for output time
%         % conversion.
%         %
%         cspice_furnsh( LSK )
%
%         %
%         % Find the set of objects in the SPK file.
%         %
%         ids = cspice_spkobj( SPK, MAXIV );
%
%         %
%         % We want to display the coverage for each object. Loop over
%         % the contents of the ID code set, find the coverage for
%         % each item in the set, and display the coverage.
%         %
%         for i=1:numel(ids)
%
%            %
%            % Extract the coverage data for object 'ids(i)'.
%            %
%            cover     = cspice_spkcov( SPK, ids(i), WINSIZ );
%            [row,col] = size(cover);
%
%            %
%            % Display a simple banner.
%            %
%            fprintf( '========================================\n')
%            fprintf( 'Coverage for object %d\n', ids(i) )
%
%            %
%            %  'cover' has dimension 2Nx1, where 'row' has the value 2N with
%            %  each window defined as a pair of endpoints such that:
%            %
%            %  window 1 = cover(1:2)
%            %  window 2 = cover(3:4)
%            %  window 3 = cover(5:6)
%            %        ...
%            %  window N = cover(2N-1,2N)
%            %
%            % Loop from 1 to 'row' with step size 2.
%            %
%            for j=1:2:row
%
%               %
%               % Convert the endpoints to TDB calendar format time strings
%               % and display them. Pass the endpoints in an array,
%               % so cspice_timout returns an array of time strings.
%               %
%               % Recall a vectorized input has dimension 1xM so transpose
%               % the 'cover' slice.
%               %
%               timstr = cspice_timout( cover(j:j+1)', ...
%                                   'YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );
%               fprintf('Interval: %d\n'  , (j+1)/2 )
%               fprintf('   Start: %s\n'  , timstr(1,:) )
%               fprintf('    Stop: %s\n\n', timstr(2,:) )
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
%      platform, with the following variable as input
%
%         SPK = { 'sat393.bsp', 'ura112.bsp' };
%
%      the output was:
%
%
%      ========================================
%      Coverage for object 3
%      Interval: 1
%         Start: 1900 JAN 01 00:00:41.183 (TDB)
%          Stop: 2099 DEC 24 00:01:07.183 (TDB)
%
%      ========================================
%      Coverage for object 6
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 7
%      Interval: 1
%         Start: 1900 JAN 01 00:00:41.183 (TDB)
%          Stop: 2099 DEC 24 00:01:07.183 (TDB)
%
%      ========================================
%      Coverage for object 10
%      Interval: 1
%         Start: 1900 JAN 01 00:00:41.183 (TDB)
%          Stop: 2099 DEC 24 00:01:07.183 (TDB)
%
%      ========================================
%      Coverage for object 399
%      Interval: 1
%         Start: 1900 JAN 01 00:00:41.183 (TDB)
%          Stop: 2099 DEC 24 00:01:07.183 (TDB)
%
%      ========================================
%      Coverage for object 610
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 611
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 612
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 613
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 614
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 615
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 616
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 617
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 632
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 633
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 634
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%          Stop: 2050 JAN 01 00:01:08.183 (TDB)
%
%      ========================================
%      Coverage for object 649
%      Interval: 1
%         Start: 1950 JAN 01 00:00:41.183 (TDB)
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 174 lines have been
%      provided.
%
%
%   2) When Example #1 was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, with the following variable as input
%
%         SPK = { 'mgs_ext12_ipng_mgs95j.bsp', 'mgs_ext26_ipng_mgs95j.bsp' };
%
%      the output was:
%
%
%      ========================================
%      Coverage for object -94
%      Interval: 1
%         Start: 2003 JUL 23 00:00:00.000 (TDB)
%          Stop: 2003 OCT 15 01:00:00.000 (TDB)
%
%      Interval: 2
%         Start: 2006 OCT 11 00:00:00.000 (TDB)
%          Stop: 2006 NOV 08 01:00:00.000 (TDB)
%
%
%-Particulars
%
%   This routine provides an API via which applications can determine
%   the coverage a specified SPK file provides for a specified
%   ephemeris object.
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
%   3)  If the input file is a binary DAF file of type other than SPK,
%       the error SPICE(INVALIDFILETYPE) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If the SPK file cannot be opened or read, an error is signaled
%       by a routine in the call tree of this routine. The output
%       window will not be modified.
%
%   5)  If the size of the output window argument `cover' is
%       insufficient to contain the actual number of intervals in the
%       coverage window for `idcode', an error is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If any of the input arguments, `spkfnm', `idcode', `room' or
%       `cover_i', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   7)  If any of the input arguments, `spkfnm', `idcode', `room' or
%       `cover_i', is not of the expected type, or it does not have
%       the expected dimensions and size, an error is signaled by the
%       Mice interface.
%
%-Files
%
%   This routine reads an SPK file.
%
%-Restrictions
%
%   1)  If an error occurs while this routine is updating the window
%       `cover_f', the window may be corrupted.
%
%-Required_Reading
%
%   DAF.REQ
%   MICE.REQ
%   SPK.REQ
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
%       Changed the argument names "spk" and "size" to "spkfnm" and "room",
%       for consistency with other routines.
%
%       Edited the header to comply with NAIF standard. Extended
%       -Index_Entries.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
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
%      Edited -I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.2.0, 03-APR-2012 (EDW)
%
%      Edits to Example code and comments. No change to Example code
%      functionality.
%
%      Renamed the argument 'size' to 'room'. "size" is a Matlab function name
%      and it's seriously dumb to use a function name word as an argument
%      name.
%
%      Edited -I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.1.0, 29-DEC-2008 (EDW)
%
%      Edited description of 'size'; 'size' now defines the maximum
%      number of intervals for the internal workspace window.
%
%      The 'cover_i' argument may now have the empty array value, [],
%      on input.
%
%      Added range restriction on size.
%
%      Corrected several typos.
%
%   -Mice Version 1.0.0, 18-JUN-2007 (EDW)
%
%-Index_Entries
%
%   get coverage window for spk_object
%   get coverage start and stop time for spk_object
%   get coverage start and stop time for ephemeris_object
%   get coverage start and stop time for body
%
%-&

function [cover] = cspice_spkcov( spkfnm, idcode, room, cover_i )

   switch nargin
      case 3

         spkfnm = zzmice_str(spkfnm);
         idcode = zzmice_int(idcode);
         room   = zzmice_int(room, [1, int32(inf)/2] );

      case 4

         spkfnm  = zzmice_str(spkfnm);
         idcode  = zzmice_int(idcode);
         room    = zzmice_int(room, [1, int32(inf)/2] );
         cover_i = zzmice_win(cover_i);

      otherwise

         error ( [ 'Usage: [cover] = cspice_spkcov( _`spkfnm`_, '          ...
                                     'idcode, room, [cover_i] )' ] )

   end

%
% The call passed either three or four arguments. Branch accordingly.
%
if ( nargin == 3 )

   %
   % Call the MEX library.
   %
   try
      [cover] = mice('spkcov_c', spkfnm, idcode, room );
   catch spiceerr
      rethrow(spiceerr)
   end

else

   %
   % Call the MEX library.
   %
   try
      cover = mice('spkcov_c', spkfnm, idcode, room );
   catch spiceerr
      rethrow(spiceerr)
   end

   cover = cspice_wnunid( cover, cover_i );

end
