%-Abstract
%
%   CSPICE_SPKOBJ finds the set of ID codes of all objects in a specified SPK
%   file.
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
%      spkfnm   the name, or cell of names, of SPICE SPK file(s).
%
%               [n,c1] = size(spkfnm); char = class(spkfnm)
%
%                  or
%
%               [1,n] = size(spkfnm); cell = class(spkfnm)
%
%      room     a parameter specifying the maximum number of elements that
%               can be accommodated by the dynamically allocated workspace
%               cell used internally by this routine.
%
%               [1,1] = size(room); int32 = class(room)
%
%               It's not necessary to compute an accurate estimate of how
%               many elements will be returned in `ids'; rather, the
%               user can pick a size considerably larger than what's
%               really required.
%
%      ids_i    an optional input describing an (Nx1) array of SPK
%               ID codes.
%
%               [m,1] = size(ids_i); int32 = class(ids_i)
%
%                  or
%
%               [0,0] = size(ids_i); int32 = class(ids_i)
%
%               Inclusion of this array results in an output array consisting
%               of a union of the data retrieved from the `spkfnm' kernels and
%               the data in `ids_i'.
%
%   the call:
%
%      [ids] = cspice_spkobj( spkfnm, room, ids_i )
%
%         or
%
%      [ids] = cspice_spkobj( spkfnm, room )
%
%   returns:
%
%      ids      an array containing the set of unique NAIF ID codes for which
%               ephemeris data exists in `spkfnm'.
%
%               [p,1] = size(ids), int32 = class(ids)
%
%               If `ids_i' exists in the argument list, `ids' returns as a
%               union of the data found in `spkfnm' and the data in
%               `ids_i'. `ids' can overwrite `ids_i'.
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
%      function spkobj_ex1( spknam )
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
%         ids = cspice_spkobj( spknam, MAXIV );
%
%         %
%         % We want to display the coverage for each object. Loop over
%         % the contents of the ID code set, find the coverage for
%         % each item in the set, and display the coverage.
%         %
%         for i=1:numel(ids)
%
%            %
%            % Extract the coverage data for object ids(i).
%            %
%            cover     = cspice_spkcov( spknam, ids(i), WINSIZ );
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
%         spknam = { 'sat393.bsp', 'ura112.bsp' };
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
%         spknam = { 'mgs_ext12_ipng_mgs95j.bsp',                          ...
%                    'mgs_ext26_ipng_mgs95j.bsp' };
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
%   the set of objects for which there are ephemeris data in a
%   specified SPK file.
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
%       by a routine in the call tree of this routine.
%
%   5)  If the size of the output set argument `ids' is insufficient to
%       contain the actual number of ID codes of objects covered by
%       the indicated SPK file, an error is signaled by a routine in
%       the call tree of this routine.
%
%   6)  If any of the input arguments, `spkfnm', `room' or `ids_i', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `spkfnm', `room' or `ids_i', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   This routine reads an SPK file.
%
%-Restrictions
%
%   1)  If an error occurs while this routine is updating the set
%       `ids', the set may be corrupted.
%
%-Required_Reading
%
%   CELLS.REQ
%   DAF.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
%   SETS.REQ
%   SPK.REQ
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
%   -Mice Version 1.4.0, 26-NOV-2021 (EDW) (JDR)
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
%   -Mice Version 1.3.0, 03-APR-2012 (EDW)
%
%      Edits to Example code and comments. No change to Example code
%      functionality.
%
%      Added error check on 'ids_i' to ensure the argument either has
%      shape [N,1] or is an empty array with shape [0,0].
%
%      Renamed the argument 'size' to 'room'. "size" is a Matlab function name
%      and it's seriously dumb to use a function name word as an argument
%      name.
%
%      Edited -I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.2.0, 13-AUG-2009 (EDW)
%
%      The union of 'ids_i'  with the interface return argument 'ids'
%      again calculated using the "unique" function, replacing "union."
%      This implementation results in the expected behavior of the
%      call in octave when 'ids_i' contains zero or one element.
%
%   -Mice Version 1.1.0, 29-DEC-2008 (EDW)
%
%      Corrected error in comment description for 'ids_i'.
%      Removed the line:
%
%         Note: 'ids_i' cannot be an empty array.
%
%      The argument can have the empty array value, [], on
%      input.
%
%      'ids_i' union with interface return call now calculated
%      using the "union" function instead of "unique."
%
%   -Mice Version 1.0.0, 18-JUN-2007 (EDW)
%
%-Index_Entries
%
%   find id codes of ephemeris objects in SPK file
%   find id codes of bodies in SPK file
%
%-&

function [ids] = cspice_spkobj( spkfnm, room, ids_i )

   switch nargin
      case 2

         spkfnm = zzmice_str(spkfnm);
         room   = zzmice_int(room);

      case 3

         spkfnm = zzmice_str(spkfnm);
         room   = zzmice_int(room);
         ids_i  = zzmice_int(ids_i);

         %
         % Check 'ids_i' has dimension Nx1 or is an empty array.
         %
         is_valid =  (  (ndims(ids_i) == 2) && (size(ids_i, 2) == 1) )     ...
                        ||                                                 ...
                        isequal( size(ids_i), [0,0] );

         if (~is_valid )
            error( 'MICE(BADARG): Argument ''ids_i'' must have size Nx1.' )
         end

      otherwise

         error( 'Usage: [ids] = cspice_spkobj( _`spkfnm`_, room, [ids_i])' )

   end

%
% The call passed either two or three arguments. Branch accordingly.
%
if ( nargin == 2 )

   %
   % Call the MEX library.
   %
   try
      ids = mice('spkobj_c', spkfnm, room );
   catch spiceerr
      rethrow(spiceerr)
   end

else

   %
   % Call the MEX library.
   %
   try
      ids = unique( [ [ids_i]; mice('spkobj_c', spkfnm, room ) ]  );
   catch spiceerr
      rethrow(spiceerr)
   end


end
