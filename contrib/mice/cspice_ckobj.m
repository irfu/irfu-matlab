%-Abstract
%
%   CSPICE_CKOBJ finds the set of ID codes of all objects in a specified CK
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
%      ckfnm    the name(s) for SPICE CKs .
%
%               [n,c1] = size(ckfnm); char = class(ckfnm)
%
%                  or
%
%               [1,n] = size(ckfnm); cell = class(ckfnm)
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
%      ids_i    an optional input describing an array of CK ID
%               codes.
%
%               [r,1] = size(ids_i); int32 = class(ids_i)
%
%                  or
%
%               [0,0] = size(ids_i); int32 = class(ids_i)
%
%               Inclusion of this array results in an output array consisting
%               of a union of the data retrieved from the `ckfnm' kernels and
%               the data in `ids_i'.
%
%
%   the call:
%
%      [ids] = cspice_ckobj( ckfnm, room, ids_i )
%
%         or
%
%      [ids] = cspice_ckobj( ckfnm, room )
%
%   returns:
%
%      ids      the set of unique CK ID codes for which pointing data exists
%               in `ckfnm'.
%
%               [p,1] = size(ids), int32 = class(ids)
%
%               If `ids_i' exists in the argument list, `ids' returns as a
%               union of the data found in `ckfnm' and the data in
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
%   1) Use a simple function to display the CK IDs found in a CK, or set of
%      CKs, and the time coverage of the data corresponding to those IDs.
%
%      This example calls both cspice_ckobj and cspice_ckcov. In practice,
%      algorithms using cspice_ckobj will also use cspice_ckcov and
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
%      function ckobj_ex1( CK, SCLK, LEVEL )
%
%         MAXIV  = 100000;
%         WINSIZ = 2 * MAXIV;
%         MAXOBJ = 1000;
%         LSK    = 'naif0012.tls';
%
%         %
%         % Load a leapseconds kernel and the SCLK corresponding to the
%         % input CK.
%         %
%         % Note, neither cspice_ckcov or cspice_ckobj require these
%         % kernels to function. We need these data for output time
%         % conversion.
%         %
%         cspice_furnsh( LSK )
%         cspice_furnsh( SCLK)
%
%         %
%         % Find the set of objects in the CK file.
%         %
%         ids = cspice_ckobj( CK, MAXOBJ );
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
%            cover    = cspice_ckcov( CK,  ids(i), 0, LEVEL,               ...
%                                     0.0, 'TDB',     WINSIZ );
%            [row,col]= size(cover);
%
%            %
%            % Display a simple banner.
%            %
%            fprintf( '========================================\n')
%            fprintf( 'Coverage for object %d\n', ids(i) )
%
%            %
%            %  `cover' has dimension 2Nx1, where 'row' has the value 2N with
%            %  each window defined as a pair of endpoints such that:
%            %
%            %  window 1 = cover(1:2)
%            %  window 2 = cover(3:4)
%            %  window 3 = cover(5:6)
%            %        ...
%            %  window N = cover(2N-1,2N)
%            %
%            % Loop from 1 to 'row' with stepsize 2.
%            %
%            for j=1:2:row
%
%               %
%               % Convert the endpoints to TDB calendar format time strings
%               % and display them. Pass the endpoints in an array,
%               % so cspice_timout returns an array of time strings.
%               %
%               % Recall a vectorized input has dimension 1xM so transpose
%               % the `cover' slice.
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
%      platform, with the following variables as inputs
%
%         CK    = { '05357_05362ra.bc', '05362_06002ra.bc' };
%         SCLK  =   'cas00101.tsc';
%         LEVEL =   'INTERVAL';
%
%      the output was:
%
%
%      ========================================
%      Coverage for object -82000
%      Interval: 1
%         Start: 2005 DEC 23 00:01:07.900 (TDB)
%          Stop: 2005 DEC 23 15:36:55.540 (TDB)
%
%      Interval: 2
%         Start: 2005 DEC 23 15:37:39.539 (TDB)
%          Stop: 2005 DEC 23 16:59:35.508 (TDB)
%
%      Interval: 3
%         Start: 2005 DEC 23 17:00:43.507 (TDB)
%          Stop: 2005 DEC 24 13:55:59.025 (TDB)
%
%      Interval: 4
%         Start: 2005 DEC 24 13:56:19.024 (TDB)
%          Stop: 2005 DEC 24 17:25:42.944 (TDB)
%
%      Interval: 5
%         Start: 2005 DEC 24 17:29:22.942 (TDB)
%          Stop: 2005 DEC 25 00:47:58.774 (TDB)
%
%      Interval: 6
%         Start: 2005 DEC 25 00:56:26.770 (TDB)
%          Stop: 2005 DEC 26 06:59:58.077 (TDB)
%
%      Interval: 7
%         Start: 2005 DEC 26 07:00:17.826 (TDB)
%          Stop: 2005 DEC 26 13:52:05.918 (TDB)
%
%      Interval: 8
%         Start: 2005 DEC 26 13:54:17.917 (TDB)
%          Stop: 2005 DEC 26 14:10:53.911 (TDB)
%
%      Interval: 9
%         Start: 2005 DEC 26 14:15:13.909 (TDB)
%          Stop: 2005 DEC 27 13:34:53.121 (TDB)
%
%      Interval: 10
%         Start: 2005 DEC 27 13:36:01.370 (TDB)
%          Stop: 2005 DEC 27 18:57:09.247 (TDB)
%
%      Interval: 11
%         Start: 2005 DEC 27 19:01:05.245 (TDB)
%          Stop: 2005 DEC 27 19:11:01.241 (TDB)
%
%      Interval: 12
%         Start: 2005 DEC 27 19:14:01.240 (TDB)
%          Stop: 2005 DEC 28 00:01:01.130 (TDB)
%
%      Interval: 13
%         Start: 2005 DEC 28 00:01:05.130 (TDB)
%          Stop: 2005 DEC 28 18:05:00.713 (TDB)
%
%      Interval: 14
%         Start: 2005 DEC 28 18:07:04.712 (TDB)
%          Stop: 2005 DEC 28 18:23:00.706 (TDB)
%
%      Interval: 15
%         Start: 2005 DEC 28 18:31:04.703 (TDB)
%          Stop: 2005 DEC 28 18:31:16.703 (TDB)
%
%      Interval: 16
%         Start: 2005 DEC 28 18:31:44.703 (TDB)
%          Stop: 2005 DEC 29 13:19:00.269 (TDB)
%
%      Interval: 17
%         Start: 2005 DEC 29 13:21:12.268 (TDB)
%          Stop: 2005 DEC 29 13:30:28.264 (TDB)
%
%      Interval: 18
%         Start: 2005 DEC 29 13:34:48.263 (TDB)
%          Stop: 2005 DEC 29 19:33:00.125 (TDB)
%
%      Interval: 19
%         Start: 2005 DEC 29 19:42:32.121 (TDB)
%          Stop: 2005 DEC 30 10:47:23.773 (TDB)
%
%      Interval: 20
%         Start: 2005 DEC 30 10:47:59.773 (TDB)
%          Stop: 2005 DEC 30 11:02:11.768 (TDB)
%
%      Interval: 21
%         Start: 2005 DEC 30 11:02:55.767 (TDB)
%          Stop: 2005 DEC 30 21:10:11.534 (TDB)
%
%      Interval: 22
%         Start: 2005 DEC 30 21:13:43.532 (TDB)
%          Stop: 2005 DEC 31 14:58:31.123 (TDB)
%
%      Interval: 23
%         Start: 2005 DEC 31 15:06:47.119 (TDB)
%          Stop: 2005 DEC 31 15:48:11.104 (TDB)
%
%      Interval: 24
%         Start: 2005 DEC 31 15:49:11.103 (TDB)
%          Stop: 2006 JAN 01 15:18:34.561 (TDB)
%
%      Interval: 25
%         Start: 2006 JAN 01 15:20:30.560 (TDB)
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 110 lines have been
%      provided.
%
%
%   2) When Example #1 was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, with the following variables as inputs
%
%         CK    = { '05357_05362ra.bc', '05362_06002ra.bc' };
%         SCLK  =   'cas00101.tsc';
%         LEVEL =   'SEGMENT';
%
%      the output was:
%
%
%      ========================================
%      Coverage for object -82000
%      Interval: 1
%         Start: 2005 DEC 23 00:01:07.900 (TDB)
%          Stop: 2005 DEC 28 00:01:01.130 (TDB)
%
%      Interval: 2
%         Start: 2005 DEC 28 00:01:05.130 (TDB)
%          Stop: 2006 JAN 02 00:01:02.360 (TDB)
%
%
%-Particulars
%
%   This routine provides an API via which applications can determine
%   the set of objects for which there are pointing data in a
%   specified CK file.
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
%       by a routine in the call tree of this routine.
%
%   5)  If the size of the output set argument `ids' is insufficient
%       to contain the actual number of ID codes of objects covered by
%       the indicated CK file, an error is signaled by a routine in
%       the call tree of this routine.
%
%   6)  If any of the input arguments, `ckfnm', `room' or `ids_i', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `ckfnm', `room' or `ids_i', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   This routine reads a C-kernel.
%
%-Restrictions
%
%   1)  If an error occurs while this routine is updating the set
%       `ids', the set may be corrupted.
%
%-Required_Reading
%
%   CELLS.REQ
%   CK.REQ
%   DAF.REQ
%   FRAMES.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
%   SETS.REQ
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
%       Changed the input argument "ck" to "ckfnm" for consistency with other
%       routines.
%
%       Edited the header to comply with NAIF standard.
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
%   -Mice Version 1.3.0, 03-APR-2012 (EDW)
%
%       Edits to Example code and comments. No change to Example code
%       functionality.
%
%       Added error check on 'ids_i' to ensure the argument either has
%       shape [N,1] or is an empty array with shape [0,0].
%
%       Renamed the argument 'size' to 'room'. "size" is a Matlab function
%       name and it's seriously dumb to use a function name word as an
%       argument name.
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Explicitly described ID variables as "CK IDs."
%
%   -Mice Version 1.2.0, 13-AUG-2009 (EDW)
%
%       The union of 'ids_i'  with the interface return argument 'ids'
%       again calculated using the "unique" function, replacing "union."
%       This implementation results in the expected behavior of the
%       call in octave when 'ids_i' contains zero or one element.
%
%       Corrected typo in previous -Version entry.
%
%   -Mice Version 1.1.0, 29-DEC-2008 (EDW)
%
%       Corrected error in comment description for 'ids_i'.
%       Removed the line:
%
%          Note: 'ids_i' cannot be an empty array.
%
%       The argument can have the empty array value, [], on
%       input.
%
%       Corrected several typos.
%
%       'ids_i' union with interface return call now calculated
%       using the 'union' function instead of 'unique'.
%
%   -Mice Version 1.0.0, 19-JUN-2007 (EDW)
%
%-Index_Entries
%
%   find id codes of objects in CK file
%
%-&

function [ids] = cspice_ckobj( ckfnm, room, ids_i )

   switch nargin
      case 2

         ckfnm = zzmice_str(ckfnm);
         room  = zzmice_int(room);

      case 3

         ckfnm = zzmice_str(ckfnm);
         room  = zzmice_int(room);
         ids_i = zzmice_int(ids_i);

         %
         % Check 'ids_i' has dimension Nx1 or is an empty array.
         %
         is_valid =  (  (ndims(ids_i) == 2) && (size(ids_i, 2) == 1) )  ...
                        ||                                              ...
                        isequal( size(ids_i), [0,0] );

         if (~is_valid )
            error( 'MICE(BADARG): Argument ''ids_i'' must have size Nx1.' )
         end

      otherwise

         error ( 'Usage: [ids] = cspice_ckobj( _`ckfnm`_, room, [ids_i])' )

   end

   %
   % The call passed either two or three arguments. Branch accordingly.
   %
   if ( nargin == 2 )

      %
      % Call the MEX library.
      %
      try
         [ids] = mice('ckobj_c', ckfnm, room );
      catch spiceerr
         rethrow(spiceerr)
      end

   else

      %
      % Call the MEX library.
      %
      try
         ids = unique( [ [ids_i]; mice('ckobj_c', ckfnm, room ) ] );
      catch spiceerr
         rethrow(spiceerr)
      end

   end
