%-Abstract
%
%   CSPICE_PCKFRM finds the set of reference frame class ID codes of all
%   frames in a specified binary PCK file.
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
%      pckfnm   the name, or cell of names, of a binary PCK file(s).
%
%               [n,c1] = size(pckfnm); char = class(pckfnm)
%
%                  or
%
%               [1,n] = size(pckfnm); cell = class(pckfnm)
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
%      ids_i    an optional input describing an (mx1) array of PCK
%               ID codes.
%
%               [m,1] = size(ids_i); int32 = class(ids_i)
%
%                  or
%
%               [0,0] = size(ids_i); int32 = class(ids_i)
%
%               Inclusion of this array results in an output array consisting
%               of a union of the data retrieved from the `pckfnm' kernels and
%               the data in `ids_i'.
%
%   the call:
%
%      [ids] = cspice_pckfrm( pckfnm, room, ids_i )
%
%         or
%
%      [ids] = cspice_pckfrm( pckfnm, room )
%
%   returns:
%
%      ids      set of unique reference frame class ID codes of each frame for
%               which data are present in the indicated PCK file.
%
%               [p,1] = size(ids); int32 = class(ids)
%
%               If `ids_i' exists in the argument list, `ids' returns as a
%               union of the IDs found in `pckfnm' and the IDs in `ids_i'.
%               `ids' can overwrite `ids_i'.
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
%   1) Use a simple routine to display the coverage for each object in a
%      specified PCK file(s). Find the set of objects in the file(s); for
%      each object, find and display the coverage.
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
%      function pckfrm_ex1( pcknam )
%
%         MAXIV  = 1000;
%         WINSIZ = 2 * MAXIV;
%         LSK    = 'naif0012.tls';
%
%         %
%         % Note, neither cspice_pckcov or cspice_pckfrm requires this
%         % kernel to function. We need the data for output time
%         % conversion.
%         %
%         cspice_furnsh( LSK )
%
%         %
%         % Find the set of frames in the PCK file.
%         %
%         ids = cspice_pckfrm( pcknam, MAXIV );
%
%         %
%         % We want to display the coverage for each frame. Loop over
%         % the contents of the ID code set, find the coverage for
%         % each item in the set, and display the coverage.
%         %
%         for i=1:numel(ids)
%
%            %
%            % Find the coverage window for the current frame.
%            %
%            cover     = cspice_pckcov( pcknam, ids(i), WINSIZ );
%            [row,col] = size(cover);
%
%            %
%            % Display a simple banner.
%            %
%            fprintf( '========================================\n')
%            fprintf( 'Coverage for frame %d\n', ids(i) )
%
%            %
%            %  `cover' has dimension 2Nx1, where `row' has the value 2N with
%            %  each window defined as a pair of endpoints such that:
%            %
%            %  window 1 = cover(1:2)
%            %  window 2 = cover(3:4)
%            %  window 3 = cover(5:6)
%            %        ...
%            %  window N = cover(2N-1,2N)
%            %
%            % Loop from 1 to `row' with step size 2.
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
%      platform, with the following variable as input
%
%         pcknam = { 'earth_720101_070426.bpc',                            ...
%                    'moon_pa_de421_1900-2050.bpc' };
%
%      the output was:
%
%
%      ========================================
%      Coverage for frame 3000
%      Interval: 1
%         Start: 1962 JAN 20 00:00:41.184 (TDB)
%          Stop: 2007 APR 26 00:01:05.185 (TDB)
%
%      ========================================
%      Coverage for frame 31006
%      Interval: 1
%         Start: 1900 JAN 01 00:00:00.000 (TDB)
%          Stop: 2051 JAN 01 00:00:00.000 (TDB)
%
%
%-Particulars
%
%   This routine provides an API via which applications can determine
%   the set of reference frames for which there are data in a
%   specified PCK file.
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
%   3)  If the input file is a binary DAF file of type other than PCK,
%       the error SPICE(INVALIDFILETYPE) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If the PCK file cannot be opened or read, an error is signaled
%       by a routine in the call tree of this routine.
%
%   5)  If the size of the output set argument `ids' is insufficient to
%       contain the actual number of ID codes of frames covered by the
%       indicated PCK file, an error is signaled by a routine in the
%       call tree of this routine.
%
%   6)  If any of the input arguments, `pckfnm', `room' or `ids_i', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `pckfnm', `room' or `ids_i', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
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
%   FRAMES.REQ
%   MICE.REQ
%   PCK.REQ
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
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Changed the argument names "pck" and "size" to "pckfnm" and "room",
%       for consistency with other routines.
%
%       Edited the header to comply with NAIF standard.
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
%   -Mice Version 1.0.0, 03-JAN-2017 (EDW)
%
%-Index_Entries
%
%   find frame class id codes of frames in binary PCK file
%
%-&

function [ids] = cspice_pckfrm( pckfnm, room, ids_i )

   switch nargin
      case 2

         pckfnm = zzmice_str(pckfnm);
         room   = zzmice_int(room);

      case 3

         pckfnm = zzmice_str(pckfnm);
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

         error( 'Usage: [ids] = cspice_pckfrm( _`pckfnm`_, room, [ids_i] )' )

   end

%
% The call passed either two or three arguments. Branch accordingly.
%
if ( nargin == 2 )

   %
   % Call the MEX library.
   %
   try
      ids = mice('pckfrm_c', pckfnm, room );
   catch spiceerr
      rethrow(spiceerr)
   end

else

   %
   % Call the MEX library.
   %
   try
      ids = unique( [ [ids_i]; mice('pckfrm_c', pckfnm, room ) ]  );
   catch spiceerr
      rethrow(spiceerr)
   end


end
