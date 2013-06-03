%-Abstract
%
%   CSPICE_CKOBJ returns the set of ID codes of all objects in a 
%   specified CK file.
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
%      ck     the string scalar or NxM character array of SPICE
%             CK file names 
%
%      size   a integer scalar defining the maximum number of IDs
%             to return from 'ck'
%
%      ids_i   an optional input describing an (Nx1) array of NAIF ID
%              codes.  Inclusion of this array results in an output
%              array consisting of a union of the data retrieved from
%              the 'ck' kernels and the data in 'ids_i'.
%
%   the call:
%
%      ids = cspice_ckobj( ck, size, ids_i)
%
%         or
%
%      ids = cspice_ckobj( ck, size)
%
%   returns:
%
%      ids   a Nx1 integer array containing the set of unique NAIF ID
%            codes for which pointing data exists in 'ck'. If 'ids_i' 
%            exists in the argument list, 'ids' returns as a union of
%            the coverage data found in 'ck' and the data in 'ids_i'.
%            'ids' can overwrite 'ids_i'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   function ckobj_t( CK, SCLK, LEVEL )
%
%      MAXIV  = 100000;
%      WINSIZ = 2 * MAXIV;
%      MAXOBJ = 1000;
%      LSK    = '/kernels/gen/lsk/naif0008.tls';
%
%      % 
%      % Load a leapseconds kernel for output time conversion.
%      % cspice_ckcov does not require a leapseconds kernel.
%      %
%      cspice_furnsh( LSK )
%      cspice_furnsh( SCLK)
%
%      %
%      % Find the set of objects in the CK file.
%      %
%      ids = cspice_ckobj( CK, MAXOBJ );
%
%      %
%      % We want to display the coverage for each object. Loop over
%      % the contents of the ID code set, find the coverage for
%      % each item in the set, and display the coverage.
%      %
%      for i=1:numel(ids)
%   
%         %
%         % Extract the coverage data for object 'ids(i)'.
%         %
%         cover     = cspice_ckcov( CK, ids(i), 0, LEVEL,  0.0, 'TDB', WINSIZ );
%         [row,col] = size(cover);
%
%         %
%         % Display a simple banner.
%         %
%         fprintf( '========================================\n')
%         fprintf( 'Coverage for object %d\n', ids(i) )
%
%         %
%         %  'cover' has dimension 2Nx1, where 'row' has the value 2N with
%         %  each window defined as a pair of endpoints such that:
%         %
%         %  window 1 = cover(1:2)
%         %  window 2 = cover(3:4)
%         %  window 3 = cover(5:6)
%         %        ...
%         %  window N = cover(2N-1,2N)
%         %
%         % Loop from 1 to 'row' with step size 2.
%         %
%         for j=1:2:row
%      
%            %
%            % Convert the endpoints to TDB calendar format time strings
%            % and display them. Pass the endpoints in an array,
%            % so cspice_timout returns an array of time strings.
%            %
%            % Recall a vectorized input has dimension 1xM so transpose
%            % the 'cover' slice.
%            %
%            timstr = cspice_timout( cover(j:j+1)', ...
%                                'YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );
%            fprintf('Interval: %d\n'  , (j+1)/2 )
%            fprintf('   Start: %s\n'  , timstr(1,:) )
%            fprintf('    Stop: %s\n\n', timstr(2,:) )
%
%         end
%
%      end
%        
%      %
%      % Empty the kernel pool.
%      %
%      cspice_kclear
%
%   Example (1):
%
%      Assign a CK kernel list as and SCLK:
%
%      >> SCLK = '/kernels/cassini/sclk/cas00101.tsc';
%      >> CK   = { '/kernels/cassini/ck/05357_05362ra.bc', ...
%                  '/kernels/cassini/ck/05362_06002ra.bc'   };
%
%      Output data using the "INTERVAL" level.
%
%      >> ckcov_t( CK, SCLK, 'INTERVAL' )
%
%   MATLAB outputs:
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
%                ... continued ...
%      
%      Interval: 24
%         Start: 2005 DEC 31 15:49:11.103 (TDB)
%          Stop: 2006 JAN 01 15:18:34.561 (TDB)
% 
%      Interval: 25
%         Start: 2006 JAN 01 15:20:30.560 (TDB)
%          Stop: 2006 JAN 01 16:43:38.528 (TDB)
%
%      Interval: 26
%         Start: 2006 JAN 01 16:45:02.528 (TDB)
%          Stop: 2006 JAN 01 22:52:10.386 (TDB)
%
%      Interval: 27
%         Start: 2006 JAN 01 22:52:38.386 (TDB)
%          Stop: 2006 JAN 02 00:01:02.360 (TDB)
%
%   Example (2):
%
%      Output data using the "SEGMENT" level.
%
%      >> ckcov_t( CK, SCLK, 'SEGMENT' ) 
%
%   MATLAB outputs:
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
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckobj_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 13-AUG-2009, EDW (JPL)
%
%      The union of 'ids_i'  with the interface return argument 'ids' 
%      again calculated using the "unique" function, replacing "union." 
%      This implementation results in the expected behavior of the 
%      call in octave when 'ids_i' contains zero or one element.
%
%      Corrected typo in previous Version entry.
%
%   -Mice Version 1.1.0, 29-DEC-2008, EDW (JPL)
%
%      Corrected error in comment description for 'ids_i'.
%      Removed the line:
%
%         Note: 'ids_i' cannot be an empty array.
%
%      The argument can have the empty array value, [], on
%      input.
%
%      Corrected several typos.
%
%      'ids_i' union with interface return call now calculated
%      using the 'union' function instead of 'unique'.
%
%   -Mice Version 1.0.0, 19-JUN-2007, EDW (JPL)
%
%-Index_Entries
% 
%   find ID codes in ck file 
% 
%-&

function [ids] = cspice_ckobj( ck, size, ids_i )

   switch nargin
      case 2

         ck   = zzmice_str(ck);
         size = zzmice_int(size);

      case 3

         ck   = zzmice_str(ck);
         size = zzmice_int(size);
         ids_i= zzmice_int(ids_i);

      otherwise

         error ( 'Usage: [ids] = cspice_ckobj( _`ck`_, size, [ids_i])' )

   end

   %
   % The call passed either two or three arguments. Branch accordingly.
   %
   if ( nargin == 2 )

      %
      % Call the MEX library.
      %
      try
         [ids] = mice('ckobj_c', ck, size );
      catch
         rethrow(lasterror)
      end

   else

      %
      % Call the MEX library.
      %
      try
         ids = unique( [ [ids_i]; mice('ckobj_c', ck, size ) ] );
      catch
         rethrow(lasterror)
      end

   end





