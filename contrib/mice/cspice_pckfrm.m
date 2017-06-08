%-Abstract
%
%   CSPICE_PCKFRM returns the set of reference frame class ID codes of 
%   all frames in a specified binary PCK file. 
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
%      pck     is the name, or cell of names, of a binary PCK file(s).
%
%              [n,c1] = size(pck); char = class(pck)
%
%                 or
%
%              [1,n] = size(pck); cell = class(pck)
%
%      room    an integer scalar defining the maximum number of PCK IDs to
%              return from 'pck'.
%
%              [1,1] = size(room); int32 = class(room)
%
%      ids_i   an optional input describing an (mx1) array of PCK
%              ID codes. Inclusion of this array results in an output
%              array consisting of a union of the data retrieved from
%              the 'pck' kernels and the data in 'ids_i'.
%
%              [m,1] = size(ids_i); int32 = class(ids_i)
%
%                 or
%
%              [0,0] = size(ids_i); int32 = class(ids_i)
%
%   the call:
%
%      ids = cspice_pckfrm( pck, room, ids_i)
%
%         or
%
%      ids = cspice_pckfrm( pck, room)
%
%   returns:
%
%      ids   set of unique reference frame class ID codes of each frame for
%            which data are present in the indicated PCK file. If 'ids_i'
%            exists in the argument list, 'ids' returns as a union of
%            the IDs found in 'pck' and the IDs in 'ids_i'.  'ids' 
%            can overwrite 'ids_i'.
%
%            [p,1] = size(ids); int32 = class(ids)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Use a simple routine to display the coverage for each object in a 
%   specified PCK file(s). Find the set of objects in the file(s); for 
%   each object, find and display the coverage. 
%
%      function pckfrm_t(pck)
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
%        cspice_furnsh( LSK )
%   
%        %
%        % Find the set of frames in the pck file.
%        %
%        ids = cspice_pckfrm( pck, MAXIV );
%   
%        %
%        % We want to display the coverage for each frame. Loop over
%        % the contents of the ID code set, find the coverage for
%        % each item in the set, and display the coverage.
%        %
%        for i=1:numel(ids)
%   
%           %
%           % Find the coverage window for the current frame. 
%           %
%           cover     = cspice_pckcov( pck, ids(i), WINSIZ );
%           [row,col] = size(cover);
%   
%           %
%           % Display a simple banner.
%           %
%           fprintf( '========================================\n')
%           fprintf( 'Coverage for frame %d\n', ids(i) )
%   
%           %
%           %  'cover' has dimension 2Nx1, where 'row' has the value 2N with
%           %  each window defined as a pair of endpoints such that:
%           %
%           %  window 1 = cover(1:2)
%           %  window 2 = cover(3:4)
%           %  window 3 = cover(5:6)
%           %        ...
%           %  window N = cover(2N-1,2N)
%           %
%           % Loop from 1 to 'row' with step size 2.
%           %
%           for j=1:2:row
%   
%              %
%              % Convert the endpoints to TDB calendar format time strings
%              % and display them. Pass the endpoints in an array,
%              % so cspice_timout returns an array of time strings.
%              %
%              % Recall a vectorized input has dimension 1xM so transpose
%              % the 'cover' slice.
%              %
%              timstr = cspice_timout( cover(j:j+1)', ...
%                                  'YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );
%              fprintf('Interval: %d\n'  , (j+1)/2 )
%              fprintf('   Start: %s\n'  , timstr(1,:) )
%              fprintf('    Stop: %s\n\n', timstr(2,:) )
%   
%           end
%   
%        end
%   
%        %
%        % Empty the kernel pool.
%        %
%        cspice_kclear
%
%   MATLAB outputs:
%
%      >> pckfrm_t( {'earth_latest_high_prec.bpc', ...
%                    'moon_pa_de421_1900-2050.bpc' })
%      ========================================
%      Coverage for frmect 3000
%      Interval: 1
%         Start: 2000 JAN 01 00:01:04.183 (TDB)
%          Stop: 2017 APR 10 00:01:09.185 (TDB)
%      
%      ========================================
%      Coverage for frmect 31006
%      Interval: 1
%         Start: 1900 JAN 01 00:00:00.000 (TDB)
%          Stop: 2051 JAN 01 00:00:00.000 (TDB)
%
%-Particulars
%
%   This routine provides an API via which applications can determine 
%   the set of reference frames for which there are data in a 
%   specified PCK file. 
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pckfrm_c.
%
%   MICE.REQ
%   CELLS.REQ 
%   DAF.REQ 
%   SETS.REQ 
%   PCK.REQ 
%
%-Version
%
%   -Mice Version 1.0.0, 03-JAN-2017, EDW (JPL)
%
%-Index_Entries
%
%   find frame class id codes of frames in binary pck file 
%
%-&

function [ids] = cspice_pckfrm( pck, room, ids_i )

   switch nargin
      case 2

         pck  = zzmice_str(pck);
         room = zzmice_int(room);

      case 3

         pck  = zzmice_str(pck);
         room = zzmice_int(room);
         ids_i= zzmice_int(ids_i);

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

         error( 'Usage: [ids] = cspice_pckfrm( _`pck`_, room, [ids_i])' )

   end

%
% The call passed either two or three arguments. Branch accordingly.
%
if ( nargin == 2 )

   %
   % Call the MEX library.
   %
   try
      ids = mice('pckfrm_c', pck, room );
   catch
      rethrow(lasterror)
   end

else

   %
   % Call the MEX library.
   %
   try
      ids = unique( [ [ids_i]; mice('pckfrm_c', pck, room ) ]  );
   catch
      rethrow(lasterror)
   end


end

