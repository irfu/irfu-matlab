%-Abstract
%
%   CSPICE_PCKCOV returns the coverage window for a specified reference 
%   frame in a specified binary PCK file. 
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
%      pck       is the name, or cell of names, of a binary PCK file(s).
%
%                [n,c1] = size(pck); char = class(pck)
%
%                  or
%
%                [1,n] = size(pck); cell = class(pck)
%
%      idcode    is the integer frame class ID code of a PCK reference
%                 frame for which data are expected to exist in the
%                 specified PCK file(s).
%
%                [1,1] = size(idcode); int32 = class(idcode)
%
%      room      the number of intervals for use as a workspace by the
%                routine. This value should equal at least the number
%                of intervals corresponding to 'idcode' in 'pck'.
%
%                [1,1] = size(room); int32 = class(room)
%
%      cover_i   an optional input describing a either an empty window or a
%                window array created from a previous cspice_pckcov call.
%                Inclusion of this window argument results in an output
%                window consisting of a union of the data retrieved from the
%                'pck' kernels and the data in 'cover_i'.
%
%                [2m,1] = size(cover_i); double = class(cover_i)
%
%                   or
%
%                [0,0] = size(cover_i); double = class(cover_i)
%
%   the call:
%
%      cover = cspice_pckcov( pck, idcode, room, cover_i )
%
%         or
%
%      cover = cspice_pckcov( pck, idcode, room )
%
%   returns:
%
%      cover   the window containing the coverage for 'idcode'', i.e.
%              the set of time intervals for which 'idcode'
%              data exist in the file 'pck'. The array 'cover'
%              contains the pairs of endpoints of these intervals.
%
%              [2p,1] = size(cover); double = class(cover)
%
%                 or
%
%              [0,1] = size(cover); double = class(cover)
%
%              Each window defined as a pair of endpoints such that:
%
%                  window 1 = cover(1:2)
%                  window 2 = cover(3:4)
%                  window 3 = cover(5:6)
%                           ...
%                  window p = cover(2p-1,2p)
%
%              The interval endpoints contained in 'cover' are ephemeris
%              times, expressed as seconds past J2000 TDB.
%
%              'cover' returns an empty set if 'pck' lacks coverage for
%              'idcode'.  If 'cover_i' exists in the argument list, 'cover'
%              returns as a union of the coverage data found in 'pck' and
%              the data in 'cover_i'. 'cover' can overwrite 'cover_i'.
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
%        % Find the set of frmects in the pck file.
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
%   the coverage a specified PCK file provides for a specified 
%   PCK class reference frame. 
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pckcov_c.
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
%   get coverage window for binary pck reference frame
%   get coverage start and stop time for binary pck frame 
%
%-&

function [cover] = cspice_pckcov( pck, idcode, room, cover_i )

   switch nargin
      case 3

         pck    = zzmice_str(pck);
         idcode = zzmice_int(idcode);
         room   = zzmice_int(room, [1, int32(inf)/2] );

      case 4

         pck    = zzmice_str(pck);
         idcode = zzmice_int(idcode);
         room   = zzmice_int(room, [1, int32(inf)/2] );
         cover_i= zzmice_win(cover_i);

      otherwise

         error ( [ 'Usage: [cover] = cspice_pckcov( _`pck`_, ' ...
                                     'idcode, room, [cover_i])' ] )

   end

   %
   % The call passed either three or four arguments. Branch accordingly.
   %
   if ( nargin == 3 )
   
      %
      % Call the MEX library.
      %
      try
         [cover] = mice('pckcov_c', pck, idcode, room );
      catch
         rethrow(lasterror)
      end
   
   else
   
      %
      % Call the MEX library.
      %
      try
         cover = mice('pckcov_c', pck, idcode, room );
      catch
         rethrow(lasterror)
      end
   
      cover = cspice_wnunid( cover, cover_i );
   
   end



