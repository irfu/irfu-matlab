%-Abstract
%
%   CSPICE_DSKOBJ returns the set of body ID codes of all objects
%   for which topographic data are provided in specified DSK files.
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
%      dsk     SPICE dsk files names.
%
%              [1,c1] = size(dsk), char = class(dsk)
%
%                 or
%
%              [1,m] = size(dsk), cell = class(dsk)
%
%      room    the maximum number of DSK IDs to return from 'dsk'.
%
%              [1,1] = size(room), int32 = class(room)
%
%      ids_i   an optional input describing an (Nx1) array of DSK
%              ID codes. Inclusion of this array results in an output
%              array consisting of a union of the data retrieved from
%              the 'dsk' kernels and the data in 'ids_i'.
%
%              [n,1] = size(ids_i), int32 = class(ids_i)
%
%                 or
%
%              [0,0] = size(ids_i), int32 = class(ids_i)
%
%   the call:
%
%      ids = cspice_dskobj( dsk, room, ids_i)
%
%         or
%
%      ids = cspice_dskobj( dsk, room)
%
%   returns:
%
%      ids   the set of unique DSK ID codes of segments in the indicated DSK
%            files. If 'ids_i'  exists in the argument list, 'ids' returns
%            as a union of data found in 'dsk' and the data in 'ids_i'.
%            'ids' can overwrite 'ids_i'.
%
%            [p,1] = size(ids), int32 = class(ids)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example (1):
%
%      Examine a DSK file and identify the set of
%      central bodies associated with the segments
%      in the file. For each body, find the
%      set of surfaces associated with that body.
%
%      function dsksrf_t(DSK)
%
%            MAXID  = 1000;
%
%            %
%            % Find the set of objects in the DSK file.
%            %
%            bodids = cspice_dskobj( DSK, MAXID );
%
%            for i=1:numel(bodids)
%
%               fprintf('Body ID:     %d\n'  , bodids(i) )
%
%               %
%               % Get the surface IDs for the Ith body.
%               %
%               srfids = cspice_dsksrf( DSK, bodids(i), MAXID );
%
%               for j=1:numel(srfids)
%
%                  fprintf('   Surface ID:  %d\n'  , srfids(j) )
%
%               end
%
%            end
%
%
%     Assign a DSK kernel as:
%
%     dsk =  'hay_a_amica_5_itokawashape_v1_0_64q.bds';
%
%     >> dsksrf_t( dsk )
%     Body ID:     2025143
%        Surface ID:  2025143
%
%     The output lists the SPK IDs in the DSK and the surface IDs.
%
%  Example (2):
%
%     Assign a DSK kernel list as:
%
%     >> dsk= { 'hay_a_amica_5_itokawashape_v1_0_64q.bds',  ...
%             'megr90n000eb_LL000E00N_UR090E90N_plate.bds', ...
%             'megr90n000eb_LL000E90S_UR090E00S_plate.bds', ...
%             'megr90n000eb_LL090E00N_UR180E90N_plate.bds', ...
%             'megr90n000eb_LL090E90S_UR180E00S_plate.bds', ...
%             'megr90n000eb_LL180E00N_UR270E90N_plate.bds', ...
%             'megr90n000eb_LL180E90S_UR270E00S_plate.bds', ...
%             'megr90n000eb_LL270E00N_UR360E90N_plate.bds', ...
%             'megr90n000eb_LL270E90S_UR360E00S_plate.bds', ...
%             'phobos_3_3.bds' };
%
%     >> dsksrf_t_t(dsk)
%     Body ID:     401
%        Surface ID:  401
%     Body ID:     499
%        Surface ID:  499001
%     Body ID:     2025143
%        Surface ID:  2025143
%
%     The output lists all SPK IDs in the DSK set and all corresponding
%     surface IDs.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dskobj_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%   SETS.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-FEB-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   find id codes of ephemeris objects in dsk file
%   find id codes of bodies in dsk file
%
%-&

function [ids] = cspice_dskobj( dsk, room, ids_i )

   switch nargin
      case 2

         dsk  = zzmice_str(dsk);
         room = zzmice_int(room, [1, int32(inf)/2] );

      case 3

         dsk  = zzmice_str(dsk);
         room = zzmice_int(room, [1, int32(inf)/2] );
         ids_i= zzmice_cell( ids_i, 'int32' );

      otherwise

         error( 'Usage: [ids] = cspice_dskobj( _`dsk`_, room, [ids_i])' )

   end

%
% The call passed either two or three arguments. Branch accordingly.
%
if ( nargin == 2 )

   %
   % Call the MEX library.
   %
   try
      ids = mice('dskobj_c', dsk, room );
   catch
      rethrow(lasterror)
   end

else

   %
   % Call the MEX library.
   %
   try
      ids = unique( [ [ids_i]; mice('dskobj_c', dsk, room ) ]  );
   catch
      rethrow(lasterror)
   end


end

