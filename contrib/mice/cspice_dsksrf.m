%-Abstract
%
%   CSPICE_DSKSRF finds the set of surface ID codes for all surfaces
%   associated with a given body in a specified DSK file.
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
%      dsk       SPICE DSK files names.
%
%                [1,c] = size(dsk), char = class(dsk)
%
%                   or
%
%                [1,m] = size(dsk), cell = class(dsk)
%
%      idcode    the SPK ID code of a body for which topographic
%                data are present in the specified DSK files.
%
%                [1,1] = size(idcode); int32 = class(idcode)
%
%      room      the maximum number of NAIF IDs to return from 'dsk'.
%
%                [1,1] = size(room), int32 = class(room)
%
%      srfids_i  an optional input describing an array of DSK
%                ID codes. Inclusion of this array results in an output
%                array consisting of a union of the data retrieved from
%                the 'dsk' kernels and the data in 'srfids_i'.
%
%                [n,1] = size(srfids_i), int32 = class(srfids_i)
%
%                   or
%
%                [0,0] = size(srfids_i), int32 = class(srfids_i)
%
%   the call:
%
%      srfids = cspice_dsksrf( dsk, idcode, room, srfids_i )
%
%         or
%
%      srfids = cspice_dsksrf( dsk, idcode, room )
%
%   returns:
%
%      srfids   an array containing the union of `srfids_i'
%               and the ID codes of the surfaces associated with the
%               body designated by `bodyid', for which segments were
%               found in the indicated DSK file.
%
%               [p,1] = size(srfids), int32 = class(srfids)
%
%               The elements of `srfids' are unique; each ID
%               code in `srfids' appears only once, even if the DSK
%               file contains multiple segments for that ID code.
%
%               See the Examples section below for a complete
%               example program showing how to retrieve body and
%               surface ID codes from a DSK file.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Use a simple function to display the DSK IDs found in a DSK or set of
%   DSKs, and the time coverage of the data corresponding to those IDs.
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
%   Matlab outputs:
%
%     >> dsksrf_t( 'hay_a_amica_5_itokawashape_v1_0_64q.bds' )
%     Body ID:     2025143
%        Surface ID:  2025143
%
%     The output lists the SPK IDs in the DSK and the surface IDs.
%
%   Example (2):
%
%   Matlab outputs:
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
%     >> dsksrf_t(dsk)
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
%   This routine provides an API via which applications can determine
%   the set of surfaces associated with a given body in a specified
%   DSK file. This routine is normally used together with cspice_dskobj.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dsksrf_c.
%
%   MICE.REQ
%   DSK.REQ
%   DAS.REQ
%   SETS.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 18-FEB-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   get coverage window for dsk object
%
%-&

function [srfids] = cspice_dsksrf( dsk, bodyid, room, srfids_i )

   switch nargin
      case 3

         dsk    = zzmice_str(dsk);
         bodyid = zzmice_int(bodyid);
         room   = zzmice_int(room, [1, int32(inf)/2] );

      case 4

         dsk    = zzmice_str(dsk);
         bodyid = zzmice_int(bodyid);
         room   = zzmice_int(room, [1, int32(inf)/2] );
         srfids_i= zzmice_cell(srfids_i);

      otherwise

         error ( [ 'Usage: [srfids] = cspice_dsksrf( _`dsk`_, ' ...
                                     'bodyid, room, srfids_i])' ] )

   end

%
% The call passed either three or four arguments. Branch accordingly.
%
if ( nargin == 3 )

   %
   % Call the MEX library.
   %
   try
      [srfids] = mice('dsksrf_c', dsk, bodyid, room );
   catch
      rethrow(lasterror)
   end

else

   %
   % Call the MEX library.
   %
   try
      [srfids] = unique( [ [srfids_i];  mice('dsksrf_c', dsk, bodyid, room) ]);
   catch
      rethrow(lasterror)
   end

end


