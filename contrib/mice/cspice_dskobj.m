%-Abstract
%
%   CSPICE_DSKOBJ finds the set of body ID codes of all objects for which
%   topographic data are provided in a specified DSK file.
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
%      dskfnm     SPICE DSK file name(s).
%
%                 [n,c1] = size(dskfnm); char = class(dskfnm)
%
%                    or
%
%                 [1,n] = size(dskfnm); cell = class(dskfnm)
%
%                 This file will be opened for read access by this routine.
%
%      room       a parameter specifying the maximum number of elements that
%                 can be accommodated by the dynamically allocated workspace
%                 cell used internally by this routine.
%
%                 [1,1] = size(room); int32 = class(room)
%
%                 It's not necessary to compute an accurate estimate of how
%                 many elements will be returned in `bodids'; rather, the
%                 user can pick a size considerably larger than what's
%                 really required.
%
%      bodids_i   an optional input describing an (Nx1) array of DSK
%                 ID codes.
%
%                 [r,1] = size(bodids_i); int32 = class(bodids_i)
%
%                 Inclusion of this array results in an output array
%                 consisting of a union of the data retrieved from the
%                 `dskfnm' kernels and the data in `bodids_i'.
%
%   the call:
%
%      [bodids] = cspice_dskobj( dskfnm, room, bodids_i )
%
%         or
%
%      [bodids] = cspice_dskobj( dskfnm, room )
%
%   returns:
%
%      bodids      the set of unique DSK ID codes of bodies in the indicated
%                  DSK files.
%
%                  [p,1] = size(bodids); int32 = class(bodids)
%
%                  If `bodids_i' exists in the argument list, `bodids'
%                  returns as a union of data found in `dskfnm' and the data
%                  in `bodids_i'. `bodids' can overwrite `bodids_i'.
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
%   1) Display the coverage for each object in a specified DSK file.
%      Find the set of objects in the file. Loop over the contents
%      of the ID code set: find the surface ID for each item in the
%      set and display the surface ID.
%
%
%      Example code begins here.
%
%
%      %
%      % Examine a DSK file and identify the set of
%      % central bodies associated with the segments
%      % in the file. For each body, find the
%      % set of surfaces associated with that body.
%      %
%      function dskobj_ex1( dsknam )
%
%         %
%         % Local constants
%         %
%         MAXID  = 1000;
%
%         %
%         % Find the set of objects in the DSK file.
%         %
%         bodids = cspice_dskobj( dsknam, MAXID );
%
%         for i=1:numel(bodids)
%
%            fprintf('Body ID:     %d\n'  , bodids(i) )
%
%            %
%            % Get the surface IDs for the Ith body.
%            %
%            srfids = cspice_dsksrf( dsknam, bodids(i), MAXID );
%
%            for j=1:numel(srfids)
%
%               fprintf('   Surface ID:  %d\n'  , srfids(j) )
%
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, with the following variable as input
%
%         dsknam = 'hay_a_amica_5_itokawashape_v1_0_64q.bds'
%
%      the output was:
%
%
%      Body ID:     2025143
%         Surface ID:  2025143
%
%
%   2) When Example #1 was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, with the following variable as input
%
%         dsknam = { 'hay_a_amica_5_itokawashape_v1_0_64q.bds',            ...
%                    'megr90n000eb_LL000E00N_UR090E90N_plate.bds',         ...
%                    'megr90n000eb_LL000E90S_UR090E00S_plate.bds',         ...
%                    'megr90n000eb_LL090E00N_UR180E90N_plate.bds',         ...
%                    'megr90n000eb_LL090E90S_UR180E00S_plate.bds',         ...
%                    'megr90n000eb_LL180E00N_UR270E90N_plate.bds',         ...
%                    'megr90n000eb_LL180E90S_UR270E00S_plate.bds',         ...
%                    'megr90n000eb_LL270E00N_UR360E90N_plate.bds',         ...
%                    'megr90n000eb_LL270E90S_UR360E00S_plate.bds',         ...
%                    'phobos_3_3.bds' };
%
%      the output was:
%
%
%      Body ID:     401
%         Surface ID:  401
%      Body ID:     499
%         Surface ID:  499001
%      Body ID:     2025143
%         Surface ID:  2025143
%
%
%      Note that the output lists all SPK IDs in the DSK set and all
%      corresponding surface IDs.
%
%-Particulars
%
%   This routine provides an API via which applications can determine
%   the set of objects for which there are topographic data in a
%   specified DSK file.
%
%-Exceptions
%
%   1)  If the input file has transfer format, the error
%       SPICE(INVALIDFORMAT) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the input file is not a transfer file but has architecture
%       other than DAS, the error SPICE(INVALIDARCHTYPE) is signaled
%       by a routine in the call tree of this routine.
%
%   3)  If the input file is a binary DAS file of type other than DSK,
%       the error SPICE(INVALIDFILETYPE) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If the DSK file cannot be opened or read, an error is signaled
%       by a routine in the call tree of this routine.
%
%   5)  If the size of the output set argument `bodids' is
%       insufficient to contain the actual number of ID codes of
%       objects covered by the indicated DSK file, the error
%       SPICE(CELLTOOSMALL) is signaled by a routine in the call tree
%       of this routine.
%
%   6)  If any of the input arguments, `dskfnm', `room' or `bodids_i',
%       is undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `dskfnm', `room' or `bodids_i',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See the description of the argument `dskfnm' above.
%
%-Restrictions
%
%   1)  If an error occurs while this routine is updating the set
%       `bodids', the set may be corrupted.
%
%-Required_Reading
%
%   CELLS.REQ
%   DAS.REQ
%   DSK.REQ
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
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 23-AUG-2021 (EDW) (JDR)
%
%       Changed the argument names "dsk", "size", "ids_i" and "ids" to
%       "dskfnm", "room", "bodids_i" and "bodids" for consistency with other
%       routines.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Edited the header to comply with NAIF standard.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%       Updated description of argument "room".
%
%   -Mice Version 1.0.0, 12-FEB-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   find id codes of ephemeris objects in DSK file
%   find id codes of bodies in DSK file
%
%-&

function [bodids] = cspice_dskobj( dskfnm, room, bodids_i )

   switch nargin
      case 2

         dskfnm = zzmice_str(dskfnm);
         room   = zzmice_int(room, [1, int32(inf)/2] );

      case 3

         dskfnm   = zzmice_str(dskfnm);
         room     = zzmice_int(room, [1, int32(inf)/2] );
         bodids_i = zzmice_cell( bodids_i, 'int32' );

      otherwise

         error( ['Usage: [bodids] = '                                      ...
                        'cspice_dskobj( _`dskfnm`_, room, [bodids_i])'] )

   end

%
% The call passed either two or three arguments. Branch accordingly.
%
if ( nargin == 2 )

   %
   % Call the MEX library.
   %
   try
      bodids = mice('dskobj_c', dskfnm, room );
   catch spiceerr
      rethrow(spiceerr)
   end

else

   %
   % Call the MEX library.
   %
   try
      bodids = unique( [ [bodids_i]; mice('dskobj_c', dskfnm, room ) ] );
   catch spiceerr
      rethrow(spiceerr)
   end


end
