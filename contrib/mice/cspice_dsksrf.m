%-Abstract
%
%   CSPICE_DSKSRF finds the set of surface ID codes for all surfaces
%   associated with a given body in a specified DSK file.
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
%      dskfnm     the SPICE DSK file name(s).
%
%                 [n,c1] = size(dskfnm); char = class(dskfnm)
%
%                   or
%
%                 [1,n] = size(dskfnm); cell = class(dskfnm)
%
%                 This file will be opened for read access by this routine.
%
%      bodyid     the integer ID code of a body for which topographic data
%                 are present in the specified DSK file.
%
%                 [1,1] = size(bodyid); int32 = class(bodyid)
%
%      room       a parameter specifying the maximum number of elements that
%                 can be accommodated by the dynamically allocated workspace
%                 cell used internally by this routine.
%
%                 [1,1] = size(room); int32 = class(room)
%
%                 It's not necessary to compute an accurate estimate of how
%                 many elements will be returned in `srfids'; rather, the
%                 user can pick a size considerably larger than what's
%                 really required.
%
%      srfids_i   an optional input describing an array of DSK ID codes.
%
%                 [r,1] = size(srfids_i); int32 = class(srfids_i)
%
%                 `srfids_i' optionally may contain a set of surface ID codes
%                 on input; on output, the ID codes already present in
%                 `srfids_i' will be combined with surface ID code set found
%                 for the body designated by `bodyid' in the file `dskfnm'.
%
%   the call:
%
%      [srfids] = cspice_dsksrf( dskfnm, bodyid, room, srfids_i )
%
%         or
%
%      [srfids] = cspice_dsksrf( dskfnm, bodyid, room )
%
%   returns:
%
%      srfids     an array containing the union of `srfids_i'
%                 and the ID codes of the surfaces associated with the
%                 body designated by `bodyid', for which segments were
%                 found in the indicated DSK file.
%
%                 [p,1] = size(srfids), int32 = class(srfids)
%
%                 The elements of `srfids' are unique; each ID
%                 code in `srfids' appears only once, even if the DSK
%                 file contains multiple segments for that ID code.
%
%                 See the -Examples section below for a complete
%                 example program showing how to retrieve body and
%                 surface ID codes from a DSK file.
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
%   1) Examine a DSK file and identify the set of central bodies
%      associated with the segments in the file. For each body, find
%      the set of surfaces associated with that body.
%
%
%      Example code begins here.
%
%
%      function dsksrf_ex1()
%
%         MAXID  = 1000;
%
%         %
%         % Prompt for the name of the file to search.
%         %
%         dskfnm = input( 'Name of DSK file > ', 's' );
%
%         %
%         % Find the set of objects in the DSK file.
%         %
%         bodids = cspice_dskobj( dskfnm, MAXID );
%
%         for i=1:numel(bodids)
%
%            fprintf('Body ID:     %d\n'  , bodids(i) )
%
%            %
%            % Get the surface IDs for the Ith body.
%            %
%            srfids = cspice_dsksrf( dskfnm, bodids(i), MAXID );
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
%      platform, using the DSK file named phobos512.bds, the output
%      was:
%
%
%      Name of DSK file > phobos512.bds
%      Body ID:     401
%         Surface ID:  401
%
%
%   2) Update the previous example to assign a DSK kernel list as
%      input (instead of prompting for it), containing shape data
%      for Itokawa (NAIF ID 2025143), Mars (499) and Phobos (401).
%
%
%      Use the DSK kernel below to provide data for Itokawa.
%
%         hay_a_amica_5_itokawashape_v1_0_64q.bds
%
%      Use the DSK kernel below to provide data for Mars.
%
%         megr90n000eb_LL000E00N_UR090E90N_plate.bds
%
%      Use the DSK kernel below to provide data for Phobos.
%
%         phobos_3_3.bds
%
%
%      Example code begins here.
%
%
%      function dsksrf_ex2()
%
%         MAXID  = 1000;
%
%         %
%         % Assing the DSK kernel list.
%         %
%         dskfnm = { 'hay_a_amica_5_itokawashape_v1_0_64q.bds',    ...
%                    'megr90n000eb_LL000E00N_UR090E90N_plate.bds', ...
%                    'phobos_3_3.bds' };
%
%         %
%         % Find the set of objects in the DSK file.
%         %
%         bodids = cspice_dskobj( dskfnm, MAXID );
%
%         for i=1:numel(bodids)
%
%            fprintf('Body ID:     %d\n'  , bodids(i) )
%
%            %
%            % Get the surface IDs for the Ith body.
%            %
%            srfids = cspice_dsksrf( dskfnm, bodids(i), MAXID );
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
%      platform, the output was:
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
%-Particulars
%
%   This routine provides an API via which applications can determine
%   the set of surfaces associated with a given body in a specified
%   DSK file. This routine is normally used together with cspice_dskobj.
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
%   5)  If the size of the output set argument `srfids' is
%       insufficient to contain the actual number of ID codes of
%       surfaces covered by the indicated DSK file, the error
%       SPICE(CELLTOOSMALL) is signaled by a routine in the call tree
%       of this routine.
%
%   6)  If any of the input arguments, `dskfnm', `bodyid', `room' or
%       `srfids_i', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   7)  If any of the input arguments, `dskfnm', `bodyid', `room' or
%       `srfids_i', is not of the expected type, or it does not have
%       the expected dimensions and size, an error is signaled by the
%       Mice interface.
%
%-Files
%
%   See the description of the argument `dsk' above.
%
%-Restrictions
%
%   1)  If an error occurs while this routine is updating the set
%       `srfids', the set may be corrupted.
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Fixed bug: modified the zzmice_cell call to include the cell type
%       identifier "int32".
%
%       Changed input argument names "dsk" and "idcode" to "dskfnm" and
%       "bodyid" for consistency with other routines.
%
%       Edited the header to comply with NAIF standard. Updated
%       code example #1 to prompt for the input DSK file and code example
%       #2 to hardcode the input DSK kernel list.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Corrected minor typo in Usage error message.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%       Updated description of argument "room".
%
%   -Mice Version 1.0.0, 18-FEB-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   find id codes of surfaces in DSK file
%
%-&

function [srfids] = cspice_dsksrf( dskfnm, bodyid, room, srfids_i )

   switch nargin
      case 3

         dskfnm   = zzmice_str(dskfnm);
         bodyid   = zzmice_int(bodyid);
         room     = zzmice_int(room, [1, int32(inf)/2] );

      case 4

         dskfnm   = zzmice_str(dskfnm);
         bodyid   = zzmice_int(bodyid);
         room     = zzmice_int(room, [1, int32(inf)/2] );
         srfids_i = zzmice_cell(srfids_i, 'int32');

      otherwise

         error ( [ 'Usage: [srfids] = cspice_dsksrf( _`dskfnm`_, '         ...
                                     'bodyid, room, [srfids_i])' ] )

   end

%
% The call passed either three or four arguments. Branch accordingly.
%
if ( nargin == 3 )

   %
   % Call the MEX library.
   %
   try
      [srfids] = mice('dsksrf_c', dskfnm, bodyid, room );
   catch spiceerr
      rethrow(spiceerr)
   end

else

   %
   % Call the MEX library.
   %
   try
      [srfids] = unique( [ [srfids_i];                                     ...
                         mice('dsksrf_c', dskfnm, bodyid, room) ]);
   catch spiceerr
      rethrow(spiceerr)
   end

end
