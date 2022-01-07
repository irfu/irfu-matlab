%-Abstract
%
%   CSPICE_DASCLS closes a DAS file.
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
%      handle   the file handle of an open DAS file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_dascls( handle )
%
%   returns:
%
%      None.
%
%      See -Particulars for a description of the effect of this routine.
%
%-Parameters
%
%   All parameters described here are declared in the Mice
%   include file MiceDAS.m. See that file for parameter values.
%
%   SPICE_DAS_FTSIZE
%
%               is the maximum number of DAS files that can be
%               open at any one time.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Dump several parameters from the first DLA segment of
%      a DSK file. The segment is assumed to be of type 2.
%
%      Example code begins here.
%
%
%      function dascls_ex1( )
%
%         %
%         % Prompt for the name of the file to search.
%         %
%         fname = input( 'Name of DSK file > ', 's' );
%
%         %
%         % Open the DSK file for read access. We use the
%         % DAS-level interface for this function.
%         %
%         handle = cspice_dasopr( fname );
%
%         %
%         % Begin a forward search through the kernel. In
%         % this example, it's a very short search.
%         %
%         [ dladsc, found] = cspice_dlabfs( handle );
%
%         if ~found
%
%            txt = sprintf( ...
%               'SPICE(NOSEGMENT): No segment found in file %s.', ...
%               dsk);
%
%            error(txt)
%
%         end
%
%         %
%         % Loop over each segment.
%         %
%         while found
%
%            %
%            % If we made it this far, DLADSC is the
%            % DLA descriptor of the first segment.
%            % Read and display type 2 bookkeeping data.
%            %
%            [ nv, np, nvxtot, vtxbds, voxsiz, voxori, vgrext, ...
%            cgscal, vtxnpl, voxnpt, voxnpl] = cspice_dskb02( handle, dladsc);
%
%            fprintf( ['\n'                                              ...
%                     'Number of vertices:                 %ld\n'        ...
%                     'Number of plates:                   %ld\n'        ...
%                     'Number of voxels:                   %ld\n'        ...
%                     'Vertex bounds in X direction (km):  %f : %f\n'    ...
%                     'Vertex bounds in Y direction (km):  %f : %f\n'    ...
%                     'Vertex bounds in Z direction (km):  %f : %f\n'    ...
%                     'Voxel edge length (km):             %f\n'         ...
%                     'Voxel grid origin (km):           ( %f %f %f )\n' ...
%                     'Voxel grid extents:                 %ld %ld %ld\n'...
%                     'Coarse voxel grid scale:            %ld\n'        ...
%                     'Size of vertex-plate list:          %ld\n'        ...
%                     'Size of voxel-plate pointer array:  %ld\n'        ...
%                     'Size of voxel-plate list:           %ld\n'],      ...
%                     nv,                              ...
%                     np,                              ...
%                     nvxtot,                          ...
%                     vtxbds(1,1), vtxbds(2,1),        ...
%                     vtxbds(1,2), vtxbds(2,2),        ...
%                     vtxbds(1,3), vtxbds(2,3),        ...
%                     voxsiz,                          ...
%                     voxori(1), voxori(2), voxori(3), ...
%                     vgrext(1), vgrext(2), vgrext(3), ...
%                     cgscal,                          ...
%                     vtxnpl,                          ...
%                     voxnpt,                          ...
%                     voxnpl )
%
%               %
%               % Search for the segment after that described by `dladsc'.
%               % `found' returns as false if no such segment located.
%               %
%               [nxtdsc, found] = cspice_dlafns( handle, dladsc);
%
%               dladsc = nxtdsc;
%
%            end
%
%         %
%         % Close the kernel. This frees program and system resources.
%         %
%         cspice_dascls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the DSK file named phobos512.bds, the output
%      was:
%
%
%      Name of DSK file > phobos512.bds
%
%      Number of vertices:                 1579014
%      Number of plates:                   3145728
%      Number of voxels:                   11914500
%      Vertex bounds in X direction (km):  -13.440030 : 12.762800
%      Vertex bounds in Y direction (km):  -11.520650 : 12.061140
%      Vertex bounds in Z direction (km):  -9.570780 : 10.055000
%      Voxel edge length (km):             0.104248
%      Voxel grid origin (km):           ( -14.073520 -11.988554 -9.903588 )
%      Voxel grid extents:                 260 235 195
%      Coarse voxel grid scale:            5
%      Size of vertex-plate list:          11010050
%      Size of voxel-plate pointer array:  1151500
%      Size of voxel-plate list:           6419540
%
%
%-Particulars
%
%   This routine provides the primary recommended method of closing an
%   open DAS file.
%
%-Exceptions
%
%   1)  If `handle' is not the handle of an open DAS file, no error
%       is signaled.
%
%   2)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `handle' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   See the description of input argument `handle' in -I/O.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       complete example code based on example in cspice_dskb02.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%   -Mice Version 1.0.0, 28-APR-2016 (NJB) (EDW)
%
%-Index_Entries
%
%   close a DAS file
%
%-&

function cspice_dascls( handle )

   switch nargin
      case 1

         handle  = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dascls( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dascls_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end
