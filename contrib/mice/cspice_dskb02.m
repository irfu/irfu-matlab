%-Abstract
%
%   CSPICE_DSKB02 return bookkeeping data from a DSK type 2 segment.
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
%      handle    the handle of a DSK file containing a type 2
%                segment from which data are to be fetched.
%
%                [1,1] = size(handle); int32 = class(handle)
%
%      dladsc    the DLA descriptor associated with the segment
%                from which data are to be fetched.
%
%                [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                int32 = class(dladsc)
%
%   the call:
%
%
%   returns:
%
%      nv         the number of vertices belonging to the specified plate
%                 model.
%
%                 [1,1] = size(nv); int32 = class(nv)
%
%      np         the number of plates belonging to the specified plate
%                 model.
%
%                 [1,1] = size(np); int32 = class(np)
%
%      nvxtot     Total number of voxels in fine grid.
%
%                 [1,1] = size(nvxtot); int32 = class(nvxtot)
%
%      vtxbds     Vertex bounds. This is an array of six values giving
%                 the minimum and maximum values of each component of
%                 the vertex set. Units are km.
%
%                 [2,3] = size(vtxbds); double = class(vtxbds)
%
%      voxsiz     Fine grid voxel size.  DSK voxels are cubes; the edge
%                 length of each cube is given by the voxel size.  This
%                 size applies to the fine voxel grid.  Units are km.
%
%                 [1,1] = size(voxsiz); double = class(voxsiz)
%
%      voxori     Voxel grid origin.  This is the location of the voxel
%                 grid origin in the body-fixed frame associated with
%                 the target body.  Units are km.
%
%                 [3,1] = size(voxori); double = class(voxori)
%
%      vgrext     Voxel grid extent.  This extent is an array of three
%                 integers indicating the number of voxels in the X, Y, and
%                 Z directions in the fine voxel grid.
%
%                 [3,1] = size(vgrext); int32 = class(vgrext)
%
%      cgscal     Coarse voxel grid scale.  The extent of the fine voxel
%                 grid is related to the extent of the coarse voxel grid by
%                 this scale factor.
%
%                 [1,1] = size(cgscal); int32 = class(cgscal)
%
%      vtxnpl     Vertex-plate correspondence list size.
%
%                 [1,1] = size(vtxnpl); int32 = class(vtxnpl)
%
%      voxnpt     Size of the voxel-to-plate pointer list.
%
%                 [1,1] = size(voxnpt); int32 = class(voxnpt)
%
%      voxnpl     Voxel-plate correspondence list size.
%
%                 [1,1] = size(voxnpl); int32 = class(voxnpl)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Dump several parameters from the first DLA segment of
%      a DSK file. The segment is assumed to be of type 2.
%
%      function dskb02_t( dsk )
%
%         %
%         % Open the DSK file for read access. We use the
%         % DAS-level interface for this function.
%         %
%         handle = cspice_dasopr( dsk );
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
%               % Search for the segment after that described by 'dladsc'.
%               % 'found' returns as false if no such segment located.
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
%   MATLAB outputs:
%
%      >> dskb02_t( 'phobos512.bds' )
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
%-Particulars
%
%   This routine supports computations involving bookkeeping information
%   stored in DSK type 2 segments. User applications typically will not
%   need to call this routine.
%
%   DSK files are built using the DLA low-level format and the DAS
%   architecture; DLA files are a specialized type of DAS file in which
%   data are organized as a doubly linked list of segments.  Each
%   segment's data belong to contiguous components of character, double
%   precision, and integer type.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine dskb02_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 27-APR-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   fetch parameters from a type 2 dsk segment
%
%-&

function [ nv, np, nvxtot, vtxbds, voxsiz, voxori,    ...
           vgrext, cgscal, vtxnpl, voxnpt, voxnpl ] = ...
           cspice_dskb02( handle, dladsc )

   switch nargin
      case 2

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );

      otherwise

         error ( [ '[ nv, np, '                 ...
                   'nvxtot, vtxbds(2,3), '      ...
                   'voxsiz, voxori(3,1), '      ...
                   'vgrext(3,1), cgscal, '      ...
                   'vtxnpl, voxnpt, voxnpl ] = '...
                   'cspice_dskb02( handle, dladsc )' ] )
   end

   %
   % Call the MEX library.
   %
   try

      [ nv, np, nvxtot, vtxbds, voxsiz, voxori,    ...
      vgrext, cgscal, vtxnpl, voxnpt, voxnpl ] = ...
      mice( 'dskb02_c', handle, dladsc  );
   catch
      rethrow(lasterror)
   end
