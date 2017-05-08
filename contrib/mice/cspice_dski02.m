%-Abstract
%
%   CSPICE_DSKI02 returns integer data from a type 2 DSK segment.
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
%      handle    the handle of a DSK file containing a type 2 segment
%                from which data are to be fetched.
%
%                [1,1] = size(handle); int32 = class(handle)
%
%      dladsc    the DLA descriptor associated with the segment from
%                which data are to be fetched.
%
%                [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                int32 = class(dladsc)
%
%      item      an integer "keyword" parameter designating the integer
%                data item to fetch.
%
%                [1,1] = size(item); int32 = class(item)
%
%                Names, values, and meanings of keyword parameters
%                supported by this routine are shown below.
%
%                Use of the names shown here is enabled by calling
%                the DSKMice parameter definition routine as shown:
%
%                   DSKMiceUser
%
%                This call must be made before the parameter names
%                are referenced. See the example program below.
%
%
%             Name               Value   Description
%             ----               -----   ----------
%
%             SPICE_DSK02_KWNV     1     Number of vertices in model.
%
%             SPICE_DSK02_KWNP     2     Number of plates in model.
%
%             SPICE_DSK02_KWNVXT   3     Total number of voxels in fine grid.
%
%             SPICE_DSK02_KWVGRX   4     Voxel grid extent.  This extent is
%                                        an array of three integers
%                                        indicating the number of
%                                        voxels in the X, Y, and Z
%                                        directions in the fine voxel
%                                        grid.
%
%             SPICE_DSK02_KWCGSC   5     Coarse voxel grid scale.  The extent
%                                        of the fine voxel grid is related to
%                                        the extent of the coarse voxel grid
%                                        by this scale factor.
%
%             SPICE_DSK02_KWVXPS   6     Size of the voxel-to-plate pointer
%                                        list.
%
%             SPICE_DSK02_KWVXLS   7     Voxel-plate correspondence list size.
%
%             SPICE_DSK02_KWVTLS   8     Vertex-plate correspondence list
%                                        size.
%
%             SPICE_DSK02_KWPLAT   9     Plate array.  For each plate, this
%                                        array contains the indices of the
%                                        plate's three vertices.  The ordering
%                                        of the array members is:
%
%                                           Plate 1 vertex index 1
%                                           Plate 1 vertex index 2
%                                           Plate 1 vertex index 3
%                                           Plate 2 vertex index 1
%                                           ...
%
%                                        The vertex indices in this
%                                        array start at 1 and end at
%                                        NV, the number of vertices
%                                        in the model.
%
%             SPICE_DSK02_KWVXPT   10    Voxel-plate pointer list. This list
%                                        contains pointers that map fine
%                                        voxels to lists of plates that
%                                        intersect those voxels. Note that
%                                        only fine voxels belonging to
%                                        non-empty coarse voxels are in the
%                                        domain of this mapping.
%
%             SPICE_DSK02_KWVXPL   11    Voxel-plate correspondence list.
%                                        This list contains lists of plates
%                                        that intersect fine voxels. (This
%                                        list is the data structure into
%                                        which the voxel-to-plate pointers
%                                        point.) This list can contain
%                                        empty lists. Plate IDs in this
%                                        list start at 1 and end at NP,
%                                        the number of plates in the model.
%
%             SPICE_DSK02_KWVTPT   12    Vertex-plate pointer list. This list
%                                        contains pointers that map vertices
%                                        to lists of plates to which those
%                                        vertices belong.
%
%                                        Note that the size of this list is
%                                        always NV, the number of vertices.
%                                        Hence there's no need for a separate
%                                        keyword for the size of this list.
%
%             SPICE_DSK02_KWVTPL   13    Vertex-plate correspondence list.
%                                        This list contains, for each vertex,
%                                        the indices of the plates to which
%                                        that  vertex belongs. Plate IDs in
%                                        this list start at 1 and end at NP,
%                                        the number of plates in the model.
%
%             SPICE_DSK02_KWCGPT   14    Coarse voxel grid pointers.  This is
%                                        an array of pointers mapping coarse
%                                        voxels to lists of pointers in the
%                                        voxel-plate pointer list.  Each
%                                        non-empty coarse voxel maps to a
%                                        list of pointers; every fine voxel
%                                        contained in a non-empty coarse voxel
%                                        has its own pointers. Grid elements
%                                        corresponding to empty coarse voxels
%                                        contain non-positive values.
%
%      start     the start index within specified data item from which
%                data are to be fetched. The index of the first element
%                of each data item is 1. This convention applies
%                uniformly to all data. For example, the plate
%                ID range starts at 1 (this fact is language-independent),
%                but a caller would use a `start' value of 1 to fetch
%                the vertex indices of the first plate.
%
%                [1,1] = size(start); int32 = class(start)
%
%      room      the amount of room in the output array. It is
%                permissible to provide an output array that has too
%                little room to fetch an item in one call.
%
%                [1,1] = size(room); int32 = class(room)
%
%   the call:
%
%      [values] = cspice_dski02( handle, dladsc, item, start, room )
%
%   returns:
%
%      values    a contiguous set of elements of the item designated by
%                `item'. The correspondence of `values' with the
%                elements of the data item is:
%
%                   values(1)      item(start)
%                     ...             ...
%                   values(n)      item(start+n-1)
%
%                If an error occurs on the call, `values' is undefined.
%
%                [1,N] = size(values); int32 = class(values)
%
%                Note, `room' >= N.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Look up all the vertices associated with each plate
%      of the model contained in a specified type 2 segment.
%      For each plate, display the plate's vertices.
%
%      For this example, we'll show the context of this look-up:
%      opening the DSK file for read access, traversing a trivial,
%      one-segment list to obtain the segment of interest.
%
%
%      function dski02_t( dsk )
%
%         %
%         % MiceUser globally defines DSK parameters.
%         % For more information, please see DSKMiceUser.m and
%         % DSKMice02.m.
%         %
%         MiceUser
%
%         %
%         % Set the dimensions of the array `vrtces', which
%         % will be used later.
%         %
%         vrtces = zeros(3,3);
%
%         %
%         % Open the DSK file for read access.
%         % We use the DAS-level interface for
%         % this function.
%         %
%         handle  = cspice_dasopr( dsk );
%
%         %
%         % Begin a forward search through the
%         % kernel, treating the file as a DLA.
%         % In this example, it's a very short
%         % search.
%         %
%         [dladsc, found] = cspice_dlabfs( handle );
%
%         if ~found
%
%            %
%            % We arrive here only if the kernel
%            % contains no segments. This is
%            % unexpected, but we're prepared for it.
%            %
%            fprintf( 'No segments found in DSK file %s\n', dsk )
%            return
%
%         end
%
%         %
%         % If we made it this far, `dladsc' is the
%         % DLA descriptor of the first segment.
%         %
%
%         %
%         % Find the number of plates in the model.
%         %
%         ival = cspice_dski02( handle, dladsc, SPICE_DSK02_KWNP, 1, 1 );
%         np = ival(1);
%
%         %
%         % For each plate, look up the desired data.
%         % Note that plate numbers range from 1 to np.
%         %
%
%         for  i = 1:np
%
%            %
%            % For the Ith plate, find the associated
%            % vertex IDs.  We must take into account
%            % the fact that each plate has three
%            % vertices when we compute the start
%            % index.
%            %
%
%            start = 3*(i-1) + 1;
%
%            %
%            % Fetch the ith plate.
%            %
%            vrtids = cspice_dski02( handle, dladsc, SPICE_DSK02_KWPLAT, ...
%                                    start,  3 );
%
%            for  j = 1:3
%
%               %
%               % Fetch the jth vertex of the ith plate.
%               %
%               start = (vrtids(j)-1) * 3 +1;
%
%               vtemp = cspice_dskd02( handle, dladsc, SPICE_DSK02_KWVERT, ...
%                                      start,  3 );
%
%               vrtces(j,:) = vtemp;
%
%            end
%
%            %
%            % Display the vertices of the ith plate:
%            %
%            fprintf( '\n' )
%            fprintf( 'Plate number: %d\n', i )
%
%            for  j = 1:3
%                 fprintf( 'Vertex %d: (%14.6e %14.6e %14.6e)\n', ...
%                                                j, vrtces(j,:) )
%            end
%
%         end
%
%         %
%         % Close the DSK.
%         %
%         cspice_dascls( handle );
%
%
%   MATLAB outputs:
%
%      dskd02_t( 'solid.bds' )
%
%            [Only the first and last few rows are shown]
%
%      Plate number: 1
%      Vertex 1: (  0.000000e+00   0.000000e+00   1.175570e+00)
%      Vertex 2: (  1.051460e+00   0.000000e+00   5.257310e-01)
%      Vertex 3: (  3.249200e-01   1.000000e+00   5.257310e-01)
%
%      Plate number: 2
%      Vertex 1: (  0.000000e+00   0.000000e+00   1.175570e+00)
%      Vertex 2: (  3.249200e-01   1.000000e+00   5.257310e-01)
%      Vertex 3: ( -8.506510e-01   6.180340e-01   5.257310e-01)
%
%                 ...
%
%      Plate number: 19
%      Vertex 1: ( -3.249200e-01  -1.000000e+00  -5.257310e-01)
%      Vertex 2: (  0.000000e+00   0.000000e+00  -1.175570e+00)
%      Vertex 3: (  8.506510e-01  -6.180340e-01  -5.257310e-01)
%
%      Plate number: 20
%      Vertex 1: (  8.506510e-01  -6.180340e-01  -5.257310e-01)
%      Vertex 2: (  0.000000e+00   0.000000e+00  -1.175570e+00)
%      Vertex 3: (  8.506510e-01   6.180340e-01  -5.257310e-01)
%
%-Particulars
%
%   Most user applications will not need to call this routine. The
%   routines
%
%      cspice_dskz02
%      cspice_dskp02
%      cspice_dskv02
%
%   have simpler interfaces and may be used to fetch the plates
%   and vertex counts, and the plates and vertices themselves,
%   from a type DSK segment. See the documentation of those
%   routines for code examples.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine dski02_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 28-NOV-2016, NJB (JPL), EDW (JPL)
%
%-Index_Entries
%
%   fetch integer data from a type 2 dsk segment
%
%-&

function [values] = cspice_dski02( handle, dladsc, item, start, room )

   switch nargin
      case 5

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         item   = zzmice_int( item   );
         start  = zzmice_int( start  );
         room   = zzmice_int( room, [1, int32(inf)/2] );

      otherwise

         error ( [ 'Usage: [values] = ' ...
                   'cspice_dski02( handle, dladsc, item, start, room ) ' ] )
   end

   %
   % Call the MEX library.
   %
   try

      [values] = mice( 'dski02_c', ...
                       handle, dladsc, item, start, room );
   catch
      rethrow(lasterror)
   end


