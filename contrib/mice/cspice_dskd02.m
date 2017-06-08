%-Abstract
%
%   CSPICE_DSKD02 returns double precision data from a type 2 DSK segment.
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
%      handle    is the handle of a DSK file containing a type 2
%                segment from which data are to be fetched.
%
%                [1,1] = size(handle); int32 = class(handle)
%
%      dladsc    is the DLA descriptor associated with the segment from
%                which data are to be fetched.
%
%                [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                int32 = class(dladsc)
%
%      item      is an integer "keyword" parameter designating the double
%                precision data item to fetch.
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
%             Name               Value   Description
%             ----               -----   ----------
%
%             SPICE_DSK02_KWDSC   15     Array containing contents of Fortran
%                                        DSK descriptor of segment. Note
%                                        that DSK descriptors are not to be
%                                        confused with DLA descriptors, which
%                                        contain segment component base
%                                        address and size information. The
%                                        dimension of this array is
%                                        SPICE_DSK_DSCSIZ.
%
%             SPICE_DSK02_KWVTBD  16     Vertex bounds. This is an array of
%                                        six values giving the minimum and
%                                        maximum values of each component
%                                        of the vertex set.
%
%             SPICE_DSK02_KWVXOR  17     Voxel grid origin. This is the
%                                        location of the voxel grid origin in
%                                        the body-fixed frame associated with
%                                        the target body.
%
%             SPICE_DSK02_KWVXSZ  18     Voxel size. DSK voxels are cubes;
%                                        the edge length of each cube is
%                                        given by the voxel size. This size
%                                        applies to the fine voxel grid. Units
%                                        are km.
%
%             SPICE_DSK02_KWVERT  19     Vertex coordinates.
%
%      start     is the start index within the specified data item from which
%                data are to be fetched. The index of the first element
%                of each data item is 1. This convention applies
%                uniformly to all data.
%
%                [1,1] = size(start); int32 = class(start)
%
%      room      is the amount of room in the output array. It is
%                permissible to provide an output array that has too
%                little room to fetch an item in one call.
%
%                [1,1] = size(room); int32 = class(room)
%
%   the call:
%
%      [values] = cspice_dskd02( handle, dladsc, item, start, room )
%
%   returns:
%
%      values    is a contiguous set of elements of the item designated by
%                `item'.
%
%                [1,n] = size(values); double = class(values)
%
%                The correspondence of `values' with the elements
%                of the data item is:
%
%                   values(1)      item(start)
%                     ...             ...
%                   values(n)      item(start+n-1)
%
%                If an error occurs on the call, `values' is undefined.
%
%                Note, `room' >= n.
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
%      function dskd02_t( dsk )
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
%   Most SPICE applications will not need to call this routine. The
%   routines cspice_dskv02, cspice_dskp02, and cspice_dskz02 provide
%   a higher-level interface for fetching DSK type 2 vertex and plate data.
%
%   DSK files are built using the DLA low-level format and
%   the DAS architecture; DLA files are a specialized type of DAS
%   file in which data are organized as a doubly linked list of
%   segments. Each segment's data belong to contiguous components of
%   character, double precision, and integer type.
%
%   Note that the DSK descriptor for the segment is not needed by this
%   routine; the DLA descriptor contains the base address and size
%   information for the integer, double precision, and character
%   components of the segment, and these suffice for the purpose of
%   fetching data.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine dskd02_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 04-APR-2017, NJB (JPL), EDW (JPL), ML (JPL)
%
%-Index_Entries
%
%   fetch double precision data from a type 2 dsk segment
%
%-&

function [values] = cspice_dskd02( handle, dladsc, item, start, room )

   switch nargin
      case 5

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         item   = zzmice_int( item   );
         start  = zzmice_int( start  );
         room   = zzmice_int( room, [1, int32(inf)/2] );

      otherwise

         error ( [ 'Usage: [values] = ' ...
                   'cspice_dskd02( handle, dladsc, item, start, room ) ' ] )
   end

   %
   % Call the MEX library.
   %
   try

      [values] = mice( 'dskd02_c', ...
                       handle, dladsc, item, start, room );
   catch
      rethrow(lasterror)
   end
