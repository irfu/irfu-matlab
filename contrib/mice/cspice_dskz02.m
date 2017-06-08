%-Abstract
%
%   CSPICE_DSKZ02 returns plate model size parameters---plate count and
%   vertex count---for a type 2 DSK segment.
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
%      [nv, np] = cspice_dskz02( handle, dladsc )
%
%   returns:
%
%      nv        the number of vertices belonging to the specified plate
%                model.
%
%                [1,1] = size(nv); int32 = class(nv)
%
%      np        the number of plates belonging to the specified plate
%                model.
%
%                [1,1] = size(np); int32 = class(np)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Look up all the vertices associated with each plate
%   of the model contained in a specified type 2 segment. For each
%   plate, display the plate's vertices and normal vector.
%
%   For this example, we'll show the context of this look-up:
%   opening the DSK file for read access, traversing a trivial,
%   one-segment list to obtain the segment of interest.
%
%      function dskz02_t( dsk )
%
%         %
%         % Declare DSK Mice parameters for use in API calls.
%         %
%         DSKMiceUser
%
%         %
%         % Constants
%         %
%         PBUFSIZ = 10000;
%
%         %
%         % Initial values
%         %
%         verts = zeros(3,3);
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
%         % Get segment vertex and plate counts.
%         %
%         [nv, np] = cspice_dskz02( handle, dladsc );
%
%         fprintf( '\n' )
%         fprintf( 'Number of vertices:  %d\n', nv )
%         fprintf( 'Number of plates:    %d\n', np )
%
%         %
%         %  Display the vertices of each plate.
%         %
%         remain = np;
%         start  = 1;
%
%         while (remain > 0 )
%
%            %
%            % `nread' is the number of plates we"ll read on this
%            % loop pass. Set `nread' to the minimum of PBUFSIZ
%            % and `remain'.
%            %
%            nread = min(PBUFSIZ, remain);
%
%            plates = cspice_dskp02( handle, dladsc, start, nread );
%
%            for  i = 1:(nread)
%
%               plix = start + i - 1;
%
%               %
%               %  Read the vertices of the current plate.
%               %
%               for  j = 1:3
%
%                  verts(j,:) = cspice_dskv02( handle, ...
%                                              dladsc, plates(j,i), 1 );
%
%               end
%
%
%               %
%               % Display the vertices of the ith plate:
%               %
%               fprintf( '\n' )
%               fprintf( 'Plate number: %d\n', i )
%
%               for  j = 1:3
%                  fprintf( '   Vertex %d: ( %16.8e %16.8e %16.8e )\n', ...
%                                                         j, verts(j,:) )
%               end
%
%               %
%               % Display the normal vector of the current plate:
%               %
%               normal = cspice_dskn02( handle, dladsc, plix );
%               fprintf( '   Normal:   ( %16.8e %16.8e %16.8e )\n', normal )
%
%            end
%
%            start  = start  + nread;
%            remain = remain - nread;
%
%         end
%
%         %
%         % Close file.
%         %
%         cspice_dascls( handle )
%
%   MATLAB outputs:
%
%      >> dskz02_t('solid.bds' )
%
%      Number of vertices:  12
%      Number of plates:    20
%
%      Plate number: 1
%         Vertex 1: (   0.00000000e+00   0.00000000e+00   1.17557000e+00 )
%         Vertex 2: (   1.05146000e+00   0.00000000e+00   5.25731000e-01 )
%         Vertex 3: (   3.24920000e-01   1.00000000e+00   5.25731000e-01 )
%         Normal:   (   4.91124160e-01   3.56821347e-01   7.94654382e-01 )
%
%      Plate number: 2
%         Vertex 1: (   0.00000000e+00   0.00000000e+00   1.17557000e+00 )
%         Vertex 2: (   3.24920000e-01   1.00000000e+00   5.25731000e-01 )
%         Vertex 3: (  -8.50651000e-01   6.18034000e-01   5.25731000e-01 )
%         Normal:   (  -1.87592328e-01   5.77350079e-01   7.94654645e-01 )
%
%           ...
%
%      Plate number: 19
%         Vertex 1: (  -3.24920000e-01  -1.00000000e+00  -5.25731000e-01 )
%         Vertex 2: (   0.00000000e+00   0.00000000e+00  -1.17557000e+00 )
%         Vertex 3: (   8.50651000e-01  -6.18034000e-01  -5.25731000e-01 )
%         Normal:   (   1.87592328e-01  -5.77350079e-01  -7.94654645e-01 )
%
%      Plate number: 20
%         Vertex 1: (   8.50651000e-01  -6.18034000e-01  -5.25731000e-01 )
%         Vertex 2: (   0.00000000e+00   0.00000000e+00  -1.17557000e+00 )
%         Vertex 3: (   8.50651000e-01   6.18034000e-01  -5.25731000e-01 )
%         Normal:   (   6.07061680e-01   0.00000000e+00  -7.94654715e-01 )
%
%-Particulars
%
%   This routine enables SPICE-based user applications to conveniently
%   fetch the plate and vertex counts of a type 2 DSK segment.
%
%   See the routine CSPICE_DSKB02 (DSK, fetch type 2 bookkeeping data)
%   for an interface that returns all type 2 DSK segment bookkeeping
%   data in a single call.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine dskz02_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 07-MAY-2014, NJB (JPL), EDW (JPL)
%
%-Index_Entries
%
%   fetch model size parameters from a type 2 dsk segment
%
%-&

function [nv, np] = cspice_dskz02( handle, dladsc )

   switch nargin
      case 2

         handle = zzmice_int( handle );

      otherwise

         error ( [ 'Usage: [nv, np] = ' ...
                   'cspice_dskz02( handle, dladsc ) ' ] )
   end

   %
   % Call the MEX library.
   %
   try

      [nv, np] = mice( 'dskz02_c', handle, dladsc  );
   catch
      rethrow(lasterror)
   end


