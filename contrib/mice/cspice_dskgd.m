%-Abstract
%
%   CSPICE_DSKGD returns the DSK descriptor from a DSK segment
%   identified by a DAS handle and DLA descriptor.
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
%      handle    the handle of a DSK file that is open for
%                read access.
%
%                [1,1] = size(handle); int32 = class(handle)
%
%      dladsc    the DLA segment descriptor corresponding to
%                a DSK segment.
%
%                [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                int32 = class(dladsc)
%
%   the call:
%
%      dskdsc = cspice_dskgd( handle, dladsc )
%
%   returns:
%
%      dskdsc    the DSK segment descriptor of the segment
%                designated by the input handle and DLA descriptor.
%
%                [SPICE_DSK_DSCSIZ,1]  = size(dskdsk)
%                               double = class(dskdsk)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Dump the DSK descriptor of the first segment of a DSK file.
%
%      function dskgd_t( dsk )
%
%         %
%         % Declare DSK Mice parameters for use in API calls.
%         %
%         MiceUser
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
%         % Fetch the DSK descriptor of this segment.
%         %
%         dskdsc = cspice_dskgd( handle, dladsc );
%
%         fprintf( '\n' )
%         fprintf( 'DSK descriptor contents:\n\n' )
%         fprintf( '   Surface ID:              %d\n',  ...
%                        dskdsc(SPICE_DSK_SRFIDX) )
%         fprintf( '   Center ID:               %d\n',  ...
%                        dskdsc(SPICE_DSK_CTRIDX) )
%         fprintf( '   Data class:               %d\n',  ...
%                        dskdsc(SPICE_DSK_CLSIDX) )
%         fprintf( '   Data type:                %d\n',  ...
%                        dskdsc(SPICE_DSK_TYPIDX) )
%         fprintf( '   Frame ID:                 %d\n',  ...
%                        dskdsc(SPICE_DSK_FRMIDX) )
%         fprintf( '   Coordinate system:        %d\n',  ...
%                        dskdsc(SPICE_DSK_SYSIDX) )
%         fprintf( '   Parameters:       %15.6f\n',      ...
%                        dskdsc(SPICE_DSK_PARIDX) )
%
%         for   i = 1:(SPICE_DSK_NSYPAR-1)
%
%            fprintf( '                     %15.6f\n', ...
%                           dskdsc(SPICE_DSK_PARIDX + i) )
%         end
%
%         fprintf( '   Coordinate 1 min: %15.6f\n',      ...
%                        dskdsc(SPICE_DSK_MN1IDX) )
%         fprintf( '   Coordinate 1 max: %15.6f\n',      ...
%                        dskdsc(SPICE_DSK_MX1IDX) )
%         fprintf( '   Coordinate 2 min: %15.6f\n',      ...
%                        dskdsc(SPICE_DSK_MN2IDX) )
%         fprintf( '   Coordinate 2 max: %15.6f\n',      ...
%                        dskdsc(SPICE_DSK_MX2IDX) )
%         fprintf( '   Coordinate 3 min: %15.6f\n',      ...
%                        dskdsc(SPICE_DSK_MN3IDX) )
%         fprintf( '   Coordinate 3 max: %15.6f\n',      ...
%                        dskdsc(SPICE_DSK_MX3IDX) )
%
%         %
%         % Close file.
%         %
%         cspice_dascls( handle )
%
%   MATLAB outputs:
%
%      >> dskgd_t( 'phobos512.bds' )
%
%      DSK descriptor contents:
%
%         Surface ID:              401
%         Center ID:               401
%         Data class:               1
%         Data type:                2
%         Frame ID:                 10021
%         Coordinate system:        1
%         Parameters:              0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%         Coordinate 1 min:       -3.141593
%         Coordinate 1 max:        3.141593
%         Coordinate 2 min:       -1.570796
%         Coordinate 2 max:        1.570796
%         Coordinate 3 min:        8.049632
%         Coordinate 3 max:       13.940940
%
%-Particulars
%
%   This is a convenience routine intended for use by low-level
%   routines that read DSK segments.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine dskgd_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 05-MAY-2014, NJB (JPL), EDW (JPL)
%
%-Index_Entries
%
%   return dsk segment_descriptor
%
%-&

function [dskdsc] = cspice_dskgd( handle, dladsc )

   switch nargin
      case 2

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );

      otherwise

         error ( 'Usage: [dskdsc(24)] = cspice_dskgd( handle, dladsc )' )

   end

   %
   % Call the MEX library.
   %
   try
      [dskdsc] = mice('dskgd_c', handle, dladsc );
   catch
      rethrow(lasterror)
   end


