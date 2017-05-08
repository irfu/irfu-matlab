%-Abstract
%
%   CSPICE_DLABFS begins a forward segment search in a DLA file.
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
%      handle    the DAS integer handle associated with the file to be
%                searched. This handle is used to identify the file in
%                subsequent calls to other DLA or DAS routines.
%
%                [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      [dladsc, found] = cspice_dlabfs( handle )
%
%   returns:
%
%      dladsc    the descriptor of the first DLA segment in the file
%                associated with `handle'. `dladsc' is valid only if the
%                output argument `found' is true.
%
%                [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                int32 = class(dladsc)
%
%      found     a logical flag indicating whether a segment was found.
%                `found' has the value true if a segment was found;
%                otherwise `found' is false.
%
%                [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Open a DLA file for read access, traverse the segment
%   list from front to back, and display segment address
%   and size attributes.
%
%      function dlab_t( dla )
%
%          %
%          % Constants
%          %
%          SPICE_DLA_BWDIDX = 1;
%          SPICE_DLA_FWDIDX = 2;
%          SPICE_DLA_IBSIDX = 3;
%          SPICE_DLA_ISZIDX = 4;
%          SPICE_DLA_DBSIDX = 5;
%          SPICE_DLA_DSZIDX = 6;
%          SPICE_DLA_CBSIDX = 7;
%          SPICE_DLA_CSZIDX = 8;
%
%          %
%          % Open the DSK file for read access.
%          % We use the DAS-level interface for
%          % this function.
%          %
%          handle = cspice_dasopr( dla );
%
%          %
%          % Begin a forward search through the
%          % kernel, treating the file as a DLA.
%          % In this example, it's a very short
%          % search.
%          %
%          segno = 1;
%
%          [dladsc, found] = cspice_dlabfs( handle );
%
%          while  found
%
%             %
%             % Display the contents of the current segment
%             % descriptor.
%             %
%             fprintf('\n\n')
%             fprintf('Segment number = %d\n', segno )
%             fprintf('\n')
%             fprintf('   Backward segment pointer         = %d\n', ...
%                                        dladsc(SPICE_DLA_BWDIDX) )
%             fprintf('   Forward segment pointer          = %d\n', ...
%                                        dladsc(SPICE_DLA_FWDIDX) )
%             fprintf('   Integer component base address   = %d\n', ...
%                                        dladsc(SPICE_DLA_IBSIDX) )
%             fprintf('   Integer component size           = %d\n', ...
%                                        dladsc(SPICE_DLA_ISZIDX) )
%             fprintf('   D.p. component base address      = %d\n', ...
%                                        dladsc(SPICE_DLA_DBSIDX) )
%             fprintf('   D.p. component size              = %d\n', ...
%                                        dladsc(SPICE_DLA_DSZIDX) )
%             fprintf('   Character component base address = %d\n', ...
%                                        dladsc(SPICE_DLA_CBSIDX) )
%             fprintf('   Character component size         = %d\n', ...
%                                        dladsc(SPICE_DLA_CSZIDX) )
%             fprintf('\n')
%             %
%             % Find the next segment.
%             %
%             current = dladsc;
%             segno = segno + 1;
%
%             [dladsc, found] = cspice_dlafns( handle, current );
%
%          end
%
%          %
%          % Close file.
%          %
%          cspice_dascls( handle )
%
%   MATLAB outputs:
%
%      >> dlab_t( 'phobos_3_3.bds' )
%
%
%      Segment number = 1
%
%         Backward segment pointer         = -1
%         Forward segment pointer          = -1
%         Integer component base address   = 11
%         Integer component size           = 3311271
%         D.p. component base address      = 0
%         D.p. component size              = 494554
%         Character component base address = 0
%         Character component size         = 0
%
%-Particulars
%
%   DLA files are built using the DAS low-level format; DLA files are
%   a specialized type of DAS file in which data are organized as a
%   doubly linked list of segments. Each segment's data belong to
%   contiguous components of character, double precision, and integer
%   type.
%
%   This routine supports forward traversal of a DLA file's segment
%   list.  Note that it is not necessary to call this routine to
%   conduct a forward traversal; all that is necessary is to have
%   access to the first descriptor in the file, which this routine
%   provides.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine dlabfs_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 05-MAY-2014, NJB, EDW (JPL)
%
%-Index_Entries
%
%   begin forward search in dla file
%
%-&

function [dladsc, found] = cspice_dlabfs( handle )

   switch nargin
      case 1

         handle = zzmice_int( handle );

      otherwise

         error ( 'Usage: [dladsc(8), found] = cspice_dlabfs( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      [dladsc, found] = mice('dlabfs_c', handle );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end


