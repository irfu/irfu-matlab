%-Abstract
%
%   CSPICE_DLABBS begins a backward segment search in a DLA file.
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
%      handle   the integer handle associated with the file to be searched.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               This handle is used to identify the file in subsequent
%               calls to other DLA or DAS routines.
%
%   the call:
%
%      [dladsc, found] = cspice_dlabbs( handle )
%
%   returns:
%
%      dladsc   the descriptor of the last DLA segment in the file associated
%               with `handle'.
%
%               [SPICE_DLA_DSCSIZ,1] = size(dladsc); int32 = class(dladsc)
%
%               The segment descriptor layout is:
%
%                  +---------------+
%                  | BACKWARD PTR  | Linked list backward pointer
%                  +---------------+
%                  | FORWARD PTR   | Linked list forward pointer
%                  +---------------+
%                  | BASE INT ADDR | Base DAS integer address
%                  +---------------+
%                  | INT COMP SIZE | Size of integer segment component
%                  +---------------+
%                  | BASE DP ADDR  | Base DAS d.p. address
%                  +---------------+
%                  | DP COMP SIZE  | Size of d.p. segment component
%                  +---------------+
%                  | BASE CHR ADDR | Base DAS character address
%                  +---------------+
%                  | CHR COMP SIZE | Size of character segment component
%                  +---------------+
%
%               `dladsc' is valid only if the output argument `found' is
%               true
%
%      found    a logical flag indicating whether a segment was found.
%
%               [1,1] = size(found); logical = class(found)
%
%               `found' has the value true if the file contains at least
%               one segment; otherwise `found' is false.
%
%-Parameters
%
%   SPICE_DLA_DSCSIZ
%
%               is the size of a SPICELIB DLA descriptor, defined in
%               MiceDLA.m.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Open a DLA file for read access, traverse the segment
%      list from back to front, and display segment address
%      and size attributes.
%
%      Example code begins here.
%
%
%      function dlabbs_ex1()
%
%         %
%         % MiceUser is a file that makes certain variables global.
%         % You must call MiceUser to have access to the parameters used
%         % in this example.
%         %
%         MiceUser;
%
%         %
%         % Prompt for the name of the file to search.
%         %
%         fname = input( 'Name of DLA file > ', 's' );
%
%         %
%         % Open the DLA file for read access.  Since DLA
%         % files use the DAS architecture, we can use DAS
%         % routines to open and close the file.
%         %
%         [handle] = cspice_dasopr( fname );
%
%         %
%         % Count the segments in the file; this allows us
%         % to label the segments in our display.
%         %
%         nsegs = 0;
%         [dladsc, found] = cspice_dlabbs( handle );
%
%         while found
%
%            nsegs  = nsegs + 1;
%            currnt = dladsc;
%            [dladsc, found] = cspice_dlafps( handle, currnt );
%
%         end
%
%         %
%         % Begin a backward search.  Let `dladsc' contain
%         % the descriptor of the last segment.
%         %
%         segno = nsegs + 1;
%
%         [dladsc, found] = cspice_dlabbs( handle );
%
%         while found
%
%            %
%            % Display the contents of the current segment
%            % descriptor.
%            %
%            segno = segno - 1;
%
%            fprintf( '\n' )
%            fprintf( '\n' )
%            fprintf( 'Segment number = %d\n', segno )
%            fprintf( '\n' )
%            fprintf( 'Backward segment pointer         = %d\n',           ...
%                                                 dladsc(SPICE_DLA_BWDIDX) )
%            fprintf( 'Forward segment pointer          = %d\n',           ...
%                                                 dladsc(SPICE_DLA_FWDIDX) )
%            fprintf( 'Character component base address = %d\n',           ...
%                                                 dladsc(SPICE_DLA_CBSIDX) )
%            fprintf( 'Character component size         = %d\n',           ...
%                                                 dladsc(SPICE_DLA_CSZIDX) )
%            fprintf( 'D.p. base address                = %d\n',           ...
%                                                 dladsc(SPICE_DLA_DBSIDX) )
%            fprintf( 'D.p. component size              = %d\n',           ...
%                                                 dladsc(SPICE_DLA_DSZIDX) )
%            fprintf( 'Integer base address             = %d\n',           ...
%                                                 dladsc(SPICE_DLA_IBSIDX) )
%            fprintf( 'Integer component size           = %d\n',           ...
%                                                 dladsc(SPICE_DLA_ISZIDX) )
%            fprintf( '\n' )
%
%            %
%            % Find the previous segment.
%            %
%            currnt = dladsc;
%            [dladsc, found] = cspice_dlafps( handle, currnt );
%
%         end
%
%         %
%         % Close the file using the DAS close routine.
%         %
%         cspice_dascls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the DSK file named phobos512.bds, the output
%      was:
%
%
%      Name of DLA file > phobos512.bds
%
%
%      Segment number = 1
%
%      Backward segment pointer         = -1
%      Forward segment pointer          = -1
%      Character component base address = 0
%      Character component size         = 0
%      D.p. base address                = 0
%      D.p. component size              = 4737076
%      Integer base address             = 11
%      Integer component size           = 29692614
%
%
%-Particulars
%
%   DLA files are built using the DAS low-level format; DLA files are
%   a specialized type of DAS file in which data are organized as a
%   doubly linked list of segments. Each segment's data belong to
%   contiguous components of character, double precision, and integer
%   type.
%
%   This routine supports backward traversal of a DLA file's segment
%   list. Note that it is not necessary to call this routine to
%   conduct a backward traversal; all that is necessary is to have
%   access to the last descriptor in the file, which this routine
%   provides.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is
%       signaled by a routine in the call tree of this routine.
%
%   2)  If an error occurs while reading the DLA file, the error
%       is signaled by a routine in the call tree of this routine.
%
%   3)  If the input descriptor is invalid, this routine will
%       fail in an unpredictable manner.
%
%   4)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   5)  If the input argument `handle' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   See description of input argument `handle'.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   DAS.REQ
%   DLA.REQ
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 19-JUL-2021 (JDR)
%
%-Index_Entries
%
%   begin backward search in DLA file
%
%-&
function [dladsc, found] = cspice_dlabbs( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( [ 'Usage: [dladsc(SPICE_DLA_DSCSIZ), found] = '            ...
                   'cspice_dlabbs( handle )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [dladsc, found] = mice('dlabbs_c', handle);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
