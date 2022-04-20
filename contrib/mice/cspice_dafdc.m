%-Abstract
%
%   CSPICE_DAFDC deletes the entire comment area of a specified DAF file.
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
%      handle   the file handle referring to a DAF file opened with
%               write access.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               This handle refers to the DAF file from which to delete
%               the comment section.
%
%   the call:
%
%      cspice_dafdc( handle )
%
%   removes the comment area of the DAF file referred to by 'handle'.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Delete the entire comment area of a DAF file. Note that this
%      action should only be performed if fresh new comments are to
%      be placed within the DAF file.
%
%      Use the SPK kernel below as input DAF file for the program.
%
%         earthstns_itrf93_201023.bsp
%
%
%      Example code begins here.
%
%
%      function dafdc_ex1()
%
%         %
%         % Local parameters
%         %
%         KERNEL = 'earthstns_itrf93_201023.bsp';
%         BUFFSZ = 10;
%         LINLEN = 1000;
%
%         %
%         % Open a DAF for write. Return a `handle' referring to the
%         % file.
%         %
%         [handle] = cspice_dafopw( KERNEL );
%
%         %
%         % Print the first 10 lines of comments from the DAF file.
%         %
%         fprintf( 'Comment area of input DAF file (max. 10 lines): \n' )
%         fprintf( ['--------------------------------',                    ...
%                   '--------------------------------\n'] )
%
%         [buffer, done] = cspice_dafec( handle, BUFFSZ, LINLEN );
%
%         for i=1:size(buffer,1)
%            fprintf( '%s\n', buffer(i,:) )
%         end
%
%         fprintf( ['--------------------------------',                    ...
%                   '--------------------------------\n'] )
%         fprintf( ' \n' )
%         fprintf( 'Deleting entire comment area...\n' )
%
%         %
%         % Delete all the comments from the DAF file.
%         %
%         cspice_dafdc( handle );
%
%         %
%         % Close the DAF file and re-open it for read
%         % access to work around the cspice_dafec restriction
%         % on comments not to be modified while they are
%         % being extracted.
%         %
%         cspice_dafcls( handle );
%
%         [handle] = cspice_dafopr( KERNEL );
%
%         %
%         % Check if the comments have indeed been deleted.
%         %
%         [buffer, done] = cspice_dafec( handle, BUFFSZ, LINLEN );
%
%         if ( done && numel(buffer) == 0 )
%            fprintf( ' \n' )
%            fprintf( '   Successful operation.\n' )
%         else
%            fprintf( ' \n' )
%            fprintf( '   Operation failed.\n' )
%         end
%
%         %
%         % Safely close the DAF.
%         %
%         cspice_dafcls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Comment area of input DAF file (max. 10 lines):
%      ----------------------------------------------------------------
%
%         SPK for DSN Station Locations
%         =====================================================================
%
%         Original file name:                   earthstns_itrf93_201023.bsp
%         Creation date:                        2020 October 28 12:30
%         Created by:                           Nat Bachman  (NAIF/JPL)
%
%
%         Introduction
%      ----------------------------------------------------------------
%
%      Deleting entire comment area...
%
%         Successful operation.
%
%
%-Particulars
%
%   A binary DAF contains an area which is reserved for storing
%   annotations or descriptive textual information about the data
%   contained in a file. This area is referred to as the "comment
%   area" of the file. The comment area of a DAF is a line oriented
%   medium for storing textual information. The comment area preserves
%   any leading or embedded white space in the line(s) of text which are
%   stored, so that the appearance of the of information will be
%   unchanged when it is retrieved (extracted) at some other time.
%   Trailing blanks, however, are NOT preserved, due to the way that
%   character strings are represented in standard Fortran 77.
%
%   This routine will delete the entire comment area from the binary DAF
%   attached to `handle'. The size of the binary DAF will remain
%   unchanged. The space that was used by the comment records is
%   reclaimed: the data area of the DAF is shifted toward the beginning
%
%-Exceptions
%
%   1)  If the binary DAF attached to `handle' is not open with write
%       access, an error is signaled by a routine in the call tree of
%       this routine.
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
%   See argument `handle' in -I/O.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   DAF.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 25-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Added complete code example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 10-JUL-2012 (EDW)
%
%-Index_Entries
%
%   delete DAF comment area
%
%-&

function cspice_dafdc( handle )

   switch nargin
      case 1

         handle  = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dafdc(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafdc_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end




