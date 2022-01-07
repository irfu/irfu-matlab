%-Abstract
%
%   CSPICE_DAFAC adds comments from a buffer of character strings to the
%   comment area of a binary DAF file, appending them to any comments which
%   are already present in the file's comment area.
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
%      handle   file handle referring to a DAF.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      buffer   vector containing comments which to write into
%               the comment area of the binary DAF attached to `handle'.
%
%               [n,c1] = size(buffer); char = class(buffer)
%
%                  or
%
%               [1,n] = size(buffer); cell = class(buffer)
%
%               Each element of `buffer' should contain one comment line.
%
%   the call:
%
%      cspice_dafac( handle, buffer )
%
%   adds the contents of `buffer' to the DAF referred to by `handle'.
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
%   1) This example demonstrates how to append new comments to the
%      comment area of a DAF file.
%
%      Use the SPK kernel below as input DAF file for the program.
%
%         earthstns_itrf93_201023.bsp
%
%
%      Example code begins here.
%
%
%      function dafac_ex1()
%
%         %
%         % Local parameters
%         %
%         KERNEL = 'earthstns_itrf93_201023.bsp';
%         BUFFSZ = 25;
%         CMTSIZ = 7;
%         LINLEN = 1000;
%
%         %
%         % Set the new comments to be added to the DAF file.
%         %
%         newcmt = {                                                       ...
%                  '================== NEW COMMENTS ==================',   ...
%                  '',                                                     ...
%                  '   New comments can be appended to the end of the',    ...
%                  '   comment area of a DAF file, with a single',         ...
%                  '   operation.',                                        ...
%                  '',                                                     ...
%                  '================ END NEW COMMENTS ================' };
%
%
%         %
%         % Open a DAF for write. Return a `handle' referring to the
%         % file.
%         %
%         [handle] = cspice_dafopw( KERNEL );
%
%         %
%         % Print the end of comment area from the DAF file.
%         % (Maximum 15 lines.)
%         %
%         done = false;
%
%         while ( ~ done )
%
%            [buffer, done] = cspice_dafec( handle, 15, LINLEN );
%
%            if ( done )
%
%               fprintf( [ 'End of comment area of input DAF file (max.',  ...
%                          ' 15 lines):\n' ]                              )
%               fprintf( ['--------------------------------',              ...
%                         '--------------------------------\n'] )
%
%               for i=1:size(buffer,1)
%                  fprintf( '%s\n', buffer(i,:) )
%               end
%
%               fprintf( ['--------------------------------',              ...
%                         '--------------------------------\n'] )
%
%            end
%         end
%
%         %
%         % Append the new comments to the DAF file.
%         %
%         cspice_dafac( handle, newcmt );
%
%         %
%         % Safely close the DAF.
%         %
%         cspice_dafcls( handle );
%
%         %
%         % Check if the comments have indeed appended.
%         %
%         % Open a DAF for read.
%         %
%         [handle] = cspice_dafopr( KERNEL );
%         done     = false;
%
%         while ( ~ done )
%
%            [buffer, done] = cspice_dafec( handle, BUFFSZ, LINLEN );
%
%            if ( done )
%
%               fprintf( [ 'End of comment area of input DAF file (max.',  ...
%                          ' 25 lines):\n' ]                              )
%               fprintf( ['--------------------------------',              ...
%                         '--------------------------------\n'] )
%
%               for i=1:size(buffer,1)
%                  fprintf( '%s\n', buffer(i,:) )
%               end
%
%               fprintf( ['--------------------------------',              ...
%                         '--------------------------------\n'] )
%
%            end
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
%      End of comment area of input DAF file (max. 15 lines):
%      ----------------------------------------------------------------
%         DSS-65_DXYZ       =    (    -0.0100          0.0242          0.015***
%         DSS-65_TOPO_EPOCH =       @2020-OCT-23/00:00
%         DSS-65_UP         =       'Z'
%         DSS-65_NORTH      =       'X'
%
%      \begintext
%      ----------------------------------------------------------------
%      End of comment area of input DAF file (max. 25 lines):
%      ----------------------------------------------------------------
%         DSS-65_DXYZ       =    (    -0.0100          0.0242          0.015***
%         DSS-65_TOPO_EPOCH =       @2020-OCT-23/00:00
%         DSS-65_UP         =       'Z'
%         DSS-65_NORTH      =       'X'
%
%      \begintext
%      ================== NEW COMMENTS ==================
%
%         New comments can be appended to the end of the
%         comment area of a DAF file, with a single
%         operation.
%
%      ================ END NEW COMMENTS ================
%      ----------------------------------------------------------------
%
%
%      Warning: incomplete output. 2 lines extended past the right
%      margin of the header and have been truncated. These lines are
%      marked by "***" at the end of each line.
%
%
%-Particulars
%
%   A binary DAF contains a data area which is reserved for storing
%   annotations or descriptive textual information about the data
%   contained in a file. This area is referred to as the "comment
%   area" of the file. The comment area of a DAF is a line oriented
%   medium for storing textual information. The comment area preserves
%   leading or embedded white space in the line(s) of text which are
%   stored so that the appearance of the information will be unchanged
%   when it is retrieved (extracted) at some other time. Trailing
%   blanks, however, are NOT preserved, due to the way that character
%   strings are represented in standard Fortran 77.
%
%   This routine will take a buffer of text lines and add (append) them
%   to the comment area of a binary DAF. If there are no comments in the
%   comment area of the file, then space will be allocated and the text
%   lines in `buffer' will be placed into the comment area. The text lines
%   may contain only printable ASCII characters (decimal values 32 -
%   126).
%
%   There is NO maximum length imposed on the significant portion of a
%   text line that may be placed into the comment area of a DAF. The
%   maximum length of a line stored in the comment area should be
%   reasonable, however, so that they may be easily extracted. A good
%   maximum value for this would be 255 characters, as this can easily
%   accommodate "screen width" lines as well as long lines which may
%   contain some other form of information.
%
%-Exceptions
%
%   1)  If the number of comments to be added is not positive, the
%       error SPICE(INVALIDARGUMENT) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If a non printing ASCII character is encountered in the
%       comments, the error SPICE(ILLEGALCHARACTER) is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If the binary DAF file attached to `handle' is not open with
%       write access, an error is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If the end of the comments cannot be found, i.e., the end of
%       comments marker is missing on the last comment record, the
%       error SPICE(BADCOMMENTAREA) is signaled by a routine in the
%       call tree of this routine.
%
%   5)  If any of the input arguments, `handle' or `buffer', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `handle' or `buffer', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   See argument `handle' in -I/O.
%
%-Restrictions
%
%   1)  This routine uses constants that are specific to the ASCII
%       character sequence. The results of using this routine with
%       a different character sequence are unpredictable.
%
%   2)  This routine is only used to extract records on environments
%       whose characters are a single byte in size. Updates to this
%       routine and routines in its call tree may be required to
%       properly handle other cases.
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
%   -Mice Version 1.1.0, 25-NOV-2020 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added complete
%       code example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 20-JUL-2012 (EDW)
%
%-Index_Entries
%
%   add comments to a binary DAF
%   append comments to a DAF comment area
%
%-&

function cspice_dafac( handle, buffer )

   switch nargin
      case 2

         handle  = zzmice_int(handle);
         buffer  = zzmice_str(buffer);

      otherwise

         error ( 'Usage: cspice_dafac( handle, buffer )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafac_c', handle, buffer );
   catch spiceerr
      rethrow(spiceerr)
   end




