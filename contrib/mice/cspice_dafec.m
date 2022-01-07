%-Abstract
%
%   CSPICE_DAFEC reads comment text from the comment area of a DAF.
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
%      handle   file handle referring to a DAF file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      bufsiz   the maximum number of comment lines to copy to `buffer'.
%
%               [1,1] = size(bufsiz); int32 = class(bufsiz)
%
%      buffln   allowed length of each string element of the output `buffer'.
%
%               [1,1] = size(buffln); int32 = class(buffln)
%
%               This length must be large enough to hold the longest output
%               string. The SPICE system imposes no limit on the length of
%               comment lines, so `buffln' normally should be set to a
%               "generous" value that is unlikely to be exceeded.
%
%   the call:
%
%      [buffer, done] = cspice_dafec( handle, bufsiz, buffln )
%
%   returns:
%
%      buffer   array containing the comment lines read from the DAF
%               associated with `handle'.
%
%               [n,c1] = size(buffer); char = class(buffer)
%
%               On output, `buffer' contains `bufsiz' or less strings of
%               comment text, with one comment line per string
%               (bufsiz >= n).
%
%      done     logical indicating whether or not all of the comment
%               lines from the comment area of the DAF have been read.
%
%               [1,1] = size(done); logical = class(done)
%
%               This variable has value true after the last comment line has
%               been read. It will have a false value otherwise.
%
%               If no comments exist in the comment area, this variable
%               returns as true.
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
%   1) The following example will extract the entire comment area of a
%      binary DAF, displaying the comments on the terminal screen.
%
%
%      Use the SPK kernel below as input DAF file for the example.
%
%         earthstns_itrf93_201023.bsp
%
%
%      Example code begins here.
%
%
%      function dafec_ex1()
%         %
%         % Define an SPK from which to extract the comment area.
%         %
%         SPK        = 'earthstns_itrf93_201023.bsp';
%
%         %
%         % Define the size of the comment area to read from the SPK.
%         % 15 lines, each with length 80 characters.
%         %
%         BUFSIZE    = 15;
%         LINLEN     = 80;
%
%         %
%         % Open the 'SPK' for reading, return the corresponding
%         % file handle to 'handle'.
%         %
%         handle = cspice_dafopr( SPK );
%
%         done = false;
%
%         [buf, done] = cspice_dafec( handle, BUFSIZE, LINLEN );
%         output = cellstr(buf);
%
%         for i=1:numel(output)
%            fprintf( '%s\n', char(output(i)) );
%         end
%
%         if done
%            fprintf( '\nAll comments read from file.\n' );
%         else
%            fprintf( '\nNot all comments read from file.\n' );
%         end
%
%         %
%         % SAFELY close the file.
%         %
%         cspice_dafcls( handle )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
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
%         =====================================================================
%
%         This file provides geocentric states---locations and velocities---***
%         set of DSN stations cited in the list below under "Position Data."***
%         position vectors point from the earth's barycenter to the stations***
%
%      Not all comments read from file.
%
%
%      Warning: incomplete output. 3 lines extended past the right
%      margin of the header and have been truncated. These lines are
%      marked by "***" at the end of each line.
%
%
%      Note that the program outputs BUFSIZ (15) lines from the SPK
%      comment area. Additional calls to cspice_dafec will read more
%      comment lines from the SPK in slices of BUFSIZ.
%
%      Reading all comment lines from the SPK requires a large value for
%      BUFSIZ. In this case, a BUFSIZ value of 139 will read all comment
%      lines from the SPK in a single cspice_dafec.
%
%-Particulars
%
%   A binary DAF contains an area which is reserved for storing
%   annotations or descriptive textual information describing the data
%   contained in a file. This area is referred to as the ``comment
%   area'' of the file. The comment area of a DAF is a line
%   oriented medium for storing textual information. The comment
%   area preserves any leading or embedded white space in the line(s)
%   of text which are stored, so that the appearance of the of
%   information will be unchanged when it is retrieved (extracted) at
%   some other time. Trailing blanks, however, are NOT preserved,
%   due to the way that character strings are represented in
%   standard Fortran 77.
%
%   This routine will read the comments from the comment area of
%   a binary DAF, placing them into a line buffer. If the line
%   buffer is not large enough to hold the entire comment area,
%   the portion read will be returned to the caller, and the DONE
%   flag will be set to false. This allows the comment area to be
%   read in ``chunks,'' a buffer at a time. After all of the comment
%   lines have been read, the `done' flag will be set to true.
%
%   This routine can be used to ``simultaneously'' extract comments
%   from the comment areas of multiple binary DAFs.
%
%-Exceptions
%
%   1)  If the size of the output line buffer is is not positive, the
%       error SPICE(INVALIDARGUMENT) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If a comment line in a DAF is longer than the length of a
%       character string array element of `buffer', the error
%       SPICE(COMMENTTOOLONG) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the end of the comments cannot be found, i.e., the end of
%       comments marker is missing on the last comment record, the
%       error SPICE(BADCOMMENTAREA) is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If the number of comment characters scanned exceeds the number
%       of comment characters computed, the error
%       SPICE(BADCOMMENTAREA) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the binary DAF attached to `handle' is not open for reading,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   6)  If any of the input arguments, `handle', `bufsiz' or `buffln',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   7)  If any of the input arguments, `handle', `bufsiz' or `buffln',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See argument `handle' in -I/O.
%
%-Restrictions
%
%   1)  The comment area may consist only of printing ASCII
%       characters, decimal values 32 - 126.
%
%   2)  There is NO maximum length imposed on the significant portion
%       of a text line that may be placed into the comment area of a
%       DAF. The maximum length of a line stored in the comment area
%       should be kept reasonable, so that they may be easily
%       extracted. A good value for this would be 1000 characters, as
%       this can easily accommodate ``screen width'' lines as well as
%       long lines which may contain some other form of information.
%
%   3)  This routine is only used to read records on environments
%       whose characters are a single byte in size. Updates
%       to this routine and routines in its call tree may be
%       required to properly handle other cases.
%
%   4)  This routine is intended to be used on DAF files whose comment
%       area does not change while this routine is called to extract
%       comments, between the start and end of the extraction process.
%       If the comment area does change (gets updated, reduced,
%       extended, or deleted) between calls to this routine on the
%       same DAF file, the routine's outputs are undefined and
%       subsequent calls to it are likely to trigger an exception.
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
%       Changed input argument name "lenout" to "buffln" for consistency
%       with other routines.
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and modified input kernel.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow. Fixed a bug in the
%       conversion of integer flag to MATLAB logical for return.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 20-JUL-2012 (EDW)
%
%-Index_Entries
%
%   extract comments from a DAF
%
%-&

function [buffer, done] = cspice_dafec( handle, bufsiz, buffln )

   switch nargin
      case 3

         handle  = zzmice_int(handle);
         bufsiz  = zzmice_int(bufsiz);
         buffln  = zzmice_int(buffln);

      otherwise

         error ( ['Usage: [buffer, done] = ' ...
                          'cspice_dafec( handle, bufsiz, buffln )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [buffer, done] = mice( 'dafec_c', handle, bufsiz, buffln );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      done = zzmice_logical(done);
   catch spiceerr
      rethrow(spiceerr)
   end
