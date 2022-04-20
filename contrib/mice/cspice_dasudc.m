%-Abstract
%
%   CSPICE_DASUDC updates character data in a specified range of DAS logical
%   addresses with substrings of a character array.
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
%      handle   a file handle of a DAS file opened for writing.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      first,
%      last     the first and last of a range of DAS logical addresses
%               of characters.
%
%               [1,1] = size(first); int32 = class(first)
%               [1,1] = size(last); int32 = class(last)
%
%               These addresses satisfy the inequality
%
%                  1  <=   first   <=   last   <=   lastc
%
%               where `lastc' is the last character logical address
%               in use in the DAS file designated by `handle'.
%
%      bpos,
%      epos     the begin and end character positions that define the
%               substrings in each of the elements of the input array that
%               are to replace the data in the range of DAS character
%               addresses given by `first' and `last'.
%
%               [1,1] = size(bpos); int32 = class(bpos)
%               [1,1] = size(epos); int32 = class(epos)
%
%      data     a two-dimensional array of 8-bit unsigned integers,
%               representing characters.
%
%               [n,m] = size(data); uint8 = class(data)
%
%               The contents of the specified substrings of the elements of
%               the array `data' will be written to the indicated DAS file in
%               order: data(1,bpos) will be written to character logical
%               address `first'; data(1,bpos+1) will be written to the
%               character logical address first+1, and so on; in this ordering
%               scheme, character `bpos' of data(i+1,:) is the successor of
%               character `epos' of data(i,:).
%
%               `data' must be declared at least as
%
%                  data = zeros( r, epos, 'uint8' )
%
%               with the dimension `r' being at least
%
%                  r = int32( ( last - first + sublen ) / sublen )
%
%               and `sublen', the length of each of the substrings in
%               the array to be written to the DAS file, being
%
%                  sublen  =  epos - bpos + 1
%
%   the call:
%
%      cspice_dasudc( handle, first, last, bpos, epos, data )
%
%   returns:
%
%      None.
%
%      See -Particulars for a description of the effect of this routine.
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
%   1) The following example demonstrates the capabilities of the
%      DAS character data routines. The reader should notice that
%      in these interfaces, the character data are treated not as
%      strings (or arrays of strings) but as a stream of single
%      characters: DAS character data are not limited to
%      human-readable text. For example, one can store images or
%      DEM data as DAS character data.
%
%      The example shows how to add a variable amount of character
%      data to a new DAS file, how to update some of the character
%      logical addresses within that file, and how to read that
%      data out to a different array.
%
%
%      Example code begins here.
%
%
%      function dasudc_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasudc_ex1.das';
%         TYPE  =   'TEST';
%
%         %
%         % Local variables.
%         %
%
%         cdatou = { '..............................',                     ...
%                    '..............................',                     ...
%                    '..............................',                     ...
%                    '..............................',                     ...
%                    '..............................',                     ...
%                    '..............................',                     ...
%                    '..............................',                     ...
%                    '..............................',                     ...
%                    '         1         2         3',                     ...
%                    '123456789012345678901234567890' };
%
%         %
%         % Open a new DAS file. Use the file name as the internal
%         % file name, and reserve no records for comments.
%         %
%         [handle] = cspice_dasonw( FNAME, TYPE, FNAME, 0 );
%
%         %
%         % Set the input data. Note that these data will be
%         % considered as a binary data stream: DAS character data
%         % are not limited to human-readable text. For example,
%         % one can store images or DEM data as DAS character data.
%         %
%         cdatin = { '--F-345678901234567890',                             ...
%                    '--S-345678901234567890',                             ...
%                    '--T-IRDxxxxxxxxxxxxxxx' };
%
%         %
%         % Add the last 20 characters of the first two elements
%         % of `cdatin', and the 3rd character from the third one.
%         %
%         cdatin = uint8(char(cdatin));
%         cspice_dasadc( handle, 41, 3, 22, cdatin );
%
%         %
%         % Update the 10th, 20th and 30th character in the DAS
%         % file with a vertical bar.
%         %
%         for i=1:3
%
%            cspice_dasudc( handle, i*10, i*10, 1, 1, uint8('|') );
%
%         end
%
%         %
%         % Close the file.
%         %
%         cspice_dascls( handle );
%
%         %
%         % Now verify the addition of data by opening the
%         % file for read access and retrieving the data.
%         %
%         [handle] = cspice_dasopr( FNAME );
%
%         %
%         % Read the 41 characters that we stored on the DAS
%         % file. Update the data on the `cdatou' array, placing
%         % 6 characters on each element, starting from the
%         % 10th position.
%         %
%         cdatou = uint8(char(cdatou));
%         [cdatou] = cspice_dasrdc( handle, 1, 41, 10, 15, cdatou );
%
%         %
%         % Dump the data to the screen. Note that the last
%         % three lines should remain unmodified, and that
%         % only 5 characters will be written on the 7th line.
%         %
%         fprintf( '\n' )
%         fprintf( 'Data from "%s":\n', FNAME )
%         fprintf( '\n' )
%
%         cdatou = cellstr(char(cdatou));
%         for i=1:10
%
%            fprintf( '%s\n', char(cdatou(i)) )
%
%         end
%
%         %
%         % Close the file.
%         %
%         cspice_dascls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Data from "dasudc_ex1.das":
%
%      .........F-3456...............
%      .........789|12...............
%      .........345678...............
%      .........9|S-34...............
%      .........56789|...............
%      .........123456...............
%      .........7890T................
%      ..............................
%               1         2         3
%      123456789012345678901234567890
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   DAS is a low-level format meant to store and transmit data. As
%   such, character data in DAS files are not interpreted by Mice
%   DAS input or output routines. There are no limits on which
%   character values may be placed in the virtual character array of a
%   DAS file.
%
%   This routine replaces the character data in the specified range
%   of logical addresses within a DAS file with the contents of the
%   specified substrings of the input array `data'.
%
%   The actual physical write operations that update the indicated
%   DAS file with the contents of the input array `data' may not take
%   place before this routine returns, since the DAS system buffers
%   data that are written as well as data that are read. In any case,
%   the data will be flushed to the file at the time the file is
%   closed, if not earlier. A physical write of all buffered
%   records can be forced by calling the Mice routine cspice_daswbr
%   (DAS, write buffered records).
%
%   In order to append character data to a DAS file, filling in a
%   range of character logical addresses that starts immediately
%   after the last character logical address currently in use, the
%   Mice routine cspice_dasadc (DAS add data, character) should be
%   used.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is signaled by
%       a routine in the call tree of this routine.
%
%   2)  Only logical addresses that already contain data may be
%       updated: if either `first' or `last' are outside the range
%
%          [ 1,  lastc ]
%
%       where `lastc' is the last character logical address that
%       currently contains data in the indicated DAS file, the error
%       SPICE(INVALIDADDRESS) is signaled by a routine in the call
%       tree of this routine. The DAS file will not be modified.
%
%   3)  If `epos' or `bpos' are outside of the range
%
%          [  1,  size(data,2)  ]
%
%       the error SPICE(INVALIDINDEX) is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If `bpos' is greater than `epos', the error
%       SPICE(INDICESOUTOFORDER) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If first > last but both addresses are valid, this routine
%       will not modify the indicated DAS file. No error will be
%       signaled.
%
%   6)  If an I/O error occurs during the data update attempted
%       by this routine, the error is signaled by a routine in the
%       call tree of this routine. `first' and `last' will not be
%       modified.
%
%   7)  If any of the input arguments, `handle', `first', `last',
%       `bpos', `epos' or `data', is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   8)  If any of the input arguments, `handle', `first', `last',
%       `bpos', `epos' or `data', is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%   9)  If the data provided in `data' is insufficient to update first-last+1
%       character addresses of the DAS file, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See the description of the argument `handle' in -I/O.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   DAS.REQ
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
%   -Mice Version 1.0.0, 26-NOV-2021 (JDR)
%
%-Index_Entries
%
%   update a range of DAS logical addresses using substrings
%   write substrings to a range of DAS logical addresses
%
%-&
function cspice_dasudc( handle, first, last, bpos, epos, data )

   switch nargin
      case 6

         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         bpos   = zzmice_int(bpos);
         epos   = zzmice_int(epos);

      otherwise

         error ( [ 'Usage: cspice_dasudc( handle, first, last, bpos, '      ...
                   'epos, data(n,m) )' ] )

   end

   %
   % Call the MEX library. Note that `data' is unchecked by this interface,
   % as it should be of a non-common class for Mice (uint8). Its checks will
   % be done by the gateway code. `data' needs to be transposed in order to
   % have the "lines" of (character) data in contiguous memory.
   %
   try
      mice('dasudc_c', handle, first, last, bpos, epos, data');
   catch spiceerr
      rethrow(spiceerr)
   end
