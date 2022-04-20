%-Abstract
%
%   CSPICE_DASADC adds character data to a DAS file.
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
%      n        the total number of characters to add to the specified DAS
%               file.
%
%               [1,1] = size(n); int32 = class(n)
%
%      bpos,
%      epos     the begin and end character positions that define the
%               substrings in each of the elements of the input array.
%
%               [1,1] = size(bpos); int32 = class(bpos)
%               [1,1] = size(epos); int32 = class(epos)
%
%               This routine writes the first `n' characters from the
%               specified set of substrings to the specified DAS file.
%
%      data     a two-dimensional array of 8-bit unsigned integers,
%               representing characters, some portion of whose contents are to
%               be added to the specified DAS file.
%
%               [n,m] = size(data); uint8 = class(data)
%
%               Specifically, the first `n' characters of the substrings
%
%                  data(i, bpos:epos),    i = 1, ...
%
%               are appended to the character data in the file.
%
%               `data' must be declared at least as
%
%                  data = zeros( r, epos, 'uint8' )
%
%               with the dimension `r' being at least
%
%                  r = int32( ( n + sublen - 1 ) / sublen )
%
%               and `sublen', the length of each of the substrings in
%               the array to be added to the DAS file, being
%
%                  sublen  =  epos - bpos + 1
%
%               The order of characters in the input substrings is
%               considered to increase from left to right within each
%               element of `data', and to increase with the indices of the
%               elements of `data'.
%
%   the call:
%
%      cspice_dasadc( handle, n, bpos, epos, data )
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
%      function dasadc_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasadc_ex1.das';
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
%      Data from "dasadc_ex1.das":
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
%   This routine adds character data to a DAS file by "appending" them
%   after any character data already in the file. The sense in which
%   the data are "appended" is that the data will occupy a range of
%   logical addresses for character data that immediately follow the
%   last logical address of a character that is occupied at the time
%   this routine is called. The diagram below illustrates this
%   addition:
%
%      +-------------------------+
%      |    (already in use)     |  Character logical address 1
%      +-------------------------+
%                  .
%                  .
%                  .
%      +-------------------------+  last character logical address
%      |   (already in use)      |  in use before call to cspice_dasadc
%      +-------------------------+
%      |     data(1,bpos)        |  First added character
%      +-------------------------+
%      |     data(1,bpos+1)      |
%      +-------------------------+
%                   .
%                   .
%                   .
%      +-------------------------+
%      |     data(1,epos)        |
%      +-------------------------+
%      |     data(2,bpos)        |
%      +-------------------------+
%                   .
%                   .
%                   .
%      +-------------------------+
%      |     data(r,C)           |  n'th added character---here `r' is
%      +-------------------------+
%                                      int( (n+l-1)/l )
%
%                                   where l = epos - bpos + 1, and
%                                   C is
%
%                                      bpos + ( n - (r-1)*l ) - 1
%
%
%   The logical organization of the characters in the DAS file is
%   independent of the order of addition to the file or physical
%   location of any data of integer or double precision type.
%
%   The actual physical write operations that add the input array
%   `data' to the indicated DAS file may not take place before this
%   routine returns, since the DAS system buffers data that are
%   written as well as data that are read. In any case, the data
%   will be flushed to the file at the time the file is closed, if
%   not earlier. A physical write of all buffered records can be
%   forced by calling the Mice routine cspice_daswbr (DAS, write
%   buffered records).
%
%   In order to update character logical addresses that already
%   contain data, the Mice routine cspice_dasudc (DAS, update data,
%   character) should be used.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is signaled
%       by a routine in the call tree of this routine.
%
%   2)  If `epos' or `bpos' are outside of the range
%
%          [  1,  size(data,2) ]
%
%       or if epos < bpos, the error SPICE(BADSUBSTRINGBOUNDS) is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If the input count `n' is less than 1, no data will be
%       added to the specified DAS file.
%
%   4)  If an I/O error occurs during the data addition attempted
%       by this routine, the error is signaled by a routine in the
%       call tree of this routine.
%
%   5)  If `n' is greater than the number of characters in the
%       specified set of input substrings, the results of calling
%       this routine are unpredictable. This routine cannot
%       detect this error.
%
%   6)  If any of the input arguments, `n', `handle', `bpos', `epos'
%       or `data', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   7)  If any of the input arguments, `n', `handle', `bpos', `epos'
%       or `data', is not of the expected type, or it does not have
%       the expected dimensions and size, an error is signaled by the
%       Mice interface.
%
%   8)  If the data provided in `data' is insufficient to add `n' characters
%       (of type uint8) to the DAS file, an error is signaled by the Mice
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
%   add character data to a DAS file
%
%-&
function cspice_dasadc( handle, n, bpos, epos, data )

   switch nargin
      case 5

         handle = zzmice_int(handle);
         n      = zzmice_int(n);
         bpos   = zzmice_int(bpos);
         epos   = zzmice_int(epos);

      otherwise

         error ( 'Usage: cspice_dasadc( handle, n, bpos, epos, data(n,m) )' )

   end

   %
   % Call the MEX library. Note that `data' is unchecked by this interface,
   % as it should be of a non-common class for Mice (uint8). Its checks will
   % be done by the gateway code. `data' needs to be transposed in order to
   % have the "lines" of (character) data in contiguous memory.
   %
   try
      mice('dasadc_c', handle, n, bpos, epos, data');
   catch spiceerr
      rethrow(spiceerr)
   end
