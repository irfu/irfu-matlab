%-Abstract
%
%   CSPICE_DASRDC reads character data from a range of DAS logical addresses.
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
%      handle   a file handle for an open DAS file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      first,
%      last     a range of DAS character logical addresses.
%
%               [1,1] = size(first); int32 = class(first)
%               [1,1] = size(last); int32 = class(last)
%
%               `first' and `last' must be greater than or equal to 1 and
%               less than or equal to the highest character logical address
%               in the DAS file designated by `handle'.
%
%      bpos,
%      epos     the begin and end character positions that define the
%               substrings in each of the elements of the output array
%               `data_i' into which character data is to be read.
%
%               [1,1] = size(bpos); int32 = class(bpos)
%               [1,1] = size(epos); int32 = class(epos)
%
%      data_i   a two-dimensional array of 8-bit unsigned integers,
%               representing characters.
%
%               [n,m] = size(data_i); uint8 = class(data_i)
%
%               `data_i' must be declared at least as
%
%                  data_i = zeros( r, epos, 'uint8' )
%
%               with the dimension `r' being at least
%
%                  r = int32( ( last - first + sublen ) / sublen )
%
%               and `sublen', the length of each of the substrings read
%               into the array elements from the DAS file, being
%
%                  sublen  =  epos - bpos + 1
%
%   the call:
%
%      [data] = cspice_dasrdc( handle, first, last, bpos, epos, data_i )
%
%   returns:
%
%      data     a two-dimensional array of 8-bit unsigned integers,
%               representing characters.
%
%               [n,m] = size(data); uint8 = class(data)
%
%               On output, the character words in the logical address range
%               `first' through `last' are copied into the characters
%
%                  data(1, bpos),
%                  data(1, bpos+1),
%                        .
%                        .
%                        .
%                  data(1, epos),
%                  data(2, bpos),
%                  data(2, bpos+1),
%                        .
%                        .
%                        .
%                  data(r, bpos)
%                  data(r, bpos+1)
%                        .
%                        .
%                        .
%
%               in that order. Note that the character positions of `data_i'
%               **other** than the ones shown in the diagram remain
%               unmodified.
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
%      function dasrdc_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasrdc_ex1.das';
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
%      Data from "dasrdc_ex1.das":
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
%   This routine provides random read access to the character data in
%   a DAS file. These data are logically structured as a
%   one-dimensional array of characters.
%
%   However, since Fortran programs usually use strings rather than
%   arrays of individual characters, the interface of this routine
%   provides for extraction of data from a DAS file into an array of
%   strings.
%
%   cspice_dasrdc allows the caller to control the amount of character data
%   read into each array element. This feature allows a program to
%   read character data into an array that has a different string
%   length from the one used to write the character data, without
%   losing the correspondence between input and output array elements.
%   For example, an array of strings of 32 characters can be written
%   to a DAS file and read back by cspice_dasrdc into a buffer of strings
%   having length 80 characters, mapping each 32-character string to
%   characters 1--32 of the output buffer.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is signaled
%       by a routine in the call tree of this routine.
%
%   2)  If `epos' or `bpos' are outside of the range
%
%          [  1,  size(data_i,2)  ]
%
%       or if epos < bpos, the error SPICE(BADSUBSTRINGBOUNDS) is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If `first' or `last' are out of range, an error is signaled by a
%       routine in the call tree of this routine.
%
%   4)  If `first' is greater than `last', `data' is left unchanged.
%
%   5)  If any of the input arguments, `handle', `first', `last',
%       `bpos', `epos' or `data_i', is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   6)  If any of the input arguments, `handle', `first', `last',
%       `bpos', `epos' or `data_i', is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%   7)  If the room available in `data_i' is insufficient to hold the
%       last-first+1 characters (of type uint8) read from the DAS file
%       in sub-strings of epos-bpos+1 size, an error is signaled by the
%       Mice interface.
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
%   read character data from a DAS file
%
%-&
function [data] = cspice_dasrdc( handle, first, last, bpos, epos, data_i )

   switch nargin
      case 6

         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         bpos   = zzmice_int(bpos);
         epos   = zzmice_int(epos);

      otherwise

         error ( [ 'Usage: [data(n,m)] = '                                  ...
                   'cspice_dasrdc( handle, first, last, '                   ...
                   'bpos, epos, data(n,m)' ] )

   end

   %
   % Call the MEX library. Note that `data_i' is unchecked by this interface,
   % as it should be of a non-common class for Mice (uint8). Its checks will
   % be done by the gateway code. `data_i' needs to be transposed in order to
   % have the "lines" of (character) data in contiguous memory. Conversely,
   % `data' needs to be transposed before returning to the caller, in order
   % to have the uint8 2-dimensional array in the expected size and shape.
   %
   try
      [data] = mice('dasrdc_c', handle, first, last, bpos, epos, data_i')';
   catch spiceerr
      rethrow(spiceerr)
   end
