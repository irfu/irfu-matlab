%-Abstract
%
%   CSPICE_DAFOPR opens a DAF for subsequent read requests.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      fname   the string name of a DAF to open for read (search) access.
%
%              [1,c1] = size(fname); char = class(fname)
%
%                 or
%
%              [1,1] = size(fname); cell = class(fname)
%
%   the call:
%
%      [handle] = cspice_dafopr( fname )
%
%   returns:
%
%      handle   the file handle used by other DAF routines
%               to refer to `fname'.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   Use cspice_dafcls to close files opened by this routine.
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
%   1) Use a simple routine to output the double precision and integer
%      values stored in an SPK's segments descriptors. This function
%      opens a DAF for read, performs a forwards search for the DAF
%      arrays, prints segments description for each array found, then
%      closes the DAF.
%
%      Use the SPK kernel below as input DAF file for the program.
%
%         de421.bsp
%
%
%      Example code begins here.
%
%
%      function dafopr_ex1()
%
%         %
%         % Local constants
%         %
%         kernel = 'de421.bsp';
%
%         %
%         % Open a DAF for read. Return a `handle' referring to the file.
%         %
%         handle = cspice_dafopr( kernel );
%
%         %
%         % Define the summary parameters appropriate
%         % for an SPK file.
%         %
%         ND = 2;
%         NI = 6;
%
%         %
%         % Begin a forward search on the file.
%         %
%         cspice_dafbfs( handle );
%
%         %
%         % Search until a DAF array is found.
%         %
%         found = cspice_daffna;
%
%         %
%         % Loop while the search finds subsequent DAF arrays.
%         %
%         while found
%
%            [dc, ic ] = cspice_dafgs( ND, NI );
%
%            fprintf( 'Doubles:  ' )
%            fprintf( '%f   ', dc )
%            fprintf( '\n' )
%
%            fprintf( 'Integers: ' )
%            fprintf( '%d   ', ic )
%            fprintf( '\n\n' )
%
%
%            %
%            % Check for another segment.
%            %
%            found = cspice_daffna;
%
%         end
%
%         %
%         % Safely close the DAF.
%         %
%         cspice_dafcls( handle )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 1   0   1   2   641   310404
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 2   0   1   2   310405   423048
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 3   0   1   2   423049   567372
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 4   0   1   2   567373   628976
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 5   0   1   2   628977   674740
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 6   0   1   2   674741   715224
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 7   0   1   2   715225   750428
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 8   0   1   2   750429   785632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 9   0   1   2   785633   820836
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 10   0   1   2   820837   944040
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 301   3   1   2   944041   1521324
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 399   3   1   2   1521325   2098608
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 199   1   1   2   2098609   2098620
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 299   2   1   2   2098621   2098632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 499   4   1   2   2098633   2098644
%
%
%      Note, the specific contents of `ic' and `dc' depend on the
%      type of DAF.
%
%      Note, the final entries in the integer array contain the segment
%      start/end indexes. The output indicates the search proceeded
%      from the start of the file (low value index) towards the end
%      (high value index).
%
%-Particulars
%
%   Most DAFs require only read access. If you do not need to
%   change the contents of a file, you should open it with cspice_dafopr.
%
%-Exceptions
%
%   1)  If the specified file has already been opened for read
%       access, the handle already associated with the file is
%       returned.
%
%   2)  If the specified file has already been opened for write
%       access, an error is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the specified file has already been opened by a non-DAF
%       routine, an error is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If the specified file cannot be opened without exceeding
%       the maximum number of files, the error SPICE(DAFFTFULL)
%       is signaled by a routine in the call tree of this routine.
%
%   5)  If the attempt to read the file's file record fails, the error
%       SPICE(FILEREADFAILED) is signaled by a routine in the call
%       tree of this routine.
%
%   6)  If the specified file is not a DAF file, an error is
%       signaled by a routine in the call tree of this routine.
%
%   7)  If no logical units are available, an error is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If the file does not exist, an error is signaled by a routine
%       in the call tree of this routine.
%
%   9)  If an i/o error occurs in the process of opening the file,
%       the error is signaled by a routine in the call tree of this
%       routine.
%
%   10) If the file name is blank or otherwise inappropriate,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   11) If the file was transferred improperly via FTP, an error is
%       signaled by a routine in the call tree of this routine.
%
%   12) If the file utilizes a binary file format that is not
%       currently supported on this platform, an error is signaled by
%       a routine in the call tree of this routine.
%
%   13) If the input argument `fname' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   14) If the input argument `fname' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   See argument `fname'.
%
%-Restrictions
%
%   1)  Files opened using this routine must be closed with cspice_dafcls.
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Edits to the -Examples section to comply with NAIF standard.
%       Modified code example to hardcode SPK file to be used as input.
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
%   open DAF for read
%
%-&

function [handle] = cspice_dafopr( fname )

   switch nargin
      case 1

         fname  = zzmice_str(fname);

      otherwise

         error ( 'Usage: [handle] = cspice_dafopr(`fname`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'dafopr_c', fname );
   catch spiceerr
      rethrow(spiceerr)
   end




