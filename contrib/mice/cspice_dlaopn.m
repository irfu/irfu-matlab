%-Abstract
%
%   CSPICE_DLAOPN opens a new DLA file and sets the file type.
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
%      fname    the name of a new DLA file to be created.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%               The file will be left opened for write access.
%
%      ftype    a code for type of data placed into a DLA file.
%
%               [1,c2] = size(ftype); char = class(ftype)
%
%                  or
%
%               [1,1] = size(ftype); cell = class(ftype)
%
%               The non-blank part of `ftype' is used as the "file type"
%               portion of the ID word in the DLA file.
%
%               The first nonblank character and the three, or fewer,
%               characters immediately following it, giving four
%               characters, are used to represent the type of the data
%               placed in the DLA file. This is provided as a convenience
%               for higher level software. It is an error if this string
%               is blank. Also, the file type may not contain any
%               nonprinting characters. When written to the DLA file, the
%               value for the type IS case sensitive.
%
%               NAIF has reserved for its own use file types consisting
%               of the upper case letters (A-Z) and the digits 0-9. NAIF
%               recommends lower case or mixed case file types be used by
%               all others in order to avoid any conflicts with NAIF file
%               types.
%
%      ifname   the internal file name for the new file.
%
%               [1,c3] = size(ifname); char = class(ifname)
%
%                  or
%
%               [1,1] = size(ifname); cell = class(ifname)
%
%               The name may contain as many as 60 characters. This name
%               should uniquely identify the file.
%
%      ncomch   the number of comment characters to allocate.
%
%               [1,1] = size(ncomch); int32 = class(ncomch)
%
%               `ncomch' is used to establish the number of comment records
%               that will be allocated to the new DLA file. The number of
%               comment records allocated is the minimum required to
%               store the specified number of comment characters.
%
%               Allocating comment records at file creation time may
%               reduce the likelihood of having to expand the
%               comment area later.
%
%   the call:
%
%      [handle] = cspice_dlaopn( fname, ftype, ifname, ncomch )
%
%   returns:
%
%      handle   the file handle associated with the file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               This handle is used to identify the file in subsequent
%               calls to other DLA routines.
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
%   1) Create a DLA file containing one segment; the segment
%      contains character, double precision, and integer data.
%      After writing and closing the file, open the file for
%      read access; dump the data to standard output.
%
%
%      Example code begins here.
%
%
%      function dlaopn_ex1()
%
%         %
%         % MiceUser is a file that makes certain variables global.
%         % You must call MiceUser to have access to the parameters used
%         % in this example.
%         %
%         MiceUser;
%
%         %
%         % Local parameters
%         %
%         DLA    =   'dlaopn_ex1.dla';
%         LNSIZE =   61;
%         MAXC   =   5;
%         MAXD   =   50;
%         MAXI   =   100;
%
%         %
%         % Local variables
%         %
%         dvals  = zeros( MAXD, 1 );
%         ivals  = zeros( MAXI, 1, 'int32' );
%         cvals  = zeros( MAXC, LNSIZE, 'uint8' );
%         cvals2 = zeros( MAXC, LNSIZE, 'uint8' );
%
%         %
%         % Set the internal file name. Don't reserve characters in
%         % the DAS comment area.
%         %
%         ifname = 'Example DLA file for testing';
%         ncomch = 0;
%
%         %
%         % Open a new DLA file.
%         %
%         [handle] = cspice_dlaopn( DLA, 'DLA', ifname, ncomch );
%
%         %
%         % Begin a new segment.
%         %
%         cspice_dlabns( handle );
%
%         %
%         % Add character data to the segment.
%         %
%         for i=1:MAXC
%
%            for j=1:LNSIZE
%
%               k = mod( j+i-1, 10 );
%               cvals(i,j) = uint8('0') + k;
%
%            end
%
%         end
%
%         cspice_dasadc( handle, MAXC*LNSIZE, 1, LNSIZE, cvals );
%
%         %
%         % Add integer and double precision data to the segment.
%         %
%         for i=1:MAXI
%
%            ivals(i) = i;
%
%         end
%
%         cspice_dasadi( handle, ivals );
%
%         for i=1:MAXD
%
%            dvals(i) = double(i);
%
%         end
%
%         cspice_dasadd( handle, dvals );
%
%         %
%         % End the segment.
%         %
%         cspice_dlaens( handle );
%
%         %
%         % Close the file.  The routine cspice_dascls flushes the DAS
%         % buffers and segregates the file before closing it.
%         %
%         cspice_dascls( handle );
%
%         %
%         % Now read the file and check the data.
%         %
%         [handle] = cspice_dasopr( DLA );
%
%         %
%         % Obtain the segment descriptor for the sole segment
%         % in the file. We need not check the found flag
%         % in this case because we know there is one segment
%         % in the file.
%         %
%         [descr, found] = cspice_dlabfs( handle );
%
%         %
%         % Fetch character data from the segment.  Obtain the
%         % base address of the character data and the
%         % character count from the descriptor.
%         %
%         base     = descr(SPICE_DLA_CBSIDX);
%         n        = descr(SPICE_DLA_CSZIDX);
%
%         [cvals2] = cspice_dasrdc( handle, base+1, base+n,                ...
%                                   1,      LNSIZE, cvals2  );
%
%         %
%         % Display the character data.
%         %
%         fprintf( '\n' )
%         fprintf( 'Character array:\n' )
%         cvals2 = cellstr(char(cvals2));
%         for i=1:n/LNSIZE
%
%            fprintf( '%s\n', char(cvals2(i)) )
%
%         end
%
%         %
%         % Fetch and display the integer and double precision data.
%         %
%         base     = descr(SPICE_DLA_IBSIDX);
%         n        = descr(SPICE_DLA_ISZIDX);
%
%         [ivals2] = cspice_dasrdi( handle, base+1, base+n );
%
%         fprintf( '\n' )
%         fprintf( 'Integer array:\n' )
%         for i=0:n/10-1
%
%            for j=1:10
%
%               fprintf( '%6d', ivals2(i*10+j) )
%
%            end
%            fprintf( '\n' )
%
%         end
%
%         base     = descr(SPICE_DLA_DBSIDX);
%         n        = descr(SPICE_DLA_DSZIDX);
%
%         [dvals2] = cspice_dasrdd( handle, base+1, base+n );
%
%         fprintf( '\n' )
%         fprintf( 'Double precision array:\n' )
%         for i=0:n/10-1
%
%            for j=1:10
%
%               fprintf( '%6.1f', dvals2(i*10+j) )
%
%            end
%            fprintf( '\n' )
%
%         end
%
%         %
%         % Close the file.  This step is unnecessary in this
%         % program, but is a good practice in general
%         % because closing the file frees resources.
%         %
%         cspice_dascls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Character array:
%      1234567890123456789012345678901234567890123456789012345678901
%      2345678901234567890123456789012345678901234567890123456789012
%      3456789012345678901234567890123456789012345678901234567890123
%      4567890123456789012345678901234567890123456789012345678901234
%      5678901234567890123456789012345678901234567890123456789012345
%
%      Integer array:
%           1     2     3     4     5     6     7     8     9    10
%          11    12    13    14    15    16    17    18    19    20
%          21    22    23    24    25    26    27    28    29    30
%          31    32    33    34    35    36    37    38    39    40
%          41    42    43    44    45    46    47    48    49    50
%          51    52    53    54    55    56    57    58    59    60
%          61    62    63    64    65    66    67    68    69    70
%          71    72    73    74    75    76    77    78    79    80
%          81    82    83    84    85    86    87    88    89    90
%          91    92    93    94    95    96    97    98    99   100
%
%      Double precision array:
%         1.0   2.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0  10.0
%        11.0  12.0  13.0  14.0  15.0  16.0  17.0  18.0  19.0  20.0
%        21.0  22.0  23.0  24.0  25.0  26.0  27.0  28.0  29.0  30.0
%        31.0  32.0  33.0  34.0  35.0  36.0  37.0  38.0  39.0  40.0
%        41.0  42.0  43.0  44.0  45.0  46.0  47.0  48.0  49.0  50.0
%
%
%      Note that after run completion, a new DLA file exists in the
%      output directory.
%
%-Particulars
%
%   DLA files are built using the DAS low-level format; DLA files are
%   a specialized type of DAS file in which data are organized as a
%   doubly linked list of segments. Each segment's data belong to
%   contiguous components of character, double precision, and integer
%   type.
%
%   This routine creates a new DLA file and sets the type of the
%   file to the mnemonic code passed to it.
%
%   DLA files created by this routine have initialized file records.
%   The ID word in a DLA file record has the form
%
%      DAS/xxxx
%
%   where the characters following the slash are supplied by the
%   caller of this routine.
%
%-Exceptions
%
%   1)  If the input filename is blank, an error is signaled by a
%       routine in the call tree of this routine. No file will be
%       created.
%
%   2)  If the specified file cannot be opened without exceeding the
%       maximum allowed number of open DAS files, an error is signaled
%       by a routine in the call tree of this routine. No file will be
%       created.
%
%   3)  If the file cannot be opened properly, an error is signaled by
%       a routine in the call tree of this routine. No file will be
%       created.
%
%   4)  If the initial records in the file cannot be written, an error
%       is signaled by a routine in the call tree of this routine. No
%       file will be created.
%
%   5)  If no logical units are available, an error is signaled by a
%       routine in the call tree of this routine. No file will be
%       created.
%
%   6)  If the file type is blank, an error is signaled by a routine
%       in the call tree of this routine. No file will be created.
%
%   7)  If the file type contains nonprinting characters, decimal 0-31
%       and 127-255, an error is signaled by a routine in the call
%       tree of this routine. No file will be created.
%
%   8)  If the number of comment characters allocated to be allocated,
%       `ncomch', is negative, the error SPICE(BADRECORDCOUNT) is
%       signaled by a routine in the call tree of this routine. No
%       file will be created.
%
%   9)  If any of the input arguments, `fname', `ftype', `ifname' or
%       `ncomch', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   10) If any of the input arguments, `fname', `ftype', `ifname' or
%       `ncomch', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See argument `fname'.
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
%   -Mice Version 1.0.0, 01-NOV-2021 (JDR)
%
%-Index_Entries
%
%   open a new DLA file
%   open a new DLA file with write access
%
%-&
function [handle] = cspice_dlaopn( fname, ftype, ifname, ncomch )

   switch nargin
      case 4

         fname  = zzmice_str(fname);
         ftype  = zzmice_str(ftype);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int(ncomch);

      otherwise

         error ( [ 'Usage: [handle] = '                                     ...
                   'cspice_dlaopn( `fname`, `ftype`, `ifname`, ncomch )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice('dlaopn_c', fname, ftype, ifname, ncomch);
   catch spiceerr
      rethrow(spiceerr)
   end
