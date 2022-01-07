%-Abstract
%
%   CSPICE_DASUDI updates data in a specified range of integer addresses in a
%   DAS file.
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
%               of integers to update.
%
%               [1,1] = size(first); int32 = class(first)
%               [1,1] = size(last); int32 = class(last)
%
%               These addresses satisfy the inequality
%
%                  1  <=   first   <=   last   <=   lasti
%
%               where `lasti' is the last integer logical address in
%               use in the DAS file designated by `handle'.
%
%      data     an array of integers.
%
%               [n,1] = size(data); int32 = class(data)
%
%               The array elements data(1) through data(n) will be written
%               to the indicated DAS file, where `n' is last - first + 1.
%
%   the call:
%
%      cspice_dasudi( handle, first, last, data )
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
%   1) Write to addresses 1 through 200 in a DAS file in random-access
%      fashion by updating the file. Recall that data must be present
%      in the file before it can be updated.
%
%
%      Example code begins here.
%
%
%      function dasudi_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasudi_ex1.das';
%         TYPE  =   'TEST';
%
%         %
%         % Open a new DAS file. Use the file name as the internal
%         % file name, and reserve no records for comments.
%         %
%         [handle] = cspice_dasonw( FNAME, TYPE, FNAME, 0 );
%
%         %
%         % Append 200 integers to the file; after the data are
%         % present, we're free to update it in any order we
%         % please. (zero out an integer array.)
%         %
%         data = zeros( 200, 1, 'int32' );
%         cspice_dasadi( handle, data );
%
%         %
%         % Now the integer logical addresses 1:200 can be
%         % written to in random-access fashion.  We'll fill them
%         % in reverse order.
%         %
%         for i=200:-1:1
%
%            cspice_dasudi( handle, i, i, i );
%
%         end
%
%         %
%         % Close the file.
%         %
%         cspice_dascls( handle );
%
%         %
%         % Now make sure that we updated the file properly.
%         % Open the file for reading and dump the contents
%         % of the integer logical addresses 1:200.
%         %
%         [handle] = cspice_dasopr( FNAME );
%         [data]   = cspice_dasrdi( handle, 1, 200 );
%
%         fprintf( '\n' )
%         fprintf( 'Data from "%s":\n', FNAME )
%         fprintf( '\n' )
%         for i=0:19
%
%            for j=1:10
%
%               fprintf( '%5d', data(i*10+j) )
%
%            end
%            fprintf( '\n' )
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
%      Data from "dasudi_ex1.das":
%
%          1    2    3    4    5    6    7    8    9   10
%         11   12   13   14   15   16   17   18   19   20
%         21   22   23   24   25   26   27   28   29   30
%         31   32   33   34   35   36   37   38   39   40
%         41   42   43   44   45   46   47   48   49   50
%         51   52   53   54   55   56   57   58   59   60
%         61   62   63   64   65   66   67   68   69   70
%         71   72   73   74   75   76   77   78   79   80
%         81   82   83   84   85   86   87   88   89   90
%         91   92   93   94   95   96   97   98   99  100
%        101  102  103  104  105  106  107  108  109  110
%        111  112  113  114  115  116  117  118  119  120
%        121  122  123  124  125  126  127  128  129  130
%        131  132  133  134  135  136  137  138  139  140
%        141  142  143  144  145  146  147  148  149  150
%        151  152  153  154  155  156  157  158  159  160
%        161  162  163  164  165  166  167  168  169  170
%        171  172  173  174  175  176  177  178  179  180
%        181  182  183  184  185  186  187  188  189  190
%        191  192  193  194  195  196  197  198  199  200
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   This routine replaces the integer data in the specified range of
%   logical addresses within a DAS file with the contents of the
%   input array `data'.
%
%   The actual physical write operations that update the indicated
%   DAS file with the contents of the input array `data' might not take
%   place before this routine returns, since the DAS system buffers
%   data that is written as well as data that is read. In any case,
%   the data will be flushed to the file at the time the file is
%   closed, if not earlier. A physical write of all buffered
%   records can be forced by calling the Mice routine cspice_daswbr
%   (DAS, write buffered records).
%
%   In order to append integer data to a DAS file, filling in a range
%   of integer logical addresses that starts immediately after the
%   last integer logical address currently in use, the Mice
%   routine cspice_dasadi (DAS add data, integer) should be used.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is
%       signaled by a routine in the call tree of this routine.
%
%   2)  Only logical addresses that already contain data may be
%       updated: if either `first' or `last' are outside the range
%
%          [ 1,  lasti ]
%
%       where `lasti' is the last integer logical address that currently
%       contains data in the indicated DAS file, the error
%       SPICE(INVALIDADDRESS) is signaled by a routine in the call
%       tree of this routine. The DAS file will not be modified.
%
%   3)  If first > last but both addresses are valid, this routine
%       will not modify the indicated DAS file. No error will be
%       signaled.
%
%   4)  If an I/O error occurs during the data update attempted
%       by this routine, the error is signaled by a routine in the
%       call tree of this routine.
%
%   5)  If any of the input arguments, `handle', `first', `last' or
%       `data', is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   6)  If any of the input arguments, `handle', `first', `last' or
%       `data', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%   7)  If the data provided in `data' is insufficient to update first-last+1
%       integer addresses of the DAS file, an error is signaled by the Icy
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
%   -Mice Version 1.0.0, 19-JUL-2021 (JDR)
%
%-Index_Entries
%
%   update integer data in a DAS file
%
%-&
function cspice_dasudi( handle, first, last, data )

   switch nargin
      case 4

         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         data   = zzmice_int(data);

      otherwise

         error ( 'Usage: cspice_dasudi( handle, first, last, data )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dasudi_c', handle, first, last, data);
   catch spiceerr
      rethrow(spiceerr)
   end
