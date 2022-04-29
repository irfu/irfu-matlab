%-Abstract
%
%   CSPICE_DASOPW opens a DAS file for writing.
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
%      fname    the name of a DAS file to be opened with write access.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%   the call:
%
%      [handle] = cspice_dasopw( fname )
%
%   returns:
%
%      handle   the handle that is associated with the file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               This handle is used to identify the file in subsequent
%               calls to other DAS routines.
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
%   1) Create a new DAS file containing 200 integer addresses set
%      to zero. Re-open the file for write access again, and write
%      to its addresses 1 through 200 in random-access fashion by
%      updating the file.
%
%
%      Example code begins here.
%
%
%      function dasopw_ex1()
%
%         %
%         % Local Parameters
%         %
%         DASNAM =   'dasopw_ex1.das';
%         IDATLN =   200;
%
%         %
%         % Open a new DAS file. Reserve no comment records.
%         %
%         type     = 'TEST';
%         [handle] = cspice_dasonw( DASNAM, type, 'TEST/DASOPW_EX1', 0 );
%
%         %
%         % Append 200 integers to the file; after the data are
%         % present, we're free to update it in any order we
%         % please. (native zeros out an integer array.)
%         %
%         idata = zeros( IDATLN, 1, 'int32' );
%         cspice_dasadi( handle, idata );
%
%         %
%         % Close the file.
%         %
%         cspice_dascls( handle );
%
%         %
%         % Open the file again for writing.
%         %
%         [handle] = cspice_dasopw( DASNAM );
%
%         %
%         % Read the data from DAS file.
%         %
%         [idata] = cspice_dasrdi( handle, 1, IDATLN );
%
%         %
%         % Print the contents of the file before updating it.
%         %
%         fprintf( 'Contents of %s before update:\n', DASNAM )
%         fprintf( '\n' )
%         for i=0:19
%            for j=1:10
%               fprintf( '%5d', idata(i*10+j) )
%            end
%            fprintf( '\n' )
%         end
%
%         %
%         % Now the integer logical addresses 1:200 can be
%         % written to in random-access fashion. We'll fill them
%         % in reverse order.
%         %
%         for i=IDATLN:-1:1
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
%         [handle] = cspice_dasopr( DASNAM );
%
%         [idata]  = cspice_dasrdi( handle, 1, IDATLN );
%
%         fprintf( '\n' )
%         fprintf( 'Contents of %s after update:\n', DASNAM )
%         fprintf( '\n' )
%         for i=0:19
%            for j=1:10
%               fprintf( '%5d', idata(i*10+j) )
%            end
%            fprintf( '\n' )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Contents of dasopw_ex1.das before update:
%
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%          0    0    0    0    0    0    0    0    0    0
%
%      Contents of dasopw_ex1.das after update:
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
%   Most DAS files require only read access. If you do not need to
%   change the contents of a file, you should open it with cspice_dasopr.
%
%-Exceptions
%
%   1)  If the input filename is blank, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the specified file does not exist, an error is signaled by
%       a routine in the call tree of this routine.
%
%   3)  If the specified file has already been opened, either by the
%       DAS file routines or by other code, an error is signaled by a
%       routine in the call tree of this routine. Note that this
%       response is not paralleled by cspice_dasopr, which allows you to open
%       a DAS file for reading even if it is already open for reading.
%
%   4)  If the specified file cannot be opened without exceeding the
%       maximum allowed number of open DAS files, the error
%       SPICE(DASFTFULL) is signaled by a routine in the call tree of
%       this routine.
%
%   5)  If the specified file cannot be opened properly, an error
%       is signaled by a routine in the call tree of this routine.
%
%   6)  If the file record cannot be read, an error is signaled by a
%       routine in the call tree of this routine.
%
%   7)  If the specified file is not a DAS file, as indicated by the
%       file's ID word, an error is signaled by a routine in the call
%       tree of this routine.
%
%   8)  If no logical units are available, an error is signaled
%       by a routine in the call tree of this routine.
%
%   9)  If the input argument `fname' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   10) If the input argument `fname' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
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
%   open a DAS file for writing
%   open a DAS file for write access
%
%-&
function [handle] = cspice_dasopw( fname )

   switch nargin
      case 1

         fname = zzmice_str(fname);

      otherwise

         error ( 'Usage: [handle] = cspice_dasopw( `fname` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice('dasopw_c', fname);
   catch spiceerr
      rethrow(spiceerr)
   end
