%-Abstract
%
%   CSPICE_DASADI adds an array of integers to a DAS file.
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
%      data     an array of integers to be added to the specified DAS file.
%
%               [n,1] = size(data); int32 = class(data)
%
%               Elements 1 through N are appended to the integer data in
%               the file.
%
%   the call:
%
%      cspice_dasadi( handle, data )
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
%   1) Create a new DAS file and add 200 integers to it. Close the
%      file, then re-open it and read the data back out.
%
%
%      Example code begins here.
%
%
%      function dasadi_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasadi_ex1.das';
%         TYPE  =   'TEST';
%
%         %
%         % Local variables.
%         %
%         data = zeros(100,1, 'int32');
%
%         %
%         % Open a new DAS file. Use the file name as the internal
%         % file name, and reserve no records for comments.
%         %
%         [handle] = cspice_dasonw( FNAME, TYPE, FNAME, 0 );
%
%         %
%         % Fill the array `data' with the integers 1 through
%         % 100, and add this array to the file.
%         %
%         for i=1:100
%
%            data(i) = i;
%
%         end
%
%         cspice_dasadi( handle, data );
%
%         %
%         % Now append the array `data' to the file again.
%         %
%         cspice_dasadi( handle, data );
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
%         [data]   = cspice_dasrdi( handle, 1, 200 );
%
%         %
%         % Dump the data to the screen.  We should see the
%         % sequence  1, 2, ..., 100, 1, 2, ... , 100.
%         %
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
%      Data from "dasadi_ex1.das":
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
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   This routine adds integer data to a DAS file by "appending" them
%   after any integer data already in the file. The sense in which
%   the data are "appended" is that the data will occupy a range of
%   logical addresses for integer data that immediately follow the
%   last logical address of an integer that is occupied at the time
%   this routine is called. The diagram below illustrates this
%   addition:
%
%      +-------------------------+
%      |    (already in use)     |  Integer logical address 1
%      +-------------------------+
%                  .
%                  .
%                  .
%      +-------------------------+
%      |    (already in use)     |  last integer logical address
%      +-------------------------+  in use before call to cspice_dasadi
%      |        data(1)          |
%      +-------------------------+
%                  .
%                  .
%                  .
%      +-------------------------+
%      |        data(n)          |
%      +-------------------------+
%
%
%   The logical organization of the integers in the DAS file is
%   independent of the location in the file of any data of double
%   precision or character type.
%
%   The actual physical write operations that add the input array
%   `data' to the indicated DAS file might not take place before this
%   routine returns, since the DAS system buffers data that are
%   written as well as data that are read. In any case, the data
%   will be flushed to the file at the time the file is closed, if
%   not earlier. A physical write of all buffered records can be
%   forced by calling the Mice routine cspice_daswbr (DAS, write
%   buffered records).
%
%   In order to update integer logical addresses that already contain
%   data, the Mice routine cspice_dasudi (DAS update data, integer)
%   should be used.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If an I/O error occurs during the data addition attempted by
%       this routine, the error is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If any of the input arguments, `handle' or `data', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `handle' or `data', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%   add integer data to a DAS file
%
%-&
function cspice_dasadi( handle, data )

   switch nargin
      case 2

         handle = zzmice_int(handle);
         data   = zzmice_int(data);

      otherwise

         error ( 'Usage: cspice_dasadi( handle, data(n) )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dasadi_c', handle, data);
   catch spiceerr
      rethrow(spiceerr)
   end
