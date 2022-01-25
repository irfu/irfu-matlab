%-Abstract
%
%   CSPICE_DASRDI reads integer data from a range of DAS logical addresses.
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
%      last     the lower and upper bounds of a range of DAS integer
%               logical addresses.
%
%               [1,1] = size(first); int32 = class(first)
%               [1,1] = size(last); int32 = class(last)
%
%               The range includes these bounds. `first' and `last' must be
%               greater than or equal to 1 and less than or equal to the
%               highest integer DAS address in the DAS file designated by
%               `handle'.
%
%   the call:
%
%      [data] = cspice_dasrdi( handle, first, last )
%
%   returns:
%
%      data     an array of integers.
%
%               [n,1] = size(data); int32 = class(data)
%
%               `data' has length n = last - first + 1.
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
%      function dasrdi_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasrdi_ex1.das';
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
%      Data from "dasrdi_ex1.das":
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
%   This routine provides random read access to the integer data in
%   a DAS file. This data are logically structured as a
%   one-dimensional array of integers.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is signaled
%       by a routine in the call tree of this routine.
%
%   2)  If `first' or `last' are out of range, an error is signaled
%       by a routine in the call tree of this routine.
%
%   3)  If `first' is greater than `last', `data' is empty.
%
%   4)  If `data' is declared with length less than first - last + 1,
%       the error cannot be diagnosed by this routine.
%
%   5)  If a file read error occurs, the error is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If any of the input arguments, `handle', `first' or `last', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `handle', `first' or `last', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
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
%   read integer data from a DAS file
%
%-&
function [data] = cspice_dasrdi( handle, first, last )

   switch nargin
      case 3

         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);

      otherwise

         error ( 'Usage: [data] = cspice_dasrdi( handle, first, last )' )

   end

   %
   % Call the MEX library.
   %
   try
      [data] = mice('dasrdi_c', handle, first, last);
   catch spiceerr
      rethrow(spiceerr)
   end
