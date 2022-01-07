%-Abstract
%
%   CSPICE_DASADD adds an array of double precision numbers to a DAS file.
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
%      data     an array of double precision numbers to be added to the
%               specified DAS file.
%
%               [n,1] = size(data); double = class(data)
%
%               Elements 1 through n are appended to the double precision
%               data in the file.
%
%   the call:
%
%      cspice_dasadd( handle, data )
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
%   1) Create a new DAS file and add 200 double precision numbers
%      to it. Close the file, then re-open it and read the data back
%      out.
%
%
%      Example code begins here.
%
%
%      function dasadd_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasadd_ex1.das';
%         TYPE  =   'TEST';
%
%         %
%         % Local variables.
%         %
%         data = zeros(100,1);
%
%         %
%         % Open a new DAS file. Use the file name as the internal
%         % file name, and reserve no records for comments.
%         %
%         [handle] = cspice_dasonw( FNAME, TYPE, FNAME, 0 );
%
%         %
%         % Fill the array `data' with the double precision
%         % numbers 1.0 through 100.0, and add this array
%         % to the file.
%         %
%         for i=1:100
%
%            data(i) = double(i);
%
%         end
%
%         cspice_dasadd( handle, data );
%
%         %
%         % Now append the array `data' to the file again.
%         %
%         cspice_dasadd( handle, data );
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
%         [data]   = cspice_dasrdd( handle, 1, 200 );
%
%         %
%         % Dump the data to the screen. We should see the
%         % sequence 1.0, 2.0, ..., 100.0, 1.0, 2.0, ... , 100.0.
%         % The numbers will be represented as double precision
%         % numbers in the output.
%         %
%         fprintf( '\n' )
%         fprintf( 'Data from "%s":\n', FNAME )
%         fprintf( '\n' )
%         for i=0:24
%
%            for j=1:8
%
%               fprintf( '%7.1f', data(i*8+j) )
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
%      Data from "dasadd_ex1.das":
%
%          1.0    2.0    3.0    4.0    5.0    6.0    7.0    8.0
%          9.0   10.0   11.0   12.0   13.0   14.0   15.0   16.0
%         17.0   18.0   19.0   20.0   21.0   22.0   23.0   24.0
%         25.0   26.0   27.0   28.0   29.0   30.0   31.0   32.0
%         33.0   34.0   35.0   36.0   37.0   38.0   39.0   40.0
%         41.0   42.0   43.0   44.0   45.0   46.0   47.0   48.0
%         49.0   50.0   51.0   52.0   53.0   54.0   55.0   56.0
%         57.0   58.0   59.0   60.0   61.0   62.0   63.0   64.0
%         65.0   66.0   67.0   68.0   69.0   70.0   71.0   72.0
%         73.0   74.0   75.0   76.0   77.0   78.0   79.0   80.0
%         81.0   82.0   83.0   84.0   85.0   86.0   87.0   88.0
%         89.0   90.0   91.0   92.0   93.0   94.0   95.0   96.0
%         97.0   98.0   99.0  100.0    1.0    2.0    3.0    4.0
%          5.0    6.0    7.0    8.0    9.0   10.0   11.0   12.0
%         13.0   14.0   15.0   16.0   17.0   18.0   19.0   20.0
%         21.0   22.0   23.0   24.0   25.0   26.0   27.0   28.0
%         29.0   30.0   31.0   32.0   33.0   34.0   35.0   36.0
%         37.0   38.0   39.0   40.0   41.0   42.0   43.0   44.0
%         45.0   46.0   47.0   48.0   49.0   50.0   51.0   52.0
%         53.0   54.0   55.0   56.0   57.0   58.0   59.0   60.0
%         61.0   62.0   63.0   64.0   65.0   66.0   67.0   68.0
%         69.0   70.0   71.0   72.0   73.0   74.0   75.0   76.0
%         77.0   78.0   79.0   80.0   81.0   82.0   83.0   84.0
%         85.0   86.0   87.0   88.0   89.0   90.0   91.0   92.0
%         93.0   94.0   95.0   96.0   97.0   98.0   99.0  100.0
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   This routine adds double precision data to a DAS file by
%   "appending" them after any double precision data already in the
%   file. The sense in which the data are "appended" is that the
%   data will occupy a range of logical addresses for double precision
%   data that immediately follow the last logical address of a double
%   precision number that is occupied at the time this routine is
%   called. The diagram below illustrates this addition:
%
%      +-------------------------+
%      |    (already in use)     |  D.p. logical address 1
%      +-------------------------+
%                  .
%                  .
%                  .
%      +-------------------------+
%      |    (already in use)     |  last d.p. logical address
%      +-------------------------+  in use before call to cspice_dasadd
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
%   The logical organization of the double precision numbers in the
%   DAS file is independent of the location in the file of any data
%   of integer or character type.
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
%   In order to update double precision logical addresses that
%   already contain data, the Mice routine cspice_dasudd
%   (DAS update data, double precision) should be used.
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
%   add double precision data to a DAS file
%
%-&
function cspice_dasadd( handle, data )

   switch nargin
      case 2

         handle = zzmice_int(handle);
         data   = zzmice_dp(data);

      otherwise

         error ( 'Usage: cspice_dasadd( handle, data(n) )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dasadd_c', handle, data);
   catch spiceerr
      rethrow(spiceerr)
   end
