%-Abstract
%
%   CSPICE_DASUDD updates data in a specified range of double precision
%   addresses in a DAS file.
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
%               of double precision numbers to update.
%
%               [1,1] = size(first); int32 = class(first)
%               [1,1] = size(last); int32 = class(last)
%
%               These addresses satisfy the inequality
%
%                  1  <=   first   <=   last   <=   lastd
%
%               where `lastd' is the last double precision logical
%               address in use in the DAS file designated by
%               `handle'.
%
%      data     an array of double precision numbers.
%
%               [n,1] = size(data); double = class(data)
%
%               The array elements data(1) through data(n) will be written
%               to the indicated DAS file, where `n' is last - first + 1.
%
%   the call:
%
%      cspice_dasudd( handle, first, last, data )
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
%      function dasudd_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasudd_ex1.das';
%         TYPE  =   'TEST';
%
%         %
%         % Open a new DAS file. Use the file name as the internal
%         % file name, and reserve no records for comments.
%         %
%         [handle] = cspice_dasonw( FNAME, TYPE, FNAME, 0 );
%
%         %
%         % Append 200 double precision numbers to the file;
%         % after the data are present, we're free to update it
%         % in any order we please.  (zeros out a double
%         % precision array.)
%         %
%         data = zeros( 200, 1 );
%         cspice_dasadd( handle, data );
%
%         %
%         % Now the double precision logical addresses 1:200
%         % can be written to in random-access fashion. We'll
%         % fill them in reverse order.
%         %
%         for i=200:-1:1
%
%            cspice_dasudd( handle, i, i, double(i) );
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
%         % of the double precision logical addresses 1:200.
%         %
%         [handle] = cspice_dasopr( FNAME );
%         [data]   = cspice_dasrdd( handle, 1, 200 );
%
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
%      Data from "dasudd_ex1.das":
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
%         97.0   98.0   99.0  100.0  101.0  102.0  103.0  104.0
%        105.0  106.0  107.0  108.0  109.0  110.0  111.0  112.0
%        113.0  114.0  115.0  116.0  117.0  118.0  119.0  120.0
%        121.0  122.0  123.0  124.0  125.0  126.0  127.0  128.0
%        129.0  130.0  131.0  132.0  133.0  134.0  135.0  136.0
%        137.0  138.0  139.0  140.0  141.0  142.0  143.0  144.0
%        145.0  146.0  147.0  148.0  149.0  150.0  151.0  152.0
%        153.0  154.0  155.0  156.0  157.0  158.0  159.0  160.0
%        161.0  162.0  163.0  164.0  165.0  166.0  167.0  168.0
%        169.0  170.0  171.0  172.0  173.0  174.0  175.0  176.0
%        177.0  178.0  179.0  180.0  181.0  182.0  183.0  184.0
%        185.0  186.0  187.0  188.0  189.0  190.0  191.0  192.0
%        193.0  194.0  195.0  196.0  197.0  198.0  199.0  200.0
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   This routine replaces the double precision data in the specified
%   range of logical addresses within a DAS file with the contents of
%   the input array `data'.
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
%   In order to append double precision data to a DAS file, filling
%   in a range of double precision logical addresses that starts
%   immediately after the last double precision logical address
%   currently in use, the Mice routine cspice_dasadd (DAS add data,
%   double precision) should be used.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is
%       signaled by a routine in the call tree of this routine.
%
%   2)  Only logical addresses that already contain data may be
%       updated: if either `first' or `last' are outside the range
%
%          [ 1,  lastd ]
%
%       where `lastd' is the last double precision logical address that
%       currently contains data in the indicated DAS file, the error
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
%       double precision addresses of the DAS file, an error is signaled by
%       the Mice interface.
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
%   update double precision data in a DAS file
%
%-&
function cspice_dasudd( handle, first, last, data )

   switch nargin
      case 4

         handle = zzmice_int(handle);
         first  = zzmice_int(first);
         last   = zzmice_int(last);
         data   = zzmice_dp(data);

      otherwise

         error ( 'Usage: cspice_dasudd( handle, first, last, data(N) )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dasudd_c', handle, first, last, data);
   catch spiceerr
      rethrow(spiceerr)
   end
