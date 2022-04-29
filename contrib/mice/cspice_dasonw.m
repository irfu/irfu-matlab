%-Abstract
%
%   CSPICE_DASONW opens a new DAS file and sets the file type.
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
%      fname    the name of a new DAS file to be created (and consequently
%               opened for write access).
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%      ftype    a string indicating the type of data placed into a DAS file.
%
%               [1,c2] = size(ftype); char = class(ftype)
%
%                  or
%
%               [1,1] = size(ftype); cell = class(ftype)
%
%               The first nonblank character and the three, or fewer,
%               characters immediately following it are stored as the part of
%               the file's ID word following the forward slash. It is an
%               error if `ftype' is blank.
%
%               The file type may not contain any nonprinting characters.
%               `ftype' is case sensitive.
%
%               NAIF has reserved for its own use file types
%               consisting of the upper case letters (A-Z) and the
%               digits 0-9. NAIF recommends lower case or mixed case
%               file types be used by all others in order to avoid any
%               conflicts with NAIF file types.
%
%      ifname   a string containing the internal file name for the new file.
%
%               [1,c3] = size(ifname); char = class(ifname)
%
%                  or
%
%               [1,1] = size(ifname); cell = class(ifname)
%
%               The name may contain as many as 60 characters. This should
%               uniquely identify the file.
%
%      ncomr    the number of comment records to allocate.
%
%               [1,1] = size(ncomr); int32 = class(ncomr)
%
%               Allocating comment records at file creation time may reduce
%               the likelihood of having to expand the comment area later.
%
%   the call:
%
%      [handle] = cspice_dasonw( fname, ftype, ifname, ncomr )
%
%   returns:
%
%      handle   the file handle associated with the file.
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
%   1) Create a new DAS file and add 200 integers to it. Close the
%      file, then re-open it and read the data back out.
%
%
%      Example code begins here.
%
%
%      function dasonw_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'dasonw_ex1.das';
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
%      Data from "dasonw_ex1.das":
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
%   The DAS files created by this routine have initialized file
%   records.
%
%   This routine creates a new DAS file and sets the type of the
%   file to the mnemonic code passed to it.
%
%-Exceptions
%
%   1)  If the input filename is blank, the error SPICE(BLANKFILENAME)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If the specified file cannot be opened without exceeding the
%       maximum allowed number of open DAS files, the error
%       SPICE(DASFTFULL) is signaled by a routine in the call tree of
%       this routine. No file will be created.
%
%   3)  If the file cannot be opened properly, an error is signaled
%       by a routine in the call tree of this routine. No file will
%       be created.
%
%   4)  If the initial records in the file cannot be written, an
%       error is signaled by a routine in the call tree of this
%       routine. No file will be created.
%
%   5)  If the file type is blank, the error SPICE(BLANKFILETYPE) is
%       signaled by a routine in the call tree of this routine.
%
%   6)  If the file type contains nonprinting characters---decimal
%       0-31 and 127-255---, the error SPICE(ILLEGALCHARACTER) is
%       signaled by a routine in the call tree of this routine.
%
%   7)  If the number of comment records allocated `ncomr' is negative,
%       the error SPICE(INVALIDCOUNT) is signaled by a routine in the
%       call tree of this routine.
%
%   8)  If any of the input arguments, `fname', `ftype', `ifname' or
%       `ncomr', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   9)  If any of the input arguments, `fname', `ftype', `ifname' or
%       `ncomr', is not of the expected type, or it does not have the
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
%   -Mice Version 1.0.0, 26-AUG-2021 (JDR)
%
%-Index_Entries
%
%   open a new DAS file
%   open a new DAS file with write access
%
%-&
function [handle] = cspice_dasonw( fname, ftype, ifname, ncomr )

   switch nargin
      case 4

         fname  = zzmice_str(fname);
         ftype  = zzmice_str(ftype);
         ifname = zzmice_str(ifname);
         ncomr  = zzmice_int(ncomr);

      otherwise

         error ( [ 'Usage: [handle] = '                                     ...
                   'cspice_dasonw( `fname`, `ftype`, `ifname`, ncomr )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice('dasonw_c', fname, ftype, ifname, ncomr);
   catch spiceerr
      rethrow(spiceerr)
   end
