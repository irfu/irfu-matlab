%-Abstract
%
%   CSPICE_DAFGDA reads the double precision data bounded by two addresses
%   within a DAF.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      handle   file handle referring to a DAF.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      baddr,
%      eaddr    initial and final addresses of a contiguous set of double
%               precision numbers within a DAF. Presumably, these make up
%               all or part of a particular array.
%
%               Note that DAF addresses begin at 1 as in the
%               FORTRAN version of the SPICE Toolkit.
%
%               [1,1] = size(baddr); int32 = class(baddr)
%               [1,1] = size(eaddr); int32 = class(eaddr)
%
%   the call:
%
%      data = cspice_dafgda( handle, baddr, eaddr )
%
%   returns:
%
%      data   are the double precision data contained between
%             the specified addresses within the specified file.
%
%             'data' has length = end - begin + 1.
%
%             [1,length] = size(data); double = class(data)
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
%   1) Open a type 8 SPK for read access, retrieve the data for
%      the first segment and identify the beginning and end addresses,
%      the number of data elements within, the size of the data array,
%      and print the first two records.
%
%      Use the SPK kernel below as input type 8 SPK file for the example.
%
%         mer1_ls_040128_iau2000_v1.bsp
%
%      Each segment contains only two records which provide the start
%      and end position for the MER-1 rover landing site in the IAU_MARS
%      frame. Since the landing site does not change over time, it is
%      expected that both records are equal.
%
%
%      Example code begins here.
%
%
%      function dafgda_ex1()
%
%         %
%         % Open the type 8 SPK for read access then read the
%         % data from the first segment.
%         %
%         handle = cspice_dafopr( 'mer1_ls_040128_iau2000_v1.bsp');
%
%         %
%         % Begin a forward search; find the first segment; read the
%         % segment summary.
%         %
%         cspice_dafbfs( handle )
%         found    = cspice_daffna;
%         [dc, ic] = cspice_dafgs( 2, 6 );
%
%         %
%         % Retrieve the data begin and end addresses.
%         %
%         baddr = ic(5);
%         eaddr = ic(6);
%
%         fprintf( 'Beginning address       : %d\n', baddr )
%         fprintf( 'Ending address          : %d\n', eaddr )
%         fprintf( 'Number of data elements : %d\n', eaddr - baddr + 1 )
%
%         %
%         % Extract all data bounded by the begin and end addresses.
%         %
%         data = cspice_dafgda( handle, baddr, eaddr );
%
%         %
%         % Check `data'. It should contain 2 * 6 + 4 elements.
%         %
%         fprintf( 'Size of data array      :(%d,%d)\n', size(data) )
%
%         %
%         % Check the data. Each set of 6 element records should possess the
%         % property:
%         %
%         %   record(6) = record(6)
%         %        i            i-1
%         %
%         fprintf( 'The first and second states stored in the segment:\n' );
%         fprintf( ' %9.3f ', data(1:6) )
%         fprintf('\n')
%
%         fprintf( ' %9.3f ', data(7:12) )
%         fprintf('\n')
%
%         %
%         % SAFELY close the file
%         %
%         cspice_dafcls(handle)
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Beginning address       : 897
%      Ending address          : 912
%      Number of data elements : 16
%      Size of data array      :(1,16)
%      The first and second states stored in the segment:
%        3376.422   -326.649   -115.392      0.000      0.000      0.000
%        3376.422   -326.649   -115.392      0.000      0.000      0.000
%
%
%-Particulars
%
%   The principal reason that DAFs are so easy to use is that
%   the data in each DAF are considered to be one long contiguous
%   set of double precision numbers. You can grab data from anywhere
%   within a DAF without knowing (or caring) about the physical
%   records in which they are stored.
%
%-Exceptions
%
%   1)  If `baddr' is zero or negative, the error SPICE(DAFNEGADDR)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If baddr > eaddr, the error SPICE(DAFBEGGTEND) is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If `handle' is invalid, an error is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If the range of addresses covered between `baddr' and `eaddr'
%       includes records that do not contain strictly double
%       precision data, then the values returned in `data' are
%       undefined. See the -Restrictions section below for details.
%
%   5)  If any of the input arguments, `handle', `baddr' or `eaddr',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   6)  If any of the input arguments, `handle', `baddr' or `eaddr',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  There are several types of records in a DAF. This routine
%       is only to be used to read double precision data bounded
%       between two DAF addresses. The range of addresses input
%       may not cross data and summary record boundaries.
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
%   -Mice Version 1.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Edited -Examples section to comply with NAIF standard. Added
%       example's problem statement and a reference to required SPK file.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 17-JUL-2012 (EDW)
%
%-Index_Entries
%
%   read data from DAF address
%
%-&

function [data] = cspice_dafgda( handle, baddr, eaddr)

   switch nargin
      case 3

         handle = zzmice_int(handle);
         baddr  = zzmice_int(baddr);
         eaddr  = zzmice_int(eaddr);

      otherwise

         error ( 'Usage: data = cspice_dafgda( handle, baddr, eaddr)' )

   end

   %
   % Call the MEX library.
   %
   try
      [data] = mice( 'dafgda_c', handle, baddr, eaddr );
   catch spiceerr
      rethrow(spiceerr)
   end
