%-Abstract
%
%   CSPICE_DASLLA returns last DAS logical addresses of character, double
%   precision and integer type that are currently in use in a specified DAS
%   file.
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
%      handle   the file handle of a DAS file whose active logical address
%               ranges are desired.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      [lastc, lastd, lasti] = cspice_daslla( handle )
%
%   returns:
%
%      lastc,
%      lastd,
%      lasti    respectively, the last 1-based logical addresses of character,
%               double precision, and integer type in use in the specified DAS
%               file.
%
%               [1,1] = size(lastc); int32 = class(lastc)
%               [1,1] = size(lastd); int32 = class(lastd)
%               [1,1] = size(lasti); int32 = class(lasti)
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
%   1) Create a DAS file containing 10 integers, 5 double precision
%      numbers, and 4 characters, then use cspice_daslla to find the logical
%      address ranges in use.
%
%
%      Example code begins here.
%
%
%      function daslla_ex1()
%
%         %
%         % Local parameters.
%         %
%         FNAME =   'daslla_ex1.das';
%
%         %
%         % Open a new DAS file. Use the file name as the internal
%         % file name, and reserve no records for comments.
%         %
%         type     = 'TEST';
%         ifname   = 'TEST.DAS/NAIF/NJB/11-NOV-1992-20:12:20';
%
%         [handle] = cspice_dasonw( FNAME, type, ifname, 0 );
%
%         for i=1:10
%
%            cspice_dasadi( handle, i );
%
%         end
%
%         for i=1:5
%
%            cspice_dasadd( handle, double(i) );
%
%         end
%
%         %
%         % Add character data to the file. DAS character data are
%         % treated as a character array, not as a string. The
%         % following call adds only the first 4 characters to the
%         % DAS file.
%         %
%         cspice_dasadc( handle, 4, 1, 4, uint8('SPUDWXY') );
%
%         %
%         % Now check the logical address ranges.
%         %
%         [lastc, lastd, lasti] = cspice_daslla( handle );
%
%         fprintf( 'Last character address in use: %d\n', lastc )
%         fprintf( 'Last d.p. address in use     : %d\n', lastd )
%         fprintf( 'Last integer address in use  : %d\n', lasti )
%
%         %
%         % Close the DAS file.
%         %
%         cspice_dascls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Last character address in use: 4
%      Last d.p. address in use     : 5
%      Last integer address in use  : 10
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   This routine is a utility that allows a calling program to
%   find the range of logical addresses currently in use in any
%   DAS file.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is signaled by
%       a routine in the call tree of this routine.
%
%   2)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `handle' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
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
%   -Mice Version 1.0.0, 30-JUN-2021 (JDR)
%
%-Index_Entries
%
%   return last logical addresses in DAS file
%   return logical address range of DAS file
%
%-&
function [lastc, lastd, lasti] = cspice_daslla( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( 'Usage: [lastc, lastd, lasti] = cspice_daslla( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      [lastc, lastd, lasti] = mice('daslla_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end
