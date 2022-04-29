%-Abstract
%
%   CSPICE_TSETYR sets the lower bound on the 100 year range.
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
%      year     the year associated with the lower bound on all year
%               expansions computed by the SPICELIB routine TEXPYR.
%
%               [1,1] = size(year); int32 = class(year)
%
%               For example if `year' is 1980, then the range of years that
%               can be abbreviated is from 1980 to 2079.
%
%   the call:
%
%      cspice_tsetyr( year )
%
%   sets the lower bound on the 100 year range.
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
%   1) Suppose that you need to manipulate time strings and that
%      you want to treat years components in the range from 0 to 99
%      as being abbreviations for years in the range from
%      1980 to 2079 (provided that the years are not modified by
%      an ERA substring). The example code below shows how you
%      could go about this.
%
%      Use the LSK kernel below to load the leap seconds and time
%      constants required for the conversions.
%
%         naif0012.tls
%
%
%      Example code begins here.
%
%
%      function tsetyr_ex1()
%
%         %
%         % Local parameters.
%         %
%         NTSTRS = 7;
%
%         %
%         % Assign an array of calendar dates.
%         %
%         date = { '00 JAN 21', '01 FEB 22', '48 MAR 23', '49 APR 24',     ...
%                  '79 JUL 14', '80 FEB 02', '99 DEC 31'               };
%
%         %
%         % Load the required LSK.
%         %
%         cspice_furnsh( 'naif0012.tls' );
%
%         %
%         % Set up the lower bound for the
%         % expansion of abbreviated years.
%         %
%         cspice_tsetyr( 1980 );
%
%         %
%         % Expand the years in input time strings.
%         %
%         fprintf( 'Time string    Expansion\n' )
%         fprintf( '-----------    -----------\n' )
%
%         for i=1:NTSTRS
%
%            [et]     = cspice_str2et( date(i) );
%            [timstr] = cspice_timout( et, 'YYYY MON DD' );
%
%            fprintf( '%s      %s\n', char(date(i)), timstr )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Time string    Expansion
%      -----------    -----------
%      00 JAN 21      2000 JAN 21
%      01 FEB 22      2001 FEB 22
%      48 MAR 23      2048 MAR 23
%      49 APR 24      2049 APR 24
%      79 JUL 14      2079 JUL 14
%      80 FEB 02      1980 FEB 02
%      99 DEC 31      1999 DEC 31
%
%
%-Particulars
%
%   This routine is used to set the range to which years
%   abbreviated to the last two digits will be expanded, allowing all
%   of the SPICE time subsystem routines to handle uniformly the
%   expansion those "abbreviated" years (i.e. the remainder after
%   dividing the actual year by 100.) The input supplied to this
%   routine represents the lower bound of the expansion interval. The
%   upper bound of the expansion interval is year + 99.
%
%   The default expansion interval is from 1969 to 2068.
%
%-Exceptions
%
%   1)  If `year' is less than 1, no action is taken.
%
%   2)  If the input argument `year' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `year' is not of the expected type, or
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
%   MICE.REQ
%   TIME.REQ
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Fixed bug: Added check in underlying code for "year" to be positive in
%       order to update the lower bound for the expansion.
%
%       Updated the header to comply with NAIF standard. Added
%       complete code example to -Examples section.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 13-FEB-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 12-JAN-2006 (EDW)
%
%-Index_Entries
%
%   Set the interval of expansion for abbreviated years
%
%-&
function cspice_tsetyr( year )

   switch nargin
      case 1

         year = zzmice_int(year);

      otherwise

         error ( 'Usage: cspice_tsetyr( year )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('tsetyr_c', year);
   catch spiceerr
      rethrow(spiceerr)
   end


