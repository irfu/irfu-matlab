%-Abstract
%
%   CSPICE_ET2UTC converts an input time from ephemeris seconds
%   past J2000 to Calendar, Day-of-Year, or Julian Date format, UTC.
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
%      et       the input epoch(s), ephemeris seconds past J2000.
%
%               [1,n] = size(et); double = class(et)
%
%      format   the format of the output time string.
%
%               [1,c1] = size(format); char = class(format)
%
%                  or
%
%               [1,1] = size(format); cell = class(format)
%
%               It may be any of the following:
%
%                  'C'      Calendar format, UTC.
%
%                  'D'      Day-of-Year format, UTC.
%
%                  'J'      Julian Date format, UTC.
%
%                  'ISOC'   ISO Calendar format, UTC.
%
%                  'ISOD'   ISO Day-of-Year format, UTC.
%
%      prec     the number of digits of precision to which fractional seconds
%               (for Calendar and Day-of-Year formats) or days (for Julian
%               Date format) are to be computed.
%
%               [1,1] = size(prec); int32 = class(prec)
%
%               If `prec' is zero or smaller, no decimal point is appended
%               to the output string. If `prec' is greater than 14, it is
%               treated as 14.
%
%   the call:
%
%      [utcstr] = cspice_et2utc( et, format, prec )
%
%   returns:
%
%      utcstr   the output time string(s) equivalent to the input epoch, in the
%               specified format.
%
%               [1,c2] = size(utcstr); char = class(utcstr)
%
%               Some examples are shown below.
%
%                  'C'      '1986 APR 12 16:31:09.814'
%                  'D'      '1986-102 // 16:31:12.814'
%                  'J'      'JD 2446533.18834276'
%                  'ISOC'   '1987-04-12T16:31:12.814'
%                  'ISOD'   '1987-102T16:31:12.814'
%
%               Fractional seconds, or for Julian dates, fractional
%               days, are rounded to the precision level specified
%               by the input argument `prec'.
%
%               For epochs prior to 1000 A.D. Jan 1 calendar
%               and day of year formats are returned with the
%               era (A.D. or B.C.) attached to the year. For
%               example
%
%                  '877 A.D. MAR 17 13:29:11.829'
%                  '471 B.C. Jan 01 12:00:00.000'
%                  '471 B.C. 001 // 12:00:00.000'
%
%               ISO formats do not support the inclusion of an era.
%               For years prior to 1 A.D. an error will be signaled
%               if ISO format has been requested.
%
%               `utcstr' returns with the same vectorization measure, N,
%               as `et'.
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
%   1) Convert a single ephemeris time to UTC in Julian Date format, and
%      an array of ephemeris times to UTC in Calendar format.
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
%      function et2utc_ex1()
%         %
%         % Define an arbitrary ephemeris time.
%         %
%         et     = -527644192.5403653;
%         format = 'J';
%         prec   = 6;
%         SIZE   = 5;
%
%         %
%         % Load a leapseconds kernel.
%         %
%         cspice_furnsh( 'naif0012.tls' )
%
%         disp( 'Ephemeris time          Converted UTC time' )
%         disp( '---------------------   ---------------------------' )
%
%         %
%         % Convert the ephemeris time to Julian Date
%         % 'format'. Define precision to 6 decimal
%         % places.
%         %
%         utcstr = cspice_et2utc( et, format, prec );
%         disp( 'Scalar (Julian Date format):' )
%         fprintf( '%21.8f   %s\n', et, utcstr )
%
%         %
%         % Create an array of ephemeris times beginning
%         % at -527644192.5403653 with graduations of 10000.0
%         % ephemeris seconds.
%         %
%         et     = [0:(SIZE-1)]*10000. -527644192.5403653;
%         format = 'C';
%
%         %
%         % Convert the array of ephemeris times 'et' to an
%         % array of UTC strings, 'utcstr', in calendar
%         % 'format'.
%         %
%         utcstr= cspice_et2utc( et, format, prec );
%
%         disp( 'Vector (Calendar format):' )
%
%         for n=1:SIZE
%
%            fprintf( '%21.8f   %s\n', et(n), utcstr(n,:) )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Ephemeris time          Converted UTC time
%      ---------------------   ---------------------------
%      Scalar (Julian Date format):
%        -527644192.54036528   JD 2445438.006415
%      Vector (Calendar format):
%        -527644192.54036528   1983 APR 13 12:09:14.274000
%        -527634192.54036528   1983 APR 13 14:55:54.274001
%        -527624192.54036528   1983 APR 13 17:42:34.274001
%        -527614192.54036528   1983 APR 13 20:29:14.274002
%        -527604192.54036528   1983 APR 13 23:15:54.274002
%
%
%-Particulars
%
%   This routine handles the task of converting a double precision
%   representation of an epoch to a character string suitable for
%   human consumption. The more general routine cspice_timout may also be
%   used to convert `et' to time strings.
%
%-Exceptions
%
%   1)  If the format for the output string is not recognized, the
%       error SPICE(INVALIDTIMEFORMAT) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If `prec' is less than or equal to zero, it is treated as
%       zero. If `prec' is greater than 14, it is treated as 14.
%
%   3)  If one of the ISO formats is specified (ISOC or ISOD) but the
%       year corresponding to `et' is prior to 1 A.D. on the Gregorian
%       Calendar, the error SPICE(YEAROUTOFRANGE) is signaled by a
%       routine in the call tree of this routine.
%
%   4)  Epochs prior to 15 Oct, 1582 on the Gregorian calendar (the
%       calendar commonly used in western societies) are returned in
%       the "extended" Gregorian Calendar. To convert epochs to the
%       Julian calendar see the SPICELIB routine GR2JUL.
%
%   5)  This routine does not attempt to account for variations
%       in the length of the second that were in effect prior
%       to Jan 1, 1972. For days prior to that date, we assume
%       there are exactly 86400 ephemeris seconds. Consequently
%       the UTC Gregorian calendar strings produced for epochs
%       prior to Jan 1, 1972 differ from the corresponding
%       TDB calendar strings by approximately 41.18 seconds.
%       (TDB Gregorian calendar strings are produced by the
%       routine cspice_etcal).
%
%   6)  If a leapseconds kernel has not been loaded prior to calling
%       this routine, an error is signaled by a routine in the
%       call tree of this routine.
%
%   7)  If any of the input arguments, `et', `format' or `prec', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   8)  If any of the input arguments, `et', `format' or `prec', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   A leapseconds kernel must be loaded via cspice_furnsh prior to calling
%   this routine. The kernel need be loaded only once during a program
%   run.
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
%   [1]  J. Jespersen and J. Fitz-Randolph, "From Sundials to Atomic
%        Clocks, Understanding Time and Frequency," Dover
%        Publications, Inc. New York, 1982.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and meta-kernel. Modified example's
%       output.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 05-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   ephemeris time to utc
%
%-&

function [utcstr] = cspice_et2utc(et, format, prec )

   switch nargin
      case 3

         et     = zzmice_dp(et);
         format = zzmice_str(format);
         prec   = zzmice_int(prec);

      otherwise

         error ( 'Usage: [_`utcstr`_] = cspice_et2utc(_et_, `format`, prec)' )

   end

   %
   % Call the MEX library.
   %
   try
      [utcstr] = mice('et2utc_c',et,format,prec);
   catch spiceerr
      rethrow(spiceerr)
   end
