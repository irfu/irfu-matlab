%-Abstract
%
%   CSPICE_ETCAL converts from an ephemeris epoch measured in seconds past
%   the epoch of J2000 to a calendar string format using a
%   formal calendar free of leapseconds.
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
%      et       the ephemeris time(s) expressed as ephemeris seconds
%               past J2000.
%
%               [1,n] = size(et); double = class(et)
%
%   the call:
%
%      [calstr] = cspice_etcal( et )
%
%   returns:
%
%      calstr   the array of time string(s) representing the input ephemeris
%               epoch `et'.
%
%               [n,c1] = size(calstr); char = class(calstr)
%
%               This string is based upon extending the Gregorian Calendar
%               backward and forward indefinitely keeping the same rules for
%               determining leap years. Moreover, there is no accounting for
%               leapseconds.
%
%               The string will have the following format
%
%                  year (era) mon day hr:mn:sc.sss
%
%               Where:
%
%                  year --- is the year
%                  era  --- is the chronological era associated with
%                           the date. For years after 999 A.D.
%                           the era is omitted. For years
%                           between 1 A.D. and 999 A.D. (inclusive)
%                           era is the string 'A.D.' For epochs
%                           before 1 A.D. Jan 1 00:00:00, era is
%                           given as 'B.C.' and the year is converted
%                           to years before the "Christian Era".
%                           The last B.C. epoch is
%
%                             1 B.C. DEC 31 23:59:59.999
%
%                           The first A.D. epoch (which occurs .001
%                           seconds after the last B.C. epoch) is:
%
%                              1 A.D. JAN 1 00:00:00.000
%
%                           Note: there is no year 0 A.D. or 0 B.C.
%                  mon  --- is a 3-letter abbreviation for the month
%                           in all capital letters.
%                  day  --- is the day of the month
%                  hr   --- is the hour of the day (between 0 and 23)
%                           leading zeros are added to hr if the
%                           numeric value is less than 10.
%                  mn   --- is the minute of the hour (0 to 59)
%                           leading zeros are added to mn if the
%                           numeric value is less than 10.
%                  sc.sss   is the second of the minute to 3 decimal
%                           places ( 0 to 59.999). Leading zeros
%                           are added if the numeric value is less
%                           than 10. Seconds are truncated, not
%                           rounded.
%
%               `calstr' returns with the same vectorization measure, N,
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
%   1) Convert a UTC time string to ephemeris time and then this time
%      to string. Note that the later conversion does not require
%      loading a leapseconds kernel.
%
%      Use the LSK kernel below to load the leap seconds and time
%      constants required for the initial conversion from UTC time
%      string to ephemeris time.
%
%         naif0012.tls
%
%
%      Example code begins here.
%
%
%      function etcal_ex1()
%
%         %
%         % Load a leapseconds kernel.
%         %
%         cspice_furnsh( 'naif0012.tls' )
%
%         %
%         % Define a UTC time string.
%         %
%         TIMESTR = '2013 JUN 30 00:00:00.000';
%
%         %
%         % Convert the time string to ephemeris time.
%         %
%         et  = cspice_str2et( TIMESTR );
%
%         %
%         % Convert the ephemeris time to a time string, the conversion
%         % ignoring leapseconds. Note, this evaluation does not require
%         % loading a leapsecond kernel.
%         %
%         cal = cspice_etcal( et );
%
%         %
%         % Display the two time strings.
%         %
%         disp( ['Original times string: ' TIMESTR] )
%         disp( ['ETCAL time string    : ' cal    ] )
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
%      Original times string: 2013 JUN 30 00:00:00.000
%      ETCAL time string    : 2013 JUN 30 00:01:07.184
%
%
%-Particulars
%
%   This is an error free routine for converting ephemeris epochs
%   represented as seconds past the J2000 epoch to formal
%   calendar strings based upon the Gregorian Calendar. This formal
%   time is often useful when one needs a human recognizable
%   form of an ephemeris epoch. There is no accounting for leap
%   seconds in the output times produced.
%
%   Note: The calendar epochs produced are not the same as the
%         UTC calendar epochs that correspond to `et'. The strings
%         produced by this routine may vary from the corresponding
%         UTC epochs by more than 1 minute.
%
%   This routine can be used in creating error messages or
%   in routines and programs in which one prefers to report
%   times without employing leapseconds to produce exact UTC
%   epochs.
%
%-Exceptions
%
%   1)  If the input `et' is so large that the corresponding
%       number of days since 1 A.D. Jan 1, 00:00:00 is
%       within 1 of overflowing or underflowing an integer,
%       `et' will not be converted to the correct string
%       representation rather, the string returned will
%       state that the epoch was before or after the day
%       that is cspice_intmin + 1 or cspice_intmax - 1 days after
%       1 A.D. Jan 1, 00:00:00.
%
%   2)  If the input argument `et' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   3)  If the input argument `et' is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  One must keep in mind when using this routine that
%       ancient times are not based upon the Gregorian
%       calendar. For example the 0 point of the Julian
%       Date system is 4713 B.C. Jan 1, 12:00:00 on the Julian
%       Calendar. If one formalized the Gregorian calendar
%       and extended it indefinitely, the zero point of the Julian
%       date system corresponds to 4714 B.C. NOV 24 12:00:00 on
%       the Gregorian calendar. There are several reasons for this.
%       Leap years in the Julian calendar occur every
%       4 years (including *all* centuries). Moreover,  the
%       Gregorian calendar "effectively" begins on 15 Oct, 1582 A.D.
%       which is 5 Oct, 1582 A.D. in the Julian Calendar.
%
%       Therefore you must be careful in your interpretation
%       of ancient dates produced by this routine.
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
%       Changed the output argument name "string" to "calstr" for consistency
%       with other routines.
%
%       Edited the header to comply with NAIF standard. Extended "calstr"
%       argument description. Added TIME required reading.
%
%       Added example's problem statement and a reference to the required LSK.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.2, 05-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%      Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 07-MAR-2007 (EDW)
%
%-Index_Entries
%
%   Convert ephemeris time to a formal calendar date
%
%-&

function [calstr] = cspice_etcal(et)

   switch nargin
      case 1

         et = zzmice_dp(et);

      otherwise

         error( 'Usage: [_`calstr`_] = cspice_etcal(_et_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [calstr] = mice('etcal_c', et );
   catch spiceerr
      rethrow(spiceerr)
   end



