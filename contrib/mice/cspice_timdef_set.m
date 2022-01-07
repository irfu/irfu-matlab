%-Abstract
%
%   CSPICE_TIMDEF_SET sets the default zone/calendar/system definitions
%   associated with calendar input strings.
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
%      item     the Time subsystem parameter value to set.
%
%               [1,c1] = size(item); char = class(item)
%
%                  or
%
%               [1,1] = size(item); cell = class(item)
%
%               The parameters and corresponding values:
%
%                  item        Allowed Values
%                  ---------   --------------
%                  CALENDAR    GREGORIAN
%                              JULIAN
%                              MIXED
%
%                  SYSTEM      TDB
%                              TDT
%                              TT
%                              UTC
%
%                  ZONE        EST, EDT, CST, CDT, MST, MDT, PST, PDT
%                              UTC+HR
%                              UTC-HR       ( 0 <= HR < 13 )
%                              UTC+HR:MN    ( 0 <= MN < 60 )
%                              UTC-HR:MN
%
%               The case of `item' is not significant.
%
%      value    the value to associate to `item'.
%
%               [1,1] = size(value); double = class(value)
%
%               Note that value is checked to ensure it is within the range
%               of allowed values for item. If it is not within the expected
%               range and appropriate error message signals.
%
%               The case of `value' is not significant.
%
%   the call:
%
%      cspice_timdef_set( item, value )
%
%   associates `item' to `value'.
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
%   1) Set the calendar parameter of the time system calendar string
%      processing to 'MIXED' mode.
%
%      Example code begins here.
%
%
%      function timdef_set_ex1()
%
%         ITEMS = { 'CALENDAR', 'SYSTEM' };
%
%         %
%         % Set the calendar to 'MIXED'
%         %
%         ITEM  = 'CALENDAR';
%         VALUE = 'MIXED';
%
%         cspice_timdef_set( ITEM, VALUE )
%
%         %
%         % Retrieve the time definition settings.
%         %
%         for i=1:numel(ITEMS)
%
%            value = cspice_timdef_get( ITEMS(i) );
%            fprintf( '%s -> %s\n',  char(ITEMS(i)), value )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      CALENDAR -> MIXED
%      SYSTEM -> UTC
%
%
%-Particulars
%
%   The routines cspice_timdef_get and cspice_timdef_set exist
%   to allow SPICE toolkit users to alter the default interpretation
%   of time strings made by the routine cspice_str2et.
%
%   Normally, unlabeled time strings are assumed to belong to
%   the Gregorian Calendar and are UTC times. However, you
%   may alter the default behavior by calling cspice_timdef_set.
%
%   Calendar
%   --------
%
%   You may set the calendar to be one of the following
%
%   Gregorian   --- This is the calendar used daily the
%                   Western Hemisphere. Leap years occur in this
%                   calendar every 4 years except on centuries
%                   such as 1900 that are not divisible by 400.
%
%   Julian      --- This is the calendar that was in use prior
%                   to October 15, 1582. Leap years occur every
%                   4 years on the Julian Calendar (including all
%                   centuries.) October 5, 1582 on the Julian
%                   calendar corresponds to October 15, 1582 of the
%                   Gregorian Calendar.
%
%   Mixed       --- This calendar uses the Julian calendar
%                   for days prior to October 15, 1582 and
%                   the Gregorian calendar for days on or after
%                   October 15, 1582.
%
%   To set the default calendar, select on of the above for value
%   and make the following call.
%
%      cspice_timdef_set( 'CALENDAR', value )
%
%
%   System
%   -------
%
%   You may set the system used for keeping time to be UTC (default)
%   TDB (barycentric Dynamical Time), TDT (Terrestrial Dynamical
%   Time), or TT (Terrestrial Time). TDT and TT represent the same
%   time system. Both TDB and TT (TDT) have no leapseconds. As such
%   the time elapsed between any two epochs on these calendars does
%   not depend upon when leapseconds occur.
%
%   To set the default time system, select TDT, TT, TDB or UTC for
%   `value' and make the following call.
%
%      cspice_timdef_set( 'SYSTEM', value )
%
%   Note that such a call has the side effect of setting the value
%   associated with 'ZONE' to a blank.
%
%   Zone
%   -----
%
%   You may alter the UTC system by specifying a time zone (UTC
%   offset). For example you may specify that epochs are referred
%   to Pacific Standard Time (PST --- UTC-7). The standard
%   abbreviations for U.S. time zones are recognized:
%
%      EST   UTC-5
%      EDT   UTC-4
%      CST   UTC-6
%      CDT   UTC-5
%      MST   UTC-7
%      MDT   UTC-6
%      PST   UTC-8
%      PDT   UTC-7
%
%   In addition you may specify any commercial time zone by using
%   "offset" notation. This notation starts with the letters "UTC"
%   followed by a + for time zones east of Greenwich and - for time
%   zones west of Greenwich. This is followed by the number of hours
%   to add or subtract from UTC. This is optionally followed by a
%   colon ":" and the number of minutes to add or subtract (based on
%   the sign that follows "UTC") to get the local time zone. Thus to
%   specify the time zone of Calcutta you would specify the time zone
%   to be UTC+5:30. To specify the time zone of Newfoundland use the
%   time zone UTC-3:30.
%
%   To set a default time zone, select one of the "built-in" U.S.
%   zones or construct an offset as discussed above. Then make the
%   call
%
%      cspice_timdef_set( 'ZONE', value );
%
%   If you "GET" a 'ZONE' it will either be blank, or have the
%   form "UTC+/-HR[:MN]"
%
%   Note that such a call has the side effect of setting the value
%   associated with SYSTEM to a blank.
%
%-Exceptions
%
%   1)  If the `item' specified is not one the recognized items, the
%       error SPICE(BADTIMEITEM) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the value associated with the item input is not one of
%       the recognized items, the error SPICE(BADDEFAULTVALUE) is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If any of the input arguments, `item' or `value', is undefined,
%       an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `item' or `value', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
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
%       Added time system name 'TT' (Terrestrial Time) as alternate
%       assignment of 'TDT' (Terrestrial Dynamical Time).
%
%       Updated the header to comply with NAIF standard. Added complete
%       code example based on existing fragment and removed non applicable
%       example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
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
%   -Mice Version 1.0.0, 14-NOV-2013 (EDW)
%
%-Index_Entries
%
%   Change time software defaults.
%   Time Zones
%   Gregorian and Julian Calendars
%
%-&

function cspice_timdef_set( item, value )

   switch nargin
      case 2

         item  = zzmice_str(item);
         value = zzmice_str(value);

      otherwise

         error ( 'Usage: cspice_timdef_set( `item`, `value`)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'timdef_set_c', item, value );
   catch spiceerr
      rethrow(spiceerr)
   end



