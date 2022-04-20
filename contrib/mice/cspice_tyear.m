%-Abstract
%
%   CSPICE_TYEAR returns the double precision value of the number
%   of seconds in a tropical year: 31556925.9747
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
%      None.
%
%   the call:
%
%      [tyear] = cspice_tyear
%
%   returns:
%
%      tyear   the number of seconds per tropical year. This value is
%              taken from the 1992 Explanatory Supplement to the
%              Astronomical Almanac.
%
%              [1,n] = size(tyear); double = class(tyear)
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
%   1) The following code example returns the double precision value of
%      the number of seconds in a tropical year, and prints it out.
%
%      Example code begins here.
%
%
%      function tyear_ex1()
%
%         %
%         % Print the double precision value of seconds in a tropical
%         % year.
%         %
%         fprintf( 'Seconds per tropical year: %20.8f\n', cspice_tyear )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Seconds per tropical year:    31556925.97470000
%
%
%-Particulars
%
%   The tropical year is often used as a fundamental unit
%   of time when dealing with older ephemeris data. For this
%   reason its value in terms of ephemeris seconds is
%   recorded in this function.
%
%-Exceptions
%
%   Error free.
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
%   [1]  P. Kenneth Seidelmann (Ed.), "Explanatory Supplement to the
%        Astronomical Almanac," p 80, University Science Books, 1992.
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
%       Corrected the value given for the number of seconds in the
%       tropical year in the -Abstract section.
%
%       Edited -Examples section to comply with NAIF standard. Added
%       example's problem statement.
%
%       Changed output argument name "return_val" to "tyear" to comply
%       with NAIF standard.
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
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   Number of seconds per tropical year
%
%-&

function [tyear] = cspice_tyear

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [tyear] = cspice_tyear' )

   end

   %
   % Call the MEX library.
   %
   try
      [tyear] =  mice('tyear_c' );
   catch spiceerr
      rethrow(spiceerr)
   end
