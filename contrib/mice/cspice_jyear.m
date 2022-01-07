%-Abstract
%
%   CSPICE_JYEAR returns the value for the number of seconds in a
%   Julian year.
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
%      No input required.
%
%   the call:
%
%      jyear = cspice_jyear
%
%   returns:
%
%      jyear   the number of seconds per Julian year as a double precision
%              scalar.
%
%              [1,1] = size(jyear); double = class(jyear)
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
%   1) Display the number of seconds in a Julian year.
%
%      Example code begins here.
%
%
%      function jyear_ex1()
%         %
%         % Display the number of seconds in a Julian Year, in
%         % 16.3 floating point format
%         %
%         fprintf( 'Seconds per Julian year: %16.3f\n', cspice_jyear)
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Seconds per Julian year:     31557600.000
%
%
%-Particulars
%
%   The julian year is often used as a fundamental unit of time when
%   dealing with ephemeris data. For this reason its value in terms of
%   ephemeris seconds is recorded in this function.
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
%        Astronomical Almanac," Page 8, University Science Books,
%        1992.
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
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and reformatted example's output.
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
%   -Mice Version 1.0.1, 11-JUN-2013 (EDW)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   Number of seconds per julian year
%
%-&

function [jyear] = cspice_jyear

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [jyear] = cspice_jyear' )

   end

   %
   % Call the MEX library.
   %
   try
      [jyear] =  mice('jyear_c' );
   catch spiceerr
      rethrow(spiceerr)
   end

