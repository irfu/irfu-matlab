%-Abstract
%
%   CSPICE_J2000 returns the value for the Julian Date of
%   2000 JAN 01 12:00:00 (2000 JAN 1.5).
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
%      [j2000] = cspice_j2000
%
%   returns:
%
%      j2000   the value 2451545.0, the Julian Date corresponding to
%              2000 JAN 01 12:00:00 (2000 JAN 1.5).
%
%              [1,1] = size(j2000); double = class(j2000)
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
%   1) Display the double precision value for the J2000 date.
%
%      Example code begins here.
%
%
%      function j2000_ex1()
%         %
%         % Display the J2000 date in 16.8 format
%         %
%         fprintf( 'J2000 date: %16.8f\n', cspice_j2000)
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      J2000 date: 2451545.00000000
%
%
%-Particulars
%
%   The function always returns the constant value shown above.
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
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement.
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
%   julian date of 2000 jan 1.5
%
%-&

function [j2000] = cspice_j2000

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [j2000] = cspice_j2000' )

   end

   %
   % Call the MEX library.
   %
   try
      [j2000] =  mice('j2000_c' );
   catch spiceerr
      rethrow(spiceerr)
   end


