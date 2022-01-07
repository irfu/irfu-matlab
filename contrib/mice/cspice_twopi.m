%-Abstract
%
%   CSPICE_TWOPI returns the value of the constant 2*pi.
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
%   The call:
%
%      [twopi] = cspice_twopi
%
%   returns:
%
%      twopi    twice the value of pi (the ratio of a circle's circumference
%               to its diameter), determined by the acos function.
%
%               [1,1] = size(twopi); double = class(twopi)
%
%               That is,
%
%                     twopi = acos ( -1.0 ) * 2.0;
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
%      the constant 2*PI and prints it out.
%
%      Example code begins here.
%
%
%      function twopi_ex1()
%
%         %
%         % Print the double precision value of 2*Pi
%         %
%         fprintf( '2*Pi: %25.22f\n', cspice_twopi )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      2*Pi:  6.2831853071795862319959
%
%
%-Particulars
%
%   The first time the function is referenced, the value is computed
%   as shown above. The value is saved, and returned directly upon
%   subsequent reference.
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
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel and reformatted example's output.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
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
%   twice the value of pi
%
%-&

function [twopi] = cspice_twopi

   switch nargin
      case 0
         ;
      otherwise
         error ( 'Usage: [twopi] = cspice_twopi' )
   end

   %
   % Call the MEX library.
   %
   try
      [twopi] =  mice('twopi_c');
   catch spiceerr
      rethrow(spiceerr)
   end

