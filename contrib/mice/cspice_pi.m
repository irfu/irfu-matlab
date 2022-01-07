%-Abstract
%
%   CSPICE_PI returns the value of the constant pi.
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
%      [onepi] = cspice_pi
%
%   returns:
%
%      onepi    the value of pi (the ratio of a circle's circumference to its
%               diameter), determined by the acos function.
%
%               [1,1] = size(onepi); double = class(onepi)
%
%               That is,
%
%                     onepi = acos ( -1.0 );
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) The following code example returns the double precision value of
%      the constant pi and prints it out.
%
%      Example code begins here.
%
%
%      function pi_ex1()
%
%         %
%         % Print the double precision value of Pi
%         %
%         fprintf( 'Pi: %25.22f\n', cspice_pi )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Pi:  3.1415926535897931159980
%
%
%      Note that the MATLAB system variable "pi" returns a double
%      precision value for PI that equates the value returned by
%      cspice_pi, to machine roundoff.
%
%   2) Suppose that you need to compute the vector that results from
%      off-pointing an instrument boresight vector by pi/10 radians
%      about the Y-axis of the instrument's reference frame.
%
%      Example code begins here.
%
%
%      function pi_ex2()
%
%         %
%         % Let's assume that the instrument boresight is not aligned
%         % to the Z-axis of the instrument's reference frame.
%         %
%         bsight = [ 0.2; 0.04; 1.0 ];
%
%         %
%         % A Pi/10 rotation about the Y axis.
%         %
%         rotmat = cspice_rotate( 0.1*cspice_pi, 2 );
%
%         %
%         % Apply the coordinate rotation to the boresight.
%         %
%         vec = rotmat * bsight;
%         fprintf( 'Pointing vector:\n' )
%         fprintf( '  %16.12f  %16.12f  %16.12f\n', vec );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Pointing vector:
%         -0.118805691116    0.040000000000    1.012859915170
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
%       Edited the header to comply with NAIF standard. Added second
%       example.
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
%       Corrected minor typo in header.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   value of pi
%
%-&

function [onepi] = cspice_pi

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [onepi] = cspice_pi' )

   end

   %
   % Call the MEX library.
   %
   try
      [onepi] =  mice('pi_c');
   catch spiceerr
      rethrow(spiceerr)
   end
