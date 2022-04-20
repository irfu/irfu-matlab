%-Abstract
%
%   CSPICE_WNUNID returns the window array union of two double
%   precision window arrays.
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
%      a        SPICE window containing zero or more intervals.
%
%               [2l,1] = size(a); double = class(a)
%
%      b        SPICE window containing zero or more intervals.
%
%               [2m,1] = size(b); double = class(b)
%
%   the call:
%
%      [c] = cspice_wnunid( a, b )
%
%   returns:
%
%      c        SPICE window resulting from the union (in the SPICE sense)
%               of `a' and `b'.
%
%               [2n,1] = size(c); double = class(c)
%
%               The function output can overwrite either `a' or `b'.
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
%   1) Given a set of SPK files, compute the total coverage of the
%      ephemeris data for the Earth system barycenter.
%
%      Use the SPK kernel below as input for Earth system barycenter's
%      ephemeris data covering from 01 Jan 1996 to 02 Jan 2111.
%
%         de403.bsp
%
%      Use the SPK kernel below as input for Earth system barycenter's
%      ephemeris data covering from 01 Jan 1050 to 01 Jan 2050.
%
%         de405.bsp
%
%
%      Example code begins here.
%
%
%      function wnunid_ex1()
%
%         SPK1 = 'de403.bsp';
%         SPK2 = 'de405.bsp';
%
%         %
%         % Retrieve the coverage for body 3 from SPK1
%         %
%         cov1 = cspice_spkcov( SPK1, 3, 10 );
%         fprintf( '%s coverage:\n', SPK1 )
%         fprintf( '   %25.12f   %25.12f\n', cov1)
%
%         %
%         % Retrieve the coverage for body 3 from SPK2
%         %
%         cov2 = cspice_spkcov( SPK2, 3, 10 );
%         fprintf( '%s coverage:\n', SPK2 )
%         fprintf( '   %25.12f   %25.12f\n', cov2)
%
%         %
%         % Perform a windows array union on 'cov1' and 'cov2'
%         %
%         cov3 = cspice_wnunid( cov1, cov2 );
%         fprintf( '\nCombined coverage:\n' )
%         fprintf( '   %25.12f   %25.12f\n', cov3)
%
%         %
%         % The output can overwrite the input.
%         %
%         cov1 = cspice_wnunid( cov1, cov2 );
%         fprintf( '\nCombined coverage (overwriting input):\n' )
%         fprintf( '   %25.12f   %25.12f\n', cov1)
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      de403.bsp coverage:
%           -126273537.816086068749     3502872062.183888912201
%      de405.bsp coverage:
%          -1577879958.816058635712     1577880064.183913230896
%
%      Combined coverage:
%          -1577879958.816058635712     3502872062.183888912201
%
%      Combined coverage (overwriting input):
%          -1577879958.816058635712     3502872062.183888912201
%
%
%-Particulars
%
%   The union of two windows contains every point contained in the
%   first window, or the second window, or both.
%
%-Exceptions
%
%   1)  The cardinality of the input windows must be even. Left
%       endpoints of stored intervals must be strictly greater than
%       preceding right endpoints. Right endpoints must be greater
%       than or equal to corresponding left endpoints. Invalid window
%       data are not diagnosed by this routine and may lead to
%       unpredictable results.
%
%   2)  If any of the input arguments, `a' or `b', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `a' or `b', is not of the
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
%   WINDOWS.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and modified code example to produce
%       formatted output.
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
%   -Mice Version 1.0.1, 08-NOV-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 26-JUN-2007 (EDW)
%
%-Index_Entries
%
%   union two d.p. windows
%
%-&

function [c] = cspice_wnunid( a, b )

   switch nargin

      case 2

         a    = zzmice_win(a);
         b    = zzmice_win(b);

      otherwise

         error( 'Usage: [c] = cspice_wnunid( a, b )' )

   end

   %
   % Call the windows routine, add to 'a' and 'b' the space needed for
   % the control segments.
   %
   try
      [c] = mice( 'wnunid_c', [zeros(6,1); a], [zeros(6,1); b] );
   catch spiceerr
      rethrow(spiceerr)
   end





