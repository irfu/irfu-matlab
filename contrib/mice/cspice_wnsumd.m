%-Abstract
%
%   CSPICE_WNSUMD summarizes the contents of a double precision window.
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
%      window   SPICE window containing zero or more intervals.
%
%               [2n,1] = size(window); double = class(window)
%
%   the call:
%
%      [meas, avg, stddev, idxsml, idxlon] = cspice_wnsumd( window )
%
%   returns:
%
%      meas     total measure of the intervals in the input  window. This is
%               just the sum of the measures of the individual intervals.
%
%               [1,1] = size(meas); double = class(meas)
%
%      avg      average measure of the intervals in the input window.
%
%               [1,1] = size(avg); double = class(avg)
%
%      stddev   standard deviation of the measures of the intervals in the
%               input window.
%
%               [1,1] = size(stddev); double = class(stddev)
%
%      idxsml,
%      idxlon   indices of the left endpoint of, respectively, the shortest
%               and longest intervals in the data contained in `window'.
%
%               [1,1] = size(idxlon); int32 = class(idxlon)
%               [1,1] = size(idxsml); int32 = class(idxsml)
%
%               The following algorithm describes the relation of
%               `idxsml' and `idxlon' to the window data:
%
%                  The shortest interval:
%
%                     [ window(idxsml), window(idxsml+1) ]
%
%                  The longest interval:
%
%                     [ window(idxlon), window(idxlon+1) ]
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
%   1) Define an array representing a window with six intervals, and
%      calculate the summary for that window.
%
%      Example code begins here.
%
%
%      function wnsumd_ex1()
%
%         %
%         % Define an array representing a window with six intervals.
%         % The values in `window' have correct order for a
%         % SPICE window.
%         %
%         window = [ [  1.;  3.]; ...
%                    [  7.; 11.]; ...
%                    [ 18.; 18.]; ...
%                    [ 23.; 27.]; ...
%                    [ 30.; 69.]; ...
%                    [ 72.; 80.] ];
%
%         %
%         % Calculate the summary for `window'.
%         %
%         [meas, avg, stddev, idxsml, idxlon] = cspice_wnsumd( window );
%
%         %
%         % `idxsml' and `idxlon' refer to the indices of
%         % the "cell" data array.
%         %
%         intrvl_short= (idxsml+1)/2;
%         intrvl_long = (idxlon+1)/2;
%
%         fprintf( 'Measure           : %f\n', meas        )
%         fprintf( 'Average           : %f\n', avg         )
%         fprintf( 'Standard Dev      : %f\n', stddev      )
%         fprintf( 'Index shortest    : %i\n', idxsml      )
%         fprintf( 'Index longest     : %i\n', idxlon      )
%         fprintf( 'Interval shortest : %i\n', intrvl_short)
%         fprintf( 'Interval longest  : %i\n', intrvl_long )
%
%         fprintf( 'Shortest interval : %f %f\n', window(idxsml),  ...
%                                                 window(idxsml+1) )
%
%         fprintf( 'Longest interval  : %f %f\n', window(idxlon), ...
%                                                 window(idxlon+1) )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Measure           : 57.000000
%      Average           : 9.500000
%      Standard Dev      : 13.413302
%      Index shortest    : 5
%      Index longest     : 9
%      Interval shortest : 3
%      Interval longest  : 5
%      Shortest interval : 18.000000 18.000000
%      Longest interval  : 30.000000 69.000000
%
%
%   2) Let A contain the intervals
%
%         [ 1, 3 ]  [ 7, 11 ]  [ 23, 27 ]
%
%      Let B contain the singleton intervals
%
%         [ 2, 2 ]  [ 9, 9 ]  [ 27, 27 ]
%
%      The measures of A and B are
%
%         (3-1) + (11-7) + (27-23) = 10
%
%      and
%
%         (2-2) + (9-9) + (27-27) = 0
%
%      respectively. Each window has three intervals; thus, the average
%      measures of the windows are 10/3 and 0. The standard deviations
%      are
%
%           ------------------------------------------
%          |       2         2          2
%          |  (3-1)  + (11-7)  + (27-23)           2           1/2
%          |  ---------------------------  - (10/3)     = (8/9)
%          |             3
%        \ |
%         \|
%
%      and 0. Neither window has one "shortest" interval or "longest"
%      interval; so the first ones found are returned: `idxsml' and
%      `idxlon' are 0 and 2 for A, 0 and 0 for B.
%
%-Particulars
%
%   This routine provides a summary of the input window, consisting
%   of the following items:
%
%   -  The measure of the window.
%
%   -  The average and standard deviation of the measures
%      of the individual intervals in the window.
%
%   -  The indices of the left endpoints of the shortest
%      and longest intervals in the window.
%
%   All of these quantities are zero if the window contains no
%   intervals.
%
%-Exceptions
%
%   1)  If `window' has odd cardinality, the error
%       SPICE(INVALIDCARDINALITY) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  Left endpoints of stored intervals must be strictly greater
%       than preceding right endpoints. Right endpoints must be
%       greater than or equal to corresponding left endpoints.
%       Invalid window data are not diagnosed by this routine and may
%       lead to unpredictable results.
%
%   3)  If the input argument `window' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `window' is not of the expected type, or
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
%       Changed output argument names "shortest" and "longest" to "idxsml"
%       and "idxlon".
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement, and a second example.
%
%       Corrected error message format.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 15-DEC-2008 (EDW)
%
%-Index_Entries
%
%   summary of a d.p. window
%
%-&

function  [meas, avg, stddev, idxsml, idxlon] = cspice_wnsumd( window )

   switch nargin

      case 1

         window = zzmice_win( window );

      otherwise

         error( [ 'Usage: [meas, avg, stddev, idxsml, idxlon] ' ...
                  ' = cspice_wnsumd( window )' ] )

      end

   try
      [wnsumd] = mice( 'wnsumd_s',  [zeros(6,1);window] );

      meas     = reshape( [wnsumd.meas],   1, [] );
      avg      = reshape( [wnsumd.avg],    1, [] );
      stddev   = reshape( [wnsumd.stddev], 1, [] );
      idxsml   = reshape( [wnsumd.idxsml], 1, [] );
      idxlon   = reshape( [wnsumd.idxlon], 1, [] );

   catch spiceerr
      rethrow(spiceerr)
   end
