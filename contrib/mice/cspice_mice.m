%-Abstract
%
%   CSPICE_MICE returns compile time, date, and version information
%   on the Mice shared object library.
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
%      item   a string identifying the compile value to return.
%
%             [1,c1] = size(item); char = class(item)
%
%                or
%
%             [1,1] = size(item); cell = class(item)
%
%             Allowed values for 'item':
%
%                'DATE'    - the date of the Mice interface compile
%
%                'TIME'    - the time on the date of the Mice
%                            interface compile
%
%                'VERSION' - the Mice version string
%
%   the call:
%
%      value = cspice_mice( item )
%
%   returns:
%
%      value   the string corresponding to the quantity identified by 'item'.
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
%   1) Display the Mice library version, and its compilation date and
%      time, as well as the Mice toolkit version.
%
%      Note that the Mice toolkit version is different than the Mice
%      library shared object version. The earlier is common to all
%      toolkits, while the later is specific to Mice.
%
%      Example code begins here.
%
%
%      function mice_ex1()
%
%         fprintf( 'Mice toolkit version: %s\n', cspice_tkvrsn( 'toolkit' ) )
%         fprintf( 'Mice library version: %s\n', cspice_mice  ( 'version' ) )
%         fprintf( '   Compiled on %s at %s\n', ...
%                  cspice_mice( 'date' ),         ...
%                  cspice_mice( 'time' )         )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Mice toolkit version: CSPICE_N0066
%      Mice library version: Mice 1.5.0 05-JAN-2017 (EDW) (NJB)
%         Compiled on Apr 11 2018 at 23:27:58
%
%
%-Particulars
%
%   None.
%
%-Exceptions
%
%   None.
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard.
%       Reformatted example's output and added problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%   -Mice Version 1.0.1, 13-FEB-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 26-JUN-2006 (EDW)
%
%-Index_Entries
%
%   return compile time, date, and version
%
%-&

function [value] = cspice_mice( item )

   switch nargin
      case 1

         item = zzmice_str(item);

      otherwise

         error ( 'Usage: [`value`] = cspice_mice( `item` )' )

   end


   %
   % Call the MEX library.
   %
   try
      [value] = mice( 'cspice_mice', item );
   catch spiceerr
      rethrow(spiceerr)
   end



