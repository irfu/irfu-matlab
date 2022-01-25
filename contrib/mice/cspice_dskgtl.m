%-Abstract
%
%   CSPICE_DSKGTL retrieves the value of a specified DSK tolerance
%   or margin parameter.
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
%      keywrd   an integer code specifying the parameter to
%               retrieve.
%
%               [1,1] = size(keywrd); int32 = class(keywrd)
%
%               See the include file MiceDtl.m for a description of the
%               possible keywords.
%
%   the call:
%
%      dpval = cspice_dskgtl( keywrd )
%
%   returns:
%
%      dpval    the value of the parameter specified by `keywrd'.
%
%               [1,1] = size(dpval); double = class(dpval)
%
%-Parameters
%
%   See the include file
%
%      MiceDtl.m
%
%   for descriptions and values of the tolerance or margin parameters
%   accessed by this routine, and of the keyword parameters used to
%   refer to them.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Retrieve the DSK tolerance keys and their corresponding values.
%      Alter the values of those parameters for which it is possible,
%      and confirm the change.
%
%      Example code begins here.
%
%
%      function dskgtl_ex1
%
%         %
%         % Retrieve the DSK tolerance keys and corresponding values.
%         %
%         MiceUser;
%
%         %
%         % Define an array of keys, and an array of key names.
%         %
%         dskkeys = [ SPICE_DSK_KEYXFR, SPICE_DSK_KEYSGR,                  ...
%                     SPICE_DSK_KEYSPM, SPICE_DSK_KEYPTM,                  ...
%                     SPICE_DSK_KEYAMG, SPICE_DSK_KEYLAL ];
%
%         knames = { 'SPICE_DSK_KEYXFR', 'SPICE_DSK_KEYSGR',               ...
%                    'SPICE_DSK_KEYSPM', 'SPICE_DSK_KEYPTM',               ...
%                    'SPICE_DSK_KEYAMG', 'SPICE_DSK_KEYLAL' };
%
%
%         %
%         % Output the tolerance keys and values.
%         %
%         for i=1:numel( dskkeys )
%
%            dpval = cspice_dskgtl( dskkeys(i) );
%
%            fprintf( 'Key %s, key value %d, parameter value %d\n',        ...
%                             char( knames(i) ), dskkeys(i), dpval )
%
%         end
%
%         fprintf( '\n' )
%         %
%         % Alter the values of the parameters. Note, you cannot change
%         % parameters SPICE_DSK_KEYAMG and SPICE_DSK_KEYLAL.
%         % Confirm the change.
%         %
%         for i=[1:4]
%
%            cspice_dskstl( dskkeys(i), i*10 );
%            dpval = cspice_dskgtl( dskkeys(i) );
%
%            fprintf( 'Key %s, key value %d, new parameter value %d\n',    ...
%                                char( knames(i) ), dskkeys(i), dpval )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Key SPICE_DSK_KEYXFR, key value 1, parameter value 1e-10
%      Key SPICE_DSK_KEYSGR, key value 2, parameter value 1e-08
%      Key SPICE_DSK_KEYSPM, key value 3, parameter value 1e-10
%      Key SPICE_DSK_KEYPTM, key value 4, parameter value 1e-07
%      Key SPICE_DSK_KEYAMG, key value 5, parameter value 1e-12
%      Key SPICE_DSK_KEYLAL, key value 6, parameter value 1e-12
%
%      Key SPICE_DSK_KEYXFR, key value 1, new parameter value 10
%      Key SPICE_DSK_KEYSGR, key value 2, new parameter value 20
%      Key SPICE_DSK_KEYSPM, key value 3, new parameter value 30
%      Key SPICE_DSK_KEYPTM, key value 4, new parameter value 40
%
%
%-Particulars
%
%   The DSK tolerance routines centralize numeric tolerance and margin
%   values used by the DSK subsystem. The DSK subsystem retrieves values
%   from the DSK tolerance subsystem to use at run time.
%
%   The DSK tolerance access functions are
%
%      cspice_dskgtl {DSK, get tolerance value}
%      cspice_dskstl {DSK, set tolerance value}
%
%   To minimize run time overhead, the "keywords" used by these routines
%   to identify parameters are actually integer codes.
%
%   SPICE users may override certain values maintained by this subsystem;
%   others values are fixed. It is recommended that any change to the
%   tolerance values made at run time be performed only by expert SPICE
%   users.
%
%-Exceptions
%
%   1)  If the input keyword is not recognized, the error
%       SPICE(INDEXOUTOFRANGE) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the input argument `keywrd' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `keywrd' is not of the expected type, or
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
%   DSK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 21-JUL-2020 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Added example's problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 10-MAR-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   retrieve DSK tolerance or margin parameters
%
%-&

function [dpval] = cspice_dskgtl( keywrd )

   switch nargin
      case 1

         %
         % The maximum value, 6, defined in the range argument corresponds to
         % the number of of parameters assigned in MiceDtl.m.
         %

         keywrd  = zzmice_int(keywrd, [1,6]);

      otherwise

         error ( 'Usage: [dpval] = cspice_dskgtl( keywrd )' )

   end

   %
   % Call the MEX library.
   %
   try
      [dpval] = mice( 'dskgtl_c', keywrd );
   catch spiceerr
      rethrow(spiceerr)
   end


