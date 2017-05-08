%-Abstract
%
%   CSPICE_DSKSTL sets the value of a specified DSK tolerance or 
%   margin parameter.
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
%      keywrd     is an integer code specifying the parameter to set. See
%                 the include file dsktol.inc for a description of the
%                 possible keywords.
%   
%                 [1,1] = size(keywrd); int32 = class(keywrd)
%   
%      dpval      is the new value of the parameter specified by `keywrd'. 
%   
%                 [1,1] = size(dpval); double = class(dpval)
%    
%                 <<< Use extreme caution. This routine performs no
%                 checks on `dpval'. >>>
%
%   the call:
%
%      cspice_dskstl( keywrd, dpval )
%
%   returns:
%
%      None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      function dskgtl_t
%      
%         %
%         % Retrieve the DSK tolerance keys and corresponding values.
%         %
%         MiceUser
%      
%      
%         %
%         % Define an array of keys, and an array of key names.
%         %
%         dskkeys = [ SPICE_DSK_KEYXFR, SPICE_DSK_KEYSGR,  ...
%                     SPICE_DSK_KEYSPM, SPICE_DSK_KEYPTM,  ...
%                     SPICE_DSK_KEYAMG, SPICE_DSK_KEYLAL ];
%      
%         knames = { 'SPICE_DSK_KEYXFR', 'SPICE_DSK_KEYSGR',  ...
%                    'SPICE_DSK_KEYSPM', 'SPICE_DSK_KEYPTM',  ...
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
%            fprintf( 'Key %s, key value %d, parameter value %d\n', ...
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
%            fprintf( 'Key %s, key value %d, new parameter value %d\n', ...
%                             char( knames(i) ), dskkeys(i), dpval )
%         end
%         
%   Matlab outputs:
%      
%      Key SPICE_DSK_KEYXFR, key value 1, parameter value 1.000000e-10
%      Key SPICE_DSK_KEYSGR, key value 2, parameter value 1.000000e-08
%      Key SPICE_DSK_KEYSPM, key value 3, parameter value 1.000000e-10
%      Key SPICE_DSK_KEYPTM, key value 4, parameter value 1.000000e-07
%      Key SPICE_DSK_KEYAMG, key value 5, parameter value 1.000000e-12
%      Key SPICE_DSK_KEYLAL, key value 6, parameter value 1.000000e-12
%      
%      Key SPICE_DSK_KEYXFR, key value 1, new parameter value 10
%      Key SPICE_DSK_KEYSGR, key value 2, new parameter value 20
%      Key SPICE_DSK_KEYSPM, key value 3, new parameter value 30
%      Key SPICE_DSK_KEYPTM, key value 4, new parameter value 40
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
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dskstl_c.
%
%   MICE.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 10-MAR-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   set dsk tolerance or margin parameters 
%
%-&

function cspice_dskstl( keywrd, dpval )

   switch nargin
      case 2

         %
         % The maximum value, 6, defined in the range argument corresponds to
         % the number of of parameters assigned in dsktol.inc.
         %

         keywrd  = zzmice_int(keywrd, [1,6]);
         dpval   = zzmice_dp(dpval);

      otherwise

         error ( 'Usage: cspice_dskstl( keywrd, dpval )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dskstl_c', keywrd, dpval );
   catch
      rethrow(lasterror)
   end


