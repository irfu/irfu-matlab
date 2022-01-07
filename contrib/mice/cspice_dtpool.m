%-Abstract
%
%   CSPICE_DTPOOL returns descriptive data about a kernel pool variable
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
%      name     name(s) of variables whose values are to be returned.
%
%               [n,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,n] = size(name); cell = class(name)
%
%   the call:
%
%      [found, n, type] = cspice_dtpool( name )
%
%   returns:
%
%      found    flag(s) returning as true if the variable `name' exists in the
%               pool; false if not.
%
%               [1,n] = size(found); logical = class(found)
%
%      n        the number of values associated with `name'. If `name' does
%               not exist in the pool, `n' returns with the value 0.
%
%               [1,n] = size(n); int32 = class(n)
%
%      type     indicating the variable type associated with `name'.
%
%               [n,1] = size(type); char = class(type)
%
%                   C if the data is character data
%                   N if the data is numeric
%                   X if there is no variable name in the pool
%
%               `found', `n', and `type' return with the same vectorization
%               measure, N, as `name'.
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
%   1) Check for the variables defined in the leapseconds kernel and
%      a name probably (hopefully) not in the kernel pool.
%
%      Use the LSK kernel below as test file to generate the results.
%
%         naif0009.tls
%
%
%      Example code begins here.
%
%
%      function dtpool_ex1()
%
%         %
%         % Load a leapsecond kernel.
%         %
%         cspice_furnsh('naif0009.tls' )
%
%         %
%         % Check for the variables defined in the leapseconds kernel
%         % and a name probably (hopefully) not in the kernel pool.
%         %
%         lmpoolNames  = strvcat(              ...
%                       'DELTET/DELTA_T_A',    ...
%                       'DELTET/K',            ...
%                       'DELTET/EB',           ...
%                       'DELTET/M',            ...
%                       'ECHO419',             ...
%                       'DELTET/DELTA_AT',     ...
%                       'EVERLASTING_GOBSTOPPER' );
%
%         [found, n, dtype] = cspice_dtpool( lmpoolNames );
%
%         for i = 1:size(lmpoolNames,1)
%
%            name = lmpoolNames(i,:);
%
%            if (found(i))
%               fprintf( 'Variable name : %s\n', name       )
%               fprintf( 'Variable size : %d\n', n(i)       )
%               fprintf( 'Variable type : %s\n\n', dtype(i) )
%            else
%               fprintf( 'Unable to find variable name : %s\n\n', name )
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Variable name : DELTET/DELTA_T_A
%      Variable size : 1
%      Variable type : N
%
%      Variable name : DELTET/K
%      Variable size : 1
%      Variable type : N
%
%      Variable name : DELTET/EB
%      Variable size : 1
%      Variable type : N
%
%      Variable name : DELTET/M
%      Variable size : 2
%      Variable type : N
%
%      Unable to find variable name : ECHO419
%
%      Variable name : DELTET/DELTA_AT
%      Variable size : 50
%      Variable type : N
%
%      Unable to find variable name : EVERLASTING_GOBSTOPPER
%
%
%-Particulars
%
%   A sister version of this routine exists named mice_dtpool that returns
%   the output arguments as fields in a single structure.
%
%   This routine allows you to determine whether or not a kernel
%   pool variable is present and to determine its size and type
%   if it is.
%
%-Exceptions
%
%   1)  If the name requested is not in the kernel pool `found'
%       will be set to false, `n' to zero and `type' to 'X'.
%
%   2)  If the input argument `name' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `name' is not of the expected type, or
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added -Parameters,
%       -Exceptions, -Files, -Restrictions, -Literature_References and
%       -Author_and_Institution sections. Updated -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 03-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.2, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Edits to Example section, proper description of "standard.tm"
%       meta kernel.
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 07-MAR-2007 (EDW)
%
%-Index_Entries
%
%   return summary information about a kernel pool variable
%
%-&

function [found, n, type] = cspice_dtpool(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [_found_, _n_, _`type`_] = cspice_dtpool(_`name`_)' )

   end

%
% Call the MEX library. The "_s" suffix indicates a structure type
% return argument.
%
   %
   % Call the MEX library.
   %
   try
      [dtpool] = mice('dtpool_s',name);
      found    = reshape( [dtpool.found], 1, [] );
      n        = reshape( [dtpool.n],     1, [] );
      type     = char( dtpool.type );
   catch spiceerr
      rethrow(spiceerr)
   end
