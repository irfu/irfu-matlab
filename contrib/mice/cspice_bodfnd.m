%-Abstract
%
%   CSPICE_BODFND determines whether values exist for some item for any body
%   in the kernel pool.
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
%      body   the ID code of the body for which the item is
%             requested. Bodies are numbered according to the
%             standard NAIF numbering scheme.
%
%             [1,1]   = size(body); int32 = class(body)
%
%
%      item   the item to be returned. Together, the body and
%             item name combine to form a variable name, e.g.,
%
%                "BODY599_RADII"
%                "BODY4_POLE_RA"
%
%             [1,c1] = size(item); char = class(item)
%
%                or
%
%             [1,1] = size(item); cell = class(item)
%
%   the call:
%
%      [found] = cspice_bodfnd( body, item )
%
%   returns:
%
%      found    true if the  `item` is in the kernel pool; false
%               if it is not.
%
%               [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example:
%
%      %
%      % Test if Earth's radii values exist in the
%      % kernel pool.
%      %
%      % The procedure searches for the kernel variable
%      % BODY399_RADII.
%      %
%      found = cspice_bodfnd( 399, 'RADII' );
%   
%      %
%      % If found, retrieve the values.
%      %
%      
%      if (found)
%      
%         radii = cspice_bodvcd( 399, 'RADII', 3 ) 
%
%      else
%
%         radii = [ 100; 100; 100 ]
%
%      end
%
% Matlab outputs
%
%   (If BODY399_RADII exists in the kernel pool)
%
%   radii =
%   
%        6.378136600000000e+03
%        6.378136600000000e+03
%        6.356751900000000e+03
%
%   (If BODY399_RADII does not exist in the kernel pool)
%   
%   radii =
%   
%      100
%      100
%      100
%
%-Particulars
%
%   The Mice routines cspice_bodvcd and cspice_bodvrd, which return values
%   from the kernel pool, signal an error if the specified item is not found.
%   In many cases, this is appropriate. However, sometimes the program
%   may attempt to recover, by providing default values, prompting for
%   replacements, and so on.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine bodfnd_c.
%
%   MICE.REQ
%   KERNEL.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 14-NOV-2016, EDW (JPL)
%
%-Index_Entries
%
%   find constants for a body in the kernel pool
%
%-&

function [found] = cspice_bodfnd( body, item )

   switch nargin
      case 2

         body   = zzmice_int(body);
         item   = zzmice_str(item);

      otherwise

         error ( ['Usage: [found] = cspice_bodfnd( body, `item` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [found] = mice( 'bodfnd_c', body, item );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end
