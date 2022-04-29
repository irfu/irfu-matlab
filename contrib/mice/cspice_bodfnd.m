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
%      body     the ID code of the body for which the item is
%               requested.
%
%               [1,1]   = size(body); int32 = class(body)
%
%               Bodies are numbered according to the standard NAIF numbering
%               scheme.
%
%      item     the item to be returned.
%
%               [1,c1] = size(item); char = class(item)
%
%                  or
%
%               [1,1] = size(item); cell = class(item)
%
%               Together, the body and item name combine to form a variable
%               name, e.g.,
%
%                  'BODY599_RADII'
%                  'BODY4_POLE_RA'
%
%   the call:
%
%      [bodfnd] = cspice_bodfnd( body, item )
%
%   returns:
%
%      bodfnd   true if the `item' is in the kernel pool; false
%               if it is not.
%
%               [1,1] = size(bodfnd); logical = class(bodfnd)
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
%   1) Test if the Earth's radii values exist in the kernel pool.
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model for the Earth.
%
%         pck00008.tpc
%
%
%      Example code begins here.
%
%
%      function bodfnd_ex1()
%         %
%         % Load a PCK file.
%         %
%         cspice_furnsh( 'pck00008.tpc' );
%
%         %
%         % Test if Earth's radii values exist in the
%         % kernel pool. The procedure searches for the kernel variable
%         % BODY399_RADII.
%         %
%         found = cspice_bodfnd( 399, 'RADII' );
%
%         %
%         % If found, retrieve the values.
%         %
%         if (found)
%
%            radii = cspice_bodvcd( 399, 'RADII', 3 )
%
%         else
%
%            fprintf( ['The kernel pool does not ', ...
%                      'contain Earth''s radii values'] )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      radii =
%
%         6378.1
%         6378.1
%         6356.8
%
%
%-Particulars
%
%   The Mice routines cspice_bodvcd and cspice_bodvrd, which return values
%   from the kernel pool, signal an error if the specified item is not found.
%   In many cases, this is appropriate. However, sometimes the program
%   may attempt to recover, by providing default values, prompting for
%   replacements, and so on.
%
%-Exceptions
%
%   1)  If any of the input arguments, `body' or `item', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   2)  If any of the input arguments, `body' or `item', is not of the
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
%   KERNEL.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
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
%   -Mice Version 1.1.0, 21-JUL-2020 (EDW) (JDR)
%
%       Changed output argument name "found" to "bodfnd" to comply with
%       NAIF standard.
%
%       Edited the header to comply with NAIF standard. Added -Parameters,
%       -Exceptions, -Files, -Restrictions, -Literature_References and
%       -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 14-NOV-2016 (EDW)
%
%-Index_Entries
%
%   find constants for a body in the kernel pool
%
%-&

function [bodfnd] = cspice_bodfnd( body, item )

   switch nargin
      case 2

         body   = zzmice_int(body);
         item   = zzmice_str(item);

      otherwise

         error ( ['Usage: [bodfnd] = cspice_bodfnd( body, `item` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [bodfnd] = mice( 'bodfnd_c', body, item );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      [bodfnd] = zzmice_logical(bodfnd);
   catch spiceerr
      rethrow(spiceerr)
   end
