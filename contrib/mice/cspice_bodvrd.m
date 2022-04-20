%-Abstract
%
%   CSPICE_BODVRD fetches from the kernel pool the double
%   precision values of an item associated with a body.
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
%      bodynm   the body name for which `item' is requested.
%
%               [1,c1] = size(bodynm); char = class(bodynm)
%
%                  or
%
%               [1,1] = size(bodynm); cell = class(bodynm)
%
%               `bodynm' is case-insensitive, and leading
%               and trailing  blanks in `bodynm' are not significant.
%               Optionally, you may supply the integer ID code for the
%               object as an integer string. For example both 'MOON'
%               and '301' are legitimate strings that indicate the
%               moon is the body of interest.
%
%      item     the item name to return.
%
%               [1,c2] = size(item); char = class(item)
%
%                  or
%
%               [1,1] = size(item); cell = class(item)
%
%               Together, the NAIF ID code of the body and the item name
%               combine to form a kernel variable name, e.g.,
%
%                    'BODY599_RADII'
%                    'BODY401_POLE_RA'
%
%               The values associated with the kernel variable having
%               the name constructed as shown are sought. Below
%               we'll take the shortcut of calling this kernel variable
%               the "requested kernel variable."
%
%               Note that `item' *is* case-sensitive. This attribute
%               is inherited from the case-sensitivity of kernel
%               variable names.
%
%      maxn     the maximum number of kernel pool values to returns.
%
%               [1,1] = size(maxn); int32 = class(maxn)
%
%   the call:
%
%      values = cspice_bodvrd(bodynm, item, maxn)
%
%   returns:
%
%      values   an array of at most `maxn' values associated with the
%               requested kernel variable.
%
%               [1,n] = size(values); double = class(values)
%               with n <= `maxn'.
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
%   1) Retrieve the radii of the Earth from the kernel pool, using both
%      'RADII' and 'radii' as the item name to return. Since the `item'
%      variable possesses case sensitivity, the later case should fail.
%      Trap the error and print it to the output.
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
%      function bodvrd_ex1()
%         %
%         % Load a PCK file.
%         %
%         cspice_furnsh( 'pck00008.tpc' );
%
%         %
%         % When the kernel variable
%         %
%         %    BODY399_RADII
%         %
%         % is present in the kernel pool---normally because a PCK
%         % defining this variable has been loaded (as is the case
%         % here)---the call
%         %
%         values1 = cspice_bodvrd( 'EARTH', 'RADII', 3);
%         fprintf( 'EARTH RADII: %10.3f  %10.3f  %10.3f\n', values1 )
%
%         %
%         % returns the dimension and values associated with the
%         % variable "BODY399_RADII".
%         %
%
%         %
%         % The call lacks case sensitivity in the `bodynm' variable.
%         %
%         values2 = cspice_bodvrd( 'earth', 'RADII', 3);
%         fprintf( 'earth RADII: %10.3f  %10.3f  %10.3f\n', values2 )
%
%         %
%         % The `item' variable possesses case sensitivity.
%         %
%         try
%
%            %
%            % A call with improper case in `item' will fail.
%            %
%            values3 = cspice_bodvrd( 'EARTH', 'radii', 3)
%
%         catch
%
%            %
%            % Catch the error, return the error string to the user.
%            %
%            disp( 'Expected error signaled:' )
%            disp( ' ' )
%            disp( lasterr )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      EARTH RADII:   6378.140    6378.140    6356.750
%      earth RADII:   6378.140    6378.140    6356.750
%      Expected error signaled:
%
%      mice: SPICE(KERNELVARNOTFOUND): [bodvrd_c->BODVRD] The variable
%      BODY399_radii could not be found in the kernel pool. (CSPICE_N0066)
%
%
%      Note that the SPICE(KERNELVARNOTFOUND) error is signaled
%      when the requested item is not found in the kernel pool.
%
%-Particulars
%
%   This routine simplifies looking up PCK kernel variables by
%   constructing names of requested kernel variables and by
%   performing error checking.
%
%   This routine is intended for use in cases where the maximum number of
%   values that may be returned is known at compile time. The caller fetches
%   all of the values associated with the specified kernel variable via a
%   single call to this routine. If the number of values to be fetched cannot
%   be known until run time, the lower-level routine cspice_gdpool should be
%   used instead. cspice_gdpool supports fetching arbitrary amounts of data
%   in multiple "chunks."
%
%   This routine is intended for use in cases where the requested
%   kernel variable is expected to be present in the kernel pool. If
%   the variable is not found or has the wrong data type, this
%   routine signals an error. In cases where it is appropriate to
%   indicate absence of an expected kernel variable by returning a
%   boolean "found flag" with the value false, again the routine
%   cspice_gdpool should be used.
%
%-Exceptions
%
%   1)  If the input body name cannot be translated to an ID code,
%       and if the name is not a string representation of an integer
%       (for example, '399'), the error SPICE(NOTRANSLATION) is
%       signaled by a routine in the call tree of this routine.
%
%   2)  If the requested kernel variable is not found in the kernel
%       pool, the error SPICE(KERNELVARNOTFOUND) is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If the requested kernel variable is found but the associated
%       values aren't numeric, the error SPICE(TYPEMISMATCH) is
%       signaled by a routine in the call tree of this routine.
%
%   4)  If the dimension of `values' indicated by `maxn' is too small to
%       contain the requested values, the error SPICE(ARRAYTOOSMALL)
%       is signaled by a routine in the call tree of this routine. The
%       output array `values' must be declared with sufficient size to
%       contain all of the values associated with the requested kernel
%       variable.
%
%   5)  If the input dimension `maxn' indicates there is more room in
%       `values' than there really is---for example, if `maxn' is 10 but
%       values is declared with dimension 5---and the dimension of the
%       requested kernel variable is larger than the actual dimension
%       of `values', then this routine may overwrite memory. The results
%       are unpredictable.
%
%   6)  If any of the input arguments, `bodynm', `item' or `maxn', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `bodynm', `item' or `maxn', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
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
%       Changed argument name "body" to "bodynm" to comply with NAIF
%       standard.
%
%       Fixed typo in Usage message.
%
%       Edited the header to comply with NAIF standard. Formatted
%       example's output.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 29-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   fetch constants for a body from the kernel pool
%   physical constants for a body
%
%-&

function [values] = cspice_bodvrd(bodynm, item, maxn)

   switch nargin
      case 3

         bodynm = zzmice_str(bodynm);
         item   = zzmice_str(item);
         maxn   = zzmice_int(maxn);

      otherwise

         error ( ['Usage:  [values] = cspice_bodvrd( `bodynm`, ', ...
                                                    '`item`, ',   ...
                                                    'maxn )' ] )

   end

   try
      [values] = mice( 'bodvrd_c', bodynm, item, maxn);
   catch spiceerr
      rethrow(spiceerr)
   end
