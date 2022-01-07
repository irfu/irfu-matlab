%-Abstract
%
%   CSPICE_BODC2S translates a body ID code to either the corresponding name
%   or if no name to ID code mapping exists, the string representation of
%   the body ID value.
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
%      code     SPICE code(s) for a set of bodies: planets, satellites,
%               barycenters, DSN stations, spacecraft, asteroids, comets,
%               or other ephemeris object.
%
%               [1,n] = size(code); int32 = class(code)
%
%   the call:
%
%      [name] = cspice_bodc2s( code )
%
%   returns:
%
%      name     the name(s) corresponding to `code' if a mapping between
%               `code' and a body name exists within SPICE, assigned either
%               in SPICE or by the user.
%
%               [n,c1] = size(name); char = class(name)
%
%               If `code' has more than one translation, then the most
%               recently defined `name' corresponding to `code' is returned.
%               `name' will have the exact format (case and blanks) as when
%               the name/code pair was defined.
%
%               If the input value of `code' does not map to a body name,
%               `name' returns with the string representation of `code'.
%
%               `name' returns with the same vectorization measure (N) as
%               `code'.
%
%-Parameters
%
%   MAXL        is the maximum allowable length of a body name. The
%               current value of this parameter is 36.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Apply the cspice_bodc2s call to several IDs representing codes
%      included in the default SPICE ID-name lists and codes not
%      included in the list.
%
%      Example code begins here.
%
%
%      function bodc2s_ex1()
%         %
%         % Assign an array of body ID codes. Not all the listed codes
%         % map to a body name.
%         %
%         code = [ 399, 0, 3, -77, 11, -1, 6000001 ];
%
%         name = cspice_bodc2s( code );
%
%
%         %
%         % Loop over the `code' array, call cspice_bodc2s for each
%         % element of `code'.
%         %
%         fprintf( 'Code      Name\n' )
%         fprintf( '-------   -----------------------\n' )
%
%         for i=1:7
%            fprintf( '%7d   %s\n', code(i), name(i,:) )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Code      Name
%      -------   -----------------------
%          399   EARTH
%            0   SOLAR SYSTEM BARYCENTER
%            3   EARTH BARYCENTER
%          -77   GALILEO ORBITER
%           11   11
%           -1   GEOTAIL
%      6000001   6000001
%
%
%      Note that the codes 11 and 6000001 did not map to a name so the
%      call returns as `name' the string expression of the codes.
%
%-Particulars
%
%   cspice_bodc2s is one of five related functions,
%
%      cspice_bods2c      Body string to code
%      cspice_bodc2s      Body code to string
%      cspice_bodn2c      Body name to code
%      cspice_bodc2n      Body code to name
%      cspice_boddef      Body name/code definition
%
%   cspice_bods2c, cspice_bodc2s, cspice_bodn2c, and cspice_bodc2n
%   perform translations between body names and their corresponding
%   integer ID codes which are used in SPICE files and routines.
%
%   cspice_bods2c is a slightly more general version of cspice_bodn2c:
%   support for strings containing ID codes in string format enables a caller
%   to identify a body using a string, even when no name is associated with
%   that body.
%
%   cspice_bodc2s is a general version of cspice_bodc2n; the routine returns
%   either the name assigned in the body ID to name mapping or a string
%   representation of the `code' value if no mapping exists.
%
%   cspice_boddef assigns a body name to ID mapping. The mapping has
%   priority in name-to-ID and ID-to-name translations.
%
%   Refer to naif_ids.req for the list of name/code associations built
%   into SPICE, and for details concerning adding new name/code
%   associations at run time by loading text kernels.
%
%-Exceptions
%
%   1)  If there is any problem with the body name-ID mapping kernel
%       variables present in the kernel pool, an error is signaled by
%       a routine in the call tree of this routine.
%
%   2)  If the input argument `code' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `code' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Body-name mappings may be defined at run time by loading text
%   kernels containing kernel variable assignments of the form
%
%      NAIF_BODY_NAME += ( <name 1>, ... )
%      NAIF_BODY_CODE += ( <code 1>, ... )
%
%   See naif_ids.req for details.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
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
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Corrected typos in the
%       header.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 28-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 01-JUN-2009 (EDW)
%
%-Index_Entries
%
%   body ID code to string
%
%-&

function [name] = cspice_bodc2s(code)

   switch nargin
      case 1

         code = zzmice_int(code);

      otherwise

         error ( 'Usage: [_`name`_] = cspice_bodc2s(_code_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [bodc2s] = mice('bodc2s_s', code);
      name     = char( bodc2s.name );
   catch spiceerr
      rethrow(spiceerr)
   end

