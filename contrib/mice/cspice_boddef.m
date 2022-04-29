%-Abstract
%
%   CSPICE_BODDEF defines a body name/ID code pair for later translation via
%   cspice_bodn2c or cspice_bodc2n.
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
%      name     the string defining the name to associate to the ID value
%               `code'.
%
%               [1,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,1] = size(name); cell = class(name)
%
%               The case and positions of blanks in a name are significant.
%               cspice_bodc2n returns the same string (case and space) most
%               recently mapped to a code. When `name' consists of more than
%               one word, the words require separation by at least one blank.
%
%               The kernel sub-system stores `name' as described in the
%               cspice_boddef call, but creates an equivalence class based on
%               `name' for comparisons in cspice_bodn2c. This class ignores
%               leading and trailing whitespace, compresses interior
%               whitespace to a single space, and ignores character case.
%
%               The following strings belong to the same equivalence
%               class:
%
%                  'JUPITER BARYCENTER'
%                  'Jupiter Barycenter'
%                  'JUPITER BARYCENTER   '
%                  'JUPITER    BARYCENTER'
%                  '   JUPITER BARYCENTER'
%
%               However, 'JUPITERBARYCENTER' is distinct from the names above.
%
%      code     the NAIF ID code corresponding to `name'.
%
%               [1,1] = size(code); int32 = class(code)
%
%   the call:
%
%      cspice_boddef( name, code )
%
%   performs the mapping assignment
%
%         `name' -> `code'
%
%      and
%
%         `code' -> 'name'
%
%      The `code' -> `name' assignment supersedes any other mapping for
%      `code'.
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
%   1) Map a non-existent code and name to each other, and test that
%      the translation is found when searching by code or by name.
%
%      Example code begins here.
%
%
%      function boddef_ex1()
%         %
%         % Map a non-existent code and name to each other.
%         %
%         cspice_boddef( 'spud',  -69 );
%
%         %
%         % Retrieve the code for name 'spud'.
%         %
%         [ code, found ] = cspice_bodn2c( 'spud' );
%
%         %
%         % Check we found a mapping.
%         %
%         if ( found )
%            txt = sprintf( 'ID for spud : %i', code );
%         else
%            txt = 'Found no mapping for spud.';
%         end
%
%         disp( txt )
%
%         %
%         % Retrieve the name for ID -69.
%         %
%         [ name, found ] = cspice_bodc2n( -69 );
%
%         %
%         % Check we found a mapping.
%         %
%         if (found)
%            txt = sprintf( 'Name for -69: %s', name );
%         else
%            txt = 'Found no mapping for -69.';
%         end
%
%         disp( txt )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      ID for spud : -69
%      Name for -69: spud
%
%
%-Particulars
%
%   cspice_boddef is one of five related subroutines,
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
%   1)  If improper inputs are supplied, or if there is insufficient
%       room to store the requested addition, an error is signaled by
%       a routine in the call tree of this routine.
%
%   2)  If the length of `name' exceeds the maximum allowed length for a
%       body name, the name stored in the kernel pool will be
%       truncated on the right.
%
%   3)  If a name-code definition inserted into this routine seems to
%       have no effect, it is possible that the contents of the
%       definition are masked by the higher precedence kernel pool
%       assignments. See the -Particulars section of this document
%       for more information.
%
%   4)  If any of the input arguments, `name' or `code', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `name' or `code', is not of the
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
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
%   -Mice Version 1.0.2, 28-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 16-MAY-2009 (EDW)
%
%       Edit to -Particulars section to document the cspice_bodc2s routine.
%       Extended argument descriptions in the -I/O section.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   body name/id code definition
%
%-&

function cspice_boddef(name, code)

   switch nargin
      case 2

         name = zzmice_str(name);
         code = zzmice_int(code);

      otherwise

         error ( 'Usage: cspice_boddef(`name`, code)' )

   end

   try
      mice('boddef_c', name, code);
   catch spiceerr
      rethrow(spiceerr)
   end


