%-Abstract
%
%   CSPICE_BODN2C translates the name of a body or object to the corresponding
%   SPICE integer ID code.
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
%      name     name(s) of a body or object,  such as a planet, satellite,
%               comet, asteroid, barycenter, DSN station, spacecraft, or
%               instrument, "known" to the SPICE system, whether through
%               hard-coded registration or run-time registration in the
%               SPICE kernel pool.
%
%               [n,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,n] = size(bodynm); cell = class(bodynm)
%
%               Case and leading and trailing blanks in a name are not
%               significant. However when a name is made up of more than one
%               word, they must be separated by at least one blank. That is,
%               all of the following strings are equivalent names:
%
%                  'JUPITER BARYCENTER'
%                  'Jupiter Barycenter'
%                  'JUPITER BARYCENTER   '
%                  'JUPITER    BARYCENTER'
%                  '   JUPITER BARYCENTER'
%
%               However, 'JUPITERBARYCENTER' is not equivalent to the names
%               above.
%
%   the call:
%
%      [code, found] = cspice_bodn2c( name )
%
%   returns:
%
%      code     containing the SPICE code(s) assigned either by SPICE or the
%               user to `name'.
%
%               [1,n] = size(code); int32 = class(code)
%
%      found    flag(s) indicating if the kernel subsystem translated `name'
%               to a corresponding `code'.
%
%               [1,n] = size(found); logical = class(found)
%
%               `found' and `code' return with the same vectorization
%               measure, N, as `name'.
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
%   1) Apply the cspice_bodn2c call to several body names to retrieve
%      their associated NAIF IDs included in the default SPICE ID-name
%      lists and a name not included in that list.
%
%      Example code begins here.
%
%
%      function bodn2c_ex1()
%         %
%         % Retrieve the NAIF ID associated to a body name.
%         %
%         disp( 'Scalar:' )
%         name            = 'Hyperion';
%         [code, found] = cspice_bodn2c( name );
%
%         %
%         % Output the mapping if it exists.
%         %
%         if ( found )
%            txt = sprintf( 'Body code %i maps to name %s',  code, name );
%            disp(txt)
%         end
%
%         disp(' ')
%
%         %
%         % Create an array of body names. Include one unknown name.
%         %
%         disp( 'Vector:' )
%         name           = strvcat( 'Triton',   'Mimas', 'Oberon',         ...
%                                   'Callisto', 'Halo' );
%         [code, found] = cspice_bodn2c( name );
%
%         n_elements = size(code,2);
%
%         %
%         % Loop over the output array.
%         %
%         for n=1:n_elements
%
%            %
%            % Check for a valid name/code mapping.
%            %
%            if( found(n) )
%               txt = sprintf( 'Body code %i maps to name %s', code(n),    ...
%                                                              name(n,:) );
%               disp(txt)
%            else
%               txt = sprintf( 'Unknown body name %s', name(n,:) );
%               disp(txt)
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Scalar:
%      Body code 607 maps to name Hyperion
%
%      Vector:
%      Body code 801 maps to name Triton
%      Body code 601 maps to name Mimas
%      Body code 704 maps to name Oberon
%      Body code 504 maps to name Callisto
%      Unknown body name Halo
%
%
%-Particulars
%
%   A sister version of this routine exists named mice_bodn2c that returns
%   the output arguments as fields in a single structure.
%
%   cspice_bodn2c is one of five related subroutines,
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
%   Programmers writing user interface code should consider using the
%   Mice routine cspice_bods2c. cspice_bods2c provides more flexibility
%   in handling input strings, since it accepts both body names and
%   strings representing integer ID codes, for example '399'.
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
%   2)  Body name strings are upper-cased, their leading and trailing
%       blanks removed, and embedded blanks are compressed out, after
%       which they get truncated to the maximum body name length MAXL.
%       Therefore, two body names that differ only after that maximum
%       length are considered equal.
%
%   3)  If the input argument `name' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `name' is not of the expected type, or
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
%   1)  See exception <2>.
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
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 01-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.2, 12-MAR-2012 (EDW) (SCK)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%       Corrected minor typo in header.
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
%   body name to code
%
%-&

function [code,found] = cspice_bodn2c(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [_code_, _found_] = cspice_bodn2c(_`name`_)' )

   end


   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [bodn2c] = mice('bodn2c_s',name);
      code     = reshape( [bodn2c.code],  1, [] );
      found    = reshape( [bodn2c.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end


