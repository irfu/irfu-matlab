%-Abstract
%
%   CSPICE_SRFS2C translates a surface string, together with a body string,
%   to the corresponding surface ID code. The input strings may contain
%   names or integer ID codes.
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
%      srfstr   a string designating a surface.
%
%               [1,c1] = size(srfstr); char = class(srfstr)
%
%                  or
%
%               [1,1] = size(srfstr); cell = class(srfstr)
%
%               `srfstr' may contain a surface name or a string
%               representation of the surface's integer ID code.
%
%               If, for the body specified by `bodstr', multiple surface
%               names are associated with one surface ID code, then any
%               of these names may be used as the value of `srfstr'.
%
%               Case and leading and trailing blanks in a surface name
%               are not significant. Sequences of consecutive embedded
%               blanks are considered equivalent to a single blank.
%               For example, all of the strings below are considered
%               to be equivalent:
%
%                  'MGS MOLA 128 PIXEL/DEG'
%                  'MGS MOLA 128 pixel/deg'
%                  'MGS MOLA 128 PIXEL/DEG   '
%                  'MGS MOLA 128    PIXEL/DEG'
%                  '   MGS MOLA 128 PIXEL/DEG'
%
%               However,
%
%                  'MGSMOLA 128PIXEL/DEG'
%
%               is not equivalent to the names above.
%
%
%      bodstr   a string designating the body associated with the input
%               surface string.
%
%               [1,c2] = size(bodstr); char = class(bodstr)
%
%                  or
%
%               [1,1] = size(bodstr); cell = class(bodstr)
%
%               `bodstr' may contain a body name or a string representation
%               of the body's integer ID code. For example, `bodstr' may
%               contain
%
%                  '1000012'
%
%               instead of
%
%                  '67P/CHURYUMOV-GERASIMENKO (1969 R1)'
%
%               Case and leading and trailing blanks in a name are not
%               significant. The treatment of blanks in `bodstr' is the
%               same as for `srfstr'.
%
%   the call:
%
%      [code, found] = cspice_srfs2c( srfstr, bodstr )
%
%   returns:
%
%      code     integer ID code of the surface designated by `srfstr', for
%               the body designated by `bodstr', if for this body an
%               association exists between the input surface string and a
%               surface ID code.
%
%               [1,1] = size(code); int32 = class(code)
%
%               `code' is defined if and only if the output flag `found' is
%               true
%
%      found    a logical flag that is true if a surface code corresponding
%               to the input strings was found and false otherwise.
%
%               [1,1] = size(found); logical = class(found)
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
%   1) Supposed a text kernel has been loaded that contains
%      the following assignments:
%
%         NAIF_SURFACE_NAME += ( 'MGS MOLA  64 pixel/deg',
%                                'MGS MOLA 128 pixel/deg',
%                                'PHOBOS GASKELL Q512'     )
%         NAIF_SURFACE_CODE += (   1,   2,    1 )
%         NAIF_SURFACE_BODY += ( 499, 499,  401 )
%
%      Translate each surface and body ID code pair to the
%      associated surface name. Also perform a translation
%      for a surface ID having no matching name.
%
%      Use the meta-kernel shown below to define the required SPICE
%      kernel variables.
%
%
%         KPL/MK
%
%         File: srfs2c_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The file contents shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%
%         \begindata
%
%         NAIF_SURFACE_NAME += ( 'MGS MOLA  64 pixel/deg',
%                                'MGS MOLA 128 pixel/deg',
%                                'PHOBOS GASKELL Q512'     )
%         NAIF_SURFACE_CODE += (   1,   2,    1 )
%         NAIF_SURFACE_BODY += ( 499, 499,  401 )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function srfs2c_ex1()
%
%         bodstr = { 'MARS', 'PHOBOS', 'MARS', ...
%                    '499',  'PHOBOS', '499',  ...
%                     'MARS', 'ZZZ' };
%         srfstr = { 'MGS MOLA  64 pixel/deg', ...
%                    'PHOBOS GASKELL Q512',    ...
%                    'MGS MOLA 128 pixel/deg', ...
%                    'MGS MOLA  64 pixel/deg', ...
%                    '1',   ...
%                    '2',   ...
%                    'ZZZ', ...
%                    '1' };
%
%         tf     = { 'false', 'true' };
%         meta   = 'srfs2c_ex1.tm';
%
%         cspice_furnsh( meta );
%
%         for i= 1:numel( srfstr )
%
%            [ surfid, found] =  cspice_srfs2c ( srfstr(i), bodstr(i) );
%
%            fprintf([ 'surface string   = %s\n'   ...
%                      'body string      = %s\n'   ...
%                      'surface ID found = %s\n'], ...
%                      char(srfstr(i)), ...
%                      char(bodstr(i)), ...
%                      char( tf(int32(found)+1) )  )
%
%         if ( found )
%             fprintf( 'surface ID       = %d\n', surfid );
%         end
%
%         fprintf( '\n' )
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
%      surface string   = MGS MOLA  64 pixel/deg
%      body string      = MARS
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = PHOBOS GASKELL Q512
%      body string      = PHOBOS
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = MGS MOLA 128 pixel/deg
%      body string      = MARS
%      surface ID found = true
%      surface ID       = 2
%
%      surface string   = MGS MOLA  64 pixel/deg
%      body string      = 499
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = 1
%      body string      = PHOBOS
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = 2
%      body string      = 499
%      surface ID found = true
%      surface ID       = 2
%
%      surface string   = ZZZ
%      body string      = MARS
%      surface ID found = false
%
%      surface string   = 1
%      body string      = ZZZ
%      surface ID found = false
%
%
%-Particulars
%
%   Surfaces are always associated with bodies (which usually are
%   ephemeris objects). For any given body, a mapping between surface
%   names and surface ID codes can be established.
%
%   Bodies serve to disambiguate surface names and ID codes: the set
%   of surface names and surface ID codes for a given body can be
%   thought of as belonging to a name space. A given surface ID code
%   or surface name may be used for surfaces of multiple bodies,
%   without conflict.
%
%   Associations between surface names and ID codes are always made
%   via kernel pool assignments; there are no built-in associations.
%
%   cspice_srfs2c is one of four related subroutines:
%
%      cspice_srfs2c    Surface string and body string to surface ID code
%      cspice_srfscc    Surface string and body ID code to surface ID code
%      cspice_srfc2s    Surface ID code and body ID code to surface string
%      cspice_srfcss    Surface ID code and body string to surface string
%
%   cspice_srfs2c, cspice_srfc2s, cspice_srfscc, and cspice_srfcss perform
%   translations between surface strings and their corresponding integer
%   ID codes.
%
%   Refer to naif_ids.req for details concerning adding new surface
%   name/code associations at run time by loading text kernels.
%
%-Exceptions
%
%   1)  If the input surface string does not map to an ID code
%       and does not represent an integer, the output `code' is
%       undefined and the output `found' is set to false.
%
%       This case is not treated as an error.
%
%   2)  If the input body string does not map to an ID code and does
%       not represent an integer, the output `code' is undefined and
%       the output `found' is set to false.
%
%       This case is not treated as an error.
%
%   3)  If any of the input arguments, `srfstr' or `bodstr', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `srfstr' or `bodstr', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   Surface name-to-ID mappings may be defined at run time by loading
%   text kernels containing kernel variable assignments of the form
%
%      NAIF_SURFACE_NAME += ( <surface name 1>, ... )
%      NAIF_SURFACE_CODE += ( <surface code 1>, ... )
%      NAIF_SURFACE_BODY += ( <body code 1>,    ... )
%
%   Above, the ith elements of the lists on the assignments' right
%   hand sides together define the ith surface name/ID mapping.
%
%   The same effect can be achieved using assignments formatted as
%   follows:
%
%      NAIF_SURFACE_NAME += <surface name 1>
%      NAIF_SURFACE_CODE += <surface code 1>
%      NAIF_SURFACE_BODY += <body code 1>
%
%      NAIF_SURFACE_NAME += <surface name 2>
%      NAIF_SURFACE_CODE += <surface code 2>
%      NAIF_SURFACE_BODY += <body code 2>
%
%         ...
%
%   Note the use of the
%
%      +=
%
%   operator; this operator appends to rather than overwrites the
%   kernel variable named on the left hand side of the assignment.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   DSK.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
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
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR) (NJB)
%
%       Updated description of "srfstr" to indicate that any
%       surface name alias may be used.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections. Edited
%       the header to comply with NAIF standard.
%
%       Added call to cspice_kclear in code example.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 01-MAR-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   surface string and body string to surface ID code
%
%-&

function [code,found] = cspice_srfs2c( srfstr, bodstr)

   switch nargin
      case 2

         srfstr = zzmice_str(srfstr);
         bodstr = zzmice_str(bodstr);

      otherwise

         error ( 'Usage: [ code, found] = cspice_srfs2c( `srfstr`, `bodstr`)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [srfs2c] = mice('srfs2c_s', srfstr, bodstr);
      code     = reshape( [srfs2c.code],  1, [] );
      found    = reshape( [srfs2c.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end
