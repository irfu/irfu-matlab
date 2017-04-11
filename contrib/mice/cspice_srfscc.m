%-Abstract
%
%   CSPICE_SRFSCC translates a surface string, together with a body ID code,
%   to the corresponding surface ID code. The input surface string may 
%   contain a name or an integer ID code. 
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
%      srfstr     is a string designating a surface. `srfstr' may contain
%                 a surface name or a string representation of the
%                 surface's integer ID code.
%
%                 [1,c1] = size(srfstr), char = class(srfstr)
%
%                    or
%
%                 [1,1] = size(srfstr), cell = class(srfstr)
%
%                 Case and leading and trailing blanks in a surface name
%                 are not significant. Sequences of consecutive embedded
%                 blanks are considered equivalent to a single blank. For
%                 example, all of the strings below are considered to be
%                 equivalent:
%
%                    'MGS MOLA 128 PIXEL/DEG'
%                    'MGS MOLA 128 pixel/deg'
%                    'MGS MOLA 128 PIXEL/DEG   '
%                    'MGS MOLA 128    PIXEL/DEG'
%                    '   MGS MOLA 128 PIXEL/DEG'
%
%                 However,
%
%                    'MGSMOLA 128PIXEL/DEG'
%
%                 is not equivalent to the names above.
%
%
%      bodyid     is the integer ID code of the body associated with the
%                 surface designated by `srfstr'.
%
%                 [1,1] = size(bodyid); int32 = class(bodyid)
%
%   the call:
%
%      [ code, found] = cspice_srfscc( surfce, bodyid )
%
%   returns:
%
%      code       is integer ID code of the surface designated by `srfstr',
%                 for the body designated by `bodyid', if for this body an
%                 association exists between the input surface string and a
%                 surface ID code. `code' is defined if and only if the
%                 output flag `found' is true.
%
%                 [1,1] = size(code); int32 = class(code)
%
%      found      is a logical flag that is true if a surface code
%                 corresponding to the input surface string and body ID
%                 code was found. `found' is false otherwise.
%
%                 [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%   The formatting of the results shown for this example may differ
%   across platforms.
%
%      Supposed a text kernel has been loaded that contains
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
%         File: srfc2s_t1.tm
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
%   Example(1):
%
%      function srfscc_t
%
%         bodyid = [ 499, 401, 499, 499, 401, 499, 499 ];
%         srfstr = { 'MGS MOLA  64 pixel/deg', ...
%                       'PHOBOS GASKELL Q512',    ...
%                       'MGS MOLA 128 pixel/deg', ...
%                       'MGS MOLA  64 pixel/deg', ...
%                       '1',   ...
%                       '2',   ...
%                       'ZZZ', ...
%                       '1' };
%
%         tf     = { 'false', 'true' };
%         meta   = 'srfc2s_t1.tm';
%
%         cspice_furnsh( meta );
%
%         for i= 1:numel( bodyid )
%
%            [ surfid, found] = cspice_srfscc( srfstr(i), bodyid(i));
%
%            fprintf([ 'surface string   = %s\n'   ...
%                         'body ID          = %d\n'   ...
%                         'surface ID found = %s\n'], ...
%                         char(srfstr(i)), ...
%                         bodyid(i),       ...
%                         char( tf(int32(found)+1) )  )
%
%         if ( found )
%             fprintf( 'surface ID       = %d\n', surfid );
%         end
%
%         fprintf( '\n' )
%
%         end
%
%   MATLAB outputs:
%
%      surface string   = MGS MOLA  64 pixel/deg
%      body ID          = 499
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = PHOBOS GASKELL Q512
%      body ID          = 401
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = MGS MOLA 128 pixel/deg
%      body ID          = 499
%      surface ID found = true
%      surface ID       = 2
%
%      surface string   = MGS MOLA  64 pixel/deg
%      body ID          = 499
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = 1
%      body ID          = 401
%      surface ID found = true
%      surface ID       = 1
%
%      surface string   = 2
%      body ID          = 499
%      surface ID found = true
%      surface ID       = 2
%
%      surface string   = ZZZ
%      body ID          = 499
%      surface ID found = false
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
%   cspice_srfscc is one of four related subroutines:
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
%-Required Reading
%
%   the CSPICE routine srfscc_c.
%   For important details concerning this module's function, please refer to
%
%   MICE.REQ
%   DSK.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 01-MAR-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   surface string and body ID code to surface ID code
%
%-&

function [ code, found] = cspice_srfscc( surfce, bodyid )

   switch nargin
      case 2

         surfce = zzmice_str(surfce);
         bodyid = zzmice_int(bodyid);

      otherwise

         error ( 'Usage: [ code, found] = cspice_srfscc( `surfce`, bodyid )' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [srfscc] = mice('srfscc_s', surfce, bodyid);
      code   = reshape( [srfscc.code],  1, [] );
      found    = reshape( [srfscc.found], 1, [] );
   catch
      rethrow(lasterror)
   end




