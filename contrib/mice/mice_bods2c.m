%-Abstract
%
%   MICE_BODS2C translates a string containing a body name or
%   ID code to the corresponding integer code.
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
%      name   containing the name or ID code of a body or  object, such as a
%             planet, satellite, comet, asteroid, barycenter, DSN station,
%             spacecraft, or instrument.
%
%             [n,m] = size(name); char = class(name)
%
%                  or
%
%             [1,n] = size(name); cell = class(name)
%
%             If 'name' contains the name of a body or object, that name must be
%             "known" to the SPICE system, whether through hard-coded
%             registration or run-time registration in the SPICE kernel pool.
%
%             Case and leading and trailing blanks in `name' are not
%             significant. However when a name is made up of more than one
%             word, they must be separated by at least one blank.  That is, all
%             of the following strings are equivalent names:
%
%                     'JUPITER BARYCENTER'
%                     'Jupiter Barycenter'
%                     'JUPITER BARYCENTER   '
%                     'JUPITER    BARYCENTER'
%                     '   JUPITER BARYCENTER'
%
%             However, 'JUPITERBARYCENTER' is not equivalent to the names above.
%
%             If 'name' is a string representation of an integer, for example
%
%                '399'
%
%             the string will be translated to the equivalent integer datum.
%             The input integer need not be one recognized by the SPICE system:
%             the integer need not be a built-in NAIF ID code, nor need it be
%             associated with a name via run-time registration.
%
%   the call:
%
%      ID = mice_bods2c( name )
%
%   returns:
%
%      ID   the structure(s) associating a body name with a corresponding 
%           SPICE ID. 
%
%           [1,n] = size(ID); struct = class(ID)
%
%           Each structure consists of the fields:
%
%              name   the "name" of a particular body. If a mapping
%                     does not exist, the 'name' field returns as NULL.
%
%                     [1,c1] = size(ID(i).name); char = class(ID(i).name)
%
%              code   the SPICE code assigned either
%                     by SPICE or the user to 'name'. If a mapping
%                     does not exist, the 'code' field returns as 0.
%
%                     [1,1] = size(ID(i).code); int32 = class(ID(i).code)
%
%      'ID' returns with the same vectorization measure, N, as 'name'.
%
%-Examples
%
%      %
%      % Retrieve the NAIF ID associated to a body name.
%      %
%      disp( 'Scalar:' )
%      name = 'Hyperion';
%      ID   = mice_bods2c( name );
%
%      %
%      % Output the mapping if it exists.
%      %
%      if ( ID.found )
%         txt = sprintf( 'String %s maps to ID %i', ...
%                         ID.name, ID.code );
%         disp(txt)
%      end
%
%      disp(' ')
%
%      %
%      % Create an array of strings. Include one string not an integer
%      % and unknown to the SPICE system.
%      %
%      disp( 'Vector:' )
%      name = strvcat( 'Cassini'   , '399',  ...
%                      'Planet Bob', 'MARS', ...
%                      '-123456'   , '987654' );
%      ID   = mice_bods2c( name );
%
%      n_elements = size(ID);
%
%      %
%      % Loop over the output array.
%      %
%      for i=1:n_elements(1)
%
%         %
%         % Check for a valid name/ID mapping.
%         %
%         if( ID(i).found ) )
%            txt = sprintf( 'String %s maps to ID %i', ...
%                            ID(i).name, ID(i).code );
%            disp(txt)
%         else
%            txt = sprintf( 'Unknown string ID %s', name(i,:) );
%            disp(txt)
%         end
%
%   MATLAB outputs:
%
%      Scalar:
%      String Hyperion maps to ID 607
%
%      Vector:
%      String Cassini    maps to ID -82
%      String 399        maps to ID 399
%      Unknown string ID Planet Bob
%      String MARS       maps to ID 499
%      String -123456    maps to ID -123456
%      String 987654     maps to ID 987654
%
%-Particulars
%
%   A sister version of this routine exists named cspice_bods2c that returns
%   the structure field data as separate arguments.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine bods2c_c.
%
%   MICE.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 01-DEC-2014, EDW (JPL)
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   body name to code
%
%-&

function [code] = mice_bods2c(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [_code_] = mice_bods2c(_`name`_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [code] = mice('bods2c_s',name);
   catch
      rethrow(lasterror)
   end


