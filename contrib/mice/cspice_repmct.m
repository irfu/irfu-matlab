%-Abstract
%
%   CSPICE_REPMCT replaces a marker with the text representation of a
%   cardinal number.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      in       an arbitrary character string.
%
%               [1,c1] = size(in); char = class(in)
%
%                  or
%
%               [1,1] = size(in); cell = class(in)
%
%      marker   an arbitrary character string.
%
%               [1,c2] = size(marker); char = class(marker)
%
%                  or
%
%               [1,1] = size(marker); cell = class(marker)
%
%               The first occurrence of `marker' in the input string is to
%               be replaced by the text representation of the cardinal number
%               `value'.
%
%               Leading and trailing blanks in `marker' are NOT
%               significant. In particular, no substitution is performed
%               if `marker' is blank or empty.
%
%      value    an arbitrary integer.
%
%               [1,1] = size(value); int32 = class(value)
%
%      rtcase   indicates the case of the replacement text.
%
%               [1,1] = size(rtcase); char = class(rtcase)
%
%                  or
%
%               [1,1] = size(rtcase); cell = class(rtcase)
%
%               `rtcase' may be any of the following:
%
%                  rtcase   Meaning        Example
%                  ------   -----------    -----------------------
%                  U, u     Uppercase      ONE HUNDRED FIFTY-THREE
%
%                  L, l     Lowercase      one hundred fifty-three
%
%                  C, c     Capitalized    One hundred fifty-three
%
%   the call:
%
%      [out] = cspice_repmct( in, marker, value, rtcase )
%
%   returns:
%
%      out      the string obtained by substituting the text representation
%               of the cardinal number `value' for the first occurrence of
%               `marker' in the input string.
%
%               [1,c3] = size(out); char = class(out)
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
%   1) The following example illustrate the use of cspice_repmct to
%      replace a marker within a string with the cardinal text
%      corresponding to an integer.
%
%
%      Example code begins here.
%
%
%      function repmct_ex1()
%
%         %
%         % 1. Uppercase
%         %
%         marker   = '#';
%         instr    = 'INVALID COMMAND. WORD # NOT RECOGNIZED.';
%
%         [outstr] = cspice_repmct( instr, marker, 5, 'U' );
%
%         fprintf( 'Case 1: Replacement text in uppercase.\n' )
%         fprintf( '   Input : %s\n', instr )
%         fprintf( '   Output: %s\n', outstr )
%         fprintf( '\n' )
%
%         %
%         % 2. Lowercase
%         %
%         marker   = ' XX ';
%         instr    = 'Word XX of the XX sentence was ...';
%
%         [outstr] = cspice_repmct( instr, marker, 5, 'L' );
%
%         fprintf( 'Case 2: Replacement text in lowercase.\n' )
%         fprintf( '   Input : %s\n', instr )
%         fprintf( '   Output: %s\n', outstr )
%         fprintf( '\n' )
%
%         %
%         % 2. Capitalized
%         %
%         marker   = ' XX ';
%         instr    = 'Name:  YY.  Rank:  XX.';
%
%         [outstr] = cspice_repmc( instr, 'YY', 'Moriarty' );
%         [outstr] = cspice_repmct( outstr, marker, 5, 'C' );
%
%         fprintf( 'Case 3: Replacement text capitalized.\n' )
%         fprintf( '   Input : %s\n', instr )
%         fprintf( '   Output: %s\n', outstr )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Case 1: Replacement text in uppercase.
%         Input : INVALID COMMAND. WORD # NOT RECOGNIZED.
%         Output: INVALID COMMAND. WORD FIVE NOT RECOGNIZED.
%
%      Case 2: Replacement text in lowercase.
%         Input : Word XX of the XX sentence was ...
%         Output: Word five of the XX sentence was ...
%
%      Case 3: Replacement text capitalized.
%         Input : Name:  YY.  Rank:  XX.
%         Output: Name:  Moriarty.  Rank:  Five.
%
%
%-Particulars
%
%   This is one of a family of related routines for inserting values
%   into strings. They are typically used to construct messages that
%   are partly fixed, and partly determined at run time. For example,
%   a message like
%
%      'Fifty-one pictures were found in directory [USER.DATA].'
%
%   might be constructed from the fixed string
%
%      '#1 pictures were found in directory #2.'
%
%   by the calls
%
%      [string] = cspice_repmct( string, '#1', 51, 'C' );
%      [string] = cspice_repmc( string, '#2', '[USER.DATA]' );
%
%   which substitute the cardinal text 'Fifty-one' and the character
%   string '[USER.DATA]' for the markers '#1' and '#2' respectively.
%
%   The complete list of routines is shown below.
%
%      cspice_repmc    ( Replace marker with character string value )
%      cspice_repmd    ( Replace marker with double precision value )
%      cspice_repmf    ( Replace marker with formatted d.p. value   )
%      cspice_repmi    ( Replace marker with integer value          )
%      cspice_repml    ( Replace marker with logical value          )
%      cspice_repmct   ( Replace marker with cardinal text          )
%      cspice_repmot   ( Replace marker with ordinal text           )
%
%-Exceptions
%
%   1)  If `marker' is blank or empty, or if `marker' is not a substring of
%       `in', no substitution is performed. (`out' and `in' are identical.)
%
%   2)  If the value of `rtcase' is not recognized, the error
%       SPICE(INVALIDCASE) is signaled by a routine in the call tree
%       of this routine. `out' is not changed.
%
%   3)  If any of the input arguments, `in', `marker', `value' or
%       `rtcase', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   4)  If any of the input arguments, `in', `marker', `value' or
%       `rtcase', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  `value' must be in the range accepted by the SPICELIB routine
%       INTTXT. This range is currently
%
%          ( -10**12, 10**12 )
%
%       Note that the endpoints of the interval are excluded.
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
%
%-Version
%
%   -Mice Version 1.0.0, 23-JAN-2021 (JDR)
%
%-Index_Entries
%
%   replace marker with cardinal text
%
%-&
function [out] = cspice_repmct( in, marker, value, rtcase )

   switch nargin
      case 4

         in     = zzmice_str(in, true);
         marker = zzmice_str(marker, true);
         value  = zzmice_int(value);
         rtcase = zzmice_str(rtcase);

      otherwise

         error ( [ 'Usage: [`out`] = '                                      ...
                   'cspice_repmct( `in`, `marker`, value, `rtcase` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [out] = mice('repmct_c', in, marker, value, rtcase);
   catch spiceerr
      rethrow(spiceerr)
   end
