%-Abstract
%
%   CSPICE_NEXTWD returns the next word in a given character string, and
%   left justifies the rest of the string.
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
%      string   the input string to be parsed.
%
%               [1,c1] = size(string); char = class(string)
%
%                  or
%
%               [1,1] = size(string); cell = class(string)
%
%               Each word of this string is a maximal sequence of
%               consecutive non-blank characters.
%
%   the call:
%
%      [next, rest] = cspice_nextwd( string )
%
%   returns:
%
%      next     the first word in `string'.
%
%               [1,c2] = size(next); char = class(next)
%
%               It is called the "next" word because cspice_nextwd is
%               typically called repeatedly to find the words of the input
%               string in left-to-right order. A word is a maximal sequence
%               of consecutive non-blank characters. `next' is always
%               returned left-justified.
%
%               If `string' is blank or empty, `next' is empty.
%
%      rest     the remaining part of `string', left-justified after the
%               removal of `next'.
%
%               [1,c3] = size(rest); char = class(rest)
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
%   1) Given a character string, get the sequence of words within.
%
%
%      Example code begins here.
%
%
%      function nextwd_ex1()
%
%         %
%         % Local variables.
%         %
%
%         rest = '  Now is the time,  for all good men to come.';
%
%         fprintf( 'Next   Rest of the string\n' )
%         fprintf( '-----  ------------------------------------------\n' )
%
%         while ~ isequal( rest, '' )
%
%            string = rest;
%            [next, rest] = cspice_nextwd( string );
%
%            fprintf( '%-5s  %s\n', next, rest )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Next   Rest of the string
%      -----  ------------------------------------------
%      Now    is the time,  for all good men to come.
%      is     the time,  for all good men to come.
%      the    time,  for all good men to come.
%      time,  for all good men to come.
%      for    all good men to come.
%      all    good men to come.
%      good   men to come.
%      men    to come.
%      to     come.
%      come.
%
%
%-Particulars
%
%   cspice_nextwd is used primarily for parsing input commands consisting
%   of one or more words, where a word is defined to be any sequence
%   of consecutive non-blank characters. Successive calls to cspice_nextwd,
%   each using the previous value of `rest' as the input string, allow
%   the calling routine to neatly parse and process one word at a
%   time.
%
%   cspice_nextwd cuts the input string into two pieces, and returns them
%   separately. The first piece is the first word in the string.
%   (Leading blanks are ignored. The first word, which is returned in
%   the output argument `next', runs from the first non-blank character
%   in the string up to the first blank that follows it.) The second
%   piece is whatever is left after the first word is removed. The
%   second piece is left justified, to simplify later calls to cspice_nextwd.
%
%-Exceptions
%
%   1)  If the input argument `string' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   2)  If the input argument `string' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
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
%   -Mice Version 1.0.0, 23-JUN-2021 (JDR)
%
%-Index_Entries
%
%   next word in a character_string
%
%-&
function [next, rest] = cspice_nextwd( string )

   switch nargin
      case 1

         string = zzmice_str(string, true);

      otherwise

         error ( 'Usage: [`next`, `rest`] = cspice_nextwd( `string` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [next, rest] = mice('nextwd_c', string);
   catch spiceerr
      rethrow(spiceerr)
   end
