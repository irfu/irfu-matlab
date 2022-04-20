%-Abstract
%
%   CSPICE_NTHWD returns the nth word in a character string, and its location
%   in the string.
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
%      nth      the index of the word to be returned.
%
%               [1,1] = size(nth); int32 = class(nth)
%
%               (One for the first word, two for the second, and so on.)
%
%   the call:
%
%      [word, loc] = cspice_nthwd( string, nth )
%
%   returns:
%
%      word     the `nth' word in `string'.
%
%               [1,c2] = size(word); char = class(word)
%
%               If `string' is blank or empty, or `nth' is non-positive or too
%               large, `word' is empty.
%
%               `word' may overwrite `string'.
%
%      loc      the location of `word' in `string'.
%
%               [1,1] = size(loc); int32 = class(loc)
%
%               (That is, `word' begins at string(loc)). If `string' is blank
%               or empty, or `nth' is non-positive or too large, `loc' is
%               zero.
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
%   1) Given a character string, get the N'th word within, and the
%      word's location.
%
%
%      Example code begins here.
%
%
%      function nthwd_ex1()
%
%         %
%         % Local parameters.
%         %
%         STRING =   ' Now is the time,   for all good men     to come.';
%
%         %
%         % Local variables.
%         %
%
%         for nth=-1:11
%
%            [word, loc] = cspice_nthwd( STRING, nth );
%            fprintf( 'Word # %2d  is <%s>, starting at position %2d\n',   ...
%                                           nth, word, loc )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Word # -1  is <>, starting at position  0
%      Word #  0  is <>, starting at position  0
%      Word #  1  is <Now>, starting at position  2
%      Word #  2  is <is>, starting at position  6
%      Word #  3  is <the>, starting at position  9
%      Word #  4  is <time,>, starting at position 13
%      Word #  5  is <for>, starting at position 21
%      Word #  6  is <all>, starting at position 25
%      Word #  7  is <good>, starting at position 29
%      Word #  8  is <men>, starting at position 34
%      Word #  9  is <to>, starting at position 42
%      Word # 10  is <come.>, starting at position 45
%      Word # 11  is <>, starting at position  0
%
%
%-Particulars
%
%   cspice_nthwd, like cspice_nextwd, is useful primarily for parsing input
%   commands consisting of one or more words, where a word is defined to be a
%   maximal sequence of consecutive non-blank characters. Each word is
%   bounded on both sides by a blank character, or by the start or end of the
%   input string. Successive calls to cspice_nextwd allow the calling routine
%   to neatly parse and process one word at a time.
%
%   The chief difference between the two routines is that
%   cspice_nthwd allows the calling routine to access the words making
%   up the input string in random order. (cspice_nextwd allows only
%   sequential access.)
%
%   cspice_nthwd may be more efficient than cspice_nextwd, since cspice_nthwd
%   doesn't update an output string consisting of the remaining, unparsed
%   string.
%
%-Exceptions
%
%   1)  If any of the input arguments, `string' or `nth', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `string' or `nth', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%   -Mice Version 1.0.0, 01-NOV-2021 (JDR)
%
%-Index_Entries
%
%   n'th word in a character_string
%
%-&
function [word, loc] = cspice_nthwd( string, nth )

   switch nargin
      case 2

         string = zzmice_str(string, true);
         nth    = zzmice_int(nth);

      otherwise

         error ( 'Usage: [`word`, loc] = cspice_nthwd( `string`, nth )' )

   end

   %
   % Call the MEX library.
   %
   try
      [word, loc] = mice('nthwd_c', string, nth);
   catch spiceerr
      rethrow(spiceerr)
   end
