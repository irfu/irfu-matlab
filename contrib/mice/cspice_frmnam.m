%-Abstract
%
%   CSPICE_FRMNAM retrieves the name of a reference frame associated with a
%   SPICE frame ID code.
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
%      frcode   value(s) defining a SPICE reference frame ID code(s).
%
%               [1,n] = size(frcode); int32 = class(frcode)
%
%   the call:
%
%      [frname] = cspice_frmnam( frcode )
%
%   returns:
%
%      frname   the frame name(s) corresponding to the `frcode' code(s).
%
%               [n,c1] = size(frname); char = class(frname)
%
%               If `frcode' is not recognized as the name of a known reference
%               frame, `frname' will be returned as an empty string.
%
%               `frname' returns with the same vectorization measure, N,
%               as `frcode'.
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
%   1) Given a set of SPICE frame IDs, retrieve their associated frame
%      names.
%
%      Example code begins here.
%
%
%      function frmnam_ex1()
%
%         %
%         % Retrieve frame name for a scalar code.
%         %
%         code = 13000;
%
%         %
%         % Output the frame name corresponding to `code'.
%         %
%         frname = cspice_frmnam( code );
%         disp(' ID code  Frame Name' )
%         disp('--------  ----------' )
%         disp('Scalar' )
%         fprintf('%8d  %s\n', code, frname )
%
%
%         %
%         % Retrieve frame information for a vector of codes.
%         %
%         disp('Vector' )
%         codes = [1, 2, 3, 4, 5, 245];
%
%         %
%         % Output the frame names corresponding to `codes'.
%         %
%         frname = cspice_frmnam( codes );
%
%         for i=1:numel( codes )
%
%            fprintf( '%8d  %s\n', codes(i), frname(i, :) )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%       ID code  Frame Name
%      --------  ----------
%      Scalar
%         13000  ITRF93
%      Vector
%             1  J2000
%             2  B1950
%             3  FK4
%             4  DE-118
%             5  DE-96
%           245
%
%
%      Note that 245 does not correspond to any known frame in SPICE,
%      and therefore a blank string is returned.
%
%-Particulars
%
%   This routine retrieves the name of a reference frame associated
%   with a SPICE frame ID code.
%
%   The ID codes stored locally are scanned for a match with `frcode'.
%   If a match is found, the name stored locally will be returned
%   as the name for the frame.
%
%   If `frcode' is not a member of the list of internally stored
%   ID codes, the kernel pool will be examined to see if the
%   variable
%
%      FRAME_idcode_NAME
%
%   is present (where idcode is the decimal character equivalent
%   of `frcode'). If the variable is located and it has both
%   character type and dimension 1, the string value of the
%   kernel pool variable is returned as the name of the reference
%   frame.
%
%   Note that because the local information is always examined
%   first and searches of the kernel pool are performed only
%   after exhausting local information, it is not possible to
%   override the local name for any reference frame that is
%   known by this routine.
%
%-Exceptions
%
%   1)  If `frcode' is not recognized as the name of a known reference
%       frame, `frname' will be returned as an empty string.
%
%   2)  If the input argument `frcode' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `frcode' is not of the expected type, or
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
%   FRAMES.REQ
%   MICE.REQ
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "frmname" to "frname".
%
%       Edited the header to comply with NAIF standard. Reformatted
%       example's output and added problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 12-MAR-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   Frame ID code to frame name translation
%
%-&

function [frname] = cspice_frmnam( frcode )

   switch nargin
      case 1

         frcode = zzmice_int(frcode);

      otherwise

         error( 'Usage: [_`frname`_] = cspice_frmnam(_frcode_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [frname] = mice('frmnam_c', frcode);
   catch spiceerr
      rethrow(spiceerr)
   end




