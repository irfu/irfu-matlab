%-Abstract
%
%   CSPICE_NAMFRM retrieves the SPICE frame ID code associated
%   with a frame name.
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
%      frname   the name(s) of some reference frame(s) (either inertial or
%               non-inertial).
%
%               [n,c1] = size(frname); char = class(frname)
%
%                  or
%
%               [1,n] = size(frname); cell = class(frname)
%
%               Leading blanks in `frname' are ignored as is character case.
%
%               Note that all legitimate frame names contain 32 or fewer
%               characters.
%
%   the call:
%
%      frcode = cspice_namfrm(frname)
%
%   returns:
%
%      frcode   the SPICE code(s) used for internal representation of the
%               named reference frame.
%
%               [1,n] = size(frcode); int32 = class(frcode)
%
%               If the name input through `frname' is not recognized, `frcode'
%               will be returned with a value of zero.
%
%               `frcode' returns with the same vectorization measure (N)
%               as `frname'.
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
%   1) Given a set of frame names, retrieve their associated SPICE
%      frame ID.
%
%      Example code begins here.
%
%
%      function namfrm_ex1()
%
%         %
%         % Retrieve frame information for a single frame.
%         %
%         name = 'ITRF93';
%
%         %
%         % Output the SPICE frame ID corresponding to 'name'.
%         %
%         frcode = cspice_namfrm( name );
%
%         disp('Frame Name  ID code' )
%         disp('----------  -------' )
%         disp('Scalar:' )
%         fprintf('%10s  %d\n', name, frcode )
%
%         %
%         % Retrieve frame information for a vector of names.
%         %
%         disp('Vector:' )
%         names = { 'J2000', 'IAU_MARS', 'FK4', 'ECLIPJ2000', 'MYFRAME' };
%
%         %
%         % Output the frame IDs corresponding to 'names'.
%         %
%         frcode = cspice_namfrm( names );
%
%         for i=1:numel( frcode )
%
%            if ( frcode(i) )
%               fprintf( '%10s  %d\n', char(names(i)), frcode(i) )
%            else
%               fprintf( 'No SPICE frame ID associated to the name  %s\n', ...
%                        char(names(i)) )
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Frame Name  ID code
%      ----------  -------
%      Scalar:
%          ITRF93  13000
%      Vector:
%           J2000  1
%        IAU_MARS  10014
%             FK4  3
%      ECLIPJ2000  17
%      No SPICE frame ID associated to the name  MYFRAME
%
%
%-Particulars
%
%   This is a low level interface routine intended primarily for
%   use within the SPK and CK systems to assist in the transformation
%   to user specified reference frames.
%
%   The routine first consults a stored list of reference frame
%   names in an attempt to determine the appropriate reference
%   frame code.
%
%   If this search is unsuccessful, the routine then examines the
%   kernel pool to determine whether or not a variable of the
%   form
%
%      'FRAME_' + frname
%
%      (where leading blanks of `frname' are ignored)
%
%   is present. If it is and the number of values associated with the
%   name is 1, this value is taken to be the frame ID code.
%
%   Note: It is NOT possible to override the default names and
%   ID codes stored locally in this routine by placing an
%   appropriately variable in the kernel pool with a different
%   ID code. The predefined values always take precedence.
%
%   Consult the frames.req required reading document for more details
%   about constructing your own frame definitions.
%
%-Exceptions
%
%   1)  If the input name is not recognized, `frcode' will be
%       returned with a value of 0.
%
%   2)  If the input argument `frname' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `frname' is not of the expected type, or
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
%   FRAMES.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   M. Liukis           (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Reformatted example's output and added problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 03-JAN-2016 (EDW) (ML)
%
%       Corrected minor typo, "or" rather than "of."
%
%   -Mice Version 1.0.0, 14-NOV-2014 (EDW) (SCK)
%
%-Index_Entries
%
%   frame name to frame ID code translation
%
%-&

function [frcode] = cspice_namfrm(frname)

   switch nargin
      case 1

         frname = zzmice_str(frname);

      otherwise

         error ( 'Usage: [_frcode_] = cspice_namfrm(_`frname`_)' )

   end

   try
      [frcode] = mice('namfrm_c',frname);
   catch spiceerr
      rethrow(spiceerr)
   end


