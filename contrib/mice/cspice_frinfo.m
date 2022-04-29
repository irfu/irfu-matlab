%-Abstract
%
%   CSPICE_FRINFO retrieves the minimal attributes associated with a
%   frame needed for converting transformations to and from it.
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
%      frcode   a SPICE ID(s) for some reference frame(s).
%
%               [1,n] = size(frcode); int32 = class(frcode)
%
%   the call:
%
%      [cent, frclss, clssid, found] = cspice_frinfo( frcode )
%
%   returns:
%
%      cent     the SPICE body ID(s) for the center of the reference frame(s)
%               (if such an ID is appropriate).
%
%               [1,n] = size(cent); int32 = class(cent)
%
%      frclss   the class ID(s) or type(s) of the frame. This identifies which
%               subsystem will perform frame transformations.
%
%               [1,n] = size(frclss); int32 = class(frclss)
%
%      clssid   the ID(s) used for the frame within its class. This may be
%               different from the frame ID.
%
%               [1,n] = size(clssid); int32 = class(clssid)
%
%      found    flag(s) returning true if `cent', `frclss' and `frcode' are
%               available, false if not.
%
%               [1,n] = size(found); logical = class(found)
%
%               `cent', `frclss', `clssid' and `found' return with the same
%               vectorization measure (N) as `frcode'.
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
%   1) Given a set of frame IDs, retrieve the SPICE body ID associated
%      with the frame's center, the frame class (or type of the frame),
%      and the ID used for the frame within its class.
%
%
%      Example code begins here.
%
%
%      function frinfo_ex1()
%
%         %
%         % Retrieve frame information for a scalar 'code'.
%         %
%         code = 13000;
%
%         [cent, frclss, clssid, found] = cspice_frinfo( code );
%         disp(' code   center  class  class_ID' )
%         disp('-----  -------  -----  --------' )
%         disp('Scalar:' )
%         if ( found )
%            fprintf( '%5d  %7d  %5d  %6d\n', ...
%                     code, cent, frclss, clssid   );
%         end
%
%         %
%         % Retrieve frame information for a vector of 'codes'.
%         %
%         disp('Vector:' )
%         codes = [1, 2, 3, 4, 5, 245];
%
%         [cent, frclss, clssid, found] = cspice_frinfo( codes );
%
%         for i=1:numel( codes )
%
%            if ( found(i) )
%               fprintf( '%5d  %7d  %5d  %6d\n', ...
%                       codes(i), cent(i), frclss(i), clssid(i) )
%            else
%               fprintf( 'No data available for frame ID %d\n', codes(i) )
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%       code   center  class  class_ID
%      -----  -------  -----  --------
%      Scalar:
%      13000      399      2    3000
%      Vector:
%          1        0      1       1
%          2        0      1       2
%          3        0      1       3
%          4        0      1       4
%          5        0      1       5
%      No data available for frame ID 245
%
%
%-Particulars
%
%   This is a low level routine needed by state transformation
%   software to transform states and attitudes between different
%   reference frames.
%
%   The routine first examines local "hard-coded" information about
%   reference frames to see if the requested frame belongs to this
%   set. If it does that information is returned.
%
%   If the requested information is not stored locally, the routine
%   then examines the kernel pool to see if the requested information
%   is stored there. If it is and has the expected format, the data
%   is retrieved and returned.
%
%-Exceptions
%
%   1)  If a frame definition is encountered that does not define a
%       central body for the frame, an error is signaled by a routine
%       in the call tree of this routine.
%
%   2)  If a frame definition is encountered that does not define
%       a class for the frame, an error is signaled by a routine
%       in the call tree of this routine.
%
%   3)  If a frame definition is encountered that does not define a
%       class ID for the frame, an error is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If a kernel variable defining a frame name is found, but that
%       variable has dimension greater than 1, the error
%       SPICE(INVALIDDIMENSION) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the input argument `frcode' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   6)  If the input argument `frcode' is not of the expected type, or
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Changed the output argument name "clss" to "frclss" for
%       consistency with other routines.
%
%       Edited the -Examples section to comply with NAIF standard.
%       Reformatted example's output and added problem statement.
%       Added case for demonstrating usage of `found' flag.
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
%   fetch reference frame attributes
%
%-&

function [cent, frclss, clssid, found] = cspice_frinfo( frcode )

   switch nargin
      case 1

         frcode = zzmice_int(frcode);

      otherwise

         error( ['Usage: [_cent_, _frclss_, _clssid_, _found_] = ' ...
                                           'cspice_frinfo(_frcode_)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [frinfo] = mice('frinfo_s', frcode);

      cent   = reshape( [frinfo.center],   1, [] );
      frclss = reshape( [frinfo.class],    1, [] );
      clssid = reshape( [frinfo.class_ID], 1, [] );
      found  = reshape( [frinfo.found],    1, [] );

   catch spiceerr
      rethrow(spiceerr)
   end
