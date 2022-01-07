%-Abstract
%
%   CSPICE_CCIFRM returns the frame name, frame id, and center associated with
%   a given frame class and class id.
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
%     frclss   is the class or type of the frame. This identifies which
%              subsystem will be used to perform frame transformations.
%
%              [1,1] = size(frclss); int32 = class(frclss)
%
%     clssid   is the ID code used for the frame within its class. This
%              may be different from the frame ID code.
%
%              [1,1] = size(clssid); int32 = class(clssid)
%
%   the call:
%
%      [frcode, frname, cent, found] = cspice_ccifrm( frclss, clssid )
%
%   returns:
%
%      frcode  is the frame ID code for the reference frame
%              identified by `frclss' and `clssid'.
%
%              [1,1] = size(frcode); int32 = class(frcode)
%
%      frname  is the name of the frame identified by
%              `frclss' and `clssid'.
%
%              [1,c1] = size(frname); char = class(frname)
%
%              If `frname' does not have enough room to hold
%              the full name of the frame, the name will be truncated
%              on the right.
%
%      cent    is the body ID code for the center of the reference
%              frame identified  by `frclss' and `clssid'.
%
%              [1,1] = size(cent); int32 = class(cent)
%
%      found   is true if a valid frame specification
%              corresponding to the input frame class and frame class
%              ID is available, in which case the other outputs are
%              valid. Otherwise, `found' is returned with the value
%              false.
%
%              [1,1] = size(found); logical = class(found)
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
%   1) The following code example demonstrates how to find the frame
%      information about a frame by its ID using cspice_frinfo and
%      by its class and class ID using cspice_ccifrm.
%
%      Example code begins here.
%
%
%      function ccifrm_ex1()
%
%         frcode1 = cspice_namfrm( 'ITRF93');
%         [center1, clss, clss_ID, found] = cspice_frinfo( frcode1 );
%
%         if ( ~found )
%               error('No info found for ITRF93' )
%         end
%
%         fprintf( 'Frame ITRF93 info:\n'          )
%         fprintf( '   Frame Code: %d\n', frcode1 )
%         fprintf( '   Center ID : %d\n', center1 )
%         fprintf( '   Class     : %d\n', clss    )
%         fprintf( '   Class ID  : %d\n', clss_ID )
%
%         [frcode2, frname, center2, found] = cspice_ccifrm( clss, clss_ID );
%
%         if ( ~found )
%               error('No info found for type 2 frame 3000.' )
%         end
%
%         fprintf( 'Type 2 frame 3000 info:\n'    )
%         fprintf( '   Frame name: %s\n', frname  )
%         fprintf( '   Frame Code: %d\n', frcode2 )
%         fprintf( '   Center ID : %d\n', center2 )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Frame ITRF93 info:
%         Frame Code: 13000
%         Center ID : 399
%         Class     : 2
%         Class ID  : 3000
%      Type 2 frame 3000 info:
%         Frame name: ITRF93
%         Frame Code: 13000
%         Center ID : 399
%
%
%-Particulars
%
%   This routine allows the user to determine the frame associated
%   with a given frame class and class ID code. The kernel pool is
%   searched first for a matching frame; if no match is found, then
%   the set of built-in frames is searched.
%
%   Since the neither the frame class nor the class ID are primary
%   keys, searching for matching frames is a linear (and therefore
%   typically slow) process.
%
%-Exceptions
%
%   1)  This routine assumes that the first frame found with matching
%       class and class ID is the correct one. SPICE's frame system
%       does not diagnose the situation where there are multiple,
%       distinct frames with matching classes and class ID codes, but
%       this situation could occur if such conflicting frame
%       specifications are loaded via one or more frame kernels. The
%       user is responsible for avoiding such frame specification
%       conflicts.
%
%   2)  If a frame class assignment is found that associates a string
%       (as opposed to numeric) value with a frame class keyword, the
%       error SPICE(INVALIDFRAMEDEF) is signaled by a routine in the
%       call tree of this routine.
%
%   3)  If a frame class assignment is found that matches the input
%       class, but a corresponding class ID assignment is not
%       found in the kernel pool, the error SPICE(INVALIDFRAMEDEF)
%       is signaled by a routine in the call tree of this routine.
%
%   4)  If a frame specification is found in the kernel pool with
%       matching frame class and class ID, but either the frame name
%       or frame ID code are not found, the error
%       SPICE(INVALIDFRAMEDEF) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If a frame specification is found in the kernel pool with
%       matching frame class and class ID, but the frame center
%       is not found, an error is signaled by a routine
%       in the call tree of this routine.
%
%   6)  If any of the input arguments, `frclss' or `clssid', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   7)  If any of the input arguments, `frclss' or `clssid', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  See item (1) in the -Exceptions section above.
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
%   M. Liukis           (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Changed the output argument name "center" to "cent" for
%       consistency with other routines.
%
%       Edited the -Examples section to comply with NAIF standard. Updated
%       code example and added problem's statement. Added -Parameters,
%       -Exceptions, -Files, -Restrictions, -Literature_References and
%       -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 01-MAR-2017 (ML) (EDW)
%
%-Index_Entries
%
%   Find info associated with a frame class and class id
%   Map frame class and class id to frame info
%   Map frame class and class id to frame name, id, and center
%
%-&

function [frcode, frname, cent, found] = cspice_ccifrm( frclss, clssid )

   switch nargin
      case 2

         frclss = zzmice_int(frclss);
         clssid = zzmice_int(clssid);

      otherwise

         error( ['Usage: [ frcode, `frname`, cent, found] = ' ...
                'cspice_ccifrm( frclss, clssid)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [ccifrm, cent] = mice( 'ccifrm_s', frclss, clssid );
      frcode   = reshape( [ccifrm.code],  1, [] );
      frname   = char( ccifrm.name );
      found    = reshape( [ccifrm.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end
