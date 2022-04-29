%-Abstract
%
%   CSPICE_TKFRAM finds the position rotation matrix from a Text Kernel (TK)
%   frame with the specified frame class ID to its base frame.
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
%      frcode   the unique frame class ID of the TK frame for which data is
%               being requested.
%
%               [1,1] = size(frcode); int32 = class(frcode)
%
%               For TK frames the frame class ID is always equal to the
%               frame ID.
%
%   the call:
%
%      [rot, frame, found] = cspice_tkfram( frcode )
%
%   returns:
%
%      rot      a position rotation matrix that converts positions relative
%               to the TK frame given by its frame class ID, `frcode', to
%               positions relative to the base frame given by its frame ID,
%               `frame'.
%
%               [3,3] = size(rot); double = class(rot)
%
%               Thus, if a position S has components x,y,z in the TK
%               frame, then S has components x', y', z' in the base
%               frame.
%
%                  .-  -.     .-     -. .- -.
%                  | x' |     |       | | x |
%                  | y' |  =  |  rot  | | y |
%                  | z' |     |       | | z |
%                  `-  -'     `-     -' `- -'
%
%
%      frame    the ID code of the base reference frame to which `rot' will
%               transform positions.
%
%               [1,1] = size(frame); int32 = class(frame)
%
%      found    a logical indicating whether or not a frame definition for
%               the TK frame with the frame class ID, `frcode', was
%               constructed from kernel pool data.
%
%               [1,1] = size(found); logical = class(found)
%
%               If `rot' and `frame' were constructed, `found' will be
%               returned with the value true. Otherwise it will be returned
%               with the value false.
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
%   1) Compute the rotation from the DSS-34 topocentric frame to
%      its base Earth body-fixed frame and use it to determine the
%      geodetic latitude and longitude of the DSS-34 site.
%
%
%      Use the FK kernel below to load the required topocentric
%      reference frame definition for the DSS-34 site.
%
%         earth_topo_050714.tf
%
%
%      Example code begins here.
%
%
%      function tkfram_ex1()
%
%         %
%         % Local parameters
%         %
%         MYTOPO =   'DSS-34_TOPO';
%
%         %
%         % Load the FK that contains the topocentric reference
%         % frame definition for DSS-34.
%         %
%         cspice_furnsh( 'earth_topo_050714.tf' );
%
%         %
%         % The name of the topocentric frame is MYTOPO.
%         % First we get the ID code of the topocentric frame.
%         %
%         [frcode] = cspice_namfrm( MYTOPO );
%
%         %
%         % Next get the rotation from the topocentric frame to
%         % the body-fixed frame. We can use the TK frame ID in
%         % place of the TK frame class ID in this call because
%         % for TK frames these IDs are identical.
%         %
%         [rot, frame, found] = cspice_tkfram( frcode );
%
%         %
%         % Make sure the topocentric frame is relative to one of
%         % the Earth fixed frames.
%         %
%         [frname] = cspice_frmnam( frame );
%
%         if ( ~strcmp( frname, 'IAU_EARTH'   ) &&                         ...
%              ~strcmp( frname, 'EARTH_FIXED' ) &&                         ...
%              ~strcmp( frname, 'ITRF93'      ) )
%
%            fprintf( 'The frame %s does not appear to be\n', MYTOPO )
%            fprintf( 'defined relative to an Earth fixed frame.\n' )
%            STOP(  );
%
%         end
%
%         %
%         % Things look ok. Get the location of the Z-axis in the
%         % topocentric frame.
%         %
%         z = rot(:,3);
%
%         %
%         % Convert the `z' vector to latitude, longitude and radius.
%         %
%         [rad, lat, lon] = cspice_reclat( z );
%
%         fprintf( 'The geodetic coordinates of the center\n' )
%         fprintf( 'of the topographic frame are:\n' )
%         fprintf( '\n' )
%         fprintf( '   Latitude  (deg):  %19.13f\n', lat*cspice_dpr )
%         fprintf( '   Longitude (deg):  %19.13f\n', lon*cspice_dpr )
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
%      The geodetic coordinates of the center
%      of the topographic frame are:
%
%         Latitude  (deg):    148.9819650021110
%         Longitude (deg):    -35.3984778756552
%
%
%-Particulars
%
%   This routine is used to construct the rotation from some frame
%   that is a constant rotation offset from some other reference
%   frame. This rotation is derived from data stored in the kernel
%   pool.
%
%   This routine is intended to be used as a low level routine by the
%   frame system software. However, you could use this routine to
%   directly retrieve the rotation from a fixed offset TK frame to
%   its base frame.
%
%-Exceptions
%
%   1)  If some kernel variable associated with this frame is not
%       present in the kernel pool, or does not have the proper type
%       or dimension, an error is signaled by a routine in the call
%       tree of this routine. In such a case `found' will be set to
%       false.
%
%   2)  If the input `frcode' has the value 0, the error
%       SPICE(ZEROFRAMEID) is signaled by a routine in the call tree
%       of this routine. `found' will be set to false.
%
%   3)  If the name of the frame corresponding to `frcode' cannot be
%       determined, the error SPICE(INCOMPLETEFRAME) is signaled by a
%       routine in the call tree of this routine.
%
%   4)  If the frame given by `frcode' is defined relative to a frame
%       that is unrecognized, the error SPICE(BADFRAMESPEC) is
%       signaled by a routine in the call tree of this routine. `found'
%       will be set to false.
%
%   5)  If the kernel pool specification for the frame given by `frcode'
%       is not one of 'MATRIX', 'ANGLES' or 'QUATERNION', the error
%       SPICE(UNKNOWNFRAMESPEC) is signaled by a routine in the call
%       tree of this routine. `found' will be set to false.
%
%   6)  If the frame `frcode' is equal to the relative frame ID (i.e.
%       the frame is defined relative to itself), the error
%       SPICE(BADFRAMESPEC2) is signaled by a routine in the call tree
%       of this routine. `found' will be set to false.
%
%   7)  If name-based and ID-based forms of any TKFRAME_ keyword are
%       detected in the kernel pool at the same time, the error
%       SPICE(COMPETINGFRAMESPEC) is signaled by a routine in the call
%       tree of this routine. `found' will be set to false.
%
%   8)  If the input argument `frcode' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   9)  If the input argument `frcode' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   This routine makes use of the loaded text kernels to determine
%   the rotation from a constant offset TK frame to its base frame.
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
%
%-Version
%
%   -Mice Version 1.0.0, 22-JUN-2021 (JDR)
%
%-Index_Entries
%
%   Fetch the rotation and frame of a text kernel frame
%   Fetch the rotation and frame of a constant offset frame
%
%-&
function [rot, frame, found] = cspice_tkfram( frcode )

   switch nargin
      case 1

         frcode = zzmice_int(frcode);

      otherwise

         error ( 'Usage: [rot(3,3), frame, found] = cspice_tkfram( frcode )' )

   end

   %
   % Call the MEX library.
   %
   try
      [rot, frame, found] = mice('tkfram_c', frcode);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
