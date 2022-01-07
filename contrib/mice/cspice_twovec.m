%-Abstract
%
%   CSPICE_TWOVEC calculates the transformation matrix to the
%   right-handed reference frame having an input vector as a
%   specified axis and a second input vector lying in a
%   define coordinate plane.
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
%      axdef    the principal axes of a coordinate frame.
%
%               [3,1] = size(axdef); double = class(axdef)
%
%      indexa   the index signifying which of the three coordinate
%               axes contains 'axdef' (1, 2 or 3).
%
%               [1,1]   = size(indexa); int32 = class(indexa)
%
%                  If 'indexa' is 1 then axdef defines the X axis of the
%                  coordinate frame.
%
%                  If 'indexa' is 2 then axdef defines the Y axis of the
%                  coordinate frame.
%
%                  If 'indexa' is 3 then axdef defines the Z axis of the
%                  coordinate frame.
%
%      plndef   a vector in the same plane as 'axdef'. 'axdef' and
%               'plndef' must be linearly independent.
%
%               [3,1] = size(plndef); double = class(plndef)
%
%      indexp   the index signifying the second principle axis,
%               orthogonal to 'axdef' (1, 2 or 3).
%
%               [1,1]   = size(indexp); int32 = class(indexp)
%
%                  If 'indexp' is 1, the second axis of the principal
%                  plane is the X-axis.
%
%                  If 'indexp' is 2, the second axis of the principal
%                  plane is the Y-axis.
%
%                  If 'indexp' is 3, the second axis of the principal plane
%                  is the Z-axis.
%
%   the call:
%
%      mout = cspice_twovec( axdef, indexa, plndef, indexp)
%
%   returns:
%
%      mout     a double precision 3x3 array defining a rotation matrix from
%               the frame of the original vectors to the new frame
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
%   1) Calculate the transformation matrix to the right-handed
%      reference frame having the +i unit vector as primary axis,
%      aligned to the frame's +X axis, and the -j unit vector as
%      secondary axis, aligned to the +Y axis.
%
%      Example code begins here.
%
%
%      function twovec_ex1()
%
%         %
%         % A trivial example. Define the reference vectors...
%         %
%         %  The i unit vector
%         %
%         axdef  = [ 1.; 0; 0.];
%         indexa = 1 ;
%
%         %
%         %  The -j unit vector. For this example, any vector
%         %  in the x-y plane linearly independent of 'axdef'
%         %  will suffice.
%         %
%         plndef = [ 0.; -1.; 0.];
%         indexp = 2;
%
%         %
%         % Calculate the transformation matrix. The new frame
%         % has 'axdef' as axis 'indexa', with 'plndef' in the same
%         % plane, the direction axis 'indexp' in that plane
%         % and orthogonal to 'axdef'. A third direction vector
%         % completes the right handed frame.
%         %
%         mout = cspice_twovec( axdef, indexa, plndef, indexp );
%         fprintf( '%15.7f %15.7f %15.7f\n', mout(1,:));
%         fprintf( '%15.7f %15.7f %15.7f\n', mout(2,:));
%         fprintf( '%15.7f %15.7f %15.7f\n', mout(3,:));
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%            1.0000000       0.0000000       0.0000000
%            0.0000000      -1.0000000       0.0000000
%            0.0000000       0.0000000      -1.0000000
%
%
%-Particulars
%
%   Given two linearly independent vectors there is a unique
%   right-handed coordinate frame having:
%
%      1) `axdef' lying along the `indexa' axis.
%
%      2) `plndef' lying in the indexa-indexp coordinate plane.
%
%   This routine determines the transformation matrix that transforms
%   from coordinates used to represent the input vectors to the
%   the system determined by `axdef' and `plndef'. Thus a vector
%   (x,y,z) in the input coordinate system will have coordinates
%
%                            t
%              mout * (x,y,z)
%
%   in the frame determined by `axdef' and `plndef'.
%
%-Exceptions
%
%   1)  If `indexa' or `indexp' is not in the set {1,2,3}, the error
%       SPICE(BADINDEX) is signaled by a routine in the call tree of
%       this routine.
%
%   2)  If `indexa' and `indexp' are the same, the error
%       SPICE(UNDEFINEDFRAME) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the cross product of the vectors `axdef' and `plndef' is zero,
%       the error SPICE(DEPENDENTVECTORS) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If any of the input arguments, `axdef', `indexa', `plndef' or
%       `indexp', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   5)  If any of the input arguments, `axdef', `indexa', `plndef' or
%       `indexp', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's meta-kernel and reformatted example's output.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 12-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 10-JAN-2006 (EDW)
%
%-Index_Entries
%
%   define an orthonormal frame from two vectors
%
%-&

function [mout] = cspice_twovec(axdef, indexa, plndef, indexp)

   switch nargin
      case 4

         axdef  = zzmice_dp(axdef);
         indexa = zzmice_int(indexa);
         indexp = zzmice_int(indexp);
         plndef = zzmice_dp(plndef);

      otherwise

         error ( ['Usage: [mout(3,3)] = cspice_twovec( axdef(3), ' ...
                                       'indexa, plndef(3), indexp)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [mout] = mice('twovec_c', axdef, indexa, plndef, indexp);
   catch spiceerr
      rethrow(spiceerr)
   end


