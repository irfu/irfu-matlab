%-Abstract
%
%   CSPICE_CKFROT finds the position rotation matrix from a C-kernel (CK)
%   frame with the specified frame class ID (CK ID) to the base frame of the
%   highest priority CK segment containing orientation data for this CK frame
%   at the time requested.
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
%      inst     the unique frame class ID (CK ID) of the CK frame for which
%               data is being requested.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%      et       the epoch for which the position rotation is desired.
%
%               [1,1] = size(et); double = class(et)
%
%               `et' should be given in seconds past the epoch of J2000
%               TDB.
%
%   the call:
%
%      [rotate, ref, found] = cspice_ckfrot( inst, et )
%
%   returns:
%
%      rotate   a position rotation matrix that converts positions relative
%               to the CK frame given by its frame class ID, `inst', to
%               positions relative to the base frame given by its frame ID,
%               `ref'.
%
%               [3,3] = size(rotate); double = class(rotate)
%
%               Thus, if a position S has components x,y,z in the CK
%               frame, then S has components x', y', z' in the base
%               frame.
%
%                  .-  -.     .-        -. .- -.
%                  | x' |     |          | | x |
%                  | y' |  =  |  rotate  | | y |
%                  | z' |     |          | | z |
%                  `-  -'     `-        -' `- -'
%
%
%      ref      the ID code of the base reference frame to which `rotate'
%               will transform positions.
%
%               [1,1] = size(ref); int32 = class(ref)
%
%      found    true if a record was found to satisfy the pointing request.
%
%               [1,1] = size(found); logical = class(found)
%
%               `found' will be false otherwise.
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
%   1) Use cspice_ckfrot to compute the instantaneous angular velocity
%      vector for the Mars Global Surveyor (MGS) spacecraft frame,
%      'MGS_SPACECRAFT', relative to the inertial frame used as the
%      base frame in CK files containing MGS spacecraft orientation
%      at 2003-JUL-25 13:00:00. The frame class ID (CK ID) for the
%      'MGS_SPACECRAFT' frame is -94000.
%
%
%      Suppose that r(t) is the rotation matrix whose columns
%      represent the inertial pointing vectors of the MGS spacecraft
%      axes at time `t'.
%
%      Then the angular velocity vector points along the vector given
%      by:
%
%                              T
%          limit  axis( r(t+h)r )
%          h-->0
%
%
%      And the magnitude of the angular velocity at time `t' is given
%      by:
%
%                              T
%          d angle ( r(t+h)r(t) )
%         ------------------------   at   h = 0
%                    dh
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: ckfrot_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            naif0012.tls                  Leapseconds
%            mgs_sclkscet_00061.tsc        MGS SCLK coefficients
%            mgs_sc_ext12.bc               MGS s/c bus attitude
%
%         \begindata
%
%         KERNELS_TO_LOAD = ( 'naif0012.tls',
%                             'mgs_sclkscet_00061.tsc',
%                             'mgs_sc_ext12.bc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function ckfrot_ex1()
%
%         %
%         % Local parameters.
%         %
%         EPOCH =   '2003-JUL-25 13:00:00';
%         INST  =   -94000;
%
%         %
%         % Load the required LSK, SCLK and CK. Use a
%         % meta-kernel for convenience.
%         %
%         cspice_furnsh( 'ckfrot_ex1.tm' );
%
%         %
%         % First convert the time to seconds past J2000. Set the
%         % delta time (1 ms).
%         %
%         [et] = cspice_str2et( EPOCH );
%         h    = 1.e-3;
%
%         %
%         % Now, look up the rotation from the MGS spacecraft
%         % frame specified by its frame class ID (CK ID) to a
%         % base reference frame (returned by cspice_ckfrot), at `et'
%         % and et+h.
%         %
%         [ret,  ref,  found]  = cspice_ckfrot( INST, et   );
%         [reth, refh, foundh] = cspice_ckfrot( INST, et+h );
%
%         %
%         % If both rotations were computed and if the base
%         % reference frames are the same, compute the
%         % instantaneous angular velocity vector.
%         %
%         if ( found && foundh && ref == refh )
%
%            %
%            % Compute the infinitesimal rotation r(t+h)r(t)**T.
%            %
%            infrot = reth * ret';
%
%            %
%            % Compute the `axis' and `angle' of the infinitesimal
%            % rotation.
%            %
%            [axis, angle] = cspice_raxisa( infrot );
%
%            %
%            % Scale `axis' to get the angular velocity vector.
%            %
%            angvel = angle/h * axis;
%
%            %
%            % Output the results.
%            %
%            fprintf( 'Instantaneous angular velocity vector:\n' )
%            fprintf( '%15.10f %14.10f %14.10f\n', angvel )
%            fprintf( 'Reference frame ID: %4d\n', ref )
%
%         else
%
%            fprintf( 'ERROR: data not found or frame mismatch.\n' )
%
%         end
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
%      Instantaneous angular velocity vector:
%         0.0001244121   0.0008314866   0.0003028634
%      Reference frame ID:    1
%
%
%-Particulars
%
%   cspice_ckfrot searches through loaded CK files to satisfy a pointing
%   request. Last-loaded files are searched first, and individual
%   files are searched in backwards order, giving priority to
%   segments that were added to a file later than the others.
%
%   The search ends when a segment is found that can give pointing
%   for the specified CK frame at the request time.
%
%   Segments with and without angular velocities are considered by
%   this routine.
%
%   This routine uses the cspice_ckmeta routine to determine the SCLK ID
%   used to convert the input `et' time to the encoded SCLK time used
%   to look up pointing data in loaded CK files.
%
%-Exceptions
%
%   1)  If no CK files were loaded prior to calling this routine, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   2)  If no SCLK correlation data needed to read CK files were
%       loaded prior to calling this routine, an error is signaled by
%       a routine in the call tree of this routine.
%
%   3)  If the input time `et' cannot be converted to an encoded SCLK
%       time, using SCLK data associated with `inst', an error is
%       signaled by a routine in the call tree of this routine.
%
%   4)  If any of the input arguments, `inst' or `et', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `inst' or `et', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   cspice_ckfrot searches through loaded CK files to locate a segment that
%   can satisfy the request for position rotation data for the CK
%   frame with the specified frame class ID at time `et'. You must load
%   a CK file containing such data before calling this routine. You
%   must also load SCLK and possibly LSK files needed to convert the
%   input `et' time to the encoded SCLK time with which the orientation
%   data stored inside that CK is tagged.
%
%-Restrictions
%
%   1)  A CK file must be loaded prior to calling this routine.
%
%   2)  LSK and SCLK files needed for time conversions must be loaded
%       prior to calling this routine.
%
%-Required_Reading
%
%   CK.REQ
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
%   -Mice Version 1.0.0, 13-DEC-2021 (JDR)
%
%-Index_Entries
%
%   get instrument frame rotation and reference frame
%
%-&
function [rotate, ref, found] = cspice_ckfrot( inst, et )

   switch nargin
      case 2

         inst = zzmice_int(inst);
         et   = zzmice_dp(et);

      otherwise

         error ( [ 'Usage: [rotate(3,3), ref, found] = '                    ...
                   'cspice_ckfrot( inst, et )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [rotate, ref, found] = mice('ckfrot_c', inst, et);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
