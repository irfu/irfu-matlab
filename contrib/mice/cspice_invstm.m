%-Abstract
%
%   CSPICE_INVSTM returns the inverse of a state transformation matrix.
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
%      mat      a state transformation matrix for converting states relative
%               to one frame to states relative to another.
%
%               [6,6] = size(mat); double = class(mat)
%
%               The state transformation of a state vector, `s', is
%               performed by the matrix-vector product.
%
%                  mat * s.
%
%               For `mat' to be a "true" state transformation matrix
%               it must have the form
%
%                   .-            -.
%                   |       :      |
%                   |   r   :   0  |
%                   |.......:......|
%                   |       :      |
%                   |  w*r  :   r  |
%                   |       :      |
%                   `-            -'
%
%               where `r' is a 3x3 rotation matrix, 0 is the 3x3 zero
%               matrix and `w' is a 3x3 skew-symmetric matrix.
%
%               NOTE: no checks are performed on `mat' to ensure that it
%                     does indeed have the form described above.
%
%   the call:
%
%      [invmat] = cspice_invstm( mat )
%
%   returns:
%
%      invmat   the inverse of `mat' under the operation of matrix
%               multiplication.
%
%               [6,6] = size(invmat); double = class(invmat)
%
%               If `mat' has the form described above, then `invmat' has
%               the form shown below.
%
%                  .-             -.
%                  |     t  :      |
%                  |    r   :   0  |
%                  |........:......|
%                  |      t :    t |
%                  | (w*r)  :   r  |
%                  |        :      |
%                  `-             -'
%
%               (The superscript "t" denotes the matrix transpose
%               operation.)
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
%   1) Suppose you have a geometric state of a spacecraft in Earth
%      body-fixed reference frame and wish to express this state
%      relative to an Earth centered J2000 frame. The following
%      example code illustrates how to carry out this computation.
%
%      Use the PCK kernel below to load the required high precision
%      orientation of the ITRF93 Earth body-fixed reference frame.
%      Note that the body ID code used in this file for the Earth is
%      3000.
%
%         earth_720101_070426.bpc
%
%
%      Example code begins here.
%
%
%      function invstm_ex1()
%
%         %
%         % Define the state of the spacecraft, in km and
%         % km/s, and the `et' epoch, in seconds past J2000.
%         %
%         et    = 0.0;
%         state = [ 175625246.29100420, 164189388.12540060,                ...
%                   -62935198.26067264,     11946.73372264,                ...
%                      -12771.29732556,        13.84902914 ]';
%
%         %
%         % Load the required high precision Earth PCK.
%         %
%         cspice_furnsh( 'earth_720101_070426.bpc' );
%
%         %
%         % First get the state transformation from J2000 frame
%         % to Earth body-fixed frame at the time of interest `et'.
%         % The body ID code used in high precision PCK files for
%         % the Earth is 3000; this number indicates that the
%         % terrestrial frame used is ITRF93.
%         %
%         earth = 3000;
%         [mat] = cspice_tisbod( 'J2000', earth, et );
%
%         %
%         % Get the inverse of `mat'.
%         %
%         [invmat] = cspice_invstm( mat );
%
%         %
%         % Transform from bodyfixed state to inertial state.
%         %
%         istat1 = invmat * state;
%
%         %
%         % Print the resulting state.
%         %
%         fprintf( [ 'Input state in Earth centered J2000 frame, using',   ...
%                    ' cspice_invstm:\n' ]                               )
%         fprintf( '   Position: %15.3f %15.3f %15.3f\n', istat1(1:3) )
%         fprintf( '   Velocity: %15.3f %15.3f %15.3f\n', istat1(4:6) )
%
%         %
%         % Compute the same state using cspice_sxform.
%         %
%         [xmat] = cspice_sxform( 'ITRF93', 'J2000', et );
%         istat2 = xmat * state;
%
%         fprintf( '\n' )
%         fprintf( [ 'Input state in Earth centered J2000 frame, using',   ...
%                    ' cspice_sxform:\n' ]                               )
%         fprintf( '   Position: %15.3f %15.3f %15.3f\n', istat2(1:3) )
%         fprintf( '   Velocity: %15.3f %15.3f %15.3f\n', istat2(4:6) )
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
%      Input state in Earth centered J2000 frame, using cspice_invstm:
%         Position:   192681395.921  -143792821.383   -62934296.473
%         Velocity:          30.312          32.007          13.876
%
%      Input state in Earth centered J2000 frame, using cspice_sxform:
%         Position:   192681395.921  -143792821.383   -62934296.473
%         Velocity:          30.312          32.007          13.876
%
%
%-Particulars
%
%   Given a matrix for transforming states relative frame 1 to
%   states relative frame 2, the routine produces the inverse
%   matrix. That is, it returns the matrix for transforming states
%   relative to frame 2 to states relative to frame 1.
%
%   This special routine exists because unlike the inverse of a
%   rotation matrix, the inverse of a state transformation matrix,
%   is NOT simply the transpose of the matrix.
%
%-Exceptions
%
%   1)  No checks are performed to ensure that the input matrix is
%       indeed a state transformation matrix.
%
%   2)  If the input argument `mat' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   3)  If the input argument `mat' is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
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
%   ROTATION.REQ
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
%   -Mice Version 1.0.0, 25-NOV-2021 (JDR)
%
%-Index_Entries
%
%   inverse of state transformation matrix
%
%-&
function [invmat] = cspice_invstm( mat )

   switch nargin
      case 1

         mat = zzmice_dp(mat);

      otherwise

         error ( 'Usage: [invmat(6,6)] = cspice_invstm( mat(6,6) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [invmat] = mice('invstm_c', mat);
   catch spiceerr
      rethrow(spiceerr)
   end
