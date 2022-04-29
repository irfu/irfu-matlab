%-Abstract
%
%   CSPICE_BLTFRM returns a SPICE set containing the frame IDs of all built-
%   in frames of a specified class.
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
%      frmcls   an integer code specifying the frame class or classes for
%               which built-in frame ID codes are requested.
%
%               [1,1] = size(frmcls); int32 = class(frmcls)
%
%               `frmcls' may designate a single class or "all classes."
%
%               The Mice parameter definitions file MiceFrm.m declares
%               parameters identifying frame classes. The supported values
%               and corresponding meanings of `frmcls' are
%
%                  Parameter            Value    Meaning
%                  ===================  =====    ==================
%                  SPICE_FRMTYP_ALL       -1     All frame classes.
%                  SPICE_FRMTYP_INERTL     1     Built-in inertial.
%                  SPICE_FRMTYP_PCK        2     PCK-based frame.
%                  SPICE_FRMTYP_CK         3     CK-based frame.
%                  SPICE_FRMTYP_TK         4     Fixed offset ("text
%                                                kernel") frame.
%                  SPICE_FRMTYP_DYN        5     Dynamic frame.
%                  SPICE_FRMTYP_SWTCH      6     Switch frame.
%
%      room     a parameter specifying the maximum number of elements that
%               can be accommodated by the dynamically allocated workspace
%               cell used internally by this routine.
%
%               [1,1] = size(room); int32 = class(room)
%
%               It's not necessary to compute an accurate estimate of how
%               many elements will be returned in `idset'; rather, the
%               user can pick a size considerably larger than what's
%               really required.
%
%   the call:
%
%      [idset] = cspice_bltfrm( frmcls, room )
%
%   returns:
%
%      idset    a SPICE set containing the ID codes of all built-in reference
%               frames of the specified class or classes.
%
%               [r,1] = size(idset); int32 = class(idset)
%
%-Parameters
%
%   See the Mice parameter definitions file MiceFrm.m.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Display the IDs and names of all SPICE built-in frames.
%      Group the outputs by frame class. Also fetch and display
%      the entire set of IDs and names using the parameter
%      SPICE_FRMTYP_ALL.
%
%
%      Example code begins here.
%
%
%      function bltfrm_ex1()
%
%         %
%         % MiceUser is a file that makes certain variables global.
%         % You must call MiceUser to have access to the parameters used
%         % in this example.
%         %
%         MiceUser;
%
%         %
%         % Local parameters
%         %
%         NFRAME = ( SPICE_NFRAME_NINERT + SPICE_NFRAME_NNINRT );
%
%         %
%         % Get the Toolkit version number and display it.
%         %
%         fprintf( "Toolkit version: %s\n", cspice_tkvrsn( "TOOLKIT" ) );
%
%         %
%         % Fetch and display the frames of each class.
%         %
%         for i=1:7
%
%            if ( i < 7 )
%
%               %
%               % Fetch the frames of class i.
%               %
%               [idset] = cspice_bltfrm( i, NFRAME );
%
%               outlin  = sprintf( 'Number of frames of class %d: %d',     ...
%                                                  i, size( idset )(1) );
%
%            else
%
%               %
%               % Fetch IDs of all built-in frames.
%               %
%               [idset] = cspice_bltfrm( SPICE_FRMTYP_ALL, NFRAME );
%
%               outlin  = sprintf( 'Number of built-in frames: %d',        ...
%                                                  size( idset )(1) );
%
%            end
%
%            %
%            % Display the NAIF ID and name of a maximum of 5 frames
%            % per family.
%            %
%            fprintf( '\n' )
%            fprintf( '%s\n', outlin )
%            fprintf( '   Frame IDs and names\n' )
%
%            nfrms = min( [ 5, size( idset )(1) ] );
%
%            for j=1:nfrms
%
%               [frname] = cspice_frmnam( idset(j) );
%               fprintf( '%12d   %s\n', idset(j), frname )
%
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Toolkit version: CSPICE_N0067
%
%      Number of frames of class 1: 21
%         Frame IDs and names
%                 1   J2000
%                 2   B1950
%                 3   FK4
%                 4   DE-118
%                 5   DE-96
%
%      Number of frames of class 2: 105
%         Frame IDs and names
%             10001   IAU_MERCURY_BARYCENTER
%             10002   IAU_VENUS_BARYCENTER
%             10003   IAU_EARTH_BARYCENTER
%             10004   IAU_MARS_BARYCENTER
%             10005   IAU_JUPITER_BARYCENTER
%
%      Number of frames of class 3: 0
%         Frame IDs and names
%
%      Number of frames of class 4: 1
%         Frame IDs and names
%             10081   EARTH_FIXED
%
%      Number of frames of class 5: 0
%         Frame IDs and names
%
%      Number of frames of class 6: 0
%         Frame IDs and names
%
%      Number of built-in frames: 145
%         Frame IDs and names
%                 1   J2000
%                 2   B1950
%                 3   FK4
%                 4   DE-118
%                 5   DE-96
%
%
%      Note that the set of built-in frames, particularly the
%      non-inertial ones, will grow over time, so the output
%      shown here may be out of sync with that produced by a
%      current SPICE Toolkit. Only the first 5 frames of each
%      family are presented in the output.
%
%-Particulars
%
%   This routine has a counterpart
%
%      cspice_kplfrm
%
%   which fetches the frame IDs of all frames specified in the kernel
%   pool.
%
%-Exceptions
%
%   1)  If the input frame class argument is not defined in
%       MiceFrm.m, the error SPICE(BADFRAMECLASS) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the size of `idset' is too small to hold the requested frame
%       ID set, the error SPICE(SETTOOSMALL) is signaled by a routine
%       in the call tree of this routine.
%
%   3)  If any of the input arguments, `frmcls' or `room', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `frmcls' or `room', is not of
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
%   FRAMES.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
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
%   -Mice Version 1.0.0, 26-AUG-2021 (JDR)
%
%-Index_Entries
%
%   fetch IDs of built-in reference frames
%
%-&
function [idset] = cspice_bltfrm( frmcls, room )

   switch nargin
      case 2

         frmcls = zzmice_int(frmcls);
         room = zzmice_int(room, [1, int32(inf)] );

      otherwise

         error ( 'Usage: [idset] = cspice_bltfrm( frmcls, room )' )

   end

   %
   % Call the MEX library.
   %
   try
      [idset] = mice('bltfrm_c', frmcls, room);
   catch spiceerr
      rethrow(spiceerr)
   end
