%-Abstract
%
%   CSPICE_KPLFRM returns a SPICE set containing the frame IDs of all
%   reference frames of a given class having specifications in the kernel
%   pool.
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
%               which frame ID codes are requested.
%
%               [1,1] = size(frmcls); int32 = class(frmcls)
%
%               The applicable reference frames are those having
%               specifications present in the kernel pool.
%
%               `frmcls' may designate a single class or "all
%               classes."
%
%               The Mice parameter definitions file MiceFrm.m declares
%               parameters identifying frame classes. The supported values
%               and corresponding meanings of `frmcls' are
%
%                  Parameter            Value   Meaning
%                  ===================  =====   ====================
%                  SPICE_FRMTYP_ALL       -1    All frame classes
%                                               specified in the
%                                               kernel pool. Class 1
%                                               is not included.
%
%                  SPICE_FRMTYP_INERTL     1    Built-in inertial.
%                                               No frames will be
%                                               returned in the
%                                               output set.
%
%                  SPICE_FRMTYP_PCK        2    PCK-based frame.
%
%                  SPICE_FRMTYP_CK         3    CK-based frame.
%
%                  SPICE_FRMTYP_TK         4    Fixed rotational
%                                               offset ("text
%                                               kernel") frame.
%
%                  SPICE_FRMTYP_DYN        5    Dynamic frame.
%
%                  SPICE_FRMTYP_SWTCH      6    Switch frame.
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
%      [idset] = cspice_kplfrm( frmcls, room )
%
%   returns:
%
%      idset    a SPICE set containing the ID codes of all reference frames
%               having specifications present in the kernel pool and
%               belonging to the specified class or classes.
%
%               [r,1] = size(idset); int32 = class(idset)
%
%               Note that if `frmcls' is set to SPICE_FRMTYP_INERTL, `idset'
%               will be empty on output.
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
%   1) Display the IDs and names of all reference frames having
%      specifications present in the kernel pool. Group the outputs
%      by frame class. Also fetch and display the entire set of IDs
%      and names using the parameter SPICE_FRMTYP_ALL.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: kplfrm_ex1.tm
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
%            File name            Contents
%            --------------       --------------------------
%            clem_v20.tf          Clementine FK
%            moon_060721.tf       Generic Lunar SPICE frames
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'clem_v20.tf'
%                                'moon_060721.tf' )
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function kplfrm_ex1()
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
%         META   = 'kplfrm_ex1.tm';
%         NFRAME = 1000;
%
%         %
%         % Load kernels that contain frame specifications.
%         %
%         cspice_furnsh( META );
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
%               [idset] = cspice_kplfrm( i, NFRAME );
%
%               outlin  = sprintf( 'Number of frames of class %d: %d', ...
%                                                  i, size( idset )(1) );
%
%            else
%
%               %
%               % Fetch IDs of all frames specified in the kernel pool.
%               %
%               [idset] = cspice_kplfrm( SPICE_FRMTYP_ALL, NFRAME );
%
%               outlin  = sprintf( [ 'Number of frames in the kernel', ...
%                                    ' pool: %d' ], size( idset )(1)   );
%
%            end
%
%            %
%            % Display the fetched frame IDs and corresponding names.
%            %
%            fprintf( '\n' )
%            fprintf( '%s\n', outlin )
%            fprintf( '   Frame IDs and names\n' )
%
%            for j=1:size( idset )(1)
%
%               [frname] = cspice_frmnam( idset(j) );
%
%               fprintf( '%12d   %s\n', idset(j), frname )
%
%            end
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
%      Number of frames of class 1: 0
%         Frame IDs and names
%
%      Number of frames of class 2: 1
%         Frame IDs and names
%             31002   MOON_PA_DE403
%
%      Number of frames of class 3: 1
%         Frame IDs and names
%            -40000   CLEM_SC_BUS
%
%      Number of frames of class 4: 11
%         Frame IDs and names
%            -40008   CLEM_CPT
%            -40007   CLEM_BSTAR
%            -40006   CLEM_ASTAR
%            -40005   CLEM_LIDAR
%            -40004   CLEM_LWIR
%            -40003   CLEM_NIR
%            -40002   CLEM_UVVIS
%            -40001   CLEM_HIRES
%             31000   MOON_PA
%             31001   MOON_ME
%             31003   MOON_ME_DE403
%
%      Number of frames of class 5: 0
%         Frame IDs and names
%
%      Number of frames of class 6: 0
%         Frame IDs and names
%
%      Number of frames in the kernel pool: 13
%         Frame IDs and names
%            -40008   CLEM_CPT
%            -40007   CLEM_BSTAR
%            -40006   CLEM_ASTAR
%            -40005   CLEM_LIDAR
%            -40004   CLEM_LWIR
%            -40003   CLEM_NIR
%            -40002   CLEM_UVVIS
%            -40001   CLEM_HIRES
%            -40000   CLEM_SC_BUS
%             31000   MOON_PA
%             31001   MOON_ME
%             31002   MOON_PA_DE403
%             31003   MOON_ME_DE403
%
%
%-Particulars
%
%   This routine enables SPICE-based applications to conveniently
%   find the frame ID codes of reference frames having specifications
%   present in the kernel pool. Such frame specifications are
%   introduced into the kernel pool either by loading frame kernels
%   or by means of calls to the kernel pool 'put' API routines
%
%      cspice_pcpool
%      cspice_pdpool
%      cspice_pipool
%
%   Given a reference frame's ID code, other attributes of the
%   frame can be obtained via calls to the Mice APIs
%
%      cspice_frmnam {Return a frame's name}
%      cspice_frinfo {Return a frame's center, class, and class ID}
%
%   This routine has a counterpart
%
%      cspice_bltfrm
%
%   which fetches the frame IDs of all built-in reference frames.
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
%   3)  Frames of class 1 may not be specified in the kernel pool.
%       However, for the convenience of users, this routine does not
%       signal an error if the input class is set to SPICE_FRMTYP_INERTL.
%       In this case the output set will be empty.
%
%   4)  This routine relies on the presence of just three kernel
%       variable assignments for a reference frame in order to
%       determine that that reference frame has been specified:
%
%          FRAME_<frame name>       = <ID code>
%          FRAME_<ID code>_NAME     = <frame name>
%
%       and either
%
%          FRAME_<ID code>_CLASS    = <class>
%
%       or
%
%          FRAME_<frame name>_CLASS = <class>
%
%       It is possible for the presence of an incomplete frame
%       specification to trick this routine into incorrectly
%       deciding that a frame has been specified. This routine
%       does not attempt to diagnose this problem.
%
%   5)  If any of the input arguments, `frmcls' or `room', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `frmcls' or `room', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   Reference frame specifications for frames that are not
%   built in are typically established by loading frame kernels.
%
%-Restrictions
%
%   1)  This routine will work correctly if the kernel pool
%       contains no invalid frame specifications. See the
%       description of exception 4 above. Users must ensure
%       that no invalid frame specifications are introduced
%       into the kernel pool, either by loaded kernels or
%       by means of the kernel pool "put" APIs.
%
%-Required_Reading
%
%   FRAMES.REQ
%   KERNEL.REQ
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
%   -Mice Version 1.0.0, 08-AUG-2021 (JDR)
%
%-Index_Entries
%
%   fetch IDs of reference_frames from the kernel_pool
%
%-&
function [idset] = cspice_kplfrm( frmcls, room )

   switch nargin
      case 2

         frmcls = zzmice_int(frmcls);
         room = zzmice_int(room, [1, int32(inf)] );

      otherwise

         error ( 'Usage: [idset] = cspice_kplfrm( frmcls, room )' )

   end

   %
   % Call the MEX library.
   %
   try
      [idset] = mice('kplfrm_c', frmcls, room);
   catch spiceerr
      rethrow(spiceerr)
   end
