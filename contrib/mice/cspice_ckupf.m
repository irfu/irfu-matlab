%-Abstract
%
%   CSPICE_CKUPF unloads a CK pointing file so that it will no longer be
%   searched by the readers.
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
%      handle   the integer handle assigned to the CK file upon loading.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_ckupf( handle )
%
%   returns:
%
%      None.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Unload a CK kernel specified by an integer handle, making
%      room to load another CK.
%
%         cspice_ckupf( handle );
%
%   2) Load a CASSINI CK file and obtain the position transformation
%      matrix from J2000 to the spacecraft reference frame provided
%      by the CK and its angular velocity vector at a given spacecraft
%      clock time.
%
%
%      Use the CK kernel below to load the CASSINI image navigated
%      spacecraft pointing and orientation data.
%
%         04153_04182ca_ISS.bc
%
%
%      In order to convert from spacecraft clock time to 'ticks,'
%      (units of encoded SCLK) as required by cspice_ckgpav, we will need
%      to load as well a CASSINI SCLK.
%
%      Use the SCLK kernel below to load the CASSINI spacecraft clock
%      time correlation data required for the conversion between
%      spacecraft clock string representation and double precision
%      encoding of spacecraft clock counts.
%
%         cas00071.tsc
%
%
%      Example code begins here.
%
%
%      function ckupf_ex2()
%
%         %
%         % Constants for this program:
%         %
%         % -- The code for the CASSINI spacecraft clock is -82.
%         %
%         % -- The code for CASSINI spacecraft reference frame is -82000.
%         %
%         % -- Tolerance: 1 second. It must be converted to 'ticks'
%         %    (units of encoded SCLK) for input to cspice_ckgpav.
%         %
%         % -- The reference frame we want is J2000.
%         %
%         SC   =   -82;
%         INST =   -82000;
%         REF  =   'J2000';
%         TOL  =   '1.0';
%
%         CK   =   '04153_04182ca_ISS.bc';
%         SCLK =   'cas00071.tsc';
%         sclkch = '1465644281.0';
%
%         %
%         % Load the CK for read access. This call may be replaced (as
%         % recommended by NAIF) by cspice_furnsh.
%         %
%         [handle] = cspice_cklpf( CK );
%
%         %
%         % We need to load a CASSINI SCLK kernel to convert from
%         % clock string to ticks.  Although not required for
%         % the CASSINI spacecraft clock, most modern spacecraft
%         % clocks require a leapseconds kernel to be loaded in
%         % addition to an SCLK kernel.
%         %
%         cspice_furnsh( SCLK );
%
%         %
%         % Convert tolerance from CASSINI formatted character string
%         % SCLK to ticks, which are units of encoded SCLK.
%         %
%         [toltik] = cspice_sctiks( SC, TOL );
%
%         %
%         % cspice_ckgpav requires encoded spacecraft clock time.
%         %
%         [sclkdp] = cspice_scencd( SC, sclkch );
%
%         [cmat, av, clkout, found] = cspice_ckgpav( INST,   sclkdp,       ...
%                                                    toltik, REF     );
%
%         %
%         % Display the results.
%         %
%         if ( found )
%
%            [clkch] = cspice_scdecd( SC, clkout );
%
%            fprintf( 'Requested SCLK time         : %s\n', sclkch )
%            fprintf( '   CASSINI SCLK time        : %s\n', clkch )
%            fprintf( '   J2000 to S/C frame matrix:\n' )
%            fprintf( '\n' )
%            for i=1:3
%
%               fprintf( '%20.10f %19.10f %19.10f\n', cmat(i,:) )
%
%            end
%            fprintf( '\n' )
%            fprintf( [ '   Angular velocity vector  : %10.7f %10.7f',     ...
%                       ' %10.7f\n' ], av             )
%
%         else
%
%               fprintf( 'Pointing not found for time %s\n', sclkch )
%
%         end
%
%         %
%         % Close the CK file. This call may be replaced (as
%         % recommended by NAIF) by cspice_unload, if cspice_furnsh has
%         % been used to load the file.
%         %
%         cspice_ckupf( handle );
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
%      Requested SCLK time         : 1465644281.0
%         CASSINI SCLK time        : 1/1465644281.171
%         J2000 to S/C frame matrix:
%
%             -0.3353514559        0.8643744402        0.3746948467
%             -0.9378874268       -0.3438519652       -0.0461844200
%              0.0889189272       -0.3669095980        0.9259971767
%
%         Angular velocity vector  :  0.0000000  0.0000000  0.0000000
%
%
%   3) The following example extracts the first 20 lines of the
%      comment area of a CK, displaying the comments on the terminal
%      screen.
%
%
%      Example code begins here.
%
%
%      function ckupf_ex3()
%
%         %
%         % Local parameters.
%         %
%         LINLEN =   1001;
%         BUFFSZ =   20;
%
%         %
%         % Local variables.
%         %
%
%         ckname = input( 'Enter name of CK > ', 's' );
%
%         %
%         % Open the CK for read access. This operation could have
%         % been done with cspice_dafopr.
%         %
%         [handle] = cspice_cklpf( ckname );
%
%         %
%         % Extract up to 20 lines from the comment area of the
%         % loaded CK file and display them on the terminal screen.
%         %
%         [buffer, done] = cspice_dafec( handle, BUFFSZ, LINLEN );
%
%         buffer = cellstr(buffer);
%         for i=1:numel(buffer)
%
%               fprintf( '%s\n', char(buffer(i)) )
%
%         end
%
%         %
%         % Close the CK file. This operation could have been done
%         % with cspice_dafcls.
%         %
%         cspice_ckupf( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the Cassini CK file named 04161_04164ra.bc as
%      input CK file, the output was:
%
%
%      Enter name of CK > 04161_04164ra.bc
%      \beginlabel
%      PDS_VERSION_ID               = PDS3
%      RECORD_TYPE                  = FIXED_LENGTH
%      RECORD_BYTES                 = 1024
%      ^SPICE_KERNEL                = "04161_04164ra.bc"
%      MISSION_NAME                 = "CASSINI-HUYGENS"
%      SPACECRAFT_NAME              = "CASSINI ORBITER"
%      DATA_SET_ID                  = "CO-S/J/E/V-SPICE-6-V1.0"
%      KERNEL_TYPE_ID               = CK
%      PRODUCT_ID                   = "04161_04164ra.bc"
%      PRODUCT_CREATION_TIME        = 2005-06-29T21:28:09
%      PRODUCER_ID                  = "CASSINI_AACS/JPL"
%      MISSION_PHASE_NAME           = "SCIENCE CRUISE"
%      PRODUCT_VERSION_TYPE         = ACTUAL
%      PLATFORM_OR_MOUNTING_NAME    = "N/A"
%      START_TIME                   = 2004-06-09T12:00:03.631
%      STOP_TIME                    = 2004-06-12T11:58:57.943
%      SPACECRAFT_CLOCK_START_COUNT = "1/1465475046.160"
%      SPACECRAFT_CLOCK_STOP_COUNT  = "1/1465734182.160"
%      TARGET_NAME                  = "N/A"
%
%
%-Particulars
%
%   Unloading a file with cspice_ckupf removes that file from consideration
%   by the CK readers. In doing so, it frees up space for another
%   file to be loaded.
%
%-Exceptions
%
%   1)  Unloading a file that has not been loaded is a no-op.
%       No error is signaled.
%
%   2)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `handle' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   The file referred to by `handle' is unloaded.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   CK.REQ
%   DAF.REQ
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
%   -Mice Version 1.0.0, 30-JUN-2021 (JDR)
%
%-Index_Entries
%
%   unload CK pointing file
%
%-&
function cspice_ckupf( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_ckupf( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('ckupf_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end
