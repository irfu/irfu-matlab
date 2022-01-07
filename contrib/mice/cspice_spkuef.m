%-Abstract
%
%   CSPICE_SPKUEF unloads an ephemeris file so that it will no longer be
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
%      handle   Integer handle assigned to the file upon loading.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_spkuef( handle )
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
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Load a planetary ephemeris SPK; then look up a series of
%      geometric states of the Earth relative to the solar system
%      barycenter, referenced to the J2000 frame.
%
%      Use the SPK kernel below to load the required ephemerides
%      for the Earth and the Earth Barycenter relative to the solar
%      system barycenter.
%
%         de405.bsp
%
%
%      Example code begins here.
%
%
%      function spkuef_ex1()
%
%         %
%         % Local constants.
%         %
%         MAXITR   = 5;
%         ET0      = -315576000.0;
%         STEP     = 3600.0;
%
%         ABCORR   = 'NONE';
%         FRAME    = 'J2000';
%         OBSERVER = 'SOLAR SYSTEM BARYCENTER';
%         SPK      = 'de405.bsp';
%         TARGET   = 'EARTH';
%
%         %
%         % Load the SPK file.
%         %
%         [handle] = cspice_spklef( SPK );
%
%         %
%         % Step through a series of epochs, looking up a state vector
%         % at each one.
%         %
%         for i=0:MAXITR-1
%
%            et =  ET0 + i*STEP;
%
%            [state, lt] = cspice_spkezr( TARGET,  et,     ...
%                                         FRAME,   ABCORR, ...
%                                         OBSERVER         );
%
%            fprintf( '\n' )
%            fprintf( 'et = %20.10f\n', et )
%            fprintf( '\n' )
%            fprintf( 'J2000 x-position (km):   %20.10f\n', state(1) )
%            fprintf( 'J2000 y-position (km):   %20.10f\n', state(2) )
%            fprintf( 'J2000 z-position (km):   %20.10f\n', state(3) )
%            fprintf( 'J2000 x-velocity (km/s): %20.10f\n', state(4) )
%            fprintf( 'J2000 y-velocity (km/s): %20.10f\n', state(5) )
%            fprintf( 'J2000 z-velocity (km/s): %20.10f\n', state(6) )
%            fprintf( '\n' )
%
%         end
%
%         %
%         % Unload the SPK kernel. This isn't necessary in a stand-
%         % alone program, but it's good practice in functions
%         % because it frees program and system resources.
%         %
%         cspice_spkuef( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      et = -315576000.0000000000
%
%      J2000 x-position (km):   -26772058.9514643848
%      J2000 y-position (km):   132760135.1677220613
%      J2000 z-position (km):    57557579.2735445350
%      J2000 x-velocity (km/s):       -29.7772753957
%      J2000 y-velocity (km/s):        -5.0656884328
%      J2000 z-velocity (km/s):        -2.1979102802
%
%
%      et = -315572400.0000000000
%
%      J2000 x-position (km):   -26879249.7439419106
%      J2000 y-position (km):   132741862.7243705541
%      J2000 z-position (km):    57549651.2066062242
%      J2000 x-velocity (km/s):       -29.7731620671
%      J2000 y-velocity (km/s):        -5.0856683968
%      J2000 z-velocity (km/s):        -2.2065710777
%
%
%      et = -315568800.0000000000
%
%      J2000 x-position (km):   -26986425.6981768459
%      J2000 y-position (km):   132723518.3595090210
%      J2000 z-position (km):    57541691.9637668282
%      J2000 x-velocity (km/s):       -29.7690319295
%      J2000 y-velocity (km/s):        -5.1056448242
%      J2000 z-velocity (km/s):        -2.2152302239
%
%
%      et = -315565200.0000000000
%
%      J2000 x-position (km):   -27093586.7536762133
%      J2000 y-position (km):   132705102.0859030634
%      J2000 z-position (km):    57533701.5509854183
%      J2000 x-velocity (km/s):       -29.7648849936
%      J2000 y-velocity (km/s):        -5.1256176961
%      J2000 z-velocity (km/s):        -2.2238877108
%
%
%      et = -315561600.0000000000
%
%      J2000 x-position (km):   -27200732.8499865979
%      J2000 y-position (km):   132686613.9163857996
%      J2000 z-position (km):    57525679.9742503539
%      J2000 x-velocity (km/s):       -29.7607212708
%      J2000 y-velocity (km/s):        -5.1455869940
%      J2000 z-velocity (km/s):        -2.2325435301
%
%
%-Particulars
%
%   A file is removed from consideration by the readers by a call to
%   cspice_spkuef.
%
%   The file table entry corresponding to the file referenced by
%   handle, is removed. Also any segment table entry which came from
%   the specified file is also deleted.
%
%   If the file specified by handle does not appear in the file table,
%   nothing happens.
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
%   MICE.REQ
%   SPK.REQ
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
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   unload SPK ephemeris file
%
%-&
function cspice_spkuef( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_spkuef( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('spkuef_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end
