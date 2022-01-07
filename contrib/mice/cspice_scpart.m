%-Abstract
%
%   CSPICE_SCPART gets spacecraft clock partition information from a
%   spacecraft clock kernel file.
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
%      sc       the NAIF ID for the spacecraft whose clock partition
%               information is being requested.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%   the call:
%
%      [pstart, pstop] = cspice_scpart( sc )
%
%   returns:
%
%      pstart   an array containing nparts partition start times represented
%               as double precision, encoded SCLK ("ticks").
%
%               [nparts,1] = size(pstart); double = class(pstart)
%
%               The values contained in `pstart' are whole numbers.
%
%      pstop    an array containing nparts partition end times represented as
%               double precision, encoded SCLK ("ticks").
%
%               [nparts,1] = size(pstop); double = class(pstop)
%
%               The values contained in `pstop' are whole numbers.
%
%-Parameters
%
%   MXPART      is the maximum number of spacecraft clock partitions
%               expected in the kernel file for any one spacecraft.
%               MXPART is currently set to 9999.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) The following code example finds partition start and stop
%      times for the Stardust spacecraft from a spacecraft clock
%      kernel file. Since those times are always returned in units
%      of ticks, the program uses cspice_scfmt to print the times in
%      Stardust clock format.
%
%      Use the SCLK kernel below to load the Stardust time
%      correlation data and spacecraft clock partition information.
%
%         sdu_sclkscet_00074.tsc
%
%
%      Example code begins here.
%
%
%      function scpart_ex1()
%
%         %
%         % Assign the value for the Stardust spacecraft ID.
%         %
%         sc = -29;
%
%         %
%         % Load the SCLK file.
%         %
%         cspice_furnsh( 'sdu_sclkscet_00074.tsc' );
%
%         %
%         % Retrieve the arrays for `pstart' and `pstop' and the
%         % number of partitions within the SCLK.
%         %
%         [pstart, pstop] = cspice_scpart( sc );
%
%         %
%         % Loop over each array value.
%         %
%         nparts = size(pstart,1);
%         for i=1:nparts
%
%            [start] = cspice_scfmt( sc, pstart(i) );
%            [stop]  = cspice_scfmt( sc, pstop(i) );
%
%            fprintf( '\n' )
%            fprintf( 'Partition: %d\n', i )
%            fprintf( '   Start : %s\n', start )
%            fprintf( '   Stop  : %s\n', stop )
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
%      Partition: 1
%         Start : 0000000000.000
%         Stop  : 0602741011.080
%
%      Partition: 2
%         Start : 0602741014.217
%         Stop  : 0605660648.173
%
%      Partition: 3
%         Start : 0605660649.000
%         Stop  : 0631375256.224
%
%      Partition: 4
%         Start : 0631375257.000
%         Stop  : 0633545577.218
%
%      Partition: 5
%         Start : 0633545578.000
%         Stop  : 0644853954.043
%
%      Partition: 6
%         Start : 0644853954.000
%         Stop  : 0655316480.089
%
%      Partition: 7
%         Start : 0655316480.000
%         Stop  : 0660405279.066
%
%      Partition: 8
%         Start : 0660405279.000
%         Stop  : 0670256568.229
%
%      Partition: 9
%         Start : 0670256569.000
%         Stop  : 0674564039.091
%
%      Partition: 10
%         Start : 0674564040.000
%         Stop  : 4294537252.255
%
%
%-Particulars
%
%   cspice_scpart looks for two variables in the kernel pool for each
%   spacecraft's partition information. If sc = -nn, then the names of
%   the variables are
%
%      SCLK_PARTITION_START_nn
%      SCLK_PARTITION_END_nn
%
%   The start and stop times returned are in units of "ticks".
%
%-Exceptions
%
%   1)  If the kernel variables containing the spacecraft clock
%       partition start and stop times have not been loaded in the
%       kernel pool, an error is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the number of start and stop times are different, the error
%       SPICE(NUMPARTSUNEQUAL) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the input argument `sc' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   4)  If the input argument `sc' is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   An SCLK kernel containing spacecraft clock partition start
%   and stop times for the spacecraft clock indicated by `sc' must
%   be loaded into the kernel pool.
%
%-Restrictions
%
%   1)  This routine assumes that an SCLK kernel appropriate to the
%       spacecraft identified by `sc' has been loaded into the kernel
%       pool.
%
%-Required_Reading
%
%   MICE.REQ
%   SCLK.REQ
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
%   -Mice Version 1.0.0, 26-NOV-2021 (JDR)
%
%-Index_Entries
%
%   spacecraft_clock partition information
%
%-&
function [pstart, pstop] = cspice_scpart( sc )

   switch nargin
      case 1

         sc = zzmice_int(sc);

      otherwise

         error ( [ 'Usage: [pstart(nparts), pstop(nparts)] = '              ...
                   'cspice_scpart( sc )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [pstart, pstop] = mice('scpart_c', sc);
   catch spiceerr
      rethrow(spiceerr)
   end
