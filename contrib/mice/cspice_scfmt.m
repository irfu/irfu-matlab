%-Abstract
%
%   CSPICE_SCFMT converts encoded spacecraft clock ticks to character clock
%   format.
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
%      sc       the NAIF ID number for the spacecraft whose clock's time is
%               being decoded.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%      ticks    the double precision encoding of a clock time in units of
%               ticks.
%
%               [1,1] = size(ticks); double = class(ticks)
%
%               Partition information is not reflected in this value.
%
%               An analogy may be drawn between a spacecraft clock and
%               a standard wall clock. The number of ticks
%               corresponding to the wall clock string
%
%                  hh:mm:ss
%
%               would be the number of seconds represented by that
%               time.
%
%               For example,
%
%                  Clock string    Number of ticks
%                  ------------    ---------------
%                    00:00:10             10
%                    00:01:00             60
%                    00:10:00            600
%                    01:00:00           3600
%                    01:01:00           3660
%
%               If `ticks' contains a fractional part the result is the
%               same as if `ticks' had been rounded to the nearest whole
%               number.
%
%               See the -Examples section below for examples of
%               actual spacecraft clock conversions.
%
%   the call:
%
%      [clkstr] = cspice_scfmt( sc, ticks )
%
%   returns:
%
%      clkstr   the spacecraft clock character string corresponding to
%               `ticks'.
%
%               [1,c1] = size(clkstr); char = class(clkstr)
%
%                  or
%
%               [1,1] = size(clkstr); cell = class(clkstr)
%
%               Partition information is not included in `clkstr'.
%
%               Using Galileo as an example, the full format clock
%               string is
%
%                  wwwwwwww:xx:y:z
%
%               where z is a mod-8 counter (values 0-7) which
%               increments approximately once every 8 1/3 ms., y is a
%               mod-10 counter (values 0-9) which increments once
%               every time z turns over, i.e., approximately once every
%               66 2/3 ms., xx is a mod-91 (values 0-90) counter
%               which increments once every time y turns over, i.e.,
%               once every 2/3 seconds. wwwwwwww is the Real-Time Image
%               Count (RIM), which increments once every time xx turns
%               over, i.e., once every 60 2/3 seconds. The roll-over
%               expression for the RIM is 16777215, which corresponds
%               to approximately 32 years.
%
%               wwwwwwww, xx, y, and z are referred to interchangeably
%               as the fields or components of the spacecraft clock.
%               SCLK components may be separated by any of these five
%               characters: ' '  ':'  ','  '-'  '.'
%               The delimiter used is determined by a kernel pool
%               variable and can be adjusted by the user.
%
%               Some spacecraft clock components have offset, or
%               starting, values different from zero. For example,
%               with an offset value of 1, a mod 20 counter would
%               cycle from 1 to 20 instead of from 0 to 19.
%
%               See the SCLK required reading for a detailed
%               description of the Voyager and Mars Observer clock
%               formats.
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
%      function scfmt_ex1()
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
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
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
%   2) Below are some examples illustrating various input numbers of
%      ticks and the resulting clock string outputs for the Galileo
%      spacecraft.
%
%         ticks                clkstr
%         ----------------     --------------------
%         -1                   Error: Ticks must be a positive number
%         0                    '00000000:00:0:0'
%         1                    '00000000:00:0:1'
%         1.3                  '00000000:00:0:1'
%         1.5                  '00000000:00:0:2'
%         2                    '00000000:00:0:2'
%         7                    '00000000:00:0:7'
%         8                    '00000000:00:1:0'
%         80                   '00000000:01:0:0'
%         88                   '00000000:01:1:0'
%         7279                 '00000000:90:9:7'
%         7280                 '00000001:00:0:0'
%         1234567890           '00169583:45:6:2'
%
%      The following examples are for the Voyager 2 spacecraft.
%      Note that the third component of the Voyager clock has an
%      offset value of one.
%
%         ticks                clkstr
%         ----------------     --------------------
%         -1                   Error: Ticks must be a positive number
%         0                    '00000:00:001'
%         1                    '00000:00:002'
%         1.3                  '00000:00:002'
%         1.5                  '00000:00:003'
%         2                    '00000:00:003'
%         799                  '00000:00:800'
%         800                  '00000:01:001'
%         47999                '00000:59:800'
%         48000                '00001:00:001'
%         3145727999           '65535:59:800'
%
%-Particulars
%
%   The routine cspice_sctiks performs the inverse operation to cspice_scfmt,
%   converting from clock format to number of ticks.
%
%   Note the important difference between cspice_scfmt and cspice_scdecd.
%   cspice_scdecd converts some number of ticks since the spacecraft clock
%   start time to a character string which includes a partition number.
%   cspice_scfmt, which is called by cspice_scdecd, does not make use of
%   partition information.
%
%-Exceptions
%
%   1)  If the data type for the spacecraft is not supported, the
%       error SPICE(NOTSUPPORTED) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the value for `ticks' is negative, an error is signaled
%       by a routine in the call tree of this routine.
%
%   3)  If the SCLK kernel file does not contain data for the
%       spacecraft specified by `sc', an error is signaled by a routine
%       in the call tree of this routine.
%
%   4)  If any of the input arguments, `sc' or `ticks', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `sc' or `ticks', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
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
%   -Mice Version 1.0.0, 07-SEP-2020 (JDR)
%
%-Index_Entries
%
%   convert spacecraft_clock ticks to character clock format
%
%-&
function [clkstr] = cspice_scfmt( sc, ticks )

   switch nargin
      case 2

         sc    = zzmice_int(sc);
         ticks = zzmice_dp(ticks);

      otherwise

         error ( 'Usage: [`clkstr`] = cspice_scfmt( sc, ticks )' )

   end

   %
   % Call the MEX library.
   %
   try
      [clkstr] = mice('scfmt_c', sc, ticks);
   catch spiceerr
      rethrow(spiceerr)
   end
