%-Abstract
%
%   CSPICE_SCDECD converts a double precision encoding of spacecraft clock
%   time into a character representation.
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
%      sc       the NAIF integer code of the spacecraft whose clock's time is
%               being decoded.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%      sclkdp   the double precision encoding(s) of a clock time(s) in units
%               of ticks since the spacecraft clock start time.
%
%               [1,n] = size(sclkdp); double = class(sclkdp)
%
%               This value does reflect partition information.
%
%               An analogy may be drawn between a spacecraft clock
%               and a standard wall clock. The number of ticks
%               corresponding to the wall clock string
%
%                  hh:mm:ss
%
%               would be the number of seconds represented by that
%               time.
%
%               For example:
%
%                  Clock string      Number of ticks
%                  ------------      ---------------
%                    00:00:10              10
%                    00:01:00              60
%                    00:10:00             600
%                    01:00:00            3600
%
%               If `sclkdp' contains a fractional part the result
%               is the same as if `sclkdp' had been rounded to the
%               nearest whole number.
%
%   the call:
%
%      [sclkch] = cspice_scdecd( sc, sclkdp )
%
%   returns:
%
%      sclkch   the character representation(s) of the clock count(s).
%
%               [n,c1] = size(sclkch); char = class(sclkch)
%
%               The exact form that `sclkch' takes depends on the
%               spacecraft.
%
%               Nevertheless, `sclkch' will have the following general
%               format:
%
%                  'pp/sclk_string'
%
%               'pp' is an integer greater than or equal to one and
%               represents a "partition number".
%
%               Each mission is divided into some number of partitions.
%               A new partition starts when the spacecraft clock
%               resets, either to zero, or to some other
%               value. Thus, the first partition for any mission
%               starts with launch, and ends with the first clock
%               reset. The second partition starts immediately when
%               the first stopped, and so on.
%
%               In order to be completely unambiguous about a
%               particular time, you need to specify a partition number
%               along with the standard clock string.
%
%               Information about when partitions occur for different
%               missions is contained in a spacecraft clock kernel
%               file which needs to be loaded into the kernel pool
%               before calling cspice_scdecd.
%
%               The routine cspice_scpart may be used to read the partition
%               start and stop times, in encoded units of ticks, from
%               the kernel file.
%
%               Since the end time of one partition is coincident with
%               the begin time of the next, two different time strings
%               with different partition numbers can encode into the
%               same value.
%
%               For example, if partition 1 ends at time t1, and
%               partition 2 starts at time t2, then
%
%                  '1/t1' and '2/t2'
%
%               will be encoded into the same value, say `x'. cspice_scdecd
%               always decodes such values into the latter of the
%               two partitions. In this example,
%
%                  [clkstr] = cspice_scdecd( x, sc );
%
%               will result in
%
%                  clkstr = '2/t2'.
%
%               'sclk_string' is a spacecraft specific clock string,
%               typically consisting of a number of components
%               separated by delimiters.
%
%               Using Galileo as an example, the full format is
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
%               `sclkch' returns with the same vectorization measure, N,
%               as `sclkdp'.
%
%-Parameters
%
%   MXPART      is the maximum number of spacecraft clock partitions
%               expected in the kernel file for any one spacecraft.
%               MXPART is currently set to 9999.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Double precision encodings of spacecraft clock counts are used
%      to tag pointing data in the C-kernel.
%
%      In the following example, pointing for a sequence of images
%      from the CASSINI Imaging Science Subsystem (ISS) is requested
%      from the C-kernel using an array of character spacecraft clock
%      counts as input. The clock counts attached to the output are
%      then decoded to character and compared with the input strings.
%
%      Use the CK kernel below to load the CASSINI image navigated
%      spacecraft pointing and orientation data.
%
%         04153_04182ca_ISS.bc
%
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
%      function scdecd_ex1()
%
%         %
%         % Local parameters.
%         %
%         % The instrument we want pointing for is the CASSINI
%         % spacecraft. The reference frame we want is
%         % J2000. The spacecraft is CASSINI.
%         %
%         SC     =   -82;
%         INST   =   -82000;
%         REF    =   'J2000';
%         CK     =   '04153_04182ca_ISS.bc';
%         SCLK   =   'cas00071.tsc';
%         NPICS  =   4;
%         CLKTOL =   '1.0';
%
%         %
%         % Set the input SCLK strings.
%         %
%         sclkin = { '1/1465644279.0', '1/1465644281.0',                   ...
%                    '1/1465644351.0', '1/1465644361.0' };
%
%         %
%         % Load the appropriate files. We need
%         %
%         %    1. CK file containing pointing data.
%         %    2. Spacecraft clock kernel file.
%         %
%         cspice_furnsh( CK   );
%         cspice_furnsh( SCLK );
%
%         %
%         % Convert the tolerance string to ticks.
%         %
%         [tol] = cspice_sctiks( SC, CLKTOL );
%
%         for i=1:NPICS
%
%            [timein]  = cspice_scencd( SC, sclkin(i) );
%
%            [cmat, timeout, found] = cspice_ckgp( INST, timein, tol, REF );
%
%            fprintf( '\n' )
%            fprintf( 'Input s/c clock count : %s\n', char(sclkin(i)) )
%
%            if ( found )
%
%               [sclkout] = cspice_scdecd( SC, timeout );
%
%               fprintf( 'Output s/c clock count: %s\n', sclkout )
%               fprintf( 'Output C-Matrix:\n' )
%               fprintf( '%21.15f %20.15f %20.15f\n', cmat' );
%
%            else
%
%               fprintf( 'No pointing found.\n' )
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
%      Input s/c clock count : 1/1465644279.0
%      No pointing found.
%
%      Input s/c clock count : 1/1465644281.0
%      Output s/c clock count: 1/1465644281.171
%      Output C-Matrix:
%         -0.335351455948710    0.864374440205611    0.374694846658341
%         -0.937887426812980   -0.343851965210223   -0.046184419961653
%          0.088918927227039   -0.366909598048763    0.925997176691424
%
%      Input s/c clock count : 1/1465644351.0
%      Output s/c clock count: 1/1465644351.071
%      Output C-Matrix:
%         -0.335380929397586    0.864363638262230    0.374693385378623
%         -0.937874292008090   -0.343889838107825   -0.046169163264003
%          0.088946301703530   -0.366899550417080    0.925998528787713
%
%      Input s/c clock count : 1/1465644361.0
%      No pointing found.
%
%
%   2) Convert a series of CASSINI encoded SCLK values to their
%      corresponding character representation of spacecraft clock
%      time, and convert back that SCLK string to double precision
%      form.
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
%      function scdecd_ex2()
%
%         %
%         % Assign values for the spacecraft ID (CASSINI),
%         % the SCLK kernel, and a double precision
%         % encodings of SCLK strings
%         %
%         SC     = -82;
%         SCLK   = 'cas00071.tsc';
%         timein = 197483587237.0;
%
%         %
%         % Load the kernel file.
%         %
%         cspice_furnsh( SCLK )
%
%         %
%         % Convert the CASSINI encoded SCLK to an
%         % SCLK string.
%         %
%         sclkch = cspice_scdecd( SC, timein );
%
%         %
%         % Convert the SCLK string to double precision form.
%         % The output value should match the original.
%         %
%         sclkdp = cspice_scencd( SC, sclkch );
%
%         disp( 'Scalar:' )
%
%         txt = sprintf( 'Original: %20.8f', timein );
%         disp( txt )
%
%         txt = sprintf( 'SCLKCH  : %s',     sclkch );
%         disp( txt )
%
%         txt = sprintf( 'Decoded : %20.8f', sclkdp );
%         disp( txt )
%
%         disp( ' ' )
%
%         %
%         % Convert a vector of SCLK values.
%         %
%         timein = [ 197483587237.0,                                       ...
%                    197483587250.0,                                       ...
%                    197485901583.201,                                     ...
%                    197486447183.0,                                       ...
%                    198136032015.400 ];
%
%         %
%         % Convert the SCLK double precision values to the string
%         % representation, then convert to the dp form. As before, the
%         % output value should match the original.
%         %
%         sclkch = cspice_scdecd( SC, timein );
%         sclkdp = cspice_scencd( SC, sclkch );
%
%         disp( 'Vector:' )
%         for i=1:5
%
%            txt = sprintf( 'Original: %20.8f', timein(i)   );
%            disp( txt )
%
%            txt = sprintf( 'SCLKCH  : %s',     sclkch(i,:) );
%            disp( txt )
%
%            txt = sprintf( 'Decoded : %20.8f', sclkdp(i)   );
%            disp( txt )
%
%            disp( ' ' )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Scalar:
%      Original: 197483587237.00000000
%      SCLKCH  : 1/1465644281.165
%      Decoded : 197483587237.00000000
%
%      Vector:
%      Original: 197483587237.00000000
%      SCLKCH  : 1/1465644281.165
%      Decoded : 197483587237.00000000
%
%      Original: 197483587250.00000000
%      SCLKCH  : 1/1465644281.178
%      Decoded : 197483587250.00000000
%
%      Original: 197485901583.20098877
%      SCLKCH  : 1/1465653322.015
%      Decoded : 197485901583.00000000
%
%      Original: 197486447183.00000000
%      SCLKCH  : 1/1465655453.079
%      Decoded : 197486447183.00000000
%
%      Original: 198136032015.39999390
%      SCLKCH  : 1/1468192894.015
%      Decoded : 198136032015.00000000
%
%
%-Particulars
%
%   In general, it is difficult to compare spacecraft clock counts
%   numerically since there are too many clock components for a
%   single comparison. The routine cspice_scencd provides a method of
%   assigning a single double precision number to a spacecraft's
%   clock count, given one of its character representations.
%
%   This routine performs the inverse operation to cspice_scencd, converting
%   an encoded double precision number to character format.
%
%   To convert the number of ticks since the start of the mission to
%   a clock format character string, cspice_scdecd:
%
%      1) Determines the spacecraft clock partition that TICKS falls
%         in.
%
%      2) Subtracts off the number of ticks occurring in previous
%         partitions, to get the number of ticks since the beginning
%         of the current partition.
%
%      3) Converts the resulting ticks to clock format and forms the
%         string
%
%            'partition_number/clock_string'
%
%-Exceptions
%
%   1)  If kernel variables required by this routine are unavailable,
%       an error is signaled by a routine in the call tree of this
%       routine. `sclkch' will be returned as an empty string in this
%       case.
%
%   2)  If the number of partitions in the kernel file for spacecraft
%       `sc' exceeds the parameter MXPART, the error SPICE(TOOMANYPARTS)
%       is signaled by a routine in the call tree of this routine.
%       `sclkch' will be returned as an empty string in this case.
%
%   3)  If the encoded value does not fall in the boundaries of the
%       mission, the error SPICE(VALUEOUTOFRANGE) is signaled by a
%       routine in the call tree of this routine. `sclkch' will be
%       returned as an empty string in this case.
%
%   4)  If any of the input arguments, `sc' or `sclkdp', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `sc' or `sclkdp', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   A kernel file containing spacecraft clock partition information
%   for the desired spacecraft must be loaded, using the routine
%   cspice_furnsh, before calling this routine.
%
%-Restrictions
%
%   1)  Assumes that an SCLK kernel file appropriate for the clock
%       designated by `sc' is loaded in the kernel pool at the time
%       this routine is called.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Extended -I/O
%       section to provide further description of arguments. Added
%       complete examples; second one based on existing fragment using
%       CASSINI PDS archived data.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 06-JAN-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 18-APR-2006 (EDW)
%
%-Index_Entries
%
%   decode spacecraft_clock
%
%-&

function [sclkch] = cspice_scdecd(sc, sclkdp)

   switch nargin
      case 2

        sc     = zzmice_int(sc);
        sclkdp = zzmice_dp(sclkdp);

      otherwise

         error ( 'Usage: [_`sclkch`_] = cspice_scdecd(sc, _sclkdp_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [sclkch] = mice('scdecd_c',sc,sclkdp);
   catch spiceerr
      rethrow(spiceerr)
   end
