%-Abstract
%
%   CSPICE_SCTIKS converts a spacecraft clock format string to
%   number of 'ticks'.
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
%      sc       the NAIF ID of the spacecraft clock whose time is
%               to encode.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%      clkstr   the scalar string or N-vector representation of the
%               `sc' spacecraft's clock time, WITHOUT PARTITION NUMBER.
%
%               [n,c1] = size(clkstr); char = class(clkstr)
%
%                  or
%
%               [1,n] = size(clkstr); cell = class(clkstr)
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
%               SCLK components may be separated by any of the
%               following characters: ' '  '.'  ':'  ','  '-'
%               Any number of spaces may separate the components and
%               the delimiters. The presence of the RIM component
%               is required. Successive components may be omitted, and
%               in such cases are assumed to represent zero values.
%
%               Values for the individual components may exceed the
%               maximum expected values. For instance, '0:0:0:9' is
%               an acceptable Galileo clock string, and will convert
%               to the same number of ticks as '0:0:1:1'.
%
%               Consecutive delimiters containing no intervening digits
%               are treated as if they delimit zero components.
%
%               Trailing zeros should always be included to match the
%               length of the counter. For example, a Galileo clock
%               count of '25684.90' should not be represented as
%               '25684.9'.
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
%   the call:
%
%      ticks = cspice_sctiks( sc, clkstr )
%
%   returns:
%
%      ticks    the tick values(s) represented by the spacecraft clock
%               string `clkstr'.
%
%               [1,n] = size(ticks); double = class(ticks)
%
%               A tick is defined to be the smallest time increment
%               expressible by the spacecraft clock.
%
%               An analogy may be drawn between a spacecraft clock
%               and a standard wall clock, measuring hours, minutes
%               and seconds. The number of ticks represented by the
%               wall clock string
%
%                  hh:mm:ss
%
%               would be the number of seconds represented by that
%               time.
%
%               For example:
%
%                  00:00:10  would convert to 10
%                  00:01:00  would convert to 60
%                  00:10:00  would convert to 600
%                  01:00:00  would convert to 3600
%                  01:01:00  would convert to 3660
%
%               See the -Examples section below for examples for
%               actual spacecraft clocks.
%
%               `ticks' returns with the same vectorization measure, N,
%               as `clkstr'.
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
%   1) cspice_sctiks is used as part of the process of encoding spacecraft
%      clock by cspice_scencd, though cspice_sctiks does not process any
%      partition information.
%
%      Another use of cspice_sctiks, however, is to convert a clock
%      measurement to ticks for use as a tolerance for the CK reader
%      cspice_ckgp.
%
%      In the following example, pointing for a sequence of images from
%      the CASSINI Imaging Science Subsystem (ISS) is requested from
%      the C-kernel using an array of character spacecraft clock counts
%      as input. The clock counts attached to the output are then
%      decoded to character and compared with the input strings.
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: sctiks_ex1.tm
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
%            File name              Contents
%            --------------------   -----------------------------------
%            04153_04182ca_ISS.bc   CASSINI image navigated spacecraft
%                                   orientation.
%            cas00071.tsc           CASSINI SCLK.
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( '04153_04182ca_ISS.bc',
%                                'cas00071.tsc'        )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function sctiks_ex1()
%
%         %
%         % The instrument we want pointing for is the CASSINI
%         % spacecraft.  The reference frame we want is
%         % J2000. The spacecraft is CASSINI.
%         %
%         SC   = -82;
%         INST = -82000;
%         REF = 'J2000';
%         META = 'sctiks_ex1.tm';
%
%         %
%         % Load the appropriate files. We need
%         %
%         %    1. CK file containing pointing data.
%         %    2. Spacecraft clock kernel file.
%         %
%         cspice_furnsh( META );
%
%         %
%         % The SCLK string includes a partition number. Pictures are
%         % never shuttered at intervals smaller than 1.0 seconds
%         % from each other.
%         %
%         clkstr = '1/1465644281.0';
%         tolstr = '1.0';
%
%         %
%         % Encode the clock string and the tolerance.
%         %
%         [timein] = cspice_scencd( SC, clkstr );
%         [tol] = cspice_sctiks( SC, tolstr );
%
%         %
%         % Get the pointing from the C-kernel.
%         %
%         [cmat, timeout, found] = cspice_ckgp( INST, timein, tol, REF );
%
%         if ( found )
%            [sclkout] = cspice_scdecd( SC, timeout );
%
%            fprintf( '\n' )
%            fprintf( 'Input  s/c clock count:  %s\n', clkstr  )
%            fprintf( 'Output s/c clock count:  %s\n', sclkout )
%            fprintf( 'Tolerance             :  %13.2f\n', tol )
%            fprintf( 'Output C-Matrix:          \n' )
%            fprintf( '    %18.15f   %18.15f   %18.15f\n', ...
%                              cmat(1,1), cmat(2,1), cmat(3,1) )
%            fprintf( '    %18.15f   %18.15f   %18.15f\n', ...
%                              cmat(1,2), cmat(2,2), cmat(3,2) )
%            fprintf( '    %18.15f   %18.15f   %18.15f\n', ...
%                              cmat(1,3), cmat(2,3), cmat(3,3) )
%         else
%            fprintf( '\n' )
%            fprintf( 'Input  s/c clock count:  %s\n', clkstr  )
%            fprintf( 'Tolerance             :  %13.2f\n', tol )
%            fprintf( 'No pointing found.\n' )
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
%      Input  s/c clock count:  1/1465644281.0
%      Output s/c clock count:  1/1465644281.171
%      Tolerance             :         256.00
%      Output C-Matrix:
%          -0.335351455948710   -0.937887426812980    0.088918927227039
%           0.864374440205611   -0.343851965210223   -0.366909598048763
%           0.374694846658341   -0.046184419961653    0.925997176691424
%
%
%   2) Below are some examples illustrating various clock string inputs
%      and the resulting outputs for the Galileo spacecraft. See the
%      SCLK required reading for a detailed description of the Galileo
%      clock format.
%
%         CLKSTR                TICKS
%         ----------------      --------------------
%         '0:0:0:1'             1
%         '0:0:1'               8
%         '0:1'                 80
%         '1'                   7280
%         '1 0 0 0'             7280
%         '1,0,0,0'             7280
%         '1:90'                14480
%         '1:9'                 8000
%         '1:09'                8000
%         '0-0-10'              80   |--  Third component is supposed
%         '0-1-0'               80   |    to be a mod-10 count.
%         '0/1/0'               Error: '/' is not an accepted delimiter.
%         '1: 00 : 0 : 1'       7281
%         '1:::1'               7281
%         '1.1.1.1.1'           Error: Too many components
%         '1.1.1.1.'            Error: The last delimiter signals that
%                                      a fifth component will follow.
%
%
%      The following examples are for the Voyager 2 spacecraft. Note
%      that the last component of the Voyager clock has an offset
%      value of 1.
%
%         CLKSTR                TICKS
%         ----------------      --------------------
%          '0.0.001'              0
%          '0:0:002'              1
%          '0:01'                 800
%          '1'                    48000
%          '1.0'                  48000
%          '1.0.0'                Error: The 3rd component is never 0.
%          '0.0:100'              99
%          '0-60-1'               48000
%          '1-1-1'                48800
%          '1-1-2'                48801
%
%-Particulars
%
%   Each spacecraft is assigned a clock type code in the kernel file.
%   cspice_sctiks calls the function SCTYPE to determine this value. If the
%   clock type is supported by SPICE, then the SPICELIB routine TIKSnn
%   is called to handle the actual conversion from clock format to
%   number of ticks. The nn in TIKSnn refers to the spacecraft clock
%   type code. Different spacecraft have distinct clock formats but
%   can still be of the same clock type.
%
%   The TIKSnn routines are routines to the routines SCLKnn, which
%   also contain the ticks-to-clock format conversion routines FMTnn.
%   FMTnn is called by the routine cspice_scfmt, which performs the
%   inverse operation to cspice_sctiks.
%
%   Note the important difference between cspice_scencd and cspice_sctiks.
%   cspice_scencd converts a clock string to the number of ticks it
%   represents since the beginning of the mission, and so uses partition
%   information. cspice_sctiks just converts to absolute ticks.
%
%-Exceptions
%
%   1)  If the spacecraft clock type is not supported, the error
%       SPICE(NOTSUPPORTED) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If any of the extracted clock components cannot be parsed as
%       integers, or the string has too many components, or the value
%       of one of the components is less than the offset value, then,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   3)  Invalid spacecraft ID's are not diagnosed.
%
%   4)  If any of the input arguments, `sc' or `clkstr', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `sc' or `clkstr', is not of the
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Added complete examples to -Examples section, using CASSINI PDS
%       archived data.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Extended argument detailed descriptions in the -I/O section.
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
%   -Mice Version 1.0.0, 07-JUN-2006 (EDW)
%
%-Index_Entries
%
%   convert spacecraft_clock string to ticks
%
%-&

function [ticks] = cspice_sctiks(sc, clkstr)

   switch nargin
      case 2

         sc     = zzmice_int(sc);
         clkstr = zzmice_str(clkstr);

      otherwise

         error ( 'Usage: [_ticks_] = cspice_sctiks(sc, _`clkstr`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [ticks] = mice('sctiks_c', sc, clkstr);
   catch spiceerr
      rethrow(spiceerr)
   end





