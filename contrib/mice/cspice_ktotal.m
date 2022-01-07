%-Abstract
%
%   CSPICE_KTOTAL returns the number of kernels of a specified type that are
%   currently loaded via the cspice_furnsh interface.
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
%      kind     a list of types of kernels to count when computing loaded
%               kernels.
%
%               [1,c1] = size(kind); char = class(kind)
%
%                  or
%
%               [1,1] = size(kind); cell = class(kind)
%
%               `kind' should consist of a list of words of kernels to
%               examine. Recognized types are
%
%                  SPK  --- All SPK files are counted in the total.
%                  CK   --- All CK files are counted in the total.
%                  PCK  --- All binary PCK files are counted in the
%                           total.
%                  DSK  --- All DSK files are counted in the total.
%                  EK   --- All EK files are counted in the total.
%                  TEXT --- All text kernels that are not meta-text
%                           kernels are included in the total.
%                  META --- All meta-text kernels are counted in the
%                           total.
%                  ALL  --- Every type of kernel is counted in the
%                           total.
%
%               `kind' is case insensitive. If a word appears in `kind'
%               that is not one of those listed above, it is ignored.
%
%               When `kind' consists of multiple words, the words must
%               be separated by blanks. Examples of valid lists are the
%               strings
%
%                  'SPK CK TEXT'
%                  'SPK CK text'
%                  'PCK DSK'
%                  'CK'
%                  'ALL'
%
%               See the -Examples section for illustrations of the
%               use of `kind'.
%
%   the call:
%
%      [count] = cspice_ktotal( kind )
%
%   returns:
%
%      count    the number of kernels loaded through cspice_furnsh that
%               belong to the list specified by `kind'.
%
%               [1,1] = size(count); int32 = class(count)
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
%   1) Load a meta-kernel with a PCK, an LSK and an SPK, and
%      separately, a text kernel and a binary PCK. Show the
%      total number of kernels and meta-kernels loaded. Determine the
%      number of text kernels loaded, and the number of binary
%      kernels.
%
%      Unload all kernels and clear the kernel pool using
%      cspice_kclear, and check that no kernels are loaded.
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: ktotal_ex1.tm
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
%            de421.bsp                     Planetary ephemeris
%            pck00008.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00008.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Use the PCK kernel below as the binary PCK required for the
%      example.
%
%         earth_latest_high_prec.bpc
%
%
%      Use the FK kernel below as the text kernel required for the
%      example.
%
%         RSSD0002.TF
%
%
%      Example code begins here.
%
%
%      function ktotal_ex1()
%
%         %
%         % Load several kernel files.
%         %
%         cspice_furnsh( 'ktotal_ex1.tm' )
%         cspice_furnsh( 'RSSD0002.TF' )
%         cspice_furnsh( 'earth_latest_high_prec.bpc' )
%
%         %
%         % Count the number of loaded kernel files.
%         %
%         n   = cspice_ktotal( 'ALL' );
%         txt = sprintf(['The total number of kernels after ',             ...
%                        'cspice_kclear call: %d'], n        );
%         disp( txt )
%
%         %
%         % Count the number of meta-kernels.
%         %
%         n   = cspice_ktotal( 'META' );
%         txt = sprintf(['The total number of meta-kernels  ',             ...
%                        '                  : %d'], n        );
%         disp( txt )
%
%         %
%         % Count the number of text kernels.
%         %
%         n   = cspice_ktotal( 'TEXT' );
%         txt = sprintf(['The total number of text kernels  ',             ...
%                        '                  : %d'], n        );
%         disp( txt )
%
%         %
%         % Count the number of binary kernels. These kernels
%         % are of type CK, DSK, EK, PCK or SPK.
%         %
%         n   = cspice_ktotal( 'CK DSK EK PCK SPK' );
%         txt = sprintf(['The total number of binary kernels',             ...
%                        '                  : %d'], n        );
%         disp( txt )
%
%         %
%         % Clear the KEEPER system, retrieve the number of loaded
%         % after the clear.
%         %
%         cspice_kclear
%
%         n   = cspice_ktotal( 'ALL' );
%         txt = sprintf(['The total number of kernels after ',             ...
%                        'cspice_kclear     : %d'], n        );
%         disp( txt )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      The total number of kernels after cspice_kclear call: 6
%      The total number of meta-kernels                    : 1
%      The total number of text kernels                    : 3
%      The total number of binary kernels                  : 2
%      The total number of kernels after cspice_kclear     : 0
%
%
%-Particulars
%
%   cspice_ktotal allows you to easily determine the number of kernels
%   loaded via the interface cspice_furnsh that are of a type of interest.
%
%-Exceptions
%
%   1)  If a word on the list specified by `kind' is not recognized,
%       it is ignored.
%
%   2)  If `kind' is blank, or none of the words in `kind' is on the
%       list specified above, `count' will be returned as zero.
%
%   3)  If the input argument `kind' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `kind' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
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
%   KERNEL.REQ
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 2.1.0, 01-NOV-2021 (EDW) (JDR) (NJB)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and meta-kernel. Merged the existing
%       code fragments into a complete example.
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
%       Updated -I/O description of input argument "kind" to illustrate
%       use of multi-word lists. Added kernel.req to the list of required
%       readings, and removed dsk.req. Improved -Abstract section.
%       Corrected class type description for output argument `count', to
%       int32.
%
%   -Mice Version 2.0.0, 20-JAN-2016 (EDW) (NJB)
%
%       Header update to reflect support for use of DSKs. Corrected
%       class type description for output argument `count', to double.
%
%   -Mice Version 1.0.2, 01-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 01-DEC-2006 (EDW)
%
%-Index_Entries
%
%   Number of loaded kernels of a given type
%
%-&

function [count] = cspice_ktotal( kind )

   switch nargin
      case 1

         kind = zzmice_str(kind);

      otherwise

         error ( 'Usage: count = cspice_ktotal(`kind`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [count] = mice( 'ktotal_c', kind );

      %
      % Convert the integers returned from the interface to double precision
      % in case a user includes the return arguments in a calculation
      % with other doubles.
      %
      count = zzmice_dp(count);

   catch spiceerr
      rethrow(spiceerr)
   end



