%-Abstract
%
%   CSPICE_UNITIM transforms time from one uniform scale to another. The
%   uniform time scales are TAI, GPS, TT, TDT, TDB, ET, JED, JDTDB, JDTDT.
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
%      epoch    the epoch(s) relative to the `insys' time scale.
%
%               [1,n] = size(epoch); double = class(epoch)
%
%      insys    a time scale.
%
%               [1,c1] = size(insys); char = class(insys)
%
%                  or
%
%               [1,1] = size(insys); cell = class(insys)
%
%               Acceptable values are:
%
%                  'TAI'     International Atomic Time.
%                  'TDB'     Barycentric Dynamical Time.
%                  'TDT'     Terrestrial Dynamical Time.
%                  'TT'      Terrestrial Time, identical to TDT.
%                  'ET'      Ephemeris time (in the SPICE system, this is
%                            equivalent to TDB).
%                  'JDTDB'   Julian Date relative to TDB.
%                  'JDTDT'   Julian Date relative to TDT.
%                  'JED'     Julian Ephemeris date (in the SPICE system
%                            this is equivalent to JDTDB).
%                  'GPS'     Global Positioning System Time.
%
%               The routine is not sensitive to the case of the
%               characters in `insys'; 'tai' 'Tai' and 'TAI' are all
%               equivalent from the point of view of this routine.
%
%      outsys   the time scale to which `epoch' should be converted.
%
%               [1,c2] = size(outsys); char = class(outsys)
%
%                  or
%
%               [1,1] = size(outsys); cell = class(outsys)
%
%               Acceptable values are the same as for `insys'. The routine
%               is not sensitive to the case of `outsys'.
%
%   the call:
%
%      [unitim] = cspice_unitim( epoch, insys, outsys )
%
%   returns:
%
%      unitim   the time(s) in the system specified by `outsys' that is
%               equivalent to the `epoch' in the `insys' time scale.
%
%               [1,n] = size(unitim); double = class(unitim)
%
%               `unitim' returns with the same vectorization measure, N,
%               as `epoch'.
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
%   1) Convert an input UTC string to Julian Ephemeris Date seconds.
%
%      Use the LSK kernel below to load the leap seconds and time
%      constants required for the conversions.
%
%         naif0012.tls
%
%
%      Example code begins here.
%
%
%      function unitim_ex1()
%
%         %
%         % Load a leapseconds kernel.
%         %
%         cspice_furnsh( 'naif0012.tls' )
%
%         utcstr = 'Dec 19 2003  16:48:00';
%         et     = cspice_str2et( utcstr );
%
%         converted_et = cspice_unitim(et, 'ET','JED');
%
%         fprintf( 'UTC time             : %s\n', utcstr )
%         fprintf( 'Ephemeris time       : %21.9f\n', et )
%         fprintf( 'Julian Ephemeris Date: %21.9f\n', converted_et )
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
%      UTC time             : Dec 19 2003  16:48:00
%      Ephemeris time       :   125124544.183560610
%      Julian Ephemeris Date:     2452993.200742865
%
%
%-Particulars
%
%   We use the term uniform time scale to refer to those
%   representations of time that are numeric (each epoch is
%   represented by a number) and additive. A numeric time system is
%   additive if given the representations, `e1' and `e2', of any pair of
%   successive epochs, the time elapsed between the epochs is given by
%   e2 - e1.
%
%   Given an epoch in one of the uniform time scales specified by
%   `insys', the function returns the equivalent representation in the
%   scale specified by `outsys'. A list of the recognized uniform time
%   scales is given in the detailed input for `insys'.
%
%-Exceptions
%
%   1)  The kernel pool must contain the variables:
%
%          'DELTET/DELTA_T_A'
%          'DELTET/K'
%          'DELTET/EB'
%          'DELTET/M'
%
%       If these are not present, the error SPICE(MISSINGTIMEINFO) is
%       signaled by a routine in the call tree of this routine. (These
%       variables are typically inserted into the kernel pool by
%       loading a leapseconds kernel with the SPICE routine cspice_furnsh.)
%
%   2)  If the names of either the input or output time types are
%       unrecognized, the error SPICE(BADTIMETYPE) is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If any of the input arguments, `epoch', `insys' or `outsys',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   4)  If any of the input arguments, `epoch', `insys' or `outsys',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The appropriate variable must be loaded into the SPICE kernel
%       pool (normally by loading a leapseconds kernel with cspice_furnsh)
%       prior to calling this routine.
%
%-Required_Reading
%
%   MICE.REQ
%   TIME.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "output" to "unitim" to comply with NAIF
%       standard.
%
%       Added time system name 'TT' (Terrestrial Time) as alternate
%       assignment of 'TDT' (Terrestrial Dynamical Time).
%
%       Included GPS time system mapping.
%
%       Edited the header to comply with NAIF standard. Added a reference to
%       the required LSK. Changed example's output format.
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
%   -Mice Version 1.0.0, 12-MAR-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   Transform between two uniform numeric time systems
%   Transform between two additive numeric time systems
%   Convert one uniform numeric time system to another
%   Convert one additive numeric time system to another
%
%-&

function [unitim] = cspice_unitim( epoch, insys, outsys )

   switch nargin
      case 3

         epoch  = zzmice_dp(epoch);
         insys  = zzmice_str(insys);
         outsys = zzmice_str(outsys);

      otherwise

         error( [ 'Usage: [_unitim_] = '                                    ...
                  'cspice_unitim( _epoch_, `insys`, `outsys` )' ] )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [unitim] = mice( 'unitim_c', epoch, insys, outsys);
   catch spiceerr
      rethrow(spiceerr)
   end
