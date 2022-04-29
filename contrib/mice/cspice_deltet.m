%-Abstract
%
%   CSPICE_DELTET returns value of Delta ET (ET-UTC)
%   for an input epoch.
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
%      epoch    the epoch(s) at which "delta ET" is to be computed.
%
%               [1,1] = size(epoch); double = class(epoch)
%
%               `epoch' may be either UTC or ephemeris seconds past J2000,
%               as specified by `eptype'.
%
%      eptype   the type of input epoch.
%
%               [1,c1] = size(eptype); char = class(eptype)
%
%                  or
%
%               [1,1] = size(eptype); cell = class(eptype)
%
%               It may be either of the following:
%
%                  'UTC'    UTC seconds past J2000 UTC.
%
%                  'ET'     Ephemeris seconds past J2000 TDB,
%                           also known as barycentric dynamical
%                           time (TDB).
%
%   the call:
%
%      [delta] = cspice_deltet( epoch, eptype )
%
%   returns:
%
%      delta    the value of
%
%                  "delta ET" = ET - UTC
%
%               at the input epoch.
%
%               [1,1] = size(delta); double = class(delta)
%
%               This is added to UTC to give ET, or subtracted from ET to
%               give UTC. The routine is reversible: that is, given the
%               following calls,
%
%                  [del1] = cspice_deltet( utc,      'UTC' );
%                  [del2] = cspice_deltet( utc+del1, 'ET'  );
%
%               the expression
%
%                  ( del1 == del2 )
%
%               is always true.
%
%               `delta' returns with the same vectorization measure
%               (N) as `epoch'.
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
%   1) Calculate the ET to UTC delta times in seconds, at January 1, 1997
%      and January 1, 2004, and for every Julian year in-between.
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
%      function deltet_ex1()
%
%         %
%         % Load a leapsecond file.
%         %
%         cspice_furnsh( 'naif0012.tls' )
%
%         %
%         % Define times of interest and the array size
%         % parameter.
%         %
%         SIZE     = 2004 - 1997 +1;
%         UTC_1997 = '1997 JAN 01 00:00:00.000';
%         UTC_2004 = '2004 JAN 01 00:00:00.000';
%
%         %
%         % Convert the UTC time strings to ET.
%         %
%         et_1997 = cspice_str2et( UTC_1997 );
%         et_2004 = cspice_str2et( UTC_2004 );
%
%         %
%         % Calculate the ET-UTC delta at Jan 1 1997
%         % and Jan 1 2004.
%         %
%         delt_1997 = cspice_deltet( et_1997, 'ET' );
%         delt_2004 = cspice_deltet( et_2004, 'ET' );
%
%         disp( '     UTC time             Delta ET-UTC' )
%         disp( '------------------------  ------------' )
%         disp( 'Scalar:' )
%         fprintf( '%s  %12.8f\n', UTC_1997, delt_1997 )
%         fprintf( '%s  %12.8f\n', UTC_2004, delt_2004 )
%
%         %
%         % Given an array of 'SIZE' ephemeris times
%         % starting from value 'et_1997' with steps being
%         % of the number of seconds per Julian year, return
%         % the ET-UTC delta value for each time.
%         %
%         et   = [0:SIZE-1]*cspice_jyear + et_1997;
%         delt = cspice_deltet( et, 'ET' );
%
%         %
%         % Convert 'et' to 'utc'.
%         %
%         utc = cspice_et2utc( et, 'C', 3 );
%
%         disp( 'Vector:' )
%         for n=1:SIZE
%            fprintf( '%s  %12.8f\n', utc(n,:), delt(n) );
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
%           UTC time             Delta ET-UTC
%      ------------------------  ------------
%      Scalar:
%      1997 JAN 01 00:00:00.000   62.18393536
%      2004 JAN 01 00:00:00.000   64.18391170
%      Vector:
%      1997 JAN 01 00:00:00.000   62.18393536
%      1998 JAN 01 05:59:59.000   63.18393508
%      1999 JAN 01 11:59:58.000   64.18393480
%      2000 JAN 01 17:59:58.000   64.18393452
%      2000 DEC 31 23:59:58.000   64.18393424
%      2002 JAN 01 05:59:58.000   64.18393396
%      2003 JAN 01 11:59:58.000   64.18393368
%      2004 JAN 01 17:59:58.000   64.18393341
%
%
%-Particulars
%
%   The constants necessary for computing the offset are taken
%   from the kernel pool, where they are assumed to have been
%   loaded from a kernel file.
%
%   The tables are consulted to determine the number of leap seconds
%   preceding the input epoch. Also, an approximation to the periodic
%   yearly variation (which has an amplitude of just under two
%   milliseconds) in the difference between ET and TAI (Atomic Time)
%   is computed. The final value of Delta ET is given by
%
%      Delta ET = ( et - tai ) + leap seconds
%
%-Exceptions
%
%   1)  If the input epoch is not recognized, the error
%       SPICE(INVALIDEPOCH) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the variables necessary for the computation of `delta' have
%       not been loaded into the kernel pool, the error
%       SPICE(KERNELVARNOTFOUND) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the number of leapseconds in the pool is greater than the
%       local leapseconds buffer size, the error SPICE(BUFFEROVERFLOW)
%       is signaled by a routine in the call tree of this routine.
%
%   4)  If any of the input arguments, `epoch' or `eptype', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `epoch' or `eptype', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The routines cspice_str2et and cspice_et2utc are preferred for
%       conversions between UTC and ET. This routine is provided mainly as a
%       utility for cspice_str2et and cspice_et2utc.
%
%   2)  A leapseconds kernel containing leapseconds and relativistic
%       terms MUST be loaded prior to calling this routine.
%       Examples demonstrating how to load a kernel pool are included
%       in the Required Reading file time.req and in the -Examples
%       section of this header. For more general information about
%       kernel pools, please consult the Required Reading file
%       kernel.req.
%
%-Required_Reading
%
%   MICE.REQ
%   TIME.REQ
%   KERNEL.REQ
%
%-Literature_References
%
%   [1]  "The Astronomical Almanac for the Year 1990," United States
%        Naval Observatory, U.S. Government Printing Office,
%        Washington, D.C., 1989.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and a reference to the required LSK.
%       Modified example's output.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 30-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   difference between ephemeris time and utc
%
%-&

function [delta] = cspice_deltet( epoch, eptype )

   switch nargin
      case 2

         eptype = zzmice_str(eptype);
         epoch  = zzmice_dp(epoch);

      otherwise

         error ( 'Usage: [_delta_] = cspice_deltet( _epoch_, `eptype`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [delta] = mice('deltet_c', epoch, eptype) ;
   catch spiceerr
      rethrow(spiceerr)
   end


