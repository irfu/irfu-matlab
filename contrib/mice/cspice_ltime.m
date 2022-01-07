%-Abstract
%
%   CSPICE_LTIME computes the transmission (or reception) time of a signal at
%   a specified target, given the reception (or transmission) time at a
%   specified observer. This routine also returns the elapsed time between
%   transmission and reception.
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
%      etobs    the epoch(s) expressed as ephemeris seconds past J2000 TDB.
%
%               [1,n] = size(etobs); double = class(etobs)
%
%               This is the time at which an electromagnetic signal is "at"
%               the observer.
%
%      obs      the NAIF ID of some observer.
%
%               [1,1] = size(obs); int32 = class(obs)
%
%      dir      the direction the signal travels.
%
%               [1,2] = size(dir); char = class(dir)
%
%                  or
%
%               [1,1] = size(dir); cell = class(dir)
%
%               The acceptable values are '->' and '<-'. When you read the
%               calling sequence from left to right, the "arrow" given by
%               `dir' indicates which way the electromagnetic signal is
%               traveling.
%
%               If the argument list reads as below,
%
%                  ..., `obs', '->', `targ', ...
%
%               the signal is traveling from the observer to the
%               target.
%
%               If the argument reads as
%
%                  ..., `obs', '<-', `targ'
%
%               the signal is traveling from the target to
%               the observer.
%
%      targ     the NAIF ID of the target.
%
%               [1,1] = size(targ); int32 = class(targ)
%
%   the call:
%
%      [ettarg, elapsd] = cspice_ltime( etobs, obs, dir, targ )
%
%   returns:
%
%      ettarg   the epoch(s) expressed as ephemeris seconds past J2000 TDB
%               at which the electromagnetic signal is "at" the target body.
%
%               [1,n] = size(ettarg); double = class(ettarg)
%
%               Note `ettarg' is computed using only Newtonian
%               assumptions about the propagation of light.
%
%      elapsd   the number of ephemeris seconds (TDB) between transmission
%               and receipt of the signal.
%
%               [1,n] = size(elapsd); double = class(elapsd)
%
%               `elapsd' is computed as:
%
%                  elapsd = abs( etobs - ettarg )
%
%               `ettarg' and `elapsd' return with the same vectorization
%               measure, N, as `etobs'.
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
%   1) Suppose a signal is transmitted from Earth towards the Jupiter
%      system barycenter on July 4, 2004.
%
%         signal traveling to Jupiter system barycenter
%         *  -._.-._.-._.-._.-._.-._.-._.-._.->  *
%
%         Earth (399)            Jupiter system barycenter (5)
%
%      Compute the time at which the signal arrives at Jupiter
%      and the time it took the signal to arrive there (propagation
%      time).
%
%      Suppose also that another signal is received at the Earth from
%      Jupiter system barycenter at the same time.
%
%         signal sent from Jupiter system barycenter
%         *  <-._.-._.-._.-._.-._.-._.-._.-._.-  *
%
%         Earth (399)            Jupiter system barycenter (5)
%
%      Compute the time at which the signal was transmitted from Jupiter,
%      and its propagation time.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: ltime_ex1.tm
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
%            naif0012.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'naif0012.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function ltime_ex1()
%
%         %
%         %  Load an SPK and leapseconds kernel.
%         %
%         cspice_furnsh( 'ltime_ex1.tm' )
%
%         %
%         % Suppose a signal originates from Earth towards the
%         % Jupiter system barycenter. Define the NAIF IDs
%         % for the observer, Earth (399), the target, Jupiter
%         % barycenter (5), and time of interest.
%         %
%         OBS      = 399;
%         TARGET   = 5;
%         TIME_STR = 'July 4, 2004';
%
%         %
%         %  Convert the transmission time to ET.
%         %
%         et = cspice_str2et( TIME_STR);
%
%         %
%         %  Determine the arrival time and the time for propagation.
%         %
%         [arrive, lt] = cspice_ltime( et, OBS, '->', TARGET);
%
%         %
%         %  Convert the arrival time (ET) to UTC.
%         %
%         arrive_utc = cspice_et2utc( arrive, 'C', 3 );
%
%         %
%         %  Output the results.
%         %
%         txt = sprintf( 'Transmission at (UTC)       : %s', TIME_STR );
%         disp(txt)
%
%         txt = sprintf( 'The signal arrived at (UTC) : %s', arrive_utc );
%         disp(txt)
%
%         txt = sprintf( 'Time for propagation (secs) : %16.4f', lt );
%         disp(txt)
%         disp( ' ' )
%
%         %
%         % Now assume the signal originated at Jupiter barycenter,
%         % received by Earth at TIME_STR. Determine the transmission
%         % time and the time for propagation.
%         %
%         [receive, lt] = cspice_ltime( et, OBS, '<-', TARGET);
%
%         %
%         % Convert the reception time (ET) to UTC.
%         %
%         receive_utc = cspice_et2utc( receive, 'C', 3 );
%
%         %
%         %  Output the results.
%         %
%         txt = sprintf( 'Reception at (UTC)          : %s', TIME_STR );
%         disp(txt)
%
%         txt = sprintf( 'The signal sent at (UTC)    : %s', receive_utc );
%         disp(txt)
%
%         txt = sprintf( 'Time for propagation (secs) : %16.4f', lt );
%         disp(txt)
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
%      Transmission at (UTC)       : July 4, 2004
%      The signal arrived at (UTC) : 2004 JUL 04 00:48:38.717
%      Time for propagation (secs) :        2918.7171
%
%      Reception at (UTC)          : July 4, 2004
%      The signal sent at (UTC)    : 2004 JUL 03 23:11:21.248
%      Time for propagation (secs) :        2918.7525
%
%
%-Particulars
%
%   Suppose a radio signal travels between two solar system
%   objects. Given an ephemeris for the two objects, which way
%   the signal is traveling, and the time when the signal is
%   "at" at one of the objects (the observer `obs'), this routine
%   determines when the signal is "at" the other object (the
%   target `targ'). It also returns the elapsed time between
%   transmission and receipt of the signal.
%
%-Exceptions
%
%   1)  If `dir' is not one of '->' or '<-', the error
%       SPICE(BADDIRECTION) is signaled by a routine in the call tree
%       of this routine. In this case `ettarg' and `elapsd' will not be
%       altered from their input values.
%
%   2)  If insufficient ephemeris information is available to
%       compute the outputs `ettarg' and `elapsd', or if observer
%       or target is not recognized, an error is signaled
%       by a routine in the call tree of this routine.
%
%       In this case, the value of `ettarg' will be set to `etobs'
%       and `elapsd' will be set to zero.
%
%   3)  If any of the input arguments, `etobs', `obs', `dir' or
%       `targ', is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   4)  If any of the input arguments, `etobs', `obs', `dir' or
%       `targ', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
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
%   TIME.REQ
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
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and meta-kernel.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 10-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-JAN-2006 (EDW)
%
%-Index_Entries
%
%   Compute uplink and downlink light time
%
%-&

function [ettarg, elapsd] = cspice_ltime(etobs, obs, dir, targ)

   switch nargin
      case 4

         etobs = zzmice_dp(etobs);
         obs   = zzmice_int(obs);
         targ  = zzmice_int(targ);
         dir   = zzmice_str(dir);

      otherwise

         error ( ['Usage: [_ettarg_, _elapsd_] = ' ...
                  'cspice_ltime( _etobs_, obs, `dir`, targ)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [ettarg, elapsd] = mice('ltime_c',etobs, obs, dir, targ);
   catch spiceerr
      rethrow(spiceerr)
   end
