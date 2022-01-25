%-Abstract
%
%   CSPICE_SPKSFS searches through loaded SPK files to find the
%   highest-priority segment applicable to the body and time specified.
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
%      body   the SPK ID code of an ephemeris object, typically a solar
%             system body.
%
%             [1,1] = size(dc); int32 = class(body)
%
%      et     the time, in seconds past the epoch J2000 TDB.
%
%             [1,1] = size(et); double = class(et)
%
%   the call:
%
%      [handle, descr, ident, found] = cspice_spksfs( body, et)
%
%   returns:
%
%      handle   the handle of the SPK file containing a located segment.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      descr    the descriptor of a located SPK segment.
%
%               [5,1] = size(descr); double = class(descr)
%
%      ident    the string SPK segment identifier of a located SPK segment.
%
%               [1,c1] = size(ident); char = class(ident)
%
%      found    indicates whether a requested segment was found or not.
%               The other output arguments are valid only if `found'
%               is set to true.
%
%               [1,1] = size(found); logical = class(found)
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
%   1) Find a segment for the Pluto barycenter, with coverage for
%      a specified epoch, in a JPL planetary SPK file, and display
%      the segment's information.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: spksfs_ex1.tm
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
%            naif0010.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'naif0010.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function spksfs_ex1()
%
%         %
%         % Local constants
%         %
%         META   =  'spksfs_ex1.tm';
%         ND     =  2;
%         NI     =  6;
%
%         %
%         % Load meta-kernel.
%         %
%         cspice_furnsh( META )
%
%         %
%         % Convert starting time to seconds past J2000 TDB.
%         %
%         timstr = '2012 APR 27 00:00:00.000 TDB';
%
%         et0 = cspice_str2et(timstr);
%
%         %
%         % Get the NAIF ID code for the Pluto system barycenter.
%         % This is a built-in ID code, so something's seriously
%         % wrong if we can't find the code.
%         %
%         [idcode, found] = cspice_bodn2c( 'PLUTO BARYCENTER' );
%
%         if ~found
%            cspice_kclear
%            errot( 'SPICE(BUG)' )
%         end
%
%         [handle, descr, segid, found] = cspice_spksfs( idcode, et0);
%
%
%         if ~found
%            cspice_kclear
%            txt = sprintf( 'No SPK segment found for body %d at time %s', ...
%                            body, timstr );
%            error( txt )
%         end
%
%         %
%         % Unpack the descriptor of the current segment.
%         %
%         [dc, ic] = cspice_dafus( descr, ND, NI );
%
%         frname = cspice_frmnam( ic(3) );
%
%         fprintf( 'Body        = %d\n', ic(1) )
%         fprintf( 'Center      = %d\n', ic(2) )
%         fprintf( 'Frame       = %s\n', frname)
%         fprintf( 'Data type   = %d\n', ic(4) )
%         fprintf( 'Start ET    = %f\n', dc(1) )
%         fprintf( 'Stop ET     = %f\n', dc(2) )
%         fprintf( 'Segment ID  = %s\n\n', segid )
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
%      Body        = 9
%      Center      = 0
%      Frame       = J2000
%      Data type   = 2
%      Start ET    = -3169195200.000000
%      Stop ET     = 1696852800.000000
%      Segment ID  = DE-0421LE-0421
%
%
%-Particulars
%
%   This routine finds the highest-priority segment, in any loaded
%   SPK file, such that the segment provides data for the specified
%   body and epoch.
%
%-Exceptions
%
%   1)  If an attempt is made to call cspice_spksfs when there aren't any
%       files loaded, the error SPICE(NOLOADEDFILES) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If any of the input arguments, `body' or `et', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `body' or `et', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   All files loaded by cspice_spklef are potential search targets for
%   cspice_spksfs.
%
%-Restrictions
%
%   1)  If Fortran i/o errors occur while searching a loaded SPK
%       file, the internal state of this suite of routines may
%       be corrupted. It may be possible to correct the state
%       by unloading the pertinent SPK files and then re-loading
%       them.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited -Examples section to comply with NAIF standard. Changed
%       example code to focus on cspice_spksfs example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 30-OCT-2012 (EDW)
%
%-Index_Entries
%
%   select SPK file and segment
%
%-&

function [handle, descr, ident, found] = cspice_spksfs(body, et)

   switch nargin
      case 2

         body = zzmice_int(body);
         et   = zzmice_dp(et);

      otherwise

         error ( ['Usage: [handle, descr(5), `ident`, found] =' ...
                              ' cspice_spksfs( body, et)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      spksfs = mice( 'spksfs_s', body, et );

      handle  = reshape( [spksfs.handle], 1, [] );
      descr   = reshape( [spksfs.descr ], 5, [] );
      ident   = char( spksfs.ident );
      found   = reshape( [spksfs.found ], 1, [] );

   catch spiceerr
      rethrow(spiceerr)
   end
