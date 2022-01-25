%-Abstract
%
%   CSPICE_DAFUS unpacks an array summary into its double precision and
%   integer components.
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
%      sum      an containing the DAF array summary. This identifies the
%               contents and location of a single array within a DAF.
%
%               [n,1] = size(sum); double = class(sum)
%
%      nd       the size of the return double precision array.
%
%               [1,1] = size(nd); int32 = class(nd)
%
%      ni       the size of the return integer array.
%
%               [1,1] = size(ni); int32 = class(ni)
%
%               For an SPK file, `nd' always equals 2, `ni' always equals 6.
%               The precise contents of the vectors depend on the type of DAF
%               but the final two elements of the `ic' (integer) vector always
%               contains the initial and final addresses respectively of the
%               array.
%
%   the call:
%
%      [dc, ic] = cspice_dafus( sum, nd, ni )
%
%   returns:
%
%      dc       the array of double precision components of the summary.
%
%               [1,nd] = size(dc); double = class(dc)
%
%      ic       the array of integer components of the summary.
%
%               [1,ni] = size(ic); int32 = class(ic)
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
%         File name: dafus_ex1.tm
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
%      function dafus_ex1()
%
%         %
%         % Local constants
%         %
%         META   =  'dafus_ex1.tm';
%         ND     =  2;
%         NI     =  6;
%
%         %
%         % Load a meta-kernel that specifies a planetary SPK file
%         % and leapseconds kernel. The contents of this meta-kernel
%         % are displayed above.
%         %
%         cspice_furnsh( META )
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
%            error( 'SPICE(BUG)' )
%         end
%
%         %
%         % Pick a request time; convert to seconds past J2000 TDB.
%         %
%         reqtim = '2011 FEB 18 UTC';
%
%         et = cspice_str2et( reqtim );
%
%         %
%         % Find a loaded segment for the specified body and time.
%         %
%
%         [handle, descr, segid, found] = cspice_spksfs( idcode, et );
%
%         if ~found
%            cspice_kclear
%            txt = sprintf( ['No descriptor found for the ',      ...
%                            'body %d at time %s'],  idcode, et );
%            error( txt )
%         else
%
%            %
%            % Display the DAF file handle.
%            %
%            fprintf( 'DAF handle:    %d\n', handle )
%
%            %
%            % Display the segment ID.
%            %
%            %
%            % Unpack the descriptor. Display the contents.
%            %
%            [dc, ic] = cspice_dafus( descr, ND, NI );
%
%            fprintf( 'Segment ID:       %s\n', segid )
%            fprintf( 'Body ID code:     %d\n', ic(1) )
%            fprintf( 'Center ID code:   %d\n', ic(2) )
%            fprintf( 'Frame ID code:    %d\n', ic(3) )
%            fprintf( 'SPK data type:    %d\n', ic(4) )
%            fprintf( 'Start time (TDB): %f\n', dc(1) )
%            fprintf( 'Stop time  (TDB): %f\n', dc(2) )
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
%      DAF handle:    1
%      Segment ID:       DE-0421LE-0421
%      Body ID code:     9
%      Center ID code:   0
%      Frame ID code:    1
%      SPK data type:    2
%      Start time (TDB): -3169195200.000000
%      Stop time  (TDB): 1696852800.000000
%
%
%-Particulars
%
%   The components of array summaries are packed into double
%   precision arrays.
%
%   The total size of the summary is
%
%           (ni - 1)
%      nd + -------- + 1
%               2
%
%   double precision words (where `nd', `ni' are nonnegative).
%
%-Exceptions
%
%   1)  If `nd' is zero or negative, no double precision components
%       are returned.
%
%   2)  If `ni' is zero or negative, no integer components are returned.
%
%   3)  If the total size of the summary is greater than 125 double
%       precision words, some components may not be returned.
%
%   4)  If any of the input arguments, `sum', `nd' or `ni', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `sum', `nd' or `ni', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%   DAF.REQ
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
%       Corrected size information of input argument "sum".
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and meta-kernel.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 29-OCT-2012 (EDW)
%
%-Index_Entries
%
%   get DAF summary
%
%-&

function [dc, ic] = cspice_dafus( sum, nd, ni )

   switch nargin
      case 3

         sum = zzmice_dp(sum);
         nd  = zzmice_int(nd);
         ni  = zzmice_int(ni);

      otherwise

         error ( 'Usage: [dc(nd), ic(ni)] = cspice_dafus( sum(N), nd, ni )' )

   end

   %
   % Call the MEX library.
   %
   try
      [dc, ic] = mice( 'dafus_c', sum, nd, ni );
   catch spiceerr
      rethrow(spiceerr)
   end

