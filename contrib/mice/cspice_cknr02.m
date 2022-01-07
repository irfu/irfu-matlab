%-Abstract
%
%   CSPICE_CKNR02 returns the number of pointing records in a CK type 02
%   segment. The segment is identified by a CK file handle and segments
%   descriptor.
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
%      handle   the handle of the binary CK file containing the segment.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               The file should have been opened for read or write access,
%               by cspice_cklpf, cspice_dafopr or cspice_dafopw.
%
%      descr    the packed descriptor of a data type 2 CK segment.
%
%               [5,1] = size(descr); double = class(descr)
%
%   the call:
%
%      [nrec] = cspice_cknr02( handle, descr )
%
%   returns:
%
%      nrec     the number of pointing records in the type 2 segment
%               associated with `handle' and `descr'.
%
%               [1,1] = size(nrec); int32 = class(nrec)
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
%   1) The following code example extracts the start and end SCLK
%      time, seconds per tick rate, platform's +Z axis direction,
%      and angular velocity vector for each pointing instance in
%      the first segment in a CK file that contains segments of data
%      type 2.
%
%      Use the CK kernel below, available in the Viking Orbiter PDS
%      archives, as input for the code example.
%
%         vo2_swu_ck2.bc
%
%      Example code begins here.
%
%
%      function cknr02_ex1()
%
%         %
%         % First load the file. (The file may also be opened by
%         % using cspice_cklpf).
%         %
%         [handle] = cspice_dafopr( 'vo2_swu_ck2.bc' );
%
%         %
%         % Begin forward search.  Find the first array.
%         %
%         cspice_dafbfs( handle );
%         [found] = cspice_daffna;
%
%         %
%         % Get segment descriptor.
%         %
%         % Unpack the segment descriptor into its double precision
%         % and integer components.
%         %
%         [dcd, icd] = cspice_dafgs( 2,   6   );
%         [descr]    = cspice_dafps( dcd, icd );
%
%         %
%         % The data type for a segment is located in the third
%         % integer component of the descriptor.
%         %
%         if ( icd(3) == 2 )
%
%            %
%            % How many records does this segment contain?
%            %
%            [nrec] = cspice_cknr02( handle, descr );
%
%            for i=1:nrec
%
%               %
%               % Get the ith record in the segment.
%               %
%               [record] = cspice_ckgr02( handle, descr, i );
%
%               %
%               % Unpack `record' into the start and end time, rate in
%               % TDB seconds/tick, quaternion, and av.
%               %
%               sclks = record(1);
%               sclke = record(2);
%               sclkr = record(3);
%
%               quat  = record(4:7);
%               av    = record(8:10);
%
%               %
%               % The +Z axis direction is the third row of the
%               % C-matrix.
%               %
%               [cmat] = cspice_q2m( quat );
%
%               z      = cmat(:,3);
%
%               %
%               % Write out the results.
%               %
%               fprintf( 'Record:  %1d\n', i )
%               fprintf( '   Start encoded SCLK: %20.6f\n', sclks )
%               fprintf( '   End encoded SCLK  : %20.6f\n', sclke )
%               fprintf( '   TDB Seconds/tick  : %20.6f\n', sclkr )
%               fprintf( '   +Z axis           : %12.8f %12.8f %12.8f\n', z  )
%               fprintf( '   Angular velocity  : %12.8f %12.8f %12.8f\n', av )
%               fprintf( '\n' )
%
%            end
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Record:  1
%         Start encoded SCLK:   32380393707.000015
%         End encoded SCLK  :   32380395707.000015
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.37100754  -0.58030627   0.72498141
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  2
%         Start encoded SCLK:   32380402605.999947
%         End encoded SCLK  :   32380404605.999947
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.36011140  -0.57851153   0.73187716
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  3
%         Start encoded SCLK:   32380412542.000053
%         End encoded SCLK  :   32380414542.000053
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.34654725  -0.57947579   0.73764003
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  4
%         Start encoded SCLK:   32827264875.000000
%         End encoded SCLK  :   32827266875.000000
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.52790086  -0.78142544   0.33270853
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  5
%         Start encoded SCLK:   32827403805.999992
%         End encoded SCLK  :   32827405805.999992
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.53591430  -0.75336795   0.38109395
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  6
%         Start encoded SCLK:   32827412705.000042
%         End encoded SCLK  :   32827414705.000042
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.52897927  -0.75383975   0.38975193
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  7
%         Start encoded SCLK:   32827417284.000038
%         End encoded SCLK  :   32827419284.000038
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.52739277  -0.75283638   0.39382008
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  8
%         Start encoded SCLK:   33793314593.000053
%         End encoded SCLK  :   33793316593.000053
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.54702960  -0.69406245   0.46801275
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  9
%         Start encoded SCLK:   33793332478.000046
%         End encoded SCLK  :   33793334478.000046
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.55111040  -0.67525183   0.49021658
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  10
%         Start encoded SCLK:   33793341463.000061
%         End encoded SCLK  :   33793343463.000061
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.55538431  -0.66374748   0.50098659
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  11
%         Start encoded SCLK:   33793350363.000034
%         End encoded SCLK  :   33793352363.000034
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.55830213  -0.65303671   0.51170478
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  12
%         Start encoded SCLK:   33984028250.000000
%         End encoded SCLK  :   33984030250.000000
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.56718006  -0.65737847   0.49614546
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  13
%         Start encoded SCLK:   33984046134.999992
%         End encoded SCLK  :   33984048134.999992
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.57392845  -0.63511679   0.51694564
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  14
%         Start encoded SCLK:   33984055121.000053
%         End encoded SCLK  :   33984057121.000053
%         TDB Seconds/tick  :             0.001000
%         +Z axis           :  -0.57600983  -0.62488782   0.52699894
%         Angular velocity  :   0.00000000   0.00000000   0.00000000
%
%      Record:  15
%         Start encoded SCLK:   33984220835.999966
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 875 lines have been
%      provided.
%
%
%-Particulars
%
%   For a complete description of the internal structure of a type 2
%   segment, see the CK required reading.
%
%   This routine returns the number of pointing records contained
%   in the specified segment. It is normally used in conjunction
%   with cspice_ckgr02, which returns the Ith record in the segment.
%
%-Exceptions
%
%   1)  If the segment indicated by `descr' is not a type 2 segment, the
%       error SPICE(CKWRONGDATATYPE) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If the specified handle does not belong to any file that is
%       currently known to be open, an error is signaled by a routine
%       in the call tree of this routine.
%
%   3)  If `descr' is not a valid descriptor of a segment in the CK
%       file specified by `handle', the results of this routine are
%       unpredictable.
%
%   4)  If any of the input arguments, `handle' or `descr', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `handle' or `descr', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   The CK file specified by `handle' should be open for read or write
%   access.
%
%-Restrictions
%
%   1)  The binary CK file containing the segment whose descriptor was
%       passed to this routine must be opened for read or write access
%       by cspice_cklpf, cspice_dafopr or cspice_dafopw.
%
%-Required_Reading
%
%   CK.REQ
%   DAF.REQ
%   MICE.REQ
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
%   -Mice Version 1.0.0, 01-NOV-2021 (JDR)
%
%-Index_Entries
%
%   number of CK type_2 records
%
%-&
function [nrec] = cspice_cknr02( handle, descr )

   switch nargin
      case 2

         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr);

      otherwise

         error ( 'Usage: [nrec] = cspice_cknr02( handle, descr(5) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [nrec] = mice('cknr02_c', handle, descr);
   catch spiceerr
      rethrow(spiceerr)
   end
