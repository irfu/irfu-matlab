%-Abstract
%
%   CSPICE_CKNR03 returns the number of pointing instances in a CK type 03
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
%      descr    the packed descriptor of a data type 3 CK segment.
%
%               [5,1] = size(descr); double = class(descr)
%
%   the call:
%
%      [nrec] = cspice_cknr03( handle, descr )
%
%   returns:
%
%      nrec     the number of pointing instances in the type 3 segment.
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
%   1) The following code example extracts the SCLK time, boresight
%      vector, and angular velocity vector for each pointing instance
%      in the first segment in a CK file that contains segments of
%      data type 3.
%
%      Use the CK kernel below, available in the Venus Express PDS
%      archives, as input for the code example.
%
%         VEX_BOOM_V01.BC
%
%      Example code begins here.
%
%
%      function cknr03_ex1()
%
%         %
%         % First load the file (it may also be opened by
%         % using cspice_cklpf).
%         %
%         [handle] = cspice_dafopr( 'VEX_BOOM_V01.BC' );
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
%         [dcd, icd] = cspice_dafgs( 2, 6 );
%         [descr]    = cspice_dafps( dcd, icd );
%
%         %
%         % The data type for a segment is located in the third
%         % integer component of the descriptor.
%         %
%         if ( icd(3) == 3 )
%
%            %
%            % Does the segment contain `av' data?
%            %
%            avseg =  ( icd(4) == 1 );
%
%            %
%            % How many records does this segment contain?
%            %
%            [nrec] = cspice_cknr03( handle, descr );
%
%            for i=1:nrec
%
%               %
%               % Get the ith pointing instance in the segment.
%               %
%               [record] = cspice_ckgr03( handle, descr, i );
%
%               %
%               % Unpack `record' into the time, quaternion, and av.
%               %
%               sclkdp = record(1);
%               quat   = record(2:5);
%
%               if ( avseg )
%
%                  av = record(6:8);
%
%               end
%
%               %
%               % The boresight vector is the third row of the
%               % C-matrix.
%               %
%               [cmat] = cspice_q2m( quat );
%               bore   = cmat(3,:);
%
%               %
%               % Write out the results.
%               %
%               fprintf( 'Record: %2d\n', i )
%               fprintf( '   SCLK time       : %24.6f\n', sclkdp )
%               fprintf( '   Boresight       : %13.9f %13.9f %13.9f\n',    ...
%                                                                  bore )
%
%               if ( avseg )
%
%                  fprintf( '   Angular velocity: %13.9f %13.9f %13.9f\n', ...
%                                                                      av  )
%
%               end
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
%         SCLK time       :           2162686.710986
%         Boresight       :  -0.999122830   0.000000000   0.041875654
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record:  2
%         SCLK time       :       54160369751.715164
%         Boresight       :  -0.999122830   0.000000000   0.041875654
%         Angular velocity:   0.000000000   1.176083393   0.000000000
%
%      Record:  3
%         SCLK time       :       54160454948.487686
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record:  4
%         SCLK time       :      299264885854.937805
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record:  5
%         SCLK time       :     2366007685832.532227
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record:  6
%         SCLK time       :     4432750485810.126953
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record:  7
%         SCLK time       :     6505155594828.757812
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record:  8
%         SCLK time       :     8571898394806.352539
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record:  9
%         SCLK time       :    10638641194783.947266
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record: 10
%         SCLK time       :    12705383994761.541016
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record: 11
%         SCLK time       :    14777789103780.169922
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record: 12
%         SCLK time       :    16844531903757.763672
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%      Record: 13
%         SCLK time       :    18911274703735.359375
%         Boresight       :   0.000000000   0.000000000   1.000000000
%         Angular velocity:   0.000000000   0.000000000   0.000000000
%
%
%-Particulars
%
%   For a complete description of the internal structure of a type 3
%   segment, see the CK required reading.
%
%   This routine returns the number of discrete pointing instances
%   contained in the specified segment. It is normally used in
%   conjunction with cspice_ckgr03 which returns the Ith pointing instance
%   in the segment.
%
%-Exceptions
%
%   1)  If the segment indicated by `descr' is not a type 3 segment, the
%       error SPICE(CKWRONGDATATYPE) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If the specified handle does not belong to any DAF file that
%       is currently known to be open, an error is signaled by a
%       routine in the call tree of this routine.
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
%   The CK file specified by `handle' should be open for read or
%   write access.
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
%   number of CK type_3 records
%
%-&
function [nrec] = cspice_cknr03( handle, descr )

   switch nargin
      case 2

         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr, true);

      otherwise

         error ( 'Usage: [nrec] = cspice_cknr03( handle, descr(5) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [nrec] = mice('cknr03_c', handle, descr);
   catch spiceerr
      rethrow(spiceerr)
   end
