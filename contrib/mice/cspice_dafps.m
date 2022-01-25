%-Abstract
%
%   CSPICE_DAFPS packs (assembles) an array summary from its double precision
%   and integer components.
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
%      dc       the double precision components of the summary.
%
%               [1,nd] = size(dc); double = class(dc)
%
%      ic       the integer components of the summary.
%
%               [1,ni] = size(ic); int32 = class(ic)
%
%   the call:
%
%      [sum] = cspice_dafps( dc, ic )
%
%   returns:
%
%      sum      an array summary containing the components in `dc' and `ic'.
%
%               [n,1] = size(sum); double = class(sum)
%
%               This identifies the contents and location of a single array
%               within a DAF.
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
%   The components of array summaries are packed into double precision arrays
%   for reasons outlined in the DAF required reading. Two routines,
%   cspice_dafps (pack summary) and cspice_dafus (unpack summary) are
%   provided for packing and unpacking summaries.
%
%   The total size of the summary is
%
%            (ni - 1)
%      nd + ---------- + 1
%               2
%
%   double precision words (where the sizes of `dc' and `ic', nd and ni,
%   are non-negative).
%
%-Exceptions
%
%   1)  If `nd', the number of double precision components, is zero or
%       negative, no DP components are stored.
%
%   2)  If `ni', the number of integer components, is zero or
%       negative, no integer components are stored.
%
%   3)  If the total size of the summary is greater than 125 double
%       precision words, some components may not be stored. See
%       -Particulars for details.
%
%   4)  If any of the input arguments, `dc' or `ic', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   5)  If any of the input arguments, `dc' or `ic', is not of the
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
%   -Mice Version 1.0.0, 26-NOV-2021 (JDR)
%
%-Index_Entries
%
%   pack DAF summary
%
%-&
function [sum] = cspice_dafps( dc, ic )

   switch nargin
      case 2

         dc = zzmice_dp(dc);
         ic = zzmice_int(ic);

      otherwise

         error ( 'Usage: [sum()] = cspice_dafps( dc(nd), ic(ni) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [sum] = mice('dafps_c', dc, ic);
   catch spiceerr
      rethrow(spiceerr)
   end
