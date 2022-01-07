%-Abstract
%
%   CSPICE_SPKCLS closes a SPK file opened for read or write.
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
%      handle   the file handle for an open SPK file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_spkcls( handle )
%
%   closes the file attached to `handle'.
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
%   1) This example demonstrates how to create an SPK type 8 kernel
%      containing only one segment, given a time-ordered set of
%      discrete states and epochs.
%
%      Note that after run completion, a new SPK type 8 exists in the
%      output directory.
%
%      Example code begins here.
%
%
%      function spkcls_ex1()
%
%         %
%         % Define the segment identifier parameters.
%         %
%         BODY       = 3;
%         CENTER     = 10;
%         REF        = 'J2000';
%         POLY_DEG   = 3;
%         SPK8       = 'spkcls_ex1.bsp';
%         N_DISCRETE = 9;
%
%         %
%         % A set of epochs.
%         %
%         DISCRETEEPOCHS = (1:9)*100;
%
%         %
%         % An array of discrete states to write to the SPK segment.
%         %
%         base = [ (1:6)*100 ]';
%
%         %
%         % Create the 6xN array of states.
%         %
%         DISCRETESTATES = [(base+1), (base+2), (base+3), ...
%                           (base+4), (base+5), (base+6), ...
%                           (base+7), (base+8), (base+9) ];
%
%         %
%         % Create a segment identifier.
%         %
%         segid = 'SPK type 8 test segment';
%
%         %
%         % Open a new SPK file.
%         %
%         handle = cspice_spkopn( SPK8, segid, 4 );
%
%         step   = DISCRETEEPOCHS(2) - DISCRETEEPOCHS(1);
%
%         %
%         % Create a type 8 segment.
%         %
%         cspice_spkw08( handle,                       ...
%                        BODY,                         ...
%                        CENTER,                       ...
%                        REF,                          ...
%                        DISCRETEEPOCHS(1),            ...
%                        DISCRETEEPOCHS(N_DISCRETE),   ...
%                        segid,                        ...
%                        POLY_DEG,                     ...
%                        DISCRETESTATES,               ...
%                        DISCRETEEPOCHS(1),            ...
%                        step )
%
%         %
%         % Close the SPK file.
%         %
%         cspice_spkcls( handle )
%
%
%      When this program is executed, no output is presented on
%      screen. After run completion, a new SPK type 8 exists in
%      the output directory.
%
%-Particulars
%
%   Close the SPK file attached to `handle'. The close operation tests the
%   file to ensure the presence of data segments.
%
%   A cspice_spkcls call should balance each call to cspice_spkopn.
%
%-Exceptions
%
%   1)  If there are no segments in the file, the error
%       SPICE(NOSEGMENTSFOUND) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `handle' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   See argument `handle'.
%
%-Restrictions
%
%   None.
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
%   -Mice Version 1.1.0, 20-JUL-2020 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       complete code example, based on the cspice_spkw08 example.
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
%   -Mice Version 1.0.0, 23-MAY-2012 (EDW)
%
%-Index_Entries
%
%   close an SPK file
%
%-&

function cspice_spkcls( handle)

   switch nargin
      case 1

         handle = zzmice_int( handle );

      otherwise

         error ( 'Usage: cspice_spkcls(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'spkcls_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end


