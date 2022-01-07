%-Abstract
%
%   CSPICE_SPKOPN creates a new SPK file, returning the handle of the opened
%   file.
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
%      fname    the name of the new SPK file to be created.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%      ifname   the internal filename for the SPK file that is being created.
%
%               [1,c2] = size(ifname); char = class(ifname)
%
%                  or
%
%               [1,1] = size(ifname); cell = class(ifname)
%
%               The internal filename may be up to 60 characters long. If
%               you do not have any conventions for tagging your files, an
%               internal filename of 'SPK_file' is perfectly acceptable. You
%               may also leave it blank if you like.
%
%      ncomch   the space, measured in characters, to be initially set aside
%               for the comment area when a new SPK file is opened.
%
%               [1,1] = size(ncomch); int32 = class(ncomch)
%
%               The amount of space actually set aside may be greater than
%               the amount requested, due to the manner in which comment
%               records are allocated in an SPK file. However, the amount of
%               space set aside for comments will always be at least the
%               amount that was requested.
%
%               The value of ncomch should be greater than or equal to
%               zero, i.e., 0 <= ncomch. A negative value, should one
%               occur, will be assumed to be zero.
%
%   the call:
%
%      [handle] = cspice_spkopn( fname, ifname, ncomch )
%
%   returns:
%
%      handle   the handle of the opened SPK file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               If an error occurs when opening the file, the value of this
%               variable should not be used, as it will not represent a valid
%               handle.
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
%
%      Example code begins here.
%
%
%      function spkopn_ex1()
%
%         %
%         % Define the segment identifier parameters.
%         %
%         BODY       = 3;
%         CENTER     = 10;
%         REF        = 'J2000';
%         POLY_DEG   = 3;
%         SPK8       = 'spkopn_ex1.bsp';
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
%   Open a new SPK file, reserving room for comments if requested.
%
%-Exceptions
%
%   1)  If the value of `ncomch' is negative, a value of zero (0) will
%       be used for the number of comment characters to be set aside
%       for comments.
%
%   2)  If an error occurs while attempting to open the SPK file, the
%       value of `handle' will not represent a valid file handle.
%
%   3)  If any of the input arguments, `fname', `ifname' or `ncomch',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   4)  If any of the input arguments, `fname', `ifname' or `ncomch',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See `fname' and `handle'.
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
%   -Mice Version 1.1.0, 05-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       complete code example, based on the cspice_spkw08 example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
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
%   open a new SPK file
%
%-&

function [handle] = cspice_spkopn( fname, ifname, ncomch )

   switch nargin
      case 3

         fname  = zzmice_str(fname);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int(ncomch);

      otherwise

         error ( 'Usage: [handle] = cspice_spkopn(`fname`, `ifname`, ncomch)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'spkopn_c', fname, ifname, ncomch );
   catch spiceerr
      rethrow(spiceerr)
   end




