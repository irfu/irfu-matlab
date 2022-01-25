%-Abstract
%
%   CSPICE_CKOPN opens a new CK file, returning the handle
%   of the opened file.
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
%      fname    the name of the CK file to open.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%      ifname   a descriptive internal filename for the CK.
%
%               [1,c2] = size(ifname); char = class(ifname)
%
%                  or
%
%               [1,1] = size(ifname); cell = class(ifname)
%
%      ncomch   the number of characters to reserve for comments.
%
%               [1,1] = size(ncomch); int32 = class(ncomch)
%
%   the call:
%
%      [handle] = cspice_ckopn( name, ifname, ncomch )
%
%   returns:
%
%      handle   the file handle assigned to `fname'
%
%               [1,1] = size(handle); int32 = class(handle)
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
%   1) Create a CK type 3 segment; fill with data for a simple time
%      dependent rotation and angular velocity, and reserve room in
%      the CK comments area for 5000 characters.
%
%      Example code begins here.
%
%
%      function ckopn_ex1()
%
%         INST3      = -77703;
%         NCOMCH     = 5000;
%         REF        = 'J2000';
%         CK3        = 'ckopn_ex1.bc';
%         IFNAME     = 'Test CK type 3 created by cspice_ckw03';
%         SEGID3     = 'Test type 3 segment test CK';
%         SECPERTICK = 0.001;
%         SPACING    = 10.0;
%         MAXREC     = 50;
%
%         %
%         % Note, sclkdp is a vector input, not a vectorized scalar.
%         %
%         sclkdp    = [1:MAXREC]';
%         sclkdp    = (sclkdp - 1)*SPACING;
%
%         spinrate  = [1:MAXREC]*1.e-6;
%
%         theta     = [0:MAXREC-1]*SPACING;
%         theta     = theta .* spinrate;
%
%         %
%         % Create a zero-filled array for the angular velocity
%         % vectors. This allocates the needed memory and
%         % defines a variable of the correct shape.
%         %
%         expavvs = zeros( [3 MAXREC] );
%
%         a1 = zeros( [1 MAXREC] );
%         a2 = a1;
%
%         size(theta)
%         size(a2)
%         size(a1)
%         r  = cspice_eul2m( theta, a2, a1, 3, 1 ,3 );
%         q  = cspice_m2q( r );
%
%         %
%         % Fill the z component of the expavvs vectors with the
%         % corresponding spinrate element scaled to SECPERTICK.
%         %
%         expavvs(3,:) = spinrate/SECPERTICK;
%
%         begtim = sclkdp(1);
%         endtim = sclkdp(MAXREC);
%         avflag = 1;
%
%         starts = [1:(MAXREC/2)]';
%         starts = (starts-1)*2*SPACING;
%
%         %
%         % Open a new CK, write the data, catch any errors.
%         %
%         try
%            handle = cspice_ckopn( CK3, IFNAME, NCOMCH )
%            cspice_ckw03( handle,  ...
%                          begtim,  ...
%                          endtim,  ...
%                          INST3,   ...
%                          REF,     ...
%                          avflag,  ...
%                          SEGID3,  ...
%                          sclkdp,  ...
%                          q,       ...
%                          expavvs, ...
%                          starts )
%         catch
%
%            error( [ 'Failure: ' lasterr] )
%         end
%
%         cspice_ckcls(handle)
%
%
%      When this program is executed, no output is presented on
%      screen. After run completion, a new CK file exists in the
%      output directory.
%
%-Particulars
%
%   Open a new CK file, reserving room for comments if requested.
%
%   A cspice_ckcls call should balance every cspice_ckopn
%   call.
%
%-Exceptions
%
%   1)  If the value of `ncomch' is negative, a value of zero (0) will
%       be used for the number of comment characters to be set aside
%       for comments.
%
%   2)  If an error occurs while attempting to open a CK file the
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
%   CK.REQ
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
%   -Mice Version 1.1.0, 26-NOV-2021 (EDW) (JDR)
%
%       Updated the header to comply with NAIF standard. Added
%       complete code example based on existing fragment.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       extended -Particulars.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 29-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   open a new CK file
%
%-&

function [handle] = cspice_ckopn( fname, ifname, ncomch )

   switch nargin
      case 3

         fname  = zzmice_str(fname);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int( ncomch );

      otherwise

         error ( 'Usage: [handle] = cspice_ckopn(`fname`, `ifname`, ncomch)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'ckopn_c', fname, ifname, ncomch );
   catch spiceerr
      rethrow(spiceerr)
   end

