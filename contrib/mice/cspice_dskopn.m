%-Abstract
%
%   CSPICE_DSKOPN opens a new DSK file, returning the handle
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
%      fname    the name of the DSK file to open.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%      ifname   the descriptive internal filename for the DSK.
%
%               [1,c2] = size(ifname); char = class(ifname)
%
%                  or
%
%               [1,1] = size(ifname); cell = class(ifname)
%
%      ncomch   the scalar integer number of characters to
%               reserve for comments.
%
%               [1,1] = size(ncomch); int32 = class(ncomch)
%
%   the call:
%
%      handle = cspice_dskopn( name, ifname, ncomch )
%
%   returns:
%
%      handle   the file handle assigned to 'fname'.
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
%   1) Create a three-segment DSK file using plate model data for
%      Phobos. Use latitudinal, rectangular, and planetodetic
%      coordinates in the respective segments. This is not a
%      realistic example, but it serves to demonstrate use of
%      the supported coordinate systems.
%
%      Use the DSK kernel below to provide, for simplicity, the input
%      plate and vertex data. This file has one segment only.
%
%         phobos_3_3.bds
%
%
%      Example code begins here.
%
%
%      function dskopn_ex1()
%
%         %
%         % MiceUser globally defines DSK parameters.
%         % For more information, please see MiceDSK.m.
%         %
%         MiceUser
%
%         NSEG = 3;
%
%         cornam = { 'radius', 'Z-coordinate', 'Z-coordinate', 'altitude'};
%
%         %
%         % Assign names of input and output DSK files.
%         %
%         indsk = 'phobos_3_3.bds';
%         dsk   = 'phobos_3_3_3seg.bds';
%
%         if ( exist( dsk, 'file' ) == 2 )
%            delete( dsk )
%         end
%
%
%         %
%         % Open input DSK for read access; find first segment.
%         %
%         inhan           = cspice_dasopr( indsk );
%         [dladsc, found] = cspice_dlabfs( inhan );
%
%
%         %
%         % Fetch vertices and plates from input DSK file.
%         %
%         % Note that vertex and plate indices are 1-based.
%         %
%         disp( 'Reading input data...' )
%
%         vrtces = cspice_dskv02( inhan, dladsc, 1, SPICE_DSK02_MAXVRT );
%         plates = cspice_dskp02( inhan, dladsc, 1, SPICE_DSK02_MAXPLT );
%
%         disp( 'Done.' )
%
%
%         %
%         % Set input array sizes required by cspice_dskmi2.
%         %
%         voxpsz = SPICE_DSK02_MAXVXP;
%         voxlsz = SPICE_DSK02_MXNVLS;
%         worksz = SPICE_DSK02_MAXCEL;
%         spaisz = SPICE_DSK02_SPAISZ;
%         makvtl = true;
%
%         %
%         % Set fine and coarse voxel scales. (These usually
%         % need to determined by experimentation.)
%         %
%         finscl = 5.0;
%         corscl = 4;
%
%         %
%         % Open a new DSK file.
%         %
%         handle = cspice_dskopn( dsk, dsk, 0 );
%
%         for segno=1:NSEG
%
%            %
%            % Create spatial index. We won't generate a
%            % vertex-plate mapping, so we set the flag
%            % for creating this map to "false."
%            %
%            fprintf( 'Creating segment %d\n', segno )
%            fprintf( 'Creating spatial index...\n' )
%
%            [spaixd, spaixi] = cspice_dskmi2( vrtces, plates, finscl, ...
%                                              corscl, worksz, voxpsz, ...
%                                              voxlsz, makvtl,         ...
%                                              spaisz );
%
%            fprintf( 'Done.\n')
%
%            %
%            % Set up inputs describing segment attributes:
%            %
%            % - Central body: Phobos
%            % - Surface ID code: user's choice.
%            %   We use the segment number here.
%            % - Data class: general (arbitrary) shape
%            % - Body-fixed reference frame
%            % - Time coverage bounds (TBD)
%            %
%            center = 401;
%            surfid = segno;
%            dclass = SPICE_DSK_GENCLS;
%            frame  = 'IAU_PHOBOS';
%
%            first = -50. * cspice_jyear();
%            last  =  50. * cspice_jyear();
%
%            %
%            % Set the coordinate system and coordinate system
%            % bounds based on the segment index.
%            %
%            % Zero out the coordinate parameters to start.
%            %
%            corpar = zeros(SPICE_DSK_NSYPAR,1);
%
%            switch segno
%
%               case 1
%
%                  %
%                  % Use planetocentric latitudinal coordinates. Set
%                  % the longitude and latitude bounds.
%                  %
%                  corsys = SPICE_DSK_LATSYS;
%
%                  mncor1 = -cspice_pi();
%                  mxcor1 =  cspice_pi();
%                  mncor2 = -cspice_halfpi();
%                  mxcor2 =  cspice_halfpi();
%
%               case 2
%
%                  %
%                  % Use rectangular coordinates. Set the
%                  % X and Y bounds.
%                  %
%                  % The bounds shown here were derived from
%                  % the plate data. They lie slightly outside
%                  % of the range spanned by the plates.
%                  %
%                  corsys = SPICE_DSK_RECSYS;
%
%                  mncor1 = -1.3;
%                  mxcor1 =  1.31;
%                  mncor2 = -1.21;
%                  mxcor2 =  1.2;
%
%               case 3
%
%                  %
%                  % Set the coordinate system to planetodetic.
%                  %
%                  corsys    = SPICE_DSK_PDTSYS;
%
%                  mncor1    = -cspice_pi();
%                  mxcor1    =  cspice_pi();
%                  mncor2    = -cspice_halfpi();
%                  mxcor2    =  cspice_halfpi();
%
%                  %
%                  % We'll use equatorial and polar radii from
%                  % pck00010.tpc. These normally would be fetched
%                  % at run time, but for simplicity, we'll use
%                  % hard-coded values.
%                  %
%                  re        = 13.0;
%                  rp        =  9.1;
%                  f         = ( re - rp ) / re;
%
%                  corpar = [ re, f ]';
%
%               otherwise
%
%                  error( 'Mice(BUG)' )
%
%            end
%
%            %
%            % Compute plate model radius bounds.
%            %
%            fprintf( 'Computing %s bounds of plate set...\n', ...
%                                            char(cornam(corsys)) )
%
%            [mncor3, mxcor3] = cspice_dskrb2( vrtces, plates, ...
%                                              corsys, corpar );
%
%            fprintf ( 'Done.\n' )
%
%            %
%            % Write the segment to the file.
%            %
%            fprintf( 'Writing segment...\n' )
%
%            cspice_dskw02( handle, center, surfid, dclass, frame,   ...
%                           corsys, corpar, mncor1, mxcor1, mncor2,  ...
%                           mxcor2, mncor3, mxcor3, first,  last,    ...
%                           vrtces, plates, spaixd, spaixi        )
%
%         end
%
%         %
%         % Close the input DSK.
%         %
%         cspice_dascls( inhan )
%         cspice_dskcls( handle, true )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Reading input data...
%      Done.
%      Creating segment 1
%      Creating spatial index...
%      Done.
%      Computing radius bounds of plate set...
%      Done.
%      Writing segment...
%      Creating segment 2
%      Creating spatial index...
%      Done.
%      Computing Z-coordinate bounds of plate set...
%      Done.
%      Writing segment...
%      Creating segment 3
%      Creating spatial index...
%      Done.
%      Computing altitude bounds of plate set...
%      Done.
%      Writing segment...
%
%
%      Note that after run completion, a new DSK exists in the output
%      directory.
%
%-Particulars
%
%   A cspice_dskcls call should balance every cspice_dskopn
%   call.
%
%-Exceptions
%
%   1)  If the input filename is blank, an error is signaled by a
%       routine in the call tree of this routine. No file will be
%       created.
%
%   2)  If the specified file cannot be opened without exceeding the
%       maximum allowed number of open DAS files, an error is signaled
%       by a routine in the call tree of this routine. No file will be
%       created.
%
%   3)  If the file cannot be opened properly, an error is signaled by
%       a routine in the call tree of this routine. No file will be
%       created.
%
%   4)  If the initial records in the file cannot be written, an error
%       is signaled by a routine in the call tree of this routine. No
%       file will be created.
%
%   5)  If no logical units are available, an error is signaled by a
%       routine in the call tree of this routine. No file will be
%       created.
%
%   6)  If the internal file name contains nonprinting characters
%       (ASCII codes decimal 0-31 and 127-255), an error is signaled
%       by a routine in the call tree of this routine. No file will be
%       created.
%
%   7)  If the number of comment characters allocated `ncomch' is
%       negative, an error is signaled by a routine in the call
%       tree of this routine. No file will be created.
%
%   8)  If any of the input arguments, `fname', `ifname' or `ncomch',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   9)  If any of the input arguments, `fname', `ifname' or `ncomch',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See argument `fname'.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 02-FEB-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   open a new DSK file
%
%-&

function [handle] = cspice_dskopn( fname, ifname, ncomch )

   switch nargin
      case 3

         fname  = zzmice_str(fname);
         ifname = zzmice_str(ifname);
         ncomch = zzmice_int(ncomch);

      otherwise

         error ( ['Usage: [handle] = cspice_dskopn( `fname`, ',
                                          '`ifname`, ncomch )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'dskopn_c', fname, ifname, ncomch );
   catch spiceerr
      rethrow(spiceerr)
   end




