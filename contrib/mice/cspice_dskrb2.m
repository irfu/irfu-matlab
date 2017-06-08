%-Abstract
%
%   CSPICE_DSKRB2 determine range bounds for a DSK plate set.
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
%      vrtces      is an array of coordinates of the vertices. The Ith
%                  vertex occupies elements [1:3,I] of this array.
%
%                  [3,m] = size(vrtces); double = class(vrtces)
%
%      plates      is a array representing the triangular plates of a
%                  shape model. The elements of `plates' are vertex
%                  indices; vertex indices are 1-based. The vertex
%                  indices of the Ith plate occupy elements [1:3,I] of
%                  this array.
%
%                  [3,n] = size(plates); int32 = class(plates)
%
%      corsys      is an integer parameter identifying the coordinate
%                  system in which the bounds are to be computed.
%
%                  [1,1] = size(corsys); int32 = class(corsys)
%
%                  The bounds apply to the third coordinate in each system:
%
%                     Latitudinal:           radius
%                     Planetodetic:          altitude
%                     Rectangular:           Z
%
%      corpar     is a array of parameters associated with the coordinate
%                 system.
%
%                  [2,1] = size(corpar); double = class(corpar)
%
%                 Currently the only supported system that has
%                 associated parameters is the planetodetic system. For
%                 planetodetic coordinates,
%
%                    corpar(1) is the equatorial radius
%
%                    corpar(2) is the flattening coefficient. Let `re' and
%                    `rp' represent, respectively, the equatorial and
%                    polar radii of the reference ellipsoid of the
%                    system. Then
%
%                        corpar(2) = ( re - rp ) / re
%
%   the call:
%
%      [mncor3, mxcor3] = cspice_dskrb2( vrtces, plates, corsys, corpar )
%
%   returns:
%
%      mncor3    is a lower bound on the range of the third coordinate
%                of the system identified by `corsys' and `corpar', taken
%                over all plates.
%
%                [1,1] = size(mncor3); double = class(mncor3)
%
%                For latitudinal and rectangular coordinates, `mncor3'
%                is the greatest lower bound of the third coordinate.
%
%                For planetodetic coordinates, `mncor3' is an
%                approximation: it is less than or equal to the greatest
%                lower bound.
%
%      mxcor3    is the least upper bound on the range of the third
%                coordinate of the system identified by `corsys' and
%                `corpar', taken over all plates.
%
%                [1,1] = size(mxcor3); double = class(mxcor3)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      function dskd02_t
%
%         %
%         % MiceUser globally defines DSK parameters.
%         % For more information, please see DSKMiceUser.m and
%         % DSKMice02.m.
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
%         indsk = '/kernels/gen/dsk/phobos_3_3.bds';
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
%            [mncor3, mxcor3] = cspice_dskrb2( vrtces, plates, corsys, corpar );
%
%            fprintf ( 'Done.\n' )
%
%            %
%            % Write the segment to the file.
%            %
%            fprintf( 'Writing segment...\n' )
%
%            cspice_dskw02( handle, ...
%                              center, ...
%                              surfid, ...
%                              dclass, ...
%                              frame,  ...
%                              corsys, ...
%                              corpar, ...
%                              mncor1, ...
%                              mxcor1, ...
%                              mncor2, ...
%                              mxcor2, ...
%                              mncor3, ...
%                              mxcor3, ...
%                              first,  ...
%                              last,   ...
%                              vrtces, ...
%                              plates, ...
%                              spaixd,  ...
%                              spaixi )
%
%         end
%
%         cspice_dascls( inhan )
%         cspice_dskcls( handle, true )
%
%         %
%         % Close the input DSK.
%         %
%         cspice_dascls( inhan )
%
%   MATLAB outputs:
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
%      After run completion, A DSK exists in the output directory.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dskrb2_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 04-FEB-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   compute range bounds for type 2 dsk segment
%
%-&

function [mncor3, mxcor3] = cspice_dskrb2( vrtces, plates, corsys, corpar )

   switch nargin
      case 4

         vrtces = zzmice_dp(vrtces);
         plates = zzmice_int(plates);
         corsys = zzmice_int(corsys);
         corpar = zzmice_dp(corpar);

      otherwise

         error ( ['Usage:  '] )

   end

   %
   % Call the MEX library.
   %
   try
      [mncor3, mxcor3] = mice( 'dskrb2_c', vrtces, plates, corsys, corpar );

   catch
      rethrow(lasterror)
   end




