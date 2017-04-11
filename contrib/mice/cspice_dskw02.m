%-Abstract
%
%   CSPICE_DSKW02 writes a type 2 DSK segment to a DSK file.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      handle      is the DAS file handle associated with a DSK file.
%                  The file must be open for write access.
%
%                  [1,1] = size(handle); int32 = class(handle)
%
%      center      is the ID code of the body whose surface is described
%                  by the input plate model. `center' refers to an
%                  ephemeris object.
%
%                  [1,1] = size(center); int32 = class(center)
%
%      surfid      is the ID code of the surface described by the input
%                  plate model. Multiple surfaces (for example, surfaces
%                  having different resolutions) may be associated with a
%                  given body.
%
%                  [1,1] = size(surfid); int32 = class(surfid)
%
%      dclass      is the data class of the input data set. See the
%                  header file DSKMice02.m for values and meanings.
%
%                  [1,1] = size(dclass); int32 = class(dclass)
%
%      frame       is the name of the reference frame with respect to
%                  which the input data are expressed.
%
%                  [1,c1] = size(fname); char = class(fname)
%
%                     or
%
%                  [1,1] = size(fname); cell = class(fname)
%
%      corsys      is the coordinate system in which the spatial coverage
%                  of the input data is expressed. `corsys' is an integer
%                  code. See the header file DSKMice02.m for values and
%                  meanings.
%
%                  [1,1] = size(corsys); int32 = class(corsys)
%
%      corpar      is an array of parameters associated with the input
%                  coordinate system.
%
%                  For latitudinal and rectangular coordinates, `corpar'
%                  is ignored.
%
%                  For planetodetic coordinates, the contents of `corpar'
%                  are:
%
%                    Element         Contents
%                    ---------       -----------------------------------
%                    corpar(1)       Equatorial radius of reference
%                    corpar(2)       Flattening coefficient. The polar
%                                    radius of the reference spheroid
%                                    is given by
%
%                                        corpar(1) * ( 1 - corpar(2) )
%
%                    [2,1] = size(corpar); double = class(corpar)
%
%      mncor1,
%      mxcor1,
%      mncor2,
%      mxcor2,
%      mncor3,
%      mxcor3      are, respectively, lower and upper bounds of
%                  each of the coordinates of the input data, where the
%                  coordinate system is defined by `corsys' and `corpar'.
%                  These bounds define the region for which the segment
%                  provides data.
%
%                  [1,1] = size(mncor1); double = class(mncor1)
%                  [1,1] = size(mxcor1); double = class(mxcor1)
%                  [1,1] = size(mncor2); double = class(mncor2)
%                  [1,1] = size(mxcor2); double = class(mxcor2)
%                  [1,1] = size(mncor3); double = class(mncor3)
%                  [1,1] = size(mxcor3); double = class(mxcor3)
%
%                  Distance units are km. Angular units are radians.
%
%                  The interpretation of these bounds depends on the data
%                  class; see `dclass' above.
%
%                     Single-valued surfaces
%                     ----------------------
%
%                     If the segment has data class SPICE_DSK_SVFCLS (see
%                     DSKMice02.m), the segment defines a surface as a
%                     single-valued function of its domain coordinates:
%                     for example, it may define the radius of the
%                     surface as a function of planetocentric longitude
%                     and latitude.
%
%                     In this case, the input data must cover a
%                     rectangle in dimensions 1 and 2 of the input
%                     coordinate system: the set of points
%
%                        R = { (x,y): mncor1 < x < mxcor1;
%                                     mncor2 < y < mxcor2  }
%
%                     must be completely covered by the input data. In
%                     other words, for each point (x,y) of R, there must
%                     be some plate containing a point whose first two
%                     coordinates are (x,y).
%
%                     The plate set may extend beyond the coordinate
%                     range defined by the bounds on the domain.
%
%                     Normally for single-valued surfaces, `mncor3' and
%                     `mxcor3' are the minimum and maximum values of the
%                     function attained on the domain.
%
%
%                     General surfaces
%                     ----------------
%
%                     If the segment has data class SPICE_DSK_GENCLS (see
%                     DSKMice02.m), the segment simply contains a
%                     collection of plates: no guarantees are made about the
%                     topology of the surface. The coordinate bounds simply
%                     indicate the spatial region for which the segment
%                     provides data.
%
%                     Note that shapes of small bodies such as asteroids
%                     and comet nuclei may fall into the "general
%                     surface" category. Surface features such as cliffs,
%                     caves, and arches prevent a surface from being
%                     represented as a single-valued function of latitude
%                     and longitude, for example.
%
%
%                   Longitude interpretation and restrictions
%                   -----------------------------------------
%
%                   The following rules apply to longitudes provided in
%                   the arguments
%
%                       mncor1
%                       mxcor1
%
%                   All angles have units of radians. The tolerance
%                   SPICE_DSK_ANGMRG is used for the comparisons shown
%                   below.
%
%                      1) Longitudes must be in the range
%
%                         -2*pi  :  2*pi
%
%                         Values that are out of range by no more than
%                         SPICE_DSK_ANGMRG are bracketed to be in range.
%
%
%                      2) It is acceptable for the longitude bounds to be
%                         equal or out of order. If
%
%                             mxcor1 < mncor1
%                                    -
%
%                         then either `mxcor1' is treated by the DSK
%                         subsystem as though it were mxcor1 + 2*pi, or
%                         `mncor1' is treated as mncor1 - 2*pi: whichever
%                         shift puts the bounds in the allowed range is
%                         made.
%
%                         The input longitude bounds must not be equal.
%                         If the lower bound is greater than the upper
%                         bound, the difference between the bounds must
%                         not be an integer multiple of 2*pi.
%
%                         Aside from any small changes made to move the
%                         input values of `mncor1' or `mxcor1' into range,
%                         the values are stored in the DSK segment as is.
%
%
%                      3) `mxcor1' must not exceed `mncor1' by more than 2*pi.
%                         Values that are out of range by no more than
%                         SPICE_DSK_ANGMRG are bracketed to be in range.
%
%
%      first,
%      last        are the endpoints of the time interval over which this
%                  data set is applicable. These times are expressed as
%                  seconds past J2000 TDB.
%
%                  [1,1] = size(first); double = class(first)
%                  [1,1] = size(last);  double = class(last)
%
%      vrtces      is an array of coordinates of the vertices.
%                  The ith vertex occupies elements [1:3,i] of
%                  this array.
%
%                  [3,m] = size(vrtces); double = class(vrtces)
%
%      plates      is an array representing the plates of the model.
%                  The elements of `plates' are vertex indices. The vertex
%                  indices of the ith plate occupy elements [1:3,i] of
%                  this array.
%
%                  [3,n] = size(plates); int32 = class(plates)
%
%      spaixd,
%      spaixi      are, respectively, the double precision and integer
%                  components of the spatial index of the segment.
%
%                  [p,1] = size(spaixd); double = class(spaixd)
%                  [q,1] = size(spaixi); int32 = class(spaixi)
%
%                  It is strongly recommended that the helper routine
%                  cspice_dskmi2 be used to create these arrays. See the
%                  Examples section below.
%
%   the call:
%
%      cspice_dskw02( handle, ...
%                     center, ...
%                     surfid, ...
%                     dclass, ...
%                     frame,  ...
%                     corsys, ...
%                     corpar, ...
%                     mncor1, ...
%                     mxcor1, ...
%                     mncor2, ...
%                     mxcor2, ...
%                     mncor3, ...
%                     mxcor3, ...
%                     first,  ...
%                     last,   ...
%                     vrtces, ...
%                     plates, ...
%                     spaixd,  ...
%                     spaixi )
%
%   returns:
%
%      None. This routine operates by side effects.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Create a three-segment DSK file using plate model data for
%   Phobos. Use latitudinal, rectangular, and planetodetic
%   coordinates in the respective segments. This is not a
%   realistic example, but it serves to demonstrate use of
%   the supported coordinate systems.
%
%   For simplicity, use an existing DSK file to provide the
%   input plate and vertex data. The selected input file has one
%   segment.
%
%      function dskw02_t
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
%         %
%         % Close the input DSK.
%         %
%         cspice_dascls( inhan )
%         cspice_dskcls( handle, true )
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
%   This routine writes a type 2 segment to a DSK file that has been
%   opened for write access.
%
%   Users planning to create DSK files should consider whether the
%   SPICE DSK creation utility MKDSK may be suitable for their needs.
%
%   This routine is supported by the routines cspice_dskmi2 and cspice_dskrb2
%   cspice_dskmi2 simplifies use of this routine by creating the "spatial
%   index" arrays required as inputs by this routine. cspice_dskrb2 computes
%   bounds on the third coordinate of the input plate set.
%
%   Spatial Indexes
%   ===============
%
%   A spatial index is a group of data structures that facilitates
%   rapid high-level computations involving sets of plates. The data
%   structures created by this routine are aggregated into arrays
%   of type SpiceInt and type SpiceDouble.
%
%
%   Voxel grids
%   ===========
%
%   A key geometric computation---probably the most important, as it
%   serves as a foundation for other high-level computations---is
%   finding the intersection of a ray with the plate set. DSK type 2
%   segments use data structures called "voxel grids" as part of
%   their indexing mechanism. There is a "coarse grid": a box that
%   completely encloses a DSK type 2 segment's plate set, and which
%   is composed of identically-sized cubes called "coarse voxels."
%   Each coarse voxel in composed of smaller cubes called "fine
%   voxels." When the term "voxel" is used without qualification, it
%   refers to fine voxels.
%
%   Type 2 DSK segments contain data structures that associate plates
%   with the fine voxels intersected by those plates. These
%   structures enable the type 2 DSK software to rapidly find plates
%   in a given region of space.
%
%   Voxel scales
%   ============
%
%   There are two voxel scales:
%
%      - The coarse voxel scale is the integer ratio of the
%        edge length of a coarse voxel to the edge length of
%        a fine voxel
%
%      - The fine voxel scale is the double precision ratio
%        of the edge length of a fine voxel to the average
%        extent of the plates in the input plate set. "Extents"
%        of a plate are the absolute values of the differences
%        between the respective maximum and minimum X, Y, and Z
%        coordinates of the plate's vertices.
%
%   Voxel scales determine the resolution of the voxel grid.
%   Voxel scales must be chosen to satisfy size constraints and
%   provide reasonable plate lookup performance.
%
%   The following considerations apply to spatial indexes of
%   type 2 DSK segments:
%
%      1)  The maximum number of coarse voxels is fixed at
%          SPICE_DSK02_MAXCGR (declared in SpiceDSK.h).
%
%      2)  If there are too few fine voxels, the average number of
%          plates per fine voxel will be very large. This largely
%          negates the performance improvement afforded by having an
%          index. Also, the number of plates per voxel may exceed limits
%          imposed by DSK subroutines that use static arrays.
%
%      3)  If there are too many fine voxels, the average number of
%          voxels intersected by a given plate may be too large for all
%          the plate-voxel associations to be stored. In addition, the
%          time needed to examine the plate lists for each voxel
%          (including the empty ones) may become quite large, again
%          negating the value of the index.
%
%   In many cases, voxel scales yielding optimum performance must be
%   determined by experiment. However, the following heuristics can
%   provide reasonable starting values:
%
%      Let `np' be the number of plates. Let `fs' be the fine voxel
%      scale. Then a reasonable value of `fs' may be
%
%                 (0.25)
%         fs =  np       / 8.
%
%      In general, `fs' should not smaller than 1.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dskw02_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 04-FEB-2016, EDW (JPL), NJB (JPL), ML (JPL)
%
%-Index_Entries
%
%   write a type 2 dsk segment
%
%-&

function cspice_dskw02( handle, ...
                        center, ...
                        surfid, ...
                        dclass, ...
                        frame,  ...
                        corsys, ...
                        corpar, ...
                        mncor1, ...
                        mxcor1, ...
                        mncor2, ...
                        mxcor2, ...
                        mncor3, ...
                        mxcor3, ...
                        first,  ...
                        last,   ...
                        vrtces, ...
                        plates, ...
                        spaixd,  ...
                        spaixi )

   switch nargin

      case 19

         handle = zzmice_int(handle);
         center = zzmice_int(center);
         surfid = zzmice_int(surfid);
         dclass = zzmice_int(dclass);
         frame  = zzmice_str(frame);
         corsys = zzmice_int(corsys);
         corpar = zzmice_dp(corpar);
         mncor1 = zzmice_dp(mncor1);
         mxcor1 = zzmice_dp(mxcor1);
         mncor2 = zzmice_dp(mncor2);
         mxcor2 = zzmice_dp(mxcor2);
         mncor3 = zzmice_dp(mncor3);
         mxcor3 = zzmice_dp(mxcor3);
         first  = zzmice_dp(first);
         last   = zzmice_dp(last);
         vrtces = zzmice_dp(vrtces);
         plates = zzmice_int(plates);
         spaixd = zzmice_dp(spaixd);
         spaixi = zzmice_int(spaixi);

      otherwise

         error ( [ 'Usage: cspice_dskw02( handle, center, surfid, dclass, ' ...
                   '`frame`, corsys, corpar(2), mncor1, mxcor1, mncor2, '   ...
                   'mxcor2, mncor3, mxcor3, first, last, vrtces(3,m), '     ...
                   'plates(3,n), spaixd(p), spaixi(q) )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dskw02_c', handle, center, surfid, dclass, frame,  ...
                       corsys, corpar, mncor1, mxcor1, mncor2, ...
                       mxcor2, mncor3, mxcor3, first,  last,   ...
                       vrtces, plates, spaixd, spaixi  );
   catch
      rethrow(lasterror)
   end
