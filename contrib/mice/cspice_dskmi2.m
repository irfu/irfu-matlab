%-Abstract
%
%   CSPICE_DSKMI2 makes spatial index for a DSK type 2 segment. The index is
%   returned as a pair of arrays, one of type integer and one of type
%   double-precision. These arrays are suitable for use with the DSK type 2
%   writer cspice_dskw02.
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
%      vrtces   an array of coordinates of the vertices.
%
%               [3,m] = size(vrtces); double = class(vrtces)
%
%               The Ith vertex occupies elements [1:3,I] of this array.
%
%      plates   an array representing the triangular plates of a
%               shape model.
%
%               [3,n] = size(plates); int32 = class(plates)
%
%               The vertex indices of the Ith plate occupy elements [1:3,I] of
%               this array.
%
%      finscl   the fine voxel scale.
%
%               [1,1] = size(finscl); double = class(finscl)
%
%               This scale determines the edge length of the cubical
%               voxels comprising the fine voxel grid: the edge length
%               `voxsiz' is approximately
%
%                  finscl * {average plate extent}
%
%               where the extents of a plate are the respective
%               differences between the maximum and minimum
%               coordinate values of the plate's vertices.
%
%               The relationship between `voxsiz' and the average plate
%               extent is approximate because the `voxsiz' is adjusted
%               so that each dimension of the fine voxel grid is an
%               integer multiple of the coarse voxel scale.
%
%               See the -Particulars section below for further
%               information on voxel scales.
%
%      corscl   the coarse voxel scale.
%
%               [1,1] = size(corscl); int32 = class(corscl)
%
%               This integer scale is the ratio of the edge length of coarse
%               voxels to that of fine voxels. The coarse scale must be large
%               enough so that the total number of coarse voxels does not
%               exceed SPICE_DSK02_MAXCGR (see MiceDSK.m).
%
%      worksz   the second dimension of the workspace array `work'.
%
%               [1,1] = size(worksz); int32 = class(worksz)
%
%               `worksz' must be at least as large as the greater of
%
%                  - the number of fine voxel-plate associations
%
%                    This number is equal to
%
%                       np * {average number of fine voxels
%                             intersected by each plate}
%
%                  - the number of vertex-plate associations, if
%                    the vertex-plate mapping is constructed.
%
%                    This number is equal to
%
%                       nv + ( 3 * np )
%
%      voxpsz   the size of the fine voxel-plate pointer array.
%
%               [1,1] = size(voxpsz); int32 = class(voxpsz)
%
%               This array maps fine voxels to lists of plates that
%               intersect those voxels. `voxpsz' must be at least as large as
%
%                        3
%                  corscl  * {number of non-empty coarse voxels}
%
%      voxlsz   the size of the fine voxel-plate list array.
%
%               [1,1] = size(voxlsz); int32 = class(voxlsz)
%
%               This array contains, for each non-empty fine voxel, the
%               count of plates that intersect that voxel and the IDs of
%               those plates. `voxlsz' must be at least as large as
%
%                       `np' * {average number of fine voxels
%                               intersected by each plate}
%
%                   +   {number of non-empty fine voxels}
%
%      makvtl   a logical flag that, when set to true, indicates that a
%               vertex-plate association list is to be constructed.
%
%               [1,1] = size(makvtl); logical = class(makvtl)
%
%               The amount of workspace that is needed may depend on
%               whether a vertex-plate association list is
%               constructed. When this list is constructed, the size
%               of the integer component of the spatial index is
%               increased by the size of the list and the size of a
%               vertex-plate pointer array; the total of these sizes
%               is
%
%                  ( 2 * nv ) + ( 3 * np )
%
%      spxisz   the declared size of the output array `spaixi'.
%
%               [1,1] = size(spxisz); int32 = class(spxisz)
%
%               This size must be at least as large as the sum of
%
%                  - the fixed-size part of the integer component of
%                    the index, which includes the coarse voxel grid;
%                    this value is
%
%                       SPICE_DSK02_IXIFIX
%
%                  - the size `voxpsz' of the voxel-plate pointer array
%
%                  - the size `voxlsz' of the voxel-plate association
%                    list
%
%               plus, if the vertex-plate association list is
%               constructed,
%
%                  - the size `nv' of the vertex-plate pointer array
%
%                  - the size of the vertex-plate association list;
%                    this size is
%
%                       nv + ( 3 * np )
%
%   the call:
%
%      [spaixd, spaixi] = cspice_dskmi2( vrtces, plates, finscl,           ...
%                                        corscl, worksz, voxpsz,           ...
%                                        voxlsz, makvtl, spaisz );
%
%   returns:
%
%      spaixd,
%      spaixi   respectively, the double precision and integer
%               components of the spatial index of the segment.
%
%               [p,1] = size(spaixd); double = class(spaixd)
%               [q,1] = size(spaixi); int32 = class(spaixi)
%
%               `spaixd' must be declared with size at least
%               SPICE_DSK02_IXDFIX.
%
%               `spaixi' must be declared with size at least `spxisz'.
%
%-Parameters
%
%   See the parameter definitions file
%
%      MiceDSK.m
%
%   for declarations of DSK data type 2 (plate model) parameters.
%
%   See the parameter definitions file
%
%      MiceDLA.m
%
%   for declarations of DLA descriptor sizes and documentation of the
%   contents of DLA descriptors.
%
%   See the parameter definitions file
%
%      MiceDSK.m
%
%   for declarations of DSK descriptor sizes and documentation of the
%   contents of DSK descriptors.
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
%      Use the DSK kernel below to provide, for simplicity, the
%      input plate and vertex data. The selected input file has one
%      segment.
%
%         phobos_3_3.bds
%
%
%      Example code begins here.
%
%
%      function dskmi2_ex1()
%
%         %
%         % MiceUser globally defines DSK parameters.
%         % For more information, please see MiceDSK.m.
%         %
%         MiceUser
%
%         NSEG = 3;
%
%         cornam = {'radius', 'Z-coordinate', 'Z-coordinate', 'altitude'};
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
%            [spaixd, spaixi] = cspice_dskmi2( vrtces, plates, finscl,     ...
%                                              corscl, worksz, voxpsz,     ...
%                                              voxlsz, makvtl, spaisz );
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
%            fprintf( 'Computing %s bounds of plate set...\n',             ...
%                                            char(cornam(corsys)) )
%
%            [mncor3, mxcor3] = cspice_dskrb2( vrtces, plates,             ...
%                                              corsys, corpar );
%
%            fprintf ( 'Done.\n' )
%
%            %
%            % Write the segment to the file.
%            %
%            fprintf( 'Writing segment...\n' )
%
%            cspice_dskw02( handle, center, surfid, dclass, frame,         ...
%                           corsys, corpar, mncor1, mxcor1, mncor2,        ...
%                           mxcor2, mncor3, mxcor3, first,  last,          ...
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
%   Users planning to create DSK files should consider whether the
%   SPICE DSK creation utility MKDSK may be suitable for their needs.
%
%   This routine supports use of the DSK type 2 segment writer cspice_dskw02
%   by creating the "spatial index" arrays required as inputs to that
%   routine.
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
%   -  The coarse voxel scale is the integer ratio of the
%      edge length of a coarse voxel to the edge length of
%      a fine voxel
%
%   -  The fine voxel scale is the double precision ratio
%      of the edge length of a fine voxel to the average
%      extent of the plates in the input plate set. "Extents"
%      of a plate are the absolute values of the differences
%      between the respective maximum and minimum X, Y, and Z
%      coordinates of the plate's vertices.
%
%   Voxel scales determine the resolution of the voxel grid.
%   Voxel scales must be chosen to satisfy size constraints and
%   provide reasonable plate lookup performance.
%
%   The following considerations apply to spatial indexes of
%   type 2 DSK segments:
%
%      1)  The maximum number of coarse voxels is fixed at
%          SPICE_DSK02_MAXCGR (declared in MiceDSK.m).
%
%      2)  If there are too few fine voxels, the average number of
%          plates per fine voxel will be very large. This largely
%          negates the performance improvement afforded by having an
%          index. Also, the number of plates per voxel may exceed
%          limits imposed by DSK subroutines that use static arrays.
%
%      3)  If there are too many fine voxels, the average number of
%          voxels intersected by a given plate may be too large for
%          all the plate-voxel associations to be stored. In
%          addition, the time needed to examine the plate lists for
%          each voxel (including the empty ones) may become quite
%          large, again negating the value of the index.
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
%-Exceptions
%
%   1)  If the fine voxel scale is non-positive, the error
%       SPICE(BADFINEVOXELSCALE) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the coarse voxel scale is less than 1, the error
%       SPICE(BADCOARSEVOXSCALE) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If `nv', the number of vertices, is less than 3 or greater
%       than SPICE_DSK02_MAXVRT, the error SPICE(BADVERTEXCOUNT) is
%       signaled by a routine in the call tree of this routine.
%
%   4)  If `np', the number of plates, is less than 1 or greater than
%       SPICE_DSK02_MAXPLT, the error SPICE(BADPLATECOUNT) is signaled
%       by a routine in the call tree of this routine.
%
%   5)  If the workspace size `worksz' is less than np+1, where `np' is
%       the number of plates, the error SPICE(WORKSPACETOOSMALL) is
%       signaled by a routine in the call tree of this routine. This
%       is merely a sanity check; normally the workspace will need to
%       be substantially larger than this reference value. See the
%       description of `worksz' in the header section -I/O
%       above.
%
%   6)  If the voxel-plate pointer array size `voxpsz' is less than 1,
%       the error SPICE(PTRARRAYTOOSMALL) is signaled by a routine in
%       the call tree of this routine. This is merely a sanity check;
%       normally this pointer array will need to be substantially
%       larger than this reference value. See the description of
%       `voxpsz' in the header section -I/O above.
%
%   7)  If the voxel-plate list array size `voxlsz' is less than np+1,
%       where `np' is the number of plates, the error
%       SPICE(PLATELISTTOOSMALL) is signaled by a routine in the call
%       tree of this routine. This is merely a sanity check; normally
%       this array will need to be substantially larger than this
%       reference value. See the description of `voxlsz' in the header
%       section -I/O above.
%
%   8)  If the size `spxisz' of the integer array `spaixi' is too small
%       to contain its constituent structures, where the sizes
%       of these structures are derived from the inputs
%
%           `nv', `np', `voxpsz', `voxlsz'
%
%       the error SPICE(INTINDEXTOOSMALL) is signaled by a routine in
%       the call tree of this routine.
%
%   9)  If there is insufficient room to create any of the data
%       structures contained in the spatial index, an error is
%       signaled by a routine in the call tree of this routine.
%
%   10) If any of the input arguments, `vrtces', `plates', `finscl',
%       `corscl', `worksz', `voxpsz', `voxlsz', `makvtl' or `spxisz',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   11) If any of the input arguments, `vrtces', `plates', `finscl',
%       `corscl', `worksz', `voxpsz', `voxlsz', `makvtl' or `spxisz',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
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
%   DAS.REQ
%   DSK.REQ
%   MICE.REQ
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
%       Edited the header to comply with NAIF standard.
%
%       Changed the output argument name "spaix" to "spaixi" to comply
%       with NAIF standard.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Added proper usage string. Added missing information to -I/O
%       descriptions.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 04-FEB-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   make spatial index for type 2 DSK segment
%
%-&

function [spaixd, spaixi] = cspice_dskmi2( vrtces, plates, finscl, corscl, ...
                                           worksz, voxpsz, voxlsz, makvtl, ...
                                           spxisz )


   switch nargin
      case 9

         vrtces = zzmice_dp(vrtces);
         plates = zzmice_int(plates);
         finscl = zzmice_dp(finscl);
         corscl = zzmice_int(corscl);
         worksz = zzmice_int(worksz);
         voxpsz = zzmice_int(voxpsz);
         voxlsz = zzmice_int(voxlsz);
         makvtl = zzmice_int(makvtl);
         spxisz = zzmice_int(spxisz);

      otherwise

         error ( ['Usage: [ spaixd(SPICE_DSK02_IXDFIX), spaixi(spxisz)] = ' ...
                  'cspice_dskmi2( vrtces(3,m), plates(3,n), '               ...
                  'finscl, corscl, worksz, voxpsz, voxlsz, makvtl, spxisz)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [spaixd, spaixi] = mice( 'dskmi2_c', vrtces, plates, finscl, corscl, ...
                                           worksz, voxpsz, voxlsz, makvtl, ...
                                           spxisz);

   catch spiceerr
      rethrow(spiceerr)
   end




