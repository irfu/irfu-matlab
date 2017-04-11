%-Abstract
%
%   DSKMice02.pro declares DSK type 2-specific constants.
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
%   None.
%
%-Examples
%
%   Include the DSKMice02 definitions by the call:
%
%      DSKMice02
%
%   Please see the example program in cspice_dski02.
%
%-Particulars
%
%   CSPICE_DSKI02 keyword parameters
%   ================================
%
%   The following parameters may be supplied as actual values
%   of the `item' argument of CSPICE_DSKI02:
%
%
%   Name               Value   Description
%   ----               -----   ----------
%
%   SPICE_DSK02_KWNV     1     Number of vertices in model.
%
%   SPICE_DSK02_KWNP     2     Number of plates in model.
%
%   SPICE_DSK02_KWNVXT   3     Total number of voxels in fine grid.
%
%   SPICE_DSK02_KWVGRX   4     Voxel grid extent.  This extent is
%                              an array of three integers
%                              indicating the number of
%                              voxels in the X, Y, and Z
%                              directions in the fine voxel
%                              grid.
%
%   SPICE_DSK02_KWCGSC   5     Coarse voxel grid scale.  The extent
%                              of the fine voxel grid is related to
%                              the extent of the coarse voxel grid
%                              by this scale factor.
%
%   SPICE_DSK02_KWVXPS   6     Size of the voxel-to-plate pointer
%                              list.
%
%   SPICE_DSK02_KWVXLS   7     Voxel-plate correspondence list size.
%
%   SPICE_DSK02_KWVTLS   8     Vertex-plate correspondence list
%                              size.
%
%   SPICE_DSK02_KWPLAT   9     Plate array.  For each plate, this
%                              array contains the indices of the
%                              plate's three vertices.  The ordering
%                              of the array members is:
%
%                                 Plate 1 vertex index 1
%                                 Plate 1 vertex index 2
%                                 Plate 1 vertex index 3
%                                 Plate 2 vertex index 1
%                                 ...
%
%                              The vertex indices in this
%                              array start at 1 and end at
%                              NV, the number of vertices
%                              in the model.
%
%   SPICE_DSK02_KWVXPT   10    Voxel-plate pointer list. This list
%                              contains pointers that map fine
%                              voxels to lists of plates that
%                              intersect those voxels. Note that
%                              only fine voxels belonging to
%                              non-empty coarse voxels are in the
%                              domain of this mapping.
%
%   SPICE_DSK02_KWVXPL   11    Voxel-plate correspondence list.
%                              This list contains lists of plates
%                              that intersect fine voxels. (This
%                              list is the data structure into
%                              which the voxel-to-plate pointers
%                              point.) This list can contain
%                              empty lists. Plate IDs in this
%                              list start at 1 and end at NP,
%                              the number of plates in the model.
%
%   SPICE_DSK02_KWVTPT   12    Vertex-plate pointer list. This list
%                              contains pointers that map vertices
%                              to lists of plates to which those
%                              vertices belong.
%
%                              Note that the size of this list is
%                              always NV, the number of vertices.
%                              Hence there's no need for a separate
%                              keyword for the size of this list.
%
%   SPICE_DSK02_KWVTPL   13    Vertex-plate correspondence list.
%                              This list contains, for each vertex,
%                              the indices of the plates to which
%                              that  vertex belongs. Plate IDs in
%                              this list start at 1 and end at NP,
%                              the number of plates in the model.
%
%   SPICE_DSK02_KWCGPT   14    Coarse voxel grid pointers.  This is
%                              an array of pointers mapping coarse
%                              voxels to lists of pointers in the
%                              voxel-plate pointer list.  Each
%                              non-empty coarse voxel maps to a
%                              list of pointers; every fine voxel
%                              contained in a non-empty coarse voxel
%                              has its own pointers. Grid elements
%                              corresponding to empty coarse voxels
%                              contain non-positive values.
%
%
%   CSPICE_DSKD02 keyword parameters
%   ================================
%
%   The following parameters may be supplied as actual values
%   of the `item' argument of CSPICE_DSKD02:
%
%   Name               Value   Description
%   ----               -----   ----------
%
%   SPICE_DSK02_KWDSC   15     Array containing contents of Fortran
%                              DSK descriptor of segment. Note
%                              that DSK descriptors are not to be
%                              confused with DLA descriptors, which
%                              contain segment component base address
%                              and size information.  The dimension of
%                              this array is SPICE_DSK_DSCSIZ.
%
%   SPICE_DSK02_KWVTBD  16     Vertex bounds. This is an array of
%                              six values giving the minimum and
%                              maximum values of each component of the
%                              vertex set.
%
%   SPICE_DSK02_KWVXOR  17     Voxel grid origin. This is the location
%                              of the voxel grid origin in the
%                              body-fixed frame associated with the
%                              target body.
%
%   SPICE_DSK02_KWVXSZ  18     Voxel size.  DSK voxels are cubes; the
%                              edge length of each cube is given by the
%                              voxel size.  This size applies to the
%                              fine voxel grid. Units are km.
%
%   SPICE_DSK02_KWVERT  19     Vertex coordinates.
%
%
%-Required Reading
%
%   None.
%
%-Version
%
%   -Mice Version 1.0.0 20-OCT-2016, NJB (JPL), SCK (JPL), EDW (JPL)
%
%-Index_Entries
%
%   Include dsk type 2 parameters
%
%-&


%
% General parameters
%


%
% Dimension parameters
%

%
% Size of a SPICELIB DSK descriptor (in units of d.p. numbers):
%
   SPICE_DSK_DSCSIZ = 24;

%
% Number of coordinate system parameters in DSK descriptor:
%
   SPICE_DSK_NSYPAR = 10;

%
% Index parameters
%
%
% Fortran DSK descriptor index parameters. Matlab bases array indices
% at 1, as is proper:
%
   SPICE_DSK_SRFIDX = 1;
   SPICE_DSK_CTRIDX = 2;
   SPICE_DSK_CLSIDX = 3;
   SPICE_DSK_TYPIDX = 4;
   SPICE_DSK_FRMIDX = 5;
   SPICE_DSK_SYSIDX = 6;
   SPICE_DSK_PARIDX = 7;

%
% The offset between the indices immediately above and below
% is given by the parameter SPICE_DSK_NSYPAR. Literal values
% are used below for convenience of the reader.
%
   SPICE_DSK_MN1IDX = 17;
   SPICE_DSK_MX1IDX = 18;
   SPICE_DSK_MN2IDX = 19;
   SPICE_DSK_MX2IDX = 20;
   SPICE_DSK_MN3IDX = 21;
   SPICE_DSK_MX3IDX = 22;
   SPICE_DSK_BTMIDX = 23;
   SPICE_DSK_ETMIDX = 24;

%
% Data class parameters
%

%
% Spherical data class
%
   SPICE_DSK_SPHCLS = 1;


%
% Data class parameters
%

%
% Single-valued surface data class
%
   SPICE_DSK_SVFCLS = 1;

%
% General surface data class
%
   SPICE_DSK_GENCLS = 2;


%
% Coordinate system parameters
%


%
% Latitudinal coordinate system
%
   SPICE_DSK_LATSYS = 1;

%
% Cylindrical coordinate system
%
   SPICE_DSK_CYLSYS = 2;

%
% Rectangular coordinate system
%
   SPICE_DSK_RECSYS = 3;

%
% Planetodetic coordinate system
%
   SPICE_DSK_PDTSYS = 4;

%
% Type 2 definitions
% ==================
%
% Dimension parameters
%

%
% Maximum vertex count for single segment:
%
   SPICE_DSK02_MAXVRT =   16000002;

%
% Maximum plate count for single segment:
%
   SPICE_DSK02_MAXPLT =   ( 2 * (SPICE_DSK02_MAXVRT - 2 ) );

%
% The maximum allowed number of vertices, not taking into
% account shared vertices.
%
% Note that this value is not sufficient to create a vertex-plate
% mapping for a model of maximum plate count.
%
   SPICE_DSK02_MAXNPV =   ( 3 * (SPICE_DSK02_MAXPLT/2) + 1 );

%
% Maximum number of fine voxels.
%
   SPICE_DSK02_MAXVOX =   100000000;

%
% Maximum size of the coarse voxel grid array (in units of
% integers):
%
   SPICE_DSK02_MAXCGR =   100000;

%
% Maximum allowed number of vertex or plate
% neighbors a vertex may have.
%
   SPICE_DSK02_MAXEDG =   120;


%
% DSK type 2 spatial index parameters
% ===================================
%
%    DSK type 2 spatial index integer component
%    ------------------------------------------
%
%       +-----------------+
%       | VGREXT          |  (voxel grid extents, 3 integers)
%       +-----------------+
%       | CGRSCL          |  (coarse voxel grid scale, 1 integer)
%       +-----------------+
%       | VOXNPT          |  (size of voxel-plate pointer list)
%       +-----------------+
%       | VOXNPL          |  (size of voxel-plate list)
%       +-----------------+
%       | VTXNPL          |  (size of vertex-plate list)
%       +-----------------+
%       | CGRPTR          |  (coarse grid occupancy pointers)
%       +-----------------+
%       | VOXPTR          |  (voxel-plate pointer array)
%       +-----------------+
%       | VOXPLT          |  (voxel-plate list)
%       +-----------------+
%       | VTXPTR          |  (vertex-plate pointer array)
%       +-----------------+
%       | VTXPLT          |  (vertex-plate list)
%       +-----------------+
%

%
% Index parameters
%

%
% Grid extent index:
%
   SPICE_DSK02_SIVGRX =   1;

%
% Coarse grid scale index:
%
   SPICE_DSK02_SICGSC = ( SPICE_DSK02_SIVGRX + 3 );

%
% Voxel pointer count index:
%
   SPICE_DSK02_SIVXNP = ( SPICE_DSK02_SICGSC + 1 );

%
% Voxel-plate list count index:
%
   SPICE_DSK02_SIVXNL = ( SPICE_DSK02_SIVXNP + 1 );

%
% Vertex-plate list count index:
%
   SPICE_DSK02_SIVTNL =  ( SPICE_DSK02_SIVXNL + 1 );

%
% Coarse grid pointer array index:
%
   SPICE_DSK02_SICGRD =  ( SPICE_DSK02_SIVTNL + 1 );


%
% Spatial index integer component dimensions
%

%
% Size of fixed-size portion of integer component:
%
   SPICE_DSK02_IXIFIX = ( SPICE_DSK02_MAXCGR + 7 );


%
%
% DSK type 2 spatial index double precision component
% ---------------------------------------------------
%
%    +-----------------+
%    | Vertex bounds   |  6 values (min/max for each component)
%    +-----------------+
%    | Voxel origin    |  3 elements
%    +-----------------+
%    | Voxel size      |  1 element
%    +-----------------+
%


%
% Spatial index double precision indices
%

%
% Vertex bounds index:
%
   SPICE_DSK02_SIVTBD =   1;

%
% Voxel grid origin index:
%
   SPICE_DSK02_SIVXOR = ( SPICE_DSK02_SIVTBD + 6 );

%
% Voxel size index:
%
   SPICE_DSK02_SIVXSZ = ( SPICE_DSK02_SIVXOR + 3 );


%
% Spatial index double precision component dimensions
%

%
% Size of fixed-size portion of double precision component:
%
   SPICE_DSK02_IXDFIX =   10;

%
% Size of double precision component. This is a convenience
% parameter chosen to have a name consisent with the
% integer spatial index size.
%
   SPICE_DSK02_SPADSZ =   SPICE_DSK02_IXDFIX;

%
% The limits below are used to define a suggested maximum
% size for the integer component of the spatial index.
%

%
% Maximum number of entries in voxel-plate pointer array:
%
   SPICE_DSK02_MAXVXP = ( SPICE_DSK02_MAXPLT /2 );

%
% Maximum cell size:
%
   SPICE_DSK02_MAXCEL =   60000000;

%
% Maximum number of entries in voxel-plate list:
%
   SPICE_DSK02_MXNVLS =   SPICE_DSK02_MAXCEL +   ...
                 ( SPICE_DSK02_MAXVXP / 2 );

%
% Spatial index integer component size:
%
   SPICE_DSK02_SPAISZ = ( SPICE_DSK02_IXIFIX +   ...
                   SPICE_DSK02_MAXVXP +  ...
                   SPICE_DSK02_MXNVLS +   ...
                   SPICE_DSK02_MAXVRT +   ...
                   SPICE_DSK02_MAXNPV );


%
% Keyword parameters for integer data items:
%
   SPICE_DSK02_KWNV   = 1;
   SPICE_DSK02_KWNP   = SPICE_DSK02_KWNV   + 1;
   SPICE_DSK02_KWNVXT = SPICE_DSK02_KWNP   + 1;
   SPICE_DSK02_KWVGRX = SPICE_DSK02_KWNVXT + 1;
   SPICE_DSK02_KWCGSC = SPICE_DSK02_KWVGRX + 1;
   SPICE_DSK02_KWVXPS = SPICE_DSK02_KWCGSC + 1;
   SPICE_DSK02_KWVXLS = SPICE_DSK02_KWVXPS + 1;
   SPICE_DSK02_KWVTLS = SPICE_DSK02_KWVXLS + 1;
   SPICE_DSK02_KWPLAT = SPICE_DSK02_KWVTLS + 1;
   SPICE_DSK02_KWVXPT = SPICE_DSK02_KWPLAT + 1;
   SPICE_DSK02_KWVXPL = SPICE_DSK02_KWVXPT + 1;
   SPICE_DSK02_KWVTPT = SPICE_DSK02_KWVXPL + 1;
   SPICE_DSK02_KWVTPL = SPICE_DSK02_KWVTPT + 1;
   SPICE_DSK02_KWCGPT = SPICE_DSK02_KWVTPL + 1;

%
% Keyword parameters for double precision data items:
%
   SPICE_DSK02_KWDSC  = SPICE_DSK02_KWCGPT + 1;
   SPICE_DSK02_KWVTBD = SPICE_DSK02_KWDSC  + 1;
   SPICE_DSK02_KWVXOR = SPICE_DSK02_KWVTBD + 1;
   SPICE_DSK02_KWVXSZ = SPICE_DSK02_KWVXOR + 1;
   SPICE_DSK02_KWVERT = SPICE_DSK02_KWVXSZ + 1;

%
% End of parameter declaration file DSKMice02.pro
%
