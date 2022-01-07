%-Abstract
%
%   CSPICE_GETFOV returns the field-of-view (FOV) parameters for a specified
%   instrument. The instrument is specified by its NAIF ID code.
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
%      instid   the NAIF ID of an instrument.
%
%               [1,1] = size(instid); int32 = class(instid)
%
%      room     the maximum number of 3-dimensional vectors that can be
%               returned in `bounds'.
%
%               [1,1] = size(room); int32 = class(room)
%
%   the call:
%
%      [shape, frame, bsight, bounds] = cspice_getfov( instid, room )
%
%   returns:
%
%      shape    a character string that describes the "shape" of the field of
%               view.
%
%               [1,c1] = size(shape); char = class(shape)
%
%               Possible values returned are:
%
%                  'POLYGON'
%                  'RECTANGLE'
%                  'CIRCLE'
%                  'ELLIPSE'
%
%               If the value of `shape' is 'POLYGON' the field of view
%               of the instrument is a pyramidal polyhedron. The
%               vertex of the pyramid is at the instrument focal
%               point. The rays along the edges of the pyramid are
%               parallel to the vectors returned in `bounds'.
%
%               If the value of `shape' is 'RECTANGLE' the field of view
%               of the instrument is a rectangular pyramid. The vertex
%               of the pyramid is at the instrument focal point. The
%               rays along the edges of the pyramid are parallel to
%               the vectors returned in `bounds'. Moreover, in this
%               case, the boresight points along the axis of symmetry
%               of the rectangular pyramid.
%
%               If the value of `shape' is 'CIRCLE' the field of view of
%               the instrument is a circular cone centered on the
%               boresight vector. The vertex of the cone is at the
%               instrument focal point. A single vector will be
%               returned in `bounds'. This vector will be parallel to a
%               ray that lies in the cone that makes up the boundary
%               of the field of view.
%
%               If the value of `shape' is 'ELLIPSE' the field of view
%               of the instrument is an elliptical cone with the
%               boresight vector as the axis of the cone. In this
%               case two vectors are returned in `bounds'. One of the
%               vectors returned in `bounds' points to the end of the
%               semi-major axis of a perpendicular cross section of
%               the elliptic cone. The other vector points to the end
%               of the semi-minor axis of a perpendicular cross
%               section of the cone.
%
%      frame    the name of the reference frame in which the field of view
%               boundary vectors are defined.
%
%               [1,c2] = size(frame); char = class(frame)
%
%      bsight   a vector representing the principal instrument view
%               direction that can be
%
%                  -  the central pixel view direction,
%                  -  the optical axis direction,
%                  -  the FOV geometric center view direction,
%                  -  an axis of the FOV frame,
%
%               or any other vector specified for this purpose
%               in the IK FOV definition.
%
%               [3,1] = size(bsight); double = class(bsight)
%
%               The length of `bsight' is not specified other than being
%               non-zero.
%
%      bounds   an array of vectors that point to the "corners" of the
%               instrument field of view.
%
%               [3,n] = size(bounds); double = class(bounds)
%
%               (See the discussion accompanying `shape' for an expansion
%               of the term "corner of the field of view.") Note that the
%               vectors returned in `bounds' are not necessarily unit
%               vectors. Their magnitudes will be as set in the IK (for
%               'CORNERS'-style FOV specifications) or the same as the
%               magnitude of the boresight (for 'ANGLES'-style FOV
%               specifications.)
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
%   1) Load an IK, fetch the parameters for each of the FOVs defined
%      within and print these parameters to the screen.
%
%      Use the kernel shown below, an IK defining four FOVs of
%      various shapes and sizes, to load the FOV definitions.
%
%
%         KPL/IK
%
%         File name: getfov_ex1.ti
%
%         The keywords below define a circular, 10-degree wide FOV with
%         the boresight along the +Z axis of the 'SC999_INST001' frame
%         for an instrument with ID -999001 using the "angles"-class
%         specification.
%
%         \begindata
%            INS-999001_FOV_CLASS_SPEC       = 'ANGLES'
%            INS-999001_FOV_SHAPE            = 'CIRCLE'
%            INS-999001_FOV_FRAME            = 'SC999_INST001'
%            INS-999001_BORESIGHT            = ( 0.0, 0.0, 1.0 )
%            INS-999001_FOV_REF_VECTOR       = ( 1.0, 0.0, 0.0 )
%            INS-999001_FOV_REF_ANGLE        = ( 5.0 )
%            INS-999001_FOV_ANGLE_UNITS      = ( 'DEGREES' )
%         \begintext
%
%         The keywords below define an elliptical FOV with 2- and
%         4-degree angular extents in the XZ and XY planes and the
%         boresight along the +X axis of the 'SC999_INST002' frame for
%         an instrument with ID -999002 using the "corners"-class
%         specification.
%
%         \begindata
%            INS-999002_FOV_SHAPE            = 'ELLIPSE'
%            INS-999002_FOV_FRAME            = 'SC999_INST002'
%            INS-999002_BORESIGHT            = ( 1.0, 0.0, 0.0 )
%            INS-999002_FOV_BOUNDARY_CORNERS = ( 1.0, 0.0, 0.01745506,
%                                             1.0, 0.03492077, 0.0 )
%         \begintext
%
%         The keywords below define a rectangular FOV with 1.2- and
%         0.2-degree angular extents in the ZX and ZY planes and the
%         boresight along the +Z axis of the 'SC999_INST003' frame for
%         an instrument with ID -999003 using the "angles"-class
%         specification.
%
%         \begindata
%            INS-999003_FOV_CLASS_SPEC       = 'ANGLES'
%            INS-999003_FOV_SHAPE            = 'RECTANGLE'
%            INS-999003_FOV_FRAME            = 'SC999_INST003'
%            INS-999003_BORESIGHT            = ( 0.0, 0.0, 1.0 )
%            INS-999003_FOV_REF_VECTOR       = ( 1.0, 0.0, 0.0 )
%            INS-999003_FOV_REF_ANGLE        = ( 0.6 )
%            INS-999003_FOV_CROSS_ANGLE      = ( 0.1 )
%            INS-999003_FOV_ANGLE_UNITS      = ( 'DEGREES' )
%         \begintext
%
%         The keywords below define a triangular FOV with the boresight
%         along the +Y axis of the 'SC999_INST004' frame for an
%         instrument with ID -999004 using the "corners"-class
%         specification.
%
%         \begindata
%            INS-999004_FOV_SHAPE            = 'POLYGON'
%            INS-999004_FOV_FRAME            = 'SC999_INST004'
%            INS-999004_BORESIGHT            = (  0.0,  1.0,  0.0 )
%            INS-999004_FOV_BOUNDARY_CORNERS = (  0.0,  0.8,  0.5,
%                                                 0.4,  0.8, -0.2,
%                                                -0.4,  0.8, -0.2 )
%         \begintext
%
%         End of IK
%
%
%      Example code begins here.
%
%
%      function getfov_ex1()
%
%         %
%         % Set maximum number of boundary vectors, number of
%         % instruments and instrument IDs.
%         %
%         MAXBND = 4;
%         NUMINS = 4;
%         insids = [ -999001, -999002, -999003, -999004 ];
%
%         %
%         % Load the IK file.
%         %
%         cspice_furnsh ( 'getfov_ex1.ti' );
%
%         %
%         % For each instrument ...
%         %
%         fprintf ( '--------------------------------------\n' );
%         for i = 1:NUMINS
%
%            %
%            % ... fetch FOV parameters and ...
%            %
%            [shape, frame, bsight, bounds] = ...
%                     cspice_getfov( insids(i), MAXBND );
%
%            %
%            % ... print them to the screen.
%            %
%            fprintf ( 'Instrument ID: %i\n', insids(i) );
%            fprintf ( '    FOV shape: %s\n', shape );
%            fprintf ( '    FOV frame: %s\n', frame );
%            fprintf ( 'FOV boresight: %f %f %f\n', ...
%                      bsight(1), bsight(2), bsight(3) );
%
%            fprintf ( '  FOV corners: \n' );
%            [m,n] = size(bounds);
%            for j= 1:n
%               fprintf ( '               %f %f %f\n', ...
%                         bounds(1,j), bounds(2,j), bounds(3,j) );
%            end
%
%            fprintf ( '--------------------------------------\n' );
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      --------------------------------------
%      Instrument ID: -999001
%          FOV shape: CIRCLE
%          FOV frame: SC999_INST001
%      FOV boresight: 0.000000 0.000000 1.000000
%        FOV corners:
%                     0.087156 0.000000 0.996195
%      --------------------------------------
%      Instrument ID: -999002
%          FOV shape: ELLIPSE
%          FOV frame: SC999_INST002
%      FOV boresight: 1.000000 0.000000 0.000000
%        FOV corners:
%                     1.000000 0.000000 0.017455
%                     1.000000 0.034921 0.000000
%      --------------------------------------
%      Instrument ID: -999003
%          FOV shape: RECTANGLE
%          FOV frame: SC999_INST003
%      FOV boresight: 0.000000 0.000000 1.000000
%        FOV corners:
%                     0.010472 0.001745 0.999944
%                     -0.010472 0.001745 0.999944
%                     -0.010472 -0.001745 0.999944
%                     0.010472 -0.001745 0.999944
%      --------------------------------------
%      Instrument ID: -999004
%          FOV shape: POLYGON
%          FOV frame: SC999_INST004
%      FOV boresight: 0.000000 1.000000 0.000000
%        FOV corners:
%                     0.000000 0.800000 0.500000
%                     0.400000 0.800000 -0.200000
%                     -0.400000 0.800000 -0.200000
%      --------------------------------------
%
%
%-Particulars
%
%   This routine provides a common interface for retrieving from the
%   kernel pool the geometric characteristics of an instrument field
%   of view for a wide variety of remote sensing instruments
%   across many different space missions.
%
%   Given the NAIF instrument ID, (and having "loaded" the
%   instrument field of view description via the routine cspice_furnsh)
%   this routine returns the boresight of the instrument, the
%   "shape" of the field of view, a collection of vectors
%   that point along the edges of the field of view, and the
%   name of the reference frame in which these vectors are defined.
%
%   Currently this routine supports two classes of specifications
%   for FOV definitions: "corners" and "angles".
%
%   The "corners" specification requires that the following keywords
%   defining the shape, boresight, boundary vectors, and reference
%   frame of the FOV be provided in one of the text kernel files
%   (normally an IK file) loaded into the kernel pool (in the
%   keywords below <INSTID> is replaced with the instrument ID as
%   passed into the module):
%
%      INS<INSTID>_FOV_CLASS_SPEC         must be set to 'CORNERS' or
%                                         omitted to indicate the
%                                         "corners"-class
%                                         specification.
%
%      INS<INSTID>_FOV_SHAPE              must be set to one of these
%                                         values:
%
%                                            'CIRCLE'
%                                            'ELLIPSE'
%                                            'RECTANGLE'
%                                            'POLYGON'
%
%      INS<INSTID>_FOV_FRAME              must contain the name of
%                                         the frame in which the
%                                         boresight and boundary
%                                         corner vectors are defined.
%
%      INS<INSTID>_BORESIGHT              must be set to a 3D vector
%                                         defining the boresight in
%                                         the FOV frame specified in
%                                         the FOV_FRAME keyword.
%
%      INS<INSTID>_FOV_BOUNDARY   or
%      INS<INSTID>_FOV_BOUNDARY_CORNERS   must be set to one (for
%                                         FOV_SHAPE = 'CIRCLE'), two
%                                         (for FOV_SHAPE =
%                                         'ELLIPSE'), four (for
%                                         FOV_SHAPE = 'RECTANGLE'),
%                                         or three or more (for
%                                         'POLYGON') 3D vectors
%                                         defining the corners of the
%                                         FOV in the FOV frame
%                                         specified in the FOV_FRAME
%                                         keyword. The vectors should
%                                         be listed in either
%                                         clockwise or
%                                         counterclockwise order.
%                                         This is required by some
%                                         SPICE routines that make
%                                         use of FOV specifications.
%
%   The "angles" specification requires the following keywords
%   defining the shape, boresight, reference vector, reference and
%   cross angular extents of the FOV be provided in one of the text
%   kernel files (normally an IK file) loaded into the kernel
%   pool (in the keywords below <INSTID> is replaced with the
%   instrument ID as passed into the module):
%
%      INS<INSTID>_FOV_CLASS_SPEC         must be set to 'ANGLES' to
%                                         indicate the "angles"-class
%                                         specification.
%
%      INS<INSTID>_FOV_SHAPE              must be set to one of these
%                                         values:
%
%                                            'CIRCLE'
%                                            'ELLIPSE'
%                                            'RECTANGLE'
%
%      INS<INSTID>_FOV_FRAME              must contain the name of
%                                         the frame in which the
%                                         boresight and the computed
%                                         boundary corner vectors are
%                                         defined.
%
%      INS<INSTID>_BORESIGHT              must be set to a 3D vector
%                                         defining the boresight in
%                                         the FOV frame specified in
%                                         the FOV_FRAME keyword.
%
%      INS<INSTID>_FOV_REF_VECTOR         must be set to a 3D vector
%                                         that together with the
%                                         boresight vector defines
%                                         the plane in which the
%                                         first angular extent of the
%                                         FOV specified in the
%                                         FOV_REF_ANGLE keyword is
%                                         measured.
%
%      INS<INSTID>_FOV_REF_ANGLE          must be set to the angle
%                                         that is 1/2 of the total
%                                         FOV angular extent in the
%                                         plane defined by the
%                                         boresight and the vector
%                                         specified in the
%                                         FOV_REF_VECTOR keyword. The
%                                         the FOV angular half-extents
%                                         are measured from the
%                                         boresight vector.
%
%      INS<INSTID>_FOV_CROSS_ANGLE        must be set to the angle
%                                         that is 1/2 of the total
%                                         FOV angular extent in the
%                                         plane containing the
%                                         boresight and perpendicular
%                                         to the plane defined by the
%                                         boresight and the vector
%                                         specified in the
%                                         FOV_REF_VECTOR keyword. The
%                                         the FOV angular half-extents
%                                         are measured from the
%                                         boresight vector. This
%                                         keyword is not required for
%                                         FOV_SHAPE = 'CIRCLE'.
%
%      INS<INSTID>_FOV_ANGLE_UNITS        must specify units for the
%                                         angles given in the
%                                         FOV_REF_ANGLE and
%                                         FOV_CROSS_ANGLE keywords.
%                                         Any angular units
%                                         recognized by cspice_convrt are
%                                         acceptable.
%
%   The INS<INSTID>_FOV_REF_ANGLE and INS<INSTID>_FOV_CROSS_ANGLE
%   keywords can have any values for the 'CIRCLE' and 'ELLIPSE'
%   FOV shapes but must satisfy the condition cos( angle ) > 0 for
%   the 'RECTANGLE' shape.
%
%   This routine is intended to be an intermediate level routine.
%   It is expected that users of this routine will be familiar
%   with the SPICE frames subsystem and will be comfortable writing
%   software to further manipulate the vectors retrieved by this
%   routine.
%
%-Exceptions
%
%   1)  If the frame associated with the instrument can not be found,
%       the error SPICE(FRAMEMISSING) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If the shape of the instrument field of view can not be found
%       in the kernel pool, the error SPICE(SHAPEMISSING) is signaled
%       by a routine in the call tree of this routine signaled.
%
%   3)  If the FOV_SHAPE specified by the instrument kernel is not one
%       of the four values: 'CIRCLE', 'POLYGON', 'ELLIPSE', or
%       'RECTANGLE', the error SPICE(SHAPENOTSUPPORTED) is signaled by
%       a routine in the call tree of this routine. If the 'ANGLES'
%       specification is used, FOV_SHAPE must be one of the three
%       values: 'CIRCLE', 'ELLIPSE', or 'RECTANGLE'.
%
%   4)  If the direction of the boresight cannot be located in the
%       kernel pool, the error SPICE(BORESIGHTMISSING) is signaled by
%       a routine in the call tree of this routine.
%
%   5)  If the number of components for the boresight vector in the
%       kernel pool is not 3, or they are not numeric, the error
%       SPICE(BADBORESIGHTSPEC) is signaled by a routine in the call
%       tree of this routine.
%
%   6)  If the boresight vector is the zero vector, the error 
%       SPICE(ZEROBORESIGHT) is signaled by a routine in the call
%       tree of this routine.
%
%   7)  If the 'ANGLES' specification is not present in the kernel
%       pool and the boundary vectors for the edge of the field of
%       view cannot be found in the kernel pool, the error
%       SPICE(BOUNDARYMISSING) is signaled by a routine in the call
%       tree of this routine.
%
%   8)  If there is insufficient room (as specified by the argument
%       `room') to return all of the vectors associated with the
%       boundary of the field of view, the error SPICE(BOUNDARYTOOBIG)
%       is signaled by a routine in the call tree of this routine.
%
%   9)  If the number of components of vectors making up the field of
%       view is not a multiple of 3, the error SPICE(BADBOUNDARY) is
%       signaled by a routine in the call tree of this routine.
%
%   10) If the number of components of vectors making up the field of
%       view is not compatible with the shape specified for the field
%       of view, the error SPICE(BADBOUNDARY) is signaled by a routine
%       in the call tree of this routine.
%
%   11) If the reference vector for the 'ANGLES' specification can not
%       be found in the kernel pool, the error SPICE(REFVECTORMISSING)
%       is signaled by a routine in the call tree of this routine.
%
%   12) If the reference vector stored in the kernel pool to support
%       the 'ANGLES' specification contains an incorrect number of
%       components, contains 3 character components, or is parallel to
%       the boresight, the error SPICE(BADREFVECTORSPEC) is signaled
%       by a routine in the call tree of this routine.
%
%   13) If the 'ANGLES' specification is present in the kernel pool
%       and the reference angle stored in the kernel pool to support
%       the 'ANGLES' specification is absent from the kernel pool, the
%       error SPICE(REFANGLEMISSING) is signaled by a routine in the
%       call tree of this routine.
%
%   14) If the keyword that stores the angular units for the angles
%       used in the 'ANGLES' specification is absent from the kernel
%       pool, the error SPICE(UNITSMISSING) is signaled by a routine
%       in the call tree of this routine.
%
%   15) If the value used for the units in the 'ANGLES' specification
%       is not one of the supported angular units of cspice_convrt, an error
%       is signaled by a routine in the call tree of this routine.
%
%   16) If the keyword that stores the cross angle for the 'ANGLES'
%       specification is needed and is absent from the kernel pool,
%       the error SPICE(CROSSANGLEMISSING) is signaled by a routine in
%       the call tree of this routine.
%
%   17) If the angles for the 'RECTANGLE'/'ANGLES' specification case
%       have cosines that are less than those stored in the parameter
%       MINCOS, the error SPICE(BADBOUNDARY) is signaled by a routine
%       in the call tree of this routine.
%
%   18) If the class specification contains something other than
%       'ANGLES' or 'CORNERS', the error SPICE(UNSUPPORTEDSPEC) is
%       signaled by a routine in the call tree of this routine.
%
%   19) In the event that the CLASS_SPEC keyword is absent from the
%       kernel pool for the instrument whose FOV is sought, this
%       module assumes the 'CORNERS' specification is to be utilized.
%
%   20) If any of the input arguments, `instid' or `room', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   21) If any of the input arguments, `instid' or `room', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   This routine relies upon having successfully loaded an instrument
%   kernel (IK file) via the routine cspice_furnsh prior to calling this
%   routine.
%
%-Restrictions
%
%   1)  This routine will not operate unless an I-kernel for the
%       instrument with the NAIF ID specified in `instid' have been
%       loaded via a call to cspice_furnsh prior to calling this routine and
%       this IK contains the specification for the instrument field of
%       view consistent with the expectations of this routine.
%
%-Required_Reading
%
%   MICE.REQ
%   NAIF_IDS.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   B.V. Semenov        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 17-DEC-2021 (EDW) (JDR)
%
%       Bug fix: added missing exception for the boresight vector
%       being the zero vector.
%
%       Extended description of "bsight" output argument.
%
%       Edited header to comply with NAIF standard. Updated all sections to
%       improve the interface's documentation.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 12-MAR-2012 (EDW) (SCK)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.2, 24-APR-2010 (EDW)
%
%       Minor edit to code comments eliminating typo.
%
%   -Mice Version 1.0.1, 05-FEB-2009 (BVS)
%
%       Header update: added information about required IK keywords;
%       replaced old example with a new one more focused on getfov_c and
%       IK keywords.
%
%   -Mice Version 1.0.0, 07-DEC-2005 (EDW)
%
%-Index_Entries
%
%   return instrument's FOV parameters
%
%-&

function [shape, frame, bsight, bounds] = cspice_getfov( instid, room )

   switch nargin
      case 2

         instid = zzmice_int( instid );
         room   = zzmice_int( room   );

      otherwise

         error ( ['Usage: [`shape`, `frame`, bsight(3), bounds(3,N)] = ' ...
                          'cspice_getfov( instid, room )']  )

   end

   %
   % Call the MEX library.
   %
   try
      [shape, frame, bsight, bounds] = mice( 'getfov_c', instid, room );
   catch spiceerr
      rethrow(spiceerr)
   end
