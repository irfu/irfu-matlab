%-Abstract
%
%   CSPICE_GETFVN returns the field-of-view (FOV) parameters for a specified
%   instrument. The instrument is specified by name.
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
%      inst     the name of an instrument, such as a spacecraft-mounted
%               framing camera, for which the field of view parameters are to
%               be retrieved from the kernel pool.
%
%               [1,c1] = size(inst); char = class(inst)
%
%                  or
%
%               [1,1] = size(inst); cell = class(inst)
%
%               `inst' is case-insensitive, and leading and trailing blanks
%               in `inst' are not significant. Optionally, you may supply the
%               integer ID for the instrument as an integer string. For
%               example, both 'CASSINI_ISS_NAC' and '-82360' are legitimate
%               strings that indicate the CASSINI ISS NAC camera is the
%               instrument of interest.
%
%      room     the maximum number of 3-dimensional vectors that can be
%               returned in `bounds'.
%
%               [1,1] = size(room); int32 = class(room)
%
%   the call:
%
%      [shape, frame, bsight, bounds] = cspice_getfvn( inst, room )
%
%   returns:
%
%      shape    a character string that describes the "shape" of the field of
%               view.
%
%               [1,c2] = size(shape); char = class(shape)
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
%               [1,c3] = size(frame); char = class(frame)
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
%      various shapes and sizes, to load the FOV definitions, and the
%      mapping between their NAMES and NAIF IDs.
%
%
%         KPL/IK
%
%         File name: getfvn_ex1.ti
%
%         The keywords below define a circular, 10-degree wide FOV
%         with the boresight along the +Z axis of the 'SC999_INST001'
%         frame for an instrument with ID -999001 using the
%         "angles"-class specification.
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
%            INS-999002_FOV_BOUNDARY_CORNERS = (
%                                    1.0,  0.0,        0.01745506,
%                                    1.0,  0.03492077, 0.0        )
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
%         The keywords below define a triangular FOV with the
%         boresight along the +Y axis of the 'SC999_INST004' frame
%         for an instrument with ID -999004 using the "corners"-class
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
%         The keywords below define the INSTRUMENT name to NAIF ID
%         mappings for this example. For convenience we will keep them
%         in the example's IK although they could be placed elsewhere
%         (normally on the mission's frames kernel)
%
%         \begindata
%            NAIF_BODY_NAME += ( 'SC999_INST001' )
%            NAIF_BODY_CODE += ( -999001 )
%
%            NAIF_BODY_NAME += ( 'SC999_INST002' )
%            NAIF_BODY_CODE += ( -999002 )
%
%            NAIF_BODY_NAME += ( 'SC999_INST003' )
%            NAIF_BODY_CODE += ( -999003 )
%
%            NAIF_BODY_NAME += ( 'SC999_INST004' )
%            NAIF_BODY_CODE += ( -999004 )
%         \begintext
%
%         End of IK
%
%
%      Example code begins here.
%
%
%      function getfvn_ex1()
%
%         %
%         % Local parameters
%         %
%         MAXBND =   4;
%         NUMINS =   4;
%
%         %
%         % Define the instrument IDs.
%         %
%         instnm = { 'SC999_INST001', 'SC999_INST002',                     ...
%                    'SC999_INST003', 'SC999_INST004' };
%
%         %
%         % Load the IK file.
%         %
%         cspice_furnsh( 'getfvn_ex1.ti' );
%
%         %
%         % For each instrument ...
%         %
%         fprintf( '--------------------------------------\n' )
%         for i=1:NUMINS
%
%            %
%            % ... fetch FOV parameters and ...
%            %
%            [shape,  frame,                                               ...
%             bsight, bounds] = cspice_getfvn( instnm(i), MAXBND );
%
%            %
%            % ... print them to the screen.
%            %
%            fprintf( 'Instrument NAME: %s\n', char(instnm(i)) )
%            fprintf( '    FOV shape: %s\n', shape )
%            fprintf( '    FOV frame: %s\n', frame )
%            fprintf( 'FOV boresight:  %11.8f %11.8f %11.8f\n', bsight )
%
%            fprintf( '  FOV corners:\n' )
%            [m,n] = size(bounds);
%            for j=1:n
%
%               fprintf( '                %11.8f %11.8f %11.8f\n',         ...
%                                                      bounds(:,j) )
%
%            end
%            fprintf( '--------------------------------------\n' )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      --------------------------------------
%      Instrument NAME: SC999_INST001
%          FOV shape: CIRCLE
%          FOV frame: SC999_INST001
%      FOV boresight:   0.00000000  0.00000000  1.00000000
%        FOV corners:
%                       0.08715574  0.00000000  0.99619470
%      --------------------------------------
%      Instrument NAME: SC999_INST002
%          FOV shape: ELLIPSE
%          FOV frame: SC999_INST002
%      FOV boresight:   1.00000000  0.00000000  0.00000000
%        FOV corners:
%                       1.00000000  0.00000000  0.01745506
%                       1.00000000  0.03492077  0.00000000
%      --------------------------------------
%      Instrument NAME: SC999_INST003
%          FOV shape: RECTANGLE
%          FOV frame: SC999_INST003
%      FOV boresight:   0.00000000  0.00000000  1.00000000
%        FOV corners:
%                       0.01047177  0.00174523  0.99994365
%                      -0.01047177  0.00174523  0.99994365
%                      -0.01047177 -0.00174523  0.99994365
%                       0.01047177 -0.00174523  0.99994365
%      --------------------------------------
%      Instrument NAME: SC999_INST004
%          FOV shape: POLYGON
%          FOV frame: SC999_INST004
%      FOV boresight:   0.00000000  1.00000000  0.00000000
%        FOV corners:
%                       0.00000000  0.80000000  0.50000000
%                       0.40000000  0.80000000 -0.20000000
%                      -0.40000000  0.80000000 -0.20000000
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
%   This routine is identical in function to the routine cspice_getfov except
%   that it allows you to refer to an instrument by name (via a
%   character string) instead of by its NAIF instrument ID.
%
%   Given the NAIF instrument name, (and having "loaded" the
%   instrument field of view description and instrument name to NAIF
%   ID mapping) this routine returns the boresight of the instrument,
%   the "shape" of the field of view, a collection of vectors
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
%   keywords below <INSTID> is replaced with the instrument ID which
%   corresponds to the `inst' name passed into the module):
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
%   instrument ID which corresponds to the `inst' name passed into the
%   module):
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
%   1)  If the name of the instrument cannot be translated to its NAIF
%       ID code, the error SPICE(IDCODENOTFOUND) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the frame associated with the instrument can not be found,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   3)  If the shape of the instrument field of view can not be found
%       in the kernel pool, an error is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If the FOV_SHAPE specified by the instrument kernel is not
%       one of the four values: 'CIRCLE', 'POLYGON', 'ELLIPSE', or
%       'RECTANGLE', an error is signaled by a routine in the call
%       tree of this routine. If the 'ANGLES' specification is used,
%       FOV_SHAPE must be one of the three values: 'CIRCLE',
%       'ELLIPSE', or 'RECTANGLE'.
%
%   5)  If the direction of the boresight cannot be located in the
%       kernel pool, an error is signaled by a routine in the call
%       tree of this routine.
%
%   6)  If the number of components for the boresight vector in the
%       kernel pool is not 3, or they are not numeric, an error is
%       signaled by a routine in the call tree of this routine.
%
%   7)  If the boresight vector is the zero vector, an error is 
%       is signaled by a routine in the call tree of this routine.
%
%   8)  If the 'ANGLES' specification is not present in the kernel
%       pool and the boundary vectors for the edge of the field of
%       view cannot be found in the kernel pool, an error is signaled
%       by a routine in the call tree of this routine.
%
%   9)  If there is insufficient room (as specified by the argument
%       `room') to return all of the vectors associated with the
%       boundary of the field of view, an error is signaled by a
%       routine in the call tree of this routine.
%
%   10) If the number of components of vectors making up the field of
%       view is not a multiple of 3, an error is signaled by a routine
%       in the call tree of this routine.
%
%   11) If the number of components of vectors making up the field of
%       view is not compatible with the shape specified for the field
%       of view, an error is signaled by a routine in the call tree of
%       this routine.
%
%   12) If the reference vector for the 'ANGLES' specification can not
%       be found in the kernel pool, an error is signaled by a routine
%       in the call tree of this routine.
%
%   13) If the reference vector stored in the kernel pool to support
%       the 'ANGLES' specification contains an incorrect number of
%       components, contains 3 character components, or is parallel to
%       the boresight, an error is signaled by a routine in the call
%       tree of this routine.
%
%   14) If the 'ANGLES' specification is present in the kernel pool
%       and the reference angle stored in the kernel pool to support
%       the 'ANGLES' specification is absent from the kernel pool, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   15) If the keyword that stores the angular units for the angles
%       used in the 'ANGLES' specification is absent from the kernel
%       pool, an error is signaled by a routine in the call tree of
%       this routine.
%
%   16) If the value used for the units in the 'ANGLES' specification
%       is not one of the supported angular units of cspice_convrt, an error
%       is signaled by a routine in the call tree of this routine.
%
%   17) If the keyword that stores the cross angle for the 'ANGLES'
%       specification is needed and is absent from the kernel pool, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   18) If the angles for the 'RECTANGLE'/'ANGLES' specification case
%       have cosines that are less than those stored in the parameter
%       MINCOS defined in the cspice_getfov routine, an error is signaled by
%       a routine in the call tree of this routine.
%
%   19) If the class specification contains something other than
%       'ANGLES' or 'CORNERS', an error is signaled by a routine in
%       the call tree of this routine.
%
%   20) In the event that the CLASS_SPEC keyword is absent from the
%       kernel pool for the instrument whose FOV is sought, this
%       module assumes the 'CORNERS' specification is to be utilized.
%
%   21) If any of the input arguments, `inst' or `room', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   22) If any of the input arguments, `inst' or `room', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPICE kernels must be loaded by the calling program
%   before this routine is called.
%
%   This routine relies upon having successfully loaded an instrument
%   kernel (IK file) via the routine cspice_furnsh prior to calling this
%   routine.
%
%   The name `inst' must be associated with an NAIF ID code, normally
%   through a frames kernel (FK file) or instrument kernel (IK file).
%
%   Kernel data are normally loaded via cspice_furnsh once per program run,
%   NOT every time this routine is called.
%
%-Restrictions
%
%   1)  This routine will not operate unless proper mapping between
%       the instrument name `inst' and a NAIF ID exists, an I-kernel for
%       that instrument ID has been loaded via a call to cspice_furnsh prior
%       to calling this routine and this IK contains the specification
%       for the instrument field of view consistent with the
%       expectations of this routine.
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
%
%-Version
%
%   -Mice Version 1.0.0, 17-DEC-2021 (JDR)
%
%-Index_Entries
%
%   return instrument's FOV parameters, using instrument name
%
%-&
function [shape, frame, bsight, bounds] = cspice_getfvn( inst, room )

   switch nargin
      case 2

         inst = zzmice_str(inst);
         room = zzmice_int(room);

      otherwise

         error ( [ 'Usage: [`shape`, `frame`, bsight(3), bounds(3,N)] = '   ...
                   'cspice_getfvn( `inst`, room )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [shape, frame, bsight, bounds] = mice('getfvn_c', inst, room);
   catch spiceerr
      rethrow(spiceerr)
   end
