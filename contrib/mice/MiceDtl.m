%-Abstract
%
%   MiceDtl.m declares DSK tolerance variables for use with
%   Mice DSK APIs.
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
%-Parameters
%
%   None.
%
%-Examples
%
%   Include these definitions by using the call:
%
%      MiceUser;
%
%   from the user's application code.
%
%-Particulars
%
%     This file contains declarations of tolerance and margin values
%     used by the DSK subsystem.
%
%     It is recommended that the default values defined in this file be
%     changed only by expert SPICE users.
%
%     The values declared in this file are accessible at run time
%     through the routines
%
%        cspice_dskgtl  {DSK, get tolerance value}
%        cspice_dskstl  {DSK, set tolerance value}
%
%-Exceptions
%
%   None.
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
%   None.
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
%   -Mice Version 1.1.0, 11-JUN-2020 (JDR)
%
%       Updated the header to comply with NAIF standard. Renamed the include
%       file from DSKtol.m to MiceDtl.m for consistency with other
%       include files.
%
%       Removed warning note for developers from the abstract section.
%
%       Updated the header for compliance with NAIF standard. Added
%       -Parameters, -Exceptions, -Files, -Restrictions, -Literature_References
%       and -Author_and_Institution sections.
%
%   -Mice Version 1.0.0, 10-MAR-2016 (NJB) (EDW)
%
%-Index_Entries
%
%   Include DSK Mice parameters
%
%-&


%
%     Parameter declarations
%     ======================
%
%     DSK type 2 plate expansion factor
%     ---------------------------------
%
%     The factor XFRACT is used to slightly expand plates read from DSK
%     type 2 segments in order to perform ray-plate intercept
%     computations.
%
%     This expansion is performed to prevent rays from passing through
%     a target object without any intersection being detected. Such
%     "false miss" conditions can occur due to round-off errors.
%
%     Plate expansion is done by computing the difference vectors
%     between a plate's vertices and the plate's centroid, scaling
%     those differences by (1 + XFRACT), then producing new vertices by
%     adding the scaled differences to the centroid. This process
%     doesn't affect the stored DSK data.
%
%     Plate expansion is also performed when surface points are mapped
%     to plates on which they lie, as is done for illumination angle
%     computations.
%
%     This parameter is user-adjustable.
%

      SPICE_DSK_XFRACT = 1.D-10;

%
%     The keyword for setting or retrieving this factor is
%

      SPICE_DSK_KEYXFR = 1;

%
%     Greedy segment selection factor
%     -------------------------------
%
%     The factor SEGEXP is used to slightly expand DSK segment
%     boundaries in order to select segments to consider for
%     ray-surface intercept computations. The effect of this factor is
%     to make the multi-segment intercept algorithm consider all
%     segments that are sufficiently close to the ray of interest, even
%     if the ray misses those segments.
%
%     This expansion is performed to prevent rays from passing through
%     a target object without any intersection being detected. Such
%     "false miss" conditions can occur due to round-off errors.
%
%     The exact way this parameter is used is dependent on the
%     coordinate system of the segment to which it applies, and the DSK
%     software implementation. This parameter may be changed in a
%     future version of SPICE.
%

      SPICE_DSK_SGREED = 1.D-8;

%
%     The keyword for setting or retrieving this factor is
%

      SPICE_DSK_KEYSGR = SPICE_DSK_KEYXFR + 1;

%
%     Segment pad margin
%     ------------------
%
%     The segment pad margin is a scale factor used to determine when a
%     point resulting from a ray-surface intercept computation, if
%     outside the segment's boundaries, is close enough to the segment
%     to be considered a valid result.
%
%     This margin is required in order to make DSK segment padding
%     (surface data extending slightly beyond the segment's coordinate
%     boundaries) usable: if a ray intersects the pad surface outside
%     the segment boundaries, the pad is useless if the intercept is
%     automatically rejected.
%
%     However, an excessively large value for this parameter is
%     detrimental, since a ray-surface intercept solution found "in" a
%     segment can supersede solutions in segments farther from the
%     ray's vertex. Solutions found outside of a segment thus can mask
%     solutions that are closer to the ray's vertex by as much as
%     the value of this margin, when applied to a segment's boundary
%     dimensions.
%

      SPICE_DSK_SGPADM = 1.D-10;

%
%     The keyword for setting or retrieving this factor is
%

      SPICE_DSK_KEYSPM = SPICE_DSK_KEYSGR + 1;

%
%     Surface-point membership margin
%     -------------------------------
%
%     The surface-point membership margin limits the distance
%     between a point and a surface to which the point is
%     considered to belong. The margin is a scale factor applied
%     to the size of the segment containing the surface.
%
%     This margin is used to map surface points to outward
%     normal vectors at those points.
%
%     If this margin is set to an excessively small value,
%     routines that make use of the surface-point mapping won't
%     work properly.
%

      SPICE_DSK_PTMEMM = 1.D-7;

%
%     The keyword for setting or retrieving this factor is
%

      SPICE_DSK_KEYPTM = SPICE_DSK_KEYSPM + 1;

%
%     Angular rounding margin
%     -----------------------
%
%     This margin specifies an amount by which angular values
%     may deviate from their proper ranges without a SPICE error
%     condition being signaled.
%
%     For example, if an input latitude exceeds pi/2 radians by a
%     positive amount less than this margin, the value is treated as
%     though it were pi/2 radians.
%
%     Units are radians.
%

      SPICE_DSK_ANGMRG = 1.D-12;

%
%     This parameter is not user-adjustable.
%
%     The keyword for retrieving this parameter is
%

      SPICE_DSK_KEYAMG = SPICE_DSK_KEYPTM + 1;


%
%     Longitude alias margin
%     ----------------------
%
%     This margin specifies an amount by which a longitude
%     value can be outside a given longitude range without
%     being considered eligible for transformation by
%     addition or subtraction of 2*pi radians.
%
%     A longitude value, when compared to the endpoints of
%     a longitude interval, will be considered to be equal
%     to an endpoint if the value is outside the interval
%     differs from that endpoint by a magnitude less than
%     the alias margin.
%
%
%     Units are radians.
%

      SPICE_DSK_LONALI = 1.D-12;

%
%     This parameter is not user-adjustable.
%
%     The keyword for retrieving this parameter is
%

      SPICE_DSK_KEYLAL = SPICE_DSK_KEYAMG + 1;

%
%     End of include file MiceDtl.m
%
