%-Abstract
%
%   CSPICE_VPRJPI calculates the vector in a specified plane that
%   maps under orthogonal projection to a specified vector in
%   another plane.
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
%      vin      an arbitrary 3-vector.
%
%               [3,1] = size(vin); double = class(vin)
%
%      projpl   a SPICE plane that represents the geometric plane containing
%               `vin'.
%
%               [1,1] = size(projpl); struct = class(projpl)
%
%               The structure has the fields:
%
%                 normal:   [3,1] = size(normal);   double = class(normal)
%                 constant: [1,1] = size(constant); double = class(constant)
%
%      invpl    a SPICE plane that represents the geometric plane containing
%               the inverse image of `vin' under orthogonal projection onto
%               `projpl'.
%
%               [1,1] = size(invpl); struct = class(invpl)
%
%               The structure has the fields:
%
%                 normal:   [3,1] = size(normal);   double = class(normal)
%                 constant: [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [vout, found] = cspice_vprjpi( vin, projpl, invpl )
%
%   returns:
%
%      vout     inverse orthogonal projection of `vin'.
%
%               [3,1] = size(vout); double = class(vout)
%
%               This is the vector lying in the plane `invpl' whose orthogonal
%               projection onto the plane `projpl' is `vin'. `vout' is valid
%               only when `found' is true. Otherwise, `vout' is undefined.
%
%      found    flag(s) indicating whether the inverse orthogonal projection
%               of `vin' could be computed.
%
%               [1,1] = size(found); logical = class(found)
%
%               `found' is true if so, false otherwise.
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
%   1) Suppose
%
%         vin    =  ( 0.0, 1.0, 0.0 ),
%
%      and that `projpl', the SPICE plane that represents the geometric
%      plane containing `vin', the has normal vector
%
%         projn  =  ( 0.0, 0.0, 1.0 ).
%
%
%      Also, let's suppose that `invpl' has normal vector and constant
%
%         invn   =  ( 0.0, 2.0, 2.0 )
%         invc   =    4.0.
%
%      Then `vin' lies on the y-axis in the x-y plane, and we want to
%      find the vector `vout' lying in `invpl' such that the orthogonal
%      projection of `vout' the x-y plane is `vin'. Let the notation
%      < a, b > indicate the inner product of vectors a and b.
%      Since every point x in `invpl' satisfies the equation
%
%         <  x,  (0.0, 2.0, 2.0)  >  =  4.0,
%
%      we can verify by inspection that the vector
%
%         ( 0.0, 1.0, 1.0 )
%
%      is in `invpl' and differs from `vin' by a multiple of `projn'. So
%
%         ( 0.0, 1.0, 1.0 )
%
%      must be `vout'.
%
%      The following code example is used to find this result
%      using Mice.
%
%
%      Example code begins here.
%
%
%      function vprjpi_ex1()
%
%         %
%         % Define a vector in plane1...
%         %
%         vin = [ 0., 1., 0.]';
%
%         %
%         % Construct 2 planes. Define the normal vectors for both
%         % planes and constant for the inverse plane.
%         %
%         norm1 = [ 0., 0., 1.]';
%         norm2 = [ 0., 2., 2.]';
%         con2  = 4.0;
%
%         %
%         % Create the SPICE planes
%         %
%         plane1 = cspice_nvp2pl( norm1, vin  );
%         plane2 = cspice_nvc2pl( norm2, con2 );
%
%         %
%         % Calculate the inverse projection to plane2.
%         %
%         [ vec_iproj, found] = cspice_vprjpi( vin, plane1, plane2);
%
%         if ( found )
%            disp( 'Found inverse vector:' )
%            fprintf('  %7.3f %7.3f %7.3f\n', vec_iproj )
%         else
%            disp( 'Could not find the inverse vector.' )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Found inverse vector:
%          0.000   1.000   1.000
%
%
%-Particulars
%
%   Projecting a vector orthogonally onto a plane can be thought of
%   as finding the closest vector in the plane to the original vector.
%   This "closest vector" always exists; it may be coincident with the
%   original vector. Inverting an orthogonal projection means finding
%   the vector in a specified plane whose orthogonal projection onto
%   a second specified plane is a specified vector. The vector whose
%   projection is the specified vector is the inverse projection of
%   the specified vector, also called the "inverse image under
%   orthogonal projection" of the specified vector. This routine
%   finds the inverse orthogonal projection of a vector onto a plane.
%
%   Related routines are cspice_vprjp, which projects a vector onto a plane
%   orthogonally, and cspice_vproj, which projects a vector onto another
%   vector orthogonally.
%
%-Exceptions
%
%   1)  If the normal vector of either input plane does not have unit
%       length (allowing for round-off error), the error
%       SPICE(NONUNITNORMAL) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the geometric planes defined by `projpl' and `invpl' are
%       orthogonal, or nearly so, the inverse orthogonal projection
%       of `vin' may be undefined or have magnitude too large to
%       represent with double precision numbers. In either such
%       case, `found' will be set to false.
%
%   3)  Even when `found' is true, `vout' may be a vector of extremely
%       large magnitude, perhaps so large that it is impractical to
%       compute with it. It's up to you to make sure that this
%       situation does not occur in your application of this routine.
%
%   4)  If any of the input arguments, `vin', `projpl' or `invpl', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `vin', `projpl' or `invpl', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  It is recommended that the input planes be created by one of
%       the Mice routines
%
%          cspice_nvc2pl ( Normal vector and constant to plane )
%          cspice_nvp2pl ( Normal vector and point to plane    )
%          cspice_psv2pl ( Point and spanning vectors to plane )
%
%       In any case each input plane must have a unit length normal
%       vector and a plane constant consistent with the normal
%       vector.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  G. Thomas and R. Finney, "Calculus and Analytic Geometry,"
%        7th Edition, Addison Wesley, 1988.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and modified code example to match
%       the input values used in the statement.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 12-MAR-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   vector projection onto plane
%
%-&

function [vout, found] = cspice_vprjpi( vin, projpl, invpl )

   switch nargin

      case 3

         vin    = zzmice_dp( vin );
         projpl = zzmice_pln( projpl );
         invpl  = zzmice_pln( invpl );

      otherwise

         error ( ['Usage: [vout(3), found] = ' ...
                                  'cspice_vprjpi( vin(3), projpl, invpl )'] )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [vout, found] = mice('vprjpi_c', vin, projpl, invpl );
   catch spiceerr
      rethrow(spiceerr)
   end

