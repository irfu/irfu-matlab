%-Abstract
%
%   CSPICE_VUPACK unpacks three scalar components from a vector.
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
%      v        a double precision 3-dimensional vector.
%
%               [3,1] = size(v); double = class(v)
%
%   the call:
%
%      [x, y, z] = cspice_vupack( v )
%
%   returns:
%
%      x,
%      y,
%      z        the double precision scalar components of the vector `v'.
%
%               [1,1] = size(x); double = class(x)
%               [1,1] = size(y); double = class(y)
%               [1,1] = size(z); double = class(z)
%
%               The following equalities hold:
%
%                  x = v(1)
%                  y = v(2)
%                  z = v(3)
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
%   1) Suppose that you have an instrument kernel that provides,
%      within a single keyword, the three frequencies used by the
%      instrument, and that you want to use these frequencies
%      independently within your code.
%
%      The following code example demonstrates how to use cspice_vupack
%      to get these frequencies into independent scalar variables.
%
%      Use the kernel shown below, an IK defining the three
%      frequencies used by an instrument with NAIF ID -999001.
%
%
%         KPL/IK
%
%         File name: vupack_ex1.ti
%
%         The keyword below define the three frequencies used by a
%         hypothetical instrument (NAIF ID -999001). They correspond
%         to three filters: red, green and blue. Frequencies are
%         given in micrometers.
%
%         \begindata
%
%            INS-999001_FREQ_RGB   = (  0.65,  0.55, 0.475 )
%            INS-999001_FREQ_UNITS = ( 'MICROMETERS'       )
%
%         \begintext
%
%
%         End of IK
%
%
%      Example code begins here.
%
%
%      function vupack_ex1()
%
%         %
%         % Local parameters.
%         %
%         IKNAME = 'vupack_ex1.ti';
%         KEYWRD = 'INS-999001_FREQ_RGB';
%
%         %
%         % Load the instrument kernel.
%         %
%         cspice_furnsh( IKNAME );
%
%         %
%         % Get the frequency data from the kernel pool.
%         %
%         [ddata, found] = cspice_gdpool( KEYWRD, 1, 3 );
%
%         if ( found )
%
%            [red, green, blue] = cspice_vupack( ddata );
%            fprintf( 'Blue  (nm):  %5.2f\n', blue  * 1000.0 )
%            fprintf( 'Green (nm):  %5.2f\n', green * 1000.0 )
%            fprintf( 'Red   (nm):  %5.2f\n', red   * 1000.0 )
%
%         else
%
%            fprintf( 'No data found in the kernel pool for %s\n', KEYWRD )
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
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Blue  (nm):  475.00
%      Green (nm):  550.00
%      Red   (nm):  650.00
%
%
%-Particulars
%
%   Basically, this is just shorthand notation for the common
%   sequence
%
%      x = v(1)
%      y = v(2)
%      z = v(3)
%
%   The routine is useful largely for two reasons. First, it
%   reduces the chance that the programmer will make a "cut and
%   paste" mistake, like
%
%      x = v(1)
%      y = v(1)
%      z = v(1)
%
%   Second, it makes conversions between equivalent units simpler,
%   and clearer. For instance, the sequence
%
%      x = v(1) * cspice_rpd()
%      y = v(2) * cspice_rpd()
%      z = v(3) * cspice_rpd()
%
%   can be replaced by the (nearly) equivalent sequence
%
%      v         = cspice_rpd() * v;
%      [x, y, z] = cspice_vupack( v );
%
%-Exceptions
%
%   1)  If the input argument `v' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   2)  If the input argument `v' is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
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
%   MICE.REQ
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
%   -Mice Version 1.0.0, 07-SEP-2020 (JDR)
%
%-Index_Entries
%
%   unpack three scalar components from a vector
%
%-&
function [x, y, z] = cspice_vupack( v )

   switch nargin
      case 1

         v = zzmice_dp(v);

      otherwise

         error ( 'Usage: [x, y, z] = cspice_vupack( v(3) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [x, y, z] = mice('vupack_c', v);
   catch spiceerr
      rethrow(spiceerr)
   end
