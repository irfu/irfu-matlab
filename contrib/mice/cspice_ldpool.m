%-Abstract
%
%   CSPICE_LDPOOL loads the variables contained in a NAIF ASCII kernel file
%   into the kernel pool.
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
%      fname    the name of the text kernel file whose variables will be
%               loaded into the pool.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%   the call:
%
%      cspice_ldpool( fname )
%
%   returns:
%
%      None.
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
%   1) The following program demonstrates how to load the variables
%      contained in a NAIF ASCII kernel file into the kernel pool
%      and how to determine the properties of a stored kernel
%      variable.
%
%      The program prompts for text kernel name and for the name of
%      a kernel variable. If the variable is present in the kernel
%      pool, the dimension and type of the variable are displayed.
%
%
%      Example code begins here.
%
%
%      function ldpool_ex1()
%
%         %
%         % Prompt for the name of a text-kernel file.
%         %
%         fname = input( 'Enter text-kernel name        > ', 's' );
%
%         %
%         % Load the kernel. The same operation could be done using
%         % a cspice_furnsh call.
%         %
%         cspice_ldpool( fname );
%
%         varnam = input( 'Enter name of kernel variable > ', 's' );
%
%         [found, n, vtype] = cspice_dtpool( varnam );
%
%         if ( found )
%
%            fprintf( '\n' )
%            fprintf( 'Properties of variable %s:\n', varnam )
%            fprintf( '\n' )
%            fprintf( '   Size:   %d\n', n )
%
%            if ( strcmp( vtype, 'C' ) )
%
%               fprintf( '   Type:   Character\n' )
%
%            else
%
%               fprintf( '   Type:   Numeric\n' )
%
%            end
%
%         else
%
%            fprintf( '%s is not present in the kernel pool.\n', varnam )
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
%      platform, using the PCK file gm_de431.tpc to ask for the
%      variable 'BODY000_GMLIST', the output was:
%
%
%      Enter text-kernel name        > gm_de431.tpc
%      Enter name of kernel variable > BODY000_GMLIST
%
%      Properties of variable BODY000_GMLIST:
%
%         Size:   65
%         Type:   Numeric
%
%
%-Particulars
%
%   None.
%
%-Exceptions
%
%   1)  If an I/O error occurs while opening or reading a text kernel,
%       the error is signaled by a routine in the call tree of this
%       routine.
%
%   2)  If any text kernel parsing error occurs, the error is signaled
%       by a routine in the call tree of this routine.
%
%   3)  If a kernel pool overflow is detected, an error is signaled by
%       a routine in the call tree of this routine.
%
%   4)  If the input argument `fname' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   5)  If the input argument `fname' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   See `fname' in -I/O.
%
%-Restrictions
%
%   1)  Normally SPICE applications should load kernels via the
%       cspice_furnsh routine.
%
%-Required_Reading
%
%   KERNEL.REQ
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
%   -Mice Version 1.0.0, 30-JUN-2021 (JDR)
%
%-Index_Entries
%
%   LOAD variables from a text kernel file into the pool
%
%-&
function cspice_ldpool( fname )

   switch nargin
      case 1

         fname = zzmice_str(fname);

      otherwise

         error ( 'Usage: cspice_ldpool( `fname` )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('ldpool_c', fname);
   catch spiceerr
      rethrow(spiceerr)
   end
