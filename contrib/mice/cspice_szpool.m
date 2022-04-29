%-Abstract
%
%   CSPICE_SZPOOL returns the kernel pool size limitations.
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
%      name     the name of a kernel pool size parameter.
%
%               [1,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,1] = size(name); cell = class(name)
%
%               The following parameters may be specified.
%
%                  'MAXVAR'
%                  'MAXVAL'
%                  'MAXLIN'
%                  'MAXCHR'
%                  'MXNOTE'
%                  'MAXLEN'
%                  'MAXAGT'
%
%               See the main routine for a description of the
%               meaning of these parameters. Note that the case
%               of `name' is insignificant.
%
%   the call:
%
%      [n, found] = cspice_szpool( name )
%
%   returns:
%
%      n        the value of the parameter specified by `name'.
%
%               [1,1] = size(n); int32 = class(n)
%
%               If `name' is not one of the items specified above, `n' will
%               be returned with the value 0.
%
%      found    true if the parameter is recognized false if it is not.
%
%               [1,1] = size(found); logical = class(found)
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
%   1) The following code example demonstrates how to determine the
%      size of a kernel reader parameter.
%
%
%      Example code begins here.
%
%
%      function szpool_ex1()
%
%         %
%         % Local Variables
%         %
%         varname = 'MAXLEN';
%
%         %
%         % Make the call to retrieve the value of MAXLEN
%         %
%         [n, found] = cspice_szpool( varname );
%
%         %
%         % If MAXLEN parameter was found, print it out
%         %
%         if ( found )
%            fprintf( 'Kernel parameter found.\n' )
%            fprintf( 'value:  %s = %d\n', char(varname), n )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Kernel parameter found.
%      value:  MAXLEN = 32
%
%
%-Particulars
%
%   This routine provides the a programmatic interface to the
%   parameters used to define the kernel pool. It is not
%   anticipated that most kernel pool users will need to use this
%   routine.
%
%-Exceptions
%
%   1)  If the specified parameter is not recognized, the value of `n'
%       returned will be zero and `found' will be set to false.
%
%   2)  If the input argument `name' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `name' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
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
%   return a kernel pool definition parameter
%
%-&
function [n, found] = cspice_szpool( name )

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [n, found] = cspice_szpool( `name` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [n, found] = mice('szpool_c', name);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
