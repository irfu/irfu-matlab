%-Abstract
%
%   CSPICE_EXPOOL confirms the existence of a numeric kernel variable in the
%   kernel pool.
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
%      name     the name of the numeric kernel variable whose existence in
%               the kernel pool is to be checked.
%
%               [1,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,1] = size(name); cell = class(name)
%
%   the call:
%
%      [found] = cspice_expool( name )
%
%   returns:
%
%      found    true whenever the specified variable is included in the pool.
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
%   1) The following code example demonstrates how to use cspice_expool
%      to confirm the existence of numeric kernel pool variables.
%      In the example, we will look for different variables;
%      some of them numeric, some string valued and some not
%      present in the kernel pool.
%
%      Use the kernel shown below; an IK defining two keywords
%      used to provide data for an instrument with NAIF ID -999001.
%
%
%         KPL/IK
%
%         File name: expool_ex1.ti
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
%      function expool_ex1()
%
%         %
%         % Local parameters.
%         %
%         IKNAME =   'expool_ex1.ti';
%         NKPVNM =   3;
%
%         %
%         % Define the variable names
%         %
%         keywrd = { 'INS-999001_FREQ_RGB',                                ...
%                    'NOT_IN_THE_POOL',                                    ...
%                    'INS-999001_FREQ_UNITS' };
%
%         %
%         % Load the instrument kernel.
%         %
%         cspice_furnsh( IKNAME );
%
%         for i=1:NKPVNM
%
%            %
%            % Check if the variable is numeric and present
%            % in the kernel pool.
%            %
%            [found] = cspice_expool( keywrd(i) );
%
%            fprintf( 'Variable name: %s\n', char(keywrd(i)) )
%
%            if ( found )
%
%               fprintf( [ '   It is numeric and exists in the kernel',    ...
%                          ' pool.\n' ]                                 )
%
%            else
%
%               fprintf( [ '   Either it is not numeric or it is not in',  ...
%                          ' the kernel pool.\n' ]                        )
%
%            end
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
%      Variable name: INS-999001_FREQ_RGB
%         It is numeric and exists in the kernel pool.
%      Variable name: NOT_IN_THE_POOL
%         Either it is not numeric or it is not in the kernel pool.
%      Variable name: INS-999001_FREQ_UNITS
%         Either it is not numeric or it is not in the kernel pool.
%
%
%-Particulars
%
%   This routine determines whether or not a numeric kernel pool
%   variable exists. It does not detect the existence of
%   string valued kernel pool variables.
%
%   A better routine for determining the existence of numeric kernel
%   pool variables is the routine cspice_dtpool which determines the
%   existence, size and type of kernel pool variables.
%
%-Exceptions
%
%   1)  If the input argument `name' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   2)  If the input argument `name' is not of the expected type, or
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
%   -Mice Version 1.0.0, 07-SEP-2020 (JDR)
%
%-Index_Entries
%
%   CONFIRM the existence of a pooled numeric kernel variable
%
%-&
function [found] = cspice_expool( name )

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [found] = cspice_expool( `name` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [found] = mice('expool_c', name);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
