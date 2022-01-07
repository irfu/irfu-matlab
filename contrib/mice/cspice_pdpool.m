%-Abstract
%
%   CSPICE_PDPOOL inserts double precision data into the kernel pool.
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
%      name     name of the kernel pool variable to associate with the values
%               supplied in the array `values'. `name' is restricted to a
%               length of 32 characters or less.
%
%               [1,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,1] = size(name); cell = class(name)
%
%      values   values to load into the kernel pool sub-system with the
%               assigned variable name `name'.
%
%               [n,1] = size(values); double = class(values)
%
%   the call:
%
%      cspice_pdpool( name, values )
%
%   returns:
%
%      Inserts the variable `name' into the kernel pool with values as
%      defined in `values'.
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
%   1) The following example code shows how a topocentric frame for a
%      point on the surface of the earth may be defined at run time using
%      cspice_pcpool, cspice_pdpool, and cspice_pipool. In this example,
%      the surface point is associated with the body code 300000. To
%      facilitate testing, the location of the surface point coincides
%      with that of the DSN station DSS-12; the reference frame MYTOPO
%      defined here coincides with the reference frame DSS-12_TOPO.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: pdpool_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                        Contents
%            ---------                        --------
%            earth_720101_070426.bpc          Earth historical
%                                             binary PCK
%            earth_topo_050714.tf             DSN station FK
%
%         \begindata
%
%         KERNELS_TO_LOAD = ( 'earth_720101_070426.bpc',
%                             'earth_topo_050714.tf'    )
%
%         \begintext
%
%         End of meta-kernel.
%
%
%      Example code begins here.
%
%
%      function pdpool_ex1()
%
%         angles  = [-243.1945102442646, -54.7000629043147, 180.0]';
%
%         et      = 0.0;
%
%         axes    = [3, 2, 3]';
%         center  = 300000;
%         frclass = 4;
%         frclsid = 1500000;
%         frcode  = 1500000;
%
%         %
%         % Define the MYTOPO reference frame.
%         %
%         %
%         cspice_pipool( 'FRAME_MYTOPO',             frcode   );
%         cspice_pcpool( 'FRAME_1500000_NAME',      'MYTOPO'  );
%         cspice_pipool( 'FRAME_1500000_CLASS',      frclass  );
%         cspice_pipool( 'FRAME_1500000_CLASS_ID',   frclsid  );
%         cspice_pipool( 'FRAME_1500000_CENTER',     center   );
%
%         cspice_pcpool( 'OBJECT_300000_FRAME',     'MYTOPO'  );
%
%         cspice_pcpool( 'TKFRAME_MYTOPO_RELATIVE', 'ITRF93'  );
%         cspice_pcpool( 'TKFRAME_MYTOPO_SPEC',     'ANGLES'  );
%         cspice_pcpool( 'TKFRAME_MYTOPO_UNITS',    'DEGREES' );
%         cspice_pipool( 'TKFRAME_MYTOPO_AXES',      axes     );
%         cspice_pdpool( 'TKFRAME_MYTOPO_ANGLES',    angles   );
%
%         %
%         % Load a high precision binary earth PCK. Also load a
%         % topocentric frame kernel for DSN stations. Use a meta-kernel
%         % for convenience.
%         %
%         cspice_furnsh( 'pdpool_ex1.tm' );
%
%         %
%         % Look up transformation from DSS-12_TOPO frame to MYTOPO frame.
%         % This transformation should differ by round-off error from
%         % the identity matrix.
%         %
%         [rmat] = cspice_pxform( 'DSS-12_TOPO', 'MYTOPO', et );
%
%         fprintf( '\n' )
%         fprintf(['DSS-12_TOPO to MYTOPO transformation at et', ...
%                  '  %14.5f :\n'], et )
%         fprintf( '\n' )
%         fprintf( '    %18.15f  %18.15f  %18.15f\n', ...
%                  rmat(1,1), rmat(2,1), rmat(3,1) )
%         fprintf( '    %18.15f  %18.15f  %18.15f\n', ...
%                  rmat(1,2), rmat(2,2), rmat(3,2) )
%         fprintf( '    %18.15f  %18.15f  %18.15f\n', ...
%                  rmat(1,3), rmat(2,3), rmat(3,3) )
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
%      DSS-12_TOPO to MYTOPO transformation at et         0.00000 :
%
%           1.000000000000000   0.000000000000000   0.000000000000000
%           0.000000000000000   1.000000000000000  -0.000000000000000
%           0.000000000000000  -0.000000000000000   1.000000000000000
%
%
%-Particulars
%
%   This routine provides a programmatic interface for inserting
%   data into the SPICE kernel pool without reading an external file.
%
%-Exceptions
%
%   1)  If `name' is already present in the kernel pool and there
%       is sufficient room to hold all values supplied in `values',
%       the old values associated with `name' will be overwritten.
%
%   2)  If there is not sufficient room to insert a new variable into
%       the kernel pool and `name' is not already present in the kernel
%       pool, an error is signaled by a routine in the call tree of
%       this routine.
%
%   3)  If there is not sufficient room to insert the values
%       associated with `name', the error SPICE(NOMOREROOM) is signaled
%       by a routine in the call tree of this routine.
%
%   4)  If the kernel pool variable name length exceeds its maximum
%       allowed length (see Kernel Required Reading, kernel.req), the
%       error SPICE(BADVARNAME) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If any of the input arguments, `name' or `values', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `name' or `values', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%   KERNEL.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.2.0, 26-NOV-2021 (EDW) (JDR)
%
%       Changed the input argument name "dvals" to "values" for consistency
%       with other routines.
%
%       Edited the header to comply with NAIF standard. Updated example
%       code to provide parallel version of the one in the CSPICE
%       pdpool_c header.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.1.1, 12-MAR-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Added mention of the length restriction on the kernel pool variable
%       name "name".
%
%   -Mice Version 1.1.0, 23-FEB-2009 (EDW)
%
%       Added zzmice_str call on input "name" to convert string cells to
%       character arrays if "name" has type string cells. Added proper
%       markers for usage string variable types.
%
%   -Mice Version 1.0.0, 24-JAN-2006 (EDW)
%
%-Index_Entries
%
%   Set the value of a d.p. kernel pool variable
%
%-&

function cspice_pdpool( name, values )

   switch nargin
      case 2

         name   = zzmice_str(name);
         values = zzmice_dp(values);

      otherwise

         error ( 'Usage: cspice_pdpool( `name`, values(n) )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('pdpool_c', name, values );
   catch spiceerr
      rethrow(spiceerr)
   end



