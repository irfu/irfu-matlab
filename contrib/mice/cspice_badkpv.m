%-Abstract
%
%   CSPICE_BADKPV determines if a kernel pool variable is present and if so
%   that it has the correct size and type.
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
%      caller   the name of the routine calling this routine to check
%               correctness of kernel pool variables.
%
%               [1,c1] = size(caller); char = class(caller)
%
%                  or
%
%               [1,1] = size(caller); cell = class(caller)
%
%      name     the name of a kernel pool variable that the calling program
%               expects to be present in the kernel pool.
%
%               [1,c2] = size(name); char = class(name)
%
%                  or
%
%               [1,1] = size(name); cell = class(name)
%
%      comp     the comparison operator to use when comparing the number of
%               components of the kernel pool variable specified by `name'
%               with the integer `size'.
%
%               [1,c3] = size(comp); char = class(comp)
%
%                  or
%
%               [1,1] = size(comp); cell = class(comp)
%
%               If `dim' is is the actual size of the kernel pool variable
%               then cspice_badkpv will check that the sentence
%
%                  dim comp size
%
%               is a true statement. If it is not a true statement
%               an error will be signaled.
%
%               Allowed values for `comp' and their meanings are:
%
%                  '='      dim == size
%                  '<'      dim <  size
%                  '>'      dim >  size
%                  '=>'     dim >= size
%                  '<='     dim <= size
%
%      size     an integer to compare with the actual number of components of
%               the kernel pool variable specified by `name'.
%
%               [1,1] = size(size); int32 = class(size)
%
%      divby    an integer that is one of the factors of the actual dimension
%               of the specified kernel pool variable.
%
%               [1,1] = size(divby); int32 = class(divby)
%
%               In other words, it is expected that `divby' evenly divides
%               the actual dimension of `name'. In those cases in which the
%               factors of the dimension of `name' are not important, set
%               `divby' to 1 in the calling program.
%
%      type     the expected type of the kernel pool variable.
%
%               [1,c4] = size(type); char = class(type)
%
%                  or
%
%               [1,1] = size(type); cell = class(type)
%
%               Recognized values are
%
%                  'C' for character type
%                  'N' for numeric type (integer and double precision)
%
%               The case of `type' is insignificant. If the value
%               of `type' is not one of the 2 values given above
%               no check for the type of the variable will be
%               performed.
%
%   the call:
%
%      cspice_badkpv( caller, name, comp, size, divby, type )
%
%   signals a SPICE error if the kernel pool variable lacks the described
%   properties, otherwise the call has not effect.
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
%   1) Suppose that you need to fetch a number of variables
%      from the kernel pool and want to check that the requested
%      items are in fact available prior to performing further
%      computations. The code example shows how you might use
%      this routine to handle the details of checking of
%      the various items.
%
%      Although by default the SPICE error handling system will
%      report the error and halt the execution of the program, in
%      this example we have decided to change this behavior to
%      display the error messages and continue the execution of
%      the program.
%
%      Use the kernel shown below to define some variables related
%      to the Earth.
%
%
%         KPL/PCK
%
%         File name: badkpv_ex1.tpc
%
%         The contents of this kernel are not intended for
%         real applications. Use only with this example.
%
%         \begindata
%
%            BODY_399_DATA  = ( 3.1416, 2.71828, 0.5, 12.0 )
%            BODY_399_NAMES = ( 'PI', 'E', 'HALF', 'DOZEN' )
%
%         \begintext
%
%         End of constants kernel
%
%
%      Example code begins here.
%
%
%      function badkpv_ex1()
%
%         %
%         % Local parameters.
%         %
%         CALLER =   'BADKPV_EX1';
%
%         %
%         % Load the test kernel.
%         %
%         cspice_furnsh( 'badkpv_ex1.tpc' );
%
%         %
%         % Assume that we need some data for body 399 and we expect
%         % there to be an even number of items available and at
%         % least 4 such items. Moreover we expect these items to be
%         % numeric. Note that the variable assignments below are
%         % present only to assist in understanding the calls to
%         % cspice_badkpv.
%         %
%         name  = 'BODY_399_DATA';
%         comp  = '=>';
%         size  =  4;
%         divby =  2;
%         type  = 'N';
%
%         try
%
%            cspice_badkpv( CALLER, name, comp, size, divby, type );
%
%            fprintf( [ 'Expected form of variable %s found in kernel',    ...
%                       ' pool.\n' ], name                             )
%
%         catch
%
%            %
%            % Catch the error, return the error string to the user.
%            %
%            disp( '======================================================' )
%            disp( lasterr )
%            disp( '======================================================' )
%
%         end
%
%         %
%         % In addition we need the names given to these items.
%         % Improperly indicate the array has type numeric.
%         %
%         name  = 'BODY_399_NAMES';
%         comp  = '=>';
%         size  =  4;
%         divby =  1;
%         type  = 'N';
%
%         try
%
%            cspice_badkpv( CALLER, name, comp, size, divby, type );
%
%            fprintf( [ 'Expected form of variable %s found in kernel',    ...
%                       ' pool.\n' ], name                             )
%
%         catch
%
%            %
%            % Catch the error, return the error string to the user.
%            %
%            disp( '======================================================' )
%            disp( lasterr )
%            disp( '======================================================' )
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
%      Expected form of variable BODY_399_DATA found in kernel pool.
%      ======================================================
%      mice: SPICE(BADVARIABLETYPE): [badkpv_c->BADKPV] BADKPV_EX1: The
%      kernel pool variable 'BODY_399_NAMES' must be of type "NUMERIC".
%      However, the current type is character. (CSPICE_N0066)
%      ======================================================
%
%
%      Note that, as expected, the error SPICE(BADVARIABLETYPE) is
%      signaled by the second cspice_badkpv call, since we have improperly
%      indicated that the requested array is numeric, when actually
%      it is of character type.
%
%-Particulars
%
%   This routine takes care of routine checking that often needs
%   to be done by programs and routines that rely upon kernel
%   pool variables being present and having the correct attributes.
%
%   It checks for the presence of the kernel pool variable and
%   examines the type and dimension of the variable to make sure
%   they conform to the requirements of the calling routine.
%
%-Exceptions
%
%   1)  If the kernel pool variable specified by `name' is not present
%       in the kernel pool, the error SPICE(VARIABLENOTFOUND) is
%       signaled by a routine in the call tree of this routine.
%
%   2)  If the comparison operator specified by `comp' is unrecognized,
%       the error SPICE(UNKNOWNCOMPARE) is signaled by a routine in
%       the call tree of this routine.
%
%   3)  If the expected type of the kernel pool variable `type' is not
%       one of the supported types, the error SPICE(INVALIDTYPE) is
%       signaled by a routine in the call tree of this routine.
%
%   4)  If the comparison of the actual size of the kernel pool
%       variable with `size' is not satisfied, the error
%       SPICE(BADVARIABLESIZE) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the variable does not have the expected type, the error
%       SPICE(BADVARIABLETYPE) is signaled by a routine in the call
%       tree of this routine.
%
%   6)  If any of the input arguments, `caller', `name', `comp',
%       `size', `divby' or `type', is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   7)  If any of the input arguments, `caller', `name', `comp',
%       `size', `divby' or `type', is not of the expected type, or it
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
%   ERROR.REQ
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
%   -Mice Version 1.0.0, 05-SEP-2021 (JDR)
%
%-Index_Entries
%
%   Check the properties of a kernel pool variable
%
%-&
function cspice_badkpv( caller, name, comp, size, divby, type )

   switch nargin
      case 6

         caller = zzmice_str(caller);
         name   = zzmice_str(name);
         comp   = zzmice_str(comp);
         size   = zzmice_int(size);
         divby  = zzmice_int(divby);
         type   = zzmice_str(type);

      otherwise

         error ( [ 'Usage: '                                                ...
                   'cspice_badkpv( `caller`, `name`, `comp`, size, divby, ' ...
                   '`type` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice('badkpv_c', caller, name, comp, size, divby, type);
   catch spiceerr
      rethrow(spiceerr)
   end
