%-Abstract
%
%   CSPICE_STPOOL retrieves the nth string from a kernel pool variable, where
%   the string may be continued across several components of the kernel pool
%   variable.
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
%      item     the name of a kernel pool variable for which the caller wants
%               to retrieve a full (potentially continued) string component.
%
%               [1,c1] = size(item); char = class(item)
%
%                  or
%
%               [1,1] = size(item); cell = class(item)
%
%      nth      the number of the string to retrieve from the kernel pool.
%
%               [1,1] = size(nth); int32 = class(nth)
%
%               The range of `nth' is 1 to the number of full strings that
%               are present.
%
%      contin   a sequence of characters which (if they appear as the last
%               non-blank sequence of characters in a component of a value of
%               a kernel pool variable) act as a continuation marker: the
%               marker indicates that the string associated with the
%               component is continued into the next literal component of the
%               kernel pool variable.
%
%               [1,c2] = size(contin); char = class(contin)
%
%                  or
%
%               [1,1] = size(contin); cell = class(contin)
%
%               If `contin' is blank, all of the components of `item' will be
%               retrieved as a single string.
%
%   the call:
%
%      [nthstr, found] = cspice_stpool( item, nth, contin )
%
%   returns:
%
%      nthstr   the `nth' full string associated with the kernel pool
%               variable specified by `item'.
%
%               [1,c3] = size(nthstr); char = class(nthstr)
%
%      found    a logical variable indicating success of the request to
%               retrieve the `nth' string associated with `item'.
%
%               [1,1] = size(found); logical = class(found)
%
%               If an nth string exists, `found' will be set to true;
%               otherwise `found' will be set to false.
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
%   1) Retrieve the n'th string from a kernel pool variable, where the string
%      may be continued across several components of the kernel pool variable.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: stpool_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs.
%
%         This meta-kernel contains a single variable assigned to an
%         array of two character strings that are split over several
%         components of the variable, with '//' as continuation
%         marker.
%
%         \begindata
%
%         LONG_VAL = ( 'String 1: inserted into //'
%                      'the kernel pool using //'
%                      '3 components.'
%                      'String 2: split up as 2 //'
%                      'components of a kernel pool variable.' )
%
%         \begintext
%
%         End of meta-kernel.
%
%
%      Example code begins here.
%
%
%      function stpool_ex1()
%         %
%         % Load the meta-kernel kernel containing the variable
%         % assignment.
%         %
%         cspice_furnsh( 'stpool_ex1.tm' )
%
%         %
%         % Retrieve the `nth' entry for kernel pool variable
%         % 'LONG_VAL' to `nthstr'.
%         %
%         ITEM   = 'LONG_VAL';
%         CONTIN = '//';
%
%         for nth=1:3
%
%            [nthstr, found] = cspice_stpool( ITEM, nth, CONTIN );
%
%            if ( found )
%
%               fprintf( ['Found index = %d component of kernel ' ...
%                         'variable %s. String:\n\n'], ...
%                         nth, ITEM)
%               fprintf( ' ``%s``\n\n', nthstr )
%
%            else
%
%               fprintf( ['No index = %d component of kernel ' ...
%                         'variable %s found \n'               ...
%                         'in the kernel pool.\n'], ...
%                          nth, ITEM)
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Found index = 1 component of kernel variable LONG_VAL. String:
%
%       ``String 1: inserted into the kernel pool using 3 components.``
%
%      Found index = 2 component of kernel variable LONG_VAL. String:
%
%       ``String 2: split up as 2 components of a kernel pool variable.``
%
%      No index = 3 component of kernel variable LONG_VAL found
%      in the kernel pool.
%
%
%-Particulars
%
%   The SPICE Kernel Pool provides a very convenient interface
%   for supplying both numeric and textual data to user application
%   programs. However, any particular component of a character
%   valued component of a kernel pool variable is limited to 80
%   or fewer characters in length.
%
%   This routine allows you to overcome this limitation by
%   "continuing" a character component of a kernel pool variable.
%   To do this you need to select a continuation sequence
%   of characters and then insert this sequence as the last non-blank
%   set of characters that make up the portion of the component
%   that should be continued.
%
%   For example, you may decide to use the sequence '//' to indicate
%   that a string should be continued to the next component of
%   a kernel pool variable. Then set up the
%   kernel pool variable as shown below
%
%      LONG_STRINGS = ( 'This is part of the first component //'
%                       'that needs more than one line when //'
%                       'inserting it into the kernel pool.'
%                       'This is the second string that is split //'
%                       'up as several components of a kernel pool //'
%                       'variable.' )
%
%   When loaded into the kernel pool, the variable LONG_STRINGS
%   will have six literal components:
%
%      component (1) == 'This is part of the first component //'
%      component (2) == 'that needs more than one line when //'
%      component (3) == 'inserting it into the kernel pool.'
%      component (4) == 'This is the second string that is split //'
%      component (5) == 'up as several components of a kernel pool //'
%      component (6) == 'variable.'
%
%   These are the components that would be retrieved by the call
%
%      [component, found] = cspice_gcpool( 'LONG_STRINGS', 1, 6 );
%
%   However, using the routine cspice_stpool you can view the variable
%   LONG_STRINGS as having two long components.
%
%      strgna = [ 'This is part of the first component that '              ...
%                 'needs more than one line when inserting '               ...
%                 'it into the kernel pool. ' ];
%
%      strgnb = [ 'This is the second string that is split '               ...
%                 'up as several components of a kernel pool '             ...
%                 'variable. ' ];
%
%
%   These string components would be retrieved by the following two
%   calls.
%
%      [strgna, found] = cspice_stpool( 'LONG_STRINGS, 1, '//' );
%      [strgnb, found] = cspice_stpool( 'LONG_STRINGS, 2, '//' );
%
%-Exceptions
%
%   1)  If the variable specified by `item' is not present in the kernel
%       pool or is present but is not character valued, `nthstr' will be
%       returned as a blank and `found' will be set to false. In
%       particular if `nth' is less than 1, `nthstr' will be returned as a
%       blank and `found' will be false.
%
%   2)  If the variable specified has a blank string associated with
%       its `nth' full string, `nthstr' will be blank and `found' will be
%       set to true.
%
%   3)  If the continuation character is a blank, every component
%       of the variable specified by `item' will be inserted into
%       the output string.
%
%   4)  If the continuation character is blank, then a blank component
%       of a variable is treated as a component with no letters. For
%       example:
%
%          STRINGS = ( 'This is a variable'
%                      'with a blank'
%                      ' '
%                      'component.' )
%
%       Is equivalent to
%
%          STRINGS = ( 'This is a variable'
%                      'with a blank'
%                      'component.' )
%
%       from the point of view of cspice_stpool if `contin' is set to the
%       blank character.
%
%   5)  If any of the input arguments, `item', `nth' or `contin', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `item', `nth' or `contin', is
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.2.0, 01-NOV-2021 (EDW) (JDR)
%
%       Changed the output argument name "string" to "nthstr" for consistency
%       with other routines.
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel. Reformatted example's output to fit within
%       maximum line length for SPICE headers.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.1.1, 12-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.1.0, 10-MAY-2011 (EDW)
%
%       "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.0, 26-SEP-2007 (EDW)
%
%-Index_Entries
%
%   Retrieve a continued string value from the kernel pool
%
%-&

function [nthstr, found] = cspice_stpool( item, nth, contin )

   switch nargin
      case 3

         item   = zzmice_str(item);
         nth    = zzmice_int(nth);
         contin = zzmice_str(contin);

      otherwise

         error ( ['Usage: [`nthstr`, found] = '                            ...
                   'cspice_stpool( `item`, nth, `contin` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [nthstr, found] = mice( 'stpool_c', item, nth, contin );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end


