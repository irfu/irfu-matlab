%-Abstract
%
%   CSPICE_LMPOOL loads the variables contained in a text buffer
%   into the kernel pool.
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
%      cvals   string(s) defining SPICE kernel variable assignments
%              that could serve as a SPICE text kernel.
%
%              [n,c1] = size(cvals); char = class(cvals)
%
%                 or
%
%              [1,n] = size(cvals); cell = class(cvals)
%
%   the call:
%
%       cspice_lmpool( cvals)
%
%   inserts the variable assignments defined by 'cvals' into the
%   kernel pool subsystem. Once inserted, the user can access the
%   variables using the cspice_gcpool, cspice_gipool, or cspice_gdpool
%   calls.
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
%   1) Create a kernel in a text buffer and load the variables
%      contained within the buffer into the kernel pool. Ensure the
%      loaded data exists in the kernel pool. Query the pool for
%      each expected name, and print the size of the variable with
%      that name, and the type of data for that name.
%
%      Example code begins here.
%
%
%      function lmpool_ex1()
%
%         %
%         % Kernel pool variable's names.
%         %
%         lmpoolNames  = {                   ...
%                       'DELTET/DELTA_T_A',  ...
%                       'DELTET/K',          ...
%                       'DELTET/EB',         ...
%                       'DELTET/M',          ...
%                       'DELTET/DELTA_AT'    ...
%                        };
%
%         %
%         % Create a kernel in a text buffer.
%         %
%         textbuf = {                                             ...
%                'DELTET/DELTA_T_A = 32.184',                     ...
%                'DELTET/K         = 1.657D-3',                   ...
%                'DELTET/EB        = 1.671D-2',                   ...
%                'DELTET/M         = ( 6.239996 1.99096871D-7 )', ...
%                'DELTET/DELTA_AT  = ( 10, @1972-JAN-1',          ...
%                '                     11, @1972-JUL-1',          ...
%                '                     12, @1973-JAN-1',          ...
%                '                     13, @1974-JAN-1',          ...
%                '                     14, @1975-JAN-1',          ...
%                '                     15, @1976-JAN-1',          ...
%                '                     16, @1977-JAN-1',          ...
%                '                     17, @1978-JAN-1',          ...
%                '                     18, @1979-JAN-1',          ...
%                '                     19, @1980-JAN-1',          ...
%                '                     20, @1981-JUL-1',          ...
%                '                     21, @1982-JUL-1',          ...
%                '                     22, @1983-JUL-1',          ...
%                '                     23, @1985-JUL-1',          ...
%                '                     24, @1988-JAN-1',          ...
%                '                     25, @1990-JAN-1',          ...
%                '                     26, @1991-JAN-1',          ...
%                '                     27, @1992-JUL-1',          ...
%                '                     28, @1993-JUL-1',          ...
%                '                     29, @1994-JUL-1',          ...
%                '                     30, @1996-JAN-1',          ...
%                '                     31, @1997-JUL-1',          ...
%                '                     32, @1999-JAN-1 )'         ...
%                   };
%
%         %
%         % Load the kernel data into the kernel pool.
%         %
%         cspice_lmpool( textbuf )
%
%         %
%         % Ensure the loaded data exists in the kernel pool.
%         % Query the pool for each expected name, size of the
%         % variable with that name, and the type of data
%         % for that name.
%         %
%
%        [found, n, type] = cspice_dtpool( lmpoolNames );
%
%         for i = 1:numel(lmpoolNames)
%
%            if ( found(i) )
%
%               fprintf( ['Found %s, with %i values assigned' ...
%                         ' of data type %s.\n\n'],   ...
%                         char(lmpoolNames(i)), n(i), type(i) )
%
%            end
%
%         end
%
%         %
%         %  It's always good form to unload kernels after use,
%         %  particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Found DELTET/DELTA_T_A, with 1 values assigned of data type N.
%
%      Found DELTET/K, with 1 values assigned of data type N.
%
%      Found DELTET/EB, with 1 values assigned of data type N.
%
%      Found DELTET/M, with 2 values assigned of data type N.
%
%      Found DELTET/DELTA_AT, with 46 values assigned of data type N.
%
%
%      Note that the query found the five kernel variables, returned
%      the number of elements assigned to each kernel variable, and
%      the data type associated with the variable, 'N' (numerical)
%      for all cases.
%
%-Particulars
%
%   This routine allows you to store a text kernel in an internal
%   array of your program and load this array into the kernel pool
%   without first storing its contents as a text kernel.
%
%   Kernel pool variable names are restricted to a length of 32
%   characters or less.
%
%-Exceptions
%
%   1)  If any of the kernel pool variables names or their values, as
%       provided in the input `cvals' array, cannot be parsed, an error
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If there is no room left in the kernel pool to store all
%       variables present in the input `cvals' array, an error is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If the length of any kernel pool variable name present in the
%       input `cvals' array exceeds its maximum allowed length (see
%       Kernel Required Reading, kernel.req), an error is signaled by
%       a routine in the call tree of this routine.
%
%   4)  If the input argument `cvals' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   5)  If the input argument `cvals' is not of the expected type, or
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.2, 13-FEB-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 10-FEB-2010 (EDW)
%
%       Added mention of the length restriction on kernel pool variable
%       names.
%
%   -Mice Version 1.0.0, 23-FEB-2009 (EDW)
%
%-Index_Entries
%
%   Load the kernel pool from an internal text buffer
%
%-&

function cspice_lmpool( cvals )

   switch nargin
      case 1

         cvals = zzmice_str( cvals);

      otherwise

         error ( 'Usage: cspice_lmpool( _`cvals`_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('lmpool_c', cvals );
   catch spiceerr
      rethrow(spiceerr)
   end



