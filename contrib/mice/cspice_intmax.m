%-Abstract
%
%   CSPICE_INTMAX returns the value of the largest (positive) number
%   representable in an integer variable.
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
%      None.
%
%   the call:
%
%      [intmax] = cspice_intmax
%
%   returns:
%
%      intmax   the value of the largest (positive) number that can be
%               represented in an integer variable.
%
%               [1,1] = size(intmax); int32 = class(intmax)
%
%               This varies from machine to machine.
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
%   1) Obtain the integer part of a double precision number. If the
%      integer component of the number is out of range, avoid
%      overflow by making it as large or small as possible.
%
%
%      Example code begins here.
%
%
%      function intmax_ex1()
%
%         %
%         % Define a set of three numbers, two of them having an
%         % integer component that is out of range.
%         %
%         number = [2.e40, -1.5e35, 1.6]';
%
%         for i=1:3
%
%            fprintf( 'Double precision number: %8.1e\n', number(i) )
%
%            %
%            % If the integer component is out of range, avoid
%            % overflow by making it as large as possible.
%            %
%            if ( number(i) > double( cspice_intmax ) )
%
%               fprintf( '   Overflow! Greater than cspice_intmax.\n' )
%               ivalue = cspice_intmax;
%
%            elseif ( number(i) < double( cspice_intmin ) )
%
%               fprintf( '   Overflow! Smaller than cspice_intmin.\n' )
%               ivalue = cspice_intmin;
%
%            else
%
%               ivalue = floor( number(i) );
%
%            end
%
%            fprintf( '   Integer part        :  %d\n', ivalue )
%            fprintf( '\n' )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Double precision number:  2.0e+40
%         Overflow! Greater than cspice_intmax.
%         Integer part        :  2147483647
%
%      Double precision number: -1.5e+35
%         Overflow! Smaller than cspice_intmin.
%         Integer part        :  -2147483648
%
%      Double precision number:  1.6e+00
%         Integer part        :  1
%
%
%-Particulars
%
%   None.
%
%-Exceptions
%
%   Error free.
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
%   [1]  "Programming in VAX FORTRAN", Digital Equipment Corporation,
%        September 1984, Appendix C, FORTRAN Data Representation,
%        page C-2.
%
%   [2]  "Microsoft FORTRAN Reference", Microsoft Corporation
%        1989, Section 1.3.1, page 10.
%
%   [3]  "Sun FORTRAN Programmer's Guide, Sun Microsystems,
%        Revision A of 6 May 1988, Appendix F, Manual Pages for
%        FORTRAN, page 306 (RANGE).
%
%   [4]  "Language Systems FORTRAN Reference Manual", Language
%        Systems Corporation, version 1.2.1, page 3-2.
%
%   [5]  "Lahey F77L EM/32 Programmers Reference Manual",
%        version 4.0, page 95.
%
%   [6]  "FORTRAN/9000 Reference HP 9000 Series 700 Computers",
%        First Edition, June 1991, Hewlett Packard Company, page 4-4.
%
%   [7]  "SGI Fortran 77 Programmer's Guide", Document number
%        007-0711-030, page 2-2.
%
%   [8]  "Language Reference Manual", Absoft Fortran V3.2, 1993,
%        page 3-14, section 3.6.1.5. (for the NeXT)
%
%   [9]  "Unix/VMS Compatibility Libraries", Absoft Fortran V3.2,
%        1993; Chapter 3, Support Libraries, page 3-14, inmax.
%        (for the NeXT)
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 25-AUG-2021 (JDR)
%
%-Index_Entries
%
%   largest integer number
%
%-&
function [intmax] = cspice_intmax

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [intmax] = cspice_intmax' )

   end

   %
   % Call the MEX library.
   %
   try
      [intmax] = mice('intmax_c');
   catch spiceerr
      rethrow(spiceerr)
   end
